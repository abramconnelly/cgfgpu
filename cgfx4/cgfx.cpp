
#include "cgfx.hpp"
#include "app_perf.h"

#define P_PUSH(x)
#define P_POP()

/*#define P_PUSH(x)		PERF_PUSH(x); \
						for (int t=0; t < 10000; t++) \
						{

#define P_POP()			}  \
						PERF_POP(); 
	*/					

				
bool cudaCheck(CUresult status, char* msg)
{
	if (status != CUDA_SUCCESS) {
		const char* stat = "";
		cuGetErrorString(status, &stat);
		printf("CUDA ERROR: #%d %s (in %s)\n", status, stat, msg);
		exit(-1);
		return false;
	}
	return true;
}

void cudaStart (int devsel, CUdevice& dev, CUcontext& ctx, bool verbose)
{
	// NOTES:
	// CUDA and OptiX Interop: (from Programming Guide 3.8.0)
	// - CUDA must be initialized using run-time API
	// - CUDA may call driver API, but only after context created with run-time API
	// - Once app created CUDA context, OptiX will latch onto the existing CUDA context owned by run-time API
	// - Alternatively, OptiX can create CUDA context. Then set runtime API to it. (This is how Ocean SDK sample works.)

	int version = 0;
	char name[128];

	int cnt = 0;
	CUdevice dev_id;
	cuInit(0);

	//--- List devices
	cuDeviceGetCount(&cnt);
	if (cnt == 0) {
		printf("ERROR: No CUDA devices found.\n");
		dev = NULL; ctx = NULL;
		return;
	}
	if (verbose) printf("  Device List:\n");
	for (int n = 0; n < cnt; n++) {
		cuDeviceGet(&dev_id, n);
		cuDeviceGetName(name, 128, dev_id);
		if (verbose) printf("   %d. %s\n", n, name);
	}

	devsel = 0;
	if (devsel >= 0) {
		//--- Create new context with Driver API 
		cudaCheck(cuDeviceGet(&dev, devsel), "cuDeviceGet");
		cudaCheck(cuCtxCreate(&ctx, CU_CTX_SCHED_AUTO, dev), "cuCtxCreate");
	}
	cuDeviceGetName(name, 128, dev);
	if (verbose) printf("   Using Device: %d, %s, Context: %p\n", (int)dev, name, (void*)ctx);

	cuCtxSetCurrent(NULL);
	cuCtxSetCurrent(ctx);
}

float cudaGetMemUsage(bool verbose)
{
	size_t free, total;
	
	cuMemGetInfo(&free, &total);
	float freeMB = free / (1024.0f*1024.0f);		// MB
	float totalMB = total / (1024.0f*1024.0f);
	
	if (verbose) printf("GPU memory: %5.2fMB used, %5.2MB total\n", totalMB - freeMB, totalMB);
	return totalMB - freeMB;		// used
}


void CGFX_t::Initialize(int tile_res, int tilesteps_per_genome)
{
	mbVerbose = false;
	mNumGenome = 0;
	mMaxGenome = 0;
	mTilesPerBlock = 32;
	mTagSize = 24;
	mTileRes = 201;	
	mTilesPerGenome = 10669088;
	mTilepathPerGenome = 862;
	mBlocksPerGenome = mTilesPerGenome / 32;

	// Start CUDA
	cudaStart(0, mCuDevice, mCuContext, mbVerbose);

	// Cache offset optimization table
	if (mbVerbose) printf("Building table.\n");
	BuildTable();

	// Load CUDA kernels
	if (mbVerbose) printf("Loading CUDA kernels.\n");
	cudaCheck(cuModuleLoad(&mCuModule, "kernel.ptx"), "cuModuleLoad");
	cudaCheck(cuModuleGetFunction(&mCuConcordanceKernel, mCuModule, "concordanceKernel"), "cuModuleGetFunction (concordance)");
}

void CGFX_t::Resize (int num)
{	
	if (num >= mMaxGenome) {

		int max = mMaxGenome + 4; 
		mBitvecSize = mBlocksPerGenome * sizeof(uint32_t);		// # Blocks * 32 bits/block
		mOveroffSize = mTilepathPerGenome * sizeof(uint64_t);	// # Tilepath * 64 bits/path

		mm.Resize(mCanonVec, max, mBitvecSize);
		mm.Resize(mCacheVec, max, mBitvecSize);
		mm.Resize(mLoqVec, max, mBitvecSize);
		mm.Resize(mSpanVec, max, mBitvecSize);
		mm.Resize(mOverflowOffsetVec, max, mOveroffSize );
		mm.Resize(mOverflowAllocVec, max, sizeof(uint64_t));

		mMaxGenome = num;
	}
}
void CGFX_t::InitializeGlobalData (cgf_t* cgf)
{
	// Prepare TileMap
	//
	tilemap_t tilemap;
	str2tilemap(cgf->TileMap, &tilemap);			// expand string to tilemap array	
	int tilemap_cnt = tilemap.offset.size();

	mm.Resize(mTileMap, tilemap_cnt, 3 * sizeof(uint32_t));

	uint32_t* tmap = (uint32_t*)mm.getCPU(mTileMap, 0);

	for (int n = 0; n < tilemap_cnt; n++)
	{
		*tmap++ = tilemap.offset[n];
		*tmap++ = tilemap.variant[0][n];
		*tmap++ = tilemap.variant[1][n];
	}
	mm.Commit(mTileMap);

	// Prepare Block-to-Tilepath Mapping
	//
	mm.Resize(mBlockToTilepath, mBlocksPerGenome, 2 * sizeof(uint32_t));

	uint32_t* dat = (uint32_t*) mm.getCPU(mBlockToTilepath, 0);

	int start, blk = 0;
	for (int p = 0; p < cgf->TilePathCount; p++) {
		int offsprev = (p == 0) ? 0 : cgf->StrideOffset[p - 1];
		int blkcnt = cgf->StrideOffset[p] - offsprev;
		start = offsprev;
		for (int t = 0; t < blkcnt; t++) {
			*dat++ = p;
			*dat++ = start;
			blk++;
			if (blk > mBlocksPerGenome) {
				printf("ERROR: Genome tilepath and step count does not match global number of blocks.\n");
			}
		}
	}
	mm.Commit( mBlockToTilepath );	
}

void CGFX_t::LoadCgf4 (char* fname)
{
	if (mbVerbose) printf("Loading %s.\n", fname );
	FILE* ifp;
	if ( (ifp = fopen( fname, "rb")) == NULL) { 
		printf ( "ERROR: Cannot read %s\n", fname); 
		return; 
	}
	cgf_t* cgf = NULL;

	Resize ( mNumGenome + 1 );

	// Read cgf4 genome	
	//cgf = cgf_read(ifp);
	cgf = cgf_read_hiq (ifp);
	mCgfVec.push_back( cgf );
	fclose(ifp);

	// Initialize Globals
	if (mBlockToTilepath.cpu == 0x0) {
		InitializeGlobalData ( cgf );		
	}

	// Allocate space in cohort
	int id, chk[10];
	id = mm.Alloc(mCanonVec);
	chk[1] = mm.Alloc(mCacheVec);
	chk[2] = mm.Alloc(mLoqVec);
	chk[3] = mm.Alloc(mSpanVec);
	chk[4] = mm.Alloc(mOverflowOffsetVec);
	chk[5] = mm.Alloc(mOverflowAllocVec);
	assert ( (id == chk[1]) && (id == chk[2]) && (id == chk[3]) && (id==chk[4]) && (id==chk[5]) );

	// Transfer cgf4 bit-vectors into cohort
	memcpy ( mm.getCPU(mCanonVec, id),	&cgf->Canon[0], mBitvecSize);
	memcpy ( mm.getCPU(mCacheVec, id),	&cgf->CacheOverflow[0], mBitvecSize);
	memcpy ( mm.getCPU(mSpanVec, id),	&cgf->Span[0], mBitvecSize);
	memcpy ( mm.getCPU(mLoqVec, id),	&cgf->Loq[0], mBitvecSize);

	// Transfer overflow offsets
	memcpy ( mm.getCPU(mOverflowOffsetVec, id), &cgf->OverflowOffset[0], mOveroffSize);

	// Transfer overflow	
	uint64_t sz_overflow = sizeof(uint16_t) * cgf->OverflowOffset[cgf->OverflowOffset.size() - 1];
	char* overflow = (char*) &cgf->Overflow[0];
	char* cpu_data = (char*) malloc( sz_overflow );							// alloc overflow data
	char** cpu_addr = &cpu_data;											// pointer to address
	mm.SetDataCPU (mOverflowAllocVec, id, 0, cpu_addr, sizeof(uint64_t));	// add to pointer table
	memcpy(cpu_data, overflow, sz_overflow);								// transfer data 

	CUdeviceptr gpu_data;
	cuMemAlloc (&gpu_data, sz_overflow);									// alloc overflow data
	char** gpu_addr = (char**) gpu_data;									// pointer to address
	mm.SetDataGPU(mOverflowAllocVec, id, 0, gpu_addr, sizeof(uint64_t));	// add to pointer table
	cuMemcpyHtoD(gpu_data, overflow, sz_overflow);							// transfer data

	// Commit vectors to GPU
	mm.Commit(mCanonVec);
	mm.Commit(mCacheVec);
	mm.Commit(mSpanVec);
	mm.Commit(mLoqVec);
	mm.Commit(mOverflowOffsetVec);

	mNumGenome++;
}

void CGFX_t::SetConcordanceRange(int start_path, int start_step, int end_path, int end_step)
{
	mStartPath = start_path;
	mStartStep = start_step;
	mEndPath = end_path;
	mEndStep = end_step;
}

void CGFX_t::SetConcordanceOutput(int num_concord )
{
	// Concordance output vectors
	mm.Create(mOutMatchList, num_concord, mBlocksPerGenome);
	mm.Create(mOutOverflowList, num_concord, mBlocksPerGenome);
	mm.Create(mOutTotalList, num_concord, mBlocksPerGenome);
}

const int tab32[32] = {
	0,   9, 1,  10, 13, 21,  2, 29,
	11, 14, 16, 18, 22, 25,  3, 30,
	8,  12, 20, 28, 15, 17, 24,  7,
	19, 27, 23,  6, 26,  5,  4, 31 };

inline int log2_32(uint32_t value) {
	value |= value >> 1;
	value |= value >> 2;
	value |= value >> 4;
	value |= value >> 8;
	value |= value >> 16;
	return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

void CGFX_t::BuildTable()
{
	mm.Create( mTable, 256, 9 );
	uint8_t* table = (uint8_t*) mTable.cpu;

	uint8_t k, p;
	uint8_t* entry;

	for (int i = 0; i < 256; i++) {
		entry = &( table[9 * i] );
		for (p = 0, k = 0; k < 8; k++) {
			if ((0x1u << (uint32_t)k) & (uint8_t) i) {
				entry[p++] = k;
			}
		}
		entry[8] = p;
	}
	mm.Commit(mTable);
}
#define NULLMIN		-2147483647
#define NULLMAX		 2147483646


void CGFX_t::ConcordanceOverflowOnBlock(int cgfa, int cgfb, uint64_t ii, int blk_start, int blk_end,
	uchar& blk_match, uchar& blk_tot)
{
	// input reads
	uint32_t hiq_a, hiq_b;
	uint32_t span_a, span_b;
	uint32_t canon_a, canon_b;
	uint32_t cache_a, cache_b;

	// computed
	uint32_t xspan_a, xspan_b;
	uint32_t actual_canon_a, actual_canon_b;
	uint32_t noncanon_a, noncanon_b;
	uint32_t anchor_mask_a, anchor_mask_b;
	uint32_t non_anchor_span_mask_a, non_anchor_span_mask_b;
	uint32_t tmp;

	int hexit_a[8], hexit_b[8],
		hexit_relative_step_a[32], hexit_relative_step_b[32];
	int* hout;
	uint8_t* table_entry;

	int i, j, p, k;

	// Input bit vectors
	hiq_a = ~*(((uint32_t*)mm.getCPU(mLoqVec, cgfa)) + ii);		// flip loq_a
	hiq_b = ~*(((uint32_t*)mm.getCPU(mLoqVec, cgfb)) + ii);		// flip loq_b
	span_a = *(((uint32_t*)mm.getCPU(mSpanVec, cgfa)) + ii);		// span bits
	span_b = *(((uint32_t*)mm.getCPU(mSpanVec, cgfb)) + ii);
	canon_a = *(((uint32_t*)mm.getCPU(mCanonVec, cgfa)) + ii);		// canon bits
	canon_b = *(((uint32_t*)mm.getCPU(mCanonVec, cgfb)) + ii);
	cache_a = *(((uint32_t*)mm.getCPU(mCacheVec, cgfa)) + ii);		// cache bits
	cache_b = *(((uint32_t*)mm.getCPU(mCacheVec, cgfb)) + ii);

	// Computed bit-vectors
	xspan_a = ~span_a;
	actual_canon_a = canon_a & xspan_a & hiq_a;		// get actual canon bits, masking off the span bit rules
	anchor_mask_a = span_a & hiq_a & (~canon_a);	// anchor tile bit vector for convenience.
	noncanon_a = (anchor_mask_a & hiq_a) | ((~span_a) & (~actual_canon_a) & hiq_a); // non-canonical tiles and anchor tiles, excluding spanning tiles.
	non_anchor_span_mask_a = span_a & (~anchor_mask_a);

	xspan_b = ~span_b;
	actual_canon_b = canon_b & xspan_b & hiq_b;
	anchor_mask_b = span_b & hiq_b & (~canon_b);
	noncanon_b = (anchor_mask_b & hiq_b) | ((~span_b) & (~actual_canon_b) & hiq_b);
	non_anchor_span_mask_b = span_b & (~anchor_mask_b);

	// Count the number of canonical matches.
	blk_match = NumberOfSetBits32(actual_canon_a & actual_canon_b);

	// Cache values for genome a, pulled out into array 	
	hexit_a[0] = (cache_a & ((uint32_t)0xf));
	hexit_a[1] = (cache_a & ((uint32_t)0xf << 4)) >> 4;
	hexit_a[2] = (cache_a & ((uint32_t)0xf << 8)) >> 8;
	hexit_a[3] = (cache_a & ((uint32_t)0xf << 12)) >> 12;
	hexit_a[4] = (cache_a & ((uint32_t)0xf << 16)) >> 16;
	hexit_a[5] = (cache_a & ((uint32_t)0xf << 20)) >> 20;
	hexit_a[6] = (cache_a & ((uint32_t)0xf << 24)) >> 24;
	hexit_a[7] = (cache_a & ((uint32_t)0xf << 28)) >> 28;

	// Fill out the hexit_relative_step array
	// This is a count of the canonical bits to get the cache entry positions.
	//
	tmp = noncanon_a;
	hout = &hexit_relative_step_a[0];
	for (int jmp = 0; jmp < 32; jmp += 8) {
		table_entry = (uint8_t*)mTable.cpu + 9 * (tmp & 0xff);
		hout[0] = table_entry[0] + jmp;
		hout[1] = table_entry[1] + jmp;
		hout[2] = table_entry[2] + jmp;
		hout[3] = table_entry[3] + jmp;
		hout[4] = table_entry[4] + jmp;
		hout[5] = table_entry[5] + jmp;
		hout[6] = table_entry[6] + jmp;
		hout[7] = table_entry[7] + jmp;
		hout += (int)table_entry[8];
		tmp >>= 8;
	}

	/*tmp = overflow_a;
	for (k = 0; (k<8) && (tmp != 0); k++, tmp &= (tmp - 1)) {
	hexit_relative_step_a[k] = log2_32(tmp & ~(tmp - 1));
	}*/

	// Cache values for genome b, pulled out into array 
	//
	hexit_b[0] = (cache_b & ((uint32_t)0xf));
	hexit_b[1] = (cache_b & ((uint32_t)0xf << 4)) >> 4;
	hexit_b[2] = (cache_b & ((uint32_t)0xf << 8)) >> 8;
	hexit_b[3] = (cache_b & ((uint32_t)0xf << 12)) >> 12;
	hexit_b[4] = (cache_b & ((uint32_t)0xf << 16)) >> 16;
	hexit_b[5] = (cache_b & ((uint32_t)0xf << 20)) >> 20;
	hexit_b[6] = (cache_b & ((uint32_t)0xf << 24)) >> 24;
	hexit_b[7] = (cache_b & ((uint32_t)0xf << 28)) >> 28;

	// Fill out the hexit_relative_step array
	// This is a count of the canonical bits to get the cache entry positions.
	//
	tmp = noncanon_b;
	hout = &hexit_relative_step_b[0];
	for (int jmp = 0; jmp < 32; jmp += 8) {
		table_entry = (uint8_t*)mTable.cpu + 9 * (tmp & 0xff);
		hout[0] = table_entry[0] + jmp;
		hout[1] = table_entry[1] + jmp;
		hout[2] = table_entry[2] + jmp;
		hout[3] = table_entry[3] + jmp;
		hout[4] = table_entry[4] + jmp;
		hout[5] = table_entry[5] + jmp;
		hout[6] = table_entry[6] + jmp;
		hout[7] = table_entry[7] + jmp;
		hout += (int)table_entry[8];
		tmp >>= 8;
	}

	/*tmp = overflow_b;
	for (k = 0; (k<8) && (tmp != 0); k++, tmp &= (tmp - 1)) {
	hexit_relative_step_b[k] = log2_32(tmp & ~(tmp - 1));
	}*/


	// Cache value testing
	//
	// Do a zipper match to count the number of cache overflow hits.
	// Tile variants that overflow from the cache will be picked
	// up by the overflow count.
	//
	for (i = 0, j = 0; (i<8) && (j<8); ) {

		if ((hexit_relative_step_a[i] < 0) || (hexit_relative_step_b[j] < 0)) { break; }

		if (hexit_relative_step_a[i] < hexit_relative_step_b[j]) {
			i++; continue;
		}
		if (hexit_relative_step_a[i] > hexit_relative_step_b[j]) {
			j++; continue;
		}
		if (hexit_relative_step_a[i] == hexit_relative_step_b[j]) {
			if ((hexit_a[i] > 0) && (hexit_a[i] < 0xf) &&
				(hexit_b[j] > 0) && (hexit_b[j] < 0xf) &&
				(hexit_a[i] == hexit_b[j])) {

				blk_match++;
			}
			i++;
			j++;
		}
	}
	// The total number of high quality matches is the number of high quality
	// tiles both have in common, excluding the non-anchor high quality tiles
	//
	blk_tot = NumberOfSetBits32(hiq_a & hiq_b & (~non_anchor_span_mask_a) & (~non_anchor_span_mask_b));


	// OVERFLOW
	//
	int over_aL[32], over_aR[32];
	int over_bL[32], over_bR[32];
	uint16_t *overptr;
	int L, R, M;
	bool chk;

	for (int n = 0; n < 32; n++) {
		over_aL[n] = NULLMIN;
		over_aR[n] = NULLMIN;
		over_bL[n] = NULLMAX;
		over_bR[n] = NULLMAX;
	}
	//printf("blk: %d, tilestart: %d\n", ii, ii * 32);

	// Gather cache spill overflow B
	int variant;
	uint32_t* tmap = (uint32_t*) mm.getCPU(mTileMap, 0);
	if (i == 8) {
		for (; j < 8; j++) {
			M = hexit_relative_step_b[j];
			variant = hexit_b[j];
			if (M < 0 || M > 31) break;
			if (variant <= 0 || variant >= 0xf) { continue; }
			// variants for each knot position
			for (k = tmap[variant * 3 - 3]; k < tmap[variant * 3] && M < 32; k++) {
				assert(M >= 0 && M < 32);
				over_bL[ M ] = tmap[ k*3 + 1];
				over_bR[ M ] = tmap[ k*3 + 2];
				M++;
			}
		}
	}
	// Gather cache spill overflow A
	if (j == 8) {
		for (; i < 8; i++) {
			M = hexit_relative_step_a[i];
			variant = hexit_a[i];
			if (M < 0 || M > 31) break;
			if (variant <= 0 || variant >= 0xf) { continue; }
			// variants for each knot position
			for (k = tmap[variant * 3 - 3]; k < tmap[variant * 3] && M < 32; k++) {
				assert(M >= 0 && M < 32);
				over_aL[ M ] = tmap[k * 3 + 1];
				over_aR[ M ] = tmap[k * 3 + 2];
				M++;
			}
		}
	}

	// Get tilepath and tilestart for current block (using Block-to-Tilepath map)
	uint64_t* map = (uint64_t*) mm.getCPU(mBlockToTilepath, ii);	// global lookup indexed by block id
	int tilepath = (*map) & 0xFFFFFFFF;		// low  32-bit word is path
	int tilestart = (*map) >> 32;			// high 32-bit word is starting block of path
	tilestart *= 32;		// convert block index to absolute tile index

	if (ii == 181640)
		chk = true;

	// Overflow for A
	// Find start and end of overflow for entire tilepath in cgf'a'
	map = ((uint64_t*) mm.getCPU(mOverflowOffsetVec, cgfa)) + tilepath;	
	overptr = *(uint16_t**)mm.getCPU(mOverflowAllocVec, cgfa);
	L = *(map - 1);		// tilepath_start_offset
	R = *map;			// tilepath_start_offnxt
	// Binary search for overflow starting in this block	
	while (L < R)
	{
		M = floor(((L + R) / 6)) * 3;
		if (tilestart + *(overptr + M) < ii*32) 
			L = M + 3;
		else
			R = M;
	}		
	// Gather all overflow A in this block
	if (tilestart + *(overptr + L) < (ii + 1) * 32) {
		for (R = L; tilestart + *(overptr + R) < (ii + 1) * 32 && R < *map; R += 3) {
			M = tilestart + *(overptr + R) - (ii * 32);
			assert(M >= 0 && M < 32);
			over_aL[M] = *(overptr + R + 1);
			over_aR[M] = *(overptr + R + 2);
			//printf(" overA: %d, %d:%d\n", tilestart + *(overptr + R), *(overptr + R + 1), *(overptr + R + 2));
		}
	}
	
	// Overflow for B
	// Find start and end of overflow for entire tilepath in cgf'b'
	map = ((uint64_t*) mm.getCPU(mOverflowOffsetVec, cgfb)) + tilepath;
	overptr = *(uint16_t**) mm.getCPU(mOverflowAllocVec, cgfb);
	L = *(map - 1);		// tilepath_start_offset
	R = *map;			// tilepath_start_offnxt
	// Binary search for overflow starting in this block
	while (L < R)
	{
		M = floor(((L + R) / 6)) * 3;
		if (tilestart + *(overptr + M) < ii * 32)
			L = M + 3;
		else
			R = M;
	}
	// Gather all overflow B in this block
	if (tilestart + *(overptr + L) < (ii + 1) * 32) {
		for (R = L; tilestart + *(overptr + R) < (ii + 1) * 32 && (R < *map); R += 3) {
			M = tilestart + *(overptr + R) - (ii * 32);
			assert(M >= 0 && M < 32);
			over_bL[M] = *(overptr + R + 1);
			over_bR[M] = *(overptr + R + 2);
			//printf(" overB: %d, %d:%d\n", tilestart + *(overptr + R), *(overptr + R + 1), *(overptr + R + 2));
		}
	}

	// Count overflow matches
	for (int n = 0; n < 32; n++)
		if (over_aL[n] == over_bL[n] && over_aR[n] == over_bR[n]) {
			//printf(" MATCH: %d, %d:%d\n", n, over_aL[n], over_aR[n]);
			blk_match++;
		}

	chk = true;

	//int cnt_over_a 
	//int cnt_over_b
	//blk_match = 0;
	//blk_tot = 1;
}

void CGFX_t::ConcordanceOnBlock( int cgfa, int cgfb, uint64_t ii, int blk_start, int blk_end, 
									uchar& blk_match, uchar& blk_tot )
{
	// input reads
	uint32_t hiq_a, hiq_b;
	uint32_t span_a, span_b;	
	uint32_t canon_a, canon_b;	
	uint32_t cache_a, cache_b;
	
	// computed
	uint32_t xspan_a, xspan_b;
	uint32_t actual_canon_a, actual_canon_b;
	uint32_t overflow_a, overflow_b;
	uint32_t anchor_mask_a, anchor_mask_b;
	uint32_t non_anchor_span_mask_a, non_anchor_span_mask_b;
	uint32_t tmp;

	int hexit_a[8], hexit_b[8],
		hexit_relative_step_a[32], hexit_relative_step_b[32];
	int* hout;
	uint8_t* table_entry;

	int i, j, p, k;

	// Input bit vectors
	hiq_a = ~ *(((uint32_t*) mm.getCPU(mLoqVec, cgfa)) + ii);		// flip loq_a
	hiq_b = ~ *(((uint32_t*) mm.getCPU(mLoqVec, cgfb)) + ii);		// flip loq_b
	span_a = *(((uint32_t*) mm.getCPU(mSpanVec, cgfa)) + ii);		// span bits
	span_b = *(((uint32_t*) mm.getCPU(mSpanVec, cgfb)) + ii);	
	canon_a = *(((uint32_t*) mm.getCPU(mCanonVec, cgfa)) + ii);		// canon bits
	canon_b = *(((uint32_t*) mm.getCPU(mCanonVec, cgfb)) + ii);
	cache_a = *(((uint32_t*) mm.getCPU(mCacheVec, cgfa)) + ii);		// cache bits
	cache_b = *(((uint32_t*) mm.getCPU(mCacheVec, cgfb)) + ii);
	
	// Computed bit-vectors
	xspan_a = ~span_a;
	actual_canon_a = canon_a & xspan_a & hiq_a;		// get actual canon bits, masking off the span bit rules
	anchor_mask_a = span_a & hiq_a & (~canon_a);	// anchor tile bit vector for convenience.
	overflow_a = (anchor_mask_a & hiq_a) | ((~span_a) & (~actual_canon_a) & hiq_a); // non-canonical tiles and anchor tiles, excluding spanning tiles.
	non_anchor_span_mask_a = span_a & (~anchor_mask_a);

	xspan_b = ~span_b;
	actual_canon_b = canon_b & xspan_b & hiq_b;
	anchor_mask_b = span_b & hiq_b & (~canon_b);
	overflow_b = (anchor_mask_b & hiq_b) | ((~span_b) & (~actual_canon_b) & hiq_b);
	non_anchor_span_mask_b = span_b & (~anchor_mask_b);

	// Count the number of canonical matches.
	blk_match = NumberOfSetBits32(actual_canon_a & actual_canon_b);

	// Cache values for genome a, pulled out into array 	
	hexit_a[0] = (cache_a & ((uint32_t)0xf));
	hexit_a[1] = (cache_a & ((uint32_t)0xf << 4)) >> 4;
	hexit_a[2] = (cache_a & ((uint32_t)0xf << 8)) >> 8;
	hexit_a[3] = (cache_a & ((uint32_t)0xf << 12)) >> 12;
	hexit_a[4] = (cache_a & ((uint32_t)0xf << 16)) >> 16;
	hexit_a[5] = (cache_a & ((uint32_t)0xf << 20)) >> 20;
	hexit_a[6] = (cache_a & ((uint32_t)0xf << 24)) >> 24;
	hexit_a[7] = (cache_a & ((uint32_t)0xf << 28)) >> 28;
	
	// Fill out the hexit_relative_step array
	// This is a count of the canonical bits to get the cache entry positions.
	//
	tmp = overflow_a;
	hout = &hexit_relative_step_a[0];
	for (int jmp = 0; jmp < 32; jmp += 8) {
		table_entry = (uint8_t*) mTable.cpu + 9 * (tmp & 0xff);
		hout[0] = table_entry[0] + jmp;
		hout[1] = table_entry[1] + jmp;
		hout[2] = table_entry[2] + jmp;
		hout[3] = table_entry[3] + jmp;
		hout[4] = table_entry[4] + jmp;
		hout[5] = table_entry[5] + jmp;
		hout[6] = table_entry[6] + jmp;
		hout[7] = table_entry[7] + jmp;
		hout += (int)table_entry[8];
		tmp >>= 8;
	}

	/*tmp = overflow_a;
	for (k = 0; (k<8) && (tmp != 0); k++, tmp &= (tmp - 1)) {
		hexit_relative_step_a[k] = log2_32(tmp & ~(tmp - 1));
	}*/

	// Cache values for genome b, pulled out into array 
	//
	hexit_b[0] = (cache_b & ((uint32_t)0xf));
	hexit_b[1] = (cache_b & ((uint32_t)0xf << 4)) >> 4;
	hexit_b[2] = (cache_b & ((uint32_t)0xf << 8)) >> 8;
	hexit_b[3] = (cache_b & ((uint32_t)0xf << 12)) >> 12;
	hexit_b[4] = (cache_b & ((uint32_t)0xf << 16)) >> 16;
	hexit_b[5] = (cache_b & ((uint32_t)0xf << 20)) >> 20;
	hexit_b[6] = (cache_b & ((uint32_t)0xf << 24)) >> 24;
	hexit_b[7] = (cache_b & ((uint32_t)0xf << 28)) >> 28;

	// Fill out the hexit_relative_step array
	// This is a count of the canonical bits to get the cache entry positions.
	//
	tmp = overflow_b;
	hout = &hexit_relative_step_b[0];
	for (int jmp = 0; jmp < 32; jmp += 8) {
		table_entry = (uint8_t*) mTable.cpu + 9 * (tmp & 0xff);
		hout[0] = table_entry[0] + jmp;
		hout[1] = table_entry[1] + jmp;
		hout[2] = table_entry[2] + jmp;
		hout[3] = table_entry[3] + jmp;
		hout[4] = table_entry[4] + jmp;
		hout[5] = table_entry[5] + jmp;
		hout[6] = table_entry[6] + jmp;
		hout[7] = table_entry[7] + jmp;
		hout += (int)table_entry[8];
		tmp >>= 8;
	}

	/*tmp = overflow_b;
	for (k = 0; (k<8) && (tmp != 0); k++, tmp &= (tmp - 1)) {
		hexit_relative_step_b[k] = log2_32(tmp & ~(tmp - 1));
	}*/


	// Cache value testing
	//
	// Do a zipper match to count the number of cache overflow hits.
	// Tile variants that overflow from the cache will be picked
	// up by the overflow count.
	//
	for (i = 0, j = 0; (i<8) && (j<8); ) {

		if ((hexit_relative_step_a[i] < 0) || (hexit_relative_step_b[j] < 0)) { break; }

		if (hexit_relative_step_a[i] < hexit_relative_step_b[j]) {
			i++; continue;
		}
		if (hexit_relative_step_a[i] > hexit_relative_step_b[j]) {
			j++; continue;
		}
		if (hexit_relative_step_a[i] == hexit_relative_step_b[j]) {
			if ((hexit_a[i] > 0) && (hexit_a[i] < 0xf) &&
				(hexit_b[j] > 0) && (hexit_b[j] < 0xf) &&
				(hexit_a[i] == hexit_b[j])) {

				if ((ii > blk_start) && (ii < blk_end)) {
					blk_match++;
				}

				else if ((ii == blk_start) && (ii < blk_end) &&
					((mStartTile + hexit_relative_step_a[i]) >= mStartStep) &&
					((mStartTile + hexit_relative_step_b[j]) >= mStartStep)) {
					blk_match++;
				}

				else if ((ii > blk_start) && (ii == blk_end) &&
					((mStartTile + hexit_relative_step_a[i]) <= mEndStep) &&
					((mStartTile + hexit_relative_step_b[j]) <= mEndStep)) {
					blk_match++;
				}

				else if ((ii == blk_start) && (ii == blk_end) &&
					((mStartTile + hexit_relative_step_a[i]) >= mStartStep) &&
					((mStartTile + hexit_relative_step_b[j]) >= mStartStep) &&
					((mStartTile + hexit_relative_step_a[i]) <= mEndStep) &&
					((mStartTile + hexit_relative_step_b[j]) <= mEndStep)) {
					blk_match++;
				}
				// else we skip over and don't count the match as it falls outside of the selected window
			}
			i++;
			j++;
		}
	}


	// The total number of high quality matches is the number of high quality
	// tiles both have in common, excluding the non-anchor high quality tiles
	//
	blk_tot = NumberOfSetBits32(hiq_a & hiq_b & (~non_anchor_span_mask_a) & (~non_anchor_span_mask_b));
}

unsigned int nth_bit_set(uint32_t value, unsigned int n)
{
	const uint32_t  pop2 = (value & 0x55555555u) + ((value >> 1) & 0x55555555u);
	const uint32_t  pop4 = (pop2 & 0x33333333u) + ((pop2 >> 2) & 0x33333333u);
	const uint32_t  pop8 = (pop4 & 0x0f0f0f0fu) + ((pop4 >> 4) & 0x0f0f0f0fu);
	const uint32_t  pop16 = (pop8 & 0x00ff00ffu) + ((pop8 >> 8) & 0x00ff00ffu);
	const uint32_t  pop32 = (pop16 & 0x000000ffu) + ((pop16 >> 16) & 0x000000ffu);
	unsigned int    rank = 0;
	unsigned int    temp;

	if (n++ >= pop32)
		return 32;

	temp = pop16 & 0xffu;
	/* if (n > temp) { n -= temp; rank += 16; } */
	rank += ((temp - n) & 256) >> 4;
	n -= temp & ((temp - n) >> 8);

	temp = (pop8 >> rank) & 0xffu;
	/* if (n > temp) { n -= temp; rank += 8; } */
	rank += ((temp - n) & 256) >> 5;
	n -= temp & ((temp - n) >> 8);

	temp = (pop4 >> rank) & 0x0fu;
	/* if (n > temp) { n -= temp; rank += 4; } */
	rank += ((temp - n) & 256) >> 6;
	n -= temp & ((temp - n) >> 8);

	temp = (pop2 >> rank) & 0x03u;
	/* if (n > temp) { n -= temp; rank += 2; } */
	rank += ((temp - n) & 256) >> 7;
	n -= temp & ((temp - n) >> 8);

	temp = (value >> rank) & 0x01u;
	/* if (n > temp) rank += 1; */
	rank += ((temp - n) & 256) >> 8;

	return rank;
}

int CGFX_t::Concordance (int cgfa, int cgfb, int outslot)
{
	int blk_start, blk_end;
	uchar blk_match, blk_over_match, blk_total;
	int match=0, over_match=0, total=0;
	
	// Ranges
	cgf_get_block_start_end(mCgfVec[cgfa], blk_start, blk_end, mStartPath, mStartStep, mEndPath, mEndStep);
	mStartTile = (mStartStep / (8 * mCgfVec[cgfa]->Stride)) * 32;

	uchar* out_match = (uchar*) mm.getCPU(mOutMatchList, outslot);
	uchar* out_over_match = (uchar*)mm.getCPU(mOutOverflowList, outslot);
	uchar* out_total  = (uchar*)mm.getCPU(mOutTotalList, outslot);

	if (out_match==0x0) {
		printf("ERROR: No output vector created. Call SetConcordanceOutput.\n");
		return 0;
	}

	// Concordance, no overflow
	PERF_PUSH("ConcNoOver");
	clock_t t = clock();
	for (uint64_t blk = blk_start; blk <= blk_end; blk++ ) {
		ConcordanceOnBlock ( cgfa, cgfb, blk, blk_start, blk_end, blk_match, blk_total );
		out_match[blk - blk_start] = blk_match;
		out_total[blk - blk_start] = blk_total;
		match += blk_match;
		total += blk_total;
	}	
	t = clock() - t;
	PERF_POP();
	
	if (mbVerbose) {
		printf("Match: %3.2f%% (%d / %d), t: %4.1f ms ---- CGFX cpu, no overflow\n", match*100.0f / total, match, total, (float)t * 1000.0f / CLOCKS_PER_SEC);
		for (int n = mStartView; n < mStartView + 25; n++)
			printf("%02d ", out_match[n]);
		printf("\n\n");
	}

	// Concordance, OVERFLOW
	PERF_PUSH("ConcOverflow");
	t = clock();
	for (uint64_t blk = blk_start; blk <= blk_end; blk++) {
		ConcordanceOverflowOnBlock(cgfa, cgfb, blk, blk_start, blk_end, blk_over_match, blk_total);
		out_over_match[blk - blk_start] = blk_over_match;
		over_match += blk_over_match;
	}
	t = clock() - t;
	PERF_POP();
	
	if (mbVerbose) {
		printf("Match: %3.2f%% (%d / %d), t: %4.1f ms ---- CGFX cpu, OVERFLOW\n", over_match*100.0f / total, over_match, total, (float)t * 1000.0f / CLOCKS_PER_SEC);
		for (int n = mStartView; n < mStartView + 25; n++)
			printf("%02d ", out_over_match[n] );
		printf("\n\n");
	}

	return (int) match * 100 / total;
}

int CGFX_t::ConcordanceGPU ( int cgfa )
{
	int blk_start, blk_end;
	uchar blk_match, blk_tot;

	// Ranges
	cgf_get_block_start_end(mCgfVec[cgfa], blk_start, blk_end, mStartPath, mStartStep, mEndPath, mEndStep);
	mStartTile = (mStartStep / (8 * mCgfVec[cgfa]->Stride)) * 32;

	// Determine GPU grid dims
	int num_blocks = blk_end - blk_start + 1;			// number of blocks to concord
	int total_blocks = mBlocksPerGenome;				// total blocks in cgf data
	int thread_x = 64;
	int thread_y = 8;
	int grid_x = int((num_blocks - 1) / thread_x) + 1;	// gpu grid to launch
	int grid_y = int((mNumGenome - 1) / thread_y) + 1;

	// Launch kernel on blocks	
	void* args[17] = {
		&cgfa, 
		&blk_start, &blk_end,
		&mStartPath, &mStartStep, &mEndPath, &mEndStep, &mStartTile,
		&total_blocks, &mNumGenome,
		&mCanonVec.gpu, &mCacheVec.gpu, &mSpanVec.gpu, &mLoqVec.gpu,
		&mOutMatchList.gpu, &mOutTotalList.gpu, &mTable.gpu
	};
	PERF_PUSH("ConcordanceGPU");
	clock_t t = clock();
	cudaCheck( cuLaunchKernel(mCuConcordanceKernel, grid_x, grid_y, 1, thread_x, thread_y, 1, 0, NULL, args, NULL), "cuLaunchKernel (concordance)");
	cuCtxSynchronize();
	t = clock() - t;
	PERF_POP();
	
	mm.Retrieve(mOutMatchList);
	mm.Retrieve(mOutTotalList);

	// Display match results
	if (mbVerbose) {
		for (int j = 1; j < mNumGenome; j++) {
			uchar* out_matches = (uchar*)mm.getCPU(mOutMatchList, j);	// CPU values now valid
			uchar* out_totals = (uchar*)mm.getCPU(mOutTotalList, j);

			int match = 0, total = 0;
			for (int n = 0; n < total_blocks; n++) {
				match += out_matches[n];
				total += out_totals[n];
			}
			printf("Match: %3.2f%% (%d / %d), t: %4.1f ms ---- CGFX GPU, no overflow\n", match*100.0f/total, match, total, (float)t * 1000.0f / CLOCKS_PER_SEC );
			for (int n = mStartView; n < mStartView+25; n++)
				printf("%02d ", out_matches[n] );
			printf("\n\n");

		}
	}

	cudaGetMemUsage( mbVerbose );

	return 1;
}
