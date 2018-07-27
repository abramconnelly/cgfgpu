


#include "cgf4.hpp"
#include "cgfx_dataptr.hpp"

extern int cgfx_main(int argc, char **argv);


class CGFX_t {
public:

	void Initialize(int tile_res, int tilesteps_per_genome);
	void InitializeGlobalData (cgf_t* cgf);
	void Resize (int max);
	void BuildTable();
	int getMatchSize() { return mOutMatchList.stride; }
	uchar* getMatchList(int i) { return (uchar*) mOutMatchList.cpu + i*mOutMatchList.stride; }
	uchar* getTotalList(int i) { return (uchar*) mOutTotalList.cpu + i*mOutTotalList.stride; }

	void LoadCgf4(char* fname);
	void SetConcordanceRange(int start_path, int start_step, int end_path, int end_step);
	void SetConcordanceOutput(int n);
	int Concordance(int ndx_a, int ndx_b, int outslot);
	void ConcordanceOnBlock(int cgfa, int cgfb, uint64_t ii, int blk_start, int blk_end, uchar& blk_match, uchar& blk_tot);
	void ConcordanceOverflowOnBlock(int cgfa, int cgfb, uint64_t ii, int blk_start, int blk_end, uchar& blk_match, uchar& blk_tot);
	int ConcordanceGPU ( int ndx_a );
	
	void SetVerbose(bool tf)	{ mbVerbose = tf; }
	void SetStartView(int v)	{ mStartView = v; }
	int getNumGenome()			{ return mNumGenome;  }

	int			mNumGenome;
	int			mMaxGenome;

	int			mTilesPerBlock;			//   fixed:  32 tiles per block
	int			mTagSize;				//   fixed:  24 base pairs
	int			mTileRes;				// default: 201 base pairs	
	int			mTilesPerGenome;		// default: 10,655,006 tiles		per genome
	int			mTilepathPerGenome;		// default:        862 tilepaths	per genome
	int			mBlocksPerGenome;		// default:    332,968 blocks		per genome
	int			mBitvecSize;			// default:  1,331,872 bytes		per genome
	int			mOveroffSize;			// default:	     6,896 bytes        per genome

	DataPtr		mCanonVec;
	DataPtr		mCacheVec;
	DataPtr		mSpanVec;
	DataPtr		mLoqVec;
	DataPtr		mOverflowOffsetVec;
	DataPtr		mOverflowAllocVec;

	DataPtr		mBlockToTilepath;
	DataPtr		mTileMap;

	std::vector<cgf_t*>  mCgfVec;

	DataPtr		mOutMatchList;			// Concordance outputs
	DataPtr		mOutOverflowList;
	DataPtr		mOutTotalList;

	int			mStartPath, mStartStep, mStartTile;
	int			mEndPath, mEndStep;
	
	DataPtr		mTable;

	Allocator	mm;

	bool		mbVerbose;
	int			mStartView;

	#ifdef USE_CUDA
		CUdevice	mCuDevice;				// CUDA
		CUcontext	mCuContext;
		CUmodule	mCuModule;
		CUfunction	mCuConcordanceKernel;
	#endif


};

