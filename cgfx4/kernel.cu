
#include <stdio.h>

#include "cuda_math.cuh"

__constant__ int tab32[32] = {
	0,   9, 1,  10, 13, 21,  2, 29,
	11, 14, 16, 18, 22, 25,  3, 30,
	8,  12, 20, 28, 15, 17, 24,  7,
	19, 27, 23,  6, 26,  5,  4, 31 };

__device__ inline int log2_32(uint32_t value) {
	value |= value >> 1;
	value |= value >> 2;
	value |= value >> 4;
	value |= value >> 8;
	value |= value >> 16;
	return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

__device__ inline int NumberOfSetBits32(uint32_t u)
{
  u = u - ((u >> 1) & 0x55555555);
  u = (u & 0x33333333) + ((u >> 2) & 0x33333333);
  return (((u + (u >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}


extern "C" __global__ void concordanceKernel ( 
			int cgfa, 
			int blk_start, int blk_end,
			int mStartPath, int mStartStep, int mEndPath, int mEndStep, int mStartTile,
			int mTotalBlocks, int mTotalGenomes,
			uint32_t* mCanonVec, uint32_t* mCacheVec, uint32_t* mSpanVec, uint32_t* mLoqVec,
			uchar* mMatchList, uchar* mTotalList, uint8_t* mTable )
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x + blk_start;		// block id
	int cgfb = blockIdx.y * blockDim.y + threadIdx.y;				// second genome
	if ( ii < blk_start || ii > blk_end ) return;
	if ( cgfb > mTotalGenomes ) return;

	// input reads
	uint32_t hiq_a, hiq_b;
	uint32_t span_a, span_b;	
	uint32_t canon_a, canon_b;	
	uint32_t cache_a, cache_b;
	
	// computed
	uint32_t actual_canon_a, actual_canon_b;
	uint32_t anchor_mask_a, anchor_mask_b;
	uint32_t non_anchor_span_mask_a, non_anchor_span_mask_b;
	uint32_t tmp;

	int hexit_a[8], hexit_b[8],
		hexit_relative_step_a[8], hexit_relative_step_b[8];

	int i, j, p, k;

	// Input bit vector s
	hiq_a = ~ *( (mLoqVec   + cgfa*mTotalBlocks) + ii);		// flip loq_a
	hiq_b = ~ *( (mLoqVec   + cgfb*mTotalBlocks) + ii);		// flip loq_b
	span_a =  *( (mSpanVec  + cgfa*mTotalBlocks) + ii);		// span bits
	span_b =  *( (mSpanVec  + cgfb*mTotalBlocks) + ii);	
	canon_a = *( (mCanonVec + cgfa*mTotalBlocks) + ii);		// canon bits
	canon_b = *( (mCanonVec + cgfb*mTotalBlocks) + ii);
	cache_a = *( (mCacheVec + cgfa*mTotalBlocks) + ii);		// cache bits
	cache_b = *( (mCacheVec + cgfb*mTotalBlocks) + ii);
	
	// Computed bit-vectors
	actual_canon_a = canon_a & ~span_a & hiq_a;		// get actual canon bits, masking off the span bit rules
	anchor_mask_a = span_a & hiq_a & (~canon_a);	// anchor tile bit vector for convenience.
	non_anchor_span_mask_a = span_a & (~anchor_mask_a);

	actual_canon_b = canon_b & ~span_b & hiq_b;
	anchor_mask_b = span_b & hiq_b & (~canon_b);
	non_anchor_span_mask_b = span_b & (~anchor_mask_b);

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
	tmp = (anchor_mask_a & hiq_a) | ((~span_a) & (~actual_canon_a) & hiq_a);	// overflow_a: non-canonical tiles and anchor tiles, excluding spanning tiles.
	for (k = 0; (k<8) && (tmp != 0); k++, tmp &= (tmp - 1)) {
		hexit_relative_step_a[k] = log2_32(tmp & ~(tmp - 1));
	}

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
	tmp = (anchor_mask_b & hiq_b) | ((~span_b) & (~actual_canon_b) & hiq_b);	// overflow_b
	for (k = 0; (k<8) && (tmp != 0); k++, tmp &= (tmp - 1)) {
		hexit_relative_step_b[k] = log2_32(tmp & ~(tmp - 1));
	}
	// Count the number of canonical matches.
	tmp = NumberOfSetBits32(actual_canon_a & actual_canon_b);

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
					tmp++;
				}

				else if ((ii == blk_start) && (ii < blk_end) &&
					((mStartTile + hexit_relative_step_a[i]) >= mStartStep) &&
					((mStartTile + hexit_relative_step_b[j]) >= mStartStep)) {
					tmp++;
				}

				else if ((ii > blk_start) && (ii == blk_end) &&
					((mStartTile + hexit_relative_step_a[i]) <= mEndStep) &&
					((mStartTile + hexit_relative_step_b[j]) <= mEndStep)) {
					tmp++;
				}

				else if ((ii == blk_start) && (ii == blk_end) &&
					((mStartTile + hexit_relative_step_a[i]) >= mStartStep) &&
					((mStartTile + hexit_relative_step_b[j]) >= mStartStep) &&
					((mStartTile + hexit_relative_step_a[i]) <= mEndStep) &&
					((mStartTile + hexit_relative_step_b[j]) <= mEndStep)) {
					tmp++;
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
	* (mMatchList + (cgfb*mTotalBlocks) + ii) = tmp;
	* (mTotalList + (cgfb*mTotalBlocks) + ii) = NumberOfSetBits32(hiq_a & hiq_b & (~non_anchor_span_mask_a) & (~non_anchor_span_mask_b));

}


