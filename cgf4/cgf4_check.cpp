#include "cgf4.hpp"

int cgf_sanity(cgf_t *cgf) {
  int loc_debug = 0;

  int i, j, k, p, idx, z;
  uint64_t n64, ii, jj, tilepath_idx, idx64;
  uint64_t s64, e64;

  uint32_t u32;

  unsigned char *loq,
                *span,
                *canon,
                *cache_ovf;
  uint16_t *overflow;
  uint64_t *overflow64;

  uint32_t loq_mask=0,
           hiq_mask=0,
           span_mask=0,
           xspan_mask=0,
           cache_mask=0,
           lo_cache=0,
           canon_mask=0,
           anchor_mask=0,
           nonnchor_span_mask=0,
           cache_ovf_mask=0;
           //env_mask=0;

  int tile_step_block_start;
  int anchor_tile;
  int tile_offset;

  int start_pos, end_pos, n_pos,
      n_q, n_q_end;
  int match=0, tot=0, stride=0;

  int hexit[8],
      hexit_relative_step[8];


  int start_block_tile = 0;

  uint64_t n_ovf;
  uint64_t iistart;
  uint64_t jjstart;

  uint16_t u16;

  int n_tile, n_r;

  int err_code = -1;
  tilemap_t tilemap;

  std::vector< uint16_t > cache16_ovf, spillover16, spillover16_knot;
  //std::vector< int > cache_ovf, spillover, spillover_knot;
  std::vector< int > spillover;

  str2tilemap(cgf->TileMap, &tilemap);

  loq = &(cgf->Loq[0]);
  span = &(cgf->Span[0]);
  canon = &(cgf->Canon[0]);
  cache_ovf = &(cgf->CacheOverflow[0]);
  overflow = &(cgf->Overflow[0]);
  overflow64 = &(cgf->Overflow64[0]);


  for (i=0; i<8; i++) {
    if (cgf->Magic[i] != CGF_MAGIC[i]) { return err_code; }
  }
  err_code--;

  if (cgf->CGFVersion.size()==0) { return err_code; }
  err_code--;

  if (cgf->LibraryVersion.size()==0) { return err_code; }
  err_code--;

  n64 = cgf->TilePathCount;
  stride = (int)cgf->Stride;

  if (stride<=0) { return err_code; }
  err_code--;

  if (n64==0) { return 0; }

  if (cgf->Loq.size() != cgf->Span.size()) { return err_code; }
  err_code--;
  if (cgf->Loq.size() != cgf->Canon.size()) { return err_code; }
  err_code--;
  if (cgf->Loq.size() != cgf->CacheOverflow.size()) { return err_code; }
  err_code--;

  stride = cgf->Stride;

  spillover16_knot.clear();

  for (tilepath_idx=0; tilepath_idx<n64; tilepath_idx++) {
    if ((tilepath_idx>0) && (cgf->StrideOffset[tilepath_idx-1] >= cgf->StrideOffset[tilepath_idx])) { return err_code-1; }
    if (cgf->StrideOffset[tilepath_idx] >= cgf->Loq.size()) { return err_code-2; }
    if (cgf->StrideOffset[tilepath_idx] >= cgf->Span.size()) { return err_code-3; }
    if (cgf->StrideOffset[tilepath_idx] >= cgf->Canon.size()) { return err_code-4; }
    if (cgf->StrideOffset[tilepath_idx] >= (stride*cgf->CacheOverflow.size())) { return err_code-5; }

    if (cgf->TileStepCount[tilepath_idx]==0) { continue; }

    start_pos = ( (tilepath_idx>0) ? 8*stride*(cgf->StrideOffset[tilepath_idx-1]) : 0 );
    end_pos = (8*stride*cgf->StrideOffset[tilepath_idx])-1;

    n_q = start_pos / (8*stride);
    n_q_end = end_pos / (8*stride);


    //DEBUG
    if (loc_debug) {
      printf("### TILEPATH %i (%x) n_q %i, n_q_end %i\n", (int)tilepath_idx, (int)tilepath_idx,
          (int)n_q, (int)n_q_end);
    }

    n_tile = cgf->TileStepCount[tilepath_idx];
    e64 = stride*cgf->StrideOffset[tilepath_idx];

    if ((n_tile%32)!=0) {
      u32 = 0;
      u32 |= ((uint32_t)cgf->Loq[ e64 - 4 ]) << 0;
      u32 |= ((uint32_t)cgf->Loq[ e64 - 3 ]) << 8;
      u32 |= ((uint32_t)cgf->Loq[ e64 - 2 ]) << 16;
      u32 |= ((uint32_t)cgf->Loq[ e64 - 1 ]) << 24;

      n_r = 32 - (n_tile%32);
      for (i=(n_tile%32); i<32; i++) {

        if (loc_debug) {
          printf("  tp%i, n_tile%%32 %i, i %i, u32 %x, val %x\n",
              (int)tilepath_idx,
              n_tile%32,
              i,
              (unsigned int)u32,
              (unsigned int) (u32 & ((uint32_t)1<<i)) );
        }

        if ( (u32 & ((uint32_t)1<<i)) == 0) { return err_code - 6; }

      }
    }

    start_block_tile = 0;
    spillover16_knot.clear();
    spillover.clear();
    for (ii=n_q; ii<=n_q_end; ii++, start_block_tile+=32) {


      loq_mask    = loq[4*ii] | (loq[4*ii+1]<<8) | (loq[4*ii+2]<<16) | (loq[4*ii+3]<<24);
      span_mask   = span[4*ii] | (span[4*ii+1]<<8) | (span[4*ii+2]<<16) | (span[4*ii+3]<<24);
      cache_mask  = canon[4*ii] | (canon[4*ii+1]<<8) | (canon[4*ii+2]<<16) | (canon[4*ii+3]<<24);
      lo_cache    = cache_ovf[4*ii] | (cache_ovf[4*ii+1]<<8) | (cache_ovf[4*ii+2]<<16) | (cache_ovf[4*ii+3]<<24);

      xspan_mask      = ~span_mask;
      hiq_mask        = ~loq_mask;
      //hiq_mask        &= env_mask;

      canon_mask      = cache_mask & xspan_mask & hiq_mask;
      anchor_mask     = span_mask & hiq_mask & (~cache_mask);
      cache_ovf_mask  = (anchor_mask & hiq_mask) | ((~span_mask) & (~canon_mask) & hiq_mask);
      nonnchor_span_mask = span_mask & (~anchor_mask);

      tot += NumberOfSetBits32( hiq_mask & (~nonnchor_span_mask) );
      match += NumberOfSetBits32( canon_mask );

      // record hexit values
      // and initialize relative step
      //
      for (i=0; i<8; i++) {
        u32 = ((lo_cache & ((uint32_t)0xf<<(4*i))) >> (4*i));
        hexit[i] = (int)u32;
        hexit_relative_step[i] = -1;
      }

      p = 0;
      for (i=0; i<32; i++) {
        //if (cache_ovf_mask & ((uint32_t)1<<i)) {
        u32 = (uint32_t)1<<i;
        if (cache_ovf_mask & u32) {
          hexit_relative_step[p++] = i;
          if (p>=8) { break; }
        }
      }

      //DEBUG
      //
      if (loc_debug) {
        printf("cache_ovf_mask %08x\n", (unsigned int)cache_ovf_mask);
        printf("cache_ovf:");
        for (i=0; i<8; i++) {
          printf("  [%i](%i %i)  (step:%i)", i, hexit_relative_step[i], hexit[i],
              (int)(start_block_tile+hexit_relative_step[i]) );
        }
        printf("\n");
      }

      spillover.clear();
      for (i=0; i<8; i++) {
        if (hexit_relative_step[i] < 0) {

          if (loc_debug) { printf("  # ending at spillover %i\n", i); }

          break;
        }

        if (loc_debug) {
          printf("  # spillover + (%i, %i)\n",
              (int)hexit_relative_step[i] + start_block_tile,
              hexit[i]);
        }

        if ((hexit[i]==0) || (hexit[i]>=0xf)) { continue; }

        spillover.push_back(hexit_relative_step[i] + start_block_tile);
        spillover.push_back(hexit[i]);
      }

      //DEBUG
      if (loc_debug) {
        printf("  # spillover(%i):", (int)spillover.size());
        for (i=0; i<spillover.size(); i+=2) {
          printf(" (%i,%i)", spillover[i], spillover[i+1]);
        }
        printf("\n");

        printf("spillover size %i (/2 -> %i)\n",
            (int)spillover.size(),
            (int)(spillover.size()/2));
      }

      // expand knots into spillover16_knot
      //
      for (i=0; i<spillover.size(); i+=2) {

        if (loc_debug) {
          printf("  ## spillover[%i] %i %i\n", i, spillover[i], spillover[i+1]);
          printf("  ## tilemap.offset[%i] %i - [%i] %i\n",
              spillover[i+1]-1,
              tilemap.offset[spillover[i+1]-1],
              spillover[i+1],
              tilemap.offset[spillover[i+1]]);
        }

        tile_offset=0;
        for ( j=tilemap.offset[ spillover[i+1]-1 ]; j<tilemap.offset[spillover[i+1]]; j++) {
          spillover16_knot.push_back(spillover[i] + tile_offset);
          u16 = ((tilemap.variant[0][j] < 0) ? OVF16_MAX : (uint16_t)tilemap.variant[0][j]);
          spillover16_knot.push_back( tilemap.variant[0][j] );
          u16 = ((tilemap.variant[1][j] < 0) ? OVF16_MAX : (uint16_t)tilemap.variant[1][j]);
          spillover16_knot.push_back( tilemap.variant[1][j] );

          tile_offset++;
        }
      }

    }

    if (loc_debug) {
      printf("spillover(%i):\n", (int)spillover16_knot.size());
      for (i=0; i<spillover16_knot.size(); i+=3) {
        printf("  spill[%i] %i %i\n",
            (int)spillover16_knot[i],
            (int)spillover16_knot[i+1],
            (int)spillover16_knot[i+2]);
      }
    }

    s64 = 0;
    n_ovf = cgf->OverflowOffset[tilepath_idx];
    if (tilepath_idx>0) {
      s64 = cgf->OverflowOffset[tilepath_idx-1];
      n_ovf -= cgf->OverflowOffset[tilepath_idx-1];
    }


    //printf("Overflow(%i):\n", (int)cgf->TileStepCount[tilepath_idx]);
    if (loc_debug) {
      printf("Overflow([%i+%i]):\n", (int)s64, (int)n_ovf);
      for (i=0; i<n_ovf; i+=3) {
        printf("  ovf[%i] %i %i\n",
            (int)cgf->Overflow[ s64 + i ],
            (int)cgf->Overflow[ s64 + i + 1 ],
            (int)cgf->Overflow[ s64 + i + 2 ]);
      }
    }

    if (loc_debug) {
      printf("... %i\n", (int)tilepath_idx);
      printf("SPILLOVER16(%i)\n", (int)spillover16_knot.size());
      for (ii=0; ii<spillover16_knot.size(); ii+=3) {
        printf("  spillover[%i] %i %i %i\n",
            (int)ii,
            (int)spillover16_knot[ii],
            (int)spillover16_knot[ii+1],
            (int)spillover16_knot[ii+2]);
      }
      printf("\n");
      printf("\n");

      printf("OVERFLOW[%i+%i]:\n",
          (int)s64, (int)n_ovf);
      for (ii=s64; ii<(s64+n_ovf); ii+=3) {
        printf("  Overflow[%i] %i %i %i\n",
            (int)ii,
            (int)cgf->Overflow[ii],
            (int)cgf->Overflow[ii+1],
            (int)cgf->Overflow[ii+2]);
      }

    }


    match =
      overflow_concordance16(NULL, NULL,
                             spillover16_knot, 0, (int)spillover16_knot.size(),
                             cgf->Overflow, (int)s64, (int)(s64+n_ovf), NULL);

    if (match != 0) {

      if (loc_debug) { printf("got match %i\n", match); }

      return err_code-7;
    }

  }
  err_code-=5;






  return 0;
}
