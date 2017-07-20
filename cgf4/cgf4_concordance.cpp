#include "cgf4.hpp"

// For the most part, high quality tile concordance is
// straight forward:
//
//   compare the tiles at a particluar position and increment the
//   count if they match, with the total being the number of high
//   quality tile positions.
//
// It's more complicated than this intuitive approach because the genomes
// under consideration are diploid and have spanning tiles involved.
//
// The basic idea is that only so-called "knots" are considered.
// Only high quality knots are counted towards the total high quality
// tile count if the anchor tile of the knots in both datasets are high
// quality.
// Any low quality tile in a knot should render the whole knot low quality
// and the data vectors shoudl reflect this.
// Only matching knots count towards the matching total.
//
// A knot is the smallest run of diploid tiles where each allele starts on
// a tile and a ends on a tile boundary.
//
// Here are some examples of knots:
//
//    0+1:0+1               or [ [0], [0] ]
//    0+2:0+2               or [ [0,-1], [0,-1] ]
//    0+2:4+1,5+1           or [ [0,-1], [4, 5] ]
//    0+2,f+2:3+1,4+2,7+1   or [ [0,-1,15,-1], [3,4,-1,7] ]
//
//
// Here is an example of an invalid knot:
//
//    0+2:f+1,e+2     or [ [0, -1], [f,e,-1] ]
//
// The knot is a self contained group of bi-allelic tiles that doesn't have
// any unmatched trailing spanning tiles hanging over the end.
//
// When doing pair concordance, only knots that are high quality are considered
// towards the total and only matching knots are considered towards the match count.
// Depending on the context, this could be counter-intuitive as individual high quality
// knot or tile counts could vary between the pair high quality knot counts.
//
// For example, consider the following example:
//
//                        |           |       |       |       
// Datasat0:    [ [ 0, -1, 13, -1, -1, 17, -1, 23, -1,  1 ], 
//                [ 1,  1,  2,  2,  1,  1,  1,  1,  1, 29 ] ]
//
// knots:          0+2:1,1 13+3:2,2,1  17+2:1,1 23+2:1,1 1:29
//
//                        |        |  |           |   |
// Datasat1:    [ [ 0, -1, 11, -1,  0, 19, -1, -1, 31,  1 ], 
//                [ 1,  1,  2,  2,  1,  1,  1,  1,  1, 29 ] ]
//
// knots:          0+2:1,1 11+2:2,2 0:1 19+3:1,1,1 31:1 1:29
//
// Assuming both datasets were of high quality, Dataset0 has 5
// knots whereas Dataset1 has 6 knots.
// The concordance between these two would give 4 knots that are
// of high quailty between the two, of which only 2 matched.
// That is, if the two datasets (suitably encoded) were passed
// to thi concordance function, the result would be:
//
//   2 match
//   4 total
//
// Even though each dataset considered individually would have more than
// 4 valid high quality knots.
//
// As stated previously, knots with any no-call in it renders the whole
// knot as a no-call knot an ddoesn't add to the count of the total or
// matched kntos returned.
//

static void print_bin32(uint32_t u32) {
  uint32_t i;
  for (i=0; i<32; i++) {
    if (u32 & (1<<i)) { printf("1"); }
    else printf(".");
  }

  return;

  for (i=0; i<32; i++) {
    if (u32 & (1<<(31-i))) { printf("1"); }
    else printf(".");
  }

}

int cgf_hiq_concordance(int *r_match, int *r_tot,
                        cgf_t *a, cgf_t *b,
                        int start_tile_path, int start_tile_step,
                        int end_tile_path, int end_tile_step) {
  int i, j, k, p;
  uint64_t ii, jj, end_noninc_a, end_noninc_b;
  int match=0, tot=0;

  int start_pos, end_pos, n_pos,
      n_q, n_q_end;

  unsigned char *loq_a, *loq_b,
                *span_a, *span_b,
                *canon_a, *canon_b,
                *cache_ovf_a, *cache_ovf_b;
  uint16_t *overflow_a, *overflow_b;
  uint64_t *overflow64_a, *overflow64_b;

  int hexit_a[8], hexit_b[8],
      hexit_relative_step_a[8], hexit_relative_step_b[8];

  uint32_t loq_mask_a, loq_mask_b,
           hiq_mask_a, hiq_mask_b,
           span_mask_a, span_mask_b,
           xspan_mask_a, xspan_mask_b,
           cache_mask_a, cache_mask_b,
           lo_cache_a, lo_cache_b,
           canon_mask_a, canon_mask_b,
           anchor_mask_a, anchor_mask_b,
           non_anchor_span_mask_a, non_anchor_span_mask_b,
           cache_ovf_mask_a, cache_ovf_mask_b,
           env_mask_a, env_mask_b;

  uint32_t u32;
  uint16_t prev_tile_step16_a, prev_tile_step16_b;

  int tile_step_block_start;
  int anchor_tile_a, anchor_tile_b;

  int stride;

  loq_a = &(a->Loq[0]);
  loq_b = &(b->Loq[0]);

  span_a = &(a->Span[0]);
  span_b = &(b->Span[0]);

  canon_a = &(a->Canon[0]);
  canon_b = &(b->Canon[0]);

  cache_ovf_a = &(a->CacheOverflow[0]);
  cache_ovf_b = &(b->CacheOverflow[0]);

  overflow_a = &(a->Overflow[0]);
  overflow_b = &(b->Overflow[0]);

  overflow64_a = &(a->Overflow64[0]);
  overflow64_b = &(b->Overflow64[0]);

  stride = a->Stride;

  start_pos = 0;
  if (start_tile_path>0) {
    start_pos = 8 * stride * (a->StrideOffset[start_tile_path-1]);
  }
  start_pos += start_tile_step;

  end_pos = 0;
  if (end_tile_path>0) {
    end_pos = 8 * stride * (a->StrideOffset[end_tile_path-1]);
  }
  end_pos += end_tile_step;

  n_pos = end_pos - start_pos + 1;

  n_q = start_pos / (8*stride);
  n_q_end = end_pos / (8*stride);

  printf("lens (%i,%i) (%i,%i) (%i,%i) (%i,%i)\n",
      (int)a->Loq.size(), (int)b->Loq.size(),
      (int)a->Span.size(), (int)b->Span.size(),
      (int)a->Canon.size(), (int)b->Canon.size(),
      (int)a->CacheOverflow.size(), (int)b->CacheOverflow.size());

  printf("pos [%i,%i (+%i)]\n", start_pos, end_pos, n_pos);
  printf("n_q %i, n_q_end %i\n", n_q, n_q_end);

  for (ii=n_q; ii<=n_q_end; ii++) {

    env_mask_a = 0xffffffff;
    env_mask_b = 0xffffffff;

    if (ii==n_q) {
      env_mask_a &= (0xffffffff << (start_tile_step%32));
      env_mask_b &= (0xffffffff << (start_tile_step%32));
    }

    // TODO: edge case when n_q == n_q_end
    //
    if (ii==n_q_end) {
      env_mask_a &= (0xffffffff >> (31-(end_tile_step%32)));
      env_mask_b &= (0xffffffff >> (31-(end_tile_step%32)));
    }

    //DEBUG
    //
    printf("\n\n---\n");
    printf("ii %llu (%llu of [%llu,%llu])\n",
        (unsigned long long)ii, (unsigned long long)ii,
        (unsigned long long)n_q, (unsigned long long)n_q_end);
    printf("env_mask_a %08x\n", (unsigned int)env_mask_a);
    printf("env_mask_b %08x\n", (unsigned int)env_mask_b);


    // collect the uint32_t bit vectors into a convenient form
    //
    loq_mask_a    = loq_a[4*ii] | (loq_a[4*ii+1]<<8) | (loq_a[4*ii+2]<<16) | (loq_a[4*ii+3]<<24);
    span_mask_a   = span_a[4*ii] | (span_a[4*ii+1]<<8) | (span_a[4*ii+2]<<16) | (span_a[4*ii+3]<<24);
    cache_mask_a  = canon_a[4*ii] | (canon_a[4*ii+1]<<8) | (canon_a[4*ii+2]<<16) | (canon_a[4*ii+3]<<24);
    lo_cache_a    = cache_ovf_a[4*ii] | (cache_ovf_a[4*ii+1]<<8) | (cache_ovf_a[4*ii+2]<<16) | (cache_ovf_a[4*ii+3]<<24);

    xspan_mask_a      = ~span_mask_a;
    hiq_mask_a        = ~loq_mask_a;
    hiq_mask_a        &= env_mask_a;

    // non anchor spanning tiles are indicated with a span bit set and a canon bit set
    // so make sure to account for them to get the actual canononical bits out.
    //
    canon_mask_a      = cache_mask_a & xspan_mask_a & hiq_mask_a;

    // anchor tile bit vector for convenience.
    //
    anchor_mask_a     = span_mask_a & hiq_mask_a & (~cache_mask_a);

    // convenience bit vector of anchor tile non-canonical tiles,
    // excluding spanning tiles.
    //
    cache_ovf_mask_a  = (anchor_mask_a & hiq_mask_a) | ((~span_mask_a) & (~canon_mask_a) & hiq_mask_a);

    non_anchor_span_mask_a = span_mask_a & (~anchor_mask_a);


    // --
    //
    loq_mask_b    = loq_b[4*ii] | (loq_b[4*ii+1]<<8) | (loq_b[4*ii+2]<<16) | (loq_b[4*ii+3]<<24);
    span_mask_b   = span_b[4*ii] | (span_b[4*ii+1]<<8) | (span_b[4*ii+2]<<16) | (span_b[4*ii+3]<<24);
    cache_mask_b  = canon_b[4*ii] | (canon_b[4*ii+1]<<8) | (canon_b[4*ii+2]<<16) | (canon_b[4*ii+3]<<24);
    lo_cache_b    = cache_ovf_b[4*ii] | (cache_ovf_b[4*ii+1]<<8) | (cache_ovf_b[4*ii+2]<<16) | (cache_ovf_b[4*ii+3]<<24);

    xspan_mask_b      = ~span_mask_b;
    hiq_mask_b        = ~loq_mask_b;
    hiq_mask_b        &= env_mask_b;

    canon_mask_b      = cache_mask_b & xspan_mask_b & hiq_mask_b;
    anchor_mask_b     = span_mask_b & hiq_mask_b & (~cache_mask_b);
    cache_ovf_mask_b  = (anchor_mask_b & hiq_mask_b) | ((~span_mask_b) & (~canon_mask_b) & hiq_mask_b);
    non_anchor_span_mask_b = span_mask_b & (~anchor_mask_b);


    // The total number of high quality matches is the number of high quality
    // tiles both have in common
    //
    //tot += NumberOfSetBits32( hiq_mask_a & hiq_mask_b );
    tot += NumberOfSetBits32( hiq_mask_a & hiq_mask_b & (~non_anchor_span_mask_a) & (~non_anchor_span_mask_b) );

    // Count the number of canonical matches.
    //
    match += NumberOfSetBits32( canon_mask_a & canon_mask_b );

    //DEBUG
    //
    printf("  hiq   %08x %08x\n",
        (unsigned int)hiq_mask_a,
        (unsigned int)hiq_mask_b);
    printf("  canon %08x %08x\n",
        (unsigned int)canon_mask_a,
        (unsigned int)canon_mask_b);
    printf("  anch  %08x %08x\n",
        (unsigned int)anchor_mask_a,
        (unsigned int)anchor_mask_b);

    printf("  hiq:    "); print_bin32(hiq_mask_a); printf(" "); print_bin32(hiq_mask_b); printf("\n");
    printf("  canon:  "); print_bin32(canon_mask_a); printf(" "); print_bin32(canon_mask_b); printf("\n");
    printf("  anchor: "); print_bin32(anchor_mask_a); printf(" "); print_bin32(anchor_mask_b); printf("\n");
    printf("  span:   "); print_bin32(span_mask_a); printf(" "); print_bin32(span_mask_b); printf("\n");
    printf("  xspan:  "); print_bin32(xspan_mask_a); printf(" "); print_bin32(xspan_mask_b); printf("\n");

    printf("  +match %i, +tot %i\n", match, tot);

    // record the hexit values and tile step positions
    // for both cgfs
    //
    for (i=0; i<8; i++) {
      u32 = ((lo_cache_a & ((uint32_t)0xf<<(4*i))) >> (4*i));
      hexit_a[i] = (int)u32;
      hexit_relative_step_a[i] = -1;
    }

    p = 0;
    for (i=0; i<32; i++) {
      if (cache_ovf_mask_a & (1<<i)) {
        hexit_relative_step_a[p++] = i;
        if (p>=8) { break; }
      }
    }

    for (i=0; i<8; i++) {
      u32 = ((lo_cache_b & ((uint32_t)0xf<<(4*i))) >> (4*i));
      hexit_b[i] = (int)u32;
      hexit_relative_step_b[i] = -1;
    }

    p = 0;
    for (i=0; i<32; i++) {
      if (cache_ovf_mask_b & (1<<i)) {
        hexit_relative_step_b[p++] = i;
        if (p>=8) { break; }
      }
    }

    //DEBUG
    //
    printf("  hexit:\n");
    for (i=0; i<8; i++) {
      printf("  [%i] %i, [%i] %i\n",
          hexit_relative_step_a[i], hexit_a[i],
          hexit_relative_step_b[i], hexit_b[i]);
    }


    // Do a zipper match to count the number of cache overflow hits.
    // Tile variants that overflow from the cache will be picked
    // up by the overflow count.
    //
    for (i=0, j=0; (i<8) && (j<8); ) {
      if ((hexit_relative_step_a[i] < 0) || (hexit_relative_step_b[j] < 0))  { break; }
      if (hexit_relative_step_a[i] < hexit_relative_step_b[j]) { i++; continue; }
      if (hexit_relative_step_a[i] > hexit_relative_step_b[j]) { j++; continue; }
      if (hexit_relative_step_a[i] == hexit_relative_step_b[j]) {
        if ((hexit_a[i] > 0) && (hexit_a[i] < 0xf) &&
            (hexit_b[j] > 0) && (hexit_b[j] < 0xf) &&
            (hexit_a[i]==hexit_b[j])) {
          match++;
        }
        i++;
        j++;
      }
    }

  }

  //DEBUG
  printf("  non-overflow match: %i (%i)\n", match, tot);

  // Count overflow matches
  //
  ii=jj=0;
  if (start_tile_path>0) {
    ii=a->OverflowOffset[start_tile_path-1];
    jj=b->OverflowOffset[start_tile_path-1];
  }

  end_noninc_a = a->OverflowOffset[end_tile_path];
  end_noninc_b = a->OverflowOffset[end_tile_path];

  printf("  [%llu,%llu : %llu,%llu]\n",
      (unsigned long long)ii,
      (unsigned long long)jj,
      (unsigned long long)end_noninc_a,
      (unsigned long long)end_noninc_b);

  anchor_tile_a = 1;
  anchor_tile_b = 1;

  prev_tile_step16_a = 0;
  prev_tile_step16_b = 0;

  while ((ii<end_noninc_a) && (jj<end_noninc_b)) {

    if ((a->Overflow[ii] - prev_tile_step16_a) > 1) {
      anchor_tile_a = 1;
    }

    if ((b->Overflow[ii] - prev_tile_step16_b) > 1) {
      anchor_tile_b = 1;
    }

    if (a->Overflow[ii] < b->Overflow[jj]) {

      if ((a->Overflow[ii+1] == OVF16_MAX) ||
          (a->Overflow[ii+2] == OVF16_MAX)) {
        anchor_tile_a = 0;
      }
      prev_tile_step16_a = a->Overflow[ii];

      ii+=3;
      continue;
    }

    if (a->Overflow[ii] > b->Overflow[jj]) {

      if ((b->Overflow[jj+1] == OVF16_MAX) ||
          (b->Overflow[jj+2] == OVF16_MAX)) {
        anchor_tile_b = 0;
      }
      prev_tile_step16_b = b->Overflow[jj];

      jj+=3;
      continue;
    }

    if (a->Overflow[ii] == b->Overflow[jj]) {
      if ((a->Overflow[ii+1] == b->Overflow[jj+1]) &&
          (a->Overflow[ii+2] == b->Overflow[jj+2])) {

        if ((anchor_tile_a==1) && (anchor_tile_b==1)) {
          match++;
        }

      }

      if (a->Overflow[ii+1] == OVF16_MAX) { anchor_tile_a = 0; }
      if (b->Overflow[ii+1] == OVF16_MAX) { anchor_tile_b = 0; }

      prev_tile_step16_a = a->Overflow[ii];
      prev_tile_step16_b = b->Overflow[jj];

      ii+=3;
      jj+=3;
      continue;
    }

    // should not get here...
    //
    break;
  }



  //DEBUG
  //
  printf("  fin: match %i, tot %i\n", match, tot);

  *r_match = match;
  *r_tot =tot; 

  return 0;
}
