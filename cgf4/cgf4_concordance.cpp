#include "cgf4.hpp"

// For the most part, high quality tile concordance is
// conceptually straight forward:
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
// and the data vectors should reflect this.
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
// The knot is a self contained group of bi-allelic tiles that don't have
// any unmatched trailing spanning tiles hanging over the end.
//
// When doing pair concordance, only knots that are high quality are considered
// towards the total and only matching knots are considered towards the match count.
// Depending on the context, this could be counter-intuitive as individual high quality
// knot or tile counts could vary between the pair high quality knot counts.
//
// For example, consider the following example (where '|' denote knot boundaries):
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
// to this concordance method, the result would be:
//
//   2 match
//   4 total
//
// Even though each dataset considered individually would have more than
// 4 valid high quality knots.
//
// As stated previously, knots with any no-call in it renders the whole
// knot as a no-call knot and don't add to the count of the total or
// matched knots returned.
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

// Attemp at an optimized overflow concordance.
// Runs concordance against overflow uint16_t arrays.
// "Elements" in the overflow arrays are grouped into 3
// uint16_t, with the first being the position, the second
// being allele 0 and the third being allele 1.
// `start_[ab]` starts at the real position and the
// end value ends at the real position.
//
// Overflow values are encoded as OVF16_MAX
//
// For example:
//
//   [ 20 0 0 21 OVF16_MAX 0 ... ]
//
//   start_a could be 0, end_noninc_a could be 2.
//
// r_match and r_tot are updated with the number matched
// and the total considered.
// The total is the number of 'knot' positions each overflow array
// shares.
// Though this function is relatively agnostic about the source of the
// data, since the overflow16 vectors typically only store high quality
// tiles, the normal usage is that this will compare only high quality
// overflow elements.
//
inline int overflow_concordance16_z(int *r_match, int *r_tot,
                             uint16_t *a_overflow, int start_a, int end_noninc_a,
                             uint16_t *b_overflow, int start_b, int end_noninc_b,
                             cgf_opt_t *cgf_opt) {

  int anchor_step_a, anchor_step_b;
  int ii, jj, idx;
  int knot_a_start, knot_a_n,
      knot_b_start, knot_b_n;
  uint16_t prev_tile_step16_a, prev_tile_step16_b;
  int match=0, tot=0;

  int loc_debug = 0;

  if (cgf_opt && (cgf_opt->verbose)) { loc_debug=1; }

  if ( ((end_noninc_b-start_b)==0) ||
       ((end_noninc_a-start_a)==0)) {
    return 0;
  }

  ii = start_a;
  jj = start_b;

  anchor_step_a = (int)a_overflow[ii];
  anchor_step_b = (int)b_overflow[jj];

#ifdef CONC_DEBUG
  if (loc_debug) {
    printf("# ovf_conc16 %i %i\n", anchor_step_a, anchor_step_b);
  }
#endif

  knot_a_start = ii;
  knot_a_n=0;
  do {
    prev_tile_step16_a = a_overflow[ii];
    knot_a_n += 3;

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("  a %i %i %i\n", a_overflow[ii], a_overflow[ii+1], a_overflow[ii+2]);
    }
#endif

    ii+=3;
  } while ((ii<end_noninc_a) &&
           ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
           ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );

  knot_b_start = jj;
  knot_b_n = 0;
  do {
    prev_tile_step16_b = b_overflow[jj];
    knot_b_n += 3;

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("  b %i %i %i\n", b_overflow[jj], b_overflow[jj+1], b_overflow[jj+2]);
    }

    if (loc_debug) {
      printf("  knot_b++ %i %i %i\n",
          (int)prev_tile_step16_b,
          (int)b_overflow[jj+1],
          (int)b_overflow[jj+2]);
    }
#endif

    jj+=3;
  } while ((jj<end_noninc_b) &&
           ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
           ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );

  while ((ii<end_noninc_a) && (jj<end_noninc_b)) {

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("ii %i (/%i), jj %i (/%i)\n",
          (int)ii, (int)end_noninc_a,
          (int)jj, (int)end_noninc_b);
      printf("  anchor a:%i (s%i), b:%i (s%i)\n",
          anchor_step_a, knot_a_n,
          anchor_step_b, knot_b_n);

      printf("  knot_a (step %i):", anchor_step_a);
      for (idx=0; idx<knot_a_n; idx++) { printf(" %i", a_overflow[knot_a_start+idx]); }
      printf("\n");
      printf("  knot_b (step %i):", anchor_step_b);
      for (idx=0; idx<knot_b_n; idx++) { printf(" %i", b_overflow[knot_b_start+idx]); }
      printf("\n");
      fflush(stdout);
    }
#endif

    if (anchor_step_a == anchor_step_b) {

      tot++;

#ifdef CONC_DEBUG
      if (loc_debug) { printf("  ==\n"); }
#endif

      if ( (knot_a_n>0) &&
           (knot_a_n == knot_b_n) &&
           (a_overflow[knot_a_start+1] != OVF16_MAX) &&
           (a_overflow[knot_a_start+2] != OVF16_MAX) &&
           (b_overflow[knot_b_start+1] != OVF16_MAX) &&
           (b_overflow[knot_b_start+2] != OVF16_MAX) ) {
        for (idx = 0; idx<knot_a_n; idx++) {
          if (a_overflow[knot_a_start + idx] != b_overflow[knot_b_start + idx]) { break; }
        }
        if (idx == knot_a_n) {

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("MATCH %i+%i (match count %i)\n", anchor_step_a, (int)knot_a_n/3, match);
            printf("  match! ovf16_z (a) %i+%i (knot_a_n %i)\n", anchor_step_a, (int)knot_a_n/3, (int)knot_a_n);
            fflush(stdout);
          }
#endif

          match++;
        }
      }

      knot_a_start = ii;
      knot_a_n=0;

      knot_b_start = jj;
      knot_b_n=0;

      if (ii >= end_noninc_a) { continue; }
      if (jj >= end_noninc_b) { continue; }

      anchor_step_a = (int)a_overflow[ii];
      do {
        prev_tile_step16_a = a_overflow[ii];
        knot_a_n+=3;

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_a++ %i %i %i\n",
              (int)a_overflow[ii],
              (int)a_overflow[ii+1],
              (int)a_overflow[ii+2]);
        }
#endif

        ii+=3;
      } while ((ii<end_noninc_a) &&
               ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
               ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );

      anchor_step_b = (int)b_overflow[jj];

      do {
        prev_tile_step16_b = b_overflow[jj];
        knot_b_n+=3;

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_b++ %i %i %i\n",
              (int)b_overflow[jj],
              (int)b_overflow[jj+1],
              (int)b_overflow[jj+2]);
        }
#endif


        jj+=3;
      } while ((jj<end_noninc_b) &&
               ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
               ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );

      continue;
    }

    if (anchor_step_a < anchor_step_b) {

#ifdef CONC_DEBUG
      if (loc_debug) { printf("  <\n");  }
#endif

      knot_a_start = ii;
      knot_a_n=0;

      if (ii >= end_noninc_a) { continue; }

      anchor_step_a = (int)a_overflow[ii];
      do {
        prev_tile_step16_a = a_overflow[ii];
        knot_a_n+=3;

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_a++ %i %i %i\n",
              (int)a_overflow[ii],
              (int)a_overflow[ii+1],
              (int)a_overflow[ii+2]);
        }
#endif

        ii+=3;
      } while ((ii<end_noninc_a) &&
               ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
               ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );

      continue;
    }

    if (anchor_step_a > anchor_step_b) {

#ifdef CONC_DEBUG
      if (loc_debug) { printf("  >\n"); }
#endif

      knot_b_start = jj;
      knot_b_n=0;

      if (jj >= end_noninc_b) { continue; }

      anchor_step_b = (int)b_overflow[jj];
      do {
        prev_tile_step16_b = b_overflow[jj];
        knot_b_n+=3;

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_b++ %i %i %i\n",
              (int)b_overflow[jj],
              (int)b_overflow[jj+1],
              (int)b_overflow[jj+2]);
        }
#endif


        jj+=3;
      } while ((jj<end_noninc_b) &&
               ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
               ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );

      continue;
    }


    if (anchor_step_a == anchor_step_b) {

      if ( (knot_a_n>0) &&
           (knot_a_n == knot_b_n) &&
           (a_overflow[knot_a_start + 1] != OVF16_MAX) &&
           (a_overflow[knot_a_start + 2] != OVF16_MAX) &&
           (b_overflow[knot_b_start + 1] != OVF16_MAX) &&
           (b_overflow[knot_b_start + 2] != OVF16_MAX) ) {
        for (idx = 0; idx<knot_a_n; idx++) {
          if (a_overflow[knot_a_start + idx] != b_overflow[knot_b_start + idx]) { break; }
        }
        if (idx == knot_a_n) {

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("MATCH %i+%i\n", anchor_step_a, (int)knot_a_n/3);
            printf("  match! ovf16 (b) %i+%i (knot_a.size %i)\n", anchor_step_a, (int)knot_a_n/3, (int)knot_a_n);
          }
#endif

          match++;
        }
      }
    }

  }

  // end conditions
  //

  // Skip ahead to the appropriate place in the 'a'
  // structure, recording the last knot to be considered
  //
  while ((ii<end_noninc_a) &&
         (anchor_step_a < anchor_step_b)) {

#ifdef CONC_DEBUG
    if (loc_debug) { printf("  < (fin)\n");  }
#endif

    knot_a_start = ii;
    knot_a_n = 0;

    if (ii >= end_noninc_a) { continue; }

    anchor_step_a = (int)a_overflow[ii];
    do {
      prev_tile_step16_a = a_overflow[ii];
      knot_a_n+=3;

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("  knot_a++ %i %i %i\n",
            (int)a_overflow[ii],
            (int)a_overflow[ii+1],
            (int)a_overflow[ii+2]);
      }
#endif

      ii+=3;
    } while ((ii<end_noninc_a) &&
             ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
             ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );
  }

  // Skip ahead to the appropriate place in the 'b'
  // structure, recording the last knot to be considered
  //
  while ((jj<end_noninc_b) &&
         (anchor_step_a > anchor_step_b)) {

#ifdef CONC_DEBUG
    if (loc_debug) { printf("  > (fin)\n"); }
#endif

    knot_b_start = jj;
    knot_b_n = 0;

    if (jj >= end_noninc_b) { continue; }

    anchor_step_b = (int)b_overflow[jj];
    do {
      prev_tile_step16_b = b_overflow[jj];
      knot_b_n+=3;

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("  knot_b++ %i %i %i\n",
            (int)b_overflow[jj],
            (int)b_overflow[jj+1],
            (int)b_overflow[jj+2]);
      }
#endif


      jj+=3;
    } while ((jj<end_noninc_b) &&
             ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
             ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );
  }




  if (anchor_step_a == anchor_step_b) {

    if ( (knot_a_n>0) &&
         (knot_a_n == knot_b_n) &&
         (a_overflow[knot_a_start + 1] != OVF16_MAX) &&
         (a_overflow[knot_a_start + 2] != OVF16_MAX) &&
         (b_overflow[knot_b_start + 1] != OVF16_MAX) &&
         (b_overflow[knot_b_start + 2] != OVF16_MAX) ) {
      for (idx = 0; idx<knot_a_n; idx++) {
        if (a_overflow[knot_a_start + idx] != b_overflow[knot_b_start + idx]) { break; }
      }

      if (idx == knot_a_n) {

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("MATCH %i+%i\n", anchor_step_a, (int)knot_a_n/3);
          printf("  match! ovf16 (c) %i+%i (knot_a.size %i)\n", anchor_step_a, (int)knot_a_n/3, (int)knot_a_n);
          printf("  match! ovf16 (c..) %i+%i (knot_b.size %i)\n", anchor_step_b, (int)knot_b_n/3, (int)knot_b_n);
          fflush(stdout);
        }
#endif

        match++;
      }
    }
  }

  if (r_match) { *r_match = match; }
  if (r_tot) { *r_tot = tot; }

  return match;


}

// "unoptimized" overflow concordance function.
// Not used but functional and kept as reference.
//
int overflow_concordance16(int *r_match, int *r_tot,
                           std::vector<uint16_t> &a_overflow, int start_a, int end_noninc_a,
                           std::vector<uint16_t> &b_overflow, int start_b, int end_noninc_b,
                           cgf_opt_t *cgf_opt) {
  int anchor_step_a, anchor_step_b;
  int ii, jj, idx;
  std::vector<uint16_t> knot_a, knot_b;
  uint16_t prev_tile_step16_a, prev_tile_step16_b;
  int match=0, tot=0;

  int loc_debug = 0;

  size_t knot_a_size=0, knot_b_size=0;

  knot_a.reserve(27);
  knot_b.reserve(27);

  if (cgf_opt && (cgf_opt->verbose)) { loc_debug=1; }

  if ( ((end_noninc_b-start_b)==0) ||
       ((end_noninc_a-start_a)==0)) {
    return 0;
  }

  ii = start_a;
  jj = start_b;

  anchor_step_a = (int)a_overflow[ii];
  anchor_step_b = (int)b_overflow[jj];

#ifdef CONC_DEBUG
  if (loc_debug) {
    printf("# ovf_conc16 %i %i\n", anchor_step_a, anchor_step_b);
  }
#endif

  knot_a.clear();
  do {
    prev_tile_step16_a = a_overflow[ii];
    knot_a.push_back(a_overflow[ii+1]);
    knot_a.push_back(a_overflow[ii+2]);

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("  a %i %i %i\n", a_overflow[ii], a_overflow[ii+1], a_overflow[ii+2]);
    }
#endif

    ii+=3;
  } while ((ii<end_noninc_a) &&
           ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
           ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );

  knot_b.clear();
  do {
    prev_tile_step16_b = b_overflow[jj];
    knot_b.push_back(b_overflow[jj+1]);
    knot_b.push_back(b_overflow[jj+2]);

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("  b %i %i %i\n", b_overflow[jj], b_overflow[jj+1], b_overflow[jj+2]);
    }

    if (loc_debug) {
      printf("  knot_b++ %i %i %i\n",
          (int)prev_tile_step16_b,
          (int)b_overflow[jj+1],
          (int)b_overflow[jj+2]);
    }
#endif

    jj+=3;
  } while ((jj<end_noninc_b) &&
           ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
           ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );

  while ((ii<end_noninc_a) && (jj<end_noninc_b)) {

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("ii %i (/%i), jj %i (/%i)\n",
          (int)ii, (int)end_noninc_a,
          (int)jj, (int)end_noninc_b);
      printf("  anchor a:%i (s%i), b:%i (s%i)\n",
          anchor_step_a, (int)knot_a.size(),
          anchor_step_b, (int)knot_b.size());

      printf("  knot_a (step %i):", anchor_step_a);
      for (idx=0; idx<knot_a.size(); idx++) { printf(" %i", knot_a[idx]); }
      printf("\n");
      printf("  knot_b (step %i):", anchor_step_b);
      for (idx=0; idx<knot_b.size(); idx++) { printf(" %i", knot_b[idx]); }
      printf("\n");
      fflush(stdout);
    }
#endif

    if (anchor_step_a == anchor_step_b) {

      tot++;

#ifdef CONC_DEBUG
      if (loc_debug) { printf("  ==\n"); }
#endif

      knot_a_size = knot_a.size();

      if ( (knot_a_size>0) &&
           (knot_a_size == knot_b.size()) &&
           (knot_a[0] != OVF16_MAX) &&
           (knot_a[1] != OVF16_MAX) &&
           (knot_b[0] != OVF16_MAX) &&
           (knot_b[1] != OVF16_MAX) ) {
        for (idx = 0; idx<knot_a_size; idx++) {
          if (knot_a[idx] != knot_b[idx]) { break; }
        }
        if (idx == knot_a_size) {

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("MATCH %i+%i (match count %i)\n", anchor_step_a, (int)knot_a.size()/2, match);
            printf("  match! ovf16 (a) %i+%i (knot_a.size %i)\n", anchor_step_a, (int)knot_a.size()/2, (int)knot_a.size());
            fflush(stdout);
          }
#endif

          match++;
        }
      }

      knot_a.clear();
      knot_b.clear();

      if (ii >= end_noninc_a) { continue; }
      if (jj >= end_noninc_b) { continue; }

      anchor_step_a = (int)a_overflow[ii];

      do {
        prev_tile_step16_a = a_overflow[ii];
        knot_a.push_back((int)a_overflow[ii+1]);
        knot_a.push_back((int)a_overflow[ii+2]);

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_a++ %i %i %i\n",
              (int)a_overflow[ii],
              (int)a_overflow[ii+1],
              (int)a_overflow[ii+2]);
        }
#endif

        ii+=3;
      } while ((ii<end_noninc_a) &&
               ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
               ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );

      anchor_step_b = (int)b_overflow[jj];

      do {
        prev_tile_step16_b = b_overflow[jj];
        knot_b.push_back((int)b_overflow[jj+1]);
        knot_b.push_back((int)b_overflow[jj+2]);

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_b++ %i %i %i\n",
              (int)b_overflow[jj],
              (int)b_overflow[jj+1],
              (int)b_overflow[jj+2]);
        }
#endif


        jj+=3;
      } while ((jj<end_noninc_b) &&
               ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
               ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );

      continue;
    }

    if (anchor_step_a < anchor_step_b) {

#ifdef CONC_DEBUG
      if (loc_debug) { printf("  <\n");  }
#endif

      knot_a.clear();

      if (ii >= end_noninc_a) { continue; }

      anchor_step_a = (int)a_overflow[ii];

      do {
        prev_tile_step16_a = a_overflow[ii];
        knot_a.push_back((int)a_overflow[ii+1]);
        knot_a.push_back((int)a_overflow[ii+2]);

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_a++ %i %i %i\n",
              (int)a_overflow[ii],
              (int)a_overflow[ii+1],
              (int)a_overflow[ii+2]);
        }
#endif

        ii+=3;
      } while ((ii<end_noninc_a) &&
               ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
               ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );

      continue;
    }

    if (anchor_step_a > anchor_step_b) {

#ifdef CONC_DEBUG
      if (loc_debug) { printf("  >\n"); }
#endif

      knot_b.clear();

      if (jj >= end_noninc_b) { continue; }

      anchor_step_b = (int)b_overflow[jj];

      do {
        prev_tile_step16_b = b_overflow[jj];
        knot_b.push_back((int)b_overflow[jj+1]);
        knot_b.push_back((int)b_overflow[jj+2]);

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  knot_b++ %i %i %i\n",
              (int)b_overflow[jj],
              (int)b_overflow[jj+1],
              (int)b_overflow[jj+2]);
        }
#endif


        jj+=3;
      } while ((jj<end_noninc_b) &&
               ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
               ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );

      continue;
    }


    if (anchor_step_a == anchor_step_b) {

      knot_a_size = knot_a.size();

      if ( (knot_a_size>0) &&
           (knot_a_size == knot_b.size()) &&
           (knot_a[0] != OVF16_MAX) &&
           (knot_a[1] != OVF16_MAX) &&
           (knot_b[0] != OVF16_MAX) &&
           (knot_b[1] != OVF16_MAX) ) {
        for (idx = 0; idx<knot_a_size; idx++) {
          if (knot_a[idx] != knot_b[idx]) { break; }
        }
        if (idx == knot_a_size) {

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("MATCH %i+%i\n", anchor_step_a, (int)knot_a.size()/2);
            printf("  match! ovf16 (b) %i+%i (knot_a.size %i)\n", anchor_step_a, (int)knot_a.size()/2, (int)knot_a.size());
          }
#endif

          match++;
        }
      }
    }

  }

  // end conditions
  //

  // Skip ahead to the appropriate place in the 'a'
  // structure, recording the last knot to be considered
  //
  while ((ii<end_noninc_a) &&
         (anchor_step_a < anchor_step_b)) {

#ifdef CONC_DEBUG
    if (loc_debug) { printf("  < (fin)\n");  }
#endif

    knot_a.clear();

    if (ii >= end_noninc_a) { continue; }

    anchor_step_a = (int)a_overflow[ii];

    do {
      prev_tile_step16_a = a_overflow[ii];
      knot_a.push_back((int)a_overflow[ii+1]);
      knot_a.push_back((int)a_overflow[ii+2]);

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("  knot_a++ %i %i %i\n",
            (int)a_overflow[ii],
            (int)a_overflow[ii+1],
            (int)a_overflow[ii+2]);
      }
#endif

      ii+=3;
    } while ((ii<end_noninc_a) &&
             ((a_overflow[ii] - prev_tile_step16_a) <= 1) &&
             ((a_overflow[ii+1] == OVF16_MAX) || (a_overflow[ii+2] == OVF16_MAX)) );
  }

  // Skip ahead to the appropriate place in the 'b'
  // structure, recording the last knot to be considered
  //
  while ((jj<end_noninc_b) &&
         (anchor_step_a > anchor_step_b)) {

#ifdef CONC_DEBUG
    if (loc_debug) { printf("  > (fin)\n"); }
#endif

    knot_b.clear();

    if (jj >= end_noninc_b) { continue; }

    anchor_step_b = (int)b_overflow[jj];

    do {
      prev_tile_step16_b = b_overflow[jj];
      knot_b.push_back((int)b_overflow[jj+1]);
      knot_b.push_back((int)b_overflow[jj+2]);

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("  knot_b++ %i %i %i\n",
            (int)b_overflow[jj],
            (int)b_overflow[jj+1],
            (int)b_overflow[jj+2]);
      }
#endif


      jj+=3;
    } while ((jj<end_noninc_b) &&
             ((b_overflow[jj] - prev_tile_step16_b) <= 1) &&
             ((b_overflow[jj+1] == OVF16_MAX) || (b_overflow[jj+2] == OVF16_MAX)) );
  }




  if (anchor_step_a == anchor_step_b) {

    knot_a_size = knot_a.size();

    if ( (knot_a_size>0) &&
         (knot_a_size == knot_b.size()) &&
         (knot_a[0] != OVF16_MAX) &&
         (knot_a[1] != OVF16_MAX) &&
         (knot_b[0] != OVF16_MAX) &&
         (knot_b[1] != OVF16_MAX) ) {
      for (idx = 0; idx<knot_a_size; idx++) {
        if (knot_a[idx] != knot_b[idx]) { break; }
      }

      if (idx == knot_a_size) {

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("MATCH %i+%i\n", anchor_step_a, (int)knot_a.size()/2);
          printf("  match! ovf16 (c) %i+%i (knot_a.size %i)\n", anchor_step_a, (int)knot_a.size()/2, (int)knot_a.size());
          printf("  match! ovf16 (c..) %i+%i (knot_b.size %i)\n", anchor_step_b, (int)knot_b.size()/2, (int)knot_b.size());
          fflush(stdout);
        }
#endif

        match++;
      }
    }
  }



  if (r_match) { *r_match = match; }
  if (r_tot) { *r_tot = tot; }

  return match;

}



// The high quality concordance does a few different types of counts:
// * Count the number of canonical (high quality) bits set
// * Count the number of cache overflow entities in common
// * Count the number of Overflow entities in common
// * Count the number of Overflow64 entities in common
// * Count the number of unmatched cache overflow entities
//   to both of the Overflow arrays
//
// Finding the total and mathced count of elements kept in the
// 'cache' is relatively straight forward.
// Any tiles in a cache that that overflow are saved and then
// compared against the overflow vectors.
// Finally, the overflow vectors are compared at the end.
//
// All high quality tiles have a corresponding low bit set
// in the `Loq` vector.
// All spanning tiles, including the anchor tile, have their
// corresponding bit set in the `Span` vector.
// The `Canon` bit vector is set high if it's the default
// tile or if it's a non-anchor spanning tile.  That is:
//
// Canon | Span
// ------------
//    0      1     -> anchor spanning tile
//    1      1     -> non-anchor spanning tile
//
// That is, anchor tiles are indicated by a 0 Canon bit and
// 1 spanning bit.
//
// `CacheOverflow` is segmented into blocks that can contain up
// to 8 high quality non-canonical (non-default) tiles.
// The count can be determined by the number of canonical
// high quality bits set.
// The codes for the `CacheOverflow` elements are:
//
//   * 0x0              complex
//   * 0x1 to 0xe       lookup value in tile map (0x1 maps to tile map entry 1 etc.)
//   * 0xf              overflow
//
// Testing `Overflow64` is not implemented.
//
// Different tilepaths are padded in each of the arrays with high quality
// values set to 0 for the tiles that run over.
// This allows for an easier concordance test as tilepath boundaries
// don't need special consideration.
//
int cgf_hiq_concordance(int *r_match, int *r_tot,
                        cgf_t *a, cgf_t *b,
                        int start_tile_path, int start_tile_step,
                        int end_tile_path_inc, int end_tile_step_inc,
                        cgf_opt_t *cgf_opt) {
  int loc_debug = 0;

  int i, j, k, p;
  uint64_t ii, jj, end_noninc_a, end_noninc_b;
  uint64_t ii4, jj4;
  uint16_t u16;

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
           cache_ovf_mask_a, cache_ovf_mask_b;

  uint32_t env_mask;

  uint32_t u32, t_u32;
  uint16_t prev_tile_step16_a, prev_tile_step16_b;
  int prev_tile_step_a, prev_tile_step_b;

  int tile_step_block_start;
  int tile_offset;

  int stride;
  int tilepath_idx;

  int start_block_tile = 0;
  int idx, z;

  int t_match, t_tot;

  uint64_t iistart;
  uint64_t jjstart;

  std::vector<int> knot_a, knot_b;
  int knot_len=0;

  tilemap_t *tilemap, tm;

  // interleaved step, cache val
  //
  std::vector< int > spillover_a, spillover_b;

  // 3 values interleaved like Overflow arrays
  //
  std::vector< int > spillover_knot_a, spillover_knot_b;
  std::vector< uint16_t > spillover16_knot_a, spillover16_knot_b;

  if (cgf_opt && (cgf_opt->verbose)) {
    loc_debug=1;
  }

  knot_a.clear();             knot_b.clear();
  spillover_a.clear();        spillover_b.clear();
  spillover_knot_a.clear();   spillover_knot_b.clear();
  spillover16_knot_a.clear(); spillover16_knot_b.clear();

  spillover_a.reserve(1024);
  spillover_b.reserve(1024);

  spillover16_knot_a.reserve(2048);
  spillover16_knot_b.reserve(2048);

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

  // reuse cached tilemap if available
  //
  if (a->TileMapCacheInit==0) {
    str2tilemap(a->TileMap, &tm);
    tilemap = &tm;
  }
  else {
    tilemap = &(a->TileMapCache);
  }

  stride = a->Stride;
  tilepath_idx = start_tile_path;

  start_pos = 0;
  if (start_tile_path>0) {
    start_pos = 8 * stride * (a->StrideOffset[start_tile_path-1]);
  }
  start_pos += start_tile_step;

  end_pos = 0;
  if (end_tile_path_inc>0) {
    end_pos = 8 * stride * (a->StrideOffset[end_tile_path_inc-1]);
  }
  end_pos += end_tile_step_inc;

  n_pos = end_pos - start_pos + 1;

  n_q = start_pos / (8*stride);
  n_q_end = end_pos / (8*stride);

#ifdef CONC_DEBUG
  if (loc_debug) {
    printf("lens (%i,%i) (%i,%i) (%i,%i) (%i,%i)\n",
        (int)a->Loq.size(), (int)b->Loq.size(),
        (int)a->Span.size(), (int)b->Span.size(),
        (int)a->Canon.size(), (int)b->Canon.size(),
        (int)a->CacheOverflow.size(), (int)b->CacheOverflow.size());

    printf("pos [%i,%i (+%i)]\n", start_pos, end_pos, n_pos);
    printf("n_q %i, n_q_end %i\n", n_q, n_q_end);
  }
#endif

  start_block_tile = start_tile_step / (8*stride);
  start_block_tile *= 32;

  for (ii=n_q; ii<=n_q_end; ii++, start_block_tile+=32) {

    // Reset relevant tile path and tile step information
    //
    env_mask = 0xffffffff;

    if (ii==n_q) {
      env_mask &= (0xffffffff << (start_tile_step%32));
    }

    if (ii==n_q_end) {
      env_mask &= (0xffffffff >> (31-(end_tile_step_inc%32)));
    }

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("\n\n---\n");
      printf("start_block_tile: %i (%x)\n", start_block_tile, start_block_tile);
      printf("ii %llu (%llu of [%llu,%llu])\n",
          (unsigned long long)ii, (unsigned long long)ii,
          (unsigned long long)n_q, (unsigned long long)n_q_end);
      printf("  env_mask: "); print_bin32(env_mask); printf("\n");
    }
#endif


    // collect the uint32_t bit vectors into a convenient form
    //
    ii4 = 4*ii;
    loq_mask_a    = loq_a[ii4]        | (loq_a[ii4+1]<<8)       | (loq_a[ii4+2]<<16)        | (loq_a[ii4+3]<<24);
    span_mask_a   = span_a[ii4]       | (span_a[ii4+1]<<8)      | (span_a[ii4+2]<<16)       | (span_a[ii4+3]<<24);
    cache_mask_a  = canon_a[ii4]      | (canon_a[ii4+1]<<8)     | (canon_a[ii4+2]<<16)      | (canon_a[ii4+3]<<24);
    lo_cache_a    = cache_ovf_a[ii4]  | (cache_ovf_a[ii4+1]<<8) | (cache_ovf_a[ii4+2]<<16)  | (cache_ovf_a[ii4+3]<<24);

    xspan_mask_a      = ~span_mask_a;
    hiq_mask_a        = ~loq_mask_a;

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
    loq_mask_b    = loq_b[ii4]        | (loq_b[ii4+1]<<8)       | (loq_b[ii4+2]<<16)        | (loq_b[ii4+3]<<24);
    span_mask_b   = span_b[ii4]       | (span_b[ii4+1]<<8)      | (span_b[ii4+2]<<16)       | (span_b[ii4+3]<<24);
    cache_mask_b  = canon_b[ii4]      | (canon_b[ii4+1]<<8)     | (canon_b[ii4+2]<<16)      | (canon_b[ii4+3]<<24);
    lo_cache_b    = cache_ovf_b[ii4]  | (cache_ovf_b[ii4+1]<<8) | (cache_ovf_b[ii4+2]<<16)  | (cache_ovf_b[ii4+3]<<24);

    xspan_mask_b      = ~span_mask_b;
    hiq_mask_b        = ~loq_mask_b;

    canon_mask_b      = cache_mask_b & xspan_mask_b & hiq_mask_b;
    anchor_mask_b     = span_mask_b & hiq_mask_b & (~cache_mask_b);
    cache_ovf_mask_b  = (anchor_mask_b & hiq_mask_b) | ((~span_mask_b) & (~canon_mask_b) & hiq_mask_b);
    non_anchor_span_mask_b = span_mask_b & (~anchor_mask_b);


    // The total number of high quality matches is the number of high quality
    // tiles both have in common, excluding the non-anchor high quality tiles
    //
    tot += NumberOfSetBits32( env_mask & hiq_mask_a & hiq_mask_b & (~non_anchor_span_mask_a) & (~non_anchor_span_mask_b) );

    // Count the number of canonical matches.
    //
    match += NumberOfSetBits32( env_mask & canon_mask_a & canon_mask_b );

#ifdef CONC_DEBUG
    if (loc_debug) {
      t_u32 = env_mask & canon_mask_a & canon_mask_b;
      for (i=0; i<32; i++) {
        if (t_u32 & ((uint32_t)1<<i)) {
          printf("MATCH %i+1 (match %i)\n", start_block_tile+i, match);
        }
      }
    }

    if (loc_debug) {
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
    }
#endif

    hexit_a[0] = (lo_cache_a & ((uint32_t)0xf));
    hexit_a[1] = (lo_cache_a & ((uint32_t)0xf<<4)) >> 4;
    hexit_a[2] = (lo_cache_a & ((uint32_t)0xf<<8)) >> 8;
    hexit_a[3] = (lo_cache_a & ((uint32_t)0xf<<12)) >> 12;

    hexit_a[4] = (lo_cache_a & ((uint32_t)0xf<<16)) >> 16;
    hexit_a[5] = (lo_cache_a & ((uint32_t)0xf<<20)) >> 20;
    hexit_a[6] = (lo_cache_a & ((uint32_t)0xf<<24)) >> 24;
    hexit_a[7] = (lo_cache_a & ((uint32_t)0xf<<28)) >> 28;

    hexit_relative_step_a[0] = -1;
    hexit_relative_step_a[1] = -1;
    hexit_relative_step_a[2] = -1;
    hexit_relative_step_a[3] = -1;
    hexit_relative_step_a[4] = -1;
    hexit_relative_step_a[5] = -1;
    hexit_relative_step_a[6] = -1;
    hexit_relative_step_a[7] = -1;

    p=0;
    do {

      if (cache_ovf_mask_a & 0x1u) { hexit_relative_step_a[p++] = 0; }
      if (cache_ovf_mask_a & 0x2u) { hexit_relative_step_a[p++] = 1; }
      if (cache_ovf_mask_a & 0x4u) { hexit_relative_step_a[p++] = 2; }
      if (cache_ovf_mask_a & 0x8u) { hexit_relative_step_a[p++] = 3; }

      if (cache_ovf_mask_a & 0x10u) { hexit_relative_step_a[p++] = 4; }
      if (cache_ovf_mask_a & 0x20u) { hexit_relative_step_a[p++] = 5; }
      if (cache_ovf_mask_a & 0x40u) { hexit_relative_step_a[p++] = 6; }
      if (cache_ovf_mask_a & 0x80u) { hexit_relative_step_a[p++] = 7; if (p>=8) { break; } }

      if (cache_ovf_mask_a & 0x100u) { hexit_relative_step_a[p++] = 8; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x200u) { hexit_relative_step_a[p++] = 9; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x400u) { hexit_relative_step_a[p++] = 10; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x800u) { hexit_relative_step_a[p++] = 11; if (p>=8) { break; } }

      if (cache_ovf_mask_a & 0x1000u) { hexit_relative_step_a[p++] = 12; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x2000u) { hexit_relative_step_a[p++] = 13; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x4000u) { hexit_relative_step_a[p++] = 14; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x8000u) { hexit_relative_step_a[p++] = 15; if (p>=8) { break; } }

      if (cache_ovf_mask_a & 0x10000u) { hexit_relative_step_a[p++] = 16; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x20000u) { hexit_relative_step_a[p++] = 17; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x40000u) { hexit_relative_step_a[p++] = 18; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x80000u) { hexit_relative_step_a[p++] = 19; if (p>=8) { break; } }

      if (cache_ovf_mask_a & 0x100000u) { hexit_relative_step_a[p++] = 20; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x200000u) { hexit_relative_step_a[p++] = 21; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x400000u) { hexit_relative_step_a[p++] = 22; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x800000u) { hexit_relative_step_a[p++] = 23; if (p>=8) { break; } }

      if (cache_ovf_mask_a & 0x1000000u) { hexit_relative_step_a[p++] = 24; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x2000000u) { hexit_relative_step_a[p++] = 25; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x4000000u) { hexit_relative_step_a[p++] = 26; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x8000000u) { hexit_relative_step_a[p++] = 27; if (p>=8) { break; } }

      if (cache_ovf_mask_a & 0x10000000u) { hexit_relative_step_a[p++] = 28; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x20000000u) { hexit_relative_step_a[p++] = 29; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x40000000u) { hexit_relative_step_a[p++] = 30; if (p>=8) { break; } }
      if (cache_ovf_mask_a & 0x80000000u) { hexit_relative_step_a[p++] = 31; if (p>=8) { break; } }

    } while (0);

    hexit_b[0] = (lo_cache_b & ((uint32_t)0xf));
    hexit_b[1] = (lo_cache_b & ((uint32_t)0xf<<4)) >> 4;
    hexit_b[2] = (lo_cache_b & ((uint32_t)0xf<<8)) >> 8;
    hexit_b[3] = (lo_cache_b & ((uint32_t)0xf<<12)) >> 12;

    hexit_b[4] = (lo_cache_b & ((uint32_t)0xf<<16)) >> 16;
    hexit_b[5] = (lo_cache_b & ((uint32_t)0xf<<20)) >> 20;
    hexit_b[6] = (lo_cache_b & ((uint32_t)0xf<<24)) >> 24;
    hexit_b[7] = (lo_cache_b & ((uint32_t)0xf<<28)) >> 28;

    hexit_relative_step_b[0] = -1;
    hexit_relative_step_b[1] = -1;
    hexit_relative_step_b[2] = -1;
    hexit_relative_step_b[3] = -1;
    hexit_relative_step_b[4] = -1;
    hexit_relative_step_b[5] = -1;
    hexit_relative_step_b[6] = -1;
    hexit_relative_step_b[7] = -1;

    p=0;
    do {
      if (cache_ovf_mask_b & 0x1u) { hexit_relative_step_b[p++] = 0; }
      if (cache_ovf_mask_b & 0x2u) { hexit_relative_step_b[p++] = 1; }
      if (cache_ovf_mask_b & 0x4u) { hexit_relative_step_b[p++] = 2; }
      if (cache_ovf_mask_b & 0x8u) { hexit_relative_step_b[p++] = 3; }

      if (cache_ovf_mask_b & 0x10u) { hexit_relative_step_b[p++] = 4; }
      if (cache_ovf_mask_b & 0x20u) { hexit_relative_step_b[p++] = 5; }
      if (cache_ovf_mask_b & 0x40u) { hexit_relative_step_b[p++] = 6; }
      if (cache_ovf_mask_b & 0x80u) { hexit_relative_step_b[p++] = 7; if (p>=8) { break; } }

      if (cache_ovf_mask_b & 0x100u) { hexit_relative_step_b[p++] = 8; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x200u) { hexit_relative_step_b[p++] = 9; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x400u) { hexit_relative_step_b[p++] = 10; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x800u) { hexit_relative_step_b[p++] = 11; if (p>=8) { break; } }

      if (cache_ovf_mask_b & 0x1000u) { hexit_relative_step_b[p++] = 12; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x2000u) { hexit_relative_step_b[p++] = 13; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x4000u) { hexit_relative_step_b[p++] = 14; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x8000u) { hexit_relative_step_b[p++] = 15; if (p>=8) { break; } }

      if (cache_ovf_mask_b & 0x10000u) { hexit_relative_step_b[p++] = 16; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x20000u) { hexit_relative_step_b[p++] = 17; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x40000u) { hexit_relative_step_b[p++] = 18; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x80000u) { hexit_relative_step_b[p++] = 19; if (p>=8) { break; } }

      if (cache_ovf_mask_b & 0x100000u) { hexit_relative_step_b[p++] = 20; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x200000u) { hexit_relative_step_b[p++] = 21; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x400000u) { hexit_relative_step_b[p++] = 22; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x800000u) { hexit_relative_step_b[p++] = 23; if (p>=8) { break; } }

      if (cache_ovf_mask_b & 0x1000000u) { hexit_relative_step_b[p++] = 24; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x2000000u) { hexit_relative_step_b[p++] = 25; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x4000000u) { hexit_relative_step_b[p++] = 26; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x8000000u) { hexit_relative_step_b[p++] = 27; if (p>=8) { break; } }

      if (cache_ovf_mask_b & 0x10000000u) { hexit_relative_step_b[p++] = 28; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x20000000u) { hexit_relative_step_b[p++] = 29; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x40000000u) { hexit_relative_step_b[p++] = 30; if (p>=8) { break; } }
      if (cache_ovf_mask_b & 0x80000000u) { hexit_relative_step_b[p++] = 31; if (p>=8) { break; } }

    } while (0);


#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("  hexit:\n");
      for (i=0; i<8; i++) {
        printf("  [%i] %i, [%i] %i\n",
            hexit_relative_step_a[i], hexit_a[i],
            hexit_relative_step_b[i], hexit_b[i]);
      }
    }
#endif


    // Do a zipper match to count the number of cache overflow hits.
    // Tile variants that overflow from the cache will be picked
    // up by the overflow count.
    //
    for (i=0, j=0; (i<8) && (j<8); ) {

      if ((hexit_relative_step_a[i] < 0) || (hexit_relative_step_b[j] < 0))  { break; }

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("  a[%i+%i=%i] %i (%i) ... b[%i+%i=%i] %i (%i)\n",
                start_block_tile, hexit_relative_step_a[i],
                start_block_tile + hexit_relative_step_a[i],
                hexit_a[i], i,

                start_block_tile , hexit_relative_step_b[j],
                start_block_tile + hexit_relative_step_b[j],
                hexit_b[j], j );
      }
#endif


      if (hexit_relative_step_a[i] < hexit_relative_step_b[j]) {

#ifdef CONC_DEBUG
        if (loc_debug) { printf("  >\n"); }
#endif

        i++; continue;
      }
      if (hexit_relative_step_a[i] > hexit_relative_step_b[j]) {

#ifdef CONC_DEBUG
        if (loc_debug) { printf("  <\n"); }
#endif

        j++; continue;
      }
      if (hexit_relative_step_a[i] == hexit_relative_step_b[j]) {
        if ((hexit_a[i] > 0) && (hexit_a[i] < 0xf) &&
            (hexit_b[j] > 0) && (hexit_b[j] < 0xf) &&
            (hexit_a[i]==hexit_b[j])) {

          if ((ii > n_q) && (ii < n_q_end)) {
            match++;
          }

          else if ( (ii == n_q) && (ii < n_q_end) &&
                    ((start_block_tile + hexit_relative_step_a[i]) >= start_tile_step) &&
                    ((start_block_tile + hexit_relative_step_b[j]) >= start_tile_step) ) {
            match++;
          }

          else if ( (ii > n_q) && (ii == n_q_end) &&
                    ((start_block_tile + hexit_relative_step_a[i]) <= end_tile_step_inc) &&
                    ((start_block_tile + hexit_relative_step_b[j]) <= end_tile_step_inc) ) {
            match++;
          }

          else if ( (ii == n_q) && (ii == n_q_end) &&
                    ((start_block_tile + hexit_relative_step_a[i]) >= start_tile_step) &&
                    ((start_block_tile + hexit_relative_step_b[j]) >= start_tile_step) &&
                    ((start_block_tile + hexit_relative_step_a[i]) <= end_tile_step_inc) &&
                    ((start_block_tile + hexit_relative_step_b[j]) <= end_tile_step_inc) ) {
            match++;
          }
          // else we skip over and don't count the match
          // as it falls outside of the selected window

#ifdef CONC_DEBUG
          if (loc_debug) {
            if ( ((ii > n_q) && (ii < n_q_end)) ||
                 ( (ii == n_q) && (ii < n_q_end) &&
                   ((start_block_tile + hexit_relative_step_a[i]) >= start_tile_step) &&
                   ((start_block_tile + hexit_relative_step_b[j]) >= start_tile_step)) ||
                 ( (ii > n_q) && (ii == n_q_end) &&
                   ((start_block_tile + hexit_relative_step_a[i]) <= end_tile_step_inc) &&
                   ((start_block_tile + hexit_relative_step_b[j]) <= end_tile_step_inc)) ||
                 ( (ii == n_q) && (ii == n_q_end) &&
                   ((start_block_tile + hexit_relative_step_a[i]) >= start_tile_step) &&
                   ((start_block_tile + hexit_relative_step_b[j]) >= start_tile_step) &&
                   ((start_block_tile + hexit_relative_step_a[i]) <= end_tile_step_inc) &&
                   ((start_block_tile + hexit_relative_step_b[j]) <= end_tile_step_inc) ) ) {

              printf("MATCH %i+%i (match count %i)\n", start_block_tile+hexit_relative_step_a[i], 1, match);
              printf("  a[%i] %i == b[%i] %i ++\n",
                  hexit_relative_step_a[i], hexit_a[i],
                  hexit_relative_step_b[j], hexit_b[j]);
            }
            else {
              printf("ignored-match %i+%i (match count %i) (falls outside window)\n", start_block_tile+hexit_relative_step_a[i], 1, match);
              printf("  ignored: a[%i] %i == b[%i] %i ++\n",
                  hexit_relative_step_a[i], hexit_a[i],
                  hexit_relative_step_b[j], hexit_b[j]);
            }
          }
#endif

        }

#ifdef CONC_DEBUG
        if (loc_debug) { printf("  ==\n"); }
#endif

        i++;
        j++;
      }
    }

#ifdef CONC_DEBUG
    if (loc_debug) { printf("  i%i j%i\n", i, j); }
#endif

    if (i==8) {
      for (; j<8; j++) {

        if (hexit_relative_step_b[j] < 0) { break; }
        if ((hexit_b[j] == 0) || (hexit_b[j] >= 0xf)) { continue; }

        if ( ( (ii == n_q) &&
               ((start_block_tile + hexit_relative_step_b[j]) < start_tile_step) )
             ||
             ( (ii == n_q_end) &&
               ((start_block_tile + hexit_relative_step_b[j]) > end_tile_step_inc) )
             ) {
          //do nothing

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("  skipping adding spill b %i %i (j%i, rel:%i) (outside of window)\n",
                hexit_relative_step_b[j] + start_block_tile,
                hexit_b[j], j, hexit_relative_step_b[j]);
          }
#endif

        }
        else {
          spillover_b.push_back(hexit_relative_step_b[j] + start_block_tile);
          spillover_b.push_back(hexit_b[j]);

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("  adding spill b %i %i (j%i, rel:%i)\n",
                hexit_relative_step_b[j] + start_block_tile,
                hexit_b[j], j, hexit_relative_step_b[j]);
          }
#endif

        }

      }
    }
    else if (j==8) {
      for (; i<8; i++) {

        if (hexit_relative_step_a[i] < 0) { break; }
        if ((hexit_a[i] == 0) || (hexit_a[i] >= 0xf)) { continue; }

        if ( ( (ii == n_q) &&
               ((start_block_tile + hexit_relative_step_a[i]) < start_tile_step) )
            ||
             ( (ii == n_q_end) &&
               ((start_block_tile + hexit_relative_step_a[i]) > end_tile_step_inc) )
             ) {

          //do nothing

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("  skipping adding spill a %i %i (i%i) (outside of window)\n",
                hexit_relative_step_a[i] + start_block_tile,
                hexit_a[i], i);
          }
#endif

        }
        else {

          spillover_a.push_back(hexit_relative_step_a[i] + start_block_tile);
          spillover_a.push_back(hexit_a[i]);

#ifdef CONC_DEBUG
          if (loc_debug) {
            printf("  adding spill a %i %i (i%i)\n",
                hexit_relative_step_a[i] + start_block_tile,
                hexit_a[i], i);
          }
#endif

        }

      }
    }

#ifdef CONC_DEBUG
    if (loc_debug) { printf("  +match %i, +tot %i\n", match, tot); }
#endif


		//  __ _ _ ___ ______ _____ _____ _ _
		// / _| '_/ _ (_-<_-</ _ \ V / -_) '_|
		// \__|_| \___/__/__/\___/\_/\___|_|
		//
    // crossover match calculations.
    //
    // calculate the number matched from the elements in the
    // cache to the elements in the `Overflow` vectors for
    // dataset a to b and vice versa.
    //

#ifdef CONC_DEBUG
		if (loc_debug) {
      printf("  ## ii %i, (ii+1) = %i, strideoffset[%i] %i\n",
          (int)ii, (int)((ii+1)),
          (int)tilepath_idx, (int)(a->StrideOffset[tilepath_idx]));
    }
#endif


    // If we've crossed a tilepath boundary or we're at the end
    // of our window, check the spillover to overflow matches.
    //
    if ( ((ii+1) >= a->StrideOffset[tilepath_idx]) ||
         (ii == n_q_end) ) {

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("CROSS ii %i (ii+1=%i) >= a->StrideOffset[tilepath_idx %i] %i\n",
            (int)ii, (int)(ii+1),
            (int)tilepath_idx,
            (int)a->StrideOffset[tilepath_idx]);
      }
#endif

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("SPILLOVER CALC\n");
        printf("---\n");
        for (i=0; i<spillover_a.size(); i+=2) {
          printf("  spillover_a[%i] step:%i tilemapid:%i\n", i, spillover_a[i], spillover_a[i+1]);
        }
        printf("---\n");
        for (i=0; i<spillover_b.size(); i+=2) {
          printf("  spillover_b[%i] step:%i tilemapid:%i\n", i, spillover_b[i], spillover_b[i+1]);
        }
        printf("---\n");
      }
#endif

      // Transfer cache spillover to unpacked spillover for easy comparison with
      // the Overflow arrays
      //
      spillover_knot_a.clear();
      for (i=0; i<spillover_a.size(); i+=2) {

        // if the knot passes over the end of the last tilestep in our window, ignore it
        //
        if (ii == n_q_end) {

          // We only check that the anchor step is in the window range, allowing for
          // matches for knots that begin inside the window but end outside the window.
          // Though technically 'wrong' this simplifies calculations since we don't need
          // to worry about edge cases of knots spilling over into neighboring blocks.
          //
          if (spillover_a[i] > end_tile_step_inc) {

#ifdef CONC_DEBUG
            if (loc_debug) {
              printf("  skipping ovf16_a %i %i\n", spillover_a[i], spillover_a[i+1]);
            }
#endif

            continue;
          }
        }

        // tile offset holds the relative offste from the beginning of the tile step
        // for the tile knot held in the tilemap
        //
        tile_offset=0;
        for ( j=tilemap->offset[ spillover_a[i+1]-1 ]; j<tilemap->offset[spillover_a[i+1]]; j++) {
          spillover_knot_a.push_back(spillover_a[i] + tile_offset);
          spillover_knot_a.push_back( tilemap->variant[0][j] );
          spillover_knot_a.push_back( tilemap->variant[1][j] );

          spillover16_knot_a.push_back(spillover_a[i] + tile_offset);
          u16 = ((tilemap->variant[0][j] < 0) ? OVF16_MAX : (uint16_t)tilemap->variant[0][j]);
          spillover16_knot_a.push_back( tilemap->variant[0][j] );
          u16 = ((tilemap->variant[1][j] < 0) ? OVF16_MAX : (uint16_t)tilemap->variant[1][j]);
          spillover16_knot_a.push_back( tilemap->variant[1][j] );

          tile_offset++;
        }
      }
      spillover_a.clear();

      spillover_knot_b.clear();
      for (i=0; i<spillover_b.size(); i+=2) {

        // if the knot passes over the end of the last tilestep in our window, ignore it
        //
        if (ii == n_q_end) {

          // We only check that the anchor step is in the window range, allowing for
          // matches for knots that begin inside the window but end outside the window.
          // Though technically 'wrong' this simplifies calculations since we don't need
          // to worry about edge cases of knots spilling over into neighboring blocks.
          //
          if (spillover_b[i] > end_tile_step_inc) {

#ifdef CONC_DEBUG
            if (loc_debug) {
              printf("  skipping ovf16_b %i %i\n", spillover_b[i], spillover_b[i+1]);
            }
#endif

            continue;
          }
        }

        // tile offset holds the relative offste from the beginning of the tile step
        // for the tile knot held in the tilemap
        //
        tile_offset=0;
        for ( j=tilemap->offset[ spillover_b[i+1]-1 ]; j<tilemap->offset[spillover_b[i+1]]; j++) {
          spillover_knot_b.push_back(spillover_b[i] + tile_offset);
          spillover_knot_b.push_back( tilemap->variant[0][j] );
          spillover_knot_b.push_back( tilemap->variant[1][j] );

          spillover16_knot_b.push_back((uint16_t)(spillover_b[i] + tile_offset));
          u16 = ((tilemap->variant[0][j] < 0) ? OVF16_MAX : (uint16_t)tilemap->variant[0][j]);
          spillover16_knot_b.push_back( u16 );
          u16 = ((tilemap->variant[1][j] < 0) ? OVF16_MAX : (uint16_t)tilemap->variant[1][j]);
          spillover16_knot_b.push_back( u16 );

          tile_offset++;
        }
      }
      spillover_b.clear();

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("spillover_knot_a:\n");
        for (i=0; i<spillover_knot_a.size(); i+=3) {
          printf(" [%i] %i %i %i (%i %i %i)\n",
              i,
              spillover_knot_a[i],
              spillover_knot_a[i+1],
              spillover_knot_a[i+2],

              (int)spillover16_knot_a[i],
              (int)spillover16_knot_a[i+1],
              (int)spillover16_knot_a[i+2]
              );
        }
        printf("\n");
        printf("---\n");

        printf("spillover_knot_b:\n");
        for (i=0; i<spillover_knot_b.size(); i+=3) {
          printf(" [%i] %i %i %i (%i %i %i)\n",
              i,
              spillover_knot_b[i],
              spillover_knot_b[i+1],
              spillover_knot_b[i+2],

              (int)spillover16_knot_b[i],
              (int)spillover16_knot_b[i+1],
              (int)spillover16_knot_b[i+2]);
        }
        printf("\n");
        printf("  non-overflow match: %i (%i)\n", match, tot);
      }
#endif

      // Find start and end of the tilepath for
      // each `Overflow` vector.
      //

      iistart=0;
      if (tilepath_idx>0) {
        iistart=a->OverflowOffset[tilepath_idx-1];
      }

      jjstart=0;
      if (tilepath_idx>0) {
        jjstart=b->OverflowOffset[tilepath_idx-1];
      }

      end_noninc_a = a->OverflowOffset[tilepath_idx];
      end_noninc_b = b->OverflowOffset[tilepath_idx];

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("\nA OVERFLOW[%i-%i] ** SPILLOVER B [%i-%i]\n",
            (int)iistart, (int)end_noninc_a,
            0, (int)spillover16_knot_b.size());
      }
#endif

      t_match=0; t_tot=0;

      overflow_concordance16_z( &t_match, &t_tot,
                                &(a->Overflow[0]), iistart, end_noninc_a,
                                &(spillover16_knot_b[0]), 0, (int)spillover16_knot_b.size(),
                                cgf_opt);

      // total has already been counted
      //
      match += t_match;

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("a->Overflow to spillover16_knot_b match %i\n", t_match);
      }


      if (loc_debug) {
        printf("\nB OVERFLOW[%i-%i] ** SPILLOVER A[%i-%i]\n",
            (int)jjstart, (int)end_noninc_b,
            0, (int)spillover16_knot_a.size()
            );
      }
#endif

      t_match=0; t_tot = 0;

      overflow_concordance16_z( &t_match, &t_tot,
                                &(spillover16_knot_a[0]), 0, (int)spillover16_knot_a.size(),
                                &(b->Overflow[0]), jjstart, end_noninc_b,
                                cgf_opt);

      // total has already been counted
      //
      match += t_match;

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("spillover16_knot_a to b->Overflow match %i\n", t_match);
      }
#endif


      // Reset relevant tile path and tile step information
      //
      start_block_tile=-32;
      tilepath_idx++;
      spillover16_knot_a.clear();
      spillover16_knot_b.clear();

    } // tilepath crossed boundary calculations

  } // for ii=n_q ...

	//           __            __
	//  _____ __/ _|  _____ __/ _|
	// / _ \ V /  _| / _ \ V /  _|
	// \___/\_/|_|   \___/\_/|_|
	//

  for (tilepath_idx = start_tile_path;
       tilepath_idx <= end_tile_path_inc;
       tilepath_idx++) {

    iistart=jjstart=0;
    if (tilepath_idx>0) {
      iistart=a->OverflowOffset[tilepath_idx-1];
      jjstart=b->OverflowOffset[tilepath_idx-1];
    }

    end_noninc_a = a->OverflowOffset[tilepath_idx];
    end_noninc_b = b->OverflowOffset[tilepath_idx];

    // If we've started on the beginning tile path,
    // run forward till we get to the first tile step
    // under ocnsiderationg.
    //
    if (tilepath_idx == start_tile_path) {

      while (iistart < end_noninc_a) {
        if ((int)(a->Overflow[iistart]) >= start_tile_step) { break; }
        iistart+=3;
      }

      while (jjstart < end_noninc_b) {
        if ((int)(b->Overflow[jjstart]) >= start_tile_step) { break; }
        jjstart+=3;
      }

    }

    if (tilepath_idx == end_tile_path_inc) {

      while (end_noninc_a > iistart) {
        if ((int)((a->Overflow[end_noninc_a-3])) <= end_tile_step_inc) { break; }

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  a->Overflow[%i] %i\n",
              (int)end_noninc_a,
              (int)(a->Overflow[end_noninc_a-3]));
        }
#endif

        end_noninc_a -= 3;
      }

#ifdef CONC_DEBUG
      if (loc_debug) {
        printf("---\n");
      }
#endif

      while (end_noninc_b > jjstart) {
        if ((int)((b->Overflow[end_noninc_b-3])) <= end_tile_step_inc) { break; }

#ifdef CONC_DEBUG
        if (loc_debug) {
          printf("  b->Overflow[%i] %i\n",
              (int)end_noninc_b,
              (int)(b->Overflow[end_noninc_b-3]));
        }
#endif

        end_noninc_b -= 3;
      }

    }

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("\n\nOVERFLOW OVERFLOW CONCORDANCE\n----\n\n");
      printf("  iistart %i, end_noninc_a %i\n", (int)iistart, (int)end_noninc_a);
      printf("  jjstart %i, end_noninc_b %i\n", (int)jjstart, (int)end_noninc_b);
    }
#endif

    t_match=0; t_tot=0;

    overflow_concordance16_z( &t_match, &t_tot,
                              &(a->Overflow[0]), iistart, end_noninc_a,
                              &(b->Overflow[0]), jjstart, end_noninc_b,
                              cgf_opt);
    match += t_match;

#ifdef CONC_DEBUG
    if (loc_debug) {
      printf("  ovf ovf conc t_match %i\n", t_match);
    }
#endif

  }

#ifdef CONC_DEBUG
  if (loc_debug) { printf("  fin: match %i, tot %i\n", match, tot); }
#endif

  // all done!
  //
  *r_match = match;
  *r_tot =tot;

  return 0;
}
