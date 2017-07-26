#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <vector>
#include <string>
#include <map>

typedef struct tileband_type {
  std::vector< int > v[2];
  std::vector< std::vector< int > > noc_v[2];
} tileband_t;

void tileband_knot(tileband_t &t, int knot_start, int *knot_len, int *loc_knot);

int tileband_concordance(tileband_t &a, tileband_t &b, int s, int n, int hiq_flag, int *r_match, int *r_tot) {
  int i, j, k;
  int end_noninc=0;
  int match=0, tot=0;

  int idx_a, idx_b;
  int knot_a_len, knot_b_len;
  int knot_a_loq, knot_b_loq;
  int is_match=0;

  if (s<=0) { s=0; }
  if (n<=0) { end_noninc = a.v[0].size(); }
  else { end_noninc = s+n; }
  if (end_noninc >= a.v[0].size()) { end_noninc = a.v[0].size(); }
  n = end_noninc - s;

  idx_a = s;
  idx_b = s;


  tileband_knot(a, idx_a, &knot_a_len, &knot_a_loq);
  tileband_knot(b, idx_b, &knot_b_len, &knot_b_loq);

  while ( (idx_a < (s+n)) &&
          (idx_b < (s+n)) ) {

    printf("idx_a %i+%i (%i), idx_b %i+%i (%i) (%i / %i)\n",
        idx_a, knot_a_len, knot_a_loq,
        idx_b, knot_b_len, knot_b_loq,
       match, tot );

    if (idx_a == idx_b) {

      if (hiq_flag) {

        if ((knot_a_loq==0) && (knot_b_loq==0)) {
          tot++;

          if ((knot_a_len == knot_b_len) &&
              (knot_a_loq == 0) &&
              (knot_b_loq == 0)) {
            is_match=1;
            for (i=idx_a; i<(idx_a+knot_a_len); i++) {
              if ((a.v[0][i] != b.v[0][i]) ||
                  (a.v[1][i] != b.v[1][i])) {
                is_match=0;
                break;
              }
            }
            if (is_match) { match++; }
          }
        }

      } else {
        tot++;

        if (knot_a_len == knot_b_len) {
          is_match=1;
          for (i=idx_a; i<(idx_a+knot_a_len); i++) {
            if ((a.v[0][i] != b.v[0][i]) ||
                (a.v[1][i] != b.v[1][i])) {
              is_match=0;
              break;
            }
          }
          if (is_match) { match++; }
        }

      }

      idx_a += knot_a_len;
      idx_b += knot_b_len;

      if ((idx_a < (s+n)) && (idx_b < (s+n))) {
        tileband_knot(a, idx_a, &knot_a_len, &knot_a_loq);
        tileband_knot(b, idx_b, &knot_b_len, &knot_b_loq);
      }
      continue;
    }

    if (idx_a < idx_b) {
      idx_a += knot_a_len;
      if (idx_a < (s+n)) {
        tileband_knot(a, idx_a, &knot_a_len, &knot_a_loq);
      }
      continue;
    }

    if (idx_a > idx_b) {
      idx_b += knot_b_len;
      if (idx_b < (s+n)) {
        tileband_knot(b, idx_b, &knot_b_len, &knot_b_loq);
      }
      continue;
    }

  }

  *r_match = match;
  *r_tot = tot;

}

void tileband_knot(tileband_t &t, int knot_start, int *knot_len, int *loc_knot) {
  int i, n;
  int kl=1;

  *loc_knot=0;

  if ((t.noc_v[0][knot_start].size()>0) || 
      (t.noc_v[1][knot_start].size()>0)) {
    *loc_knot = 1;
  }

  n = (int)t.v[0].size();
  for (i=knot_start+1; i<n; i++) {

    if ((t.v[0][i] < 0) || (t.v[1][i] < 0)) {
      kl++;

      if ((t.noc_v[0][i].size()>0) || 
          (t.noc_v[1][i].size()>0)) {
        *loc_knot = 1;
      }

      continue;
    }

    break;
  }

  *knot_len = kl;
}

void tileband_print(tileband_t &t) {
  int i, j, k, n;

  n=(int)t.v[0].size();

  printf("[");
  for (i=0; i<t.v[0].size(); i++) { printf(" %i", t.v[0][i]); }
  printf("]\n");

  printf("[");
  for (i=0; i<t.v[1].size(); i++) { printf(" %i", t.v[1][i]); }
  printf("]\n");

  printf("[");
  for (i=0; i<t.noc_v[0].size(); i++) {

    printf("[");
    for (j=0; j<t.noc_v[0][i].size(); j++) {
      printf(" %i", t.noc_v[0][i][j]);
    }
    printf(" ]");

  }
  printf("]\n");

  printf("[");
  for (i=0; i<t.noc_v[1].size(); i++) {

    printf("[");
    for (j=0; j<t.noc_v[1][i].size(); j++) {
      printf(" %i", t.noc_v[1][i][j]);
    }
    printf(" ]");

  }
  printf("]\n");


}

int tileband_read(tileband_t &t, FILE *fp) {
  int ch, band_state=0, paren=0;
  std::string buf;
  std::vector< int > loq_v;

  while (!feof(fp)) {
    ch = fgetc(fp);
    if ((ch==EOF) || (ch=='\n')) {
      band_state++;
      continue;
    }

    if (ch=='[') { paren++; continue; }
    if (ch==']') {
      paren--;

      if ((band_state==0) || (band_state==1)) {
        if (buf.size()>0) {
          t.v[band_state].push_back(atoi(buf.c_str()));
        }
      }

      else if ((band_state==2) || (band_state==3)) {
        if (paren==1) {
          t.noc_v[band_state-2].push_back(loq_v);
          loq_v.clear();
        }
      }

      buf.clear();
      continue;
    }

    if (ch==' ') {

      if ((band_state==0) || (band_state==1)) {
        if (buf.size()>0) {
          t.v[band_state].push_back(atoi(buf.c_str()));
        }
      }
      else if ((band_state==2) || (band_state==3)) {
        if (buf.size()>0) {
          loq_v.push_back(atoi(buf.c_str()));
        }
      }
      buf.clear();
      continue;
    }

    buf.push_back((char)ch);
  }

  return 0;
}

int main(int argc, char **argv) {
  int i, j, k;
  FILE *fp;
  tileband_t tileband_a, tileband_b;
  std::string ifn_a, ifn_b;
  int match=0, tot=0;

  if (argc!=3) {
    printf("provide two band files\n");
    exit(0);
  }

  if (!(fp = fopen(argv[1], "r"))) {
    perror(argv[1]);
    exit(-1);
  }
  tileband_read(tileband_a, fp);
  fclose(fp);

  if (!(fp = fopen(argv[2], "r"))) {
    perror(argv[2]);
    exit(-1);
  }
  tileband_read(tileband_b, fp);
  fclose(fp);


  tileband_print(tileband_a);
  printf("\n---\n\n");
  tileband_print(tileband_b);


  tileband_concordance(tileband_a, tileband_b,
      0, 0, 1,
      &match, &tot);

  printf("match: %i, total: %i\n", match, tot);

}
