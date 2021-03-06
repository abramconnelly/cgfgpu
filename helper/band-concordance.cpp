#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>

#include <vector>
#include <string>
#include <map>

#define VERSION_STR "0.1.0"

int VERBOSE_MATCH;

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

  int end_bound=0;

  if (s<=0) { s=0; }
  if (n<=0) { end_noninc = a.v[0].size(); }
  else { end_noninc = s+n; }
  if (end_noninc >= a.v[0].size()) { end_noninc = a.v[0].size(); }
  n = end_noninc - s;

  idx_a = s;
  idx_b = s;

  // skip past initial knot
  //
  for (i=s; i<(s+n); i++) {
    if ( (a.v[0][idx_a] >= 0) &&
         (a.v[1][idx_a] >= 0) &&
         (b.v[0][idx_b] >= 0) &&
         (b.v[1][idx_b] >= 0) ) {
      break;
    }
    idx_a++;
    idx_b++;
  }


  tileband_knot(a, idx_a, &knot_a_len, &knot_a_loq);
  tileband_knot(b, idx_b, &knot_b_len, &knot_b_loq);

  while ( (idx_a < (s+n)) &&
          (idx_b < (s+n)) ) {

    if (VERBOSE_MATCH) {
    printf("idx_a %i+%i (%i), idx_b %i+%i (%i) (%i / %i)\n",
        idx_a, knot_a_len, knot_a_loq,
        idx_b, knot_b_len, knot_b_loq,
       match, tot );
    }

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
            if (is_match) {

              if (VERBOSE_MATCH) { printf("MATCH %i+%i\n", idx_a, knot_a_len); }

              match++;
            }
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
          if (is_match) {

            if (VERBOSE_MATCH) { printf("MATCH %i+%i\n", idx_a, knot_a_len); }

            match++;
          }
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

static struct option long_options[] = {
  {"step",          required_argument,        NULL, 's'},
  {"endstep",       required_argument,        NULL, 'S'},
  {"help",                no_argument,        NULL, 'h'},
  {"version",             no_argument,        NULL, 'v'},
  {"verbose",             no_argument,        NULL, 'V'},
  {"all-endpoints",       no_argument,        NULL, 'A'},
  {0,0,0,0}
};

void show_version() {
  printf("%s", VERSION_STR);
}

void show_help() {
  printf("\n");
  printf("version: ");
  show_version();
  printf("\n\n");
  printf("usage:\n\n");
  printf("  band-concordance [-h] [-v] [-V] [-s s] [-S S] bandA bandB\n");
  printf("\n");
  printf("  [-s s]  tile step start\n");
  printf("  [-S S]  tile step end (inclusive)\n");
  printf("  [-n n]  tile step count\n");
  printf("  [-A]    do a concordance on all endpoints in path\n");
  printf("  [-h]    show help (this screen)\n");
  printf("  [-v]    verbose\n");
  printf("  [-V]    version\n");
  printf("\n");
}

int main(int argc, char **argv) {
  int i, j, k;
  FILE *fp;
  tileband_t tileband_a, tileband_b;
  std::string ifn_a, ifn_b;
  int match=0, tot=0;

  int opt, option_index;
  int verbose_flag = 0;

  int start_tilestep=0, end_tilestep_inc = -1;
  int n_tilestep = -1;
  int hiq_flag = 1;
  int all_endpoints=0;

  int s_offset,n_offset;

  VERBOSE_MATCH=0;

  while ((opt = getopt_long(argc, argv, "hvVs:S:n:A", long_options, &option_index))!=-1) switch (opt) {
    case 0:
      fprintf(stderr, "invalid argument");
      exit(-1);
      break;
    case 's': start_tilestep = atoi(optarg); break;
    case 'S': end_tilestep_inc = atoi(optarg); break;
    case 'n': n_tilestep = atoi(optarg); break;
    case 'v': verbose_flag = 1; break;
    case 'A': all_endpoints=1; break;
    case 'V': show_version(); exit(0); break;
    default:
    case 'h': show_help(); exit(0); break;
  }

  if (verbose_flag) {
    VERBOSE_MATCH=1;
  }

  if ((argc-optind)>=2) {
    ifn_a = argv[optind];
    ifn_b = argv[optind+1];
  }
  else {
    printf("provide input band files\n");
    show_help();
    exit(-1);
  }

  if (end_tilestep_inc >= 0) {
    n_tilestep = end_tilestep_inc - start_tilestep + 1;
  }

  /*
  if (argc!=3) {
    printf("provide two band files\n");
    exit(0);
  }
  */

  if (!(fp = fopen(ifn_a.c_str(), "r"))) {
    perror(ifn_a.c_str());
    exit(-1);
  }
  tileband_read(tileband_a, fp);
  fclose(fp);

  if (!(fp = fopen(ifn_b.c_str(), "r"))) {
    perror(ifn_b.c_str());
    exit(-1);
  }
  tileband_read(tileband_b, fp);
  fclose(fp);


  if (VERBOSE_MATCH) {
    tileband_print(tileband_a);
    printf("\n---\n\n");
    tileband_print(tileband_b);
  }

  if (all_endpoints) {

    if (n_tilestep < 0) {
      n_tilestep = (int)tileband_a.v[0].size();
    }

    for (s_offset=0; s_offset < (n_tilestep-1); s_offset++) {
      for (n_offset=1; n_offset < (n_tilestep-s_offset+1); n_offset++) {

        match=0; tot=0;
        tileband_concordance(tileband_a, tileband_b,
            start_tilestep + s_offset, n_offset,
            hiq_flag,
            &match, &tot);

        printf("[%i+%i] match: %i, total: %i\n",
            start_tilestep+s_offset, n_offset,
            match, tot);

      }
    }

  }

  else {

    tileband_concordance(tileband_a, tileband_b,
        start_tilestep, n_tilestep,
        hiq_flag,
        &match, &tot);

    printf("match: %i, total: %i\n", match, tot);

  }

}
