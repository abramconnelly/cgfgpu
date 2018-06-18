#include "cgf4.hpp"

#define cleanup_err() do { ret=1; show_help(stderr); goto cgf_cleanup; } while (0);
#define cleanup_fail() do { ret=-2; goto cgf_cleanup; } while (0);
#define cleanup_ok() do { ret=0; show_help(stderr); goto cgf_cleanup; } while (0);
#define cleanup() do { show_help(stderr); goto cgf_cleanup; } while (0);


static struct option long_options[] = {
  {"header",              no_argument,        NULL, 'H'},
  {"create-container",    no_argument,        NULL, 'C'},
  {"info",                no_argument,        NULL, 'I'},
  {"show-all",            no_argument,        NULL, 'A'},
  {"help",                no_argument,        NULL, 'h'},
  {"version",             no_argument,        NULL, 'v'},
  {"verbose",             no_argument,        NULL, 'V'},
  {"ez-print",            no_argument,        NULL, 'Z'},
  {"hiq",                 no_argument,        NULL, 'q'},
  {"match",               no_argument,        NULL, 'm'},
  {"sanity",              no_argument,        NULL, 'Y'},

  {"band",                required_argument,  NULL, 'b'},
  {"encode",              required_argument,  NULL, 'e'},

  {"tilepath",            required_argument,  NULL, 'p'},
  {"endtilepath",         required_argument,  NULL, 'P'},

  {"tilestep",            required_argument,  NULL, 's'},
  {"endtilestep",         required_argument,  NULL, 'S'},

  {"genotype",            no_argument,        NULL, '\0'},
  {"no-gtz",              no_argument,        NULL, '\0'},
  {"all-pairs",           no_argument,        NULL, '\0'},
  {"show-stats",          no_argument,        NULL, '\0'},
  {"repeat",              required_argument,  NULL, '\0'},

  {"input",               required_argument,  NULL, 'i'},
  {"output",              required_argument,  NULL, 'o'},
  {"tilemap",             required_argument,  NULL, 't'},
  {"version-opt",         required_argument,  NULL, 'T'},
  {"library-version",     required_argument,  NULL, 'L'},
  {0,0,0,0}
};

void show_help(FILE *fp) {
  fprintf(fp, "CGF Tool.  A tool used to inspect and edit Compact Genome Format (CGF) v4 files.\n");
  fprintf(fp, "Version: %s\n", CGF_VERSION);
  fprintf(fp, "\n");
  fprintf(fp, "usage:\n  cgf4 [-H] [-b tilepath] [-e tilepath] [-i ifn] [-o ofn] [-h] [-v] [-V] [ifn]\n");
  fprintf(fp, "\n");
  fprintf(fp, "  [-H|--header]               show header\n");
  fprintf(fp, "  [-C|--create-container]     create empty container\n");
  fprintf(fp, "  [-I|--info]                 print basic information about CGF file\n");
  fprintf(fp, "  [-b|--band tilepath]        output band for tilepath\n");
  fprintf(fp, "  [-F fill_level]             bit vector for band fillin (1 canon, 2 cache, 4 ovf, 8 noc, 0xff default)\n");
  fprintf(fp, "  [-m|--match]                run concordance on a pair of cgf files (must provide two cgf files)\n");
  fprintf(fp, "  [-e|--encode tilepath]      input tilepath band and add it to file, overwriting if it already exists\n");

  fprintf(fp, "  [-p|--tilepath tilepath]    tilepath (start)\n");
  fprintf(fp, "  [-P|--endtilepath tilepath] end tilepath\n");
  fprintf(fp, "  [-s|--tilestep tilestep]    tilestep (start)\n");
  fprintf(fp, "  [-S|--endtilestep tilestep] end tilestep\n");
  fprintf(fp, "  [--genotype]                genotype flag\n");
  fprintf(fp, "  [--no-gtz]                  do not use zlib compression for tile genotype position information\n");
  fprintf(fp, "  [--all-pairs]               all pairs concordance\n");
  fprintf(fp, "  [--show-stats]              show stats\n");


  fprintf(fp, "  [-i|--input ifn]            input file (CGF)\n");
  fprintf(fp, "  [-o|--output ofn]           output file (CGF)\n");
  fprintf(fp, "  [-A|--show-all]             show all tilepaths\n");
  fprintf(fp, "  [-h|--help]                 show help (this screen)\n");
  fprintf(fp, "  [-v|--version]              show version\n");
  fprintf(fp, "  [-V|--verbose]              set verbose level\n");
  fprintf(fp, "  [-t|--tilemap tfn]          use tilemap file (instead of default)\n");
  fprintf(fp, "  [-Z|--ez-print]             print \"ez\" structure information\n");
  fprintf(fp, "  [-T|--opt-version vopt]     CGF version option.  Must be one of \"default\" or \"noc-inv\"\n");
  fprintf(fp, "  [-L|--library-version lopt] CGF library version option.  Overwrites default value if specified\n");
  fprintf(fp, "  [-U|--update-header]        Update header only\n");
  fprintf(fp, "  [-Y|--sanity]               Run sanity checks\n");
  fprintf(fp, "  [-q|--hiq]                  Only output high quality information (band output)\n");
  fprintf(fp, "\n");
}

void show_version(FILE *fp) {
  fprintf(fp, "version %s\n", CGF_VERSION);
}

#ifdef SAMPLE_STACK_PROFILER
// g++ -rdynamic bt_debug_test.c
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <signal.h>

#include <execinfo.h>
#include <unistd.h>

#include <string>
#include <vector>


void timer_handler(int signum) {

  void *buffers[1024];
  char **strings;
  int i, nret, bufsz = 1024;
  int custom_debug=0;

  nret = backtrace(buffers, bufsz);

  printf("\n---\ntimer: got %i\n", nret);
  fflush(stdout);

  if (custom_debug) {

    strings = backtrace_symbols(buffers, nret);
    if (strings != NULL) {
      for (i=0; i<nret; i++) {
        printf("[%i] %s\n", i, strings[i]);
      }

      free(strings);
    }

  }
  else {
    backtrace_symbols_fd(buffers, nret, STDOUT_FILENO);
    fflush(stdout);
  }

}

void setup_timer(struct sigaction *sa, struct itimerval *timer, void (*timer_func)(int), suseconds_t usec) {
  memset(sa, 0, sizeof(struct sigaction));
  sa->sa_handler = timer_func;
  sigaction(SIGVTALRM, sa, NULL);

  timer->it_value.tv_sec = 0;
  timer->it_value.tv_usec = usec; //e.g. 250000;
  timer->it_interval.tv_sec = 0;
  timer->it_interval.tv_usec = usec; //e.g. 250000;
  setitimer(ITIMER_VIRTUAL, timer, NULL);

}

#endif


int main(int argc, char **argv) {
  int i, j, k, opt, ret=0, idx;
  int n_step=0;
  char buf[1024];

  std::string tilemap_str;
  FILE *ifp=NULL, *ofp=NULL;
  cgf_t *cgf=NULL, *cgf_b=NULL;
  int option_index=0;
  int def_or_nocinv=0;

  int match=0, tot=0;

  cgf_opt_t cgf_opt;

#ifdef SAMPLE_STACK_PROFILER
  struct sigaction sa;
  struct itimerval timer;
  suseconds_t usec;
  usec = 25000;
  setup_timer(&sa, &timer, timer_handler, usec);
#endif


  cgf_opt_init(&cgf_opt);

  while ((opt = getopt_long(argc, argv, "Hb:e:i:o:Ct:T:L:U:hvVAZRIqmp:P:s:S:YF:", long_options, &option_index))!=-1) switch (opt) {
    case 0:
      if (strcmp(long_options[option_index].name, "all-pairs")==0) {
        cgf_opt.all_pairs=1;
      }
      else if (strcmp(long_options[option_index].name, "show-stats")==0) {
        cgf_opt.print_stats=1;
      }
      else if (strcmp(long_options[option_index].name, "repeat")==0) {
        cgf_opt.repeat = atoi(optarg);
      }
      else if (strcmp(long_options[option_index].name, "genotype")==0) {
        cgf_opt.encode_genotype_flag=1;
      }
      else if (strcmp(long_options[option_index].name, "no-gtz")==0) {
        cgf_opt.gtz_flag=0;
      }
      else {
        fprintf(stderr, "invalid option, exiting\n");
        exit(-1);
      }
      break;
    case 'H': cgf_opt.show_header=1; break;
    case 'Y': cgf_opt.run_sanity=1; break;
    case 'C': cgf_opt.create_container=1; break;
    case 'I': cgf_opt.info=1; break;
    case 't': cgf_opt.tilemap=1; cgf_opt.tilemap_fn=strdup(optarg); break;
    case 'b': cgf_opt.show_band=1; cgf_opt.tilepath=atoi(optarg); break;
    case 'e': cgf_opt.encode=1; cgf_opt.tilepath=atoi(optarg); break;
    case 'd': cgf_opt.del=1; cgf_opt.tilepath=atoi(optarg); break;
    case 'i': cgf_opt.ifns.push_back(optarg); break;
    case 'o': cgf_opt.ofn=strdup(optarg); break;

    case 'h': cgf_opt.show_help=1; break;
    case 'v': cgf_opt.show_version=1; break;

    case 'A': cgf_opt.show_all=1; break;
    case 'q': cgf_opt.hiq_only=1; break;

    case 'p': cgf_opt.tilepath=atoi(optarg); break;
    case 'P': cgf_opt.endtilepath=atoi(optarg); break;

    case 's': cgf_opt.tilestep=atoi(optarg); break;
    case 'S': cgf_opt.endtilestep=atoi(optarg); break;

    case 'm': cgf_opt.match=1; break;

    case 'F': cgf_opt.fill_level = (uint32_t)strtoul(optarg,NULL, 10); break;

    case 'V': cgf_opt.verbose=1; break;
    case 'Z': cgf_opt.ez_print=1; break;
    case 'R': cgf_opt.run_test=1; break;
    case 'T': cgf_opt.cgf_version_opt=1;
              cgf_opt.cgf_version_opt_ele.push_back(optarg);
              cgf_opt.update_header = 1;
              break;
    case 'L': cgf_opt.cglf_version_opt=1;
              cgf_opt.cglf_version_opt_ele.push_back(optarg);
              cgf_opt.update_header=1;
              break;
    case 'U': cgf_opt.update_header=1; break;
    default: printf("unknown option"); cleanup_err(); break;
  }

  if (argc>optind) {
    if ((argc-optind)>2) { printf("Extra options specified\n"); cleanup_err(); }
    for (i=0; i<(argc-optind); i++) {
      cgf_opt.ifns.push_back(argv[optind+i]);
    }
  }

  if (cgf_opt.ifns.size()>0) { cgf_opt.ifn = cgf_opt.ifns[0]; }

  if (cgf_opt.show_help) { show_help(stdout); goto cgf_cleanup; }
  if (cgf_opt.show_version) { show_version(stdout); goto cgf_cleanup; }


  // We must have a command.  If not, exit.
  // The 'no-surprise' rule is to show help and exit gracefully when
  // no commands are specified.
  //
  if ((cgf_opt.create_container +
       cgf_opt.encode +
       cgf_opt.del +
       cgf_opt.show_band +
       cgf_opt.show_header +
       cgf_opt.show_all +
       cgf_opt.run_test +
       cgf_opt.info +
       cgf_opt.match +
       cgf_opt.run_sanity +
       cgf_opt.update_header) == 0) {
    cleanup_err();
  }

  // Don't allow more than one command
  //
  if ((cgf_opt.create_container +
       cgf_opt.encode +
       cgf_opt.del +
       cgf_opt.show_band +
       cgf_opt.show_header +
       cgf_opt.show_all +
       cgf_opt.run_test +
       cgf_opt.match +
       cgf_opt.run_sanity +
       cgf_opt.info) > 0) {
    cgf_opt.update_header = 0;
  }

  // Don't allow more than one command
  //
  if ((cgf_opt.create_container +
       cgf_opt.encode +
       cgf_opt.del +
       cgf_opt.show_band +
       cgf_opt.show_header +
       cgf_opt.show_all +
       cgf_opt.run_test +
       cgf_opt.info +
       cgf_opt.match +
       cgf_opt.run_sanity +
       cgf_opt.update_header) != 1) {
    printf("must specify exactly one of show header (-H), show band (-b), encode (-e), delete (-d), create empty container (-C) or update header (-U)\n");

    printf("cc %i\n", cgf_opt.create_container);
    printf("enc %i\n", cgf_opt.encode);
    printf("de %i\n", cgf_opt.del);
    printf("sb %i\n", cgf_opt.show_band);
    printf("sh %i\n", cgf_opt.show_header);
    printf("sa %i\n", cgf_opt.show_all);
    printf("rt %i\n", cgf_opt.run_test);
    printf("in %i\n", cgf_opt.info);
    printf("m %i\n", cgf_opt.match);
    printf("rs %i\n", cgf_opt.run_sanity);
    printf("uh %i\n", cgf_opt.update_header);

    cleanup_err();
  }

  // add to CGF header string
  //
  for (i=0; i<cgf_opt.cgf_version_opt_ele.size(); i++) {
    if (cgf_opt.cgf_version_opt_ele[i] == "noc-inv") {
      cgf_opt.cgf_version_str += ",noc-inv";
    }
  }


  //                  __       _
  //   _______  ___  / /____ _(_)__  ___ ____
  //  / __/ _ \/ _ \/ __/ _ `/ / _ \/ -_) __/
  //  \__/\___/_//_/\__/\_,_/_/_//_/\__/_/
  //

  if (cgf_opt.create_container) {


    if (cgf_opt.ofn.size()==0) {
      if (cgf_opt.ifn.size()>0) { cgf_opt.ofn = cgf_opt.ifn; }
      else { printf("specify output CGF file\n"); cleanup_err(); }
    }

    if (cgf_opt.tilemap_fn.size()==0) {
      tilemap_str = DEFAULT_TILEMAP;
    } else {
      if ((read_tilemap_from_file(tilemap_str, cgf_opt.tilemap_fn.c_str()))==NULL) {
        perror(cgf_opt.tilemap_fn.c_str());
        cleanup_err();
      }
    }

    if (cgf_opt.ofn=="-") { ofp=stdout; }
    else if ((ofp = fopen(cgf_opt.ofn.c_str(), "w"))==NULL) { perror(cgf_opt.ofn.c_str()); cleanup_err(); }

    if (cgf_opt.encode_genotype_flag) {
      cgf_opt.cgf_version_str += ",gt";
    }

    cgf_create_container(
        ofp,
        cgf_opt.cgf_version_str.c_str(),
        cgf_opt.cglf_version_str.c_str(),
        tilemap_str.c_str());
  }

  //                         __
  //   ___ ___  _______  ___/ /__
  //  / -_) _ \/ __/ _ \/ _  / -_)
  //  \__/_//_/\__/\___/\_,_/\__/
  //

  else if (cgf_opt.encode) {
    if (cgf_opt.tilepath<0) { printf("must specify tilepath\n"); cleanup_err(); }

    if ((cgf_opt.ifn.size()==0) && (cgf_opt.ofn.size()>0))        { cgf_opt.ifn=cgf_opt.ofn; }
    else if ((cgf_opt.ifn.size()>0) && (cgf_opt.ofn.size()==0))   { cgf_opt.ofn=cgf_opt.ifn; }
    else if ((cgf_opt.ifn.size()==0) && (cgf_opt.ofn.size()==0))  { printf("provide CGF file\n"); cleanup_err(); }

    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }
    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }
    if (ifp!=stdin) { fclose(ifp); }
    ifp = NULL;

    // Update header information if necessary
    //
    if (cgf_opt.update_header) {
      cgf->CGFVersion = cgf_opt.cgf_version_str;
      cgf->LibraryVersion= cgf_opt.cglf_version_str;
    }

    if (cgf_opt.band_ifn.size()>0) {
      if ((cgf_opt.band_ifp=fopen(cgf_opt.band_ifn.c_str(), "r"))==NULL) { perror(cgf_opt.band_ifn.c_str()); cleanup_err(); }
    } else {
      cgf_opt.band_ifp = stdin;
    }

    // Encode the band information into the CGF structure
    //

    if (cgf_opt.encode_genotype_flag) {

      k = cgf_read_genotype_band_tilepath(cgf_opt.band_ifp, cgf, cgf_opt.tilepath, cgf_opt.gtz_flag);
      if (k!=0) {
        fprintf(stderr, "error reading genotype band, got %i\n", k);
        cleanup_fail();
      }

    }
    else {

      k = cgf_read_band_tilepath(cgf_opt.band_ifp, cgf, cgf_opt.tilepath);
      if (k!=0) {
        fprintf(stderr, "error reading band, got %i\n", k);
        cleanup_fail();
      }

    }

    k = cgf_write_to_file(cgf, cgf_opt.ofn.c_str());
    if (k<0) {
      fprintf(stderr, "error writing CGF file, got %i\n", k);
      cleanup_fail();
    }

  }

  //    __                __
  //   / /  ___ ____  ___/ /
  //  / _ \/ _ `/ _ \/ _  /
  // /_.__/\_,_/_//_/\_,_/
  //

  else if (cgf_opt.show_band) {
    if (cgf_opt.tilepath<0) { printf("must specify tilepath\n"); cleanup_err(); }

    //if (!cgf_opt.ifn) { printf("provide input CGF file\n"); cleanup_err(); }
    //if ((ifp=fopen(cgf_opt.ifn, "r"))==NULL) { perror(cgf_opt.ifn); cleanup_err(); }

    if (cgf_opt.ifn.size()==0) { printf("provide input CGF file\n"); cleanup_err(); }
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }

    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }

    idx = cgf_opt.tilepath;

    if (cgf_opt.tilestep>=0) {
      n_step = cgf_opt.endtilestep - cgf_opt.tilestep + 1;
      if (n_step<0) { n_step=0; }
      if (cgf_opt.hiq_only) {
        cgf_output_band_format2(cgf, idx, stdout, cgf_opt.tilestep, n_step, (uint32_t)cgf_opt.fill_level);
      } else {

        if (cgf_opt.encode_genotype_flag) {
          cgf_output_band_format2(cgf, idx, stdout, cgf_opt.tilestep, n_step, (uint32_t)cgf_opt.fill_level, 0);
        }
        else {
          cgf_output_band_format2(cgf, idx, stdout, cgf_opt.tilestep, n_step, (uint32_t)cgf_opt.fill_level);
        }

      }
    }

    else {
      if (cgf_opt.hiq_only) {
        cgf_output_band_format(cgf, idx, stdout, 1);
      } else {
        cgf_output_band_format(cgf, idx, stdout, 0);
      }
    }

  }

  else if (cgf_opt.run_sanity) {
    if (cgf_opt.ifn.size()==0) { printf("provide input CGF file\n"); cleanup_err(); }
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }

    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }

    ret = cgf_sanity(cgf);
    if (ret==0) {
      printf("ok\n");
    } else {
      printf("sanity error %i\n", ret);
    }

  }

  //      __    __                        _      __
  //  ___/ /__ / /  __ _____ _  ___  ____(_)__  / /_
  // / _  / -_) _ \/ // / _ `/ / _ \/ __/ / _ \/ __/
  // \_,_/\__/_.__/\_,_/\_, / / .__/_/ /_/_//_/\__/
  //

  else if (cgf_opt.show_all) {
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }
    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }
    if (ifp!=stdin) { fclose(ifp); }
    ifp = NULL;

    cgf_print(cgf);
  }

  //    _      ___
  //   (_)__  / _/__
  //  / / _ \/ _/ _ \
  // /_/_//_/_/ \___/
  //

  else if (cgf_opt.info) {

    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }
    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }
    if (ifp!=stdin) { fclose(ifp); }
    ifp = NULL;


    if (cgf_opt.print_stats) {
      cgf4_print_tilepath_stats(cgf, &cgf_opt);
    }
    else {
      printf("sanity: %i\n", cgf_sanity(cgf));
    }

  }

  //                                 __
  //  _______  ___  _______  _______/ /__ ____  _______
  // / __/ _ \/ _ \/ __/ _ \/ __/ _  / _ `/ _ \/ __/ -_)
  // \__/\___/_//_/\__/\___/_/  \_,_/\_,_/_//_/\__/\__/
  //

  else if (cgf_opt.match) {

    if (cgf_opt.ifns.size()!=2) {
      printf("must provide two cgf files to match\n");
      cleanup_err();
    }

    if (!(ifp = fopen(cgf_opt.ifns[0].c_str(), "r"))) {
      perror(cgf_opt.ifns[0].c_str());
      cleanup_fail();
    }

    cgf = cgf_read_hiq(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }
    fclose(ifp);

    if (!(ifp = fopen(cgf_opt.ifns[1].c_str(), "r"))) {
      perror(cgf_opt.ifns[1].c_str());
      cleanup_fail();
    }

    cgf_b = cgf_read_hiq(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }
    fclose(ifp);

    ifp = NULL;

    // set up some default values
    //
    if (cgf_opt.tilepath<0) {
      cgf_opt.tilepath = 0;
    }

    if (cgf_opt.tilestep<0) {
      cgf_opt.tilestep = 0;
    }

    if (cgf_opt.endtilepath<0) {
      cgf_opt.endtilepath = (int)cgf->TilePathCount-1;
      if (cgf_b->TilePathCount  < cgf->TilePathCount) {
        cgf_opt.endtilepath = (int)cgf_b->TilePathCount-1;
      }
    }

    if (cgf_opt.endtilestep<0) {
      cgf_opt.endtilestep = cgf->TileStepCount[ cgf_opt.endtilepath ]-1;
      if (cgf_b->TileStepCount[ cgf_opt.endtilepath ] < cgf->TileStepCount[ cgf_opt.endtilepath ]) {
        cgf_opt.endtilestep = cgf_b->TileStepCount[ cgf_opt.endtilepath ]-1;
      }
    }


    if (cgf_opt.all_pairs==0) {

      str2tilemap(cgf->TileMap, &(cgf->TileMapCache));
      cgf->TileMapCacheInit=1;

      str2tilemap(cgf_b->TileMap, &(cgf->TileMapCache));
      cgf_b->TileMapCacheInit=1;

      if (cgf_opt.repeat==0) {
        k = cgf_hiq_concordance( &match, &tot,
            cgf, cgf_b,
            cgf_opt.tilepath, cgf_opt.tilestep,
            cgf_opt.endtilepath, cgf_opt.endtilestep,
            &cgf_opt);

        printf("match: %i, total: %i\n", match, tot);
      }
      else {

        for (i=0; i<cgf_opt.repeat; i++) {
          k = cgf_hiq_concordance( &match, &tot,
              cgf, cgf_b,
              cgf_opt.tilepath, cgf_opt.tilestep,
              cgf_opt.endtilepath, cgf_opt.endtilestep,
              &cgf_opt);
          printf("[%i] match: %i, total: %i\n", i, match, tot);
        }

      }

    }
    else {
      int cur_tilepath=0, cur_tilestep=0;
      int start_tilestep=0, end_tilestep=0, ss=0, ee=0;

      str2tilemap(cgf->TileMap, &(cgf->TileMapCache));
      cgf->TileMapCacheInit=1;

      str2tilemap(cgf_b->TileMap, &(cgf->TileMapCache));
      cgf_b->TileMapCacheInit=1;

      start_tilestep = cgf_opt.tilestep;
      for (cur_tilepath = cgf_opt.tilepath; cur_tilepath <= cgf_opt.endtilepath; cur_tilepath++) {

        ss = ((cur_tilepath == cgf_opt.tilepath) ? cgf_opt.tilestep : 0 );
        ee = ((cur_tilepath == cgf_opt.endtilepath) ? cgf_opt.endtilestep : (cgf_b->TileStepCount[ cur_tilepath ]-1) );

        for (start_tilestep = ss; start_tilestep < ee; start_tilestep++) {
          for (end_tilestep = start_tilestep; end_tilestep <= ee; end_tilestep++) {

            match=0; tot=0;
            k = cgf_hiq_concordance( &match, &tot,
                cgf, cgf_b,
                cur_tilepath, start_tilestep,
                cur_tilepath, end_tilestep,
                &cgf_opt);

            printf("[%i][%i+%i] match: %i, total: %i\n",
                cur_tilepath,
                start_tilestep, end_tilestep - start_tilestep+1,
                match, tot);

          }

        }

      }

    }

  }

  //    __               __
  //   / /  ___ ___ ____/ /__ ____
  //  / _ \/ -_) _ `/ _  / -_) __/
  // /_//_/\__/\_,_/\_,_/\__/_/
  //

  else if (cgf_opt.show_header) {
    tilemap_t tm;

    if (cgf_opt.ifn.size()==0) { printf("provide input CGF file\n"); cleanup_err(); }
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }

    //cgf = cgf_read(ifp);
    cgf = cgf_read_hiq(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }

    str2tilemap(cgf->TileMap, &(cgf->TileMapCache));
    cgf4_print_header_json(cgf, stdout);
  }

  else if (cgf_opt.print_stats) {

    if (cgf_opt.ifn.size()==0) { printf("provide input CGF file\n"); cleanup_err(); }
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }

    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }

    cgf4_print_tilepath_stats(cgf, &cgf_opt);
  }

cgf_cleanup:
  if (ifp && (ifp!=stdin)) { fclose(ifp); }
  if (ofp && (ofp!=stdout)) { fclose(ofp); }
  if (cgf_opt.band_ifp && (cgf_opt.band_ifp!=stdin)) { fclose(cgf_opt.band_ifp); }
  exit(ret);

}
