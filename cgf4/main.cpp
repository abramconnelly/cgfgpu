#include "cgf4.hpp"

#define cleanup_err() do { ret=1; show_help(); goto cgf_cleanup; } while (0);
#define cleanup_fail() do { ret=-2; goto cgf_cleanup; } while (0);
#define cleanup_ok() do { ret=0; show_help(); goto cgf_cleanup; } while (0);
#define cleanup() do { show_help(); goto cgf_cleanup; } while (0);


typedef struct cgf_opt_type {
  int show_header,
      show_band,
      encode,
      show_help,
      show_version,
      verbose,
      del,
      create_container,
      tilemap,
      show_all,
      ez_print;
  int run_test,
      info;

  //char *ifn,
  //     *ofn,
  //     *tilemap_fn;
  //char *band_ifn;

  FILE *band_ifp;

  int cgf_version_opt;
  std::vector< std::string > cgf_version_opt_ele;
  std::string cgf_version_str;
  int update_cgf_version;

  int cglf_version_opt;
  std::vector< std::string > cglf_version_opt_ele;
  std::string cglf_version_str;
  int update_cglf_version;

  std::string ifn, ofn, tilemap_fn, band_ifn;
  std::vector< std::string > ifns;

  int update_header;
  int hiq_only;

  int match;
  int run_sanity;

  int tilepath, endtilepath;
  int tilestep, endtilestep;

  uint32_t fill_level;

} cgf_opt_t;

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

  {"input",               required_argument,  NULL, 'i'},
  {"output",              required_argument,  NULL, 'o'},
  {"tilemap",             required_argument,  NULL, 't'},
  {"version-opt",         required_argument,  NULL, 'T'},
  {"library-version",     required_argument,  NULL, 'L'},
  {0,0,0,0}
};

void init_cgf_opt(cgf_opt_t *opt) {
  opt->show_header=0;
  opt->show_band=0;
  opt->encode=0;
  opt->show_help=0;
  opt->show_version=0;
  opt->verbose=0;
  opt->del=0;
  opt->create_container=0;
  opt->tilemap=0;
  opt->show_all=0;
  opt->ez_print=0;
  opt->run_test=0;
  opt->info=0;

  opt->tilepath=0;
  opt->endtilepath=-1;

  opt->tilestep=0;
  opt->endtilestep=-1;

  //opt->ifn = NULL;
  //opt->ofn = NULL;
  //opt->tilemap_fn = NULL;
  //opt->band_ifn = NULL;

  opt->band_ifp = NULL;
  opt->run_sanity = 0;

  opt->cgf_version_str = CGF_VERSION;
  opt->update_cgf_version = 0;

  opt->cglf_version_opt=0;
  opt->cglf_version_str = CGLF_VERSION;
  opt->update_cglf_version=0;
  opt->update_header=0;
  opt->hiq_only=0;

  opt->match=0;

  opt->fill_level=0xff;
}

void show_help() {
  printf("CGF Tool.  A tool used to inspect and edit Compact Genome Format (CGF) v4 files.\n");
  printf("Version: %s\n", CGF_VERSION);
  printf("\n");
  printf("usage:\n  cgf4 [-H] [-b tilepath] [-e tilepath] [-i ifn] [-o ofn] [-h] [-v] [-V] [ifn]\n");
  printf("\n");
  printf("  [-H|--header]               show header\n");
  printf("  [-C|--create-container]     create empty container\n");
  printf("  [-I|--info]                 print basic information about CGF file\n");
  printf("  [-b|--band tilepath]        output band for tilepath\n");
  printf("  [-m|--match]                run concordance on a pair of cgf files (must provide two cgf files)\n");
  printf("  [-e|--encode tilepath]      input tilepath band and add it to file, overwriting if it already exists\n");

  printf("  [-p|--tilepath tilepath]    tilepath (start)\n");
  printf("  [-P|--endtilepath tilepath] end tilepath\n");
  printf("  [-s|--tilestep tilestep]    tilestep (start)\n");
  printf("  [-S|--endtilestep tilestep] end tilestep\n");

  printf("  [-i|--input ifn]            input file (CGF)\n");
  printf("  [-o|--output ofn]           output file (CGF)\n");
  printf("  [-A|--show-all]             show all tilepaths\n");
  printf("  [-h|--help]                 show help (this screen)\n");
  printf("  [-v|--version]              show version\n");
  printf("  [-V|--verbose]              set verbose level\n");
  printf("  [-t|--tilemap tfn]          use tilemap file (instead of default)\n");
  printf("  [-Z|--ez-print]             print \"ez\" structure information\n");
  printf("  [-T|--opt-version vopt]     CGF version option.  Must be one of \"default\" or \"noc-inv\"\n");
  printf("  [-L|--library-version lopt] CGF library version option.  Overwrites default value if specified\n");
  printf("  [-U|--update-header]        Update header only\n");
  printf("  [-Y|--sanity]               Run sanity checks\n");
  printf("  [-q|--hiq]                  Only output high quality information (band output)\n");
  printf("\n");
}

void show_version() {
  printf("version %s\n", CGF_VERSION);
}

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

  init_cgf_opt(&cgf_opt);

  while ((opt = getopt_long(argc, argv, "Hb:e:i:o:Ct:T:L:U:hvVAZRIqmp:P:s:S:YF:", long_options, &option_index))!=-1) switch (opt) {
    case 0:
      fprintf(stderr, "sanity error, invalid option to parse, exiting\n");
      exit(-1);
      break;
    case 'H': cgf_opt.show_header=1; break;
    case 'Y': cgf_opt.run_sanity=1; break;
    case 'C': cgf_opt.create_container=1; break;
    case 'I': cgf_opt.info=1; break;
    case 't': cgf_opt.tilemap=1; cgf_opt.tilemap_fn=strdup(optarg); break;
    case 'b': cgf_opt.show_band=1; cgf_opt.tilepath=atoi(optarg); break;
    case 'e': cgf_opt.encode=1; cgf_opt.tilepath=atoi(optarg); break;
    case 'd': cgf_opt.del=1; cgf_opt.tilepath=atoi(optarg); break;
    //case 'i': cgf_opt.ifn=strdup(optarg); break;
    case 'i': cgf_opt.ifns.push_back(optarg); break;
    case 'o': cgf_opt.ofn=strdup(optarg); break;
    case 'h': cgf_opt.show_help=1; break;
    case 'A': cgf_opt.show_all=1; break;
    case 'v': cgf_opt.show_version=1; break;
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
    default: printf("unknown option"); show_help(); cleanup_ok(); break;
  }

  if (argc>optind) {
    if ((argc-optind)>2) { printf("Extra options specified\n"); cleanup_err(); }
    for (i=0; i<(argc-optind); i++) {
      cgf_opt.ifns.push_back(argv[optind+i]);
    }

    //if ((argc-optind)>1) { printf("Extra options specified\n"); cleanup_err(); }
    //if (cgf_opt.ifn) { printf("Input CGF already specified.\n"); cleanup_err(); }
    //cgf_opt.ifn = strdup(argv[optind]);
  }

  if (cgf_opt.ifns.size()>0) { cgf_opt.ifn = cgf_opt.ifns[0]; }

  if (cgf_opt.show_help) { show_help(); goto cgf_cleanup; }
  if (cgf_opt.show_version) { show_version(); goto cgf_cleanup; }


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
    cleanup_ok();
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

	//                  __       _             
	//   _______  ___  / /____ _(_)__  ___ ____
	//  / __/ _ \/ _ \/ __/ _ `/ / _ \/ -_) __/
	//  \__/\___/_//_/\__/\_,_/_/_//_/\__/_/   
	//                                         

  if (cgf_opt.create_container) {


    if (cgf_opt.ofn.size()==0) {
    //if (!cgf_opt.ofn) {
      //if (cgf_opt.ifn) { cgf_opt.ofn=strdup(cgf_opt.ifn); }
      if (cgf_opt.ifn.size()>0) { cgf_opt.ofn = cgf_opt.ifn; }
      else { printf("specify output CGF file\n"); cleanup_err(); }
    }

    //if (!cgf_opt.tilemap_fn) {
    if (cgf_opt.tilemap_fn.size()==0) {
      tilemap_str = DEFAULT_TILEMAP;
    } else {
      //if ((read_tilemap_from_file(tilemap_str, cgf_opt.tilemap_fn))==NULL) {
      if ((read_tilemap_from_file(tilemap_str, cgf_opt.tilemap_fn.c_str()))==NULL) {
        //perror(cgf_opt.tilemap_fn);
        perror(cgf_opt.tilemap_fn.c_str());
        cleanup_err();
      }
    }

    if (cgf_opt.ofn=="-") { ofp=stdout; }
    else if ((ofp = fopen(cgf_opt.ofn.c_str(), "w"))==NULL) { perror(cgf_opt.ofn.c_str()); cleanup_err(); }

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

    //if (!cgf_opt.ifn && cgf_opt.ofn)        { cgf_opt.ifn=strdup(cgf_opt.ofn); }
    //else if (cgf_opt.ifn && !cgf_opt.ofn)   { cgf_opt.ofn=strdup(cgf_opt.ifn); }
    //else if ((!cgf_opt.ifn) && (!cgf_opt.ofn)) { printf("provide CGF file\n"); cleanup_err(); }

    if ((cgf_opt.ifn.size()==0) && (cgf_opt.ofn.size()>0))        { cgf_opt.ifn=cgf_opt.ofn; }
    else if ((cgf_opt.ifn.size()>0) && (cgf_opt.ofn.size()==0))   { cgf_opt.ofn=cgf_opt.ifn; }
    else if ((cgf_opt.ifn.size()==0) && (cgf_opt.ofn.size()==0))  { printf("provide CGF file\n"); cleanup_err(); }

    //if ((ifp=fopen(cgf_opt.ifn, "r"))==NULL) { perror(cgf_opt.ifn); cleanup_err(); }
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

    //if (cgf_opt.band_ifn) {
    if (cgf_opt.band_ifn.size()>0) {
      //if ((cgf_opt.band_ifp=fopen(cgf_opt.band_ifn, "r"))==NULL) { perror(cgf_opt.band_ifn); cleanup_err(); }
      if ((cgf_opt.band_ifp=fopen(cgf_opt.band_ifn.c_str(), "r"))==NULL) { perror(cgf_opt.band_ifn.c_str()); cleanup_err(); }
    } else {
      cgf_opt.band_ifp = stdin;
    }

    //encode...
    //
    k = cgf_read_band_tilepath(cgf, cgf_opt.tilepath, cgf_opt.band_ifp);

    /*
    // Do some rudimentary sanity checks
    //
    k = cgf_sanity(cgf);
    if (k<0) {
      fprintf(stderr, "SANITY FAIL: %i\n", k);
      cleanup_fail();
    }
    */

    //k = cgf_write_to_file(cgf, cgf_opt.ofn);
    k = cgf_write_to_file(cgf, cgf_opt.ofn.c_str());

  }

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

    /*
    for (idx=0; idx<cgf->Path.size(); idx++) {
      if (cgf->Path[idx].TilePath == (uint64_t)cgf_opt.tilepath) { break; }
    }

    if ((uint64_t)idx==cgf->Path.size()) {
      printf("Tile Path %i not found\n", cgf_opt.tilepath);
      cleanup_ok();
    }
    */

    idx = cgf_opt.tilepath;

    if (cgf_opt.tilestep>=0) {
      n_step = cgf_opt.endtilestep - cgf_opt.tilestep + 1;
      if (n_step<0) { n_step=0; }
      if (cgf_opt.hiq_only) {
        cgf_output_band_format2(cgf, idx, stdout, cgf_opt.tilestep, n_step, cgf_opt.fill_level);
      } else {
        cgf_output_band_format2(cgf, idx, stdout, cgf_opt.tilestep, n_step, cgf_opt.fill_level);
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

  else if (cgf_opt.show_all) {
    //if ((ifp=fopen(cgf_opt.ifn, "r"))==NULL) { perror(cgf_opt.ifn); cleanup_err(); }
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

  else if (cgf_opt.info) {

    //if ((ifp=fopen(cgf_opt.ifn, "r"))==NULL) { perror(cgf_opt.ifn); cleanup_err(); }
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }
    cgf = cgf_read(ifp);
    if (!cgf) {
      //printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn);
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }
    if (ifp!=stdin) { fclose(ifp); }
    ifp = NULL;

    printf("sanity: %i\n", cgf_sanity(cgf));

  }

  else if (cgf_opt.match) {

    if (cgf_opt.ifns.size()!=2) {
      printf("must provide two cgf files to match\n");
      cleanup_err();
    }

    if (!(ifp = fopen(cgf_opt.ifns[0].c_str(), "r"))) {
      perror(cgf_opt.ifns[0].c_str());
      cleanup_fail();
    }
    //cgf = cgf_read(ifp);
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
    //cgf_b = cgf_read(ifp);
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

    printf("## [%i.%i,%i.%i]\n",
        cgf_opt.tilepath, cgf_opt.tilestep,
        cgf_opt.endtilepath, cgf_opt.endtilestep);

    //k = cgf_hiq_concordance(&match, &tot, cgf, cgf_b, 0x035e, 0, 0x35e, 34);
    k = cgf_hiq_concordance( &match, &tot,
        cgf, cgf_b,
        cgf_opt.tilepath, cgf_opt.tilestep,
        cgf_opt.endtilepath, cgf_opt.endtilestep);

    //printf("got(%i): %i / %i\n", k, match, tot);
    printf("match: %i, total: %i\n", match, tot);

  }

  else if (cgf_opt.show_header) {
    tilemap_t tm;

    if (cgf_opt.ifn.size()==0) { printf("provide input CGF file\n"); cleanup_err(); }
    if ((ifp=fopen(cgf_opt.ifn.c_str(), "r"))==NULL) { perror(cgf_opt.ifn.c_str()); cleanup_err(); }

    cgf = cgf_read(ifp);
    if (!cgf) {
      printf("CGF read error.  Is %s a valid CGFv3 file?\n", cgf_opt.ifn.c_str());
      cleanup_fail();
    }

    str2tilemap(cgf->TileMap, tm);
  }






cgf_cleanup:
  //if (cgf_opt.ifn) { free(cgf_opt.ifn); }
  //if (cgf_opt.ofn) { free(cgf_opt.ofn); }
  //if (cgf_opt.band_ifn) { free(cgf_opt.band_ifn); }
  //if (cgf_opt.tilemap_fn) { free(cgf_opt.tilemap_fn); }
  if (ifp && (ifp!=stdin)) { fclose(ifp); }
  if (ofp && (ofp!=stdout)) { fclose(ofp); }
  if (cgf_opt.band_ifp && (cgf_opt.band_ifp!=stdin)) { fclose(cgf_opt.band_ifp); }
  exit(ret);

}
