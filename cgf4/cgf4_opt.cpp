#include "cgf4.hpp"

void cgf_opt_init(cgf_opt_t *opt) {
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

  opt->encode_genotype_flag = 0;

  opt->tilepath=0;
  opt->endtilepath=-1;

  opt->tilestep=0;
  opt->endtilestep=-1;

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

  opt->all_pairs=0;
  opt->repeat = 0;

  opt->print_stats=0;
}


