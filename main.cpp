
#include "app_perf.h"
#include "cgf4.hpp"
#include "cgfx.hpp"


int main(int argc, char **argv) 
{
  PERF_INIT(64, true, true, true, 0, "");

  int mode = 0;		// default to cgf4

  for (int n = 0; n < argc; n++) {
	  if (strcmp(argv[n], "-x") == 0) { mode = 1; break; }
  }
  switch (mode) {
  case 0: default:	cgf4_main(argc, argv);	break;		// Run cgf4 commands  
  case 1:			cgfx_main(argc, argv);  break;		// Run cgfx commands
  };

}