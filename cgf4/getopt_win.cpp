/*
POSIX getopt for Windows

AT&T Public License

Code given out at the 1985 UNIFORUM conference in Dallas.  
*/

#ifndef __GNUC__

#include "getopt_win.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define	EMSG	""
#define	BADCH	(int)'?'
#define	BADARG	(int)':'

#define NULL	0
#define EOF	(-1)
#define ERR(s, c)	if(opterr){\
	char errbuf[2];\
	errbuf[0] = c; errbuf[1] = '\n';\
	fputs(argv[0], stderr);\
	fputs(s, stderr);\
	fputc(c, stderr);}
	//(void) write(2, argv[0], (unsigned)strlen(argv[0]));\
	//(void) write(2, s, (unsigned)strlen(s));\
	//(void) write(2, errbuf, 2);}

int	opterr = 1;
int	optind = 1;
int	optopt;
int	optreset;
char	*optarg;

#define optional_argument		0
#define required_argument		1

int getopt(int argc, char** argv, char* opts)
{
	static int sp = 1;
	register int c;
	register char *cp;

	if(sp == 1)
		if(optind >= argc ||
		   argv[optind][0] != '-' || argv[optind][1] == '\0')
			return(EOF);
		else if(strcmp(argv[optind], "--") == NULL) {
			optind++;
			return(EOF);
		}
	optopt = c = argv[optind][sp];
	if(c == ':' || (cp=strchr(opts, c)) == NULL) {
		ERR(": illegal option -- ", c);
		if(argv[optind][++sp] == '\0') {
			optind++;
			sp = 1;
		}
		return('?');
	}
	if(*++cp == ':') {
		if(argv[optind][sp+1] != '\0')
			optarg = &argv[optind++][sp+1];
		else if(++optind >= argc) {
			ERR(": option requires an argument -- ", c);
			sp = 1;
			return('?');
		} else
			optarg = argv[optind++];
		sp = 1;
	} else {
		if(argv[optind][++sp] == '\0') {
			sp = 1;
			optind++;
		}
		optarg = NULL;
	}
	return(c);
}

int getopt_internal(int nargc, char** nargv, char* ostr)
{
	static char *place = EMSG;		/* option letter processing */
	char *oli;				/* option letter list index */

	assert(nargv != NULL);
	assert(ostr != NULL);

	if (optreset || !*place) {		/* update scanning pointer */
		optreset = 0;
		if (optind >= nargc || *(place = nargv[optind]) != '-') {
			place = EMSG;
			return (-1);
		}
		if (place[1] && *++place == '-') {	/* found "--" */
											/* ++optind; */
			place = EMSG;
			return (-2);
		}
	}					/* option letter okay? */
	if ((optopt = (int)*place++) == (int)':' ||
		!(oli = strchr(ostr, optopt))) {
		/*
		* if the user didn't specify '-' as an option,
		* assume it means -1.
		*/
		if (optopt == (int)'-')
			return (-1);
		if (!*place)
			++optind;
		if (opterr && *ostr != ':')
			(void)fprintf(stderr, "%s: illegal option -- %c\n", nargv[0], optopt);
		return (BADCH);
	}
	if (*++oli != ':') {			/* don't need argument */
		optarg = NULL;
		if (!*place)
			++optind;
	}
	else {				/* need an argument */
		if (*place)			/* no white space */
			optarg = place;
		else if (nargc <= ++optind) {	/* no arg */
			place = EMSG;
			if ((opterr) && (*ostr != ':'))
				(void)fprintf(stderr,
					"%s: option requires an argument -- %c\n",
					nargv[0], optopt);
			return (BADARG);
		}
		else				/* white space */
			optarg = nargv[optind];
		place = EMSG;
		++optind;
	}
	return (optopt);			/* dump back option letter */
}

/*
* getopt_long --
*	Parse argc/argv argument vector.
*/
int getopt_long(int nargc, char** nargv, char* options, struct option* long_options, int* index)

{
	int retval;

	assert(nargv != NULL);
	assert(options != NULL);
	assert(long_options != NULL);
	/* index may be NULL */

	if ((retval = getopt_internal(nargc, nargv, options)) == -2) {
		char *current_argv = nargv[optind++] + 2, *has_equal;
		int i, current_argv_len, match = -1;

		if (*current_argv == '\0') {
			return(-1);
		}
		if ((has_equal = strchr(current_argv, '=')) != NULL) {
			current_argv_len = has_equal - current_argv;
			has_equal++;
		}
		else
			current_argv_len = strlen(current_argv);

		for (i = 0; long_options[i].name; i++) {
			if (strncmp(current_argv, long_options[i].name, current_argv_len))
				continue;

			if (strlen(long_options[i].name) == (unsigned)current_argv_len) {
				match = i;
				break;
			}
			if (match == -1)
				match = i;
		}
		if (match != -1) {
			if (long_options[match].has_arg == required_argument ||
				long_options[match].has_arg == optional_argument) {
				if (has_equal)
					optarg = has_equal;
				else
					optarg = nargv[optind++];
			}
			if ((long_options[match].has_arg == required_argument)
				&& (optarg == NULL)) {
				/*
				* Missing argument, leading :
				* indicates no error should be generated
				*/
				if ((opterr) && (*options != ':'))
					(void)fprintf(stderr,
						"%s: option requires an argument -- %s\n",
						nargv[0], current_argv);
				return (BADARG);
			}
		}
		else { /* No matching argument */
			if ((opterr) && (*options != ':'))
				(void)fprintf(stderr,
					"%s: illegal option -- %s\n", nargv[0], current_argv);
			return (BADCH);
		}
		if (long_options[match].flag) {
			*long_options[match].flag = long_options[match].val;
			retval = 0;
		}
		else
			retval = long_options[match].val;
		if (index)
			*index = match;
	}
	return(retval);
}

#endif  /* __GNUC__ */