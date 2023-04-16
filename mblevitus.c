/*--------------------------------------------------------------------
 *    The MB-system:  mblevitus.c  4/15/93
 *
 *    Copyright (c) 1993-2019 by
 *    David W. Caress (caress@mbari.org)
 *      Monterey Bay Aquarium Research Institute
 *      Moss Landing, CA 95039
 *    and Dale N. Chayes (dale@ldeo.columbia.edu)
 *      Lamont-Doherty Earth Observatory
 *      Palisades, NY 10964
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * MBLEVITUS generates an average water velocity profile for a
 * specified location from the Levitus temperature and salinity
 * database.
 *
 * The calculation of water sound velocity from salinity and
 * temperature observations proceeds in two steps. The first
 * step is to calculate the pressure as a function of depth
 * and latitude. We use equations from a 1989 book by Coates:
 * *
 * The second step is to calculate the water sound velocity.
 * We use the DelGrosso equation because of the results presented in
 *    Dusha, Brian D. Worcester, Peter F., Cornuelle, Bruce D.,
 *      Howe, Bruce. M. "On equations for the speed of sound
 *      in seawater", J. Acoust. Soc. Am., Vol 93, No 1,
 *      January 1993, pp 255-275.
 *
 * Author:  D. W. Caress
 * Date:  April 15, 1993
 *
 * GMTfied by	J. Luis. This version can be called by Matlab and Julia
 * Date:	May, 2020
 *
 */

#define THIS_MODULE_CLASSIC_NAME "mblevitus"
#define THIS_MODULE_MODERN_NAME "mblevitus"
#define THIS_MODULE_LIB "mbgmt"
#define THIS_MODULE_PURPOSE "Create a mean annual sound velocity profile for a specified 1 X 1 degree region."
#define THIS_MODULE_KEYS	">D}"
#define THIS_MODULE_NEEDS ""
#define THIS_MODULE_OPTIONS "-:>Vho"

// include gmt_def.h but first undefine PACKAGE variables to prevent
// warnings about name collision between GDAL's cpl_port.h and mb_config.h
#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif
#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif
#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif
#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif
#ifdef PACKAGE_URL
#undef PACKAGE_URL
#endif
#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

/* GMT header file */
#include "gmt_dev.h"

/* MBIO include files */
#ifndef _WIN32
#include "levitus.h"
#endif
#include "mb_define.h"
#include "mb_status.h"

EXTERN_MSC int GMT_mblevitus(void *API, int mode, void *args);
EXTERN_MSC char *gmt_runtime_bindir_win32 (char *result);

const double MBLEVITUS_NO_DATA = -1000000000.0;
#define NDEPTH_MAX 46
#define NLEVITUS_MAX 33

// TODO(schwehr): warning: excess elements in array initializer
double depth[48 /* NDEPTH_MAX + 2 */] =
	{0.0, 10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 125.0,
	 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0,
	 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0,
	 1750.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,
	 5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0,
	 9500.0, 10000.0, 10500.0, 11000.0, 11500.0, 12000.0};

/*--------------------------------------------------------------------*/
static struct MBLEVITUS_CTRL {
	struct mbsvplist_X {	/*  */
		int  dir;
		bool help, all4;
		char *read_file;
		char *out_file;
		double lon, lat;
	} X;
};

static void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct  MBLEVITUS_CTRL *Ctrl;

	Ctrl = gmt_M_memory (GMT, NULL, 1, struct MBLEVITUS_CTRL);
	Ctrl->X.read_file = NULL;
	Ctrl->X.dir = 1;
	return (Ctrl);
}

static void Free_Ctrl(struct GMT_CTRL *GMT, struct MBLEVITUS_CTRL *Ctrl) {	/* Deallocate control structure */
	if (!Ctrl) return;
	if (Ctrl->X.read_file) free (Ctrl->X.read_file);
	if (Ctrl->X.out_file)  free(Ctrl->X.out_file);
	gmt_M_free (GMT, Ctrl);
}

static int usage(struct GMTAPI_CTRL *API, int level) {
	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message (API, GMT_TIME_NONE, "mblevitus -Llon/lat [<levitus_file>] [-A] [-O[<outfile>]] [-H] [-V] [-z]\n");

	if (level == GMT_SYNOPSIS) return (EXIT_FAILURE);

	GMT_Message (API, GMT_TIME_NONE, "\t-L Sets the longitude and latitude of the location of the sound velocity profile.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   For backward compatibility -R is also accepted.\n");

	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t<levitus_file> Sets the location of the LevitusAnnual82.dat file [Default is builtin].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-A By default mblevitus prints (depth, velocity) but this option extends it\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   to print also (depth, velocity, temperature, salinity)\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-H print out a description of its operation and then exit immediately.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-O Write the SVP to <outfile>. Default (no -O) is to write to standard output.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   But -O alone will write to outfile = \"velocity\".\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-z Make z positive up (original MB has it positive down).\n");
	GMT_Option (API, "V,:,o");

	return (EXIT_FAILURE);
}

static int parse (struct GMT_CTRL *GMT, struct MBLEVITUS_CTRL *Ctrl, struct GMT_OPTION *options) {
	/* This parses the options provided to mbsvplist and sets parameters in Ctrl.
	 * Note Ctrl has already been initialized and non-zero default values set.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	unsigned int n_errors = 0, n_files = 0;
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

	for (opt = options; opt; opt = opt->next) {	/* Process all the options given */

		switch (opt->option) {
			case '<':	/* Input file (only one or three is accepted) */
				if (gmt_check_filearg (GMT, '<', opt->arg, GMT_IN, GMT_IS_DATASET)) {
					Ctrl->X.read_file = strdup(opt->arg);
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error: only one input file is allowed.\n");
					n_errors++;
				}
				break;

			case '>':	/* Write to named output file instead of stdout */
				if (opt->arg[0] && !Ctrl->X.out_file) Ctrl->X.out_file = strdup (opt->arg);
				break;

			/* Processes program-specific parameters */

			case 'A':	/*  Output D, T, S, V instead og only D,V */
				Ctrl->X.all4 = true;
				break;

			case 'H':	/*  Help */
				Ctrl->X.help = true;
				break;

			case 'L':
			case 'R':	/*  */
				const char *lonptr = strtok(opt->arg, "/");
				const char *latptr = strtok(NULL, "/");
				if (lonptr != NULL && latptr != NULL) {
					Ctrl->X.lon = mb_ddmmss_to_degree(lonptr);
					Ctrl->X.lat = mb_ddmmss_to_degree(latptr);
				}
				break;

			case 'O':	/*  */
				if (opt->arg)
					Ctrl->X.out_file = strdup (opt->arg);
				else
					Ctrl->X.out_file = strdup ("velocity");	/* Only to try to maintain backward compatibily */
				break;

			case 'z':	/*  Z positive up */
				Ctrl->X.dir = -1;
				break;

			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

/*--------------------------------------------------------------------*/
int GMT_mblevitus(void *V_API, int mode, void *args) {
	static const char program_name[] = "MBLEVITUS";
	static const char help_message[] =
		"MBLEVITUS generates an average water velocity profile for a\n"
		"specified location from the Levitus temperature and salinity database.";
	static const char usage_message[] = "mblevitus -Llon/lat [-A] [-I<levitus_file>] [-O[<outfile>]] [-V] [-H]";

	int verbose = 0, error = 0;
	unsigned int ix, iy;

	double longitude = 0.0;
	double latitude = 0.0;

	char *ofile, *read_file;
	char header[GMT_BUFSIZ * 2] = {""}; 
	char buff[GMT_LEN128] = {""};
	char dir[PATH_MAX] = {""};
	char levitus_file[PATH_MAX] = {""};

	/* process argument list */
	bool help = false;

	uint64_t dim[4] = {1, 1, 0, 2};

	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL; /* General GMT interal parameters */
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr(V_API); /* Cast from void to GMTAPI_CTRL pointer */
	struct GMT_RECORD *Out = NULL;
	struct GMT_DATASET *D = NULL;
	struct GMT_DATASEGMENT *S = NULL;
	struct GMT_DATASEGMENT_HIDDEN *SH = NULL;
	struct MBLEVITUS_CTRL *Ctrl = NULL;

	/*----------------------- Standard module initialization and parsing ----------------------*/
	if (API == NULL) Return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) Return (usage(API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) Return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE) bailout(usage (API, GMT_USAGE));	/* Return the usage message */
	if (options->option == GMT_OPT_SYNOPSIS) bailout(usage (API, GMT_SYNOPSIS));	/* Return the synopsis */

	/* Parse the command-line arguments */

#if GMT_MAJOR_VERSION >= 6
	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, NULL, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
#else
	GMT = gmt_begin_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, &GMT_cpy); /* Save current state */
#endif
	if (GMT_Parse_Common(API, THIS_MODULE_OPTIONS, options)) Return (API->error);
	Ctrl = (struct MBLEVITUS_CTRL *)New_Ctrl(GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse(GMT, Ctrl, options))) {
		GMT_Report (API, GMT_MSG_NORMAL, "usage: %s\n", usage_message);
		GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
		Return(MB_ERROR_BAD_USAGE);
	}
	/*---------------------------- This is the mbsvplist main code ----------------------------*/

	read_file = Ctrl->X.read_file;
	latitude = Ctrl->X.lat;
	longitude = Ctrl->X.lon;

	/* Try to combine the MB & GMT verbose systems and make it also compatible with GMT < 6.1 */
	if (GMT->current.setting.verbose > 3)
		verbose = 1;
	if (GMT->current.setting.verbose == GMT_MSG_DEBUG)
		verbose = 2;

	if (!Ctrl->X.out_file) 				/* Not external and no output name (-O option) provided */
		ofile = NULL;
	else if (Ctrl->X.out_file[0] == '@') 	/* External and no output name (-O option) provided */
		ofile = Ctrl->X.out_file;	
	else								/* Got one output explicitly name */
		ofile = Ctrl->X.out_file;

	if (verbose == 1 || Ctrl->X.help) {
		GMT_Message (API, GMT_TIME_NONE, "Program %s\n", program_name);
		GMT_Message (API, GMT_TIME_NONE, "MB-system Version %s\n", MB_VERSION);
	}

	if (Ctrl->X.help) {
		GMT_Message (API, GMT_TIME_NONE, "\n%s\n", help_message);
		GMT_Message (API, GMT_TIME_NONE, "usage: %s\n", usage_message);
		Return(MB_ERROR_NO_ERROR);
	}

	if (Ctrl->X.read_file && access (Ctrl->X.read_file, R_OK) == 0) 	/* Input levitus file can be read */
		strncpy(levitus_file, Ctrl->X.read_file, PATH_MAX-1);
#ifndef _WIN32
	else if (levitusfile)		/* On Unix this is set at compile time */
		strncpy(levitus_file, levitus_file, PATH_MAX-1);
#endif
	else {						/* We are probably on windows or using the GMT+MB package */
		snprintf (dir, PATH_MAX, "%s/%s", GMT->session.SHAREDIR, "mbsystem/LevitusAnnual82.dat");
		if (access (dir, R_OK) == 0) 	/* File can be read */
			strncpy(levitus_file, dir, PATH_MAX-1);
		else if (!API->external) {
#ifdef _WIN32
			/* Find the path to the bin directory and from it, the location of the Levitus file */
			char dir[PATH_MAX + 1];
			gmt_runtime_bindir_win32(dir);
			char *pch = strrchr(dir, '\\');
			pch[0] = '\0';
			strcat(dir, "\\share\\mbsystem\\LevitusAnnual82.dat");
			strncpy(levitus_file, dir, PATH_MAX-1);
#endif
		}
	}

	if (!levitus_file) {
		GMT_Report (API, GMT_MSG_NORMAL, "Could not find the Levitus database file LevitusAnnual82.dat\n");
		Return(MB_ERROR_BAD_USAGE);
	}

	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  Program <%s>\n", program_name);
		fprintf(stderr, "dbg2  MB-system Version %s\n", MB_VERSION);
		fprintf(stderr, "dbg2  Control Parameters:\n");
		fprintf(stderr, "dbg2       verbose:          %d\n", verbose);
		fprintf(stderr, "dbg2       help:             %d\n", help);
		fprintf(stderr, "dbg2       levitusfile:      %s\n", levitus_file);
		fprintf(stderr, "dbg2       ofile:            %s\n", ofile);
		fprintf(stderr, "dbg2       longitude:        %f\n", longitude);
		fprintf(stderr, "dbg2       latitude:         %f\n", latitude);
	}

	FILE *ifp = fopen(levitus_file, "rb");
	if (ifp == NULL) {
		GMT_Report (API, GMT_MSG_NORMAL, "Unable to Open Levitus database file <%s> for reading\n", levitus_file);
		GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
		Return(MB_ERROR_OPEN_FAIL);
	}

	if (longitude < -360.0 || longitude > 360.0 || latitude < -90.0 || latitude > 90.0) {
		GMT_Report (API, GMT_MSG_NORMAL, "Invalid location specified:  longitude: %f  latitude: %f\n", longitude, latitude);
		GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
		Return(MB_ERROR_BAD_PARAMETER);
	}

	/* get the longitude and latitude indices */
	int ilon;
	if (longitude < 0.0)
		ilon = (int)(longitude + 360.0);
	else if (longitude >= 360.0)
		ilon = (int)(longitude - 360.0);
	else
		ilon = (int)longitude;

	const double lon_actual = ilon + 0.5;
	const int ilat = (int)(latitude + 90.0);
	const double lat_actual = ilat - 89.5;
	GMT_Report (API, GMT_MSG_VERBOSE, "Location for mean annual water velocity profile:\n");
	GMT_Report (API, GMT_MSG_VERBOSE, "  Requested:  %6.4f longitude   %6.4f latitude\n", longitude, latitude);
	GMT_Report (API, GMT_MSG_VERBOSE, "  Used:       %6.1f longitude   %6.1f latitude\n", lon_actual, lat_actual);

	int status = MB_SUCCESS;
	error = MB_ERROR_NO_ERROR;

	/* read the temperature */
	const int record_size = sizeof(float) * NLEVITUS_MAX * 180;
	long location = ilon * record_size;
	/* int status = */ fseek(ifp, location, 0);
	float temperature[NLEVITUS_MAX][180];
	if (fread(&temperature[0][0], 1, record_size, ifp) != record_size) {
		status = MB_FAILURE;
		error = MB_ERROR_EOF;
		GMT_Report (API, GMT_MSG_NORMAL, "ERROR: EOF reading temperature\n");
	}

	/* read the salinity */
	location = location + 360 * record_size;
	/* status = */ fseek(ifp, location, 0);
	float salinity[NLEVITUS_MAX][180];
	if (fread(&salinity[0][0], 1, record_size, ifp) != record_size) {
		status = MB_FAILURE;
		error = MB_ERROR_EOF;
		GMT_Report (API, GMT_MSG_NORMAL, "ERROR: EOF reading salinity\n");
	}

	fclose(ifp);		/* Close the LevitusAnnual82.dat data file */

#ifdef BYTESWAPPED
	for (int i = 0; i < NLEVITUS_MAX; i++) {
		mb_swap_float(&temperature[i][ilat]);
		mb_swap_float(&salinity[i][ilat]);
	}
#endif

	/* calculate velocity from temperature and salinity */
	int nvelocity = 0;
	int nvelocity_tot = 0;
	int last_good = -1;
	double velocity[NDEPTH_MAX];
	for (int i = 0; i < NDEPTH_MAX; i++) {
		if (i < NLEVITUS_MAX)
			if (salinity[i][ilat] > MBLEVITUS_NO_DATA) {
				last_good = i;
				nvelocity++;
			}
		if (last_good >= 0) {		/* set counter */
			nvelocity_tot++;

			/* get pressure for a given depth as a function of latitude */
			double pressure = 1.0052405 * depth[i] * (1.0 + 0.00528 * sin(DTR * latitude) * sin(DTR * latitude)) +
							  0.00000236 * depth[i] * depth[i];
			/* calculate water sound speed using DelGrosso equations */
			/* convert decibar to kg/cm**2 */
			pressure = pressure * 0.1019716;
			const double c0 = 1402.392;
			const double dltact = temperature[last_good][ilat] * (5.01109398873 + temperature[last_good][ilat] *
			                      (-0.0550946843172 + temperature[last_good][ilat] * 0.000221535969240));
			const double dltacs = salinity[last_good][ilat] * (1.32952290781 + salinity[last_good][ilat] * 0.000128955756844);
			const double dltacp = pressure * (0.156059257041 + pressure * (2.4499868841E-5 + pressure * -0.883392332513E-8));
			const double dcstp = temperature[last_good][ilat] * (-0.0127562783426 * salinity[last_good][ilat] +
			                     pressure * (0.00635191613389 + pressure * (0.265484716608E-7 * temperature[last_good][ilat] -
								 0.159349479045E-5 + 0.522116437235E-9 * pressure) - 0.438031096213E-6 *
								 temperature[last_good][ilat] * temperature[last_good][ilat])) + salinity[last_good][ilat] *
								 (-0.161674495909E-8 * salinity[last_good][ilat] * pressure * pressure +
								 temperature[last_good][ilat] * (0.968403156410E-4 * temperature[last_good][ilat] + pressure *
								 (0.485639620015E-5 * salinity[last_good][ilat] - 0.340597039004E-3)));
			velocity[i] = rint((c0 + dltact + dltacs + dltacp + dcstp) * 1000) / 1000;	/* 3 decimal places is good enough */
		}
		else
			velocity[i] = salinity[i][ilat];
	}

	/* check for existence of water velocity profile */
	if (nvelocity < 1) {
		GMT_Report (API, GMT_MSG_NORMAL, "No water velocity profile available for specified location.\n");
		GMT_Report (API, GMT_MSG_NORMAL, "\tThis place is probably subaerial!. No output file created.\n");
		GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
		Return(MB_ERROR_BAD_PARAMETER);
	}

	if (Ctrl->X.all4) dim[3] = 4;		/* Want all 4: Depth, Velocity, Temperature, Salinity  */
	if ((D = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_LINE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
		Return (API->error);			/* An empty dataset */

	snprintf(header, GMT_BUFSIZ * 2, "Sound velocity profile created by %s MB-System Version %s\n", program_name, MB_VERSION);
#ifndef _WIN32							/* This is so unix (and not really important info) */
	{
		const time_t right_now = time((time_t *)0);
		char date[32] = {""};
		strcpy(date, ctime(&right_now));
		date[strlen(date) - 1] = '\0';
		const char *user_ptr = getenv("USER");
		if (user_ptr == NULL)
			user_ptr = getenv("LOGNAME");
		char user[128] = {""};
		(user_ptr != NULL) ? strcpy(user, user_ptr) : strcpy(user, "unknown");
		char host[128] = {""};
		gethostname(host, 128);
		snprintf(buff, GMT_LEN128, "# Run by user <%s> on cpu <%s> at <%s>\n#", user, host, date); strcat(header, buff);
	}
#endif

	snprintf(buff, GMT_LEN128, "# Water velocity profile derived from Levitus temperature and salinity database.");
		strcat(header, buff);
	snprintf(buff, GMT_LEN128, "\n# This profile represents the annual average water velocity structure for a 1 X 1");
		strcat(header, buff);
	snprintf(buff, GMT_LEN128, "\n# area centered at %g longitude and %g latitude.", lon_actual, lat_actual);
		strcat(header, buff);
	snprintf(buff, GMT_LEN128, "\n# This water velocity profile is in the form of discrete");
		strcat(header,buff);
	if (Ctrl->X.all4)
		{snprintf(buff, GMT_LEN128, "\n# (depth, velocity, temperatue, salinity) points where"); strcat(header,buff);}
	else
		{snprintf(buff, GMT_LEN128, "\n# (depth, velocity) points where");	strcat(header,buff);}
	snprintf(buff, GMT_LEN128, "\n# the depth is in meters and the velocity in meters/second. The first %d", nvelocity);
		strcat(header,buff);
	snprintf(buff, GMT_LEN128, "\n# velocity values are defined using the salinity and temperature values available in the");
		strcat(header,buff);
	snprintf(buff, GMT_LEN128, "\n# Levitus database. The remaining %d velocity values are", nvelocity_tot - nvelocity);
		strcat(header,buff);
	snprintf(buff, GMT_LEN128, "\n# calculated using the deepest temperature and salinity value available.");
		strcat(header,buff);

	S = GMT_Alloc_Segment (API, GMT_NO_STRINGS, nvelocity_tot, dim[3], NULL, NULL);
	S->header = strdup(header);

	if (Ctrl->X.dir == -1)
		for (int i = 0; i < NDEPTH_MAX; i++) depth[i] = -depth[i];

	ix = (GMT->current.setting.io_lonlat_toggle[GMT_IN]);	iy = 1 - ix;
	gmt_M_memcpy (S->data[ix], depth, nvelocity_tot, double);
	gmt_M_memcpy (S->data[iy], velocity, nvelocity_tot, double);

	if (Ctrl->X.all4) {
		int i;
		for (i = 0; i < nvelocity; i++) {
			S->data[2][i] = temperature[i][ilat];
			S->data[3][i] = salinity[i][ilat];
		}
		for (i = nvelocity; i < nvelocity_tot; i++) {
			S->data[2][i] = S->data[3][i] = 0;
		}
	}

	S->n_rows = nvelocity_tot;
	D->table[0]->segment[0] = S;
	D->table[0]->n_records += nvelocity_tot;	D->n_records += nvelocity_tot;

	GMT_Report (API, GMT_MSG_VERBOSE, "Values defined directly by Levitus database:      %2d\n", nvelocity);
	GMT_Report (API, GMT_MSG_VERBOSE, "Values assuming deepest salinity and temperature: %2d\n", nvelocity_tot - nvelocity);
	GMT_Report (API, GMT_MSG_VERBOSE, "Velocity points written:                          %2d\n", nvelocity_tot);

	/* Since we only have a single segment make the segment marker a # to not distinguish from a comment */
	char sm[2] = {""};
	sm[0] = GMT->current.setting.io_seg_marker[GMT_OUT];
	GMT->current.setting.io_seg_marker[GMT_OUT] = '#';
	if (GMT->current.setting.io_lonlat_toggle[GMT_IN])		/* Also toggles in output */
		GMT->current.setting.io_lonlat_toggle[GMT_OUT] = false;

	if (GMT_Write_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_LINE, GMT_WRITE_SET, NULL, ofile, D) != GMT_NOERROR)
		Return (API->error);

	GMT->current.setting.io_seg_marker[GMT_OUT] = sm[0];

	Return(error);
}
