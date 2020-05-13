/*--------------------------------------------------------------------
 *    The MB-system:	mbsvplist.c	1/3/2001
 *
 *    Copyright (c) 2001-2019 by
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
 * This program, mbsvplist, lists all water sound velocity
 * profiles (SVPs) within swath data files. Swath bathymetry is
 * calculated from raw angles and travel times by raytracing
 * through a model of the speed of sound in water. Many swath
 * data formats allow SVPs to be embedded in the data, and
 * often the SVPs used to calculate the data will be included.
 * By default, all unique SVPs encountered are listed to
 * stdout. The SVPs may instead be written to individual files
 * with names FILE_XXX.svp, where FILE is the swath data
 * filename and XXX is the SVP count within the file.  The -D
 * option causes duplicate SVPs to be output. The -P option
 * implies -O, and also causes the parameter file to be modified
 * so that the first svp output for each file becomes the
 * svp used for recalculating bathymetry for that swath file.
 *
 * Author:	D. W. Caress
 * Date:	January 3,  2001
 */

#define THIS_MODULE_CLASSIC_NAME "mbsvplist"
#define THIS_MODULE_MODERN_NAME "mbsvplist"
#define THIS_MODULE_LIB "mbgmt"
#define THIS_MODULE_PURPOSE "Lists all water sound velocity profiles (SVPs) within swath data files."
#define THIS_MODULE_KEYS	"<D{,>D}"
#define THIS_MODULE_NEEDS ""
#define THIS_MODULE_OPTIONS "-:>RVhi"

/* GMT header file */
#include "gmt_dev.h"

EXTERN_MSC int GMT_mbsvplist(void *API, int mode, void *args);


/* MBIO include files */
#include "mb_status.h"
#include "mb_format.h"
#include "mb_define.h"
#include "mb_process.h"

#define MBSVPLIST_SVP_NUM_ALLOC 24;
typedef enum {
	MBSVPLIST_PRINTMODE_CHANGE = 0,
	MBSVPLIST_PRINTMODE_UNIQUE = 1,
	MBSVPLIST_PRINTMODE_ALL = 2,
} printmode_t;

struct mbsvplist_svp_struct {
	bool time_set;		  /* time stamp known */
	bool position_set;	  /* position known */
	bool repeat_in_file;  /* repeats a previous svp in the same file */
	bool match_last;	  /* repeats the last svp in the same file or the previous file */
	bool depthzero_reset; /* uppermost SVP value set to zero depth */
	double time_d;
	double longitude;
	double latitude;
	double depthzero;
	int n;
	double depth[MB_SVP_MAX];
	double velocity[MB_SVP_MAX];
};

//#define PRINT(C,condition,Out,buf,fp,...) ((condition) ? fprintf(fp, __VA_ARGS__) : sprintf(buf, __VA_ARGS__), GMT_Put_Record(C, GMT_WRITE_DATA, Out))

static struct MBSVPLIST_CTRL {
	struct mbsvplist_F {	/*  */
		bool active;
	} F;
	struct mbsvplist_X {	/*  */
		bool help;
		bool ssv_output;
		bool svp_file_output;
		bool svp_setprocess;
		bool output_as_table;
		bool svp_force_zero;
		bool output_counts;
		printmode_t	svp_printmode;
		int format;
		int min_num_pairs;
		int verbose;
		char *read_file;
		char *out_file;
	} X;
};

static void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct  MBSVPLIST_CTRL *Ctrl;

	Ctrl = gmt_M_memory (GMT, NULL, 1, struct MBSVPLIST_CTRL);
	Ctrl->X.read_file = NULL;
	return (Ctrl);
}

static void Free_Ctrl(struct GMT_CTRL *GMT, struct MBSVPLIST_CTRL *Ctrl) {	/* Deallocate control structure */
	if (!Ctrl) return;
	if (Ctrl->X.read_file) free (Ctrl->X.read_file);
	if (Ctrl->X.out_file)  free(Ctrl->X.out_file);
	gmt_M_free (GMT, Ctrl);
}

static int usage(struct GMTAPI_CTRL *API, int level) {
	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message (API, GMT_TIME_NONE, "usage: mbinfo -I<inputfile>\n");

	return (EXIT_FAILURE);
}

static int parse (struct GMT_CTRL *GMT, struct MBSVPLIST_CTRL *Ctrl, struct GMT_OPTION *options) {
	/* This parses the options provided to mbsvplist and sets parameters in Ctrl.
	 * Note Ctrl has already been initialized and non-zero default values set.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	unsigned int n_errors = 0, n_files = 0, n;
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
				if (opt->arg[0]) Ctrl->X.out_file = strdup (opt->arg);
				break;

			/* Processes program-specific parameters */

			case 'C':	/*  */
				Ctrl->X.output_counts = true;
				Ctrl->X.ssv_output = false;
				break;

			case 'D':	/*  */
				Ctrl->X.svp_printmode = MBSVPLIST_PRINTMODE_ALL;
				break;

			case 'F':	/* format */
				n = sscanf(opt->arg, "%d", &(Ctrl->X.format));
				if (n == 1)
					Ctrl->F.active = true;
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -F option: \n");
					n_errors++;
				}
				break;

			case 'H':	/*  Help */
				Ctrl->X.help = true;
				break;

			case 'I':	/* -I<inputfile> */
				if (!gmt_access (GMT, opt->arg, R_OK)) {	/* Got a file */
					Ctrl->X.read_file = strdup(opt->arg);
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -I: Requires a valid file\n");
					n_errors++;
				}
				break;

			case 'M':	/*  */
				Ctrl->X.svp_printmode = (printmode_t)atoi(opt->arg);
				break;

			case 'N':	/*  */
				Ctrl->X.min_num_pairs = atoi(opt->arg);
				break;

			case 'O':	/*  */
				Ctrl->X.svp_file_output = true;
				Ctrl->X.ssv_output = false;
				break;

			case 'P':	/*  */
				Ctrl->X.svp_file_output = true;
				Ctrl->X.svp_setprocess = true;
				Ctrl->X.ssv_output = false;
				break;

			case 'S':	/*  */
				Ctrl->X.ssv_output = true;
				Ctrl->X.svp_file_output = true;
				Ctrl->X.svp_setprocess = true;
				break;

			case 'T':	/*  */
				Ctrl->X.output_as_table = true;
				Ctrl->X.ssv_output = false;
				break;

			case 'Z':	/*  */
				Ctrl->X.svp_force_zero = true;
				break;

			case 'V':	/*  Verbose */
				Ctrl->X.verbose++;
				break;

			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	if (!Ctrl->X.read_file)			/* I don't think this is a good default but won't break compat. JL */
		Ctrl->X.read_file = strdup("datalist.mb-1");

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

/*--------------------------------------------------------------------*/
int GMT_mbsvplist(void *V_API, int mode, void *args) {

	static const char program_name[] = "mbsvplist";
	static const char help_message[] =
		"mbsvplist lists all water sound velocity profiles (SVPs) within swath data files.\n"
		"Swath bathymetry is calculated from raw angles and travel times by raytracing\n"
		"through a model of the speed of sound in water. Many swath data formats allow SVPs\n"
		"to be embedded in the data, and often the SVPs used to calculate the data will be included.\n"
		"By default, all unique SVPs encountered are listed to stdout. The SVPs may instead be written\n"
		"to individual files with names FILE_XXX.svp, where FILE is the swath data filename and XXX\n"
		"is the SVP count within the file. The -D option causes duplicate SVPs to be output.\n"
		"The -T option will output a CSV table of svp#, time, longitude, latitude and number of points for SVPs.\n"
		"When the -Nmin_num_pairs option is used, only svps that have at least min_num_pairs svp values will\n"
		"be output.(This is particularly useful for .xse data where the svp is entered as a single values svp.)\n";
	static const char usage_message[] = "mbsvplist [-C -D -Fformat -H -Ifile -Mmode -O -Nmin_num_pairs -P -T -V -Z]";

	int verbose = 0;
	int format;
	int pings = 1;
	int lonflip;
	int btime_i[7];
	int etime_i[7];
	int error = MB_ERROR_NO_ERROR;
	double bounds[4];
	double speedmin;
	double timegap;
	int status = mb_defaults(verbose, &format, &pings, &lonflip, bounds, btime_i, etime_i, &speedmin, &timegap);
	bounds[0] = -360.0;
	bounds[1] = 360.0;
	bounds[2] = -90.0;
	bounds[3] = 90.0;

	printmode_t svp_printmode = MBSVPLIST_PRINTMODE_CHANGE;
	char *read_file = NULL;
	bool output_counts = false;
	bool ssv_output = false;
	bool svp_file_output = false;
	bool svp_setprocess = false;
	bool ssv_bounds_set = false;
	bool output_as_table = false;
	bool svp_force_zero = false;
	int min_num_pairs = 0;
	double ssv_bounds[4] = {-360.0, 360.0, -90.0, 90.0};

	/* determine whether to read one file or a list of files */
	bool read_datalist;
	bool read_data;
	void *datalist;
	char file[MB_PATH_MAXLINE] = {""};
	char dfile[MB_PATH_MAXLINE] = {""};
	double file_weight;
	/* MBIO read control parameters */
	double btime_d;
	double etime_d;
	int beams_bath;
	int beams_amp;
	int pixels_ss;

	/* MBIO read values */
	void *mbio_ptr = NULL;
	void *store_ptr;
	int kind;
	int time_i[7];
	double time_d;
	double navlon;
	double navlat;
	double speed;
	double heading;
	double distance;
	double altitude;
	double sonardepth;
	char   *beamflag = NULL;
	double *bath = NULL;
	double *bathacrosstrack = NULL;
	double *bathalongtrack = NULL;
	double *amp = NULL;
	double *ss = NULL;
	double *ssacrosstrack = NULL;
	double *ssalongtrack = NULL;
	char comment[MB_COMMENT_MAXLINE];

	/* save time stamp and position of last survey data */
	double last_time_d = 0.0;
	double last_navlon = 0.0;
	double last_navlat = 0.0;

	/* data record source types */
	int platform_source;
	int nav_source;
	int sensordepth_source;
	int heading_source;
	int attitude_source;
	int svp_source;

	/* output mode settings */

	/* SVP values */
	struct mbsvplist_svp_struct svp;
	struct mbsvplist_svp_struct svp_last;
	int svp_save_alloc = 0;
	struct mbsvplist_svp_struct *svp_save = NULL;
	char svp_file[MB_PATH_MAXLINE] = {""};
	int svp_read_tot = 0;
	int svp_written_tot = 0;
	int svp_repeat_in_file;
	int out_cnt = 0;
	int svp_time_i[7];

	/* ttimes values */
	int nbeams;
	double *ttimes = NULL;
	double *angles = NULL;
	double *angles_forward = NULL;
	double *angles_null = NULL;
	double *heave = NULL;
	double *alongtrack_offset = NULL;
	double ssv;

	bool svp_match_last = false;
	int svp_unique_tot = 0;

	bool svp_loaded = false;
	int svp_save_count = 0;
	int svp_read = 0;
	int svp_written = 0;
	int svp_unique = 0;

	char header[GMT_BUFSIZ * 5] = {""}; 
	char buf512[GMT_LEN512] = {""};
	unsigned int ix, iy, tbl = 0;
	size_t n_save = 0, n_alloc = 0, n_seg_alloc = 0;
	uint64_t n_seg = 0;
	uint64_t dim[4] = {1, 0, 0, 2};
	double ssv_out[2];

	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL; /* General GMT interal parameters */
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr(V_API); /* Cast from void to GMTAPI_CTRL pointer */
	struct GMT_RECORD *Out = NULL;
	struct GMT_DATASET *D = NULL;
	struct GMT_DATASEGMENT *S = NULL;
	struct MBSVPLIST_CTRL *Ctrl = NULL;
	
	svp_last.n = 0;

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage(API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE) bailout(usage (API, GMT_USAGE));	/* Return the usage message */
	if (options->option == GMT_OPT_SYNOPSIS) bailout(usage (API, GMT_SYNOPSIS));	/* Return the synopsis */

	/* Parse the command-line arguments */

	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, NULL, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
	if (GMT_Parse_Common(API, THIS_MODULE_OPTIONS, options)) Return (API->error);
	Ctrl = (struct MBSVPLIST_CTRL *)New_Ctrl(GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse(GMT, Ctrl, options))) {
		GMT_Report (API, GMT_MSG_NORMAL, "usage: %s\n", usage_message);
		GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
		Return(MB_ERROR_BAD_USAGE);
	}
	/*---------------------------- This is the mbsvplist main code ----------------------------*/

	if (verbose == 1 || Ctrl->X.help) {
		GMT_Message (API, GMT_TIME_NONE, "Program %s\n", program_name);
		GMT_Message (API, GMT_TIME_NONE, "MB-system Version %s\n", MB_VERSION);
	}

	if (Ctrl->X.help) {
		GMT_Message (API, GMT_TIME_NONE, "\n%s\n", help_message);
		GMT_Message (API, GMT_TIME_NONE, "usage: %s\n", usage_message);
		Return(MB_ERROR_NO_ERROR);
	}

	read_file = Ctrl->X.read_file;
	output_counts = Ctrl->X.output_counts;
	ssv_output = Ctrl->X.ssv_output;
	svp_file_output = Ctrl->X.svp_file_output;
	svp_setprocess = Ctrl->X.svp_setprocess;
	output_as_table = Ctrl->X.output_as_table;
	svp_force_zero = Ctrl->X.svp_force_zero;
	svp_printmode = Ctrl->X.svp_printmode;
	format = Ctrl->X.format;
	min_num_pairs = Ctrl->X.min_num_pairs;
	verbose = Ctrl->X.verbose;

	if (GMT->common.R.active[RSET]) {
		ssv_bounds[0] = GMT->common.R.wesn[XLO];	ssv_bounds[1] = GMT->common.R.wesn[XHI];
		ssv_bounds[2] = GMT->common.R.wesn[YLO];	ssv_bounds[3] = GMT->common.R.wesn[YHI];
		ssv_bounds_set = true;
	}

	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  Program <%s>\n", program_name);
		fprintf(stderr, "dbg2  MB-system Version %s\n", MB_VERSION);
		fprintf(stderr, "dbg2  Control Parameters:\n");
		fprintf(stderr, "dbg2       verbose:           %d\n", verbose);
		fprintf(stderr, "dbg2       help:              %d\n", Ctrl->X.help);
		fprintf(stderr, "dbg2       format:            %d\n", format);
		fprintf(stderr, "dbg2       pings:             %d\n", pings);
		fprintf(stderr, "dbg2       lonflip:           %d\n", lonflip);
		fprintf(stderr, "dbg2       bounds[0]:         %f\n", bounds[0]);
		fprintf(stderr, "dbg2       bounds[1]:         %f\n", bounds[1]);
		fprintf(stderr, "dbg2       bounds[2]:         %f\n", bounds[2]);
		fprintf(stderr, "dbg2       bounds[3]:         %f\n", bounds[3]);
		fprintf(stderr, "dbg2       btime_i[0]:        %d\n", btime_i[0]);
		fprintf(stderr, "dbg2       btime_i[1]:        %d\n", btime_i[1]);
		fprintf(stderr, "dbg2       btime_i[2]:        %d\n", btime_i[2]);
		fprintf(stderr, "dbg2       btime_i[3]:        %d\n", btime_i[3]);
		fprintf(stderr, "dbg2       btime_i[4]:        %d\n", btime_i[4]);
		fprintf(stderr, "dbg2       btime_i[5]:        %d\n", btime_i[5]);
		fprintf(stderr, "dbg2       btime_i[6]:        %d\n", btime_i[6]);
		fprintf(stderr, "dbg2       etime_i[0]:        %d\n", etime_i[0]);
		fprintf(stderr, "dbg2       etime_i[1]:        %d\n", etime_i[1]);
		fprintf(stderr, "dbg2       etime_i[2]:        %d\n", etime_i[2]);
		fprintf(stderr, "dbg2       etime_i[3]:        %d\n", etime_i[3]);
		fprintf(stderr, "dbg2       etime_i[4]:        %d\n", etime_i[4]);
		fprintf(stderr, "dbg2       etime_i[5]:        %d\n", etime_i[5]);
		fprintf(stderr, "dbg2       etime_i[6]:        %d\n", etime_i[6]);
		fprintf(stderr, "dbg2       speedmin:          %f\n", speedmin);
		fprintf(stderr, "dbg2       timegap:           %f\n", timegap);
		fprintf(stderr, "dbg2       svp_printmode:     %d\n", svp_printmode);
		fprintf(stderr, "dbg2       svp_file_output:   %d\n", svp_file_output);
		fprintf(stderr, "dbg2       svp_setprocess:    %d\n", svp_setprocess);
		fprintf(stderr, "dbg2       svp_force_zero:    %d\n", svp_force_zero);
		fprintf(stderr, "dbg2       ssv_output:        %d\n", ssv_output);
		fprintf(stderr, "dbg2       ssv_bounds_set:    %d\n", ssv_bounds_set);
		fprintf(stderr, "dbg2       ssv_bounds[0]:     %f\n", ssv_bounds[0]);
		fprintf(stderr, "dbg2       ssv_bounds[1]:     %f\n", ssv_bounds[1]);
		fprintf(stderr, "dbg2       ssv_bounds[2]:     %f\n", ssv_bounds[2]);
		fprintf(stderr, "dbg2       ssv_bounds[3]:     %f\n", ssv_bounds[3]);
	}

	if (format == 0) {
		mb_get_format(verbose, read_file, NULL, &format, &error);
		read_datalist = format < 0;
	}

	/* open file list */
	if (read_datalist) {
		char *pch;
		const int look_processed = MB_DATALIST_LOOK_UNSET;
		if (mb_datalist_open(verbose, &datalist, read_file, look_processed, &error) != MB_SUCCESS) {
			GMT_Report (API, GMT_MSG_NORMAL, "Unable to open data list file: %s\n", read_file);
			GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
			Return(MB_ERROR_OPEN_FAIL);
		}
		read_data = mb_datalist_read(verbose, datalist, file, dfile, &format, &file_weight, &error) == MB_SUCCESS;
#ifdef _WIN32
		pch = strchr(file, ':');
		if (error && pch && (pch - file) > 1) {
			int k, off = pch - file;
			for (k = 0; k < strlen(file) - off + 1; k++)
				file[k] = file[k+off-1];
			file[k] = '\0';
		}
#endif
	}
	else {		/* else copy single filename to be read */
		strcpy(file, read_file);
		read_data = true;
	}

	if (ssv_output || output_counts) {
		unsigned int n_col = output_counts ? 1 : 2;
		if ((error = GMT_Set_Columns (GMT->parent, GMT_OUT, n_col, GMT_COL_FIX_NO_TEXT)) != GMT_NOERROR)
			Return (error);
		if (GMT_Init_IO (API, GMT_IS_DATASET, GMT_IS_NONE, GMT_OUT, GMT_ADD_DEFAULT, 0, options) != GMT_NOERROR) 	/* Establishes data output */
			Return (API->error);
		if (GMT_Begin_IO (API, GMT_IS_DATASET, GMT_OUT, GMT_HEADER_OFF) != GMT_NOERROR)
			Return (API->error);
		Out = gmt_new_record(GMT, ssv_out, NULL);
	}
	else {
		if ((D = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_LINE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
			Return (API->error);	/* An empty dataset */
	}

	ix = (GMT->current.setting.io_lonlat_toggle[GMT_IN]);	iy = 1 - ix;

	/* loop over all files to be read */
	while (read_data) {
		/* check format and get data sources */
		if ((status = mb_format_source(verbose, &format, &platform_source, &nav_source, &sensordepth_source, &heading_source, &attitude_source, &svp_source, &error)) == MB_FAILURE) {
			char *message;
			mb_error(verbose, error, &message);
			GMT_Report (API, GMT_MSG_NORMAL, "MBIO Error returned from function <mb_format_source>:\n%s\n", message);
			GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
			Return (error);
		}

		/* initialize reading the swath file */
		if (mb_read_init(verbose, file, format, pings, lonflip, bounds, btime_i, etime_i, speedmin, timegap,
						 &mbio_ptr, &btime_d, &etime_d, &beams_bath, &beams_amp, &pixels_ss, &error) != MB_SUCCESS) {
			char *message;
			mb_error(verbose, error, &message);
			GMT_Report (API, GMT_MSG_NORMAL, "MBIO Error returned from function <mb_read_init>:\n%s\n", message);
			GMT_Report (API, GMT_MSG_NORMAL, "Multibeam File <%s> not initialized for reading\n", file);
			GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
			Return (error);
		}

		/* allocate memory for data arrays */
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), (void **)&beamflag, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&bath, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(double), (void **)&amp, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&bathacrosstrack, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&bathalongtrack, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&ss, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&ssacrosstrack, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&ssalongtrack, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&ttimes, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&angles, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&angles_forward, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&angles_null, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&heave, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&alongtrack_offset, &error);

		/* if error initializing memory then quit */
		if (error != MB_ERROR_NO_ERROR) {
			char *message;
			mb_error(verbose, error, &message);
			GMT_Report (API, GMT_MSG_NORMAL, "MBIO Error allocating data arrays:\n%s\n", message);
			GMT_Report (API, GMT_MSG_NORMAL, "Program Terminated\n");
			Return (error);
		}

		/* output info */
		if (verbose >= 1) {
			if (ssv_output)
				GMT_Report (API, GMT_MSG_INFORMATION, "Searching %s for SSV records\n", file);
			else
				GMT_Report (API, GMT_MSG_INFORMATION, "Searching %s for SVP records\n", file);
		}

		/* read and print data */
		svp.n = 0;
		while (error <= MB_ERROR_NO_ERROR) {		/* read a data record */
			status = mb_get_all(verbose, mbio_ptr, &store_ptr, &kind, time_i, &time_d, &navlon, &navlat, &speed, &heading,
								&distance, &altitude, &sonardepth, &beams_bath, &beams_amp, &pixels_ss, beamflag, bath, amp,
								bathacrosstrack, bathalongtrack, ss, ssacrosstrack, ssalongtrack, comment, &error);

			if (verbose >= 2) {
				fprintf(stderr, "\ndbg2  Ping read in program <%s>\n", program_name);
				fprintf(stderr, "dbg2       kind:           %d\n", kind);
				fprintf(stderr, "dbg2       error:          %d\n", error);
				fprintf(stderr, "dbg2       status:         %d\n", status);
			}

			/* if svp then extract data */
			if (error <= MB_ERROR_NO_ERROR && kind == svp_source && svp_source != MB_DATA_NONE) {
				/* extract svp */
				status = mb_extract_svp(verbose, mbio_ptr, store_ptr, &kind, &svp.n, svp.depth, svp.velocity, &error);
				if (status == MB_SUCCESS) {
					svp_read++;
					svp_loaded = true;
					svp.match_last = false;
					svp.repeat_in_file = false;
					if (last_time_d != 0.0) {
						svp.time_set = true;
						svp.time_d = last_time_d;
					}
					else {
						svp.time_set = false;
						svp.time_d = 0.0;
					}
					if (navlon != 0.0 || navlat != 0.0) {
						svp.position_set = true;
						svp.longitude = navlon;
						svp.latitude = navlat;
					}
					else if (last_navlon != 0.0 || last_navlat != 0.0) {
						svp.position_set = true;
						svp.longitude = last_navlon;
						svp.latitude = last_navlat;
					}
					else {
						svp.position_set = false;
						svp.longitude = 0.0;
						svp.latitude = 0.0;
					}
					svp.depthzero_reset = false;
					svp.depthzero = 0.0;
				}
				else {
					svp_loaded = false;
				}

				/* force zero depth if requested */
				if (svp_loaded && svp.n > 0 && svp_force_zero && svp.depth[0] != 0.0) {
					svp.depthzero = svp.depth[0];
					svp.depth[0] = 0.0;
					svp.depthzero_reset = true;
				}

				/* check if the svp is a duplicate to a previous svp in the same file */
				if (svp_loaded) {
					svp_match_last = false;
					for (int j = 0; j < svp_save_count && svp_match_last; j++) {
						if (svp.n == svp_save[j].n && memcmp(svp.depth, svp_save[j].depth, svp.n) == 0 &&
							memcmp(svp.velocity, svp_save[j].velocity, svp.n) == 0) {
							svp_match_last = true;
						}
					}
					svp.match_last = svp_match_last;

					/* check if svp is the same as the previous */
					if (svp.n == svp_last.n && memcmp(svp.depth, svp_last.depth, svp.n) == 0 &&
						memcmp(svp.velocity, svp_last.velocity, svp.n) == 0) {
						svp_repeat_in_file = true;
					}
					else {
						svp_repeat_in_file = false;
					}
					svp.repeat_in_file = svp_repeat_in_file;

					/* save the svp */
					svp_last.time_set = false;
					svp_last.position_set = false;
					svp_last.n = svp.n;
					for (int i = 0; i < svp.n; i++) {
						svp_last.depth[i] = svp.depth[i];
						svp_last.velocity[i] = svp.velocity[i];
					}
				}

				/* if the svp is unique so far, save it in memory */
				if (svp_loaded && !svp_match_last && svp.n >= min_num_pairs) {
					if (svp_save_count >= svp_save_alloc) {		/* allocate memory as needed */
						svp_save_alloc += MBSVPLIST_SVP_NUM_ALLOC;
						status = mb_reallocd(verbose, __FILE__, __LINE__, svp_save_alloc * sizeof(struct mbsvplist_svp_struct),
											 (void **)&svp_save, &error);
					}

					/* save the svp */
					svp_save[svp_save_count].time_set = svp.time_set;
					svp_save[svp_save_count].position_set = svp.position_set;
					svp_save[svp_save_count].repeat_in_file = svp.repeat_in_file;
					svp_save[svp_save_count].match_last = svp.match_last;
					svp_save[svp_save_count].time_d = svp.time_d;
					svp_save[svp_save_count].longitude = svp.longitude;
					svp_save[svp_save_count].latitude = svp.latitude;
					svp_save[svp_save_count].n = svp.n;
					for (int i = 0; i < svp.n; i++) {
						svp_save[svp_save_count].depth[i] = svp.depth[i];
						svp_save[svp_save_count].velocity[i] = svp.velocity[i];
					}
					svp_save_count++;
					svp_unique++;
				}
			}

			/* else if survey data save most recent ping time and if ssv output desired call mb_ttimes() and output ssv */
			else if (error <= MB_ERROR_NO_ERROR && kind == MB_DATA_DATA) {
				/* save most recent survey time stamp and position */
				last_time_d = time_d;
				last_navlon = navlon;
				last_navlat = navlat;

				/* check if any saved svps need time tags and position */
				if (time_d != 0.0 && (navlon != 0.0 || navlat != 0.0)) {
					for (int isvp = 0; isvp < svp_save_count; isvp++) {
						if (!svp_save[isvp].time_set) {
							svp_save[isvp].time_set = true;
							svp_save[isvp].time_d = time_d;
						}
						if (!svp_save[isvp].position_set) {
							svp_save[isvp].position_set = true;
							svp_save[isvp].longitude = navlon;
							svp_save[isvp].latitude = navlat;
						}
					}
				}

				/* if desired output ssv_output */
				if (ssv_output) {
					/* extract ttimes */
					status = mb_ttimes(verbose, mbio_ptr, store_ptr, &kind, &nbeams, ttimes, angles, angles_forward, angles_null,
									   heave, alongtrack_offset, &sonardepth, &ssv, &error);

					/* output ssv */
					if (status == MB_SUCCESS) {
						if (!ssv_bounds_set || (navlon >= ssv_bounds[0] && navlon <= ssv_bounds[1] &&
												navlat >= ssv_bounds[2] && navlat <= ssv_bounds[3]))
							Out->data[0] = sonardepth;
							Out->data[1] = ssv;
							GMT_Put_Record(API, GMT_WRITE_DATA, Out);
					}
				}
			}
		}

		status &= mb_close(verbose, &mbio_ptr, &error);

		/* output svps from this file if there are any and ssv_output and output_counts are false */
		if (svp_save_count > 0 && !ssv_output && !output_counts) {
			for (int isvp = 0; isvp < svp_save_count; isvp++) {
				if (svp_save[isvp].n >= min_num_pairs &&
					((svp_printmode == MBSVPLIST_PRINTMODE_CHANGE && (svp_written == 0 || !svp_save[isvp].repeat_in_file)) ||
					 (svp_printmode == MBSVPLIST_PRINTMODE_UNIQUE && !svp_save[isvp].match_last) ||
					 (svp_printmode == MBSVPLIST_PRINTMODE_ALL))) {
					/* set the output */
					if (svp_file_output) {
						sprintf(svp_file, "%s_%3.3d.svp", file, isvp);
						if (Ctrl->X.out_file) free(Ctrl->X.out_file);	/* First free any previous fname */
						Ctrl->X.out_file = strdup(svp_file);
					}

					/* get time as date */
					mb_get_date(verbose, svp_save[isvp].time_d, svp_time_i);

					/* print out the svp */
					if (output_as_table) {		/* output csv table to stdout */
						if (out_cnt == 0) {		/* output header records */
							printf("#mbsvplist CSV table output\n#navigation information is "
							       "approximate\n#SVP_cnt,date_time,longitude,latitude,num_data_points\n");
						}
						out_cnt++;
						printf("%d,%4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d.%6.6d,%.6f,%.6f,%d\n", out_cnt, svp_time_i[0],
						        svp_time_i[1], svp_time_i[2], svp_time_i[3], svp_time_i[4], svp_time_i[5], svp_time_i[6],
						        svp_save[isvp].longitude, svp_save[isvp].latitude, svp_save[isvp].n);
					}
					else {
						const time_t right_now = time((time_t *)0);
						char date[32] = {""};
						char *user_ptr = getenv("USER");
						char user[MB_PATH_MAXLINE] = {""};
						char host[MB_PATH_MAXLINE] = {""};
						struct GMT_DATASEGMENT_HIDDEN *SH = NULL;

						/* output info */
						if (verbose >= 1)
							GMT_Report (API, GMT_MSG_NORMAL, "Outputting SVP to file: %s (# svp pairs=%d)\n", svp_file, svp_save[isvp].n);

						/* write it out */
						S = GMT_Alloc_Segment (API, GMT_NO_STRINGS, svp_save[isvp].n, 2, NULL, NULL);
						snprintf (header, GMT_BUFSIZ, "\n## MB-SVP %4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d.%6.6d %.9f %.9f",
						          svp_time_i[0], svp_time_i[1], svp_time_i[2], svp_time_i[3], svp_time_i[4], svp_time_i[5], svp_time_i[6], svp_save[isvp].longitude, svp_save[isvp].latitude);
						snprintf(buf512, GMT_LEN512, "\n## Water Sound Velocity Profile (SVP)"); strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## Output by Program %s", program_name); strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## MB-System Version %s", MB_VERSION); strcat(header, buf512);
						strcpy(date, ctime(&right_now));
						date[strlen(date) - 1] = '\0';
						if (user_ptr == NULL)
							user_ptr = getenv("LOGNAME");
						if (user_ptr != NULL)
							strcpy(user, user_ptr);
						else
							strcpy(user, "unknown");
						gethostname(host, MB_PATH_MAXLINE);
						snprintf(buf512, GMT_LEN512, "\n## Run by user <%s> on cpu <%s> at <%s>", user, host, date); strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## Swath File: %s", file); strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## Start Time: %4.4d/%2.2d/%2.2d %2.2d:%2.2d:%2.2d.%6.6d", svp_time_i[0],
								svp_time_i[1], svp_time_i[2], svp_time_i[3], svp_time_i[4], svp_time_i[5], svp_time_i[6]);
						strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## SVP Longitude: %f", svp_save[isvp].longitude); strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## SVP Latitude:  %f", svp_save[isvp].latitude); strcat(header, buf512);
						snprintf(buf512, GMT_LEN512, "\n## SVP Count: %d", svp_save_count); strcat(header, buf512);
						if (fabs(svp_save[isvp].depthzero) > 0.01) {
							snprintf(buf512, GMT_LEN512, "\n## Initial depth reset from %f to 0.0 meters", svp_save[isvp].depthzero);
							strcat(header, buf512);
						}
						snprintf(buf512, GMT_LEN512, "\n## Number of SVP Points: %d", svp_save[isvp].n);	strcat(header, buf512);

						S->header = strdup(header);
						gmt_M_memcpy (S->data[ix], svp_save[isvp].depth, svp_save[isvp].n, double);
						gmt_M_memcpy (S->data[iy], svp_save[isvp].velocity, svp_save[isvp].n, double);
						S->n_rows = svp_save[isvp].n;

						if (n_seg == n_seg_alloc) {
							size_t old_n_alloc = n_seg_alloc;
							D->table[tbl]->segment = gmt_M_memory (GMT, D->table[tbl]->segment, (n_seg_alloc += 1), struct GMT_DATASEGMENT *);
							gmt_M_memset (&(D->table[tbl]->segment[old_n_alloc]), n_seg_alloc - old_n_alloc, struct GMT_DATASEGMENT *);	/* Set to NULL */
						}
						D->table[tbl]->segment[n_seg++] = S;
						D->table[tbl]->n_segments++;		D->n_segments++;
						D->table[tbl]->n_records += svp_save[isvp].n;	D->n_records += svp_save[isvp].n;

						svp_written++;
					}

					/* if desired, set first svp output to be used for recalculating bathymetry */
					if (svp_file_output && svp_setprocess && svp_save_count == 1)
						status = mb_pr_update_svp(verbose, file, true, svp_file, MBP_ANGLES_OK, true, &error);
				}
			}
		}

		/* update total counts */
		svp_read_tot += svp_read;
		svp_unique_tot += svp_unique;
		svp_written_tot += svp_written;

		/* output info */
		if (verbose >= 1) {
			GMT_Report (API, GMT_MSG_INFORMATION, "%d SVP records read\n", svp_read);
			GMT_Report (API, GMT_MSG_INFORMATION, "%d SVP unique records read\n", svp_unique);
			GMT_Report (API, GMT_MSG_INFORMATION, "%d SVP records written\n", svp_written);
		}

		/* figure out whether and what to read next */
		if (read_datalist) {
			read_data = mb_datalist_read(verbose, datalist, file, dfile, &format, &file_weight, &error) == MB_SUCCESS;
			tbl++;
			D->table = gmt_M_memory (GMT, D->table, tbl+1, struct GMT_DATATABLE *);
			if ((D->table[tbl] = gmt_create_table (GMT, 0, 0, dim[3], 0, false)) == NULL)
				Return (API->error);
			D->n_tables++;
			n_seg = n_seg_alloc = 0;		/* Reset these since a new table is on the way */
		}
		else
			read_data = false;

	}		/* end loop over files in list */


	if (D) {		/* Not all paths created a DATASET */
		if (GMT_Write_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_LINE, GMT_WRITE_SET, NULL, Ctrl->X.out_file, D) != GMT_NOERROR)
			Return (API->error);
	}
	else if (output_counts) {
		Out->data[0] = svp_unique_tot;
		GMT_Put_Record(API, GMT_WRITE_DATA, Out);
	}

	if (read_datalist)
		mb_datalist_close(verbose, &datalist, &error);

	/* output info */
	if (verbose >= 1) {
		GMT_Report (API, GMT_MSG_NORMAL, "Total %d SVP records read\n", svp_read_tot);
		GMT_Report (API, GMT_MSG_NORMAL, "Total %d SVP unique records found\n", svp_unique_tot);
		GMT_Report (API, GMT_MSG_NORMAL, "Total %d SVP records written\n", svp_written_tot);
	}

	/* deallocate memory */
	status &= mb_freed(verbose, __FILE__, __LINE__, (void **)&svp_save, &error);

	/* check memory */
	if (verbose >= 4)
		status &= mb_memory_list(verbose, &error);

	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  Program <%s> completed\n", program_name);
		fprintf(stderr, "dbg2  Ending status:\n");
		fprintf(stderr, "dbg2       status:  %d\n", status);
	}

	if (ssv_output || output_counts) {
		gmt_M_free (GMT, Out);
		if (GMT_End_IO (API, GMT_OUT, 0) != GMT_NOERROR) 	/* Disables further data output */
			Return (API->error);
	}
	
	Return (GMT_NOERROR);
}
/*--------------------------------------------------------------------*/
