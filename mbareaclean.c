/*--------------------------------------------------------------------
 *    The MB-system:	mbareaclean.c	2/27/2003
 *    $Id: mbareaclean.c 2298 2017-04-10 07:57:48Z caress $
 *
 *    Copyright (c) 2003-2016 by
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
 * mbareaclean identifies and flags artifacts in swath sonar bathymetry data.
 * The edit events are output to edit save files which can be applied
 * to the data by the program mbprocess. These are the same edit save
 * files created and/or modified by mbclean and mbedit.
 * The input data are one swath file or a datalist referencing multiple
 * swath files. An area is specified in longitude and latitude bounds,
 * along with a bin size in meters. The area is divided into a grid with
 * square cells of the specified bin size. As the data are read, each of
 * the soundings that fall within one of the bins is stored. Once all of
 * data are read, one or more statistical tests are performed on the soundings
 * within each bin, providing there are a sufficient number of soundings.
 * The user may specify one or both of the following actions:
 *   1) Previously unflagged soundings that fail a test are flagged as bad.
 *   2) Previously flagged soundings that pass all tests are unflagged.
 * If a sounding's flag status is changed, that flagging action is output
 * to the edit save file of the swath file containing that sounding. This
 * program will create edit save files if necessary, or append to those that
 * already exist.
 *
 * Author:	D. W. Caress
 * Date:	February 27, 2003
 *		Amsterdam Airport
 */

#define THIS_MODULE_NAME	"mbareaclean"
#define THIS_MODULE_LIB		"mbgmt"
#define THIS_MODULE_PURPOSE	"Automatically flag bad beams in swath sonar bathymetry data"
#define THIS_MODULE_KEYS	""
#define THIS_MODULE_NEEDS	""
#define THIS_MODULE_OPTIONS "->V"

/* GMT5 header file */
#include "gmt_dev.h"

EXTERN_MSC int GMT_mbareaclean(void *API, int mode, void *args);

#define GMT_PROG_OPTIONS "->V"

/* mbio include files */
#include "mb_status.h"
#include "mb_format.h"
#include "mb_define.h"
#include "mb_io.h"
#include "mb_swap.h"
#include "mb_process.h"
#include "mb_info.h"

/* allocation */
#define FILEALLOCNUM	16
#define PINGALLOCNUM	128
#define SNDGALLOCNUM	128

/* Redefines to use the GMT print mechanism that also works in Matlab */
#define PR GMT_Report(API, GMT_MSG_NORMAL, output); 
#define print sprintf
#define stderr_ output
/* Replacing the above defines by these should get back the original fprintf(stderr,...) behaviour
#define stderr_ stderr
#define PR
#define print fprintf
*/

struct mbareaclean_file_struct {
	char	filelist[MB_PATH_MAXLINE];
	int	file_format;
	int	nping;
	int	nping_alloc;
	int	nnull;
	int	nflag;
	int	ngood;
	int	nunflagged;
	int	nflagged;
	double	*ping_time_d;
	int	*pingmultiplicity;
	double	*ping_altitude;
	int	nsndg;
	int	nsndg_alloc;
	int	sndg_countstart;
	int	beams_bath;
	struct mbareaclean_sndg_struct *sndg;
};
struct mbareaclean_sndg_struct {
	int	sndg_file;
	int	sndg_ping;
	int	sndg_beam;
	double	sndg_depth;
	double	sndg_x;
	double	sndg_y;
	char	sndg_beamflag_org;
	char	sndg_beamflag_esf;
	char	sndg_beamflag;
	char	sndg_edit;
};

/* Control structure for mbgetdata */
GMT_LOCAL struct MBAREACLEAN_CTRL {
	int     read_datalist;
	void	*datalist;

	struct mbareaclean_B {	/*  */
		bool active;
	} B;
	struct mbareaclean_D {	/*  */
		bool active;
		int  std_dev_nmin;
		double std_dev_threshold;
	} D;
	struct mbareaclean_F {	/* -F<format> */
		bool active;
		int format;
	} F;
	struct mbareaclean_G {	/*  */
		bool active;
	} G;
	struct mbareaclean_I {	/* -I<inputfile> */
		bool active;
		char *file;
	} I;
	struct mbareaclean_L {	/* -L */
		bool active;
	} L;
	struct mbareaclean_M {	/*  */
		bool active;
		bool density_filter;
		int median_filter_nmin, density_filter_nmax;
		double median_filter_threshold;
	} M;
	struct mbareaclean_N {	/* -N */
		bool active;
		int  min_beam, beam_in, max_beam_no, max_beam;
	} N;
	struct mbareaclean_P {	/*  */
		bool   active;
		int	plane_fit_nmin;
		double	plane_fit_threshold;
	} P;
	struct mbareaclean_R {	/*  */
		bool active;
		double	areabounds[4];
	} R;
	struct mbareaclean_S {	/*  */
		bool active;
		double binsize;
	} S;
	struct mbareaclean_T {	/*  */
		bool active;
		int flag_detect;
	} T;
};

GMT_LOCAL void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct  MBAREACLEAN_CTRL *Ctrl;

	Ctrl = gmt_M_memory (GMT, NULL, 1, struct MBAREACLEAN_CTRL);

	Ctrl->I.file = NULL;
		
	/* mbswath variables */
	Ctrl->read_datalist = false;
	Ctrl->datalist = NULL;
	Ctrl->D.std_dev_threshold = 2;
	Ctrl->D.std_dev_nmin = 10;
	Ctrl->F.format = 0;
	Ctrl->M.median_filter_threshold = 0.25;
	Ctrl->M.median_filter_nmin = 10;
	Ctrl->M.density_filter_nmax = 0;
	Ctrl->N.min_beam = 0;
	Ctrl->N.max_beam = 0;
	Ctrl->N.max_beam_no = 0;
	Ctrl->N.beam_in = true;
	Ctrl->P.plane_fit_threshold = 0.05;
	Ctrl->P.plane_fit_nmin = 10;
	Ctrl->S.binsize = 0;
	Ctrl->T.flag_detect = MB_DETECT_AMPLITUDE;

	return (Ctrl);
}

GMT_LOCAL void Free_Ctrl (struct GMT_CTRL *GMT, struct MBAREACLEAN_CTRL *Ctrl) {	/* Deallocate control structure */
	if (!Ctrl) return;
	if (Ctrl->I.file) free (Ctrl->I.file);
	gmt_M_free (GMT, Ctrl);
}

GMT_LOCAL int usage (struct GMTAPI_CTRL *API, int level) {
	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message(API, GMT_TIME_NONE, "usage: mbareaclean -Iinfile -PPARAMETER:value [-E -L -N -V]\n");
	GMT_Message(API, GMT_TIME_NONE, "MB-System Version %s\n",MB_VERSION);

	if (level == GMT_SYNOPSIS) return (EXIT_FAILURE);

	return (EXIT_FAILURE);
}

/* sounding storage values and arrays */
/* THIS IS THE REMAINING GLOBAL. GLOBALS ARE AWFULL AS THEY SURVIVE BETWEEN CALLS FROM EXTERNAL INTERFACES (e.g. MATLAB) */
struct mbareaclean_file_struct 	*files = NULL;

/* sounding pointer resolving function */
int getsoundingptr(struct GMTAPI_CTRL *API, int verbose, int soundingid, int nfile,
                   struct mbareaclean_sndg_struct **sndgptr, int *error);
int flag_sounding(struct GMTAPI_CTRL *API, int verbose, int flag, int output_bad, int output_good,
                  struct mbareaclean_sndg_struct *sndg, int *error);

static char rcs_id[] = "$Id: mbareaclean.c 2298 2017-04-10 07:57:48Z caress $";

GMT_LOCAL int parse (struct GMT_CTRL *GMT, struct MBAREACLEAN_CTRL *Ctrl, struct GMT_OPTION *options) {
	/* This parses the options provided to mbswath and sets parameters in Ctrl.
	 * Note Ctrl has already been initialized and non-zero default values set.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	unsigned int n_errors = 0, n_files = 0;
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;
	int    n, i1, i2;
	double	d1, d2;

	for (opt = options; opt; opt = opt->next) {	/* Process all the options given */

		switch (opt->option) {
			case '<':	/* Input file (only one or three is accepted) */
				Ctrl->I.active = true;
				if (gmt_check_filearg (GMT, '<', opt->arg, GMT_IN, GMT_IS_DATASET)) {
					Ctrl->I.file = strdup (opt->arg);
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error: only one input file is allowed.\n");
					n_errors++;
				}
				break;

			/* Processes program-specific parameters */

			case 'B':
				Ctrl->B.active = true;
				break;
			case 'D':
				Ctrl->D.active = true;
				sscanf (opt->arg,"%lf/%d", &Ctrl->D.std_dev_threshold,&Ctrl->D.std_dev_nmin);
				break;
			case 'F':	/* format */
				n = sscanf(opt->arg, "%d", &(Ctrl->F.format));
				if (n == 1)
					Ctrl->F.active = true;
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -F option: \n");
					n_errors++;
				}
				break;
			case 'I':	/* -I<inputfile> */
				Ctrl->I.active = true;
				if (!gmt_access(GMT, opt->arg, R_OK)) {	/* Got a file */
					Ctrl->I.file = strdup(opt->arg);
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -I: Requires a valid file\n");
					n_errors++;
				}
				break;
			case 'L':	/* -L */
				Ctrl->L.active = true;
				break;
			case 'M':	/*  */
				Ctrl->M.active = true;
				n = sscanf (opt->arg,"%lf/%d/%d", &d1,&i1,&i2);
				if (n > 0) Ctrl->M.median_filter_threshold = d1;
				if (n > 1) Ctrl->M.median_filter_nmin = i1;
				if (n > 2) {
					Ctrl->M.density_filter = true;
					Ctrl->M.density_filter_nmax = i2;
				}
				break;
			case 'N':	/*  */
				Ctrl->N.active = true;
				sscanf (opt->arg,"%d/%d", &Ctrl->N.min_beam,&Ctrl->N.max_beam_no);
				if (opt->arg[0] == '-') {
					Ctrl->N.min_beam = -Ctrl->N.min_beam;
					Ctrl->N.beam_in = false;
				}
				if (Ctrl->N.max_beam_no < 0)
					Ctrl->N.max_beam_no = -Ctrl->N.max_beam_no;
				Ctrl->N.max_beam = Ctrl->N.max_beam_no;
				if (Ctrl->N.max_beam < Ctrl->N.min_beam)
					Ctrl->N.max_beam = Ctrl->N.min_beam;
				break;
			case 'P':	/*  */
				Ctrl->P.active = true;
				sscanf (opt->arg,"%lf", &Ctrl->P.plane_fit_threshold);
				n = sscanf (opt->arg,"%lf/%d/%lf", &d1,&i1,&d2);
				if (n > 0) Ctrl->P.plane_fit_threshold = d1;
				if (n > 1) Ctrl->P.plane_fit_nmin = i1;
				break;
			case 'R':	/*  */
				Ctrl->R.active = true;
				mb_get_bounds(opt->arg, Ctrl->R.areabounds);
				break;
			case 'S':	/*  */
				Ctrl->S.active = true;
				sscanf(opt->arg,"%lf", &Ctrl->S.binsize);
				break;
			case 'T':	/*  */
				Ctrl->T.active = true;
				sscanf(opt->arg,"%d", &Ctrl->T.flag_detect);
				break;
			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, n_files != 1, "Syntax error: Must specify one input file(s)\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->I.active && !Ctrl->I.file,
	                                   "Syntax error -I option: Must specify input file\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

/*--------------------------------------------------------------------*/
int GMT_mbareaclean (void *V_API, int mode, void *args) {

	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;	/* General GMT interal parameters */
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */
	struct MBAREACLEAN_CTRL *Ctrl = NULL;

	char program_name[] = "MBAREACLEAN";
	char help_message[] =  "MBAREACLEAN identifies and flags artifacts in swath bathymetry data";
	char usage_message[] = "mbareaclean [-Fformat -Iinfile -Rwest/east/south/north -B -G -Sbinsize	\n\t -Mthreshold/nmin -Dthreshold[/nmin[/nmax]] -Ttype -N[-]minbeam/maxbeam]";
	extern char *optarg;
	int	errflg = 0, help = 0;

	/* sounding storage values and arrays (they used to be globals) */
	int	nfile = 0, nfile_alloc = 0, nsndg = 0;
	int	**gsndg = NULL;
	int	*gsndgnum = NULL;
	int	*gsndgnum_alloc = NULL;
	struct mbareaclean_sndg_struct *sndg = NULL;

	/* MBIO status variables */
	int	status;
	int	verbose = 0;
	int	error = MB_ERROR_NO_ERROR;
	char	*message = NULL;

	/* MBIO read control parameters */
	void	*mbio_ptr = NULL;
	void	*store_ptr = NULL;
	int	kind;
	char	swathfile[MB_PATH_MAXLINE];
	char	swathfileread[MB_PATH_MAXLINE];
	char	dfile[MB_PATH_MAXLINE];
	void	*datalist;
	int	look_processed = MB_DATALIST_LOOK_UNSET;
	bool	read_data;
	double	file_weight;
	int	formatread;
	int	variable_beams;
	int	traveltime;
	int	beam_flagging;
	int	pings;
	int	lonflip;
	double	bounds[4];
	int	btime_i[7], etime_i[7];
	double	btime_d, etime_d, speedmin, timegap;
	struct mb_info_struct mb_info;

	int	time_i[7];
	double	time_d;
	int	pingsread;
	double	navlon;
	double	navlat;
	double	speed;
	double	heading;
	double	distance;
	double	altitude;
	double	sonardepth;
	int	beams_bath;
	int	beams_amp;
	int	pixels_ss;
	char	*beamflag;
	char	*beamflagorg;
	int	*detect;
	double	*bath;
	double	*amp;
	double	*bathlon;
	double	*bathlat;
	double	*ss;
	double	*sslon;
	double	*sslat;
	char	comment[MB_COMMENT_MAXLINE];

	/* mbareaclean control parameters */
	int	nx, ny, detect_status, detect_error;
	double	dx, dy, mtodeglon, mtodeglat, mean, std_dev;

	/* median filter parameters */
	int	binnum, binnummax;
	double	*bindepths;
	double	threshold, median_depth;

	/* counting parameters */
	int	files_tot = 0;
	int	pings_tot = 0;
	int	beams_tot = 0;
	int	beams_good_org_tot = 0;
	int	beams_flag_org_tot = 0;
	int	beams_null_org_tot = 0;
	int	pings_file = 0;
	int	beams_file = 0;
	int	beams_good_org_file = 0;
	int	beams_flag_org_file = 0;
	int	beams_null_org_file = 0;

	/* save file control variables */
	int	esffile_open = MB_NO, action;
	char	esffile[MB_PATH_MAXLINE];
	char	output[MB_PATH_MAXLINE];
	struct mb_esf_struct esf;

	int	flagsounding;
	int	done;
	int	i, j, ix, iy, ib, kgrid;
	double	xx, yy, median_depth_low, median_depth_high;

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage (API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE) bailout (usage (API, GMT_USAGE));	/* Return the usage message */
	if (options->option == GMT_OPT_SYNOPSIS) bailout (usage (API, GMT_SYNOPSIS));	/* Return the synopsis */

	/* Parse the command-line arguments */

#if GMT_MAJOR_VERSION >= 6
	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
#else
	GMT = gmt_begin_module (API, THIS_MODULE_LIB, THIS_MODULE_NAME, &GMT_cpy); /* Save current state */
#endif
	if (GMT_Parse_Common (API, GMT_PROG_OPTIONS, options)) Return (API->error);
	Ctrl = (struct MBAREACLEAN_CTRL *) New_Ctrl (GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse (GMT, Ctrl, options))) Return (error);

	/*---------------------------- This is the mbgetdata main code ----------------------------*/
	
	/* set verbosity */
	verbose = GMT->common.V.active;
	
	/* get current default values */
	status = mb_defaults(verbose,&Ctrl->F.format,&pings,&lonflip,bounds, btime_i,etime_i,&speedmin,&timegap);

	/* reset all defaults but the format and lonflip */
	pings = 1;
	bounds[0] = -360.;
	bounds[1] = 360.;
	bounds[2] = -90.;
	bounds[3] = 90.;
	btime_i[0] = 1962;
	btime_i[1] = 2;
	btime_i[2] = 21;
	btime_i[3] = 10;
	btime_i[4] = 30;
	btime_i[5] = 0;
	btime_i[6] = 0;
	etime_i[0] = 2062;
	etime_i[1] = 2;
	etime_i[2] = 21;
	etime_i[3] = 10;
	etime_i[4] = 30;
	etime_i[5] = 0;
	etime_i[6] = 0;
	speedmin = 0.0;
	timegap = 1000000000.0;

	/* if error flagged then print it and exit */
	if (errflg) {
		print(stderr_,"usage: %s\n", usage_message);PR
		print(stderr_,"Program <%s> Terminated\n", program_name);PR
		error = MB_ERROR_BAD_USAGE;
		Return(error);
	}

	/* turn on median filter if nothing specified */
	if (!Ctrl->M.active && Ctrl->P.active && !Ctrl->D.active)
		Ctrl->M.active = true;

	/* turn on output bad if nothing specified */
	if (!Ctrl->B.active && !Ctrl->G.active) Ctrl->B.active = true;

	/* print starting message */
	if (verbose == 1 || help) {
		print(stderr_,"Program %s\n",program_name);PR
		print(stderr_,"Version %s\n",rcs_id);PR
		print(stderr_,"MB-system Version %s\n",MB_VERSION);PR
	}

	/* print starting debug statements */
	if (verbose >= 2) {
		print(stderr_,"\ndbg2  Program <%s>\n",program_name);PR
		print(stderr_,"dbg2  Version %s\n",rcs_id);PR
		print(stderr_,"dbg2  MB-system Version %s\n",MB_VERSION);PR
		print(stderr_,"dbg2  Control Parameters:\n");PR
		print(stderr_,"dbg2       verbose:        %d\n",verbose);PR
		print(stderr_,"dbg2       help:           %d\n",help);PR
		print(stderr_,"dbg2       pings:          %d\n",pings);PR
		print(stderr_,"dbg2       lonflip:        %d\n",lonflip);PR
		print(stderr_,"dbg2       bounds[0]:      %f\n",bounds[0]);PR
		print(stderr_,"dbg2       bounds[1]:      %f\n",bounds[1]);PR
		print(stderr_,"dbg2       bounds[2]:      %f\n",bounds[2]);PR
		print(stderr_,"dbg2       bounds[3]:      %f\n",bounds[3]);PR
		print(stderr_,"dbg2       btime_i[0]:     %d\n",btime_i[0]);PR
		print(stderr_,"dbg2       btime_i[1]:     %d\n",btime_i[1]);PR
		print(stderr_,"dbg2       btime_i[2]:     %d\n",btime_i[2]);PR
		print(stderr_,"dbg2       btime_i[3]:     %d\n",btime_i[3]);PR
		print(stderr_,"dbg2       btime_i[4]:     %d\n",btime_i[4]);PR
		print(stderr_,"dbg2       btime_i[5]:     %d\n",btime_i[5]);PR
		print(stderr_,"dbg2       btime_i[6]:     %d\n",btime_i[6]);PR
		print(stderr_,"dbg2       etime_i[0]:     %d\n",etime_i[0]);PR
		print(stderr_,"dbg2       etime_i[1]:     %d\n",etime_i[1]);PR
		print(stderr_,"dbg2       etime_i[2]:     %d\n",etime_i[2]);PR
		print(stderr_,"dbg2       etime_i[3]:     %d\n",etime_i[3]);PR
		print(stderr_,"dbg2       etime_i[4]:     %d\n",etime_i[4]);PR
		print(stderr_,"dbg2       etime_i[5]:     %d\n",etime_i[5]);PR
		print(stderr_,"dbg2       etime_i[6]:     %d\n",etime_i[6]);PR
		print(stderr_,"dbg2       speedmin:       %f\n",speedmin);PR
		print(stderr_,"dbg2       timegap:        %f\n",timegap);PR
		print(stderr_,"dbg2       data format:    %d\n",Ctrl->F.format);PR
		print(stderr_,"dbg2       input file:     %s\n",Ctrl->I.file);PR
		print(stderr_,"dbg2       median_filter:             %d\n",Ctrl->M.active);PR
		print(stderr_,"dbg2       median_filter_threshold:   %f\n",Ctrl->M.median_filter_threshold);PR
		print(stderr_,"dbg2       median_filter_nmin:        %d\n",Ctrl->M.median_filter_nmin);PR
		print(stderr_,"dbg2       density_filter:            %d\n",Ctrl->M.density_filter);PR
		print(stderr_,"dbg2       density_filter_nmax:       %d\n",Ctrl->M.density_filter_nmax);PR
		print(stderr_,"dbg2       plane_fit:                 %d\n",Ctrl->P.active);PR
		print(stderr_,"dbg2       plane_fit_threshold:       %f\n",Ctrl->P.plane_fit_threshold);PR
		print(stderr_,"dbg2       plane_fit_nmin:            %d\n",Ctrl->P.plane_fit_nmin);PR
		print(stderr_,"dbg2       std_dev_filter:            %d\n",Ctrl->D.active);PR
		print(stderr_,"dbg2       std_dev_threshold:         %f\n",Ctrl->D.std_dev_threshold);PR
		print(stderr_,"dbg2       std_dev_nmin:              %d\n",Ctrl->D.std_dev_nmin);PR
		print(stderr_,"dbg2       use_detect:                %d\n",Ctrl->T.active);PR
		print(stderr_,"dbg2       flag_detect:               %d\n",Ctrl->T.flag_detect);PR
		print(stderr_,"dbg2       limit_beams:               %d\n",Ctrl->N.active);PR
		print(stderr_,"dbg2       beam_in:                   %d\n",Ctrl->N.beam_in);PR
		print(stderr_,"dbg2       min_beam:                  %d\n",Ctrl->N.min_beam);PR
		print(stderr_,"dbg2       max_beam_no                %d\n",Ctrl->N.max_beam_no);PR
		print(stderr_,"dbg2       areaboundsset:  %d\n",Ctrl->R.active);PR
		print(stderr_,"dbg2       areabounds[0]:  %f\n",Ctrl->R.areabounds[0]);PR
		print(stderr_,"dbg2       areabounds[1]:  %f\n",Ctrl->R.areabounds[1]);PR
		print(stderr_,"dbg2       areabounds[2]:  %f\n",Ctrl->R.areabounds[2]);PR
		print(stderr_,"dbg2       areabounds[3]:  %f\n",Ctrl->R.areabounds[3]);PR
		print(stderr_,"dbg2       binsizeset:     %d\n",Ctrl->S.active);PR
		print(stderr_,"dbg2       binsize:        %f\n",Ctrl->S.binsize);PR
	}

	/* if help desired then print it and Return */
	if (help) {
		print(stderr_,"\n%s\n",help_message);PR
		print(stderr_,"\nusage: %s\n", usage_message);PR
		Return(error);
	}

	/* if bounds not set get bounds of input data */
	if (!Ctrl->R.active) {
		formatread = Ctrl->F.format;
		status = mb_get_info_datalist(verbose, Ctrl->I.file, &formatread, &mb_info, lonflip, &error);

		Ctrl->R.areabounds[0] = mb_info.lon_min;
		Ctrl->R.areabounds[1] = mb_info.lon_max;
		Ctrl->R.areabounds[2] = mb_info.lat_min;
		Ctrl->R.areabounds[3] = mb_info.lat_max;

		if (!Ctrl->S.active)
			Ctrl->S.binsize = 0.2 * mb_info.altitude_max;
	}

	/* calculate grid properties */
	mb_coor_scale(verbose,0.5*(Ctrl->R.areabounds[2]+Ctrl->R.areabounds[3]),&mtodeglon,&mtodeglat);
	if (Ctrl->S.binsize <= 0.0)
		Ctrl->S.binsize = (Ctrl->R.areabounds[1] - Ctrl->R.areabounds[0]) / 101 / mtodeglon;
	dx = Ctrl->S.binsize * mtodeglon;
	dy = Ctrl->S.binsize * mtodeglat;
	nx = 1 + (int)((Ctrl->R.areabounds[1] - Ctrl->R.areabounds[0]) / dx);
	ny = 1 + (int)((Ctrl->R.areabounds[3] - Ctrl->R.areabounds[2]) / dy);
	if (nx > 1 && ny > 1) {
		dx = (Ctrl->R.areabounds[1] - Ctrl->R.areabounds[0]) / (nx - 1);
		dy = (Ctrl->R.areabounds[3] - Ctrl->R.areabounds[2]) / (ny - 1);
	}

	/* allocate grid arrays */
	nsndg = 0;
	status = mb_mallocd(verbose,__FILE__,__LINE__, nx * ny * sizeof(int *), (void **)&gsndg, &error);
	if (status == MB_SUCCESS)
		status = mb_mallocd(verbose,__FILE__,__LINE__, nx * ny * sizeof(int), &gsndgnum, &error);
	if (status == MB_SUCCESS)
		status = mb_mallocd(verbose,__FILE__,__LINE__, nx * ny * sizeof(int), &gsndgnum_alloc, &error);

	/* if error initializing memory then quit */
	if (error != MB_ERROR_NO_ERROR) {
		mb_error(verbose,error,&message);
		print(stderr_,"MBIO Error allocating data arrays:\n%s\n",message);PR
		print(stderr_,"Program <%s> Terminated\n", program_name);PR
		Return(error);
	}

	/* if error initializing memory then quit */
	for (i=0;i<nx*ny;i++) {
		gsndg[i] = NULL;
		gsndgnum[i] = 0;
		gsndgnum_alloc[i] = 0;
	}

	/* give the statistics */
	if (verbose >= 0) {
		print(stderr_,"Area of interest:\n");PR
		print(stderr_,"     Minimum Longitude: %.6f Maximum Longitude: %.6f\n", Ctrl->R.areabounds[0],Ctrl->R.areabounds[1]);PR
		print(stderr_,"     Minimum Latitude:  %.6f Maximum Latitude:  %.6f\n", Ctrl->R.areabounds[2],Ctrl->R.areabounds[3]);PR
		print(stderr_,"     Bin Size:   %f\n", Ctrl->S.binsize);PR
		print(stderr_,"     Dimensions: %d %d\n", nx, ny);PR
		print(stderr_,"Cleaning algorithms:\n");PR
		if (Ctrl->M.active) {
			print(stderr_,"     Median filter: ON\n");PR
			print(stderr_,"     Median filter threshold:    %f\n", Ctrl->M.median_filter_threshold);PR
			print(stderr_,"     Median filter minimum N:    %d\n", Ctrl->M.median_filter_nmin);PR
		}
		else
			print(stderr_,"     Median filter: OFF\n");PR
		if (Ctrl->M.density_filter) {
			print(stderr_,"     Density filter: ON\n");PR
			print(stderr_,"     Density filter maximum N:    %d\n", Ctrl->M.density_filter_nmax);PR
		}
		else
			print(stderr_,"     Density filter: OFF\n");PR
		if (Ctrl->P.active) {
			print(stderr_,"     Plane fit:     ON\n");PR
			print(stderr_,"     Plane fit threshold:        %f\n", Ctrl->M.median_filter_threshold);PR
			print(stderr_,"     Plane fit minimum N:        %d\n", Ctrl->M.median_filter_nmin);PR
		}
		else
			print(stderr_,"     Plane fit:     OFF\n");PR
		if (Ctrl->D.active) {
			print(stderr_,"     Standard deviation filter: ON\n");PR
			print(stderr_,"     Standard deviation filter threshold:    %f\n", Ctrl->D.std_dev_threshold);PR
			print(stderr_,"     Standard deviation filter minimum N:    %d\n", Ctrl->D.std_dev_nmin);PR
		}
		else
			print(stderr_,"     Standard deviation filter: OFF\n");PR

		print(stderr_,"Restrictions:\n");PR
		if (Ctrl->T.active) {
			print(stderr_,"     Only flag if bottom detection algorithn is: ");PR
			if (Ctrl->T.flag_detect == MB_DETECT_UNKNOWN)
				{print(stderr_,"UNKNOWN\n");PR}
			else if (Ctrl->T.flag_detect == MB_DETECT_AMPLITUDE)
				{print(stderr_,"AMPLITUDE\n");PR}
			else if (Ctrl->T.flag_detect == MB_DETECT_PHASE)
				{print(stderr_,"PHASE\n");PR}
			else
				{print(stderr_,"%d\n", Ctrl->T.flag_detect);PR}
		}
		if (Ctrl->N.active) {
			print(stderr_,"     Only flag if beams ");PR
			if (Ctrl->N.beam_in)
				{print(stderr_,"between");PR}
			else
				{print(stderr_,"outside");PR}
			print(stderr_," beams %d - %d\n", Ctrl->N.min_beam, Ctrl->N.max_beam_no);PR
		}
		else
			{print(stderr_,"     Flag all beams\n");PR}
		print(stderr_,"Output:\n");PR
		if (Ctrl->B.active)
			{print(stderr_,"     Flag unflagged soundings identified as bad:  ON\n");PR}
		else
			{print(stderr_,"     Flag unflagged soundings identified as bad:  OFF\n");PR}
		if (Ctrl->G.active)
			{print(stderr_,"     Unflag flagged soundings identified as good: ON\n");PR}
		else
			{print(stderr_,"     Unflag flagged soundings identified as good: OFF\n");PR}
	}

	/* get format if required */
	if (Ctrl->F.format == 0)
		mb_get_format(verbose,Ctrl->I.file,NULL,&Ctrl->F.format,&error);

	/* determine whether to read one file or a list of files */
	if (Ctrl->F.format < 0)
		Ctrl->read_datalist  = true;

	/* open file list */
	if (Ctrl->read_datalist) {
		if ((status = mb_datalist_open(verbose,&datalist, Ctrl->I.file,look_processed,&error)) != MB_SUCCESS) {
			error = MB_ERROR_OPEN_FAIL;
			print(stderr_,"Unable to open data list file: %s\n", Ctrl->I.file);PR
			print(stderr_,"Program <%s> Terminated\n", program_name);PR
			Return(error);
		}
		if ((status = mb_datalist_read(verbose,datalist, swathfile,dfile,&Ctrl->F.format,&file_weight,&error)) == MB_SUCCESS)
			read_data = true;
		else
			read_data = false;
	}
	/* else copy single filename to be read */
	else {
		strcpy(swathfile, Ctrl->I.file);
		read_data = true;
	}

	/* loop over all files to be read */
	while (read_data) {

		/* check format and get format flags */
		if ((status = mb_format_flags(verbose,&Ctrl->F.format, &variable_beams, &traveltime, &beam_flagging, &error)) != MB_SUCCESS) {
			mb_error(verbose,error,&message);
			print(stderr_,"\nMBIO Error returned from function <mb_format_flags> regarding input format %d:\n%s\n",Ctrl->F.format,message);PR
			print(stderr_,"\nProgram <%s> Terminated\n", program_name);PR
			Return(error);
		}

		/* check for "fast bathymetry" or "fbt" file */
		strcpy(swathfileread, swathfile);
		formatread = Ctrl->F.format;
		if (Ctrl->T.active == MB_NO)
			mb_get_fbt(verbose, swathfileread, &formatread, &error);

		/* initialize reading the input swath sonar file */
		if ((status = mb_read_init(verbose,swathfileread,formatread,pings,lonflip,bounds, btime_i,
		                           etime_i,speedmin,timegap, &mbio_ptr,&btime_d,&etime_d, &beams_bath,
		                           &beams_amp,&pixels_ss,&error)) != MB_SUCCESS) {
			mb_error(verbose,error,&message);
			print(stderr_,"MBIO Error returned from function <mb_read_init>:\n%s\n",message);PR
			print(stderr_,"Multibeam File <%s> not initialized for reading\n",swathfileread);PR
			print(stderr_,"Program <%s> Terminated\n", program_name);PR
			Return(error);
		}

		/* initialize and increment counting variables */
		pings_file = 0;
		beams_file = 0;

		/* give the statistics */
		if (verbose >= 0) {
			print(stderr_,"Processing %s\n",swathfileread);PR
		}

		/* allocate memory for data arrays */
		beamflag = NULL;
		beamflagorg = NULL;
		detect = NULL;
		bath = NULL;
		amp = NULL;
		bathlon = NULL;
		bathlat = NULL;
		ss = NULL;
		sslon = NULL;
		sslat = NULL;
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), &beamflag, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), &detect, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), &bath, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(double), &amp, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), &bathlon, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), &bathlat, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), &ss, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), &sslon, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), &sslat, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), &beamflagorg, &error);

		/* if error initializing memory then quit */
		if (error != MB_ERROR_NO_ERROR) {
			mb_error(verbose,error,&message);
			print(stderr_,"MBIO Error allocating data arrays:\n%s\n",message);PR
			print(stderr_,"Program <%s> Terminated\n", program_name);PR
			Return(error);
		}

		/* update memory for files */
		if (nfile >= nfile_alloc) {
			nfile_alloc += FILEALLOCNUM;
			status = mb_reallocd(verbose, __FILE__, __LINE__, nfile_alloc * sizeof(struct mbareaclean_file_struct),
								(void **)&files, &error);

			/* if error initializing memory then quit */
			if (error != MB_ERROR_NO_ERROR) {
				mb_error(verbose,error,&message);
				print(stderr_,"MBIO Error allocating data arrays:\n%s\n",message);PR
				print(stderr_,"Program <%s> Terminated\n", program_name);PR
				Return(error);
			}
		}

		/* initialize current file */
		strcpy(files[nfile].filelist, swathfile);
		files[nfile].file_format = Ctrl->F.format;
		files[nfile].nping = 0;
		files[nfile].nping_alloc = PINGALLOCNUM;
		files[nfile].nnull = 0;
		files[nfile].nflag = 0;
		files[nfile].ngood = 0;
		files[nfile].nflagged = 0;
		files[nfile].nunflagged = 0;
		files[nfile].ping_time_d = NULL;
		files[nfile].pingmultiplicity = NULL;
		files[nfile].ping_altitude = NULL;
		files[nfile].nsndg = 0;
		files[nfile].nsndg_alloc = SNDGALLOCNUM;
		files[nfile].sndg_countstart = nsndg;
		files[nfile].beams_bath = beams_bath;
		files[nfile].sndg = NULL;
		status = mb_mallocd(verbose,__FILE__,__LINE__, files[nfile].nping_alloc * sizeof(double),
							&(files[nfile].ping_time_d), &error);
		if (status == MB_SUCCESS)
			status = mb_mallocd(verbose,__FILE__,__LINE__, files[nfile].nping_alloc * sizeof(int),
								&(files[nfile].pingmultiplicity), &error);
		if (status == MB_SUCCESS)
			status = mb_mallocd(verbose,__FILE__,__LINE__, files[nfile].nping_alloc * sizeof(double),
								&(files[nfile].ping_altitude), &error);
		if (status == MB_SUCCESS)
			status = mb_mallocd(verbose,__FILE__,__LINE__, files[nfile].nsndg_alloc * sizeof(struct mbareaclean_sndg_struct),
								&(files[nfile].sndg), &error);
		if (error != MB_ERROR_NO_ERROR) {
			mb_error(verbose,error,&message);
			print(stderr_,"MBIO Error allocating data arrays:\n%s\n",message);PR
			print(stderr_,"Program <%s> Terminated\n", program_name);PR
			Return(error);
		}
		nfile++;

		/* now deal with old edit save file */
		if (status == MB_SUCCESS)		/* handle esf edits */
			status = mb_esf_load(verbose, program_name, swathfile, MB_YES, MB_NO, esffile, &esf, &error);

		/* read */
		done = MB_NO;
		files_tot++;
		pings_file = 0;
		beams_file = 0;
		beams_good_org_file = 0;
		beams_flag_org_file = 0;
		beams_null_org_file = 0;
		while (done == MB_NO) {
			if (verbose > 1) {print(stderr_,"\n");PR}

			/* read next record */
			error = MB_ERROR_NO_ERROR;
			status = mb_read(verbose,mbio_ptr,&kind, &pingsread,time_i,&time_d, &navlon,&navlat, &speed,&heading,
							&distance,&altitude,&sonardepth, &beams_bath,&beams_amp,&pixels_ss, beamflag,bath,amp,
							bathlon,bathlat, ss,sslon,sslat, comment,&error);
			if (verbose >= 2) {
				print(stderr_,"\ndbg2  current data status:\n");PR
				print(stderr_,"dbg2    kind:       %d\n",kind);PR
				print(stderr_,"dbg2    status:     %d\n",status);PR
			}
			if (status == MB_SUCCESS && kind == MB_DATA_DATA) {
				for (i=0;i<beams_bath;i++)
					beamflagorg[i] = beamflag[i];
				status = mb_esf_apply(verbose, &esf, time_d, 1, beams_bath, beamflagorg, &error);

				/* get detection */
				if (Ctrl->T.active) {
					status = mb_get_store(verbose,mbio_ptr,&store_ptr,&error);
					detect_status = mb_detects(verbose,mbio_ptr,store_ptr, &kind,&beams_bath,detect,&detect_error);

					if (detect_status != MB_SUCCESS) {
						status = MB_SUCCESS;
						for (i=0;i<beams_bath;i++)
							detect[i] = MB_DETECT_UNKNOWN;
					}
				}

				/* update counters */
				pings_tot++;
				pings_file++;
				for (i=0;i<beams_bath;i++) {
					if (mb_beam_ok(beamflagorg[i])) {
						beams_tot++;
						beams_file++;
						beams_good_org_tot++;
						beams_good_org_file++;
						files[nfile-1].ngood++;
					}
					else if (beamflagorg[i] == MB_FLAG_NULL) {
						beams_null_org_tot++;
						beams_null_org_file++;
						files[nfile-1].nnull++;
					}
					else {
						beams_tot++;
						beams_file++;
						beams_flag_org_tot++;
						beams_flag_org_file++;
						files[nfile-1].nflag++;
					}
				}

				/* allocate memory if necessary */
				if (files[nfile-1].nping >= files[nfile-1].nping_alloc) {
					files[nfile-1].nping_alloc += PINGALLOCNUM;
					status = mb_reallocd(verbose, __FILE__, __LINE__, files[nfile-1].nping_alloc * sizeof(double), (void **)&(files[nfile-1].ping_time_d), &error);
					if (status == MB_SUCCESS)
						status = mb_reallocd(verbose, __FILE__, __LINE__, files[nfile-1].nping_alloc * sizeof(int), (void **)&(files[nfile-1].pingmultiplicity), &error);
					if (status == MB_SUCCESS)
						status = mb_reallocd(verbose, __FILE__, __LINE__, files[nfile-1].nping_alloc * sizeof(double), (void **)&(files[nfile-1].ping_altitude), &error);
					if (error != MB_ERROR_NO_ERROR) {
						mb_error(verbose,error,&message);
						print(stderr_,"MBIO Error allocating data arrays:\n%s\n",message);PR
						print(stderr_,"Program <%s> Terminated\n", program_name);PR
						Return(error);
					}
				}

				/* store the ping data */
				files[nfile-1].ping_time_d[files[nfile-1].nping] = time_d;
				if (files[nfile-1].nping > 0 && files[nfile-1].ping_time_d[files[nfile-1].nping] ==
					files[nfile-1].ping_time_d[files[nfile-1].nping - 1]) {
					files[nfile-1].pingmultiplicity[files[nfile-1].nping] =
					files[nfile-1].pingmultiplicity[files[nfile-1].nping - 1] + 1;
				}
				else
					files[nfile-1].pingmultiplicity[files[nfile-1].nping] = 0;

				files[nfile-1].ping_altitude[files[nfile-1].nping] = altitude;
				files[nfile-1].nping++;

				/* check beam range */
				if (Ctrl->N.active && Ctrl->N.max_beam_no == 0)
					Ctrl->N.max_beam = beams_bath - Ctrl->N.min_beam;

				/* now loop over the beams and store the soundings in the grid bins */
				for (ib=0;ib<beams_bath;ib++) {
					if (beamflagorg[ib] != MB_FLAG_NULL) { /* get bin for current beam */
						ix = (int)((bathlon[ib] - Ctrl->R.areabounds[0] - 0.5 * dx) / dx);
						iy = (int)((bathlat[ib] - Ctrl->R.areabounds[2] - 0.5 * dy) / dy);
						kgrid = ix*ny + iy;

						/* add sounding */
						if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
							if (files[nfile-1].nsndg >= files[nfile-1].nsndg_alloc) {
								files[nfile-1].nsndg_alloc += SNDGALLOCNUM;
								status = mb_reallocd(verbose, __FILE__, __LINE__, files[nfile-1].nsndg_alloc * sizeof(struct mbareaclean_sndg_struct), (void **)&files[nfile-1].sndg, &error);
								if (error != MB_ERROR_NO_ERROR) {
									mb_error(verbose,error,&message);
									print(stderr_,"MBIO Error allocating sounding arrays:\n%s\n",message);PR
									print(stderr_,"Program <%s> Terminated\n", program_name);PR
									Return(error);
								}
							}
							/* allocate space for sounding if needed */
							if (gsndgnum[kgrid] >= gsndgnum_alloc[kgrid]) {
								gsndgnum_alloc[kgrid] += SNDGALLOCNUM;
								status = mb_reallocd(verbose, __FILE__, __LINE__, gsndgnum_alloc[kgrid] * sizeof(int), (void **)&gsndg[kgrid], &error);
								if (error != MB_ERROR_NO_ERROR) {
									mb_error(verbose,error,&message);
									print(stderr_,"MBIO Error allocating sounding arrays:\n%s\n",message);PR
									print(stderr_,"Program <%s> Terminated\n", program_name);PR
									Return(error);
								}
							}

							/* store sounding data */
							sndg = &(files[nfile-1].sndg[files[nfile-1].nsndg]);
							sndg->sndg_file = nfile - 1;
							sndg->sndg_ping = files[nfile - 1].nping - 1;
							sndg->sndg_beam = ib;
							sndg->sndg_depth = bath[ib];
							sndg->sndg_x = bathlon[ib];
							sndg->sndg_y = bathlat[ib];
							sndg->sndg_beamflag_org = beamflag[ib];
							sndg->sndg_beamflag_esf = beamflagorg[ib];
							sndg->sndg_beamflag = beamflagorg[ib];
							sndg->sndg_edit = MB_YES;
							if (Ctrl->T.active && detect[ib] != Ctrl->T.flag_detect)
								sndg->sndg_edit = MB_NO;
							if (Ctrl->N.active) {
								if (Ctrl->N.min_beam <= ib && ib <= Ctrl->N.max_beam) {
									if (!Ctrl->N.beam_in)
										sndg->sndg_edit = MB_NO;
								}
								else {
									if (Ctrl->N.beam_in)
										sndg->sndg_edit = MB_NO;
								}
							}
							files[nfile-1].nsndg++;
							nsndg++;
							gsndg[kgrid][gsndgnum[kgrid]] = files[nfile-1].sndg_countstart + files[nfile-1].nsndg - 1;
							gsndgnum[kgrid]++;
						}
					}
				}

			}
			else if (error > MB_ERROR_NO_ERROR)
				done = MB_YES;

			/* process a record */

			/* reset counters and data */
		}

		/* close the files */
		status = mb_close(verbose,&mbio_ptr,&error);
		mb_esf_close(verbose, &esf, &error);

		/* check memory */
		if (verbose >= 4)
			status = mb_memory_list(verbose,&error);

		/* give the statistics */
		if (verbose >= 0) {
			print(stderr_,"pings:%4d  beams: %7d good %7d flagged %7d null \n",
			               pings_file,beams_good_org_file ,beams_flag_org_file,beams_null_org_file);PR
		}

		/* figure out whether and what to read next */
		if (Ctrl->read_datalist) {
			if ((status = mb_datalist_read(verbose,datalist, swathfile,dfile,&Ctrl->F.format,&file_weight,&error)) == MB_SUCCESS)
				read_data = true;
			else
				read_data = false;
		}
		else
			read_data = false;

	}	/* end loop over files in list */

	if (Ctrl->read_datalist)
		mb_datalist_close(verbose,&datalist,&error);

	/* loop over grid cells to find maximum number of soundings */
	binnummax = 0;
	for (ix=0;ix<nx;ix++) {
		for (iy=0;iy<ny;iy++) {
			/* get cell id */
			kgrid = ix*ny + iy;
			xx = Ctrl->R.areabounds[0] + 0.5 * dx + ix * dx;
			yy = Ctrl->R.areabounds[3] + 0.5 * dy + iy * dy;
			binnummax = MAX(binnummax, gsndgnum[kgrid]);
		}
		status = mb_mallocd(verbose,__FILE__,__LINE__, binnummax * sizeof(double), &(bindepths), &error);
		if (error != MB_ERROR_NO_ERROR) {
			mb_error(verbose,error,&message);
			print(stderr_,"MBIO Error allocating sounding sorting array:\n%s\n",message);PR
			print(stderr_,"Program <%s> Terminated\n", program_name);PR
			Return(error);
		}
	}

	/* deal with median filter */
	if (Ctrl->M.active) {
	/* loop over grid cells applying median filter test */
	for (ix=0;ix<nx;ix++) {
		for (iy=0;iy<ny;iy++) {
			/* get cell id */
			kgrid = ix*ny + iy;
			xx = Ctrl->R.areabounds[0] + 0.5 * dx + ix * dx;
			yy = Ctrl->R.areabounds[3] + 0.5 * dy + iy * dy;

			/* load up array */
			binnum = 0;
			for (i=0;i<gsndgnum[kgrid];i++) {
				getsoundingptr(API, verbose, gsndg[kgrid][i], nfile, &sndg, &error);
				if (mb_beam_ok(sndg->sndg_beamflag)) {
					bindepths[binnum] = sndg->sndg_depth;
					binnum++;
				}
			}

			/* apply median filter only if there are enough soundings */
			if (binnum >= Ctrl->M.median_filter_nmin) {
				/* run qsort */
				qsort((char *)bindepths,binnum,sizeof(double), (void *)mb_double_compare);
				median_depth = bindepths[binnum / 2];
				if (Ctrl->M.density_filter && binnum / 2 - Ctrl->M.density_filter_nmax / 2 >= 0)
					median_depth_low = bindepths[binnum / 2 + Ctrl->M.density_filter_nmax / 2];
				else
					median_depth_low = bindepths[0];
				if (Ctrl->M.density_filter && binnum / 2 + Ctrl->M.density_filter_nmax / 2 < binnum)
					median_depth_high = bindepths[binnum / 2 + Ctrl->M.density_filter_nmax / 2];
				else
					median_depth_high = bindepths[binnum-1];

				/* process the soundings */
				for (i=0;i<gsndgnum[kgrid];i++) {
					getsoundingptr(API, verbose, gsndg[kgrid][i], nfile, &sndg, &error);
					threshold = fabs(Ctrl->M.median_filter_threshold * files[sndg->sndg_file].ping_altitude[sndg->sndg_ping]);
					flagsounding = MB_NO;
					if (fabs(sndg->sndg_depth - median_depth) > threshold)
						flagsounding = MB_YES;
					if (Ctrl->M.density_filter && (sndg->sndg_depth > median_depth_high || sndg->sndg_depth < median_depth_low))
						flagsounding = MB_YES;
					flag_sounding(API, verbose, flagsounding, Ctrl->B.active, Ctrl->G.active, sndg, &error);
					}
				}
			}
		}
	}

	/* deal with standard deviation filter */
	if (Ctrl->D.active) {
		/* loop over grid cells applying std dev filter test */
		for (ix=0;ix<nx;ix++) {
			for (iy=0;iy<ny;iy++) {		/* get cell id */
				kgrid = ix*ny + iy;
				xx = Ctrl->R.areabounds[0] + 0.5 * dx + ix * dx;
				yy = Ctrl->R.areabounds[3] + 0.5 * dy + iy * dy;

				/* get mean */
				mean = 0.0;
				binnum = 0;
				for (i=0;i<gsndgnum[kgrid];i++) {
					getsoundingptr(API, verbose, gsndg[kgrid][i], nfile, &sndg, &error);
					if (mb_beam_ok(sndg->sndg_beamflag)) {
						mean += sndg->sndg_depth;
						binnum++;
					}
				}
				mean /= binnum;

				/* get standard deviation */
				std_dev = 0.0;
				for (i=0;i<gsndgnum[kgrid];i++) {
					getsoundingptr(API, verbose, gsndg[kgrid][i], nfile, &sndg, &error);
					if (mb_beam_ok(sndg->sndg_beamflag))
						std_dev += (sndg->sndg_depth - mean) * (sndg->sndg_depth - mean);
				}
				std_dev = sqrt(std_dev / binnum);
				threshold = std_dev * Ctrl->D.std_dev_threshold;

				if (binnum > 0 && verbose > 0) {
					print(stderr_,"bin: %d %d %d  pos: %f %f  nsoundings:%d / %d mean:%f std_dev:%f\n",
					              ix,iy,kgrid,xx,yy,binnum,gsndgnum[kgrid],mean,std_dev);PR
				}


				/* apply standard deviation threshold only if there are enough soundings */
				if (binnum >= Ctrl->D.std_dev_nmin) {		/* process the soundings */
					for (i=0;i<gsndgnum[kgrid];i++) {
						getsoundingptr(API, verbose, gsndg[kgrid][i], nfile, &sndg, &error);
						flag_sounding(API, verbose, fabs(sndg->sndg_depth - mean) > threshold, Ctrl->B.active,
						              Ctrl->G.active, sndg, &error);
					}
				}
			}
		}	
	}

	/* loop over files checking for changed soundings */
	for (i=0; i < nfile; i++) {		/* open esf file */
		status = mb_esf_load(verbose, program_name, files[i].filelist, MB_NO, MB_YES, esffile, &esf, &error);
		if (status == MB_SUCCESS && esf.esffp != NULL)
			esffile_open = MB_YES;
		if (status == MB_FAILURE && error == MB_ERROR_OPEN_FAIL) {
			esffile_open = MB_NO;
			print(stderr_, "\nUnable to open new edit save file %s\n", esf.esffile);PR
		}

		/* loop over all of the soundings */
		for (j=0;j<files[i].nsndg;j++) {
			sndg = &(files[i].sndg[j]);
			if (sndg->sndg_beamflag != sndg->sndg_beamflag_org) {
				if (mb_beam_ok(sndg->sndg_beamflag))
					action = MBP_EDIT_UNFLAG;
				else if (mb_beam_check_flag_manual(sndg->sndg_beamflag))
					action = MBP_EDIT_FLAG;
				else if (mb_beam_check_flag_filter(sndg->sndg_beamflag))
					action = MBP_EDIT_FILTER;
				mb_esf_save(verbose, &esf, files[i].ping_time_d[sndg->sndg_ping], sndg->sndg_beam +
							files[i].pingmultiplicity[sndg->sndg_ping] * MB_ESF_MULTIPLICITY_FACTOR, action, &error);
			}
		}

		/* close esf file */
		mb_esf_close(verbose, &esf, &error);

		/* update mbprocess parameter file */
		if (esffile_open == MB_YES) {
			/* update mbprocess parameter file */
			status = mb_pr_update_format(verbose, files[i].filelist, MB_YES, files[i].file_format, &error);
			status = mb_pr_update_edit(verbose, files[i].filelist, MBP_EDIT_ON, esffile, &error);
		}
	}

	/* give the total statistics */
	if (verbose >= 0) {
		print(stderr_,"\nMBareaclean Processing Totals:\n");PR
		print(stderr_,"-------------------------\n");PR
		print(stderr_,"%d total swath data files processed\n",files_tot);PR
		print(stderr_,"%d total pings processed\n",pings_tot);PR
		print(stderr_,"%d total soundings processed\n",beams_tot);PR
		print(stderr_,"-------------------------\n");PR
		for (i=0;i<nfile;i++) {
			print(stderr_,"%7d flagged:%7d unflagged:%7d  file:%s\n",
			              files[i].ngood + files[i].nflag, files[i].nflagged, files[i].nunflagged, files[i].filelist);PR
		}
	}

	/* free arrays */
	mb_freed(verbose,__FILE__, __LINE__, &bindepths,&error);
	for (i=0;i<nx*ny;i++)
		if (gsndg[i] != NULL)
			mb_freed(verbose,__FILE__, __LINE__, &gsndg[i],&error);
	mb_freed(verbose,__FILE__, __LINE__, (void **)&gsndg,&error);
	mb_freed(verbose,__FILE__, __LINE__, (void **)&gsndgnum,&error);
	mb_freed(verbose,__FILE__, __LINE__, (void **)&gsndgnum_alloc,&error);

	for (i=0;i<nfile;i++) {
		mb_freed(verbose,__FILE__, __LINE__, &(files[nfile-1].ping_time_d),&error);
		mb_freed(verbose,__FILE__, __LINE__, &(files[nfile-1].pingmultiplicity),&error);
		mb_freed(verbose,__FILE__, __LINE__, &(files[nfile-1].ping_altitude),&error);
	}
	mb_freed(verbose,__FILE__,__LINE__,(void **)&files,&error);
	files = NULL;

	/* check memory */
	if (verbose >= 4)
		status = mb_memory_list(verbose,&error);

	/* print output debug statements */
	if (verbose >= 2) {
		print(stderr_,"\ndbg2  Program <%s> completed\n", program_name);PR
		print(stderr_,"dbg2  Ending status:\n");PR
		print(stderr_,"dbg2       status:  %d\n",status);PR
	}

	Return (EXIT_SUCCESS);
}

/*--------------------------------------------------------------------*/
int getsoundingptr(struct GMTAPI_CTRL *API, int verbose, int soundingid, int nfile, struct mbareaclean_sndg_struct **sndgptr, int *error) {
	/* local variables */
	char *function_name = "getsoundingptr";
	char  output[MB_PATH_MAXLINE];
	int   i, j;

	/* print input debug statements */
	if (verbose >= 2) {
		print(stderr_,"\ndbg2  MBIO function <%s> called\n", function_name);PR
		print(stderr_,"dbg2  Input arguments:\n");PR
		print(stderr_,"dbg2       verbose:         %d\n",verbose);PR
		print(stderr_,"dbg2       soundingid:      %d\n",soundingid);PR
		print(stderr_,"dbg2       sndgptr:         %p\n",(void *)sndgptr);PR
	}

	/* loop over the files until the sounding is found */
	*sndgptr = NULL;
	for (i=0; i < nfile && *sndgptr == NULL; i++) {
		if (soundingid >= files[i].sndg_countstart && soundingid < files[i].sndg_countstart + files[i].nsndg) {
			j = soundingid - files[i].sndg_countstart;
			*sndgptr = &(files[i].sndg[j]);
		}
	}

	/* print output debug statements */
	if (verbose >= 2) {
		print(stderr_,"\ndbg2  MBIO function <%s> completed\n", function_name);PR
		print(stderr_,"dbg2  Return values:\n");PR
		print(stderr_,"dbg2       *sndgptr:        %p\n",(void *)sndgptr);PR
		print(stderr_,"dbg2       error:           %d\n",*error);PR
	}

	return MB_SUCCESS;
}
/*--------------------------------------------------------------------*/

int flag_sounding(struct GMTAPI_CTRL *API, int verbose, int flag, int output_bad, int output_good,
	              struct mbareaclean_sndg_struct *sndg, int *error) {
	/* local variables */
	char *function_name = "flag_sounding";
	char  output[MB_PATH_MAXLINE];

	/* print input debug statements */
	if (verbose >= 2) {
		print(stderr_,"\ndbg2  MBIO function <%s> called\n", function_name);PR
		print(stderr_,"dbg2  Input arguments:\n");PR
		print(stderr_,"dbg2       verbose:       %d\n",verbose);PR
		print(stderr_,"dbg2       flag:          %d\n",flag);PR
		print(stderr_,"dbg2       output_bad:    %d\n",output_bad);PR
		print(stderr_,"dbg2       output_good:   %d\n",output_good);PR
		print(stderr_,"dbg2       sndg->sndg_edit:     %d\n",sndg->sndg_edit);PR
		print(stderr_,"dbg2       sndg->sndg_beam:     %d\n",sndg->sndg_beam);PR
		print(stderr_,"dbg2       sndg->sndg_beamflag: %d\n",sndg->sndg_beamflag);PR
	}

	if (sndg->sndg_edit == MB_YES) {
		if (output_bad == MB_YES && mb_beam_ok(sndg->sndg_beamflag) && flag ) {
			sndg->sndg_beamflag = MB_FLAG_FLAG + MB_FLAG_FILTER;
			files[sndg->sndg_file].nflagged++;
		}
		else if (output_good == MB_YES && !mb_beam_ok(sndg->sndg_beamflag) &&
		         sndg->sndg_beamflag != MB_FLAG_NULL && !flag) {
			sndg->sndg_beamflag = MB_FLAG_NONE;
			files[sndg->sndg_file].nunflagged++;
		}
		else if (output_good == MB_YES && !mb_beam_ok(sndg->sndg_beamflag) &&
		         sndg->sndg_beamflag != MB_FLAG_NULL && flag) {
			sndg->sndg_edit = MB_NO;
		}
	}

	/* print output debug statements */
	if (verbose >= 2) {
		print(stderr_,"\ndbg2  MBIO function <%s> completed\n", function_name);PR
		print(stderr_,"dbg2  Return values:\n");PR
		print(stderr_,"dbg2       sndg->sndg_edit:     %d\n",sndg->sndg_edit);PR
		print(stderr_,"dbg2       sndg->sndg_beamflag: %d\n",sndg->sndg_beamflag);PR
		print(stderr_,"dbg2       error:           %d\n",*error);PR
	}

	return MB_SUCCESS;
}
