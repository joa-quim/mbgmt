/*--------------------------------------------------------------------
 *    The MB-system:	mbinfo.c	2/1/93
 *    $Id: mbinfo.c 2298 2017-04-10 07:57:48Z caress $
 *
 *    Copyright (c) 1993-2016 by
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
 * MBINFO reads a swath sonar data file and outputs
 * some basic statistics.  If pings are averaged (pings > 2)
 * MBINFO estimates the variance for each of the swath
 * bathymetry beams by reading a set number of pings (>2) and then finding
 * the variance of the detrended values for each beam. The variances
 * for the amplitude beams and sidescan values are
 * calculated without detrending.
 * The results are dumped to stdout.
 *
 * Author:	D. W. Caress
 * Date:	February 1, 1993
 *
 */

#define THIS_MODULE_NAME	"mbinfo"
#define THIS_MODULE_LIB		"mbgmt"
#define THIS_MODULE_PURPOSE	"Read a swath sonar data file and outputs some basic statistics."
#define THIS_MODULE_KEYS	"<D{,>T},>DC" 

/* GMT5 header file */
#include "gmt_dev.h"

EXTERN_MSC int GMT_mbinfo(void *API, int mode, void *args);

#define GMT_PROG_OPTIONS "->RUV"

/* MBIO include files */
#include "mb_status.h"
#include "mb_define.h"
#include "mb_io.h"

#define MBINFO_MAXPINGS 50
GMT_LOCAL struct ping {
	char	*beamflag;
	double	*bath;
	double	*bathlon;
	double	*bathlat;
	double	*amp;
	double	*ss;
	double	*sslon;
	double	*sslat;
};

/* output formats */
#define FREE_TEXT 	0
#define JSON 		1
#define XML			2
#define MAX_OUTPUT_FORMAT 2

#define PR GMT_Put_Record(API, GMT_WRITE_TEXT, output); 
#define print sprintf

/* Control structure for mbgetdata */
GMT_LOCAL struct MBINFO_CTRL {
	struct mbinfo_I {	/* -I<inputfile> */
		bool active;
		char *inputfile;
	} I;
};

GMT_LOCAL void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct  MBINFO_CTRL *Ctrl;
	int     verbose = 0;

	Ctrl = gmt_M_memory (GMT, NULL, 1, struct MBINFO_CTRL);

	Ctrl->I.inputfile = NULL;
	
	return (Ctrl);
}

GMT_LOCAL void Free_Ctrl(struct GMT_CTRL *GMT, struct MBINFO_CTRL *Ctrl) {	/* Deallocate control structure */
	if (!Ctrl) return;
	if (Ctrl->I.inputfile) free (Ctrl->I.inputfile);
	gmt_M_free (GMT, Ctrl);
}

GMT_LOCAL int usage(struct GMTAPI_CTRL *API, int level) {
	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message (API, GMT_TIME_NONE, "usage: mbinfo -I<inputfile>\n");

	return (EXIT_FAILURE);
}

GMT_LOCAL int parse (struct GMT_CTRL *GMT, struct MBINFO_CTRL *Ctrl, struct GMT_OPTION *options) {
	/* This parses the options provided to mbswath and sets parameters in Ctrl.
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
				Ctrl->I.active = true;
				if (gmt_check_filearg (GMT, '<', opt->arg, GMT_IN, GMT_IS_DATASET)) {
					Ctrl->I.inputfile = strdup (opt->arg);
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error: only one input file is allowed.\n");
					n_errors++;
				}
				break;

			/* Processes program-specific parameters */
			case 'I':	/* -I<inputfile> */
				Ctrl->I.active = true;
				if (!gmt_access (GMT, opt->arg, R_OK)) {	/* Got a file */
					Ctrl->I.inputfile = strdup (opt->arg);
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -I: Requires a valid file\n");
					n_errors++;
				}
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, Ctrl->I.active && !Ctrl->I.inputfile,
	                                   "Syntax error -I option: Must specify input file\n");
	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}


static char rcs_id[] = "$Id: mbinfo.c 2298 2017-04-10 07:57:48Z caress $";

/*--------------------------------------------------------------------*/

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

int GMT_mbinfo(void *V_API, int mode, void *args) {
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;	/* General GMT interal parameters */
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */
	struct MBINFO_CTRL *Ctrl = NULL;
	char program_name[] = "MBINFO";
	char help_message[] =  "MBINFO reads a swath sonar data file and outputs\n"
		"some basic statistics.  If pings are averaged (pings > 2)\n"
		"MBINFO estimates the variance for each of the swath\n"
		"beams by reading a set number of pings (>2) and then finding\n"
		"the variance of the detrended values for each beam.\n"
		"The results are dumped to stdout.";
	char usage_message[] = "mbinfo [-Byr/mo/da/hr/mn/sc -C "
		"-Eyr/mo/da/hr/mn/sc -Fformat -G -Ifile -Llonflip -Mnx/ny "
		"-N -O -Ppings -Rw/e/s/n -Sspeed -W -V -H -XinfFormat]";
	extern char *optarg;
	int	errflg = 0;
	int	c;
	int	help = 0;
	int	flag = 0;

	/* MBIO status variables */
	int	status = MB_SUCCESS;
	int	verbose = 0;
	int	error = MB_ERROR_NO_ERROR;
	char	*message;
	char	format_description[MB_DESCRIPTION_LENGTH];

	/* MBIO read control parameters */
	int	read_datalist = MB_NO;
	char	read_file[MB_PATH_MAXLINE];
	void	*datalist;
	int	look_processed = MB_DATALIST_LOOK_UNSET;
	double	file_weight;
	int	format, pings, lonflip;
	double	bounds[4];
	int	btime_i[7], etime_i[7];
	double	btime_d, etime_d, speedmin, timegap;
	char	file[MB_PATH_MAXLINE], dfile[MB_PATH_MAXLINE];
	int	pings_get = 1;
	int	pings_read = 1;
	int	beams_bath_alloc = 0;
	int	beams_amp_alloc = 0;
	int	pixels_ss_alloc = 0;
	int	beams_bath_max = 0;
	int	beams_amp_max = 0;
	int	pixels_ss_max = 0;
	int	beams_bath = 0;
	int	beams_amp = 0;
	int	pixels_ss = 0;

	/* MBIO read values */
	void	*mbio_ptr = NULL;
	struct mb_io_struct *mb_io_ptr;
	int	kind;
	struct ping *data[MBINFO_MAXPINGS];
	struct ping *datacur;
	int	time_i[7];
	double	time_d, navlon, navlat, speed, heading, distance, altitude, sonardepth;
	char	*beamflag = NULL;
	double	*bath = NULL;
	double	*bathlon = NULL;
	double	*bathlat = NULL;
	double	*amp = NULL;
	double	*ss = NULL;
	double	*sslon = NULL;
	double	*sslat = NULL;
	char	comment[MB_COMMENT_MAXLINE];
	int	icomment = 0;

	/* metadata controls */
	int	imetadata = 0;
	int	meta_vessel = 0;
	int	meta_institution = 0;
	int	meta_platform = 0;
	int	meta_sonar = 0;
	int	meta_sonarversion = 0;
	int	meta_cruiseid = 0;
	int	meta_cruisename = 0;
	int	meta_pi = 0;
	int	meta_piinstitution = 0;
	int	meta_client = 0;
	int	meta_svcorrected = 0;
	int	meta_tidecorrected = 0;
	int	meta_batheditmanual = 0;
	int	meta_batheditauto = 0;
	int	meta_rollbias = 0;
	int	meta_pitchbias = 0;
	int	meta_headingbias = 0;
	int	meta_draft = 0;

	/* mbinfo control parameters */
	int	comments = MB_NO;
	int	good_nav_only = MB_NO;
	int	good_nav;
	double	speed_threshold = 50.0;
	int	bathy_in_feet = MB_NO;
	double	bathy_scale;
	int	lonflip_use = 0;
	int	lonflip_set = MB_NO;

	/* limit variables */
	double	lonmin = 0.0;
	double	lonmax = 0.0;
	double	latmin = 0.0;
	double	latmax = 0.0;
	double	sdpmin = 0.0;
	double	sdpmax = 0.0;
	double	altmin = 0.0;
	double	altmax = 0.0;
	double	bathmin = 0.0;
	double	bathmax = 0.0;
	double	ampmin = 0.0;
	double	ampmax = 0.0;
	double	ssmin = 0.0;
	double	ssmax = 0.0;
	double	bathbeg = 0.0;
	double	bathend = 0.0;
	double	lonbeg = 0.0;
	double	latbeg = 0.0;
	double	lonend = 0.0;
	double	latend = 0.0;
	double	spdbeg = 0.0;
	double	hdgbeg = 0.0;
	double	sdpbeg = 0.0;
	double	altbeg = 0.0;
	double	spdend = 0.0;
	double	hdgend = 0.0;
	double	sdpend = 0.0;
	double	altend = 0.0;
	double	timbeg = 0.0;
	double	timend = 0.0;
	int	timbeg_i[7] ={0}, timend_i[7] ={0}, timbeg_j[5] ={0}, timend_j[5] ={0};
	double	distot = 0.0;
	double	timtot = 0.0;
	double	spdavg = 0.0;
	int	irec = 0;
	int	isbtmrec = 0;
	double	timbegfile = 0.0, timendfile = 0.0, distotfile = 0.0, timtotfile = 0.0, spdavgfile = 0.0;
	int	irecfile = 0;
	int	ntdbeams = 0;
	int	ngdbeams = 0;
	int	nzdbeams = 0;
	int	nfdbeams = 0;
	int	ntabeams = 0;
	int	ngabeams = 0;
	int	nzabeams = 0;
	int	nfabeams = 0;
	int	ntsbeams = 0;
	int	ngsbeams = 0;
	int	nzsbeams = 0;
	int	nfsbeams = 0;
	double	ngd_percent, nzd_percent, nfd_percent, nga_percent, nza_percent, nfa_percent, ngs_percent, nzs_percent, nfs_percent;
	int	beginnav = MB_NO;
	int	beginsdp = MB_NO;
	int	beginalt = MB_NO;
	int	beginbath = MB_NO;
	int	beginamp = MB_NO;
	int	beginss = MB_NO;
	int	nread = 0;

	/* variance finding variables */
	int	nbath, namp, nss;
	double	sumx, sumxx, sumy, sumxy, delta;
	double	a, b, dev, mean, variance;
	double	*bathmean = NULL, *bathvar = NULL;
	int	*nbathvar = NULL;
	double	*ampmean = NULL;
	double	*ampvar = NULL;
	int	*nampvar = NULL;
	double	*ssmean = NULL;
	double	*ssvar = NULL;
	int	*nssvar = NULL;
	int	nbathtot_alloc = 0;
	int	namptot_alloc = 0;
	int	nsstot_alloc = 0;
	double	*bathmeantot = NULL, *bathvartot = NULL;
	int	*nbathvartot = NULL;
	double	*ampmeantot = NULL, *ampvartot = NULL;
	int	*nampvartot = NULL;
	double	*ssmeantot = NULL, *ssvartot = NULL;
	int	*nssvartot = NULL;

	/* coverage mask variables */
	int	coverage_mask = MB_NO;
	int	pass, done;
	int	mask_nx = 0, mask_ny = 0;
	double	mask_dx = 0.0, mask_dy = 0.0;
	int	*mask = NULL;

	/* notice variables */
	int	print_notices = MB_NO;
	int	notice_list[MB_NOTICE_MAX] = {0};
	int	notice_list_tot[MB_NOTICE_MAX] = {0};
	int	notice_total;
	char	*notice_msg;

	/* output stream for basic stuff (stdout if verbose <= 1, output if verbose > 1) */
	FILE	*stream = NULL;
	//FILE	*output = NULL;
	FILE	*output_ = NULL;
	int		output_usefile = MB_NO;
	char	output_file[MB_PATH_MAXLINE];
	char	outputf[MB_PATH_MAXLINE];
	char	*fileprint;
	int output_format = FREE_TEXT;
	int len1,len2;
	char    string[500];

	int	val_int, ix, iy, i, j, k, read_data;
	double	speed_apparent, val_double, sigma;
	double	time_d_last = 0.0;

	char	*getenv();

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage(API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE) bailout(usage (API, GMT_USAGE));	/* Return the usage message */
	if (options->option == GMT_OPT_SYNOPSIS) bailout(usage (API, GMT_SYNOPSIS));	/* Return the synopsis */

	/* Parse the command-line arguments */

	GMT = gmt_begin_module(API, THIS_MODULE_LIB, THIS_MODULE_NAME, &GMT_cpy); /* Save current state */
	if (GMT_Parse_Common(API, GMT_PROG_OPTIONS, options)) Return (API->error);
	Ctrl = (struct MBINFO_CTRL *)New_Ctrl(GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse(GMT, Ctrl, options))) Return (error);

	/*---------------------------- This is the mbswath main code ----------------------------*/

	if (GMT_Init_IO (API, GMT_IS_TEXTSET, GMT_IS_NONE, GMT_IN, GMT_ADD_DEFAULT, 0, options) != GMT_NOERROR) {
		Return (API->error);	/* Establishes data files or stdin */
	}
	if (GMT_Init_IO (API, GMT_IS_TEXTSET, GMT_IS_NONE, GMT_OUT, GMT_ADD_DEFAULT, 0, options) != GMT_NOERROR) {	/* Establishes data output */
		Return (API->error);
	}
	if (GMT_Begin_IO (API, GMT_IS_TEXTSET, GMT_OUT, GMT_HEADER_OFF) != GMT_NOERROR) {
		Return (API->error);
	}

	/* get current default values */
	status = mb_defaults(verbose,&format,&pings_get,&lonflip,bounds, btime_i,etime_i,&speedmin,&timegap);

	/* set default input to stdin */
	strcpy (read_file, "stdin");


	/* set output stream */
	if (verbose <= 1)
		stream = stdout;
	else
		stream = stderr;

	/* if error flagged then print it and exit */
	if (errflg) {
		fprintf(stream,"usage: %s\n", usage_message);
		fprintf(stream,"\nProgram <%s> Terminated\n", program_name);
		error = MB_ERROR_BAD_USAGE;
		return(error);
	}

	/* print starting message */
	if (verbose == 1 || help) {
		fprintf(stream,"\nProgram %s\n",program_name);
		fprintf(stream,"Version %s\n",rcs_id);
		fprintf(stream,"MB-system Version %s\n",MB_VERSION);
	}

	/* if help desired then print it and exit */
	if (help) {
		fprintf(stream,"\n%s\n",help_message);
		fprintf(stream,"\nusage: %s\n", usage_message);
		return(error);
	}

	/* get format if required */
	if (format == 0)
		mb_get_format(verbose,Ctrl->I.inputfile,NULL,&format,&error);

	/* set bathymetry scaling */
	if (bathy_in_feet == MB_YES)
		bathy_scale = 1.0 / 0.3048;
	else
		bathy_scale = 1.0;

	/* determine whether to read one file or a list of files */
	if (format < 0)
		read_datalist = MB_YES;

	/* if reading from datalist then variance calculations
		are disabled */
	if (read_datalist == MB_YES)
		pings_read = 1;

	/* Open out file if requested */
	if (output_usefile == MB_YES) {
		strcpy(output_file, Ctrl->I.inputfile);
		switch (output_format) {
			case FREE_TEXT:
				strcat(output_file, ".inf");
				break;
			case JSON:
				strcat(output_file,"_inf.json");
				break;
			case XML:
				strcat(output_file,"_inf.xml");
				break;
			case '?':
				break;
		}
		if ((output_ = fopen(output_file, "wt")) == NULL)
			output_ = stream;
	}
	else
		output_ = stream;

#define output outputf

	switch (output_format) {
		case FREE_TEXT:
			break;
		case JSON:
			print(output,"{\n");PR
			break;
		case XML:
			print(output,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");PR
			print(output,"<mbinfo>\n");PR
			break;
		case '?':
			break;
	}

	/* read only once unless coverage mask requested */
	pass = 0;
	done = MB_NO;
	while (done == MB_NO) {
	/* open file list */
	if (read_datalist == MB_YES) {
		if ((status = mb_datalist_open(verbose,&datalist, Ctrl->I.inputfile,look_processed,&error)) != MB_SUCCESS) {
			error = MB_ERROR_OPEN_FAIL;
			fprintf(stderr,"\nUnable to open data list file: %s\n", Ctrl->I.inputfile);
			fprintf(stderr,"\nProgram <%s> Terminated\n", program_name);
			return(error);
		}
		if ((status = mb_datalist_read(verbose,datalist, file,dfile,&format,&file_weight,&error)) == MB_SUCCESS)
			read_data = MB_YES;
		else
			read_data = MB_NO;
	}
	/* else copy single filename to be read */
	else {
		strcpy(file, Ctrl->I.inputfile);
		read_data = MB_YES;
	}

	/* loop over all files to be read */
	while (read_data == MB_YES) {

	/* initialize reading the swath file */
	if ((status = mb_read_init(verbose,file,format,pings_get,lonflip,bounds, btime_i,etime_i,speedmin,timegap,
	                           &mbio_ptr,&btime_d,&etime_d, &beams_bath_alloc, &beams_amp_alloc, &pixels_ss_alloc,
	                           &error)) != MB_SUCCESS) {
		mb_error(verbose,error,&message);
		fprintf(stream,"\nMBIO Error returned from function <mb_read_init>:\n%s\n",message);
		fprintf(stream,"\nSwath File <%s> not initialized for reading\n",file);
		fprintf(stream,"\nProgram <%s> Terminated\n", program_name);
		return(error);
	}

	/* allocate memory for data arrays */
	for (i=0;i<pings_read;i++) {
		data[i] = NULL;
		status = mb_mallocd(verbose,__FILE__,__LINE__,pings_read*sizeof(struct ping), (void **)&data[i],&error);
		if (error == MB_ERROR_NO_ERROR) {
			datacur = data[i];
			datacur->beamflag = NULL;
			datacur->bath = NULL;
			datacur->amp = NULL;
			datacur->bathlon = NULL;
			datacur->bathlat = NULL;
			datacur->ss = NULL;
			datacur->sslon = NULL;
			datacur->sslat = NULL;
		}
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), (void **)&datacur->beamflag, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&datacur->bath, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(double), (void **)&datacur->amp, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&datacur->bathlon, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&datacur->bathlat, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&datacur->ss, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&datacur->sslon, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&datacur->sslat, &error);
	}
	if (pings_read > 1 && pass == 0) {
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&bathmean, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), (void **)&bathvar, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(int), (void **)&nbathvar, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(double),(void **) &ampmean, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(double), (void **)&ampvar, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(int), (void **)&nampvar, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&ssmean, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), (void **)&ssvar, &error);
		if (error == MB_ERROR_NO_ERROR)
			status = mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(int), (void **)&nssvar, &error);
	}

	/* if coverage mask requested get cell sizes */
	if (pass == 1 && coverage_mask == MB_YES) {
		if (mask_nx > 1 && mask_ny <= 0) {
			if ((lonmax - lonmin) > (latmax - latmin))
				mask_ny = (int)(mask_nx * (latmax - latmin) / (lonmax - lonmin));
			else {
				mask_ny = mask_nx;
				mask_nx = (int)(mask_ny * (lonmax - lonmin) / (latmax - latmin));
				if (mask_ny < 2)
					mask_ny = 2;
			}
		}
		if (mask_nx < 2)
			mask_nx = 2;
		if (mask_ny < 2)
			mask_ny = 2;
		mask_dx = (lonmax - lonmin) / mask_nx;
		mask_dy = (latmax - latmin) / mask_ny;

		/* allocate mask */
		status = mb_mallocd(verbose,__FILE__,__LINE__,mask_nx*mask_ny*sizeof(int), (void **)&mask,&error);
	}

	/* if error initializing memory then quit */
	if (error != MB_ERROR_NO_ERROR) {
		mb_error(verbose,error,&message);
		fprintf(stream,"MBIO Error allocating data arrays:\n%s\n",message);
		fprintf(stream,"Program <%s> Terminated\n", program_name);
		return(error);
	}

	/* initialize data arrays */
	irecfile = 0;
	distotfile = 0.0;
	timtotfile = 0.0;
	spdavgfile = 0.0;
	if (pass == 0 && pings_read > 1) {
		for (i=0;i<beams_bath_alloc;i++) {
			bathmean[i] = 0.0;
			bathvar[i] = 0.0;
			nbathvar[i] = 0;
		}
		for (i=0;i<beams_amp_alloc;i++) {
			ampmean[i] = 0.0;
			ampvar[i] = 0.0;
			nampvar[i] = 0;
		}
		for (i=0;i<pixels_ss_alloc;i++) {
			ssmean[i] = 0.0;
			ssvar[i] = 0.0;
			nssvar[i] = 0;
		}
	}
	if (pass == 1 && coverage_mask == MB_YES) {
		for (i=0;i<mask_nx*mask_ny;i++)
			mask[i] = MB_NO;
	}

	/* initialize metadata counters */
	meta_vessel = 0;
	meta_institution = 0;
	meta_platform = 0;
	meta_sonar = 0;
	meta_sonarversion = 0;
	meta_cruiseid = 0;
	meta_cruisename = 0;
	meta_pi = 0;
	meta_piinstitution = 0;
	meta_client = 0;
	meta_svcorrected = 0;
	meta_tidecorrected = 0;
	meta_batheditmanual = 0;
	meta_batheditauto = 0;
	meta_rollbias = 0;
	meta_pitchbias = 0;
	meta_headingbias = 0;
	meta_draft = 0;

	/* printf out file and format */
	if (pass == 0) {
		if (strrchr(file, '/') == NULL)
			fileprint = file;
		else
			fileprint = strrchr(file, '/') + 1;
		mb_format_description(verbose,&format,format_description,&error);
		switch (output_format) {
			case FREE_TEXT:
				print(output,"\nSwath Data File:      %s\n",fileprint);PR
				print(output,"MBIO Data Format ID:  %d\n",format);PR
				print(output,"%s",format_description);PR
				break;
			case JSON:
				print(output,"\"file_info\":{\n");PR
				print(output,"\"swath_data_file\":\"%s\",\n",fileprint);PR
				print(output,"\"mbio_data_format_id\":\"%d\",\n",format);PR
				len1=(int)strspn(format_description,"Formatname: ");
				len2=(int)strcspn(&format_description[len1],"\n");
				strncpy(string,&format_description[len1],len2);
				print(output,"\"format_name\": \"%s\",\n",string);PR
				len1+=len2+1;
				len1+=(int)strspn(&format_description[len1],"InformalDescription: ");
				len2=(int)strcspn(&format_description[len1],"\n");
				strncpy(string,&format_description[len1],len2);
				string[len2]='\0';
				print(output,"\"informal_description\": \"%s\",\n",string);PR
				len1+=len2+1;
				len1+=(int)strspn(&format_description[len1],"Attributes: ");
				len2=(int)strlen(format_description);
				format_description[strlen(format_description)-1]='\0';
				for (len2=len1;len2<=strlen(format_description);len2++)
					if (format_description[len2]==10)format_description[len2]=';';
				print(output,"\"attributes\": \"%s\"\n",&format_description[len1]);PR
				print(output,"},\n");PR
				break;
			case XML:
				print(output,"\t<file_info>\n");PR
				print(output,"\t\t<swath_data_file>%s</swath_data_file>\n",fileprint);PR
				print(output,"\t\t<mbio_data_format_id>%d</mbio_data_format_id>\n",format);PR
				len1=(int)strspn(format_description,"Formatname: ");
				len2=(int)strcspn(&format_description[len1],"\n");
				strncpy(string,&format_description[len1],len2);
				print(output,"\t\t<format_name>%s</format_name>\n",string);PR
				len1+=len2+1;
				len1+=(int)strspn(&format_description[len1],"InformalDescription: ");
				len2=(int)strcspn(&format_description[len1],"\n");
				strncpy(string,&format_description[len1],len2);
				print(output,"\t\t<informal_description>%s</informal_description>\n",string);PR
				len1+=len2+1;
				len1+=(int)strspn(&format_description[len1],"Attributes: ");
				len2=(int)strlen(format_description);
				format_description[strlen(format_description)-1]='\0';
				for (len2=len1;len2<=strlen(format_description);len2++)
					if (format_description[len2]==10)format_description[len2]=' ';
				print(output,"\t\t<attributes>%s</attributes>\n",&format_description[len1]);PR
				print(output,"\t</file_info>\n");PR
				break;
			case '?':
				errflg++;
			}
		}

	/* read and process data */
	while (error <= MB_ERROR_NO_ERROR) {
		nread = 0;
		error = MB_ERROR_NO_ERROR;
		while (nread < pings_read && error == MB_ERROR_NO_ERROR) {

			/* read a ping of data */
			datacur = data[nread];
			status = mb_read(verbose,mbio_ptr,&kind,&pings, time_i,&time_d, &navlon,&navlat, &speed,&heading,
			                 &distance,&altitude,&sonardepth, &beams_bath,&beams_amp,&pixels_ss,
			                 datacur->beamflag,datacur->bath,datacur->amp, datacur->bathlon,
			                 datacur->bathlat, datacur->ss,datacur->sslon,datacur->sslat, comment,&error);

			/* use local pointers for convenience - do not set these before the
				mb_read call because registered arrays can be dynamically
				reallocated during mb_read, mb_get, and mb_get_all calls */
			beamflag = datacur->beamflag;
			bath = datacur->bath;
			amp = datacur->amp;
			bathlon = datacur->bathlon;
			bathlat = datacur->bathlat;
			ss = datacur->ss;
			sslon = datacur->sslon;
			sslat = datacur->sslat;

			/* increment counters */
			if (pass == 0 && (error == MB_ERROR_NO_ERROR || error == MB_ERROR_TIME_GAP)) {
				irec++;
				irecfile++;
				nread++;
			}

			/* print comment records */
			if (pass == 0 && error == MB_ERROR_COMMENT && comments == MB_YES) {
				if (strncmp(comment,"META",4) != 0) {
					if (icomment == 0) {
						switch (output_format) {
							case FREE_TEXT:
								print(output,"\nComments in file %s:",file);PR
								icomment++;
								break;
							case '?':
								break;
						}
					}
					switch (output_format) {
						case FREE_TEXT:
							print(output,"  %s",comment);PR
							break;
						case JSON:
							print(output,"\"comment\":\"%s\",",comment);PR
							break;
						case XML:
							print(output,"\t<comment>%s</comment>",comment);PR
							break;
						case '?':
							break;
					}
				}
			}

			/* print metadata */
			if (pass == 0 && error == MB_ERROR_COMMENT && strncmp(comment,"META",4) == 0) {
					switch (output_format) {
						case FREE_TEXT:
							if (imetadata == 0) {
								print(output,"\nMetadata:");PR
								imetadata++;
							}
							if (strncmp(comment, "METAVESSEL:", 11) == 0) {
								if (meta_vessel == 0)
									print(output,"Vessel:                 %s", &comment[11]);PR
								meta_vessel++;
							}
							else if (strncmp(comment, "METAINSTITUTION:", 16) == 0) {
								if (meta_institution == 0)
									print(output,"Institution:            %s", &comment[16]);PR
								meta_institution++;
							}
							else if (strncmp(comment, "METAPLATFORM:", 13) == 0) {
								if (meta_platform == 0)
									print(output,"Platform:               %s", &comment[13]);PR
								meta_platform++;
							}
							else if (strncmp(comment, "METASONARVERSION:", 17) == 0) {
								if (meta_sonarversion == 0)
									print(output,"Sonar Version:          %s", &comment[17]);PR
								meta_sonarversion++;
							}
							else if (strncmp(comment, "METASONAR:", 10) == 0) {
								if (meta_sonar == 0)
									print(output,"Sonar:                  %s", &comment[10]);PR
								meta_sonar++;
							}
							else if (strncmp(comment, "METACRUISEID:", 13) == 0) {
								if (meta_cruiseid == 0)
									print(output,"Cruise ID:              %s", &comment[13]);PR
								meta_cruiseid++;
							}
							else if (strncmp(comment, "METACRUISENAME:", 15) == 0) {
								if (meta_cruisename == 0)
									print(output,"Cruise Name:            %s", &comment[15]);PR
								meta_cruisename++;
							}
							else if (strncmp(comment, "METAPI:", 7) == 0) {
								if (meta_pi == 0)
									print(output,"PI:                     %s", &comment[7]);PR
								meta_pi++;
							}
							else if (strncmp(comment, "METAPIINSTITUTION:", 18) == 0) {
								if (meta_piinstitution == 0)
									print(output,"PI Institution:         %s", &comment[18]);PR
								meta_piinstitution++;
							}
							else if (strncmp(comment, "METACLIENT:", 11) == 0) {
								if (meta_client == 0)
									print(output,"Client:                 %s", &comment[11]);PR
								meta_client++;
							}
							else if (strncmp(comment, "METASVCORRECTED:", 16) == 0) {
								if (meta_svcorrected == 0) {
									sscanf(comment, "METASVCORRECTED:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"Corrected Depths:       YES");PR}
									else
										{print(output,"Corrected Depths:       NO");PR}
								}
								meta_svcorrected++;
							}
							else if (strncmp(comment, "METATIDECORRECTED:", 18) == 0) {
								if (meta_tidecorrected == 0) {
									sscanf(comment, "METATIDECORRECTED:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"Tide Corrected:         YES");PR}
									else
										{print(output,"Tide Corrected:         NO");PR}
								}
								meta_tidecorrected++;
							}
							else if (strncmp(comment, "METABATHEDITMANUAL:", 19) == 0) {
								if (meta_batheditmanual == 0) {
									sscanf(comment, "METABATHEDITMANUAL:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"Depths Manually Edited: YES");PR}
									else
										{print(output,"Depths Manually Edited: NO");PR}
								}
								meta_batheditmanual++;
							}
							else if (strncmp(comment, "METABATHEDITAUTO:", 17) == 0) {
								if (meta_batheditauto == 0) {
									sscanf(comment, "METABATHEDITAUTO:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"Depths Auto-Edited:     YES");PR}
									else
										{print(output,"Depths Auto-Edited:     NO");PR}
								}
								meta_batheditauto++;
							}
							else if (strncmp(comment, "METAROLLBIAS:", 13) == 0) {
								if (meta_rollbias == 0) {
									sscanf(comment, "METAROLLBIAS:%lf", &val_double);
									print(output,"Roll Bias:              %f degrees", val_double);PR
								}
								meta_rollbias++;
							}
							else if (strncmp(comment, "METAPITCHBIAS:", 14) == 0) {
								if (meta_pitchbias == 0) {
									sscanf(comment, "METAPITCHBIAS:%lf", &val_double);
									print(output,"Pitch Bias:             %f degrees", val_double);PR
								}
								meta_pitchbias++;
							}
							else if (strncmp(comment, "METAHEADINGBIAS:", 16) == 0) {
								if (meta_headingbias == 0) {
									sscanf(comment, "METAHEADINGBIAS:%lf", &val_double);
									print(output,"Heading Bias:           %f degrees", val_double);PR
								}
								meta_headingbias++;
							}
							else if (strncmp(comment, "METADRAFT:", 10) == 0) {
								if (meta_draft == 0) {
									sscanf(comment, "METADRAFT:%lf", &val_double);
									print(output,"Draft:                  %f m", val_double);PR
								}
								meta_draft++;
							}
							break;
						case JSON:
							if (strncmp(comment, "METAVESSEL:", 11) == 0) {
								if (meta_vessel == 0)
									{print(output,"\"vessel\":\"%s\",", &comment[11]);PR}
								meta_vessel++;
							}
							else if (strncmp(comment, "METAINSTITUTION:", 16) == 0) {
								if (meta_institution == 0)
									{print(output,"\"institution\":\"%s\",", &comment[16]);PR}
								meta_institution++;
							}
							else if (strncmp(comment, "METAPLATFORM:", 13) == 0) {
								if (meta_platform == 0)
									{print(output,"\"platform\": \"%s \",", &comment[13]);PR}
								meta_platform++;
							}
							else if (strncmp(comment, "METASONARVERSION:", 17) == 0) {
								if (meta_sonarversion == 0)
									{print(output,"\"sonar_version\": \"%s\",", &comment[17]);PR}
								meta_sonarversion++;
							}
							else if (strncmp(comment, "METASONAR:", 10) == 0) {
								if (meta_sonar == 0)
									{print(output,"\"sonar\": \"%s\",", &comment[10]);PR}
								meta_sonar++;
							}
							else if (strncmp(comment, "METACRUISEID:", 13) == 0) {
								if (meta_cruiseid == 0)
									{print(output,"\"cruise_id\": \"%s\",", &comment[13]);PR}
								meta_cruiseid++;
							}
							else if (strncmp(comment, "METACRUISENAME:", 15) == 0) {
								if (meta_cruisename == 0)
									{print(output,"\"cruise_name\": \"%s\",", &comment[15]);PR}
								meta_cruisename++;
							}
							else if (strncmp(comment, "METAPI:", 7) == 0) {
								if (meta_pi == 0)
									{print(output,"\"pi\": \"%s\",", &comment[7]);PR}
								meta_pi++;
							}
							else if (strncmp(comment, "METAPIINSTITUTION:", 18) == 0) {
								if (meta_piinstitution == 0)
									{print(output,"\"pi_institution\": \"%s\",", &comment[18]);PR}
								meta_piinstitution++;
							}
							else if (strncmp(comment, "METACLIENT:", 11) == 0) {
								if (meta_client == 0)
									{print(output,"\"client\": \"%s\",", &comment[11]);PR}
								meta_client++;
							}
							else if (strncmp(comment, "METASVCORRECTED:", 16) == 0) {
								if (meta_svcorrected == 0) {
									sscanf(comment, "METASVCORRECTED:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\"corrected_depths\": \"YES\",");PR}
									else
										{print(output,"\"corrected_depths\": \"NO\",");PR}
								}
								meta_svcorrected++;
							}
							else if (strncmp(comment, "METATIDECORRECTED:", 18) == 0) {
								if (meta_tidecorrected == 0) {
									sscanf(comment, "METATIDECORRECTED:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\"tide_corrected\": \"YES\",");PR}
									else
										{print(output,"\"tide_corrected\": \"NO\",");PR}
								}
								meta_tidecorrected++;
							}
							else if (strncmp(comment, "METABATHEDITMANUAL:", 19) == 0) {
								if (meta_batheditmanual == 0) {
									sscanf(comment, "METABATHEDITMANUAL:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\"depths_manually_edited\": \"YES\",");PR}
									else
										{print(output,"\"depths_manually_edited\": \"NO\",");PR}
								}
								meta_batheditmanual++;
							}
							else if (strncmp(comment, "METABATHEDITAUTO:", 17) == 0) {
								if (meta_batheditauto == 0) {
									sscanf(comment, "METABATHEDITAUTO:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\"depths_auto-edited\": \"YES\",");PR}
									else
										{print(output,"\"depths_auto-edited\": \"NO\",");PR}
								}
								meta_batheditauto++;
							}
							else if (strncmp(comment, "METAROLLBIAS:", 13) == 0) {
								if (meta_rollbias == 0) {
									sscanf(comment, "METAROLLBIAS:%lf", &val_double);
									print(output,"\"roll_bias\": \"%f\",", val_double);PR
								}
								meta_rollbias++;
							}
							else if (strncmp(comment, "METAPITCHBIAS:", 14) == 0) {
								if (meta_pitchbias == 0) {
									sscanf(comment, "METAPITCHBIAS:%lf", &val_double);
									print(output,"\"pitch_bias\": \"%f\",", val_double);PR
								}
								meta_pitchbias++;
							}
							else if (strncmp(comment, "METAHEADINGBIAS:", 16) == 0) {
								if (meta_headingbias == 0) {
									sscanf(comment, "METAHEADINGBIAS:%lf", &val_double);
									print(output,"\"heading_bias\": \"%f\",", val_double);PR
								}
								meta_headingbias++;
							}
							else if (strncmp(comment, "METADRAFT:", 10) == 0) {
								if (meta_draft == 0) {
									sscanf(comment, "METADRAFT:%lf", &val_double);
									print(output,"\"draft\": \"%f\",", val_double);PR
								}
								meta_draft++;
							}
							break;
						case XML:
							if (imetadata == 0) {
								print(output,"\t<metadata>");PR
								imetadata++;
							}
							if (strncmp(comment, "METAVESSEL:", 11) == 0) {
								if (meta_vessel == 0)
									{print(output,"\t\t<vessel>%s</vessel>", &comment[11]);PR}
								meta_vessel++;
							}
							else if (strncmp(comment, "METAINSTITUTION:", 16) == 0) {
								if (meta_institution == 0)
									{print(output,"\t\t<institution>%s</institution>", &comment[16]);PR}
								meta_institution++;
							}
							else if (strncmp(comment, "METAPLATFORM:", 13) == 0) {
								if (meta_platform == 0)
									{print(output,"\t\t<platform>%s</platform>", &comment[13]);PR}
								meta_platform++;
							}
							else if (strncmp(comment, "METASONARVERSION:", 17) == 0) {
								if (meta_sonarversion == 0)
									{print(output,"\t\t<sonar_version>%s</sonar_version>", &comment[17]);PR}
								meta_sonarversion++;
							}
							else if (strncmp(comment, "METASONAR:", 10) == 0) {
								if (meta_sonar == 0)
									{print(output,"\t\t<sonar>%s</sonar>", &comment[10]);PR}
								meta_sonar++;
							}
							else if (strncmp(comment, "METACRUISEID:", 13) == 0) {
								if (meta_cruiseid == 0)
									{print(output,"\t\t<cruise_id>%s</cruise_id>", &comment[13]);PR}
								meta_cruiseid++;
							}
							else if (strncmp(comment, "METACRUISENAME:", 15) == 0) {
								if (meta_cruisename == 0)
									{print(output,"\t\t<cruise_name>%s</cruise_name>", &comment[15]);PR}
								meta_cruisename++;
							}
							else if (strncmp(comment, "METAPI:", 7) == 0) {
								if (meta_pi == 0)
									{print(output,"\t\t<pi>%s</pi>", &comment[7]);PR}
								meta_pi++;
							}
							else if (strncmp(comment, "METAPIINSTITUTION:", 18) == 0) {
								if (meta_piinstitution == 0)
									{print(output,"\t\t<pi_institution>%s</pi_institution>", &comment[18]);PR}
								meta_piinstitution++;
							}
							else if (strncmp(comment, "METACLIENT:", 11) == 0) {
								if (meta_client == 0)
									{print(output,"\t\t<client>%s</client>", &comment[11]);PR}
								meta_client++;
							}
							else if (strncmp(comment, "METASVCORRECTED:", 16) == 0) {
								if (meta_svcorrected == 0) {
									sscanf(comment, "METASVCORRECTED:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\t\t<corrected_depths>YES</corrected_depths>");PR}
									else
										{print(output,"\t\t<corrected_depths>NO</corrected_depths>");PR}
								}
								meta_svcorrected++;
							}
							else if (strncmp(comment, "METATIDECORRECTED:", 18) == 0) {
								if (meta_tidecorrected == 0) {
									sscanf(comment, "METATIDECORRECTED:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\t\t<tide_corrected>YES</tide_corrected>");PR}
									else
										{print(output,"\t\t<tide_corrected>NO</tide_corrected>");PR}
								}
								meta_tidecorrected++;
							}
							else if (strncmp(comment, "METABATHEDITMANUAL:", 19) == 0) {
								if (meta_batheditmanual == 0) {
									sscanf(comment, "METABATHEDITMANUAL:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\t\t<depths_manually_edited>YES</depths_manually_edited>");PR}
									else
										{print(output,"\t\t<depths_manually_edited>NO</depths_manually_edited>");PR}
								}
								meta_batheditmanual++;
							}
							else if (strncmp(comment, "METABATHEDITAUTO:", 17) == 0) {
								if (meta_batheditauto == 0) {
									sscanf(comment, "METABATHEDITAUTO:%d", &val_int);
									if (val_int == MB_YES)
										{print(output,"\t\t<depths_auto_edited>YES</depths_auto_edited>");PR}
									else
										{print(output,"\t\t<depths_auto_edited>NO</depths_auto_edited>");PR}
								}
								meta_batheditauto++;
							}
							else if (strncmp(comment, "METAROLLBIAS:", 13) == 0) {
								if (meta_rollbias == 0) {
									sscanf(comment, "METAROLLBIAS:%lf", &val_double);
									print(output,"\t\t<roll_bias>%f</roll_bias>", val_double);PR
								}
								meta_rollbias++;
							}
							else if (strncmp(comment, "METAPITCHBIAS:", 14) == 0) {
								if (meta_pitchbias == 0) {
									sscanf(comment, "METAPITCHBIAS:%lf", &val_double);
									print(output,"\t\t<pitch_bias>%f</pitch_bias>", val_double);PR
								}
								meta_pitchbias++;
							}
							else if (strncmp(comment, "METAHEADINGBIAS:", 16) == 0) {
								if (meta_headingbias == 0) {
									sscanf(comment, "METAHEADINGBIAS:%lf", &val_double);
									print(output,"\t\t<heading_bias>%f</heading_bias>", val_double);PR
								}
								meta_headingbias++;
							}
							else if (strncmp(comment, "METADRAFT:", 10) == 0) {
								if (meta_draft == 0) {
									sscanf(comment, "METADRAFT:%lf", &val_double);
									print(output,"\t\t<draft>%fm</draft>\n\t</metadata>", val_double);PR
								}
								meta_draft++;
							}
							break;
							case '?':
								break;
						}
					}

			/* output error messages */
			if (pass != 0 || error == MB_ERROR_COMMENT) { /* do nothing */
			}
			else if (error == MB_ERROR_SUBBOTTOM) { /* do nothing */
			}
			else if (verbose >= 1 && error < MB_ERROR_NO_ERROR && error >= MB_ERROR_OTHER) {
				mb_error(verbose,error,&message);
				fprintf(stream,"\nNonfatal MBIO Error:\n%s", message);
				fprintf(stream,"Time: %d %d %d %d %d %d %d", time_i[0],time_i[1],time_i[2],
				               time_i[3],time_i[4],time_i[5], time_i[6]);
			}
			else if (verbose >= 1 && error < MB_ERROR_NO_ERROR) {
				mb_error(verbose,error,&message);
				fprintf(stream,"\nNonfatal MBIO Error:\n%s", message);
				fprintf(stream,"Number of good records so far: %d",irecfile);
			}
			else if (verbose >= 1 && error > MB_ERROR_NO_ERROR && error != MB_ERROR_EOF) {
				mb_error(verbose,error,&message);
				fprintf(stream,"\nFatal MBIO Error:\n%s", message);
				fprintf(stream,"Last Good Time: %d %d %d %d %d %d %d", time_i[0],time_i[1],time_i[2],
				               time_i[3],time_i[4],time_i[5], time_i[6]);
			}

			/* take note of min and maxes */
			beams_bath_max = MAX(beams_bath_max, beams_bath);
			beams_amp_max = MAX(beams_amp_max, beams_amp);
			pixels_ss_max = MAX(pixels_ss_max, pixels_ss);
			if (pass == 0 && (error == MB_ERROR_NO_ERROR || error == MB_ERROR_TIME_GAP)) {
				/* update data counts */
				ntdbeams += beams_bath;
				ntabeams += beams_amp;
				ntsbeams += pixels_ss;

				/* set lonflip if needed */
				if (lonflip_set == MB_NO && (navlon != 0.0 || navlat != 0.0)) {
					lonflip_set = MB_YES;
					if (navlon < -270.0)
						lonflip_use = 0;
					else if (navlon >= -270.0 && navlon < -90.0)
						lonflip_use = -1;
					else if (navlon >= -90.0 && navlon < 90.0)
						lonflip_use = 0;
					else if (navlon >= 90.0 && navlon < 270.0)
						lonflip_use = 1;
					else if (navlon >= 270.0)
						lonflip_use = 0;

					/* change and apply lonflip if needed */
					if (lonflip_use != lonflip) {
						/* change lonflip used in reading */
						mb_io_ptr = (struct mb_io_struct *)mbio_ptr;
						mb_io_ptr->lonflip = lonflip_use;
						lonflip = lonflip_use;

						/* apply lonflip to data already read */
						if (lonflip_use == -1) {
							if (navlon > 0.0)
								navlon -= 360.0;
							for (i=0;i<beams_bath;i++) {
								if (bathlon[i] > 0.0) bathlon[i] -= 360.0;
							}
							for (i=0;i<pixels_ss;i++) {
								if (sslon[i] > 0.0) sslon[i] -= 360.0;
							}
						}
						else if (lonflip_use == 1) {
							if (navlon < 0.0)
								navlon += 360.0;
							for (i=0;i<beams_bath;i++) {
								if (bathlon[i] < 0.0) bathlon[i] += 360.0;
							}
							for (i=0;i<pixels_ss;i++) {
								if (sslon[i] < 0.0) sslon[i] += 360.0;
							}
						}
						else if (lonflip_use == 0) {
							if (navlon < -180.0) navlon += 360.0;
							else if (navlon > 180.0) navlon -= 360.0;
							for (i=0;i<beams_bath;i++) {
								if (bathlon[i] < -180.0) bathlon[i] += 360.0;
								else if (bathlon[i] > 180.0) bathlon[i] -= 360.0;
							}
							for (i=0;i<pixels_ss;i++) {
								if (sslon[i] < -180.0) sslon[i] += 360.0;
								else if (sslon[i] > 180.0) sslon[i] -= 360.0;
							}
						}
					}
				}

				/* get beginning values */
				if (irec == 1) {
					if (beams_bath > 0) {
						if (mb_beam_ok(beamflag[beams_bath/2]))
							bathbeg = bath[beams_bath/2];
						else
							bathbeg = altitude + sonardepth;
					}
					lonbeg = navlon;
					latbeg = navlat;
					timbeg = time_d;
					timbegfile = time_d;
					for (i=0;i<7;i++)
						timbeg_i[i] = time_i[i];
					spdbeg = speed;
					hdgbeg = heading;
					sdpbeg = sonardepth;
					altbeg = altitude;
				}
				else if (good_nav_only == MB_YES) {
					if (lonbeg == 0.0 && latbeg == 0.0 && navlon != 0.0 && navlat != 0.0) {
						lonbeg = navlon;
						if (beams_bath > 0) {
							if (mb_beam_ok(beamflag[beams_bath/2]))
								bathbeg = bath[beams_bath/2];
							else
								bathbeg = altitude + sonardepth;
						}
						latbeg = navlat;
						if (spdbeg == 0.0 && speed != 0.0)
							spdbeg = speed;
						if (hdgbeg == 0.0 && heading != 0.0)
							hdgbeg = heading;
						if (sdpbeg == 0.0 && sonardepth != 0.0)
							sdpbeg = sonardepth;
						if (altbeg == 0.0 && altitude != 0.0)
							altbeg = altitude;
					}
				}

				/* reset ending values each time */
				if (beams_bath > 0) {
					if (mb_beam_ok(beamflag[beams_bath/2]))
						bathend = bath[beams_bath/2];
					else
						bathend = altitude + sonardepth;
				}
				lonend = navlon;
				latend = navlat;
				spdend = speed;
				hdgend = heading;
				sdpend = sonardepth;
				altend = altitude;
				timend = time_d;
				timendfile = time_d;
				for (i=0;i<7;i++)
					timend_i[i] = time_i[i];

				/* check for good nav */
				speed_apparent = 3600.0*distance
					/(time_d - time_d_last);
				if (good_nav_only == MB_YES) {
					if (navlon == 0.0 || navlat == 0.0) {
						good_nav = MB_NO;
					}
					else if (beginnav == MB_YES && speed_apparent >= speed_threshold) {
						good_nav = MB_NO;
					}
					else {
						good_nav = MB_YES;
					}
				}
				else
					good_nav = MB_YES;

				/* get total distance */
				if (good_nav_only == MB_NO || (good_nav == MB_YES && speed_apparent < speed_threshold)) {
					distot+= distance;
					distotfile += distance;
				}

				/* get starting mins and maxs */
				if (beginnav == MB_NO && good_nav == MB_YES) {
					lonmin = navlon;
					lonmax = navlon;
					latmin = navlat;
					latmax = navlat;
					beginnav = MB_YES;
				}
				if (beginsdp == MB_NO && sonardepth > 0.0) {
					sdpmin = sonardepth;
					sdpmax = sonardepth;
					beginsdp = MB_YES;
				}
				if (beginalt == MB_NO && altitude > 0.0) {
					altmin = altitude;
					altmax = altitude;
					beginalt = MB_YES;
				}
				if (beginbath == MB_NO && beams_bath > 0)
					for (i=0;i<beams_bath;i++)
						if (mb_beam_ok(beamflag[i])) {
							bathmin = bath[i];
							bathmax = bath[i];
							beginbath = MB_YES;
						}
				if (beginamp == MB_NO && beams_amp > 0)
					for (i=0;i<beams_amp;i++)
						if (mb_beam_ok(beamflag[i])) {
							ampmin = amp[i];
							ampmax = amp[i];
							beginamp = MB_YES;
						}
				if (beginss == MB_NO && pixels_ss > 0)
					for (i=0;i<pixels_ss;i++)
						if (ss[i] > MB_SIDESCAN_NULL) {
							ssmin = ss[i];
							ssmax = ss[i];
							beginss = MB_YES;
						}

				/* get mins and maxs */
				if (good_nav == MB_YES && beginnav == MB_YES) {
					lonmin = MIN(lonmin, navlon);
					lonmax = MAX(lonmax, navlon);
					latmin = MIN(latmin, navlat);
					latmax = MAX(latmax, navlat);
				}
				if (beginsdp == MB_YES) {
					sdpmin = MIN(sdpmin, sonardepth);
					sdpmax = MAX(sdpmax, sonardepth);
				}
				if (beginalt == MB_YES) {
					altmin = MIN(altmin, altitude);
					altmax = MAX(altmax, altitude);
				}
				for (i=0;i<beams_bath;i++) {
					if (mb_beam_ok(beamflag[i])) {
						if (good_nav == MB_YES && beginnav == MB_YES) {
							lonmin = MIN(lonmin, bathlon[i]);
							lonmax = MAX(lonmax, bathlon[i]);
							latmin = MIN(latmin, bathlat[i]);
							latmax = MAX(latmax, bathlat[i]);
						}
						bathmin = MIN(bathmin, bath[i]);
						bathmax = MAX(bathmax, bath[i]);
						ngdbeams++;
					}
					else if (beamflag[i] == MB_FLAG_NULL)
						nzdbeams++;
					else
						nfdbeams++;
				}
				for (i=0;i<beams_amp;i++) {
					if (mb_beam_ok(beamflag[i])) {
						ampmin = MIN(ampmin, amp[i]);
						ampmax = MAX(ampmax, amp[i]);
						ngabeams++;
					}
					else if (beamflag[i] == MB_FLAG_NULL)
						nzabeams++;
					else
						nfabeams++;
				}
				for (i=0;i<pixels_ss;i++) {
					if (ss[i] > MB_SIDESCAN_NULL) {
						if (good_nav == MB_YES && beginnav == MB_YES) {
							lonmin = MIN(lonmin, sslon[i]);
							lonmax = MAX(lonmax, sslon[i]);
							latmin = MIN(latmin, sslat[i]);
							latmax = MAX(latmax, sslat[i]);
						}
						ssmin = MIN(ssmin, ss[i]);
						ssmax = MAX(ssmax, ss[i]);
						ngsbeams++;
					}
					else if (ss[i] == 0.0)
						nzsbeams++;
					else
						nfsbeams++;
				}

				/* reset time of last ping */
				time_d_last = time_d;
			}

			/* update coverage mask */
			if (pass == 1 && coverage_mask == MB_YES && (error == MB_ERROR_NO_ERROR || error == MB_ERROR_TIME_GAP)) {
				ix = (int)((navlon - lonmin) / mask_dx);
				iy = (int)((navlat - latmin) / mask_dy);
				if (ix >= 0 && ix < mask_nx && iy >= 0 && iy < mask_ny) {
					mask[ix+iy*mask_nx] = MB_YES;
				}
				for (i=0;i<beams_bath;i++) {
					if (mb_beam_ok(beamflag[i])) {
						ix = (int)((bathlon[i] - lonmin) / mask_dx);
						iy = (int)((bathlat[i] - latmin) / mask_dy);
						if (ix >= 0 && ix < mask_nx && iy >= 0 && iy < mask_ny)
							mask[ix+iy*mask_nx] = MB_YES;
					}
				}
				for (i=0;i<pixels_ss;i++) {
					if (ss[i] > MB_SIDESCAN_NULL) {
						ix = (int)((sslon[i] - lonmin) / mask_dx);
						iy = (int)((sslat[i] - latmin) / mask_dy);
						if (ix >= 0 && ix < mask_nx && iy >= 0 && iy < mask_ny)
							mask[ix+iy*mask_nx] = MB_YES;
					}
				}
			}

			/* look for problems */
			if (pass == 0 && (error == MB_ERROR_NO_ERROR || error == MB_ERROR_TIME_GAP)) {
				if (navlon == 0.0 || navlat == 0.0)
				mb_notice_log_problem(verbose, mbio_ptr, MB_PROBLEM_ZERO_NAV);
				else if (beginnav == MB_YES && speed_apparent >= speed_threshold)
				mb_notice_log_problem(verbose, mbio_ptr, MB_PROBLEM_TOO_FAST);
				for (i=0;i<beams_bath;i++) {
					if (mb_beam_ok(beamflag[i])) {
						if (bath[i] > 11000.0)
							mb_notice_log_problem(verbose, mbio_ptr, MB_PROBLEM_TOO_DEEP);
					}
				}
			}
		}

		/* process the pings */
		if (pass == 0 && pings_read > 2 && nread == pings_read && (error == MB_ERROR_NO_ERROR || error == MB_ERROR_TIME_GAP)) {

			/* do the bathymetry */
			for (i=0;i<beams_bath;i++) {

				/* fit line to depths */
				nbath  = 0;
				sumx  = 0.0;
				sumxx = 0.0;
				sumy  = 0.0;
				sumxy = 0.0;
				variance = 0.0;
				for (j=0;j<nread;j++) {
					datacur = data[j];
					bath = datacur->bath;
					beamflag = datacur->beamflag;
					if (mb_beam_ok(beamflag[i])) {
						nbath++;
						sumx  = sumx + j;
						sumxx = sumxx + j*j;
						sumy  = sumy + bath[i];
						sumxy = sumxy + j*bath[i];
					}
				}
				if (nbath == pings_read) {
					delta = nbath*sumxx - sumx*sumx;
					a = (sumxx*sumy - sumx*sumxy)/delta;
					b = (nbath*sumxy - sumx*sumy)/delta;
					for (j=0;j<nread;j++) {
						datacur = data[j];
						bath = datacur->bath;
						beamflag = datacur->beamflag;
						if (mb_beam_ok(beamflag[i])) {
							dev = bath[i] - a - b*j;
							variance = variance + dev*dev;
						}
					}
					bathmean[i] = bathmean[i] + sumy;
					bathvar[i] = bathvar[i] + variance;
					nbathvar[i] = nbathvar[i] + nbath;
				}
			}

			/* do the amplitude */
			for (i=0;i<beams_amp;i++) {		/* get mean amplitude */
				namp  = 0;
				mean  = 0.0;
				variance = 0.0;
				for (j=0;j<nread;j++) {
					datacur = data[j];
					amp = datacur->amp;
					beamflag = datacur->beamflag;
					if (mb_beam_ok(beamflag[i])) {
						namp++;
						mean  = mean + amp[i];
					}
				}
				if (namp == pings_read) {
					mean = mean/namp;
					for (j=0;j<nread;j++) {
						datacur = data[j];
						amp = datacur->amp;
						if (mb_beam_ok(beamflag[i])) {
							dev = amp[i] - mean;
							variance = variance + dev*dev;
						}
					}
					ampmean[i] = ampmean[i] + namp*mean;
					ampvar[i] = ampvar[i] + variance;
					nampvar[i] = nampvar[i] + namp;
				}
			}

			/* do the sidescan */
			for (i=0;i<pixels_ss;i++) {		/* get mean sidescan */
				nss  = 0;
				mean  = 0.0;
				variance = 0.0;
				for (j=0;j<nread;j++) {
					datacur = data[j];
					ss = datacur->ss;
					if (ss[i] > MB_SIDESCAN_NULL) {
						nss++;
						mean  = mean + ss[i];
					}
				}
				if (nss == pings_read) {
					mean = mean/nss;
					for (j=0;j<nread;j++) {
						datacur = data[j];
						ss = datacur->ss;
						if (ss[i] > MB_SIDESCAN_NULL) {
							dev = ss[i] - mean;
							variance = variance + dev*dev;
						}
					}
					ssmean[i] = ssmean[i] + nss*mean;
					ssvar[i] = ssvar[i] + variance;
					nssvar[i] = nssvar[i] + nss;
				}
			}
		}
	}

	/* look for problems */
	timtotfile = (timendfile - timbegfile)/3600.0;
	if (timtotfile > 0.0)
		spdavgfile = distotfile/timtotfile;
	if (irecfile <= 0)
		mb_notice_log_problem(verbose, mbio_ptr, MB_PROBLEM_NO_DATA);
	else if (timtotfile > 0.0 && spdavgfile >= speed_threshold)
		mb_notice_log_problem(verbose, mbio_ptr, MB_PROBLEM_AVG_TOO_FAST);

	/* get notices if desired */
	if (print_notices == MB_YES && pass == 0) {
		status = mb_notice_get_list(verbose, mbio_ptr, notice_list);
		for (i=0;i<MB_NOTICE_MAX;i++)
			notice_list_tot[i] += notice_list[i];
	}

	/* deal with statistics */
	if (pings_read > 2) {
		/* allocate total statistics arrays if needed */
		if (nbathtot_alloc < beams_bath_max) {
			status = mb_reallocd(verbose,__FILE__,__LINE__,beams_bath_max*sizeof(double), (void **)&bathmeantot,&error);
			status = mb_reallocd(verbose,__FILE__,__LINE__,beams_bath_max*sizeof(double), (void **)&bathvartot,&error);
			status = mb_reallocd(verbose,__FILE__,__LINE__,beams_bath_max*sizeof(int), (void **)&nbathvartot,&error);
			if (error != MB_ERROR_NO_ERROR) {
				mb_error(verbose,error,&message);
				fprintf(stream,"\nMBIO Error allocating data arrays:\n%s",message);
				fprintf(stream,"\nProgram <%s> Terminated", program_name);
				return(error);
			}
			else {
				for (i=nbathtot_alloc;i<beams_bath_max;i++) {
					bathmeantot[i] = 0.0;
					bathvartot[i] = 0.0;
					nbathvartot[i] = 0;
				}
				nbathtot_alloc = beams_bath_max;
			}
		}
		if (namptot_alloc < beams_amp_max) {
			status = mb_reallocd(verbose,__FILE__,__LINE__,beams_amp_max*sizeof(double), (void **)&ampmeantot,&error);
			status = mb_reallocd(verbose,__FILE__,__LINE__,beams_amp_max*sizeof(double), (void **)&ampvartot,&error);
			status = mb_reallocd(verbose,__FILE__,__LINE__,beams_amp_max*sizeof(int), (void **)&nampvartot,&error);
			if (error != MB_ERROR_NO_ERROR) {
				mb_error(verbose,error,&message);
				fprintf(stream,"\nMBIO Error allocating data arrays:\n%s",message);
				fprintf(stream,"\nProgram <%s> Terminated", program_name);
				return(error);
			}
			else {
				for (i=namptot_alloc;i<beams_amp_max;i++) {
					ampmeantot[i] = 0.0;
					ampvartot[i] = 0.0;
					nampvartot[i] = 0;
				}
				namptot_alloc = beams_amp_max;
			}
		}
		if (nsstot_alloc < pixels_ss_max) {
			status = mb_reallocd(verbose,__FILE__,__LINE__,pixels_ss_max*sizeof(double), (void **)&ssmeantot,&error);
			status = mb_reallocd(verbose,__FILE__,__LINE__,pixels_ss_max*sizeof(double), (void **)&ssvartot,&error);
			status = mb_reallocd(verbose,__FILE__,__LINE__,pixels_ss_max*sizeof(int), (void **)&nssvartot,&error);
			if (error != MB_ERROR_NO_ERROR) {
				mb_error(verbose,error,&message);
				fprintf(stream,"\nMBIO Error allocating data arrays:\n%s",message);
				fprintf(stream,"\nProgram <%s> Terminated", program_name);
				return(error);
			}
			else {
				for (i=nsstot_alloc;i<pixels_ss_max;i++) {
					ssmeantot[i] = 0.0;
					ssvartot[i] = 0.0;
					nssvartot[i] = 0;
				}
				nsstot_alloc = pixels_ss_max;
			}
		}

		/* copy statistics to total statistics */
		for (i=0;i<beams_bath;i++) {
			bathmeantot[i] += bathmean[i];
			bathvartot[i] += bathvar[i];
			nbathvartot[i] += nbathvar[i];
		}
		for (i=0;i<beams_amp;i++) {
			ampmeantot[i] += ampmean[i];
			ampvartot[i] += ampvar[i];
			nampvartot[i] += nampvar[i];
		}
		for (i=0;i<pixels_ss;i++) {
			ssmeantot[i] += ssmean[i];
			ssvartot[i] += ssvar[i];
			nssvartot[i] += nssvar[i];
		}
	}

		/* close the swath file */
		status = mb_close(verbose,&mbio_ptr,&error);

		/* figure out whether and what to read next */
		if (read_datalist == MB_YES) {
			if ((status = mb_datalist_read(verbose,datalist, file,dfile,&format,&file_weight,&error)) == MB_SUCCESS)
				read_data = MB_YES;
			else
				read_data = MB_NO;
		}
		else
			read_data = MB_NO;

	} /* end loop over files in list */

		if (read_datalist == MB_YES)
			mb_datalist_close(verbose,&datalist,&error);

		if (pass > 0 || coverage_mask == MB_NO)		/* figure out if done */
			done = MB_YES;
		pass++;

	} /* end loop over reading passes */

	/* calculate final variances */
	if (pings_read > 2) {
		for (i=0;i<nbathtot_alloc;i++)
			if (nbathvartot[i] > 0) {
				bathmeantot[i] = bathmeantot[i]/nbathvartot[i];
				bathvartot[i] = bathvartot[i]/nbathvartot[i];
			}
		for (i=0;i<namptot_alloc;i++)
			if (nampvartot[i] > 0) {
				ampmeantot[i] = ampmeantot[i]/nampvartot[i];
				ampvartot[i] = ampvartot[i]/nampvartot[i];
			}
		for (i=0;i<nsstot_alloc;i++)
			if (nssvartot[i] > 0) {
				ssmeantot[i] = ssmeantot[i]/nssvartot[i];
				ssvartot[i] = ssvartot[i]/nssvartot[i];
			}
	}

	/* calculate percentages of data */
	if (ntdbeams > 0) {
		ngd_percent = 100.0*ngdbeams/ntdbeams;
		nzd_percent = 100.0*nzdbeams/ntdbeams;
		nfd_percent = 100.0*nfdbeams/ntdbeams;
	}
	else {
		ngd_percent = 0.0;
		nzd_percent = 0.0;
		nfd_percent = 0.0;
	}
	if (ntabeams > 0) {
		nga_percent = 100.0*ngabeams/ntabeams;
		nza_percent = 100.0*nzabeams/ntabeams;
		nfa_percent = 100.0*nfabeams/ntabeams;
	}
	else {
		nga_percent = 0.0;
		nza_percent = 0.0;
		nfa_percent = 0.0;
	}
	if (ntsbeams > 0) {
		ngs_percent = 100.0*ngsbeams/ntsbeams;
		nzs_percent = 100.0*nzsbeams/ntsbeams;
		nfs_percent = 100.0*nfsbeams/ntsbeams;
	}
	else {
		ngs_percent = 0.0;
		nzs_percent = 0.0;
		nfs_percent = 0.0;
	}

	/* now print out the results */
	timtot = (timend - timbeg)/3600.0;
	if (timtot > 0.0)
		spdavg = distot/timtot;
	mb_get_jtime(verbose,timbeg_i,timbeg_j);
	mb_get_jtime(verbose,timend_i,timend_j);

	switch (output_format) {
		case FREE_TEXT:
			print(output,"\nData Totals:");PR
			print(output,"Number of Records:                    %8d",irec);PR
			isbtmrec = notice_list_tot[MB_DATA_SUBBOTTOM_MCS]
					+ notice_list_tot[MB_DATA_SUBBOTTOM_CNTRBEAM]
					+ notice_list_tot[MB_DATA_SUBBOTTOM_SUBBOTTOM];
			if (isbtmrec > 0)
				{print(output,"Number of Subbottom Records:          %8d",isbtmrec);PR}
			if (notice_list_tot[MB_DATA_SIDESCAN2] > 0)
				{print(output,"Number of Secondary Sidescan Records: %8d",notice_list_tot[MB_DATA_SIDESCAN2]);PR}
			if (notice_list_tot[MB_DATA_SIDESCAN3] > 0)
				{print(output,"Number of Tertiary Sidescan Records:  %8d",notice_list_tot[MB_DATA_SIDESCAN3]);PR}
			if (notice_list_tot[MB_DATA_WATER_COLUMN] > 0)
				{print(output,"Number of Water Column Records:       %8d",notice_list_tot[MB_DATA_WATER_COLUMN]);PR}

			print(output,"Bathymetry Data (%d beams):",beams_bath_max);PR
			print(output,"  Number of Beams:         %8d", ntdbeams);PR
			print(output,"  Number of Good Beams:    %8d     %5.2f%%", ngdbeams, ngd_percent);PR
			print(output,"  Number of Zero Beams:    %8d     %5.2f%%", nzdbeams, nzd_percent);PR
			print(output,"  Number of Flagged Beams: %8d     %5.2f%%", nfdbeams, nfd_percent);PR
			print(output,"Amplitude Data (%d beams):",beams_amp_max);PR
			print(output,"  Number of Beams:         %8d", ntabeams);PR
			print(output,"  Number of Good Beams:    %8d     %5.2f%%", ngabeams, nga_percent);PR
			print(output,"  Number of Zero Beams:    %8d     %5.2f%%", nzabeams, nza_percent);PR
			print(output,"  Number of Flagged Beams: %8d     %5.2f%%", nfabeams, nfa_percent);PR
			print(output,"Sidescan Data (%d pixels):",pixels_ss_max);PR
			print(output,"  Number of Pixels:        %8d", ntsbeams);PR
			print(output,"  Number of Good Pixels:   %8d     %5.2f%%", ngsbeams, ngs_percent);PR
			print(output,"  Number of Zero Pixels:   %8d     %5.2f%%", nzsbeams, nzs_percent);PR
			print(output,"  Number of Flagged Pixels:%8d     %5.2f%%", nfsbeams, nfs_percent);PR
			print(output,"\nNavigation Totals:");PR
			print(output,"Total Time:         %10.4f hours",timtot);PR
			print(output,"Total Track Length: %10.4f km",distot);PR
			print(output,"Average Speed:      %10.4f km/hr (%7.4f knots)", spdavg,spdavg/1.85);PR
			print(output,"\nStart of Data:");PR
			print(output,"Time:  %2.2d %2.2d %4.4d %2.2d:%2.2d:%2.2d.%6.6d  JD%d (%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%6.6d)",
			             timbeg_i[1],timbeg_i[2],timbeg_i[0],timbeg_i[3], timbeg_i[4],timbeg_i[5],timbeg_i[6],timbeg_j[1],
			             timbeg_i[0],timbeg_i[1],timbeg_i[2],timbeg_i[3],timbeg_i[4],timbeg_i[5],timbeg_i[6]);PR
			if (bathy_in_feet == MB_NO)
				{print(output,"Lon: %15.9f     Lat: %15.9f     Depth: %10.4f meters", lonbeg,latbeg,bathbeg);PR}
			else
				{print(output,"Lon: %15.9f     Lat: %15.9f     Depth: %10.4f feet", lonbeg,latbeg,bathy_scale*bathbeg);PR}
			print(output,"Speed: %7.4f km/hr (%7.4f knots)  Heading:%9.4f degrees", spdbeg,spdbeg/1.85,hdgbeg);PR
			print(output,"Sonar Depth:%10.4f m  Sonar Altitude:%10.4f m", sdpbeg,altbeg);PR
			print(output,"\nEnd of Data:");PR
			print(output,"Time:  %2.2d %2.2d %4.4d %2.2d:%2.2d:%2.2d.%6.6d  JD%d (%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%6.6d)",
			             timend_i[1],timend_i[2],timend_i[0],timend_i[3], timend_i[4],timend_i[5],timend_i[6],timend_j[1],
			             timend_i[0],timend_i[1],timend_i[2],timend_i[3],timend_i[4],timend_i[5],timend_i[6]);PR
			if (bathy_in_feet == MB_NO)
				{print(output,"Lon: %15.9f     Lat: %15.9f     Depth: %10.4f meters", lonend,latend,bathend);PR}
			else
				{print(output,"Lon: %15.9f     Lat: %15.9f     Depth: %10.4f feet", lonend,latend,bathy_scale*bathend);PR}
			print(output,"Speed: %7.4f km/hr (%7.4f knots)  Heading:%9.4f degrees", spdend,spdend/1.85,hdgend);PR
			print(output,"Sonar Depth:%10.4f m  Sonar Altitude:%10.4f m", sdpend,altend);PR
			print(output,"\nLimits:");PR
			print(output,"Minimum Longitude:   %15.9f   Maximum Longitude:   %15.9f",lonmin,lonmax);PR
			print(output,"Minimum Latitude:    %15.9f   Maximum Latitude:    %15.9f",latmin,latmax);PR
			print(output,"Minimum Sonar Depth: %10.4f   Maximum Sonar Depth: %10.4f",sdpmin,sdpmax);PR
			print(output,"Minimum Altitude:    %10.4f   Maximum Altitude:    %10.4f",altmin,altmax);PR
			if (ngdbeams > 0 || verbose >= 1) {
				print(output,"Minimum Depth:       %10.4f   Maximum Depth:       %10.4f",
				             bathy_scale*bathmin,bathy_scale*bathmax);PR
			}
			if (ngabeams > 0 || verbose >= 1)
				{print(output,"Minimum Amplitude:   %10.4f   Maximum Amplitude:   %10.4f", ampmin,ampmax);PR}
			if (ngsbeams > 0 || verbose >= 1)
				{print(output,"Minimum Sidescan:    %10.4f   Maximum Sidescan:    %10.4f", ssmin,ssmax);PR}
			break;
		case JSON:
			print(output,"\"data_totals\": {");PR
			print(output,"\"number_of_records\":\"%d\"",irec);PR
			isbtmrec = notice_list_tot[MB_DATA_SUBBOTTOM_MCS] + notice_list_tot[MB_DATA_SUBBOTTOM_CNTRBEAM] +
			           notice_list_tot[MB_DATA_SUBBOTTOM_SUBBOTTOM];
			if (isbtmrec > 0)
				{print(output,",\n\"number_of_subbottom_records\":\"%d\"",isbtmrec);PR}
			if (notice_list_tot[MB_DATA_SIDESCAN2] > 0)
				{print(output,",\n\"number_of_secondary_sidescan_records\": \"%d\"",notice_list_tot[MB_DATA_SIDESCAN2]);PR}
			if (notice_list_tot[MB_DATA_SIDESCAN3] > 0)
				{print(output,",\n\"number_of_tertiary_sidescan_records\": \"%d\"",notice_list_tot[MB_DATA_SIDESCAN3]);PR}
			if (notice_list_tot[MB_DATA_WATER_COLUMN] > 0)
				{print(output,",\n\"number_of_water_column_records\": \"%d\"",notice_list_tot[MB_DATA_WATER_COLUMN]);PR}
			print(output,"\n},");PR

			print(output,"\"bathymetry_data\": {\n\"max_beams_per_ping\": \"%d\",",beams_bath_max);PR
			print(output,"\"number_beams\": \"%d\",", ntdbeams);PR
			print(output,"\"number_good_beams\": \"%d\",\n\"percent_good_beams\": \"%5.2f\",", ngdbeams, ngd_percent);PR
			print(output,"\"number_zero_beams\": \"%d\",\n\"percent_zero_beams\": \"%5.2f\",", nzdbeams, nzd_percent);PR
			print(output,"\"number_flagged_beams\": \"%d\",\n\"percent_flagged_beams\": \"%5.2f\"", nfdbeams, nfd_percent);PR
			print(output,"},");PR

			print(output,"\"amplitude_data\": {\n\"max_beams_per_ping\": \"%d\",",beams_amp_max);PR
			print(output,"\"number_beams\": \"%d\",", ntabeams);PR
			print(output,"\"number_good_beams\": \"%d\",\n\"percent_good_beams\": \" %5.2f\",", ngabeams, nga_percent);PR
			print(output,"\"number_zero_beams\": \"%d\",\n\"percent_zero_beams\": \"%5.2f\",", nzabeams, nza_percent);PR
			print(output,"\"number_flagged_beams\": \"%d\",\n\"percent_flagged_beams\": \"%5.2f\"", nfabeams, nfa_percent);PR
			print(output,"},");PR

			print(output,"\"sidescan_data\": {\n\"max_pixels_per_ping\": \"%d\",",pixels_ss_max);PR
			print(output,"\"number_of_pixels\": \"%d\",", ntsbeams);PR
			print(output,"\"number_good_pixels\": \"%d\",\n\"percent_good_pixels\": \"%5.2f\",", ngsbeams, ngs_percent);PR
			print(output,"\"number_zero_pixels\": \"%d\",\n\"percent_zero_pixels\": \"%5.2f\",", nzsbeams, nzs_percent);PR
			print(output,"\"number_flagged_pixels\": \"%d\",\n\"percent_flagged_pixels\": \"%5.2f\"", nfsbeams, nfs_percent);PR
			print(output,"},");PR

			print(output,"\"navigation_totals\": {");PR
			print(output,"\"total_time_hours\": \"%.4f\",",timtot);PR
			print(output,"\"total_track_length_km\": \"%.4f\",",distot);PR
			print(output,"\"average_speed_km_per_hr\": \"%.4f\",\n\"average_speed_knots\": \"%.4f\"", spdavg,spdavg/1.85);PR
			print(output,"},");PR

			print(output,"\"start_of_data\": {");PR
			print(output,"\"time\": \"%2.2d %2.2d %4.4d %2.2d:%2.2d:%2.2d.%6.6d  JD%d\",",
			             timbeg_i[1],timbeg_i[2],timbeg_i[0],timbeg_i[3],
			             timbeg_i[4],timbeg_i[5],timbeg_i[6],timbeg_j[1]);PR
			print(output,"\"time_iso\": \"%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%6.6d\",",
			             timbeg_i[0],timbeg_i[1],timbeg_i[2],timbeg_i[3],timbeg_i[4],timbeg_i[5],timbeg_i[6]);PR
			if (bathy_in_feet == MB_NO) {
				print(output,"\"longitude\": \"%.9f\",\n\"latitude\": \"%.9f\",\n\"depth_meters\": \"%.4f\",",
				             lonbeg,latbeg,bathbeg);PR
			}
			else
				{print(output,"\"longitude\": \"%.9f\",\n\"latitude\": \"%.9f\",\n\"depth_feet\": \"%.4f\",",
					lonbeg,latbeg,bathy_scale*bathbeg);PR}
			print(output,"\"speed_km_per_hour\": \"%.4f\",\n\"speed_knots\": \"%.4f\",\n\"heading_degrees\": \"%.4f\",",
				spdbeg,spdbeg/1.85,hdgbeg);PR
			print(output,"\"sonar_depth_meters\": \"%.4f\",\n\"sonar_altitude_meters\": \"%.4f\"", sdpbeg,altbeg);PR
			print(output,"},");PR

			print(output,"\"end_of_data\": {");PR
			print(output,"\"time\": \"%2.2d %2.2d %4.4d %2.2d:%2.2d:%2.2d.%6.6d  JD%d\",",
			             timend_i[1],timend_i[2],timend_i[0],timend_i[3],
			             timend_i[4],timend_i[5],timend_i[6],timend_j[1]);PR
			print(output,"\"time_iso\": \"%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%6.6d\",",
			             timend_i[0],timend_i[1],timend_i[2],timend_i[3],timend_i[4],timend_i[5],timend_i[6]);PR
			if (bathy_in_feet == MB_NO) {
				print(output,"\"longitude\": \"%.9f\",\n\"latitude\": \"%.9f\",\n\"depth_meters\": \"%.4f\",",
				             lonend,latend,bathend);PR
			}
			else
				{print(output,"\"longitude\": \"%.9f\",\n\"latitude\": \"%.9f\",\n\"depth_feet\": \"%.4f\",",
				              lonend,latend,bathy_scale*bathend);PR}
			print(output,"\"speed_km_per_hour\": \"%.4f\",\n\"speed_knots\": \"%.4f\",\n\"heading_degrees\": \"%.4f\",",
			             spdend,spdend/1.85,hdgend);PR
			print(output,"\"sonar_depth_meters\": \"%.4f\",\n\"sonar_altitude_meters\": \"%.4f\"", sdpend,altend);PR
			print(output,"},");PR

			print(output,"\"limits\": {");PR
			print(output,"\"minimum_longitude\": \"%.9f\",\n\"maximum_longitude\": \"%.9f\",",lonmin,lonmax);PR
			print(output,"\"minimum_latitude\": \"%.9f\",\n\"maximum_latitude\": \"%.9f\",",latmin,latmax);PR
			print(output,"\"minimum_sonar_depth\": \"%.4f\",\n\"maximum_sonar_depth\": \"%.4f\",",sdpmin,sdpmax);PR
			print(output,"\"minimum_altitude\": \"%.4f\",\n\"maximum_altitude\": \"%.4f\"",altmin,altmax);PR
			if (ngdbeams > 0 || verbose >= 1)
				{print(output,",\n\"minimum_depth\": \"%.4f\",\n\"maximum_depth\": \"%.4f\"",
				              bathy_scale*bathmin,bathy_scale*bathmax);PR}
			if (ngabeams > 0 || verbose >= 1)
				{print(output,",\n\"minimum_amplitude\": \"%.4f\",\n\"maximum_amplitude\": \"%.4f\"", ampmin,ampmax);PR}
			if (ngsbeams > 0 || verbose >= 1)
				{print(output,",\n\"minimum_sidescan\": \"%.4f\",\n\"maximum_sidescan\": \"%.4f\"", ssmin,ssmax);PR}
			print(output,"\n}");PR
			break;
		case XML:
			print(output,"\t<data_totals>");PR
			print(output,"\t\t<number_of_records>%d</number_of_records>",irec);PR
			isbtmrec = notice_list_tot[MB_DATA_SUBBOTTOM_MCS] + notice_list_tot[MB_DATA_SUBBOTTOM_CNTRBEAM] +
					   notice_list_tot[MB_DATA_SUBBOTTOM_SUBBOTTOM];
			if (isbtmrec > 0)
				{print(output,"\t\t<number_of_subbottom_records>%d</number_of_subbottom_records>",isbtmrec);PR}
			if (notice_list_tot[MB_DATA_SIDESCAN2] > 0)
				{print(output,"\t\t<number_of_secondary_sidescan_records>%d</number_of_secondary_sidescan_records>",
				                notice_list_tot[MB_DATA_SIDESCAN2]);PR}
			if (notice_list_tot[MB_DATA_SIDESCAN3] > 0)
				{print(output,"\t\t<number_of_tertiary_sidescan_records>%d</number_of_tertiary_sidescan_records>",
				                notice_list_tot[MB_DATA_SIDESCAN3]);PR}
			if (notice_list_tot[MB_DATA_WATER_COLUMN] > 0)
				{print(output,"\t\t<number_of_water_column_records>%d</number_of_water_column_records>",
				                notice_list_tot[MB_DATA_WATER_COLUMN]);PR}
			print(output,"\t</data_totals>");PR

			print(output,"\t<bathymetry_data>");PR
			print(output,"\t\t<max_beams_per_ping>%d</max_beams_per_ping>",beams_bath_max);PR
			print(output,"\t\t<number_beams>%d</number_beams>",ntdbeams);PR
			print(output,"\t\t<number_good_beams>%d</number_good_beams>",ngdbeams);PR
			print(output,"\t\t<percent_good_beams>%.2f</percent_good_beams>",ngd_percent);PR
			print(output,"\t\t<number_zero_beams>%d</number_zero_beams>",nzdbeams);PR
			print(output,"\t\t<percent_zero_beams>%.2f</percent_zero_beams>",nzd_percent);PR
			print(output,"\t\t<number_flagged_beams>%d</number_flagged_beams>",nfdbeams);PR
			print(output,"\t\t<percent_flagged_beams>%.2f</percent_flagged_beams>",nfd_percent);PR
			print(output,"\t</bathymetry_data>");PR
			print(output,"\t<amplitude_data>");PR
			print(output,"\t\t<max_beams_per_ping>%d</max_beams_per_ping>",beams_bath_max);PR
			print(output,"\t\t<number_beams>%d</number_beams>",ntabeams);PR
			print(output,"\t\t<number_good_beams>%d</number_good_beams>",ngabeams);PR
			print(output,"\t\t<percent_good_beams>%.2f</percent_good_beams>",nga_percent);PR
			print(output,"\t\t<number_zero_beams>%d</number_zero_beams>",nzabeams);PR
			print(output,"\t\t<percent_zero_beams>%.2f</percent_zero_beams>",nza_percent);PR
			print(output,"\t\t<number_flagged_beams>%d</number_flagged_beams>",nfabeams);PR
			print(output,"\t\t<percent_flagged_beams>%.2f</percent_flagged_beams>",nfa_percent);PR
			print(output,"\t</amplitude_data>");PR
			print(output,"\t<sidescan_data>");PR
			print(output,"\t\t<max_pixels_per_ping>%d</max_pixels_per_ping>",pixels_ss_max);PR
			print(output,"\t\t<number_pixels>%d</number_pixels>",ntsbeams);PR
			print(output,"\t\t<number_good_pixels>%d</number_good_pixels>",ngsbeams);PR
			print(output,"\t\t<percent_good_pixels>%.2f</percent_good_pixels>",ngs_percent);PR
			print(output,"\t\t<number_zero_pixels>%d</number_zero_pixels>",nzsbeams);PR
			print(output,"\t\t<percent_zero_pixels>%.2f</percent_zero_pixels>",nzs_percent);PR
			print(output,"\t\t<number_flagged_pixels>%d</number_flagged_pixels>",nfsbeams);PR
			print(output,"\t\t<percent_flagged_pixels>%.2f</percent_flagged_pixels>",nfs_percent);PR
			print(output,"\t</sidescan_data>");PR

			print(output,"\t<tnavigation_totals>");PR
			print(output,"\t\t<total_time_hours>%.4f</total_time_hours>",timtot);PR
			print(output,"\t\t<total_track_length_km>%.4f</total_track_length_km>",distot);PR
			print(output,"\t\t<average_speed_km_per_hr>%.4f</average_speed_km_per_hr>",spdavg);PR
			print(output,"\t\t<average_speed_knots>%.4f</average_speed_knots>",spdavg/1.85);PR
			print(output,"\t</tnavigation_totals>");PR

			print(output,"\t<start_of_data>");PR
			print(output,"\t\t<time>%2.2d %2.2d %4.4d %2.2d:%2.2d:%2.2d.%6.6d  JD%d</time>",
			             timbeg_i[1],timbeg_i[2],timbeg_i[0],timbeg_i[3],
			             timbeg_i[4],timbeg_i[5],timbeg_i[6],timbeg_j[1]);PR
			print(output,"\t\t<time_iso>%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%6.6d</time_iso>",
			             timbeg_i[0],timbeg_i[1],timbeg_i[2],timbeg_i[3],timbeg_i[4],timbeg_i[5],timbeg_i[6]);PR
			if (bathy_in_feet == MB_NO) {
				print(output,"\t\t<longitude>%.9f</longitude>",lonbeg);PR
				print(output,"\t\t<latitude>%.9f</latitude>",latbeg);PR
				print(output,"\t\t<depth_meters>%.4f</depth_meters>",bathbeg);PR
			}
			else {
				print(output,"\t\t<longitude>%.9f</longitude>",lonbeg);PR
				print(output,"\t\t<latitude>%.9f</latitude>",latbeg);PR
				print(output,"\t\t<depth_meters>%.4f</depth_meters>",bathy_scale*bathbeg);PR
			}
			print(output,"\t\t<speed_km_per_hour>%.4f</speed_km_per_hour>",spdbeg);PR
			print(output,"\t\t<speed_knots>%.4f</speed_knots>",spdbeg/1.85);PR
			print(output,"\t\t<heading_degrees>%.4f</heading_degrees>",hdgbeg);PR
			print(output,"\t\t<sonar_depth_meters>%.4f</sonar_depth_meters>",sdpbeg);PR
			print(output,"\t\t<sonar_altitude_meters>%.4f</sonar_altitude_meters>",altbeg);PR
			print(output,"\t</start_of_data>");PR

			print(output,"\t<end_of_data>");PR
			print(output,"\t\t<time>%2.2d %2.2d %4.4d %2.2d:%2.2d:%2.2d.%6.6d  JD%d</time>",
			               timend_i[1],timend_i[2],timend_i[0],timend_i[3],
			               timend_i[4],timend_i[5],timend_i[6],timend_j[1]);PR
			print(output,"\t\t<time_iso>%4.4d-%2.2d-%2.2dT%2.2d:%2.2d:%2.2d.%6.6d</time_iso>",
			               timend_i[0],timend_i[1],timend_i[2],timend_i[3],timend_i[4],timend_i[5],timend_i[6]);PR
			if (bathy_in_feet == MB_NO) {
				print(output,"\t\t<longitude>%.9f</longitude>",lonend);PR
				print(output,"\t\t<latitude>%.9f</latitude>",latend);PR
				print(output,"\t\t<depth_meters>%.4f</depth_meters>",bathend);PR
			}
			else {
				print(output,"\t\t<longitude>%.9f</longitude>",lonend);PR
				print(output,"\t\t<latitude>%.9f</latitude>",latend);PR
				print(output,"\t\t<depth_meters>%.4f</depth_meters>",bathy_scale*bathend);PR
			}
			print(output,"\t\t<speed_km_per_hour>%.4f</speed_km_per_hour>",spdend);PR
			print(output,"\t\t<speed_knots>%.4f</speed_knots>",spdend/1.85);PR
			print(output,"\t\t<heading_degrees>%.4f</heading_degrees>",hdgend);PR
			print(output,"\t\t<sonar_depth_meters>%.4f</sonar_depth_meters>",sdpend);PR
			print(output,"\t\t<sonar_altitude_meters>%.4f</sonar_altitude_meters>",altend);PR
			print(output,"\t</end_of_data>");PR

			print(output,"\t<limits>");PR
			print(output,"\t\t<minimum_longitude>%.9f</minimum_longitude>",lonmin);PR
			print(output,"\t\t<maximum_longitude>%.9f</maximum_longitude>",lonmax);PR
			print(output,"\t\t<minimum_latitude>%.9f</minimum_latitude>",latmin);PR
			print(output,"\t\t<maximum_latitude>%.9f</maximum_latitude>",latmax);PR
			print(output,"\t\t<minimum_sonar_depth>%.4f</minimum_sonar_depth>",sdpmin);PR
			print(output,"\t\t<maximum_sonar_depth>%.4f</maximum_sonar_depth>",sdpmax);PR
			print(output,"\t\t<minimum_altitude>%.4f</minimum_altitude>",altmin);PR
			print(output,"\t\t<maximum_altitude>%.4f</maximum_altitude>",altmax);PR
			if (ngdbeams > 0 || verbose >= 1) {
				print(output,"\t\t<minimum_depth>%.4f</minimum_depth>",bathy_scale*bathmin);PR
				print(output,"\t\t<maximum_depth>%.4f</maximum_depth>",bathy_scale*bathmax);PR
			}
			if (ngabeams > 0 || verbose >= 1) {
				print(output,"\t\t<minimum_amplitude>%.4f</minimum_amplitude>",ampmin);PR
				print(output,"\t\t<maximum_amplitude>%.4f</maximum_amplitude>",ampmax);PR
			}
			if (ngsbeams > 0 || verbose >= 1) {
				print(output,"\t\t<minimum_sidescan>%.4f</minimum_sidescan>",ssmin);PR
				print(output,"\t\t<maximum_sidescan>%.4f</maximum_sidescan>",ssmax);PR
			}
			print(output,"\t</limits>");PR
			break;
		case '?':
			break;
	}
	if (pings_read > 2 && beams_bath_max > 0 && (ngdbeams > 0 || verbose >= 1)) {
		switch (output_format) {
			case FREE_TEXT:
				print(output,"\nBeam Bathymetry Variances:");PR
				print(output,"Pings Averaged: %d",pings_read);PR
				print(output," Beam     N      Mean     Variance    Sigma");PR
				print(output," ----     -      ----     --------    -----");PR
				for (i=0;i<beams_bath_max;i++) {
					print(output,"%4d  %5d   %8.2f   %8.2f  %8.2f", i,nbathvartot[i],bathy_scale*bathmeantot[i],
					             bathy_scale*bathy_scale*bathvartot[i], bathy_scale*sqrt(bathvartot[i]));PR
				}
				print(output,"");PR
				break;
			case JSON:
				print(output,",\n\"beam_bathymetry_variances\":{");PR
				print(output,"\"pings_averaged\": \"%d\",",pings_read);PR
				print(output,"\"columns\" : \"#beam,N,mean,variance,sigma\",");PR
				print(output,"\"values\": [");PR
				for (i=0;i<beams_bath_max;i++) {
					if(i>0) {print(output,",");PR}
					sigma=bathy_scale*sqrt(bathvartot[i]);
					if (isnan(sigma)) sigma=0;
					print(output,"{\"row\":\"%d,%d,%.2f,%.2f,%.2f\"}", i,nbathvartot[i],bathy_scale*bathmeantot[i],
					               bathy_scale*bathy_scale*bathvartot[i], sigma);PR
				}
				print(output,"]}");PR
				break;
			case XML:
				print(output,"\t<beam_bathymetry_variances>");PR
				print(output,"\t\t<pings_averaged>%d</pings_averaged>",pings_read);PR
				print(output,"\t\t<columns>pixel,N,mean,variance,sigma</columns>");PR
				print(output,"\t\t<values>");PR
				for (i=0;i<beams_bath_max;i++) {
					if(i>0)
						sigma=bathy_scale*sqrt(bathvartot[i]);
					if(isnan(sigma)) sigma=0;
					print(output,"\t\t\t<row>%d,%d,%.2f,%.2f,%.2f</row>", i,nbathvartot[i],bathy_scale*bathmeantot[i],
					             bathy_scale*bathy_scale*bathvartot[i], sigma);PR
				}
				print(output,"\t\t</values>");PR
				print(output,"\t</beam_bathymetry_variances>");PR
				break;
			case '?':
				break;
		}
	}
	if (pings_read > 2 && beams_amp_max > 0 && (ngabeams > 0 || verbose >= 1)) {
		switch (output_format) {
			case FREE_TEXT:
				print(output,"\nBeam Amplitude Variances:");PR
				print(output,"Pings Averaged: %d",pings_read);PR
				print(output," Beam     N      Mean     Variance    Sigma");PR
				print(output," ----     -      ----     --------    -----");PR
				for (i=0;i<beams_amp_max;i++) {
					print(output,"%4d  %5d   %8.2f   %8.2f  %8.2f", i,nampvartot[i],ampmeantot[i],
					             ampvartot[i],sqrt(ampvartot[i]));PR
				}
				print(output,"");PR
				break;
			case JSON:
				print(output,",\n\"beam_amplitude_variances\":{");PR
				print(output,"\"pings_averaged\": \"%d\",",pings_read);PR
				print(output,"\"columns\":\"beam,N,mean,variance,sigma\",");PR
				print(output,"\"values\": [");PR
				for (i=0;i<beams_amp_max;i++) {
					if (i > 0) {print(output,",");PR}
					sigma=sqrt(ampvartot[i]);
					if (isnan(sigma)) sigma=0;
					print(output,"{\"row\" : \"%d,%d,%.2f,%.2f,%.2f\"}", i,nampvartot[i],ampmeantot[i], ampvartot[i],sigma);PR
				}
				print(output,"\n]}");PR
				break;
			case XML:
				print(output,"\t<beam_amplitude_variances>");PR
				print(output,"\t\t<pings_averaged>%d</pings_averaged>",pings_read);PR
				print(output,"\t\t<columns>pixel,N,mean,variance,sigma</columns>");PR
				print(output,"\t\t<values>");PR
				for (i=0;i<beams_amp_max;i++) {
					if (i > 0)
						sigma=sqrt(ampvartot[i]);
					if(isnan(sigma)) sigma=0;
					print(output,"\t\t\t<row>%d,%d,%.2f,%.2f,%.2f</row>", i,nampvartot[i],ampmeantot[i], ampvartot[i],sigma);PR
				}
				print(output,"\t\t</values>");PR
				print(output,"\t</beam_amplitude_variances>");PR
				break;
			case '?':
				break;
		}
	}
	if (pings_read > 2 && pixels_ss_max > 0 && (ngsbeams > 0 || verbose >= 1)) {
		switch (output_format) {
			case FREE_TEXT:
				print(output,"\nPixel Sidescan Variances:");PR
				print(output,"Pings Averaged: %d",pings_read);PR
				print(output," Beam     N      Mean     Variance    Sigma");PR
				print(output," ----     -      ----     --------    -----");PR
				for (i=0;i<pixels_ss_max;i++) {
					print(output,"%4d  %5d   %8.2f   %8.2f  %8.2f", i,nssvartot[i],ssmeantot[i], ssvartot[i],sqrt(ssvartot[i]));PR
				}
				print(output,"");PR
				break;
			case JSON:
				print(output,",\n\"pixel_sidescan_variances\":{");PR
				print(output,"\"pings_averaged\": \"%d\",",pings_read);PR
				print(output,"\"columns\":\"pixel,N,mean,variance,sigma\",");PR
				print(output,"\"values\": [");PR
				for (i=0;i<pixels_ss_max;i++) {
					if (i>0) {print(output,",");PR}
					sigma=sqrt(ssvartot[i]);
					if (isnan(sigma)) sigma=0;
					print(output,"{\"row\":\"%d,%d,%.2f,%.2f,%.2f\"}", i,nssvartot[i],ssmeantot[i], ssvartot[i],sigma);PR
				}
				print(output,"\n]\n}");PR
				break;
			case XML:
				print(output,"\t<pixel_sidescan_variances>");PR
				print(output,"\t\t<pings_averaged>%d</pings_averaged>",pings_read);PR
				print(output,"\t\t<columns>pixel,N,mean,variance,sigma</columns>");PR
				print(output,"\t\t<values>");PR
				for (i=0;i<pixels_ss_max;i++) {
					if (i>0)
						sigma=sqrt(ssvartot[i]);
					if (isnan(sigma)) sigma=0;
					print(output,"\t\t\t<row>%d,%d,%.2f,%.2f,%.2f</row>", i,nssvartot[i],ssmeantot[i], ssvartot[i],sigma);PR
				}
				print(output,"\t\t</values>");PR
				print(output,"\t</pixel_sidescan_variances>");PR
				break;
			case '?':
				break;
		}
	}
	if (print_notices == MB_YES) {
		switch (output_format) {
			case FREE_TEXT:
				print(output,"\nData Record Type Notices:");PR
				for (i=0;i<=MB_DATA_KINDS;i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						print(output, "DN: %d %s", notice_list_tot[i], notice_msg);PR
					}
				}
				print(output,"\nNonfatal Error Notices:");PR
				for (i=MB_DATA_KINDS+1;i<=MB_DATA_KINDS-(MB_ERROR_MIN);i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						print(output, "EN: %d %s", notice_list_tot[i], notice_msg);PR
					}
				}
				print(output,"\nProblem Notices:");PR
				for (i=MB_DATA_KINDS-(MB_ERROR_MIN)+1;i<MB_NOTICE_MAX;i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						print(output, "PN: %d %s", notice_list_tot[i], notice_msg);PR
					}
				}
				break;
			case JSON:
				print(output,",\n\"notices\": {");PR
				notice_total=0;
				print(output,"\"data_record_type_notices\": [");PR
				for (i=0;i<=MB_DATA_KINDS;i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						if (notice_total>0) {print(output,",");PR}
						print(output, "{\"notice\": {\n\"notice_number\": \"%d\",\n\"notice_message\":\"%s\"\n}}",
						              notice_list_tot[i], notice_msg);PR
						notice_total++;
					}
				}
				if (notice_total>0) {print(output,"");PR}
				print(output,"]");PR
				notice_total=0;
				print(output,",\n\"nonfatal_error_notices\": [");PR
				for (i=MB_DATA_KINDS+1;i<=MB_DATA_KINDS-(MB_ERROR_MIN);i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						if (notice_total>0) {print(output,",");PR}
						print(output, "{\"notice\": {\n\"notice_number\": \"%d\",\n\"notice_message\":\"%s\"\n}}",
						              notice_list_tot[i], notice_msg);PR
						notice_total++;
					}
				}
				if (notice_total>0) {print(output,"");PR}
				print(output,"]");PR
				notice_total=0;
				print(output,",\n\"problem_notices\": [");PR
				for (i=MB_DATA_KINDS-(MB_ERROR_MIN)+1;i<MB_NOTICE_MAX;i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						if (notice_total>0) {print(output,",");PR}
						print(output, "{\"notice\": {\n\"notice_number\": \"%d\",\n\"notice_message\":\"%s\"\n}}",
						              notice_list_tot[i], notice_msg);PR
						notice_total++;
					}
				}
				if (notice_total>0) {print(output,"");PR}
				print(output,"]");PR
				print(output,"}");PR
				break;
			case XML:
				print(output,"\t<data_record_type_notices>");PR
				for (i=0;i<=MB_DATA_KINDS;i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						print(output, "\t\t<notice_number>%d</notice_number>",notice_list_tot[i]);PR
						print(output, "\t\t<notice_messsage>%s</notice_messsage>", notice_msg);PR
					}
				}
				print(output,"\t</data_record_type_notices>");PR
				print(output,"\t<nonfatal_error_notices>");PR
				for (i=MB_DATA_KINDS+1;i<=MB_DATA_KINDS-(MB_ERROR_MIN);i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						print(output, "\t\t<notice_number>%d</notice_number>",notice_list_tot[i]);PR
						print(output, "\t\t<notice_messsage>%s</notice_messsage>", notice_msg);PR
					}
				}
				print(output,"\t</nonfatal_error_notices>");PR
				print(output,"\t<problem_notices>");PR
				for (i=MB_DATA_KINDS-(MB_ERROR_MIN)+1;i<MB_NOTICE_MAX;i++) {
					if (notice_list_tot[i] > 0) {
						mb_notice_message(verbose, i, &notice_msg);
						print(output, "\t\t<notice_number>%d</notice_number>",notice_list_tot[i]);PR
						print(output, "\t\t<notice_messsage>%s</notice_messsage>", notice_msg);PR
					}
				}
				print(output,"\t</problem_notices>");PR
				break;
			case '?':
				break;
		}
	}
	if (coverage_mask == MB_YES) {
		switch (output_format) {
			case FREE_TEXT:
				print(output,"\nCoverage Mask:\nCM dimensions: %d %d", mask_nx, mask_ny);PR
				for (j=mask_ny-1;j>=0;j--) {
					print(output, "CM:  ");PR
					for (i=0;i<mask_nx;i++) {
						k = i + j * mask_nx;
						print(output, " %1d", mask[k]);PR
					}
					print(output, "");PR
				}
				break;
			case JSON:
				print(output,",\n\"coverage_mask\": {");PR
				print(output,"\"dimensions_nx\": \"%d\",\n\"dimensions_ny\": \"%d\",", mask_nx, mask_ny);PR
				print(output, "\"mask\": \" ");PR
				for (j=mask_ny-1;j>=0;j--) {
					for (i=0;i<mask_nx;i++) {
						k = i + j * mask_nx;
						if (i>0) {print(output,",");PR}
						print(output, "%1d", mask[k]);PR
					}
					print(output,"");PR
				}
				print(output, "\"}");PR
				break;
			case '?':
				break;
		}
	}

	/* close root element for XML export */
	switch (output_format) {
		case FREE_TEXT:
			break;
		case JSON:
			print(output,"}");PR
			break;
		case XML:
			print(output,"</mbinfo>");PR
			break;
		case '?':
			break;
	}

	/* close output file */
	//if (output_usefile == MB_YES && output != NULL)
		//fclose(output);

	/* deallocate memory used for data arrays */
	mb_freed(verbose,__FILE__,__LINE__,(void **)&bathmeantot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&bathvartot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&nbathvartot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&ampmeantot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&ampvartot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&nampvartot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&ssmeantot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&ssvartot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&nssvartot,&error);
	mb_freed(verbose,__FILE__,__LINE__,(void **)&mask,&error);

	/* set program status */
	status = MB_SUCCESS;

	/* check memory */
	if (verbose >= 4)
		status = mb_memory_list(verbose,&error);

	if (GMT_End_IO (API, GMT_IN,  0) != GMT_NOERROR) {	/* Disables further data input */
		Return (API->error);
	}
	if (GMT_End_IO (API, GMT_OUT, 0) != GMT_NOERROR) {	/* Disables further data output */
		Return (API->error);
	}
	
	return error;
}
