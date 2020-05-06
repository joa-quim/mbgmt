/*--------------------------------------------------------------------
 *    The MB-system:	mbflags.c	2/26/93
 *    $Id:  $
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * Author:	J. Luis
 * Date:	May, 2017
 */

#define THIS_MODULE_MODERN_NAME	"mbflags"
#define THIS_MODULE_CLASSIC_NAME	"mbflags"
#define THIS_MODULE_LIB		"mbgmt"
#define THIS_MODULE_PURPOSE	"Update the flags .esf file from external file"
#define THIS_MODULE_KEYS	""
#define THIS_MODULE_NEEDS	""
#define THIS_MODULE_OPTIONS "->V"

/* GMT5 header file */
#include "gmt_dev.h"

#define GMT_PROG_OPTIONS "->V"

/* mbio include files */
#include "mb_status.h"
#include "mb_format.h"
#include "mb_define.h"
#include "mb_io.h"
#include "mb_swap.h"
#include "mb_process.h"

/* ping structure definition */
struct mbclean_ping_struct {
	int	time_i[7];
	double	time_d;
	double	navlon;
	double	navlat;
	double	speed;
	double	heading;
	int		beams_bath;
	char	*beamflag;
	char	*beamflagorg;
	double	*bath;
	double	*bathacrosstrack;
	double	*bathalongtrack;
};

EXTERN_MSC int GMT_mbflags(void *API, int mode, void *args);

/* Control structure for mbcontour */
struct MBFLAGS_CTRL {

	struct mbflg_E {	/* -E */
		bool   active;
		char   *file;
	} E;
	struct mbflg_F {	/* -F */
		bool   active;
		int    format;
	} F;
	struct mbflg_I {	/* -I */
		bool   active;
		char  *file;
	} I;
	struct mbflg_L {	/* -L*/
		bool   active;
		int    lonflip;
	} L;

	struct mbflg_tranf {
		int    fix_edit_timestamps;
		int    zap_beams;
		int    check_num_good_min;
	} transf;
};

static char rcs_id[] = "$Id: mbflags.c  $";

/*-----------------------------------------------------------------------------------------*/
void *New_MBFLAGS_CTRL (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct MBFLAGS_CTRL *Ctrl;

	Ctrl = gmt_M_memory (GMT, NULL, 1, struct MBFLAGS_CTRL);

	Ctrl->F.format = 0;

	return (Ctrl);
}

void Free_MBFLAGS_CTRL (struct GMT_CTRL *GMT, struct MBFLAGS_CTRL *Ctrl) {	/* Deallocate control structure */
	if (!Ctrl) return;

	if (Ctrl->I.file) free (Ctrl->I.file);
	if (Ctrl->E.file) free (Ctrl->E.file);
	gmt_M_free (GMT, Ctrl);
}

GMT_LOCAL int usage (struct GMTAPI_CTRL *API, int level) {

	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);
	GMT_Message (API, GMT_TIME_NONE, "usage: mbflags -I<inputfile> -Fflags_file -V\\n");

	if (level == GMT_SYNOPSIS) return (EXIT_FAILURE);

	GMT_Message (API, GMT_TIME_NONE, "\t<inputfile> is an MB-System datalist referencing the swath data to be cleaned.\n");
	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	return (EXIT_FAILURE);
}

GMT_LOCAL int parse (struct GMT_CTRL *GMT, struct MBFLAGS_CTRL *Ctrl, struct GMT_OPTION *options) {

	unsigned int n_errors = 0, n_files = 0;
	int    n;
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

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

			case 'I':	/* */
				if (!gmt_access (GMT, opt->arg, R_OK)) {	/* Got a file */
					Ctrl->I.file = strdup (opt->arg);
					Ctrl->I.active = true;
					n_files = 1;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -I option: \n");
					n_errors++;
				}
 				break;


			case 'L':	/* */
				n = sscanf(opt->arg, "%d", &(Ctrl->L.lonflip));
				if (n > 0) {
					Ctrl->L.active = true;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -L option: \n");
					n_errors++;
				}
 				break;

			case 'E':	/* Read external file with info to be sent to the .esf file */
				Ctrl->E.active = true;
				Ctrl->E.file = strdup(opt->arg);
 				break;

			case 'F':	/* */
				n = sscanf(opt->arg, "%d", &(Ctrl->F.format));
				if (n > 0) {
					Ctrl->F.active = true;
				}
				else {
					GMT_Report (API, GMT_MSG_NORMAL, "Syntax error -F option: \n");
					n_errors++;
				}
 				break;

			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, n_files != 1, "Syntax error: Must specify one input file(s)\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_MBFLAGS_CTRL (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

/*-----------------------------------------------------------------------------------------*/
int GMT_mbflags (void *V_API, int mode, void *args) {

	char program_name[] = "MBCLEAN";

	struct MBFLAGS_CTRL *Ctrl = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_GRID *G = NULL;
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */

	/* MBIO status variables */
	int	status;
	int	verbose = 0;
	int	error = MB_ERROR_NO_ERROR;
	char	*message = NULL;

	/* swath file locking variables */
	int	uselockfiles;
	int	lock_status;
	int	locked;
	int	lock_purpose;
	mb_path	lock_program;
	mb_path lock_cpu;
	mb_path lock_user;
	char	lock_date[25];

	/* MBIO read control parameters */
	int	look_processed = MB_DATALIST_LOOK_UNSET;
	int	format;
	int	formatread;
	int	variable_beams;
	int	traveltime;
	int	beam_flagging;
	int	pings;
	int	ndata = 0;
	int	nzeropostot = 0;
	int nlong_acrosstot=0;
	int	nmin = 0;
	int	nbad = 0;
	int	nflag = 0;
	int	nunflag = 0;
	int	nflagesf = 0;
	int	nunflagesf = 0;
	int	nzeroesf = 0;
	int	nzeropos = 0;
	int	ndatatot = 0;
	int	nflagesftot = 0;
	int	nunflagesftot = 0;
	int	nzeroesftot = 0;
	int	nflagtot = 0;
	int	nunflagtot = 0;
	double	*list = NULL;
	double	bounds[4];
	int	btime_i[7];
	int	etime_i[7];
	double	btime_d;
	double	etime_d;
	double	speedmin;
	double	timegap;
	double	distance;
	double	altitude;
	double	sonardepth;
	int     beams_bath, beams_amp, pixels_ss;
	double *amp = NULL, *ss = NULL, *ssacrosstrack = NULL, *ssalongtrack = NULL;

	/* mbio read and write values */
	void	*mbio_ptr = NULL;
	int	kind;
	struct mbclean_ping_struct ping[1];
	int	nrec;
	int	pingsread;
	char	comment[MB_COMMENT_MAXLINE];
	int	action;

	/* fix_edit_timestamps variables */
	int	fix_edit_timestamps = MB_NO;
	double	tolerance = 0.0;

	/* save file control variables */
	int	esffile_open = MB_NO;
	char	esffile[MB_PATH_MAXLINE];
	struct mb_esf_struct esf;

	/* processing variables */
	int	i, k, start, done;

	char *row_flags;
	int dims[2] = {0};
	FILE *fid;

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) 		/* Return the purpose of program */
		return (usage (API, GMT_MODULE_PURPOSE));
	options = GMT_Create_Options (API, mode, args);
	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE)	/* Return the usage message */
		bailout (usage (API, GMT_USAGE));
	if (options->option == GMT_OPT_SYNOPSIS) 			/* Return the synopsis */
		bailout (usage (API, GMT_SYNOPSIS));

	/* Parse the command-line arguments */

#if GMT_MAJOR_VERSION >= 6
	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, NULL, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
#else
	GMT = gmt_begin_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, &GMT_cpy); /* Save current state */
#endif
	if (GMT_Parse_Common (API, GMT_PROG_OPTIONS, options)) Return(API->error);
               
	Ctrl = New_MBFLAGS_CTRL (GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse (GMT, Ctrl, options))) Return(error);

	/* Copy the parsing variables stored in the Ctrl struct to original names */
	verbose = GMT->common.V.active;

	/* get current default values */
	status = mb_defaults(verbose,&format,&pings,&Ctrl->L.lonflip,bounds, btime_i,etime_i,&speedmin,&timegap);
	status = mb_uselockfiles(verbose,&uselockfiles);

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
	timegap = 1000000000.0;

	/*---------------------------- This is the mbclean main code ----------------------------*/

	if ((fid = fopen(Ctrl->E.file, "rb")) == NULL) {
		GMT_Report (API, GMT_MSG_NORMAL, "ERROR: failed to open file %s\n", Ctrl->E.file);
		Return(-1);
	}
	fread(dims, sizeof(int), 2, fid);
	row_flags = gmt_M_memory(GMT, NULL, dims[1], char);

	GMT_Report (API, GMT_MSG_VERBOSE, "Processing input files\nVersion %s\nMB-system Version %s\n",rcs_id,MB_VERSION);

	/* get format if required */
	if (format == 0)
		mb_get_format(verbose,Ctrl->I.file,NULL,&format,&error);

	/* determine whether to read one file or a list of files */
	if (format < 0) {
		GMT_Report(API, GMT_MSG_NORMAL, "This program does not support reading feom datalist file\n");
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(-1);
	}

	/* check format and get format flags */
	if ((status = mb_format_flags(verbose,&format, &variable_beams, &traveltime, &beam_flagging, &error)) != MB_SUCCESS) {
		mb_error(verbose,error,&message);
		GMT_Report(API, GMT_MSG_NORMAL, "MBIO Error returned from function <mb_format_flags> regarding input format %d:\n%s\n",
					format,message);
		GMT_Report(API, GMT_MSG_NORMAL, "File <%s> skipped by program <%s>\n", Ctrl->I.file,program_name);
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(error);
	}

	/* warn if beam flagging not supported for the current data format */
	if (beam_flagging == MB_NO) {
		GMT_Report(API, GMT_MSG_NORMAL, "Warning:MBIO format %d does not allow flagging of bad bathymetry data.\n",format);
		GMT_Report(API, GMT_MSG_NORMAL, "When mbprocess applies edits to file:\n\t%s\nthe soundings will be nulled (zeroed) rather than flagged.\n", Ctrl->I.file);
	}

	/* try to lock file */
	if (uselockfiles == MB_YES)
		status = mb_pr_lockswathfile(verbose, Ctrl->I.file, MBP_LOCK_EDITBATHY, program_name, &error);
	else {
		lock_status = mb_pr_lockinfo(verbose, Ctrl->I.file, &locked, &lock_purpose,
										lock_program, lock_user, lock_cpu, lock_date, &error);

		/* if locked get lock info */
		if (error == MB_ERROR_FILE_LOCKED) {
			GMT_Report(API, GMT_MSG_NORMAL, "File %s locked but lock ignored\n", Ctrl->I.file);
			GMT_Report(API, GMT_MSG_NORMAL, "File locked by <%s> running <%s>\n", lock_user, lock_program);
			GMT_Report(API, GMT_MSG_NORMAL, "on cpu <%s> at <%s>\n", lock_cpu, lock_date);
			error = MB_ERROR_NO_ERROR;
		}
	}

	/* if locked let the user know file can't be opened */
	if (status == MB_FAILURE) {
		/* if locked get lock info */
		if (error == MB_ERROR_FILE_LOCKED) {
			lock_status = mb_pr_lockinfo(verbose, Ctrl->I.file, &locked, &lock_purpose, lock_program,
			                             lock_user, lock_cpu, lock_date, &error);
			GMT_Report(API, GMT_MSG_NORMAL, "Unable to open input file:\n");
			GMT_Report(API, GMT_MSG_NORMAL, "  %s\n", Ctrl->I.file);
			GMT_Report(API, GMT_MSG_NORMAL, "File locked by <%s> running <%s>\n", lock_user, lock_program);
			GMT_Report(API, GMT_MSG_NORMAL, "on cpu <%s> at <%s>\n", lock_cpu, lock_date);
		}
		else if (error == MB_ERROR_OPEN_FAIL) {
			/* else if unable to create lock file there is a permissions problem */
			GMT_Report(API, GMT_MSG_NORMAL, "Unable to create lock for intended input file:\n");
			GMT_Report(API, GMT_MSG_NORMAL, "  %s\n", Ctrl->I.file);
			GMT_Report(API, GMT_MSG_NORMAL, "-Likely permissions issue\n");
		}
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(-1);
	}

	formatread = format;

	/* initialize reading the input swath sonar file */
	if ((status = mb_read_init(verbose,Ctrl->I.file,formatread,pings,Ctrl->L.lonflip,bounds, btime_i,etime_i,
		speedmin,timegap, &mbio_ptr,&btime_d,&etime_d, &beams_bath,&beams_amp,&pixels_ss,&error)) != MB_SUCCESS) {
		mb_error(verbose,error,&message);
		GMT_Report(API, GMT_MSG_NORMAL, "MBIO Error returned from function <mb_read_init>:\n%s\n",message);
		GMT_Report(API, GMT_MSG_NORMAL, "Multibeam File <%s> not initialized for reading\n",Ctrl->I.file);
		GMT_Report(API, GMT_MSG_NORMAL, "Program <%s> Terminated\n", program_name);
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(error);
	}

	/* initialize and increment counting variables */
	ndata = nzeropos = 0;
	nmin = nbad = nflag = nunflag = nflagesf = nunflagesf = nzeroesf = 0;

	/* give the statistics */
	if (verbose >= 0)
		GMT_Report(API, GMT_MSG_NORMAL, "\nProcessing %s\n",Ctrl->I.file);

	/* allocate memory for data arrays */		// PROBABLY ONLY NEED TO ALLOCATE beamflag & beamflagorig
	ping[0].beamflag = NULL;
	ping[0].beamflagorg = NULL;
	ping[0].bath = NULL;
	ping[0].bathacrosstrack = NULL;
	ping[0].bathalongtrack = NULL;
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), &ping[0].beamflag, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(char), &ping[0].beamflagorg, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), &ping[0].bath, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), &ping[0].bathacrosstrack, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, sizeof(double), &ping[0].bathalongtrack, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_AMPLITUDE, sizeof(double), &amp, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), &ss, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), &ssacrosstrack, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_SIDESCAN, sizeof(double), &ssalongtrack, &error);
	if (error == MB_ERROR_NO_ERROR)
		mb_register_array(verbose, mbio_ptr, MB_MEM_TYPE_BATHYMETRY, 4 * sizeof(double), &list, &error);

	/* if error initializing memory then quit */
	if (error != MB_ERROR_NO_ERROR) {
		mb_error(verbose,error,&message);
		GMT_Report(API, GMT_MSG_NORMAL, "\nMBIO Error allocating data arrays:\n%s\n",message);
		GMT_Report(API, GMT_MSG_NORMAL, "\nProgram <%s> Terminated\n", program_name);
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(error);
	}

	/* now deal with old edit save file */
	GMT_Report(API, GMT_MSG_VERBOSE, "Sorting old edits...\n");
	/* handle esf edits */
	status = mb_esf_load(verbose, program_name, Ctrl->I.file, MB_YES, MB_YES, esffile, &esf, &error);
	if (status == MB_SUCCESS && esf.esffp != NULL)
		esffile_open = MB_YES;
	if (status == MB_FAILURE && error == MB_ERROR_OPEN_FAIL) {
		esffile_open = MB_NO;
		GMT_Report(API, GMT_MSG_NORMAL, "Unable to open new edit save file %s\n", esf.esffile);
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(error);
	}
	else if (status == MB_FAILURE && error == MB_ERROR_MEMORY_FAIL) {
		esffile_open = MB_NO;
		GMT_Report(API, GMT_MSG_NORMAL, "Unable to allocate memory for edits in esf file %s\n", esf.esffile);
		gmt_M_free(GMT, row_flags);		fclose(fid);
		Return(error);
	}
	GMT_Report(API, GMT_MSG_VERBOSE, "%d old edits sorted...\n",esf.nedit);

	/* read */
	done = MB_NO;
	start = 0;
	nrec = 0;
	GMT_Report(API, GMT_MSG_VERBOSE, "Processing data...\n");
	while (done == MB_NO) {
		/* read next record */
		error = MB_ERROR_NO_ERROR;
		status = mb_get(verbose, mbio_ptr,&kind,&pingsread, ping[0].time_i,&ping[0].time_d,
						&ping[0].navlon,&ping[0].navlat, &ping[0].speed,&ping[0].heading,
						&distance,&altitude,&sonardepth, &ping[0].beams_bath,&beams_amp,&pixels_ss,
						ping[0].beamflag,ping[0].bath,amp, ping[0].bathacrosstrack,
						ping[0].bathalongtrack, ss, ssacrosstrack,ssalongtrack, comment, &error);

		if (status == MB_SUCCESS && kind == MB_DATA_DATA) {

			for (i = 0; i < ping[0].beams_bath; i++)		/* save original beamflags */
				ping[0].beamflagorg[i] = ping[0].beamflag[i];

			/* apply saved edits */
			status = mb_esf_apply(verbose, &esf, ping[0].time_d, 0, ping[0].beams_bath, ping[0].beamflag, &error);

			/* update counters */
			for (i = 0; i < ping[0].beams_bath; i++) {
				if (mb_beam_ok(ping[0].beamflag[i]))
					nunflagesf++;
				else
					nflagesf++;
			}
			ndata++;
			nrec++;
		}
		else if (error > MB_ERROR_NO_ERROR)
			done = MB_YES;

		/* process a record */
		if (nrec > 0) {
			if (ping[0].beams_bath != dims[1]) {
				GMT_Report (API, GMT_MSG_NORMAL, "ERROR: Number of beams in file (%d) not equal to n columns (%d) in external file\n", ping[0].beams_bath, dims[1]);
				gmt_M_free(GMT, row_flags);		fclose(fid);
				Return(-1);
			}
			k = (int)fread(row_flags, sizeof(char), (size_t)dims[1], fid);
			for (i = 0; i < dims[1]; i++) {
				if (row_flags[i]) {
					ping[0].beamflag[i] = MB_FLAG_FLAG + MB_FLAG_FILTER;
					mb_ess_save(verbose, &esf, ping[0].time_d, i, MBP_EDIT_FILTER, &error);
					nflagtot++;
				}
			}
		}

		/* write out edits */

		for (i=0;i<ping[0].beams_bath;i++) {
			if (ping[0].beamflag[i] != ping[0].beamflagorg[i]) {
				if (mb_beam_ok(ping[0].beamflag[i]))
					action = MBP_EDIT_UNFLAG;
				else if (mb_beam_check_flag_filter2(ping[0].beamflag[i]))
					action = MBP_EDIT_FILTER;
				else if (mb_beam_check_flag_filter(ping[0].beamflag[i]))
					action = MBP_EDIT_FILTER;
				else if (ping[0].beamflag[i] != MB_FLAG_NULL)
					action = MBP_EDIT_FLAG;
				else
					action = MBP_EDIT_ZERO;
				mb_esf_save(verbose, &esf, ping[0].time_d, i, action, &error);
			}
		}
	}

	gmt_M_free(GMT, row_flags);
	fclose(fid);

	status = mb_close(verbose,&mbio_ptr,&error);	/* close the file */
	status = mb_esf_close(verbose, &esf, &error);	/* close edit save file */

	if (esffile_open == MB_YES) {		/* update mbprocess parameter file */
		status = mb_pr_update_format(verbose, Ctrl->I.file, MB_YES, format, &error);
		status = mb_pr_update_edit(verbose, Ctrl->I.file, MBP_EDIT_ON, esffile, &error);
	}

	/* unlock the raw swath file */
	if (uselockfiles == MB_YES)
		status = mb_pr_unlockswathfile(verbose, Ctrl->I.file, MBP_LOCK_EDITBATHY, program_name, &error);

	/* increment the total counting variables */
	ndatatot += ndata;
	nflagesftot += nflagesf;
	nunflagesftot += nunflagesf;
	nunflagtot += nunflag;

	/* give the total statistics */
	if (verbose >= 0) {
		GMT_Report(API, GMT_MSG_NORMAL, "\nMBclean Processing Totals:\n");
		GMT_Report(API, GMT_MSG_NORMAL, "-------------------------\n");
		GMT_Report(API, GMT_MSG_NORMAL, "%d total bathymetry data records processed\n",ndatatot);
		GMT_Report(API, GMT_MSG_NORMAL, "%d total beams flagged in old esf files\n",nflagesftot);
		GMT_Report(API, GMT_MSG_NORMAL, "%d total beams unflagged in old esf files\n",nunflagesftot);
		GMT_Report(API, GMT_MSG_NORMAL, "%d total beams flagged\n",nflagtot);
	}

	Return (EXIT_SUCCESS);
}
