/*--------------------------------------------------------------------
 *
 *    Coffeeright (c) 2006-2018 by J. Luis
 *
 *--------------------------------------------------------------------
 *
 * GMTMBGRID is an adaptation of the MB-SYSTEM MBGRID program to work with
 * the GMT package. This program uses a gaussian weighted mean algorithm to
 * grid regions covered by data and than fills gaps using a minimum curvature
 * algorithm. Two algorithms are available. The minimum curvature algorithm
 * from the GMT program surface and the IGPP/SIO zgrid routine for thin plate
 * spline interpolation.
 * LOTS of the MBGRID original code was removed (for example, it now works
 * only with xyz data), but it can now use both of the interpolation algorithms
 * on running time. It is also much more memory efficient than the original.
 *
 * With high density of data points the zgrid algorithm seams to be much
 * faster than surface, but extrapolation of large areas can result in
 * bizzare results.
 *
 *--------------------------------------------------------------------
 * MBGRID Author:	D. W. Caress
 * Date:	February 22, 1993
 * Rewrite:	May 2, 1994
 * Rerewrite:	April 25, 1995
 * Rererewrite:	January 2, 1996
 */

/*--------------------------------------------------------------------
 *
 * GMTMBGRID Author:	Joaquim Luis
 * Created:		04-JUN-2006 
 *	Add breakline	08-JAN-2010
 *
 */

/* GMT include files */
#include "gmt_dev.h"

#define THIS_MODULE_CLASSIC_NAME	"mbgmtgrid"
#define THIS_MODULE_MODERN_NAME	"mbgmtgrid"
#define THIS_MODULE_LIB		"mbgmt"
#define THIS_MODULE_PURPOSE	"Grid table data using adjustable tension continuous curvature splines"
#define THIS_MODULE_KEYS	"<D{,DD(,LG(,GG}"
#define THIS_MODULE_NEEDS	"R"
#define THIS_MODULE_OPTIONS "-:RVabdefhirs" GMT_OPT("FH")

EXTERN_MSC int GMT_gmtmbgrid(void *API, int mode, void *args);

struct GMTMBGRID_CTRL {
	struct GRVF_In {
		bool active;
		unsigned int n_files;	/* 1 or 2 */
		char *file[2];
	} In;
	struct GMTMBG_A {	/* -A<aspect_ratio> */
		bool active;
		double value;
	} A;
	struct GMTMBG_C {	/* -C<converge_limit> */
		bool active;
		int clipmode;
		int clip;
		double value;
	} C;
	struct GMTMBG_D {	/* -D<line.xyz> */
		bool active;
		char *breakfile;	/* Name of file with breaklines */
	} D;
	struct GMTMBG_E {	/* -E?? */
		bool active;
		double extend;
	} E;
	struct GMTMBG_F {	/* -E?? */
		bool active;
		char *backgroundfile;
		int grdrasterid;
	} F;
	struct GMTMBG_G {	/* -G<file> */
		bool active;
		char *file;
	} G;
	struct GMTMBG_I {	/* -Idx[/dy] */
		bool active;
		double inc[2];
	} I;
	struct GMTMBG_L {	/* -Ll|u<limit> */
		bool active;
		char *low, *high;
		double min, max;
		int64_t lmode, hmode;
	} L;
	struct GMTMBG_M {	/* -M?? */
		bool active;
		int interp_method;
	} M;
	struct GMTMBG_N {	/* -N<max_iterations> */
		bool active;
		int value;
	} N;
	struct GMTMBG_Q {	/* -Q */
		bool active;
	} Q;
	struct GMTMBG_S {	/* -S<radius>[m|c] */
		bool active;
		double radius;
		char unit;
	} S;
	struct GMTMBG_T {	/* -T<tension>[i][b] */
		bool active;
		double tension;
	} T;
	struct GMTMBG_W {	/* -W<scale> */
		bool active;
		double scale;
	} W;
	struct GMTMBG_Z {	/* -Z<over_relaxation_parameter> */
		bool active;
		double value;
	} Z;
};

struct GMTMBGRID_DATA {	/* Data point and index to node it currently constrains  */
	float x;
	float y;
	float z;
	int64_t index;
};

struct SURFACE_GLOBAL {		/* Things needed inside compare function must be global for now */
	int block_nx;		/* Number of nodes in x-dir for a given grid factor */
	int block_ny;		/* Number of nodes in y-dir for a given grid factor */
	double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
	double x_min, y_min;		/* Lower left corner of grid */
} GMT_Surface_Global;

struct GMTMBGRID_INFO {     /* Control structure for surface setup and execution */
	char *iu;               /* Pointer to grid info array */
	char mode_type[2];		/* D means include data points when iterating
                               I means just interpolate from larger grid */
	char format[128];
	char *low_file, *high_file;	/* Pointers to grids with low and high limits, if selected */
	char *breakfile;        /* Pointer to file with breakline */
	int64_t n_fact;         /* Number of factors in common (ny-1, nx-1) */
	int factors[32];        /* Array of common factors */
	int set_low;            /* 0 unconstrained,1 = by min data value, 2 = by user value */
	int set_high;           /* 0 unconstrained,1 = by max data value, 2 = by user value */
	uint64_t npoints;       /* Number of data points */
	int64_t ij_sw_corner, ij_se_corner,ij_nw_corner, ij_ne_corner;
	uint64_t n_empty;        /* No of unconstrained nodes at initialization  */
	uint64_t nbinset;
	uint64_t *index;	/* EXPERIENCIA */
	int nx;             /* Number of nodes in x-dir. */
	int ny;             /* Number of nodes in y-dir. (Final grid) */
	uint64_t nxny;			/* Total number of grid nodes without boundaries  */
	int mx;
	int my;
	uint64_t mxmy;         /* Total number of grid nodes with boundaries  */
	int block_nx;          /* Number of nodes in x-dir for a given grid factor */
	int block_ny;          /* Number of nodes in y-dir for a given grid factor */
	int max_iterations;    /* Max iter per call to iterate */
	int total_iterations;
	int interp_method;     /* Either 'surface' or 'zgrid' */
	unsigned short int *num;  /* */
	int offx, offy;        /* Number of cell sizes extensible with -E option */
	int gxdim, gydim;      /* Extended dimensions (use offx,offy) */
	int grid, old_grid;    /* Node spacings for the 'surface' branch */
	int64_t grid_east;
	int64_t offset[25][12];	/* Indices of 12 nearby points in 25 cases of edge conditions  */
	int constrained;		/* TRUE if set_low or set_high is TRUE */
	float *lower, *upper;		/* arrays for minmax values, if set */
	float *xcoords, *ycoords;
	float *u;
	float *grid_aux;
	double extend, scale, clipvalue;
	double low_limit, high_limit;	/* Constrains on range of solution */
	double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
	double r_grid_xinc, r_grid_yinc;	/* Reciprocals  */
	double converge_limit;		/* Convergence limit */
	double radius;			/* Search radius for initializing grid  */
	double tension;         /* For the zgrid algo */
	double boundary_tension;
	double interior_tension;
	double a0_const_1, a0_const_2;	/* Constants for off grid point equation  */
	double e_2, e_m2, one_plus_e2;
	double eps_p2, eps_m2, two_plus_ep2, two_plus_em2;
	double x_edge_const, y_edge_const;
	double l_epsilon;
	double z_mean;
	double z_scale;			/* Root mean square range of z after removing planar trend  */
	double r_z_scale;		/* reciprocal of z_scale  */
	double plane_c0, plane_c1, plane_c2;	/* Coefficients of best fitting plane to data  */
	double small;			/* Let data point coincide with node if distance < C->small */
	double coeff[2][12];		/* Coefficients for 12 nearby points, constrained and unconstrained  */
	double relax_old, relax_new;	/* Coefficients for relaxation factor to speed up convergence */
	double wesn_orig[4];		/* Original -R domain as we might have shifted it due to -r */
	struct SURFACE_DATA  *data;
	struct SURFACE_BRIGGS *briggs;
	struct GMT_GRID *Grid;			/* The final grid */
	struct GMT_GRID *Low, *High;		/* arrays for minmax values, if set */
};

#define GMT_WITH_NO_PS

/* flag for no data in grid */
#define	NO_DATA_FLAG	99999

/* interpolation mode */
#define MBGRID_INTERP_NONE	0
#define MBGRID_INTERP_GAP	1
#define MBGRID_INTERP_NEAR	2
#define MBGRID_INTERP_ALL	3

#define INTERP_SURFACE		1
#define INTERP_ZGRID		0

#define OUTSIDE 2000000000	/* Index number indicating data is outside useable area */

#define SURFACE_OUTSIDE LONG_MAX	/* Index number indicating data is outside usable area */

/* Variables for surface */
struct SURFACE_DATA {
	float x;
	float y;
	float z;
	int index;
};

struct SURFACE_BRIGGS {
	double b[6];
};
 

int mb_zgrid(struct GMTMBGRID_INFO *C, int verbose, float *z, uint64_t n, float *zpij, int *knxt, int *imnew, int nrng);
int mb_surface(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, uint64_t ndat, float *sgrid);
int interp_breakline(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, float *norm, unsigned short int *cnt, double *wbnd, 
		struct GMT_DATATABLE *xyzline, double factor, int xtradim, int gxdim, int gydim);
int read_data(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, int ndat, float *xdat, float *ydat, float *zdat);
int	gcd_euclid(int a, int b);	/* Finds the greatest common divisor  */
int	get_prime_factors(int n, int *f);
int iterate(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, int mode);
int compare_points (const void *point_1v, const void *point_2v);
void set_grid_parameters(struct GMTMBGRID_INFO *C);
void throw_away_unusables(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C);
void remove_planar_trend(struct GMTMBGRID_INFO *C);
int  rescale_z_values(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C);
void load_constraints(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C);
void smart_divide(struct GMTMBGRID_INFO *C);
void set_offset(struct GMTMBGRID_INFO *C);
void set_index(struct GMTMBGRID_INFO *C);
void initialize_grid(struct GMTMBGRID_INFO *C);
void set_coefficients(struct GMTMBGRID_INFO *C);
void find_nearest_point(struct GMTMBGRID_INFO *C);
void fill_in_forecast(struct GMTMBGRID_INFO *C);
void check_errors(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C);
void replace_planar_trend(struct GMTMBGRID_INFO *C);
void new_initialize_grid(struct GMTMBGRID_INFO *C);
void get_output(struct GMTMBGRID_INFO *C, float *sgrid);
int  read_data_surface (struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, struct GMT_OPTION *options);
double guess_surface_time(struct GMTMBGRID_INFO *C);

static void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct GMTMBGRID_CTRL *C = NULL;
	
	C = gmt_M_memory (GMT, NULL, 1, struct GMTMBGRID_CTRL);

	C->C.clipmode = MBGRID_INTERP_ALL;
	C->M.interp_method = INTERP_SURFACE;

	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->N.value = 250;
	C->A.value = 1.0;
	C->W.scale = 1.0;
	C->Z.value = 1.4;
		
	return (C);
}

static void Free_Ctrl (struct GMT_CTRL *GMT, struct GMTMBGRID_CTRL *C) {	/* Deallocate control structure */
	if (!C) return;
	if (C->G.file) free (C->G.file);	
	if (C->D.breakfile) free (C->D.breakfile);	
	if (C->L.low)  free (C->L.low);	
	if (C->L.high) free (C->L.high);	
	free (C);	
}

static int usage (struct GMTAPI_CTRL *API, int level) {

	gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);

	GMT_Message (API, GMT_TIME_NONE, "usage: gmtmbgrid [<table>] -G<outgrid> %s\n", GMT_I_OPT);
	GMT_Option (API, "C,3,f,h,i,F,:,.");
	GMT_Message (API, GMT_TIME_NONE, "\t%s [-C<clip>[g|n|o]] [-D<breakline>]\n", GMT_Rgeo_OPT);
	GMT_Message (API, GMT_TIME_NONE, "\t[-E<extend>] [-F<background_grid>] [-Ll<limit>] [-Lu<limit>] [-T<tension>]\n");
	GMT_Message (API, GMT_TIME_NONE, "\t[-W<scale>] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT);

	if (level == GMT_SYNOPSIS) return (GMT_MODULE_SYNOPSIS);

	GMT_Message (API, GMT_TIME_NONE, "\tgmtmbgrid will read from standard input or <xyz-file[s]>.\n");
	GMT_Message (API, GMT_TIME_NONE, "\tNotice that ONLY the first file is used in the main interpolation>.\n");
	GMT_Message (API, GMT_TIME_NONE, "\tA second file, if provided, will be used as 'background'.\n\n");
	GMT_Message (API, GMT_TIME_NONE, "\tRequired arguments to gmtmbgrid:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-G sets output grd File name\n");
	GMT_Option (API, "I,R");
	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-C<clip> clip value sets the distance from data (in grid cells) that\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tthe spline interpolation may be applied.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tAppend n to fill all undefined cells within a distance of clip cells from data.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tAppend g to fill data gaps up to two times clip grid cells in size.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tUse -Co to not apply spline interpolation.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tIf no letter is appended it defaults to n.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tDefault (no -C) is to interpolate all grid cells.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-D<breakline> uses XYZ data (can be multiseg) in breakline file as a 'soft breakline'.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-E<extend> uses data extending extend*nx|ny in the interpolation.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-F<background_grid> Enables filling in all undefined grid cells with\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tdata extracted from the <background_grid>. Alternatively, give a\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tsecond xyz-file in the command line. This will than be interpolated before use.\n");
	GMT_Option (API, "H");
	GMT_Message (API, GMT_TIME_NONE, "\t-M<s|z> Select interpolation algorithm. -Mz -> zgrid. -Ms -> surface [Default]\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-T adds Tension to the gridding equation; use a value between 0 and 1.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tOr a value between  0 and Inf if the -Mz was used.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tdefault = 0 gives minimum curvature (smoothest; bicubic) solution.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-W<scale> Sets the width of the gaussian weighting function in terms of\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tthe grid spacing [Default: scale = 1.0].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t\tDefault is 3 input columns.\n\n");
	GMT_Option (API, "a,bi3,f,h,i,r,s,:,.");

	return (EXIT_FAILURE);
}

static int parse (struct GMT_CTRL *GMT, struct GMTMBGRID_CTRL *Ctrl, struct GMT_OPTION *options) {
	/* This parses the options provided to surface and sets parameters in CTRL.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	int n_errors = 0;
	size_t len;
	char fill[16];
	struct GMT_OPTION *opt = NULL;

	for (opt = options; opt; opt = opt->next) {
		switch (opt->option) {

			case '<':	/* Input file (if) */
				Ctrl->In.active = true;
				if (Ctrl->In.n_files < 2) 
					Ctrl->In.file[Ctrl->In.n_files++] = strdup (opt->arg);
				else {
					n_errors++;
					GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Error: A maximum of two input files may be processed\n");
				}
				break;

			/* Processes program-specific parameters */
			case 'C':
				Ctrl->C.active = true;
				len = strlen(opt->arg);
				strcpy (fill, opt->arg);
				if (opt->arg[len-1] == 'O' || opt->arg[len-1] == 'o')
					Ctrl->C.clipmode = MBGRID_INTERP_NONE;
				else if (opt->arg[len-1] == 'N' || opt->arg[len-1] == 'n') {
					fill[len-1] = '\0';
					Ctrl->C.clip = atoi(fill);
					Ctrl->C.clipmode = MBGRID_INTERP_NEAR;
				}
				else if (opt->arg[len-1] == 'G' || opt->arg[len-1] == 'g') {
					fill[len-1] = '\0';
					Ctrl->C.clip = atoi(fill);
					Ctrl->C.clipmode = MBGRID_INTERP_GAP;
				}
				else {
					Ctrl->C.clip = atoi(&opt->arg[1]);
					Ctrl->C.clipmode = MBGRID_INTERP_GAP;
				}
				break;
			case 'D':
				Ctrl->D.active = true;
				Ctrl->D.breakfile = strdup (opt->arg);
				break;
			case 'E':
				Ctrl->E.active = true;
				Ctrl->E.extend = atof (opt->arg);
				break;
			case 'F':
				Ctrl->F.active = true;
				Ctrl->F.backgroundfile = strdup (opt->arg);
				break;
			case 'G':
				Ctrl->G.active = true;
				Ctrl->G.file = strdup (opt->arg);
				break;
			case 'I':
				Ctrl->I.active = true;
				n_errors += gmt_parse_inc_option (GMT, 'I', opt->arg);
				if (gmt_getinc (GMT, opt->arg, Ctrl->I.inc)) {
					gmt_inc_syntax (GMT, 'I', 1);
					n_errors++;
				}
				break;
			case 'M':
				if (opt->arg[0] == 'Z' || opt->arg[0] == 'z')
					Ctrl->M.interp_method = INTERP_ZGRID;
				else if (opt->arg[0] == 'S' || opt->arg[0] == 's')
					Ctrl->M.interp_method = INTERP_SURFACE;
				else
					fprintf(stderr,"gmtmbgrid: GMT SYNTAX ERROR -M option: Trash given - Ignored\n");
				break;

			case 'L':	/* Set limits */
				Ctrl->L.active = true;
				switch (opt->arg[0]) {
					case 'l':	/* Lower limit  */
						n_errors += gmt_M_check_condition (GMT, opt->arg[1] == 0, "Syntax error -Ll option: No argument given\n");
						Ctrl->L.low = strdup (&opt->arg[1]);
						if (!gmt_access (GMT, Ctrl->L.low, R_OK))	/* file exists */
							Ctrl->L.lmode = 3;
						else if (Ctrl->L.low[0] == 'd')		/* Use data minimum */
							Ctrl->L.lmode = 1;
						else {
							Ctrl->L.lmode = 2;		/* Use given value */
							Ctrl->L.min = atof (&opt->arg[1]);
						}
						break;
					case 'u':	/* Upper limit  */
						n_errors += gmt_M_check_condition (GMT, opt->arg[1] == 0, "Syntax error -Lu option: No argument given\n");
						Ctrl->L.high = strdup (&opt->arg[1]);
						if (!gmt_access (GMT, Ctrl->L.high, R_OK))	/* file exists */
							Ctrl->L.hmode = 3;
						else if (Ctrl->L.high[0] == 'd')	/* Use data maximum */
							Ctrl->L.hmode = 1;
						else {
							Ctrl->L.hmode = 2;		/* Use given value */
							Ctrl->L.max = atof (&opt->arg[1]);
						}
						break;
					default:
						n_errors++;
						break;
				}
				break;

			case 'T':
				Ctrl->T.active = true;
				Ctrl->T.tension = atof (opt->arg);
				break;
			case 'W':
				Ctrl->W.active = true;
				Ctrl->W.scale = atof (opt->arg);
				break;

			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, !GMT->common.R.active, "Syntax error: Must specify -R option\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->D.active && Ctrl->D.breakfile && gmt_access (GMT, Ctrl->D.breakfile, R_OK), "Syntax error -D: Cannot read file %s!\n", Ctrl->D.breakfile);
	n_errors += gmt_M_check_condition (GMT, Ctrl->I.inc[GMT_X] <= 0.0 || Ctrl->I.inc[GMT_Y] <= 0.0,
	                                   "Syntax error -I option: Must specify positive increment(s)\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->N.value < 1, "Syntax error -N option: Max iterations must be nonzero\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->Z.value < 1.0 || Ctrl->Z.value > 2.0,
	                                   "Syntax error -Z option: Relaxation value must be 1 <= z <= 2\n");
	n_errors += gmt_M_check_condition (GMT, !Ctrl->G.file, "Syntax error option -G: Must specify output file\n");
	n_errors += gmt_check_binary_io (GMT, 3);

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}


/* Convenience macros to free memory before exiting due to error or completion */
#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

int GMT_gmtmbgrid (void *V_API, int mode, void *args) {

	uint64_t i, j, ii, jj, iii, jjj, kkk, ir, i1, i2, j1, j2;
	uint64_t kgrid, kout, kint;
	int	dmask[9], n_files = 0;
	int error = 0, one = 1;
	int	nbinzero, nbinspline, nbinbackground = 0;
	int	use_NaN = true;
	int	do_break = FALSE;
	int	*work2 = NULL, *work3 = NULL;
	size_t n_alloc = GMT_CHUNK;
	char *backgroundfile = NULL, *grdfile = NULL, *breakfile = NULL;

	float	outclipvalue = NO_DATA_FLAG;
	float	*work1 = NULL, *output = NULL, *sgrid = NULL;

	double wesn[4];
	double clipvalue = NO_DATA_FLAG;
	double r;
	double zflag = 5.0e34;
	double zmin, zmax, zclip;
	
	struct GMT_TABLE *xyzline = NULL;
	struct GMT_DATASET *Lin = NULL;
	struct GMTMBGRID_INFO C;
	struct GMTMBGRID_CTRL *Ctrl = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_OPTION *options = NULL;
	struct GMT_GRID_HEADER *h = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */
	/* =========================================================================---------- */

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage (API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMT_OPT_USAGE) bailout (usage (API, GMT_USAGE));	/* Return the usage message */
	if (options->option == GMT_OPT_SYNOPSIS) bailout (usage (API, GMT_SYNOPSIS));	/* Return the synopsis */

	/* Parse the commont GMT command-line options */
	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, NULL, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
	if (GMT_Parse_Common (API, THIS_MODULE_OPTIONS, options)) Return (API->error);
	Ctrl = New_Ctrl (GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse (GMT, Ctrl, options)) != 0) Return (error);

	/*----------------------- Standard module initialization and parsing ----------------------*/

	
	/*---------------------------- This is the gmtmbgrid main code ----------------------------*/

	gmt_M_memset (&C, 1, struct GMTMBGRID_INFO);
	gmt_M_memset (&GMT_Surface_Global, 1, struct SURFACE_GLOBAL);

	C.converge_limit = Ctrl->C.value;
	C.extend = Ctrl->E.extend;
	C.interp_method = Ctrl->M.interp_method;
	C.max_iterations = Ctrl->N.value;
	C.radius = Ctrl->S.radius;
	C.tension = Ctrl->T.tension;
	C.scale = Ctrl->W.scale;
	C.relax_new = Ctrl->Z.value;
	C.clipvalue = clipvalue;		/* HARDWIRED TO ---  NO_DATA_FLAG --- */
	C.mode_type[0] = 'I';
	C.mode_type[1] = 'D';	/* D means include data points when iterating */
	C.breakfile = Ctrl->D.breakfile;
	C.z_scale = C.r_z_scale = 1;

	/* Currently in 'surface.c'
	C.low_file = Ctrl->L.low;
	C.high_file = Ctrl->L.high;
	C.set_low = Ctrl->L.lmode;
	C.low_limit = Ctrl->L.min;
	C.set_high = Ctrl->L.hmode;
	C.high_limit = Ctrl->L.max;
	C.boundary_tension = Ctrl->T.b_tension;
	C.interior_tension = Ctrl->T.i_tension;
	C.l_epsilon = Ctrl->A.value;
	*/

	GMT = API->GMT;

	gmt_M_memcpy (wesn, GMT->common.R.wesn, 4, double);         /* Specified region */
	gmt_M_memcpy (C.wesn_orig, GMT->common.R.wesn, 4, double);  /* Save original region in case of -r */

	if (GMT->common.R.active[GSET]) {		/* Pixel registration request. Use the trick of offsetting area by x_inc(y_inc) / 2 */
		wesn[XLO] += Ctrl->I.inc[GMT_X] / 2.0;	wesn[XHI] += Ctrl->I.inc[GMT_X] / 2.0;
		wesn[YLO] += Ctrl->I.inc[GMT_Y] / 2.0;	wesn[YHI] += Ctrl->I.inc[GMT_Y] / 2.0;
		one++;	/* Just so we can report correct nx,ny for the grid; internally it is the same until output */
		/* nx,ny remains the same for now but nodes are in "pixel" position.  Must reduce nx,ny by 1 when we write result */
	}

	gmt_set_pad (GMT, 0U);	/* Turn off padding (and thus BC setting) */
	if ((C.Grid = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, wesn, NULL, 
		GMT_GRID_NODE_REG, GMT_NOTSET, NULL)) == NULL) Return (API->error);

	gmt_RI_prepare (GMT, C.Grid->header);	/* Ensure -R -I consistency and set nx, ny */
	gmt_M_err_fail (GMT, gmt_grd_RI_verify (GMT, C.Grid->header, 1), Ctrl->G.file);
	gmt_set_grddim (GMT, C.Grid->header);

	if (C.Grid->header->n_columns < 4 || C.Grid->header->n_rows < 4) {
		GMT_Report (API, GMT_MSG_NORMAL, "Error: Grid must have at least 4 nodes in each direction (you have %d by %d) - abort.\n", C.Grid->header->n_columns, C.Grid->header->n_rows);
		Return (EXIT_FAILURE);
	}

	C.nx = C.Grid->header->n_columns;
	C.ny = C.Grid->header->n_rows;
	C.nxny = C.Grid->header->nm;

	C.mx = C.nx + 4;
	C.my = C.ny + 4;
	C.mxmy = C.Grid->header->size;
	GMT_Surface_Global.x_min = C.Grid->header->wesn[XLO];
	GMT_Surface_Global.y_min = C.Grid->header->wesn[YLO];

	h = C.Grid->header;

	/* Since the mbgrid version used as base for gmtmbgrid was bugged regarding the
	   option of interpolating all nodes in the grid, I use this as a patch */
	if (Ctrl->C.clipmode == MBGRID_INTERP_ALL)
		Ctrl->C.clip = MAX(C.nx, C.ny) + 1;

	/* User selected a background file/grid but didn't specify -C option */
	if (Ctrl->C.clipmode == MBGRID_INTERP_ALL && (Ctrl->F.backgroundfile != NULL))
		Ctrl->C.clipmode = MBGRID_INTERP_NEAR;

	if (use_NaN == true)
		outclipvalue = GMT->session.f_NaN;

	/* check interpolation parameters */
	if ((Ctrl->C.clipmode == MBGRID_INTERP_GAP || Ctrl->C.clipmode == MBGRID_INTERP_NEAR) && 
		 Ctrl->C.clip > C.nx && Ctrl->C.clip > C.ny)
		Ctrl->C.clipmode = MBGRID_INTERP_ALL;

	if (read_data_surface (GMT, &C, options)) Return (API->error);

	nbinspline = 0;

	/* if clip set do smooth interpolation */
	if (Ctrl->C.clipmode != MBGRID_INTERP_NONE && Ctrl->C.clip > 0 && C.nbinset > 0) {
		/* set up data vector */
		sgrid = gmt_M_memory (GMT, NULL, (size_t)(C.gxdim * C.gydim), float);

		if (Ctrl->M.interp_method == INTERP_SURFACE) {
			/* do the interpolation */
			GMT_Report (API, GMT_MSG_VERBOSE, "\nDoing Surface spline interpolation with %d data points...\n",C.npoints);
			mb_surface(GMT, &C, C.nbinset, sgrid);
			gmt_M_free (GMT, C.xcoords);  gmt_M_free (GMT, C.ycoords);
		}
		else {		/* ZGRID */
			work1 = gmt_M_memory (GMT, NULL, (size_t)(C.nbinset), float);
			work2 = gmt_M_memory (GMT, NULL, (size_t)(C.nbinset), int);
			work3 = gmt_M_memory (GMT, NULL, (size_t)(C.gxdim + C.gydim), int);

			/* do the interpolation */
			GMT_Report (API, GMT_MSG_VERBOSE, "\nDoing Zgrid spline interpolation with %d data points...\n",C.nbinset);
			mb_zgrid(&C,(int)gmt_M_is_verbose (GMT, GMT_MSG_VERBOSE), sgrid, C.nbinset, work1, work2, work3, Ctrl->C.clip);
			gmt_M_free (GMT, work1);		gmt_M_free (GMT, work2);	gmt_M_free (GMT, work3);
			gmt_M_free (GMT, C.index);
		}

		if (Ctrl->C.clipmode == MBGRID_INTERP_GAP)
			GMT_Report (API, GMT_MSG_VERBOSE, "Applying spline interpolation to fill gaps of %d cells or less...\n",Ctrl->C.clip);
		else if (Ctrl->C.clipmode == MBGRID_INTERP_NEAR)
			GMT_Report (API, GMT_MSG_VERBOSE, "Applying spline interpolation to fill %d cells from data...\n",Ctrl->C.clip);
		else if (Ctrl->C.clipmode == MBGRID_INTERP_ALL)
			GMT_Report (API, GMT_MSG_VERBOSE, "Applying spline interpolation to fill all undefined cells in the grid...\n");

		/* translate the interpolation into the grid array filling only data gaps */
		if (Ctrl->C.clipmode == MBGRID_INTERP_GAP) {
			for (i = 0; i < C.gxdim; i++) {
				for (j = 0; j < C.gydim; j++) {
					kgrid = i * C.gydim + j;
					if (Ctrl->M.interp_method == INTERP_SURFACE)
						kint = i + (C.gydim -j - 1) * C.gxdim;
					else
						kint = i + j*C.gxdim;

					C.num[kgrid] = 0;
					if (C.grid_aux[kgrid] >= clipvalue && sgrid[kint] < zflag) {
						/* initialize direction mask of search */
						for (ii = 0; ii < 9; ii++)
							dmask[ii] = 0;

						/* loop over rings around point, starting close */
						for (ir = 0; ir <= Ctrl->C.clip && C.num[kgrid] == 0; ir++) {
							/* set bounds of search */
							i1 = MAX(0, i - ir);
							i2 = MIN(C.gxdim - 1, i + ir);
							j1 = MAX(0, j - ir);
							j2 = MIN(C.gydim - 1, j + ir);
							jj = j1;
							for (ii = i1; ii <= i2 && C.num[kgrid] == 0;ii++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = lrint((ii - i)/r) + 1;
									jjj = lrint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = 1;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										C.num[kgrid] = 1;
								}
							}
				      
							jj = j2;
							for (ii = i1; ii <= i2 && C.num[kgrid] == 0;ii++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = lrint((ii - i)/r) + 1;
									jjj = lrint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = 1;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										C.num[kgrid] = 1;
								}
							}
				      
							ii = i1;
							for (jj = j1; jj <= j2 && C.num[kgrid] == 0; jj++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = lrint((ii - i)/r) + 1;
									jjj = lrint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = 1;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										C.num[kgrid] = 1;
								}
							}
				      
							ii = i2;
							for (jj = j1; jj <= j2 && C.num[kgrid] == 0; jj++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = lrint((ii - i)/r) + 1;
									jjj = lrint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = 1;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										C.num[kgrid] = 1;
								}
							}
						}
					}	/* end if */
				}
			}
			for (i = 0; i < C.gxdim; i++) {
				for (j = 0; j < C.gydim; j++) {
					kgrid = i * C.gydim + j;
					if (Ctrl->M.interp_method == INTERP_SURFACE)
						kint = i + (C.gydim -j - 1) * C.gxdim;
					else
						kint = i + j*C.gxdim;
					if (C.num[kgrid] == 1) {
						C.grid_aux[kgrid] = sgrid[kint];
						nbinspline++;
					}
				}
			}
		}

		/* translate the interpolation into the grid array filling by proximity */
		else if (Ctrl->C.clipmode == MBGRID_INTERP_NEAR) {
			for (i = 0; i < C.gxdim; i++) {
				for (j = 0; j < C.gydim; j++) {
					kgrid = i * C.gydim + j;
					if (Ctrl->M.interp_method == INTERP_SURFACE)
						kint = i + (C.gydim -j - 1) * C.gxdim;
					else
						kint = i + j*C.gxdim;

					C.num[kgrid] = 0;
					if (C.grid_aux[kgrid] >= clipvalue && sgrid[kint] < zflag) {
						/* loop over rings around point, starting close */
						for (ir = 0; ir <= Ctrl->C.clip && C.num[kgrid] == 0; ir++) {
							/* set bounds of search */
							i1 = MAX(0, i - ir);
							i2 = MIN(C.gxdim - 1, i + ir);
							j1 = MAX(0, j - ir);
							j2 = MIN(C.gydim - 1, j + ir);
				      
							jj = j1;
							for (ii = i1; ii <= i2 && C.num[kgrid] == 0;ii++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue)
									C.num[kgrid] = 1;
							}
				      
							jj = j2;
							for (ii = i1; ii <= i2 && C.num[kgrid] == 0;ii++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue)
									C.num[kgrid] = 1;
							}
				      
							ii = i1;
							for (jj = j1; jj <= j2 && C.num[kgrid] == 0;jj++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue)
									C.num[kgrid] = 1;
							}
				      
							ii = i2;
							for (jj = j1; jj <= j2 && C.num[kgrid] == 0;jj++) {
								if (C.grid_aux[ii*C.gydim+jj] < clipvalue)
									C.num[kgrid] = 1;
							}
						}
					}
				}
			}

			for (i = 0; i < C.gxdim; i++) {
				for (j = 0; j < C.gydim; j++) {
					kgrid = i * C.gydim + j;
					if (Ctrl->M.interp_method == INTERP_SURFACE)
						kint = i + (C.gydim -j - 1) * C.gxdim;
					else
						kint = i + j*C.gxdim;
					if (C.num[kgrid] == 1) {
						C.grid_aux[kgrid] = sgrid[kint];
						nbinspline++;
					}
				}
			}
		}

		/* translate the interpolation into the grid array filling all empty bins */
		else {
			for (i = 0; i < C.gxdim; i++) {
				for (j = 0; j < C.gydim; j++) {
					kgrid = i * C.gydim + j;
					if (Ctrl->M.interp_method == INTERP_SURFACE)
						kint = i + (C.gydim -j - 1) * C.gxdim;
					else
						kint = i + j*C.gxdim;

					if (C.grid_aux[kgrid] >= clipvalue && sgrid[kint] < zflag) {
						C.grid_aux[kgrid] = sgrid[kint];
						nbinspline++;
					}
				}
			}
		}

		/* deallocate the interpolation arrays */
		gmt_M_free (GMT, sgrid);
	}
	gmt_M_free (GMT, C.num);

	/* Extract background data and interpolate it later onto internal grid */
	if (Ctrl->F.backgroundfile) {

		int  object_ID;
		char in_string[GMT_LEN128], out_string[GMT_LEN128], cmd[GMT_LEN128];
		struct GMT_GRID *bg_grid = NULL;
		struct GMT_GRID *G2 = NULL;

		GMT_Report (GMT, GMT_MSG_VERBOSE, "\nExtracting background from grid file %s...\n",backgroundfile);

		if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, wesn, Ctrl->F.backgroundfile, bg_grid) == NULL)
			Return (API->error);	/* Get grid data */

		/* Create option list, register Background grid as input source via reference */
		if ((object_ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_REFERENCE|GMT_IO_RESET, GMT_IS_SURFACE, GMT_IN, NULL, bg_grid)) == GMT_NOTSET) 
			Return (API->error);

		/* MUST USE Virtual files instead. This was broken by GMT6.1 */
		//if (api_encode_id (API, in_string, object_ID) != GMT_OK)
			//Return (API->error);	/* Make filename with embedded object ID for grid G */

		if ((object_ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_REFERENCE, GMT_IS_SURFACE, GMT_OUT, NULL, NULL)) == GMT_NOTSET)
			Return (API->error);

		//if (api_encode_id (API, out_string, object_ID) != GMT_OK)
			//Return (API->error);	/* Make filename with embedded object ID for result grid G2 */

		sprintf (cmd, "%s -G%s -I%d+/%d+", in_string, out_string, C.nx, C.ny);
		if (GMT_grdsample (API, 0, cmd) != GMT_OK) return (API->error);	/* Do the resampling */
		//if ((G2 = api_retrieve_data (API, object_ID)) == NULL)
			//Return (API->error);

		if (GMT_Destroy_Data (API, &bg_grid) != GMT_OK)
			Return (API->error);

		bg_grid = G2;

		/* translate the interpolation into the grid array (only into the data gaps) */
		for (i = 0; i < C.gxdim; i++) {
			for (j = 0; j < C.gydim; j++) {
				kgrid = i * C.gydim + j;
				if (C.grid_aux[kgrid] >= clipvalue) {
					C.grid_aux[kgrid] = G2->data[kgrid];
					nbinbackground++;
				}
			}
		}
			
		if (GMT_Destroy_Data (API, &bg_grid) != GMT_OK)
			Return (API->error);
	}

	/* if ... */
	if (Ctrl->In.n_files == 2 ) {

		struct GMTMBGRID_INFO C2;
	
		gmt_M_memset (&C2, 1, struct GMTMBGRID_INFO);
		gmt_M_memcpy (&C2, &C, 1, struct GMTMBGRID_INFO);		/* Make a copy of this struct */

		GMT_Report (API, GMT_MSG_NORMAL, "Extracting background from file ...\n");

		options->arg = Ctrl->In.file[1];
		if (read_data_surface (GMT, &C2, options)) Return (API->error);

		/* allocate grid and work arrays */
		sgrid = gmt_M_memory (GMT, NULL, (size_t)(C.gxdim * C.gydim), float);
		if (Ctrl->M.interp_method == INTERP_ZGRID) {
			work1 = gmt_M_memory (GMT, NULL, (size_t)(C2.nbinset), float);
			work2 = gmt_M_memory (GMT, NULL, (size_t)(C2.nbinset), int);
			work3 = gmt_M_memory (GMT, NULL, (size_t)(C.gxdim + C.gydim), int);
		}

		/* do the interpolation */
		GMT_Report (API, GMT_MSG_VERBOSE, "\nDoing spline interpolation with %d data points from background...\n",C2.nbinset);
		if (Ctrl->M.interp_method == INTERP_SURFACE) {
			mb_surface(GMT, &C2, C2.nbinset, sgrid);
			gmt_M_free (GMT, C2.xcoords);  gmt_M_free (GMT, C2.ycoords);
			gmt_M_free (GMT, C2.grid_aux);
		}
		else {
			Ctrl->C.clip = MAX(C.gxdim,C.gydim);
			mb_zgrid(&C2,(int)gmt_M_is_verbose (GMT, GMT_MSG_VERBOSE), sgrid, C2.nbinset, work1, work2, work3, Ctrl->C.clip);
			gmt_M_free (GMT, work1);		gmt_M_free (GMT, work2);	gmt_M_free (GMT, work3);
			gmt_M_free (GMT, C2.index);
		}

		gmt_M_free (GMT, C2.num);

		/* translate the interpolation into the grid array - interpolate only to fill a data gap */
		for (i = 0; i < C.gxdim; i++) {
			for (j = 0; j < C.gydim; j++) {
				kgrid = i * C.gydim + j;
				if (Ctrl->M.interp_method == INTERP_SURFACE)
					kint = i + (C.gydim -j - 1) * C.gxdim;
				else
					kint = i + j*C.gxdim;
				if (C.grid_aux[kgrid] >= clipvalue && sgrid[kint] < zflag) {
					C.grid_aux[kgrid] = sgrid[kint];
					nbinbackground++;
				}
			}
		}
		gmt_M_free (GMT, sgrid);
	}

	/* get min max of data */
	zclip = clipvalue;
	zmin = zmax = zclip;
	for (i = 0; i < C.gxdim; i++)
		for (j = 0; j < C.gydim; j++) {
			kgrid = i * C.gydim + j;
			if (zmin == zclip && C.grid_aux[kgrid] < zclip)
				zmin = C.grid_aux[kgrid];
			if (zmax == zclip && C.grid_aux[kgrid] < zclip)
				zmax = C.grid_aux[kgrid];
			if (C.grid_aux[kgrid] < zmin && C.grid_aux[kgrid] < zclip)
				zmin = C.grid_aux[kgrid];
			if (C.grid_aux[kgrid] > zmax && C.grid_aux[kgrid] < zclip)
				zmax = C.grid_aux[kgrid];
		}
	if (zmin == zclip) zmin = 0.0;
	if (zmax == zclip) zmax = 0.0;

	nbinzero = (int)(C.gxdim*C.gydim - C.nbinset - nbinspline - nbinbackground);
	GMT_Report (API, GMT_MSG_VERBOSE, "\nTotal number of bins:          %d\n",C.gxdim * C.gydim);
	GMT_Report (API, GMT_MSG_VERBOSE, "Bins set using data:             %d\n",C.nbinset);
	GMT_Report (API, GMT_MSG_VERBOSE, "Bins set using interpolation:    %d\n",nbinspline);
	GMT_Report (API, GMT_MSG_VERBOSE, "Bins set using background:       %d\n",nbinbackground);
	GMT_Report (API, GMT_MSG_VERBOSE, "Bins not set:                    %d\n",nbinzero);
	GMT_Report (API, GMT_MSG_VERBOSE, "Minimum value: %10.2f   Maximum value: %10.2f\n", zmin,zmax);

	/* write output file */
	output = gmt_M_memory (GMT, NULL, (size_t)(h->n_columns * h->n_rows), float);
	for (i = 0; i < h->n_columns; i++) {
		for (j = 0; j < h->n_rows; j++) {
			kgrid = (i + C.offx) * C.gydim + (j + C.offy);
			kout = (h->n_rows - 1 - j) * h->n_columns + i;
			output[kout] = (float) C.grid_aux[kgrid];
			if (C.grid_aux[kgrid] == clipvalue)
				output[kout] = outclipvalue;
		}
	}

	C.Grid->data = output;
	gmt_set_pad (GMT, 0U);	/* Turn off padding (and thus BC setting) */

	if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_AND_DATA, NULL, Ctrl->G.file, C.Grid) != GMT_NOERROR)
		Return (API->error);

	/* deallocate arrays */
	gmt_M_free (GMT, C.grid_aux);

	GMT_Report (API, GMT_MSG_VERBOSE, "done!\n");

	Return (EXIT_SUCCESS);
}

/* --------------------------------------------------------------------------------------------- */
int read_data_surface (struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, struct GMT_OPTION *options) {

	int i, j, ii, jj, ix, iy, ix1, iy1, ix2, iy2, xtradim, gxdim, gydim, error;
	unsigned short int *cnt;
	uint64_t k = 0, kmax = 0, kmin = 0, kgrid;
	float *norm;
	double factor, xx, yy, xx2, weight;
	double *in, zmin = DBL_MAX, zmax = -DBL_MAX, bounds[4], wbnd[4];
	struct GMT_RECORD *In = NULL;

	/* calculate other grid properties */
	factor = 4.0 / (C->scale * C->scale * C->Grid->header->inc[GMT_X] * C->Grid->header->inc[GMT_Y]);
	if (C->extend > 0) {		/* Otherwise offx|y are 0 */
		C->offx = (int) (C->extend * C->Grid->header->n_columns);
		C->offy = (int) (C->extend * C->Grid->header->n_rows);
	}
	xtradim = (int)C->scale + 2;
	gxdim   = C->Grid->header->n_columns + 2 * C->offx;
	gydim   = C->Grid->header->n_rows + 2 * C->offy;
	bounds[0] = wbnd[0] = C->Grid->header->wesn[XLO] - C->offx * C->Grid->header->inc[GMT_X];
	bounds[1] = wbnd[1] = C->Grid->header->wesn[XHI] + C->offx * C->Grid->header->inc[GMT_X];
	bounds[2] = wbnd[2] = C->Grid->header->wesn[YLO] - C->offy * C->Grid->header->inc[GMT_Y];
	bounds[3] = wbnd[3] = C->Grid->header->wesn[YHI] + C->offy * C->Grid->header->inc[GMT_Y];
		
	/* extend the bounds slightly to be sure no data gets missed */
	xx = MIN(0.05*(bounds[1] - bounds[0]), 0.1);
	yy = MIN(0.05*(bounds[3] - bounds[2]), 0.1);
	bounds[0] -= xx;
	bounds[1] += xx;
	bounds[2] -= yy;
	bounds[3] += yy;

	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Processing input table data\n");
	/* allocate memory for arrays */
	C->grid_aux = gmt_M_memory (GMT, NULL, (size_t)(gxdim * gydim), float);
	C->num  = gmt_M_memory (GMT, NULL, (size_t)(gxdim * gydim), unsigned short int);
	norm = gmt_M_memory (GMT, NULL, (size_t)(gxdim * gydim), float);
	cnt  = gmt_M_memory (GMT, NULL, (size_t)(gxdim * gydim), unsigned short int);

	/* Read in xyz data and computes index no and store it in a structure */
	if ((error = GMT_Set_Columns (GMT->parent, GMT_IN, 3, GMT_COL_FIX_NO_TEXT)) != GMT_NOERROR)
		return (error);
	if (GMT_Init_IO (GMT->parent, GMT_IS_DATASET, GMT_IS_POINT, GMT_IN, GMT_ADD_DEFAULT, 0, options) != GMT_NOERROR)	/* Establishes data input */
		return (GMT->parent->error);
	
	if (GMT_Begin_IO (GMT->parent, GMT_IS_DATASET, GMT_IN, GMT_HEADER_ON) != GMT_NOERROR)	/* Enables data input and sets access mode */
		return (GMT->parent->error);

	do {	/* Keep returning records until we reach EOF */
		if ((In = GMT_Get_Record (GMT->parent, GMT_READ_DATA, NULL)) == NULL) {	/* Read next record, get NULL if special case */
			if (gmt_M_rec_is_error (GMT)) 		/* Bail if there are any read errors */
				return (GMT_RUNTIME_ERROR);
			else if (gmt_M_rec_is_eof (GMT)) 		/* Reached end of file */
				break;
			continue;	/* Go back and read the next record */
		}

		/* Data record to process */
		in = In->data;	/* Only need to process numerical part here */
	
		if (gmt_M_is_dnan (in[GMT_Z])) continue;	/* Cannot use NaN values */

		/* get position in grid */
		ix = (int)((in[GMT_X] - wbnd[0] + 0.5 * C->Grid->header->inc[GMT_X]) / C->Grid->header->inc[GMT_X]);
		iy = (int)((in[GMT_Y] - wbnd[2] + 0.5 * C->Grid->header->inc[GMT_Y]) / C->Grid->header->inc[GMT_Y]);

		/* process the data */
		if (ix >= -xtradim && ix < gxdim + xtradim && iy >= -xtradim && iy < gydim + xtradim) {
			ix1 = MAX (ix - xtradim, 0);
			ix2 = MIN (ix + xtradim, gxdim - 1);
			iy1 = MAX (iy - xtradim, 0);
			iy2 = MIN (iy + xtradim, gydim - 1);
			for (ii = ix1; ii <= ix2; ii++) {
				xx = wbnd[0] + ii * C->Grid->header->inc[GMT_X] - in[GMT_X];
				xx2 = xx * xx;
				for (jj = iy1; jj <= iy2; jj++) {
					kgrid = ii * gydim + jj;
					yy = wbnd[2] + jj * C->Grid->header->inc[GMT_Y] - in[GMT_Y];
					weight = exp(-(xx2 + yy*yy)*factor);
					norm[kgrid] += (float)weight;
					C->grid_aux[kgrid] += (float)(weight * in[GMT_Z]);
					C->num[kgrid]++;
					if (ii == ix && jj == iy)
						cnt[kgrid]++;
				}
			}
			k++;
		}
		else if (ix >= 0 && ix < gxdim && iy >= 0 && iy < gydim) {
			kgrid = ix*gydim + iy;
			if (C->num[kgrid] <= 0) {
				norm[kgrid] = 1.0;
				C->grid_aux[kgrid] = (float)in[GMT_Z];
				C->num[kgrid] = 1;
				cnt[kgrid] = 1;
			}
			k++;
		}

	} while (true);
	
	if (GMT_End_IO (GMT->parent, GMT_IN, 0) != GMT_OK)	/* Disables further data input */
		return (EXIT_FAILURE);

	C->npoints = k;

	if (C->npoints == 0) {
		GMT_Report (GMT->parent, GMT_MSG_NORMAL, "No datapoints inside region, aborts\n");
		return (EXIT_FAILURE);
	}

	if (C->breakfile) {		/* A breakline */
		struct GMT_DATASET *Lin = NULL;
		struct GMT_DATATABLE *xyzline = NULL;

		if ((Lin = GMT_Read_Data (GMT->parent, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_LINE, GMT_READ_NORMAL, NULL, C->breakfile, NULL)) == NULL)
			return (EXIT_FAILURE);
		xyzline = Lin->table[0];			/* Can only be one table since we read a single file */
		C->npoints += interp_breakline (GMT, C, norm, cnt, wbnd, xyzline, factor, xtradim, gxdim, gydim);
		free (C->breakfile);
		C->breakfile = NULL;	/* So this code chunk is not executed again if this function gets called twice */
	}

	C->xcoords = gmt_M_memory (GMT, NULL, (size_t)gxdim, float);
	C->ycoords = gmt_M_memory (GMT, NULL, (size_t)gydim, float);
			
	if (C->interp_method == INTERP_SURFACE) {
		C->data = gmt_M_memory (GMT, NULL, (size_t)C->npoints, struct SURFACE_DATA);

		/* now loop over all points in the output grid */
		for (k = i = 0; i < gxdim; i++) {
			C->xcoords[i] = (float)(C->Grid->header->wesn[XLO] + C->Grid->header->inc[GMT_X]*i);
			for (j = 0; j < gydim; j++) {
				C->ycoords[j] = (float)(C->Grid->header->wesn[YLO] + C->Grid->header->inc[GMT_Y]*j);
				kgrid = i * gydim + j;
				if (cnt[kgrid] > 0) {
					C->grid_aux[kgrid] = C->grid_aux[kgrid] / norm[kgrid];
					C->nbinset++;
					C->data[k].index = i * C->ny + j;
					C->data[k].x = C->xcoords[i];
					C->data[k].y = C->ycoords[j];
					C->data[k].z = C->grid_aux[kgrid];
					k++;
				}
				else
					C->grid_aux[kgrid] = (float)C->clipvalue;
			}
		}
		C->data = gmt_M_memory (GMT, C->data, (size_t)C->nbinset, struct SURFACE_DATA); /* freed in mb_surface */
	}
	else {
		C->index = gmt_M_memory (GMT, NULL, (size_t)C->npoints, uint64_t);

		for (k = i = 0; i < gxdim; i++) {
			C->xcoords[i] = (float)(C->Grid->header->wesn[XLO] + C->Grid->header->inc[GMT_X]*i);
			for (j = 0; j < gydim; j++) {
				C->ycoords[j] = (float)(C->Grid->header->wesn[YLO] + C->Grid->header->inc[GMT_Y]*j);
				kgrid = i * gydim + j;
				if (cnt[kgrid] > 0) {
					C->grid_aux[kgrid] = C->grid_aux[kgrid] / norm[kgrid];
					C->index[k] = i * C->ny + j;
					C->nbinset++;
					k++;
				}
				else
					C->grid_aux[kgrid] = (float)C->clipvalue;
			}
		}
		C->index = gmt_M_memory (GMT, C->index, (size_t)C->nbinset, uint64_t);
	}

	C->npoints = C->nbinset;    /* From now on only nbinset counts (so we should only have one var) */
	C->gxdim = gxdim;           /* Make copies to later use in main */
	C->gydim = gydim;

	gmt_M_free (GMT, norm);       /* Not needed anymore */
	gmt_M_free (GMT, cnt);

	return (0);
}

/*--------------------------------------------------------------------
 *    The MB-system:	mb_zgrid.c	    4/25/95
 *    $Id: mb_zgrid.c,v 5.0 2000/12/01 22:53:59 caress Exp $
 *
 *    Copyright (c) 1993, 1994, 1995, 2000 by
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
 * This is a function to generate thin plate spline interpolation
 * of a data field. This code originated as fortran in the
 * 1960's and was used routinely at the Institute of 
 * Geophysics and Planetary Physics at the Scripps Institution
 * of Oceanography through the 1970's and 1980's. The Fortran
 * code was obtained from Professory Robert L. Parker at
 * IGPP in 1989.
 * It was converted to C in 1995 for use with the MB-System
 * software package.
 * 
 * The nature of the interpolation is controlled by the
 * parameters cay and nrng: cay sets the tension of the
 * interpolation such that cay=0.0 yields a pure Laplace
 * (minimum curvature) solution and cay=infinity yields
 * a pure thin plate spline solution. A cay=1e10 value
 * has commonly been used to yield spline solutions.
 * The nrng value sets the number of grid spaces from
 * data that will be interpolated; if nrng exceeds the
 * maximum dimension of the grid then the entire grid
 * will be interpolated.
 * 
 * The input parameters are:
 *     nx,ny = max subscripts of z in x and y directions . 
 *     x1,y1 = coordinates of z(1,1) 
 *     dx,dy = x and y increments . 
 *     xyz(3,*) = array giving x-y position and hgt of each data point. 
 *     n = length of xyz series. 
 *     zpij[n] = float work array
 *     knxt[n] = int work array
 *     imnew[MAX(nx, ny)+1] = int work array
 *     cay = k = amount of spline eqn (between 0 and inf.) 
 *     nrng...grid points more than nrng grid spaces from the nearest 
 *            data point are set to undefined. 
 *
 * Author:	Unknown, but "jdt", "ian crain",  and "dr t murty"
 *              obviously contributed.
 * Hacker:	D. W. Caress
 * Date:	April 25, 1995
 *
 *     The following are the original comments from the Fortran code:
 *
 *     sets up square grid for contouring , given arbitrarily placed 
 *     data points. laplace interpolation is used. 
 *     the method used here was lifted directly from notes left by 
 *     mr ian crain formerly with the comp.science div. 
 *     info on relaxation soln of laplace eqn supplied by dr t murty. 
 *     fortran ii   oceanography/emr   dec/68   jdt 
 *
 *     z = 2-d array of hgts to be set up. points outside region to be 
 *     contoured should be initialized to 10**35 . the rest should be 0.0 
 *
 *     modification feb/69   to get smoother results a portion of the 
 *     beam eqn  was added to the laplace eqn giving 
 *     delta2x(z)+delta2y(z) - k(delta4x(z)+delta4y(z)) = 0 . 
 *     k=0 gives pure laplace solution.  k=infinity gives pure spline 
 *     solution. 
 *
 *     nx,ny = max subscripts of z in x and y directions . 
 *     x1,y1 = coordinates of z(1,1) 
 *     dx,dy = x and y increments . 
 *     xyz(3,*) = array giving x-y position and hgt of each data point. 
 *     n = length of xyz series. 
 *     cay = k = amount of spline eqn (between 0 and inf.) 
 *     nrng...grid points more than nrng grid spaces from the nearest 
 *            data point are set to undefined. 
 *
 *     modification dec23/69   data pts no longer moved to grid pts. 
 *
 *     modification feb/85  common blocks work1 and work2 replaced by 
 *     dimension statement and parameters nwork, mwork introduced. 
 *
 *     modification feb/90  nwork and mwork replaced by maxdat and maxdim 
 *     for compatibility with command driven interface program 
 *     David W. Caress
 *
 *     zgrid.c -- translated by f2c (version 19950314) from zgrid.f.
 *     Work arrays zpij[n], knxt[n], imnew[MAX(nx, ny)+1] are now
 *     passed into the function.
 *     David W. Caress
 *     April 25,  1995
 *--------------------------------------------------------------------*/

/*----------------------------------------------------------------------- */
int mb_zgrid(struct GMTMBGRID_INFO *C, int verbose, float *z, uint64_t n, float *zpij, int *knxt, int *imnew, int nrng) {
    /* System generated locals */
    int z_offset, i__2;

    /* Local variables */
	uint64_t i, j, k, iter, nnew, itmax, jmnew = 0;
	uint64_t kk, im, jm, npg, npt;
	int nx = C->gxdim, ny = C->gydim;
	float *xyz = C->grid_aux;
	double delz, zijn, zmin, zimm, zmax, zjmm, zipp, zjpp;
	double root, zsum, zpxy, a, b, c, d;
	double x, y, zbase, relax, delzm, derzm, dzmax, dzrms;
	double dzrms8 = 0, z00, dz, ze, hrange, zn, zs, zw, zrange, dzmaxf, relaxn, rootgs, dzrmsp, big, abz;
	double eps, zim = 0, zjm = 0, wgt, zip, zjp, tpy, zxy;

	double x1 = C->Grid->header->wesn[XLO];
	double y1 = C->Grid->header->wesn[YLO];
	double dx = C->Grid->header->inc[GMT_X];
	double dy = C->Grid->header->inc[GMT_Y];
	double cay = C->tension;

    /* Parameter adjustments */
    z_offset = nx + 1;
    z -= z_offset;

    /* Function Body */
    itmax = 100;
    eps = 0.002;
    big = 9e29;

/*     get zbase which will make all zp values positive by 20*(zmax-zmin) 
 * ********************************************************************** 
*/

	zmax = -C->clipvalue;
	zmin =  C->clipvalue;
	for (k = 1; k < n; k++) {
		if (xyz[C->index[k]] > zmax)
			zmax = xyz[C->index[k]];
		else if (xyz[C->index[k]] < zmin)
			zmin = xyz[C->index[k]];
	}
    zrange = zmax - zmin;
    zbase = (zrange * 20 - zmin);
    hrange = MIN(dx * (nx - 1), dy * (ny - 1));
    derzm = (zrange * 2 / hrange);

/*     set pointer array knxt */
/* ********************************************************************** 
*/

	for (kk = 1; kk <= n; kk++) {
		k = n + 1 - kk;
		knxt[k - 1] = 0;
		i = C->index[k-1] / C->gydim;
		i = (int)((C->xcoords[i] - x1) / dx + 1.5);
		if (i * (nx + 1 - i) <= 0) continue;

		j = C->index[k-1] % C->gydim;
		j = (int)((C->ycoords[j] - y1) / dy + 1.5);
		if (j * (ny + 1 - j) <= 0) continue;

		if (z[i + j * nx] >= big) continue;

		knxt[k - 1] = (int)n + 1;

		if (z[i + j * nx] > 0.)
			knxt[k - 1] = (int)(z[i + j * nx] + .5);

		z[i + j * nx] = (float) k;
    }

/*     affix each data point zp to its nearby grid point.  take avg zp if 
 *     more than one zp nearby the grid point. add zbase and complement. 
 * ********************************************************************** */

	for (k = 1; k <= n; k++) {
		if (knxt[k - 1] <= 0)
			continue;

		npt = 0;
		zsum = 0.;
		i = (int)((C->xcoords[C->index[k-1] / C->gydim] - x1) / dx + 1.5);
		j = (int)((C->ycoords[C->index[k-1] % C->gydim] - y1) / dy + 1.5);
		kk = k;
L70:
		npt++;
		zsum += xyz[C->index[k-1]];
		knxt[kk - 1] = -knxt[kk - 1];
		kk = -knxt[kk - 1];
		if (kk <= n) goto L70;
		z[i + j * nx] = (float)(-zsum / npt - zbase);
    }

/*     initially set each unset grid point to value of nearest known pt. 
 * ********************************************************************** 
*/

	for (i = 1; i <= nx; ++i) {
		for (j = 1; j <= ny; ++j) {
			if (z[i + j * nx] != 0.)
				continue;
			z[i + j * nx] = (float)-1e35;
		}
	}

	for (iter = 1; iter <= nrng; ++iter) {
		nnew = 0;
		for (i = 1; i <= nx; ++i) {
			for (j = 1; j <= ny; ++j) {
				if (z[i + j * nx] + big >= 0.)
					goto L192;

				if (j - 1 <= 0) goto L162;

				if (jmnew <= 0) goto L154;
				else
					goto L162;

L154:
				zijn = (float)fabs(z[i + (j - 1) * nx]);
				if (zijn >= big) goto L162;
				else
					goto L195;

L162:
				if (i <= 1) goto L172;

				if (imnew[j - 1] <= 0) goto L164;
				else
					goto L172;

L164:
				zijn = (float)fabs(z[i - 1 + j * nx]);
				if (zijn >= big) goto L172;
				else
					goto L195;

L172:
				if (j >= ny) goto L182;

				zijn = (float)fabs(z[i + (j + 1) * nx]);
				if (zijn >= big) goto L182;
				else
					goto L195;

L182:
				if (i >= nx) goto L192;

				zijn = fabs(z[i + 1 + j * nx]);
				if (zijn >= big) goto L192;
				else
					goto L195;

L192:
				imnew[j - 1] = 0;
				jmnew = 0;
				goto L197;
L195:
				imnew[j - 1] = 1;
				jmnew = 1;
				z[i + j * nx] = (float)zijn;
				++nnew;
L197:
		;
			}
		}
		if (nnew <= 0) goto L200;
	}

L200:
	for (i = 1; i <= nx; ++i) {
		for (j = 1; j <= ny; ++j) {
			abz = fabs(z[i + j * nx]);
			if (abz >= big)
				z[i + j * nx] = (float)abz;

		}
	}

/*     improve the non-data points by applying point over-relaxation */
/*     using the laplace-spline equation  (carres method is used) */
/* ********************************************************************** */

    dzrmsp = zrange;
    relax = 1.;
    for (iter = 1; iter <= itmax; ++iter) {
	dzrms = 0.;
	dzmax = 0.;
	npg = 0;
	for (i = 1; i <= nx; ++i) {
		for (j = 1; j <= ny; ++j) {
			z00 = z[i + j * nx];
			if (z00 >= big) continue;
			if (z00 < 0) continue;

			wgt = 0.;
			zsum = 0.;

			im = 0;
			if (i <= 1) goto L570;

			zim = fabs(z[i - 1 + j * nx]);
			if (zim >= big) goto L570;

			im = 1;
			wgt += 1.;
			zsum += zim;
			if (i <= 2) goto L570;

			zimm = fabs(z[i - 2 + j * nx]);
			if (zimm >= big) goto L570;

			wgt += cay;
			zsum -= cay * (zimm - zim * 2);
L570:
			if (nx <= i) goto L700;

			zip = fabs(z[i + 1 + j * nx]);
			if (zip >= big) goto L700;

			wgt += 1.;
			zsum += zip;
			if (im <= 0) goto L620;

			wgt += cay * 4;
			zsum += cay * 2 * (zim + zip);
L620:
			if (nx - 1 - i <= 0)
				goto L700;

			zipp = fabs(z[i + 2 + j * nx]);
			if (zipp >= big) goto L700;

			wgt += cay;
			zsum -= cay * (zipp - zip * 2);
L700:

			jm = 0;
			if (j <= 1) goto L1570;

			zjm = fabs(z[i + (j - 1) * nx]);
			if (zjm >= big) goto L1570;

			jm = 1;
			wgt += 1.;
			zsum += zjm;
			if (j <= 2) goto L1570;

			zjmm = fabs(z[i + (j - 2) * nx]);
			if (zjmm >= big) goto L1570;

			wgt += cay;
			zsum -= (cay * (zjmm - zjm * 2));
L1570:
			if (ny <= j) goto L1700;

			zjp = fabs(z[i + (j + 1) * nx]);
			if (zjp >= big) goto L1700;

			wgt += 1.;
			zsum += zjp;
			if (jm <= 0) goto L1620;

			wgt += (cay * 4);
			zsum += (cay * 2 * (zjm + zjp));
L1620:
			if (ny - 1 - j <= 0) goto L1700;

			zjpp = fabs(z[i + (j + 2) * nx]);
			if (zjpp >= big) goto L1700;

			wgt += cay;
			zsum -= (cay * (zjpp - zjp * 2));
L1700:

			dz = zsum / wgt - z00;
			npg++;
			dzrms += dz * dz;
			dzmax = MAX(fabs(dz), dzmax);
			z[i + j * nx] = (float)(z00 + dz * relax);
		}
	}


/*     shift data points zp progressively back to their proper places as */
/*     the shape of surface z becomes evident. */
/* ********************************************************************** */

	if (iter - iter / 10 * 10 != 0)
		goto L3600;

	for (k = 1; k <= n; ++k) {
		knxt[k - 1] = (i__2 = knxt[k - 1], abs(i__2));
		if (knxt[k - 1] <= 0)
			continue;

		x = C->xcoords[C->index[k-1] / C->gydim];
		y = C->ycoords[C->index[k-1] % C->gydim];

		i = (int)(x + 1.5);
		x = (x + 1. - i);
		j = (int)(y + 1.5);
		y = (y + 1. - j);
		zpxy = xyz[C->index[k-1]] + zbase;
		z00 = fabs(z[i + j * nx]);

		zw = 1e35;
		if (i <= 1) goto L3120;

		zw = fabs(z[i - 1 + j * nx]);
L3120:
		ze = 1e35;
		if (i >= nx) goto L3140;

		ze = fabs(z[i + 1 + j * nx]);
L3140:
		if (ze >= big)
			goto L3150;
		else
			goto L3160;

L3150:
		if (zw >= big)
			goto L3170;
		else
			goto L3180;

L3160:
		if (zw >= big)
			goto L3190;
		else
			goto L3200;

L3170:
		ze = z00;
		zw = z00;
		goto L3200;
L3180:
		ze = (z00 * 2. - zw);
		goto L3200;
L3190:
		zw = (z00 * 2. - ze);

L3200:
		zs = 1e35;
		if (j <= 1) goto L3220;

		zs = fabs(z[i + (j - 1) * nx]);
L3220:
		zn = 1e35;
		if (j >= ny) goto L3240;

		zn = fabs(z[i + (j + 1) * nx]);
L3240:
		if (zn >= big)
			goto L3250;
		else
			goto L3260;

L3250:
		if (zs >= big) 
			goto L3270;
		else
			goto L3280;

L3260:
		if (zs >= big)
			goto L3290;
		else
			goto L3300;

L3270:
		zn = z00;
		zs = z00;
		goto L3300;
L3280:
		zn = (z00 * 2. - zs);
		goto L3300;
L3290:
		zs = (z00 * 2. - zn);

L3300:
		a = ((ze - zw) * .5);
		b = ((zn - zs) * .5);
		c = ((ze + zw) * .5 - z00);
		d = ((zn + zs) * .5 - z00);
		zxy = z00 + a * x + b * y + c * x * x + d * y * y;
		delz = z00 - zxy;
		delzm = (derzm * (fabs(x) * dx + fabs(y) * dy) * 0.8);
		if (delz > delzm)
			delz = delzm;
		if (delz + delzm >= 0.)
			goto L3365;

		delz = -delzm;
L3365:
		zpij[k - 1] = (float)(zpxy + delz);
	}

	for (k = 1; k <= n; ++k) {
		if (knxt[k - 1] > 0) {
			npt = 0;
			zsum = 0.;
			i = (int)((C->xcoords[C->index[k-1] / C->gydim] - x1) / dx + 1.5);
			j = (int)((C->ycoords[C->index[k-1] % C->gydim] - y1) / dy + 1.5);
			kk = k;
L3420:
			npt++;
			zsum += zpij[kk - 1];
			knxt[kk - 1] = -knxt[kk - 1];
			kk = -knxt[kk - 1];
			if (kk <= n) goto L3420;
			z[i + j * nx] = (float)(-zsum / npt);
	    }
	}
L3600:

/*     test for convergence */
/* ****************************************************************** */
/* all grid points assigned */
	if (npg <= 1) goto L4010;

	dzrms = sqrt(dzrms / npg);
	root = dzrms / dzrmsp;
	dzrmsp = dzrms;
	dzmaxf = dzmax / zrange;
	if (iter - iter / 10 * 10 - 2 != 0)
		goto L3715;

	dzrms8 = dzrms;
L3715:
	if (iter - iter / 10 * 10 != 0)
	    goto L4000;

	root = sqrt(sqrt(sqrt(dzrms / dzrms8)));
	if (root >= 0.9999) goto L4000;

	if (dzmaxf / (1. - root) <= eps)
		goto L4010;

/*     improve the relaxation factor. */
/* ******************************************************************** */

	if ((iter - 20) * (iter - 40) * (iter - 60) != 0)
	    goto L4000;

	if (relax - 1. - root >= 0.)
	    goto L4000;

	tpy = ((root + relax - 1.) / relax);
	rootgs = tpy * tpy / root;
	relaxn = (2. / (sqrt(1. - rootgs) + 1.));
	if (iter != 60)
	    goto L3780;
	else
	    goto L3785;

L3780:
	relaxn -= ((2. - relaxn) * 0.25);
L3785:
	relax = MAX(relax,relaxn);
L4000:
	;
	if (verbose)
		if (iter % 10 == 0)
			fprintf(stderr,"Iteration %lld of a maximum of %lld\r", iter, itmax); 
    }
L4010:
	if (verbose) fprintf(stderr,"\n");

/*     remove zbase from array z and return. */
/* ********************************************************************** 
*/

	for (i = 1; i <= nx; ++i) {
		for (j = 1; j <= ny; ++j) {
			if (z[i + j * nx] < big)
				z[i + j * nx] = (float)(fabs(z[i + j * nx]) - zbase);
		}
	}
    return 0;
}


/*--------------------------------------------------------------------
 *    The MB-system:	mb_surface.c	5/2/94
 *    $Id: mb_surface.c,v 5.2 2005/03/25 04:09:53 caress Exp $
 *
 *    Inclusion in MB-System:
 *    Copyright (c) 1994, 2003 by
 *    David W. Caress (caress@mbari.org)
 *      Monterey Bay Aquarium Research Institute
 *      Moss Landing, CA 95039
 *    and Dale N. Chayes (dale@ldeo.columbia.edu)
 *      Lamont-Doherty Earth Observatory
 *      Palisades, NY 10964
 *
 *    Algorithm and original code:
 *    Copyright (c) 1991 by P. Wessel and W. H. F. Smith
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * SURFUNC is a function for gridding data using a minimum curvature
 * algorithm developed by W.H.F. Smith and P. Wessel.  The source
 * code below is almost entirely taken from the source code
 * "surface.c" distributed as part of the GMT-system by Wessel
 * and Smith. The MB-System Copyright notice above applies only
 * to the inclusion of this code in MB-System and changes made to
 * the code as part of that inclusion.  The algorithm and the
 * bulk of the code remains Copyright (c) by P. Wessel and W. H. F. Smith.
 *
 *--------------------------------------------------------------------*
 *
 * The original Smith and Wessel comments follow:
 *
 * surface.c:  a gridding program.
 * reads xyz triples and fits a surface to the data.
 * surface satisfies (1 - T) D4 z - T D2 z = 0,
 * where D4 is the biharmonic operator,
 * D2 is the Laplacian,
 * and T is a "tension factor" between 0 and 1.
 * End member T = 0 is the classical minimum curvature
 * surface.  T = 1 gives a harmonic surface.  Use T = 0.25
 * or so for potential data; something more for topography.
 *
 * Program includes overrelaxation for fast convergence and
 * automatic optimal grid factorization.
 *
 * See Smith & Wessel (Geophysics, 3, 293-305, 1990) for details.
 *
 * Authors: Walter H. F. Smith and Paul Wessel
 * Date: April, 1988.
 *
 *--------------------------------------------------------------------*
 *
 * Particulars regarding turning the program "surface" version 4.3
 * (revision of 26 February, 1992) into a function "surfunc" follow:
 *
 * Author:	D. W. Caress
 * Date:	May 2, 1994
 *
 * $Log: mb_surface.c,v $
 * Revision 5.2  2005/03/25 04:09:53  caress
 * Problems with global variables in mb_surface.c stomping on similarly named global variables in some programs has been fixed.
 *
 * Revision 5.1  2003/04/29 20:27:48  caress
 * Fixed multiple definitions of "error".
 *
 * Revision 5.0  2003/03/22 03:10:36  caress
 * Reinserting this code into MB-System for first time in years.
 *
 * Revision 4.2  1994/10/21  13:02:31  caress
 * Release V4.0
 *
 * Revision 4.1  1994/06/04  02:02:01  caress
 * Fixed several bugs and made some stylistic changes to
 * the output.  Changed the data input bounds to be much
 * larger than the working grid bounds.
 *
 * Revision 4.0  1994/05/05  20:30:06  caress
 * First cut. This code derived from GMT program surface by
 * Walter Smith and Paul Wessel.
 *
 */

int mb_surface(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, uint64_t ndat, float *sgrid) {

	/* int	size_query = FALSE; */
	int	 serror = FALSE;
	
	if (C->Grid->header->wesn[XLO] >= C->Grid->header->wesn[XHI] || 
		C->Grid->header->wesn[YLO] >= C->Grid->header->wesn[YHI]) serror = TRUE;
	if (C->Grid->header->inc[GMT_X] <= 0.0 || C->Grid->header->inc[GMT_Y] <= 0.0) serror = TRUE;

	if (C->tension != 0.0) {
		C->boundary_tension = C->tension;
		C->interior_tension = C->tension;
	}
	C->relax_old = 1.0 - C->relax_new;

	/* New stuff here for v4.3:  Check out the grid dimensions:  */
	C->grid = gcd_euclid(C->nx-1, C->ny-1);

	/*
	if (C->grid == 1) fprintf(stderr,"surface:  WARNING:  Your C->grid dimensions are mutually prime.\n");
	if (C->grid == 1 || size_query) suggest_sizes_for_surface(C->nx-1, C->ny-1);
	if (size_query) exit(0);
	*/

	/* New idea: set C->grid = 1, read data, setting index.  Then throw
		away data that can't be used in end game, constraining
		size of briggs->b[6] structure.  */
	
	C->grid = 1;
	set_grid_parameters(C);
	//read_data(GMT, C, ndat, xdat, ydat, zdat);
	throw_away_unusables(GMT, C);
	remove_planar_trend(C);
	rescale_z_values(GMT, C);
	load_constraints(GMT, C);
	
	/* Set up factors and reset C->grid to first value  */
	
	C->grid = gcd_euclid(C->nx-1, C->ny-1);
	C->n_fact = get_prime_factors(C->grid, C->factors);
	set_grid_parameters(C);
	while ( C->block_nx < 4 || C->block_ny < 4 ) {
		smart_divide(C);
		set_grid_parameters(C);
	}
	set_offset(C);
	set_index(C);
	/* Now the data are ready to go for the first iteration.  */

	/* Allocate more space  */
	
	C->briggs = gmt_M_memory (GMT, NULL, (size_t)C->npoints, struct SURFACE_BRIGGS);
	C->iu = gmt_M_memory (GMT, NULL, (size_t)(C->mx * C->my), char);
	C->u = gmt_M_memory (GMT, NULL, (size_t)(C->mx * C->my), float);

	if (C->radius > 0) initialize_grid(C); /* Fill in nodes with a weighted avg in a search radius  */

	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Grid\tMode\tIteration\tMax Change\tConv Limit\tTotal Iterations\n");
	
	set_coefficients(C);
	
	C->old_grid = C->grid;
	find_nearest_point (C);
	iterate (GMT, C, 1);
	 
	while (C->grid > 1) {
		smart_divide (C);
		set_grid_parameters(C);
		set_offset(C);
		set_index (C);
		fill_in_forecast (C);
		iterate(GMT, C, 0);
		C->old_grid = C->grid;
		find_nearest_point (C);
		iterate (GMT, C, 1);
	}
	
	if (gmt_M_is_verbose (GMT, GMT_MSG_VERBOSE)) check_errors (GMT, C);

	replace_planar_trend(C);
	get_output(C, sgrid);		/* put the result in 'sgrid' */

	gmt_M_free (GMT, C->data);
	gmt_M_free (GMT, C->briggs);
	gmt_M_free (GMT, C->iu);
	gmt_M_free (GMT, C->u);
	if (C->set_low) gmt_M_free (GMT, C->lower);
	if (C->set_high) gmt_M_free (GMT, C->upper);

	return(0);
}

void set_coefficients(struct GMTMBGRID_INFO *C) {
	double	e_4, loose, a0;
	
	loose = 1.0 - C->interior_tension;
	C->e_2 = 1.0;
	e_4 = C->e_2 * C->e_2;
	C->eps_p2 = C->e_2;
	C->eps_m2 = 1.0/C->e_2;
	C->one_plus_e2 = 1.0 + C->e_2;
	C->two_plus_ep2 = 2.0 + 2.0*C->eps_p2;
	C->two_plus_em2 = 2.0 + 2.0*C->eps_m2;

	C->x_edge_const = 4 * C->one_plus_e2 - 2 * (C->interior_tension / loose);
	C->e_m2 = 1.0 / C->e_2;
	C->y_edge_const = 4 * (1.0 + C->e_m2) - 2 * (C->interior_tension * C->e_m2 / loose);

	a0 = 1.0 / ( (6 * e_4 * loose + 10 * C->e_2 * loose + 8 * loose - 2 * C->one_plus_e2) + 
		4 * C->interior_tension * C->one_plus_e2);
	C->a0_const_1 = 2 * loose * (1.0 + e_4);
	C->a0_const_2 = 2.0 - C->interior_tension + 2 * loose * C->e_2;

	C->coeff[1][4] = C->coeff[1][7] = -loose;
	C->coeff[1][0] = C->coeff[1][11] = -loose * e_4;
	C->coeff[0][4] = C->coeff[0][7] = -loose * a0;
	C->coeff[0][0] = C->coeff[0][11] = -loose * e_4 * a0;
	C->coeff[1][5] = C->coeff[1][6] = 2 * loose * C->one_plus_e2;
	C->coeff[0][5] = C->coeff[0][6] = (2 * C->coeff[1][5] + C->interior_tension) * a0;
	C->coeff[1][2] = C->coeff[1][9] = C->coeff[1][5] * C->e_2;
	C->coeff[0][2] = C->coeff[0][9] = C->coeff[0][5] * C->e_2;
	C->coeff[1][1] = C->coeff[1][3] = C->coeff[1][8] = C->coeff[1][10] = -2 * loose * C->e_2;
	C->coeff[0][1] = C->coeff[0][3] = C->coeff[0][8] = C->coeff[0][10] = C->coeff[1][1] * a0;
	
	C->e_2 *= 2;		/* We will need these in boundary conditions  */
	C->e_m2 *= 2;
	
	C->ij_sw_corner = 2 * C->my + 2;			/*  Corners of array of actual data  */
	C->ij_se_corner = C->ij_sw_corner + (C->nx - 1) * C->my;
	C->ij_nw_corner = C->ij_sw_corner + (C->ny - 1);
	C->ij_ne_corner = C->ij_se_corner + (C->ny - 1);

}

void set_offset(struct GMTMBGRID_INFO *C) {
	uint64_t	add_w[5], add_e[5], add_s[5], add_n[5], add_w2[5], add_e2[5], add_s2[5], add_n2[5];
	int	i, j, kase;
	
	add_w[0] = -C->my; add_w[1] = add_w[2] = add_w[3] = add_w[4] = -C->grid_east;
	add_w2[0] = -2 * C->my;  add_w2[1] = -C->my - C->grid_east;  add_w2[2] = add_w2[3] = add_w2[4] = -2 * C->grid_east;
	add_e[4] = C->my; add_e[0] = add_e[1] = add_e[2] = add_e[3] = C->grid_east;
	add_e2[4] = 2 * C->my;  add_e2[3] = C->my + C->grid_east;  add_e2[2] = add_e2[1] = add_e2[0] = 2 * C->grid_east;

	add_n[4] = 1; add_n[3] = add_n[2] = add_n[1] = add_n[0] = C->grid;
	add_n2[4] = 2;  add_n2[3] = C->grid + 1;  add_n2[2] = add_n2[1] = add_n2[0] = 2 * C->grid;
	add_s[0] = -1; add_s[1] = add_s[2] = add_s[3] = add_s[4] = -C->grid;
	add_s2[0] = -2;  add_s2[1] = -C->grid - 1;  add_s2[2] = add_s2[3] = add_s2[4] = -2 * C->grid;

	for (i = 0, kase = 0; i < 5; i++) {
		for (j = 0; j < 5; j++, kase++) {
			C->offset[kase][0] = add_n2[j];
			C->offset[kase][1] = add_n[j] + add_w[i];
			C->offset[kase][2] = add_n[j];
			C->offset[kase][3] = add_n[j] + add_e[i];
			C->offset[kase][4] = add_w2[i];
			C->offset[kase][5] = add_w[i];
			C->offset[kase][6] = add_e[i];
			C->offset[kase][7] = add_e2[i];
			C->offset[kase][8] = add_s[j] + add_w[i];
			C->offset[kase][9] = add_s[j];
			C->offset[kase][10] = add_s[j] + add_e[i];
			C->offset[kase][11] = add_s2[j];
		}
	}
}

void fill_in_forecast (struct GMTMBGRID_INFO *C) {

	/* Fills in bilinear estimates into new node locations after grid is divided.   */

	uint64_t i, j, ii, jj, index_0, index_1, index_2, index_3;
	uint64_t index_new;
	double delta_x, delta_y, a0, a1, a2, a3;
	double old_size;
	
	old_size = 1.0 / (double)C->old_grid;

	/* first do from southwest corner */
	
	for (i = 0; i < C->nx-1; i += C->old_grid) {
		for (j = 0; j < C->ny-1; j += C->old_grid) {
			
			/* get indices of bilinear square */
			index_0 = C->ij_sw_corner + i * C->my + j;
			index_1 = index_0 + C->old_grid * C->my;
			index_2 = index_1 + C->old_grid;
			index_3 = index_0 + C->old_grid;
			
			/* get coefficients */
			a0 = C->u[index_0];
			a1 = C->u[index_1] - a0;
			a2 = C->u[index_3] - a0;
			a3 = C->u[index_2] - a0 - a1 - a2;
			
			/* find all possible new fill ins */
			
			for (ii = i;  ii < i + C->old_grid; ii += C->grid) {
				delta_x = (ii - i) * old_size;
				for (jj = j;  jj < j + C->old_grid; jj += C->grid) {
					index_new = C->ij_sw_corner + ii * C->my + jj;
					if (index_new == index_0) continue;
					delta_y = (jj - j) * old_size;
					C->u[index_new] = (float)(a0 + a1 * delta_x + delta_y * ( a2 + a3 * delta_x));
					C->iu[index_new] = 0;
				}
			}
			C->iu[index_0] = 5;
		}
	}
	
	/* now do linear guess along east edge */
	
	for (j = 0; j < (C->ny-1); j += C->old_grid) {
		index_0 = C->ij_se_corner + j;
		index_3 = index_0 + C->old_grid;
		for (jj = j;  jj < j + C->old_grid; jj += C->grid) {
			index_new = C->ij_se_corner + jj;
			delta_y = (jj - j) * old_size;
			C->u[index_new] = C->u[index_0] + (float)(delta_y * (C->u[index_3] - C->u[index_0]));
			C->iu[index_new] = 0;
		}
		C->iu[index_0] = 5;
	}
	/* now do linear guess along north edge */
	for (i = 0; i < (C->nx-1); i += C->old_grid) {
		index_0 = C->ij_nw_corner + i * C->my;
		index_1 = index_0 + C->old_grid * C->my;
		for (ii = i;  ii < i + C->old_grid; ii += C->grid) {
			index_new = C->ij_nw_corner + ii * C->my;
			delta_x = (ii - i) * old_size;
			C->u[index_new] = C->u[index_0] + (float)(delta_x * (C->u[index_1] - C->u[index_0]));
			C->iu[index_new] = 0;
		}
		C->iu[index_0] = 5;
	}
	/* now set northeast corner to fixed and we're done */
	C->iu[C->ij_ne_corner] = 5;
}

int compare_points (const void *point_1v, const void *point_2v) {
		/*  Routine for qsort to sort data structure for fast access to data by node location.
		    Sorts on index first, then on radius to node corresponding to index, so that index
		    goes from low to high, and so does radius.
		*/
	uint64_t block_i, block_j, index_1, index_2;
	double x0, y0, dist_1, dist_2;
	const struct SURFACE_DATA *point_1 = point_1v, *point_2 = point_2v;
	
	index_1 = point_1->index;
	index_2 = point_2->index;
	if (index_1 < index_2) return (-1);
	if (index_1 > index_2) return (1);
	if (index_1 == SURFACE_OUTSIDE) return (0);
	/* Points are in same grid cell, find the one who is nearest to grid point */
	block_i = point_1->index/GMT_Surface_Global.block_ny;
	block_j = point_1->index%GMT_Surface_Global.block_ny;
	x0 = GMT_Surface_Global.x_min + block_i * GMT_Surface_Global.grid_xinc;
	y0 = GMT_Surface_Global.y_min + block_j * GMT_Surface_Global.grid_yinc;
	dist_1 = (point_1->x - x0) * (point_1->x - x0) + (point_1->y - y0) * (point_1->y - y0);
	dist_2 = (point_2->x - x0) * (point_2->x - x0) + (point_2->y - y0) * (point_2->y - y0);
	if (dist_1 < dist_2) return (-1);
	if (dist_1 > dist_2) return (1);
	return (0);
}

void smart_divide (struct GMTMBGRID_INFO *C) {		/* Divide grid by its largest prime factor */
	C->grid /= C->factors[C->n_fact - 1];
	C->n_fact--;
}

void set_index (struct GMTMBGRID_INFO *C) {
		/* recomputes data[k].index for new value of grid,
		   sorts data on index and radii, and throws away
		   data which are now outside the useable limits. */
	int i, j, k, k_skipped = 0;

	for (k = 0; k < C->npoints; k++) {
		i = (int)floor(((C->data[k].x- C->Grid->header->wesn[XLO])*C->r_grid_xinc) + 0.5);
		j = (int)floor(((C->data[k].y- C->Grid->header->wesn[YLO])*C->r_grid_yinc) + 0.5);
		if (i < 0 || i >= C->block_nx || j < 0 || j >= C->block_ny) {
			C->data[k].index = OUTSIDE;
			k_skipped++;
		}
		else
			C->data[k].index = i * C->block_ny + j;
	}
	
	qsort ((char *)C->data, C->npoints, sizeof (struct SURFACE_DATA), compare_points);
	
	C->npoints -= k_skipped;
	
}

void find_nearest_point(struct GMTMBGRID_INFO *C) {
	int i, j, k, last_index, block_i, block_j, briggs_index;
	uint64_t iu_index;
	double x0, y0, dx, dy, xys, xy1, btemp;
	double b0, b1, b2, b3, b4, b5;
	
	last_index = -1;
	C->small = 0.05 * ((C->grid_xinc < C->grid_yinc) ? C->grid_xinc : C->grid_yinc);

	for (i = 0; i < C->nx; i += C->grid)	/* Reset grid info */
		for (j = 0; j < C->ny; j += C->grid)
			C->iu[C->ij_sw_corner + i*C->my + j] = 0;
	
	briggs_index = 0;
	for (k = 0; k < C->npoints; k++) {	/* Find constraining value  */
		if (C->data[k].index != last_index) {
			block_i = C->data[k].index/C->block_ny;
			block_j = C->data[k].index % C->block_ny;
			last_index = C->data[k].index;
	 		iu_index = C->ij_sw_corner + (block_i * C->my + block_j) * C->grid;
	 		x0 = C->Grid->header->wesn[XLO] + block_i*C->grid_xinc;
	 		y0 = C->Grid->header->wesn[YLO] + block_j*C->grid_yinc;
	 		dx = (C->data[k].x - x0)*C->r_grid_xinc;
	 		dy = (C->data[k].y - y0)*C->r_grid_yinc;
	 		if (fabs(dx) < C->small && fabs(dy) < C->small) {
	 			C->iu[iu_index] = 5;
	 			C->u[iu_index] = C->data[k].z;
	 		}
	 		else {
	 			if (dx >= 0.0) {
	 				if (dy >= 0.0)
	 					C->iu[iu_index] = 1;
	 				else
	 					C->iu[iu_index] = 4;
	 			}
	 			else {
	 				if (dy >= 0.0)
	 					C->iu[iu_index] = 2;
	 				else
	 					C->iu[iu_index] = 3;
	 			}
	 			dx = fabs(dx);
	 			dy = fabs(dy);
	 			btemp = 2 * C->one_plus_e2 / ( (dx + dy) * (1.0 + dx + dy) );
	 			b0 = 1.0 - 0.5 * (dx + (dx * dx)) * btemp;
	 			b3 = 0.5 * (C->e_2 - (dy + (dy * dy)) * btemp);
	 			xys = 1.0 + dx + dy;
	 			xy1 = 1.0 / xys;
	 			b1 = (C->e_2 * xys - 4 * dy) * xy1;
	 			b2 = 2 * (dy - dx + 1.0) * xy1;
	 			b4 = b0 + b1 + b2 + b3 + btemp;
	 			b5 = btemp * C->data[k].z;
	 			C->briggs[briggs_index].b[0] = b0;
	 			C->briggs[briggs_index].b[1] = b1;
	 			C->briggs[briggs_index].b[2] = b2;
	 			C->briggs[briggs_index].b[3] = b3;
	 			C->briggs[briggs_index].b[4] = b4;
	 			C->briggs[briggs_index].b[5] = b5;
	 			briggs_index++;
	 		}
	 	}
	 }
}

void set_grid_parameters(struct GMTMBGRID_INFO *C) {			
	GMT_Surface_Global.block_ny = C->block_ny = (C->ny - 1) / C->grid + 1;
	C->block_nx = (C->nx - 1) / C->grid + 1;
	GMT_Surface_Global.grid_xinc = C->grid_xinc = C->grid * C->Grid->header->inc[GMT_X];
	GMT_Surface_Global.grid_yinc = C->grid_yinc = C->grid * C->Grid->header->inc[GMT_Y];
	C->grid_east = C->grid * C->my;
	C->r_grid_xinc = 1.0 / C->grid_xinc;
	C->r_grid_yinc = 1.0 / C->grid_yinc;
}

void initialize_grid(struct GMTMBGRID_INFO *C) {
	/*
	 * For the initial gridsize, compute weighted averages of data inside the search radius
	 * and assign the values to u[i,j] where i,j are multiples of gridsize.
	 */
	 int	irad, jrad, i, j, imin, imax, jmin, jmax, index_1, index_2, k, ki, kj, k_index;
	 double	r, rfact, sum_w, sum_zw, weight, x0, y0;

	 irad = (int)ceil(C->radius/C->grid_xinc);
	 jrad = (int)ceil(C->radius/C->grid_yinc);
	 rfact = -4.5/(C->radius*C->radius);
	 
	 for (i = 0; i < C->block_nx; i ++ ) {
	 	x0 = C->Grid->header->wesn[XLO] + i*C->grid_xinc;
	 	for (j = 0; j < C->block_ny; j ++ ) {
	 		y0 = C->Grid->header->wesn[YLO] + j*C->grid_yinc;
	 		imin = i - irad;
	 		if (imin < 0) imin = 0;
	 		imax = i + irad;
	 		if (imax >= C->block_nx) imax = C->block_nx - 1;
	 		jmin = j - jrad;
	 		if (jmin < 0) jmin = 0;
	 		jmax = j + jrad;
	 		if (jmax >= C->block_ny) jmax = C->block_ny - 1;
	 		index_1 = imin * C->block_ny + jmin;
	 		index_2 = imax * C->block_ny + jmax + 1;
	 		sum_w = sum_zw = 0.0;
	 		k = 0;
	 		while (k < C->npoints && C->data[k].index < index_1) k++;
	 		for (ki = imin; k < C->npoints && ki <= imax && C->data[k].index < index_2; ki++) {
	 			for (kj = jmin; k < C->npoints && kj <= jmax && C->data[k].index < index_2; kj++) {
	 				k_index = ki*C->block_ny + kj;
	 				while (k < C->npoints && C->data[k].index < k_index) k++;
	 				while (k < C->npoints && C->data[k].index == k_index) {
	 					r = (C->data[k].x-x0)*(C->data[k].x-x0) + (C->data[k].y-y0)*(C->data[k].y-y0);
	 					weight = exp (rfact*r);
	 					sum_w += weight;
	 					sum_zw += weight*C->data[k].z;
	 					k++;
	 				}
	 			}
	 		}
	 		if (sum_w == 0.0) {
	 			fprintf (stderr, "surface: Warning: no data inside search radius at: %.8g %.8g\n", x0, y0);
	 			C->u[C->ij_sw_corner + (i * C->my + j) * C->grid] = (float)C->z_mean;
	 		}
	 		else {
	 			C->u[C->ij_sw_corner + (i*C->my+j)*C->grid] = (float)(sum_zw/sum_w);
	 		}
		}
	}
}


void new_initialize_grid(struct GMTMBGRID_INFO *C) {
	/*
	 * For the initial gridsize, load constrained nodes with weighted avg of their data;
	 * and then do something with the unconstrained ones.
	 */
	 int	k, k_index, block_i, block_j;
	 uint64_t u_index;
	 double	sum_w, sum_zw, weight, x0, y0, dx, dy, dx_scale, dy_scale;

	dx_scale = 4.0 / C->grid_xinc;
	dy_scale = 4.0 / C->grid_yinc;
	C->n_empty = C->block_ny * C->block_nx;
	k = 0;
	while (k < C->npoints) {
		block_i = C->data[k].index / C->block_ny;
		block_j = C->data[k].index % C->block_ny;
		x0 = C->Grid->header->wesn[XLO] + block_i*C->grid_xinc;
		y0 = C->Grid->header->wesn[YLO] + block_j*C->grid_yinc;
		u_index = C->ij_sw_corner + (block_i*C->my + block_j) * C->grid;
		k_index = C->data[k].index;
		
		dy = (C->data[k].y - y0) * dy_scale;
		dx = (C->data[k].x - x0) * dx_scale;
		sum_w = 1.0 / (1.0 + dx*dx + dy*dy);
		sum_zw = C->data[k].z * sum_w;
		k++;

		while (k < C->npoints && C->data[k].index == k_index) {
			dy = (C->data[k].y - y0) * dy_scale;
			dx = (C->data[k].x - x0) * dx_scale;
			weight = 1.0 / (1.0 + dx*dx + dy*dy);
			sum_zw += C->data[k].z * weight;
			sum_w += weight;
			sum_zw += weight*C->data[k].z;
			k++;
	 	}
	 	C->u[u_index] = (float)(sum_zw/sum_w);
	 	C->iu[u_index] = 5;
	 	C->n_empty--;
	 }
}

/* This function rewritten by D.W. Caress 5/3/94 */
int read_data(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, int ndat, float *xdat, float *ydat, float *zdat) {

	int	i, j, k, kmax = 0, kmin = 0, idat;
	double	zmin = 1.0e38, zmax = -1.0e38;

	C->data = gmt_M_memory (GMT, NULL, (size_t)ndat, struct SURFACE_DATA);

	/* Read in xyz data and computes index no and store it in a structure */
	k = 0;
	C->z_mean = 0;
	for (idat = 0; idat < ndat; idat++) {		/* Here, C->block_nx|y is still equal to C->nx|y */
		i = (int)floor(((xdat[idat] - C->Grid->header->wesn[XLO])*C->r_grid_xinc) + 0.5);
		j = (int)floor(((ydat[idat] - C->Grid->header->wesn[YLO])*C->r_grid_yinc) + 0.5);
		if (i >= 0 && i < C->block_nx && j >= 0 && j < C->block_ny) {
			C->data[k].index = i * C->block_ny + j;
			C->data[k].x = C->xcoords[i];
			C->data[k].y = C->ycoords[j];
			C->data[k].z = zdat[idat];
			if (zmin > zdat[idat]) {
				zmin = zdat[idat];
				kmin = k;
			}
			if (zmax < zdat[idat]) {
				zmax = zdat[idat];
				kmax = k;
			}
			k++;
			C->z_mean += zdat[idat];
		}
	}

	C->npoints = k;
	C->z_mean /= k;
	if (C->converge_limit == 0.0)
		C->converge_limit = 0.001 * C->z_scale; /* c_l = 1 ppt of L2 scale */

	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "surface: Minimum value of your dataset x,y,z at: %g %g %g\n",
			C->data[kmin].x, C->data[kmin].y, C->data[kmin].z);
	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "surface: Maximum value of your dataset x,y,z at: %g %g %g\n",
			C->data[kmax].x, C->data[kmax].y, C->data[kmax].z);

	if (C->set_low == 1)
		C->low_limit = C->data[kmin].z;
	else if (C->set_low == 2 && C->low_limit > C->data[kmin].z) {
	/*	C->low_limit = data[kmin].z;	*/
		/*
		fprintf (stderr, "surface: Warning:  Your lower value is > than min data value.\n");
		*/
	}
	if (C->set_high == 1)
		C->high_limit = C->data[kmax].z;
	else if (C->set_high == 2 && C->high_limit < C->data[kmax].z) {
	/*	C->high_limit = data[kmax].z;	*/
		/*
		fprintf (stderr, "surface: Warning:  Your upper value is < than max data value.\n");
		*/
	}
	return(0);
}

/* this function rewritten from write_output() by D.W. Caress 5/3/94 */
void get_output(struct GMTMBGRID_INFO *C, float *sgrid) {
	uint64_t index, i, j;

	index = C->ij_sw_corner;
	for(i = 0; i < C->nx; i++, index += C->my) 
		for (j = 0; j < C->ny; j++) 
			sgrid[j*C->nx+i] = C->u[index + C->ny - j - 1];
}
	
int	iterate(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, int mode) {

	uint64_t i, j, k, ij, kase, briggs_index, ij_v2;
	int	x_case, y_case, x_w_case, x_e_case, y_s_case, y_n_case;
	int	iteration_count = 0;
	float *u = C->u;

	double	current_limit = C->converge_limit / C->grid;
	double	change, max_change = 0.0, busum, sum_ij;
	double	b0, b1, b2, b3, b4, b5;
	
	double	x_0_const = 4.0 * (1.0 - C->boundary_tension) / (2.0 - C->boundary_tension);
	double	x_1_const = (3 * C->boundary_tension - 2.0) / (2.0 - C->boundary_tension);
	double	y_denom = 2 * (1.0 - C->boundary_tension) + C->boundary_tension;
	double	y_0_const = 4 * (1.0 - C->boundary_tension) / y_denom;
	double	y_1_const = (C->boundary_tension - 2 * (1.0 - C->boundary_tension) ) / y_denom;

	do {
		briggs_index = 0;	/* Reset the constraint table stack pointer  */
		
		max_change = -1.0;
		
		/* Fill in auxiliary boundary values (in new way) */
		
		/* First set d2[]/dn2 = 0 along edges:  */
		/* New experiment : (1-T)d2[]/dn2 + Td[]/dn = 0  */
		
		
		
		for (i = 0; i < C->nx; i += C->grid) {
			/* set d2[]/dy2 = 0 on south side:  */
			ij = C->ij_sw_corner + i * C->my;
			/* u[ij - 1] = 2 * u[ij] - u[ij + C->grid];  */
			u[ij - 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij + C->grid]);
			/* set d2[]/dy2 = 0 on north side:  */
			ij = C->ij_nw_corner + i * C->my;
			/* u[ij + 1] = 2 * u[ij] - u[ij - C->grid];  */
			u[ij + 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij - C->grid]);
			
		}
		
		for (j = 0; j < C->ny; j += C->grid) {
			/* set d2[]/dx2 = 0 on west side:  */
			ij = C->ij_sw_corner + j;
			/* u[ij - C->my] = 2 * u[ij] - u[ij + grid_east];  */
			u[ij - C->my] = (float)(x_1_const * u[ij + C->grid_east] + x_0_const * u[ij]);
			/* set d2[]/dx2 = 0 on east side:  */
			ij = C->ij_se_corner + j;
			/* u[ij + C->my] = 2 * u[ij] - u[ij - grid_east];  */
			u[ij + C->my] = (float)(x_1_const * u[ij - C->grid_east] + x_0_const * u[ij]);
		}
			
		/* Now set d2[]/dxdy = 0 at each corner:  */
		
		ij = C->ij_sw_corner;
		u[ij - C->my - 1] = u[ij + C->grid_east - 1] + u[ij - C->my + C->grid] - u[ij + C->grid_east + C->grid];
				
		ij = C->ij_nw_corner;
		u[ij - C->my + 1] = u[ij + C->grid_east + 1] + u[ij - C->my - C->grid] - u[ij + C->grid_east - C->grid];
				
		ij = C->ij_se_corner;
		u[ij + C->my - 1] = u[ij - C->grid_east - 1] + u[ij + C->my + C->grid] - u[ij - C->grid_east + C->grid];
				
		ij = C->ij_ne_corner;
		u[ij + C->my + 1] = u[ij - C->grid_east + 1] + u[ij + C->my - C->grid] - u[ij - C->grid_east - C->grid];
		
		/* Now set (1-T)dC/dn + Tdu/dn = 0 at each edge :  */
		/* New experiment:  only dC/dn = 0  */
		
		x_w_case = 0;
		x_e_case = C->block_nx - 1;
		for (i = 0; i < C->nx; i += C->grid, x_w_case++, x_e_case--) {
		
			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;
				
			/* South side :  */
			kase = x_case * 5;
			ij = C->ij_sw_corner + i * C->my;
			u[ij + C->offset[kase][11]] = 
				(float)(u[ij + C->offset[kase][0]] + C->eps_m2*(u[ij + C->offset[kase][1]] + u[ij + C->offset[kase][3]]
					- u[ij + C->offset[kase][8]] - u[ij + C->offset[kase][10]])
					+ C->two_plus_em2 * (u[ij + C->offset[kase][9]] - u[ij + C->offset[kase][2]]) );
				/*  + tense * eps_m2 * (u[ij + offset[kase][2]] - u[ij + offset[kase][9]]) / (1.0 - tense);  */
			/* North side :  */
			kase = x_case * 5 + 4;
			ij = C->ij_nw_corner + i * C->my;
			u[ij + C->offset[kase][0]] = 
				-(float)(-u[ij + C->offset[kase][11]] + C->eps_m2 * (u[ij + C->offset[kase][1]] + u[ij + C->offset[kase][3]]
					- u[ij + C->offset[kase][8]] - u[ij + C->offset[kase][10]])
					+ C->two_plus_em2 * (u[ij + C->offset[kase][9]] - u[ij + C->offset[kase][2]]) );
				/*  - tense * eps_m2 * (u[ij + offset[kase][2]] - u[ij + offset[kase][9]]) / (1.0 - tense);  */
		}
		
		y_s_case = 0;
		y_n_case = C->block_ny - 1;
		for (j = 0; j < C->ny; j += C->grid, y_s_case++, y_n_case--) {
				
			if(y_s_case < 2)
				y_case = y_s_case;
			else if(y_n_case < 2)
				y_case = 4 - y_n_case;
			else
				y_case = 2;
			
			/* West side :  */
			kase = y_case;
			ij = C->ij_sw_corner + j;
			u[ij+C->offset[kase][4]] = 
				u[ij + C->offset[kase][7]] + (float)(C->eps_p2 * (u[ij + C->offset[kase][3]] + u[ij + C->offset[kase][10]]
				-u[ij + C->offset[kase][1]] - u[ij + C->offset[kase][8]])
				+ C->two_plus_ep2 * (u[ij + C->offset[kase][5]] - u[ij + C->offset[kase][6]]));
				/*  + tense * (u[ij + offset[kase][6]] - u[ij + offset[kase][5]]) / (1.0 - tense);  */
			/* East side :  */
			kase = 20 + y_case;
			ij = C->ij_se_corner + j;
			u[ij + C->offset[kase][7]] = 
				- (float)(-u[ij + C->offset[kase][4]] + C->eps_p2 * (u[ij + C->offset[kase][3]] + u[ij + C->offset[kase][10]]
				- u[ij + C->offset[kase][1]] - u[ij + C->offset[kase][8]])
				+ C->two_plus_ep2 * (u[ij + C->offset[kase][5]] - u[ij + C->offset[kase][6]]) );
				/*  - tense * (u[ij + offset[kase][6]] - u[ij + offset[kase][5]]) / (1.0 - tense);  */
		}

		/* That's it for the boundary points.  Now loop over all data  */
		
		x_w_case = 0;
		x_e_case = C->block_nx - 1;
		for (i = 0; i < C->nx; i += C->grid, x_w_case++, x_e_case--) {
		
			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;
			
			y_s_case = 0;
			y_n_case = C->block_ny - 1;
			
			ij = C->ij_sw_corner + i * C->my;
			
			for (j = 0; j < C->ny; j += C->grid, ij += C->grid, y_s_case++, y_n_case--) {
	
				if (C->iu[ij] == 5) continue;	/* Point is fixed  */
				
				if(y_s_case < 2)
					y_case = y_s_case;
				else if(y_n_case < 2)
					y_case = 4 - y_n_case;
				else
					y_case = 2;
				
				kase = x_case * 5 + y_case;
				sum_ij = 0.0;

				if (C->iu[ij] == 0) {		/* Point is unconstrained  */
					for (k = 0; k < 12; k++)
						sum_ij += (u[ij + C->offset[kase][k]] * C->coeff[0][k]);
				}
				else {				/* Point is constrained  */
					b0 = C->briggs[briggs_index].b[0];
					b1 = C->briggs[briggs_index].b[1];
					b2 = C->briggs[briggs_index].b[2];
					b3 = C->briggs[briggs_index].b[3];
					b4 = C->briggs[briggs_index].b[4];
					b5 = C->briggs[briggs_index].b[5];
					briggs_index++;
					if (C->iu[ij] < 3) {
						if (C->iu[ij] == 1) {	/* Point is in quadrant 1  */
							busum = b0 * u[ij + C->offset[kase][10]]
								+ b1 * u[ij + C->offset[kase][9]]
								+ b2 * u[ij + C->offset[kase][5]]
								+ b3 * u[ij + C->offset[kase][1]];
						}
						else {			/* Point is in quadrant 2  */
							busum = b0 * u[ij + C->offset[kase][8]]
								+ b1 * u[ij + C->offset[kase][9]]
								+ b2 * u[ij + C->offset[kase][6]]
								+ b3 * u[ij + C->offset[kase][3]];
						}
					}
					else {
						if (C->iu[ij] == 3) {	/* Point is in quadrant 3  */
							busum = b0 * u[ij + C->offset[kase][1]]
								+ b1 * u[ij + C->offset[kase][2]]
								+ b2 * u[ij + C->offset[kase][6]]
								+ b3 * u[ij + C->offset[kase][10]];
						}
						else {		/* Point is in quadrant 4  */
							busum = b0 * u[ij + C->offset[kase][3]]
								+ b1 * u[ij + C->offset[kase][2]]
								+ b2 * u[ij + C->offset[kase][5]]
								+ b3 * u[ij + C->offset[kase][8]];
						}
					}
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + C->offset[kase][k]] * C->coeff[1][k]);
					}
					sum_ij = (sum_ij + C->a0_const_2 * (busum + b5))
						/ (C->a0_const_1 + C->a0_const_2 * b4);
				}
				
				/* New relaxation here  */
				sum_ij = u[ij] * C->relax_old + sum_ij * C->relax_new;
				
				if (C->constrained) {	/* Must check limits.  Note lower/upper is v2 format and need ij_v2! */
					ij_v2 = (C->ny - j - 1) * C->nx + i;
					if (C->set_low /*&& !GMT_is_fnan((double)C->lower[ij_v2])*/ && sum_ij < C->lower[ij_v2])
						sum_ij = C->lower[ij_v2];
					else if (C->set_high /*&& !GMT_is_fnan((double)upper[ij_v2])*/ && sum_ij > C->upper[ij_v2])
						sum_ij = C->upper[ij_v2];
				}
					
				change = fabs(sum_ij - u[ij]);
				u[ij] = (float)sum_ij;
				if (change > max_change) max_change = change;
			}
		}
		iteration_count++;
		C->total_iterations++;
		max_change *= C->z_scale;	/* Put max_change into z units  */
		GMT_Report (GMT->parent, GMT_MSG_LONG_VERBOSE, "%4d\t%c\t%8d\t%10g\t%10g\t%10d\n", C->grid, 
				C->mode_type[mode], iteration_count, max_change, current_limit, C->total_iterations);

	} while (max_change > current_limit && iteration_count < C->max_iterations);
	
	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "%4d\t%c\t%8d\t%10g\t%10g\t%10d\n",
		C->grid, C->mode_type[mode], iteration_count, max_change, current_limit, C->total_iterations);

	return(iteration_count);
}

void check_errors (struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C) {

	uint64_t i, j, k, ij, n_nodes, move_over[12];	/* move_over = offset[kase][12], but C->grid = 1 so move_over is easy  */
	
	double	x0, y0, dx, dy, mean_error, mean_squared_error, z_est, z_err, curvature, c;
	double	du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2, d3u_dx3, d3u_dx2dy, d3u_dxdy2, d3u_dy3;
	
	double	x_0_const = 4.0 * (1.0 - C->boundary_tension) / (2.0 - C->boundary_tension);
	double	x_1_const = (3 * C->boundary_tension - 2.0) / (2.0 - C->boundary_tension);
	double	y_denom = 2  * (1.0 - C->boundary_tension) + C->boundary_tension;
	double	y_0_const = 4 * (1.0 - C->boundary_tension) / y_denom;
	double	y_1_const = (C->boundary_tension - 2 * (1.0 - C->boundary_tension) ) / y_denom;
	
	
	move_over[0] = 2;
	move_over[1] = 1 - C->my;
	move_over[2] = 1;
	move_over[3] = 1 + C->my;
	move_over[4] = -2 * C->my;
	move_over[5] = -C->my;
	move_over[6] = C->my;
	move_over[7] = 2 * C->my;
	move_over[8] = -1 - C->my;
	move_over[9] = -1;
	move_over[10] = -1 + C->my;
	move_over[11] = -2;

	mean_error = 0;
	mean_squared_error = 0;
	
	/* First update the boundary values  */

	for (i = 0; i < C->nx; i ++) {
		ij = C->ij_sw_corner + i * C->my;
		C->u[ij - 1] = (float)(y_0_const * C->u[ij] + y_1_const * C->u[ij + 1]);
		ij = C->ij_nw_corner + i * C->my;
		C->u[ij + 1] = (float)(y_0_const * C->u[ij] + y_1_const * C->u[ij - 1]);
	}

	for (j = 0; j < C->ny; j ++) {
		ij = C->ij_sw_corner + j;
		C->u[ij - C->my] = (float)(x_1_const * C->u[ij + C->my] + x_0_const * C->u[ij]);
		ij = C->ij_se_corner + j;
		C->u[ij + C->my] = (float)(x_1_const * C->u[ij - C->my] + x_0_const * C->u[ij]);
	}

	ij = C->ij_sw_corner;
	C->u[ij - C->my - 1] = C->u[ij + C->my - 1] + C->u[ij - C->my + 1] - C->u[ij + C->my + 1];
	ij = C->ij_nw_corner;
	C->u[ij - C->my + 1] = C->u[ij + C->my + 1] + C->u[ij - C->my - 1] - C->u[ij + C->my - 1];
	ij = C->ij_se_corner;
	C->u[ij + C->my - 1] = C->u[ij - C->my - 1] + C->u[ij + C->my + 1] - C->u[ij - C->my + 1];
	ij = C->ij_ne_corner;
	C->u[ij + C->my + 1] = C->u[ij - C->my + 1] + C->u[ij + C->my - 1] - C->u[ij - C->my - 1];

	for (i = 0; i < C->nx; i ++) {
				
		ij = C->ij_sw_corner + i * C->my;
		C->u[ij + move_over[11]] = 
			(float)(C->u[ij + move_over[0]] + C->eps_m2*(C->u[ij + move_over[1]] + C->u[ij + move_over[3]]
				- C->u[ij + move_over[8]] - C->u[ij + move_over[10]])
				+ C->two_plus_em2 * (C->u[ij + move_over[9]] - C->u[ij + move_over[2]]) );
					
		ij = C->ij_nw_corner + i * C->my;
		C->u[ij + move_over[0]] = 
			-(float)(-C->u[ij + move_over[11]] + C->eps_m2 * (C->u[ij + move_over[1]] + C->u[ij + move_over[3]]
				- C->u[ij + move_over[8]] - C->u[ij + move_over[10]])
				+ C->two_plus_em2 * (C->u[ij + move_over[9]] - C->u[ij + move_over[2]]) );
	}
		
	for (j = 0; j < C->ny; j ++) {
			
		ij = C->ij_sw_corner + j;
		C->u[ij+move_over[4]] = 
			C->u[ij + move_over[7]] + (float)(C->eps_p2 * (C->u[ij + move_over[3]] + C->u[ij + move_over[10]]
			-C->u[ij + move_over[1]] - C->u[ij + move_over[8]])
			+ C->two_plus_ep2 * (C->u[ij + move_over[5]] - C->u[ij + move_over[6]]));
				
		ij = C->ij_se_corner + j;
		C->u[ij + move_over[7]] = 
			- (float)(-C->u[ij + move_over[4]] + C->eps_p2 * (C->u[ij + move_over[3]] + C->u[ij + move_over[10]]
			- C->u[ij + move_over[1]] - C->u[ij + move_over[8]])
			+ C->two_plus_ep2 * (C->u[ij + move_over[5]] - C->u[ij + move_over[6]]) );
	}

	/* That resets the boundary values.  Now we can test all data.  
		Note that this loop checks all values, even though only nearest were used.  */
	
	for (k = 0; k < C->npoints; k++) {
		i = C->data[k].index/C->ny;
		j = C->data[k].index%C->ny;
	 	ij = C->ij_sw_corner + i * C->my + j;
	 	if ( C->iu[ij] == 5 ) continue;
	 	x0 = C->Grid->header->wesn[XLO] + i*C->Grid->header->inc[GMT_X];
	 	y0 = C->Grid->header->wesn[YLO] + j*C->Grid->header->inc[GMT_Y];
	 	dx = (C->data[k].x - x0)*C->r_grid_xinc;
	 	dy = (C->data[k].y - y0)*C->r_grid_yinc;
 
	 	du_dx = 0.5 * (C->u[ij + move_over[6]] - C->u[ij + move_over[5]]);
	 	du_dy = 0.5 * (C->u[ij + move_over[2]] - C->u[ij + move_over[9]]);
	 	d2u_dx2 = C->u[ij + move_over[6]] + C->u[ij + move_over[5]] - 2 * C->u[ij];
	 	d2u_dy2 = C->u[ij + move_over[2]] + C->u[ij + move_over[9]] - 2 * C->u[ij];
	 	d2u_dxdy = 0.25 * (C->u[ij + move_over[3]] - C->u[ij + move_over[1]]
	 			- C->u[ij + move_over[10]] + C->u[ij + move_over[8]]);
	 	d3u_dx3 = 0.5 * ( C->u[ij + move_over[7]] - 2 * C->u[ij + move_over[6]]
	 				+ 2 * C->u[ij + move_over[5]] - C->u[ij + move_over[4]]);
	 	d3u_dy3 = 0.5 * ( C->u[ij + move_over[0]] - 2 * C->u[ij + move_over[2]]
	 				+ 2 * C->u[ij + move_over[9]] - C->u[ij + move_over[11]]);
	 	d3u_dx2dy = 0.5 * ( ( C->u[ij + move_over[3]] + C->u[ij + move_over[1]] - 2 * C->u[ij + move_over[2]] )
	 				- ( C->u[ij + move_over[10]] + C->u[ij + move_over[8]] - 2 * C->u[ij + move_over[9]] ) );
	 	d3u_dxdy2 = 0.5 * ( ( C->u[ij + move_over[3]] + C->u[ij + move_over[10]] - 2 * C->u[ij + move_over[6]] )
	 				- ( C->u[ij + move_over[1]] + C->u[ij + move_over[8]] - 2 * C->u[ij + move_over[5]] ) );

	 	/* 3rd order Taylor approx:  */
	 		
	 	z_est = C->u[ij] + dx * (du_dx +  dx * ( (0.5 * d2u_dx2) + dx * (d3u_dx3 / 6.0) ) )
				+ dy * (du_dy +  dy * ( (0.5 * d2u_dy2) + dy * (d3u_dy3 / 6.0) ) )
	 			+ dx * dy * (d2u_dxdy) + (0.5 * dx * d3u_dx2dy) + (0.5 * dy * d3u_dxdy2);
	 		
	 	z_err = z_est - C->data[k].z;
	 	mean_error += z_err;
	 	mean_squared_error += (z_err * z_err);
	 }
	 mean_error /= C->npoints;
	 mean_squared_error = sqrt( mean_squared_error / C->npoints);
	 
	 curvature = 0.0;
	 n_nodes = C->nx * C->ny;
	 
	 for (i = 0; i < C->nx; i++) {
	 	for (j = 0; j < C->ny; j++) {
	 		ij = C->ij_sw_corner + i * C->my + j;
	 		c = C->u[ij + move_over[6]] + C->u[ij + move_over[5]]
	 			+ C->u[ij + move_over[2]] + C->u[ij + move_over[9]] - 4.0 * C->u[ij + move_over[6]];
			curvature += (c * c);
		}
	}

	 /*
	 fprintf(stderr, "Fit info: N data points  N nodes\tmean error\trms error\tcurvature\n");
	 fprintf (stderr,"\t%8d\t%8d\t%.8lg\t%.8lg\t%.8lg\n", npoints, n_nodes, mean_error, mean_squared_error,
	 	curvature);
	*/
	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "\nSpline interpolation fit information:\n");
	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "Data points   nodes    mean error     rms error     curvature\n");
	GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "%9d %9d   %10g   %10g  %10g\n",
			C->npoints, n_nodes, mean_error, mean_squared_error, curvature);
 }

void remove_planar_trend(struct GMTMBGRID_INFO *C) {

	int	i;
	double	a, b, c, d, xx, yy, zz;
	double	sx, sy, sz, sxx, sxy, sxz, syy, syz;
	struct GMT_GRID_HEADER *h = C->Grid->header;
	
	sx = sy = sz = sxx = sxy = sxz = syy = syz = 0.0;
	
	for (i = 0; i < C->npoints; i++) {
		xx = (C->data[i].x - h->wesn[XLO]) / h->inc[GMT_X];
		yy = (C->data[i].y - h->wesn[YLO]) / h->inc[GMT_Y];
		zz = C->data[i].z;

		sx += xx;
		sy += yy;
		sz += zz;
		sxx +=(xx * xx);
		sxy +=(xx * yy);
		sxz +=(xx * zz);
		syy +=(yy * yy);
		syz +=(yy * zz);
	}
	
	d = C->npoints*sxx*syy + 2*sx*sy*sxy - C->npoints*sxy*sxy - sx*sx*syy - sy*sy*sxx;
	
	if (d == 0.0) {
		C->plane_c0 = C->plane_c1 = C->plane_c2 = 0.0;
		return;
	}
	
	a = sz*sxx*syy + sx*sxy*syz + sy*sxy*sxz - sz*sxy*sxy - sx*sxz*syy - sy*syz*sxx;
	b = C->npoints*sxz*syy + sz*sy*sxy + sy*sx*syz - C->npoints*sxy*syz - sz*sx*syy - sy*sy*sxz;
	c = C->npoints*sxx*syz + sx*sy*sxz + sz*sx*sxy - C->npoints*sxy*sxz - sx*sx*syz - sz*sy*sxx;

	C->plane_c0 = a / d;
	C->plane_c1 = b / d;
	C->plane_c2 = c / d;

	for (i = 0; i < C->npoints; i++) {

		xx = (C->data[i].x - C->Grid->header->wesn[XLO]) * C->r_grid_xinc;
		yy = (C->data[i].y - C->Grid->header->wesn[YLO]) * C->r_grid_yinc;
		
		C->data[i].z -= (float)(C->plane_c0 + C->plane_c1 * xx + C->plane_c2 * yy);
	}

}

void replace_planar_trend(struct GMTMBGRID_INFO *C) {
	uint64_t i, j, ij;

	 for (i = 0; i < C->nx; i++) {
	 	for (j = 0; j < C->ny; j++) {
	 		ij = C->ij_sw_corner + i * C->my + j;
	 		C->u[ij] = (float)((C->u[ij] * C->z_scale) + (C->plane_c0 + C->plane_c1 * i + C->plane_c2 * j));
		}
	}
}

void throw_away_unusables(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C) {
	/* This is a new routine to eliminate data which will become
		unusable on the final iteration, when C->grid = 1.
		It assumes C->grid = 1 and set_grid_parameters has been
		called.  We sort, mark redundant data as OUTSIDE, and
		sort again, chopping off the excess.
		
		Experimental modification 5 Dec 1988 by Smith, as part
		of a new implementation using core memory for b[6]
		coefficients, eliminating calls to temp file.
	*/
	
	int	last_index, n_outside, k;
	
	/* Sort the data  */
	
	qsort ((char *)C->data, C->npoints, sizeof (struct SURFACE_DATA), compare_points);
	
	/* If more than one datum is indexed to same node, only the first should be kept.
		Mark the additional ones as OUTSIDE
	*/
	last_index = -1;
	n_outside = 0;
	for (k = 0; k < C->npoints; k++) {
		if (C->data[k].index == last_index) {
			C->data[k].index = OUTSIDE;
			n_outside++;
		}
		else {
			last_index = C->data[k].index;
		}
	}
	/* Sort again; this time the OUTSIDE points will be thrown away  */
	
	qsort ((char *)C->data, C->npoints, sizeof (struct SURFACE_DATA), compare_points);
	C->npoints -= n_outside;
	C->data = gmt_M_memory (GMT, C->data, (size_t)C->npoints, struct SURFACE_DATA);
	if (n_outside) {
		GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "surface: %d unusable points were supplied; these will be ignored.\n", n_outside);
		GMT_Report (GMT->parent, GMT_MSG_VERBOSE, "\tYou should have pre-processed the data with blockmean or blockmedian.\n");
	}

}

int rescale_z_values(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C) {
	int	i;
	double	ssz = 0.0;

	for (i = 0; i < C->npoints; i++)
		ssz += (C->data[i].z * C->data[i].z);
	
	/* Set z_scale = rms(z):  */
	
	C->z_scale = sqrt(ssz / C->npoints);

	if (C->z_scale < GMT_CONV8_LIMIT) {
		GMT_Report (GMT->parent, GMT_MSG_NORMAL, "Warning: Input data lie exactly on a plane.\n");
		C->r_z_scale = C->z_scale = 1.0;
		return (1);	/* Flag to tell the main to just write out the plane */
	}
	else
		C->r_z_scale = 1.0 / C->z_scale;

	/* Set z_scale = rms(z) */
	for (i = 0; i < C->npoints; i++)
		C->data[i].z *= (float)C->r_z_scale;

	if (C->converge_limit == 0.0)
		C->converge_limit = 0.001 * C->z_scale; /* i.e., 1 ppt of L2 scale */

	return(0);
}

void load_constraints (struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C) {
	
	/* Load C->lower/upper limits, verify range, deplane, and rescale */
}

double	guess_surface_time(struct GMTMBGRID_INFO *C) {
	/* Routine to guess a number proportional to the operations
	 * required by surface working on a user-desired C->grid of
	 * size nx by ny, where nx = (xmax - xmin)/dx, and same for
	 * ny.  (That is, one less than actually used in routine.)
	 *
	 * This is based on the following untested conjecture:
	 * 	The operations are proportional to T = nxg*nyg*L,
	 *	where L is a measure of the distance that data
	 *	constraints must propagate, and nxg, nyg are the
	 * 	current size of the grid.
	 *	For nx,ny relatively prime, we will go through only
	 * 	one grid cycle, L = max(nx,ny), and T = nx*ny*L.
	 *	But for nx,ny whose greatest common divisor is a highly
	 * 	composite number, we will have L equal to the division
	 * 	step made at each new C->grid cycle, and nxg,nyg will
	 * 	also be smaller than nx,ny.  Thus we can hope to find
	 *	some nx,ny for which the total value of T is small.
	 *
	 * The above is pure speculation and has not been derived
	 * empirically.  In actual practice, the distribution of the
	 * data, both spatially and in terms of their values, will
	 * have a strong effect on convergence.
	 *
	 * W. H. F. Smith, 26 Feb 1992.  */

	int	gcd;		/* Current value of the gcd  */
	int	nxg, nyg;	/* Current value of the C->grid dimensions  */
	int	nfactors;	/* Number of prime factors of current gcd  */
	int	factor;		/* Currently used factor  */
	/* Doubles are used below, even though the values will be integers,
		because the multiplications might reach sizes of O(n**3)  */
	double	t_sum;		/* Sum of values of T at each C->grid cycle  */
	double	length;		/* Current propagation distance.  */


	gcd = gcd_euclid(C->nx, C->ny);
	if (gcd > 1) {
		nfactors = get_prime_factors(gcd, C->factors);
		nxg = C->nx/gcd;
		nyg = C->ny/gcd;
		if (nxg < 3 || nyg < 3) {
			factor = C->factors[nfactors - 1];
			nfactors--;
			gcd /= factor;
			nxg *= factor;
			nyg *= factor;
		}
	}
	else {
		nxg = C->nx;
		nyg = C->ny;
	}
	length = MAX(nxg, nyg);
	t_sum = nxg * (nyg * length);	/* Make it double at each multiply  */

	/* Are there more C->grid cycles ?  */
	while (gcd > 1) {
		factor = C->factors[nfactors - 1];
		nfactors--;
		gcd /= factor;
		nxg *= factor;
		nyg *= factor;
		length = factor;
		t_sum += nxg * (nyg * length);
	}
	return(t_sum);
}


int	get_prime_factors(int n, int f[]) {
	/* Fills the integer array f with the prime factors of n.
	 * Returns the number of locations filled in f, which is
	 * one if n is prime.
	 *
	 * f[] should have been malloc'ed to enough space before
	 * calling prime_factors().  We can be certain that f[32]
	 * is enough space, for if n fits in a long, then n < 2**32,
	 * and so it must have fewer than 32 prime factors.  I think
	 * that in general, ceil(log2((double)n)) is enough storage
	 * space for f[].
	 *
	 * Tries 2,3,5 explicitly; then alternately adds 2 or 4
	 * to the previously tried factor to obtain the next trial
	 * factor.  This is done with the variable two_four_toggle.
	 * With this method we try 7,11,13,17,19,23,25,29,31,35,...
	 * up to a maximum of sqrt(n).  This shortened list results
	 * in 1/3 fewer divisions than if we simply tried all integers
	 * between 5 and sqrt(n).  We can reduce the size of the list
	 * of trials by an additional 20% by removing the multiples
	 * of 5, which are equal to 30m +/- 5, where m >= 1.  Starting
	 * from 25, these are found by alternately adding 10 or 20.
	 * To do this, we use the variable ten_twenty_toggle.
	 *
	 * W. H. F. Smith, 26 Feb 1992, after D.E. Knuth, vol. II  */

	int	current_factor;	/* The factor currently being tried  */
	int	max_factor;	/* Don't try any factors bigger than this  */
	int	n_factors = 0;	/* Returned; one if n is prime  */
	int	two_four_toggle = 0;	/* Used to add 2 or 4 to get next trial factor  */
	int	ten_twenty_toggle = 0;	/* Used to add 10 or 20 to skip_five  */
	int	skip_five = 25;	/* Used to skip multiples of 5 in the list  */
	int	m;	/* Used to keep a working copy of n  */


	/* Initialize m and max_factor  */
	m = abs(n);
	if (m < 2) return(0);
	max_factor = (int)floor(sqrt((double)m));

	/* First find the 2s  */
	current_factor = 2;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 3s  */
	current_factor = 3;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 5s  */
	current_factor = 5;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Now try all the rest  */

	while (m > 1 && current_factor <= max_factor) {

		/* Current factor is either 2 or 4 more than previous value  */

		if (two_four_toggle) {
			current_factor += 4;
			two_four_toggle = 0;
		}
		else {
			current_factor += 2;
			two_four_toggle = 1;
		}

		/* If current factor is a multiple of 5, skip it.  But first,
			set next value of skip_five according to 10/20 toggle:  */

		if (current_factor == skip_five) {
			if (ten_twenty_toggle) {
				skip_five += 20;
				ten_twenty_toggle = 0;
			}
			else {
				skip_five += 10;
				ten_twenty_toggle = 1;
			}
			continue;
		}

		/* Get here when current_factor is not a multiple of 2,3 or 5:  */

		while(!(m%current_factor)) {
			m /= current_factor;
			f[n_factors] = current_factor;
			n_factors++;
		}
	}

	/* Get here when all factors up to floor(sqrt(n)) have been tried.  */

	if (m > 1) {
		/* m is an additional prime factor of n  */
		f[n_factors] = m;
		n_factors++;
	}
	return (n_factors);
}
/* gcd_euclid.c  Greatest common divisor routine  */

#define IABS(i)	(((i) < 0) ? -(i) : (i))

int	gcd_euclid(int a, int b) {
	/* Returns the greatest common divisor of u and v by Euclid's method.
	 * I have experimented also with Stein's method, which involves only
	 * subtraction and left/right shifting; Euclid is faster, both for
	 * integers of size 0 - 1024 and also for random integers of a size
	 * which fits in a long integer.  Stein's algorithm might be better
	 * when the integers are HUGE, but for our purposes, Euclid is fine.
	 *
	 * Walter H. F. Smith, 25 Feb 1992, after D. E. Knuth, vol. II  */

	int	u,v,r;

	u = MAX(IABS(a), IABS(b));
	v = MIN(IABS(a), IABS(b));

	while (v > 0) {
		r = u%v;	/* Knuth notes that u < 2v 40% of the time;  */
		u = v;		/* thus we could have tried a subtraction  */
		v = r;		/* followed by an if test to do r = u%v  */
	}
	return(u);
}

int interp_breakline(struct GMT_CTRL *GMT, struct GMTMBGRID_INFO *C, float *norm,
                     unsigned short int *cnt, double *wbnd, struct GMT_DATATABLE *xyzline,
                     double factor, int xtradim, int gxdim, int gydim) {

	int n_tot = 0, this_ini = 0, this_end = 0, n_int = 0, i, j, k = 0, n;
	int kgrid, ix, iy, ix1, ix2, iy1, iy2, ii, jj, ndata = 0;
	size_t n_alloc;
	double *x, *y, *z = NULL, dx, dy, dz, r_dx, r_dy, weight, xx, xx2, yy;
	struct GMT_GRID_HEADER *h = NULL;

	n_alloc = GMT_CHUNK;
	x = gmt_M_memory (GMT, NULL, (size_t)n_alloc, double);
	y = gmt_M_memory (GMT, NULL, (size_t)n_alloc, double);
	z = gmt_M_memory (GMT, NULL, (size_t)n_alloc, double);

	r_dx = 1. / h->inc[GMT_X]; 
	r_dy = 1. / h->inc[GMT_Y]; 
	for (i = 0; i < xyzline->n_segments; i++) {
		for (j = 0; j < xyzline->segment[i]->n_rows - 1; j++) {
			dx = xyzline->segment[i]->data[GMT_X][j+1] - xyzline->segment[i]->data[GMT_X][j];
			dy = xyzline->segment[i]->data[GMT_Y][j+1] - xyzline->segment[i]->data[GMT_Y][j];
			dz = xyzline->segment[i]->data[GMT_Z][j+1] - xyzline->segment[i]->data[GMT_Z][j];
			n_int = (int)(lrint( MAX( fabs(dx) * r_dx, fabs(dy) * r_dy ) ) + 1);
			this_end += n_int;

			if (n_alloc >= this_end) {
				n_alloc += MAX (GMT_CHUNK, n_int);
				x = gmt_M_memory (GMT, x, n_alloc, double);
				y = gmt_M_memory (GMT, y, n_alloc, double);
				z = gmt_M_memory (GMT, z, n_alloc, double);
			}

			dx /= (floor(n_int) - 1);
			dy /= (floor(n_int) - 1);
			dz /= (floor(n_int) - 1);
			for (k = this_ini, n = 0; k < this_end - 1; k++, n++) {
				x[k] = xyzline->segment[i]->data[GMT_X][j] + n * dx;
				y[k] = xyzline->segment[i]->data[GMT_Y][j] + n * dy;
				z[k] = xyzline->segment[i]->data[GMT_Z][j] + n * dz;
			}
			x[this_end - 1] = xyzline->segment[i]->data[GMT_X][j+1];
			y[this_end - 1] = xyzline->segment[i]->data[GMT_Y][j+1];
			z[this_end - 1] = xyzline->segment[i]->data[GMT_Z][j+1];

			this_ini += n_int;
		}

		n_tot += this_end;
	}

	/* Now add the interpolated breakline to the C structure */

	for (n = 0; n < n_tot; n++) {

		if (gmt_M_is_dnan (z[n])) continue;

		/* get position in C->grid */
		ix = (int)((x[n] - wbnd[0] + 0.5*h->inc[GMT_X])/h->inc[GMT_X]);
		iy = (int)((y[n] - wbnd[2] + 0.5*h->inc[GMT_Y])/h->inc[GMT_Y]);

		/* process the data */
		if (ix >= -xtradim && ix < gxdim + xtradim && iy >= -xtradim && iy < gydim + xtradim) {
			ix1 = MAX(ix - xtradim, 0);
			ix2 = MIN(ix + xtradim, gxdim - 1);
			iy1 = MAX(iy - xtradim, 0);
			iy2 = MIN(iy + xtradim, gydim - 1);
			for (ii = ix1; ii <= ix2; ii++) {
				xx = wbnd[0] + ii*h->inc[GMT_X] - x[n];
				xx2 = xx * xx;
				for (jj = iy1; jj <= iy2; jj++) {
					kgrid = ii * gydim + jj;
					yy = wbnd[2] + jj*h->inc[GMT_Y] - y[n];
					weight = exp(-(xx2 + yy*yy)*factor);
					norm[kgrid] += (float)weight;
					C->grid_aux[kgrid] += (float)(weight * z[n]);
					C->num[kgrid]++;
					if (ii == ix && jj == iy) {	/* In this bin we don't want anyone else */
						norm[kgrid] = 1.0;
						C->grid_aux[kgrid] = (float)z[n];
						C->num[kgrid] = 1;
						cnt[kgrid] = 1;
					}
				}
			}
			ndata++;
		}
		else if (ix >= 0 && ix < gxdim && iy >= 0 && iy < gydim) {
			kgrid = ix*gydim + iy;
			if (C->num[kgrid] <= 0) {
				norm[kgrid] = 1.0;
				C->grid_aux[kgrid] = (float)z[n];
				C->num[kgrid] = 1;
				cnt[kgrid] = 1;
			}
			ndata++;
		}
	}

	gmt_M_free(GMT, x);
	gmt_M_free(GMT, y);
	gmt_M_free(GMT, z);

	return(ndata);
}
