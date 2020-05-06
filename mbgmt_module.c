/*
 * $Id: mbgmt_module.c 352 2015-02-01 04:34:58Z pwessel $
 *
 * Copyright (c) 2015 by P. Wessel
 * See LICENSE.TXT file for copying and redistribution conditions.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; version 3 or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * Contact info: http://www.soest.hawaii.edu/PT/GSFML
 *--------------------------------------------------------------------
 *
 * mbgmt_module.c populates the local array g_mbgmt_module
 * with parameters such as name, group and purpose strings.
 * This file also contains the following convenience function to
 * display all module purposes:
 *
 *   void gmt_mbgmt_module_show_all (struct GMTAPI_CTRL *API);
 *
 * This function will be called by gmt --help
 */

#include "gmt_dev.h"
#include "gmt_internals.h"
#include "mbgmt_module.h"

/* Sorted array with information for all mbgmt modules */

static struct GMT_MODULEINFO modules[] = {
	{"mbareaclean", "mbareaclean", "mbgmt", "Automatically flag bad beams in swath sonar bathymetry data", ""},
	{"mbcontour",   "mbcontour", "mbgmt",  "Plot swath bathymetry, amplitude, or backscatter",  "CC(,>X}"},
	{"mbclean_j",   "mbclean_j","mbgmt", "Identifies and flags artifacts in swath sonar bathymetry data", ""},
	{"mbflags",     "mbflags","mbgmt", "Update the flags .esf file from external file", ""},
	{"mbgetdata",   "mbgetdata","mbgmt", "Get swath bathymetry, amplitude, or backscatter data into Matlab", "<D{,ND),MD}"},
	{"mbgrdtiff",   "mbgrdtiff", "mbgmt",  "Project grids or images and plot them on maps",  "<G{+,CC(,IG("},
	{"mbgrid_j",    "mbgrid_j","mbgmt", "grid bathymetry, amplitude, or sidescan data of a set of swath sonar data files", ""},
	{"mbimport",    "mbimport","mbgmt", "Get swath bathymetry, amplitude, or backscatter Image into Matlab", "<D{,CC(,MI}"},
	{"mbinfo",      "mbinfo","mbgmt", "Read a swath sonar data file and outputs some basic statistics.", "<D{,>D},>DC,>DO"},
	{"mbset",       "mbset","mbgmt", "Set values in an mbprocess parameter file", ""},
	{"mbprocess",   "mbprocess","mbgmt", "Tool for processing swath sonar bathymetry data", "ID{,OD("},
	{"mbswath",     "mbswath", "mbgmt",  "Plot swath bathymetry, amplitude, or backscatter",  "CC(,NC(,>X}"},
	{"gmtmbgrid",   "gmtmbgrid","mbgmt", "Grid table data using adjustable tension continuous curvature splines", "<D{,DD(,LG(,GG}"},
	{NULL, NULL, NULL, NULL, NULL} /* last element == NULL detects end of array */
};

#if 0
/* Pretty print all GMT mbgmt module names and their purposes */
EXTERN_MSC void gmtlib_mbgmt_module_show_all (void *API) {
	gmtlib_module_show_all (API, modules, GMT_MBGMT_STRING);
}

/* Produce single list on stdout of all MBGMT supplements module names for gmt --show-modules */
EXTERN_MSC void gmtlib_mbgmt_module_list_all (void *API) {
	gmtlib_module_list_all (API, modules);
}

/* Produce single list on stdout of all MBGMT supplements module names for gmt --show-classic [i.e., classic mode names] */
EXTERN_MSC void gmtlib_mbgmt_module_classic_all (void *API) {
	gmtlib_module_classic_all (API, modules);
}

/* Lookup module id by name, return option keys pointer (for external API developers) */
EXTERN_MSC const char *gmtlib_mbgmt_module_keys (void *API, char *candidate) {
	return (gmtlib_module_keys (API, modules, candidate));
}

/* Lookup module id by name, return group char name (for external API developers) */
EXTERN_MSC const char *gmtlib_mbgmt_module_group (void *API, char *candidate) {
	return (gmtlib_module_group (API, modules, candidate));
}
#endif
