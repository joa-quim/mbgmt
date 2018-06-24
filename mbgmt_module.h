/*
 * $Id $
 * Copyright (c) 2018 by J. Luis
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
 */

/* mbgmt_module.h declares the prototypes for mbgmt module functions
 * and the array that contains MBGMT module parameters such as name
 * and purpose strings.
 */

#pragma once
#ifndef _GMT_MBGMT_MODULE_H
#define _GMT_MBGMT_MODULE_H

#ifdef __cplusplus /* Basic C++ support */
extern "C" {
#endif

/* Declaration modifiers for DLL support (MSC et al) */
#include "declspec.h"

/* Prototypes of all modules in the GMT core library */
EXTERN_MSC int GMT_mbareaclean (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbcontour (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbclean_j (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbflags   (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbgetdata (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbgrdtiff (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbimport  (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbinfo    (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbset     (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbprocess (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbswath   (void *API, int mode, void *args);
EXTERN_MSC int GMT_gmtmbgrid (void *API, int mode, void *args);

/* Pretty print all modules in the MBGMT library and their purposes */
EXTERN_MSC void gmt_mbgmt_module_show_all(void *API);

/* List all modules in the MBGMT library to stdout */
EXTERN_MSC void gmt_mbgmt_module_list_all(void *API);

/* Undocumented API function for developers to get information about a module */
EXTERN_MSC const char *gmt_mbgmt_module_info(void *API, char *candidate);

/* Lookup module id by name, return option keys pointer (for external API developers) */
EXTERN_MSC const char *gmt_mbgmt_module_keys(void *API, char *candidate);

/* Lookup module id by name, return group char name (for external API developers) */
EXTERN_MSC const char *gmt_mbgmt_module_group(void *API, char *candidate);

#ifdef __cplusplus
}
#endif

#endif
