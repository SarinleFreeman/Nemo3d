/*****************************************************************************
The JPL NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Library General Public License
along with this library.  If not, write to the Free Software Foundation,
Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

Written by Chris Bowen, Gerhard Klimeck, and Fabiano A. Oyafuso

*****************************************************************************
$Header: /repo/nemo3d/src/util/util_stringops.h,v 1.3 2004/10/27 19:43:10 gekco
Exp $
*****************************************************************************/

#ifndef UTIL_STRINGOPS_H
#define UTIL_STRINGOPS_H 1

#include <ctype.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define str_calloc calloc
#define str_malloc malloc
#define str_free free

#include "system.h"

char *n3d_strdup_n(const char *);
char *n3d_copy_str(const char *str);
char *n3d_strcatz(char *, ...);
int n3d_strcmp_null_check(const char *, const char *);
char *n3d_FileNameGet(char *);
void n3d_FileTypeSet(char **, char *, int);
int n3d_StringListNum(char **);
int n3d_StringListNum(char **);
char **n3d_StringListAdd(char ***, const char *);

#endif
