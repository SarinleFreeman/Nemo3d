/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/util/util_argv_B.h,v 1.4 2004/10/27 19:43:10 gekco Exp
$
*****************************************************************************/

/*
 * Defines for a generic argv and argc processor...
 *
 * Copyright 1995 by Gray Watson
 *
 * This file is part of the argv library.
 *
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the
 * above copyright notice and this permission notice appear in all
 * copies, and that the name of Gray Watson not be used in advertising
 * or publicity pertaining to distribution of the document or software
 * without specific, written prior permission.
 *
 * Gray Watson makes no representations about the suitability of the
 * software described herein for any purpose.  It is provided "as is"
 * without express or implied warranty.
 *
 * The author may be contacted at gray.watson@letters.com
 *
 * $Id: util_argv_B.h,v 1.4 2004/10/27 19:43:10 gekco Exp $
 */

#ifndef __UTIL_ARGV_B_H__
#define __UTIL_ARGV_B_H__

/*
 * Version string for the library
 *
 * NOTE to gray: whenever this is changed, a corresponding Changlog
 * entry *must* be entered and 2 entries in argv.texi must be updated.
 *
 * ARGV LIBRARY VERSION -- 2.1.0
 */

/* NOTE: start of $Id: util_argv_B.h,v 1.4 2004/10/27 19:43:10 gekco Exp $ */

/*
 * Generic and standardized argument processor.  You describe the arguments
 * that you are looking for along with their types and these routines do the
 * work to convert them into values.
 *
 * These routines also provide standardized error and usage messages as well
 * as good usage documentation and long and short options.
 */

#include <fcntl.h>
#include <stdio.h> /* have to for FILE * below */
#include <stdlib.h>
#include <string.h> /* Linux warnings #include <strings.h> */
#include <termios.h>

/* this defines what type the standard void memory-pointer is */
/*#ifndef AIX_SPECIFIC
#if (defined(__STDC__) && __STDC__ == 1)
#else
#if (defined(__STDC__) && __STDC__ == 1) || (defined __cplusplus)
#endif*/

#if (defined(__STDC__) && __STDC__ == 1) || (defined __cplusplus)
#define ARGV_PNT void *
#else
#define ARGV_PNT char *
#endif

/*
 * argument information structure.  this specifies the allowable options
 * and some information about each one.
 *
 * { 'O',  "optimize",  ARGV_BOOL,  &optimize,  NULL,  "turn on optimization" }
 * { 'c',  "config",  ARGV_CHARP,  &config,  "file",  "configuration file" }
 */
typedef struct {
  char ar_short_arg;    /* the char of the arg, 'd' if '-d' */
  char *ar_long_arg;    /* long version of arg, 'delete' */
  short ar_type;        /* type of option, see values below */
  ARGV_PNT ar_variable; /* address of variable that is arg */
  char *ar_var_label;   /* label for variable descriptions */
  char *ar_comment;     /* comment for usage message */
} argv_t;

/*
 * argument array type.  when ARGV_ARRAY is |'d with the ar_type in the above
 * structure then multiple instances of the option are allowed and each
 * instance is stored into the following structure that MUST be in ar_variable
 * in the above arg_t structure.
 * NOTE: after the arguments have been processed, if aa_entryn is > 0 then
 * aa_entries needs to be free'd by user. argv_cleanup() can be used for this
 */
typedef struct {
  int aa_entryn;       /* number of elements in aa_entrees */
  ARGV_PNT aa_entries; /* entry list specified */
} argv_array_t;

/*  extract the count of the elements from an argv ARRAY */
#define ARGV_ARRAY_COUNT(array) ((array).aa_entryn)

/* extract WHICH entry of TYPE from an argv ARRAY */
#define ARGV_ARRAY_ENTRY(array, type, which)                                   \
  (((type *)(array).aa_entries)[which])

/* special ar_short_arg value to mark the last entry in the argument array */
#define ARGV_LAST ((char)255)

/*
 * special ar_short_arg value to mark mandatory arguments (i.e. arguments that
 * *must* be specified.  for arguments that are not optional like [-b].
 * to have a variable number of mandatory args then make the last MAND
 * entry be a ARG_ARRAY type.
 */
#define ARGV_MAND ((char)254)

/*
 * special ar_short_arg value to mark that there is the possibility of
 * a mandatory argument here if one is specified.
 */
#define ARGV_MAYBE ((char)253)

/*
 * special ar_short_arg value to say that the previous and next arguments in
 * the list should not be used together.
 * {'a'...}, {ARG_OR}, {'b'...}, {ARG_OR}, {'c'...} means
 * the user should only specific -a or -b or -c but not 2 or more.
 */
#define ARGV_OR ((char)252)

/*
 * special ar_short_arg value that is the same as ARGV_OR but one of the args
 * must be used.
 * {'a'...}, {ARG_OR}, {'b'...}, {ARG_OR}, {'c'...} means
 * the user must specify one of -a or -b or -c but not 2 or more.
 */
#define ARGV_XOR ((char)251)

/*
 * ar_type values of arg_t
 * NOTE: if this list is changed, some defines in argv_loc need to be changed
 */
#define ARGV_BOOL 1     /* boolean type, sets to ARGV_TRUE */
#define ARGV_BOOL_NEG 2 /* like bool but sets to ARGV_FALSE */
#define ARGV_BOOL_ARG 3 /* like bool but takes a yes/no arg */
#define ARGV_CHAR 4     /* single character */
#define ARGV_CHARP 5    /* same as STRING */
#define ARGV_STRING 5   /* character string */
#define ARGV_FLOAT 6    /* floating pointer number */
#define ARGV_SHORT 7    /* integer number */
#define ARGV_INT 8      /* integer number */
#define ARGV_U_INT 9    /* unsigned integer number */
#define ARGV_LONG 10    /* long integer number */
#define ARGV_U_LONG 11  /* unsinged long integer number */
#define ARGV_BIN 12     /* binary number (0s and 1s) */
#define ARGV_OCT 13     /* octal number, (base 8) */
#define ARGV_HEX 14     /* hexadecimal number, (base 16) */
#define ARGV_INCR 15    /* int arg which gets ++ each time */

#define ARGV_TYPE(t) ((t) & 0x3F) /* strip off all but the var type */
#define ARGV_ARRAY (1 << 14)      /* OR with type to indicate array */

/* argv_usage which argument values */
#define ARGV_USAGE_SHORT 1   /* print short usage messages */
#define ARGV_USAGE_LONG 2    /* print long-format usage messages */
#define ARGV_USAGE_DEFAULT 3 /* default usage messages */

/* boolean type settings */
#define ARGV_FALSE 0
#define ARGV_TRUE 1

/* global variable and procedure scoping for code readability */
#undef IMPORT
#define IMPORT extern

/*<<<<<<<<<<  The below prototypes are auto-generated by fillproto */

/* this is a processed version of argv[0], pre-path removed: /bin/ls -> ls */
IMPORT char argv_program[/* PROGRAM_NAME + 1 */];

/* a global value of argv from main after argv_process has been called */
IMPORT char **argv_argv;

/* a global value of argc from main after argv_process has been called */
IMPORT int argv_argc;

/* this should be set externally to provide general program help to user */
IMPORT char *argv_help_string;

/* this should be set externally to provide version information to the user */
IMPORT char *argv_version_string;

/*
 * are we running interactively?  this will exit on errors.  set to
 * false to return error codes instead.
 */
IMPORT char argv_interactive;

/*
 * the FILE stream that argv outputs all its errors.  set to NULL to
 * not dump any error messages.
 */
/* IMPORT	FILE 	*argv_error_stream; */

/*
 * processes ARGC number of arguments from ARGV depending on argument
 * info array ARGS (if null then an empty array is used).  this
 * routine will not modify the argv array in any way.  NOTE: it will
 * modify the args array by setting various flags in the type field.
 * returns 0 if no error else -1.
 */
IMPORT int argv_process(argv_t *args, int argc, char **argv);

/*
 * processes arguments sent in via the STRING that a web-server might
 * send to program in ARG0.  returns 0 on noerror else -1.
 */
IMPORT int argv_web_process_string(argv_t *args, char *arg0, char *string);

/*
 * processes arguments sent in via the QUERY_STRING environmental
 * variable that a web-server might send to program in ARG0.  returns
 * 0 on noerror else -1.
 */
IMPORT int argv_web_process(argv_t *args, char *arg0);

/*
 * print the standard usage messages for argument array ARGS (if null
 * then an empty array is used).  WHICH chooses between long or short
 * messages (see argv.h).
 * NOTE: if this is called before argv_process then the program name
 * may be messed up.
 */
IMPORT void argv_usage(argv_t *args, int which);

/*
 * see if ARG argument was used in a previous call to argv_process on
 * ARGS, returns 1 if yes else 0
 */
IMPORT int argv_was_used(argv_t *args, char arg);

/*
 * frees up any allocations in ARGS that may have been done by
 * argv_process.  This should be done at the end of the program or
 * after all the arguments have been referenced.
 */
IMPORT void argv_cleanup(argv_t *args);

/*
 * copy all the args (after the 0th), one after the other, into BUF of
 * MAX_SIZE.  NOTE: you can get the 0th argument from argv_argv[0].
 */
IMPORT void argv_copy_args(char *buf, int max_size);

/*<<<<<<<<<<   This is end of the auto-generated output from fillproto. */

#endif /* ! __UTIL_ARGV_B_H__ */
