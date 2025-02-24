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
$Header: /repo/nemo3d/src/util/util_argv_B.c,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

/*
 * Generic argv processor...
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
 */

#include <ctype.h>
#include <stdio.h>

#if STDC_HEADERS
# include <string.h>
# include <stdlib.h>
#endif

#include "util_argv_B.h"
#include "util_argv_loc_B.h"


/* internal routines */
EXPORT	void	argv_usage(argv_t * args, int which);

/*
 * exported variables
 */
/* this is a processed version of argv[0], pre-path removed: /bin/ls -> ls */
EXPORT	char	argv_program[PROGRAM_NAME + 1] = "Unknown";

/* a global value of argv from main after argv_process has been called */
EXPORT	char	**argv_argv		= NULL;
/* a global value of argc from main after argv_process has been called */
EXPORT	int	argv_argc		= 0;

/* this should be set externally to provide general program help to user */
EXPORT	char	*argv_help_string	= NULL;
/* this should be set externally to provide version information to the user */
EXPORT	char	*argv_version_string	= NULL;

/*
 * are we running interactively?  this will exit on errors.  set to
 * false to return error codes instead.
 */
EXPORT	char 	argv_interactive	= ARGV_TRUE;

/*
 * the FILE stream that argv outputs all its errors.  set to NULL to
 * not dump any error messages.
 */
/* EXPORT	FILE 	*argv_error_stream	= stdout; */
#define argv_error_stream stderr

/* local variables */
LOCAL_QUEUE_DECLARE(argv_t *);			/* args waiting for values */

LOCAL	argv_t	empty[] = { { ARGV_LAST }};	/* empty argument array */

/* global settings */
LOCAL	int	global_close	= GLOBAL_CLOSE_ENABLE;	/* close processing */
LOCAL	int	global_env	= GLOBAL_ENV_BEFORE;	/* env processing */
LOCAL	int	global_error	= GLOBAL_ERROR_SEE;	/* error processing */
LOCAL	int	global_multi	= GLOBAL_MULTI_ACCEPT;	/* multi processing */
LOCAL	int	global_usage	= GLOBAL_USAGE_LONG;	/* usage processing */

/***************************** general utilities *****************************/

/*
 * binary STR to integer translation
 */
LOCAL	int	btoi(char * str)
{
  int		ret = 0;
  
  /* strip off spaces */
  for (; isspace(*str); str++);
  
  for (; *str == '0' || *str == '1'; str++) {
    ret *= 2;
    ret += *str - '0';
  }
  
  return ret;
}

/*
 * octal STR to integer translation
 */
LOCAL	int	otoi(char *str)
{
  int		ret = 0;
  
  /* strip off spaces */
  for (; isspace(*str); str++);
  
  for (; *str >= '0' && *str <= '7'; str++) {
    ret *= 8;
    ret += *str - '0';
  }
  
  return ret;
}

/*
 * hexadecimal STR to integer translation
 */
LOCAL	int	htoi(char * str)
{
  int		ret = 0;
  
  /* strip off spaces */
  for (; isspace(*str); str++);
  
  /* skip a leading 0[xX] */
  if (*str == '0' && (*(str + 1) == 'x' || *(str + 1) == 'X'))
    str += 2;
  
  for (; isdigit(*str) ||
       (*str >= 'a' && *str <= 'f') || (*str >= 'A' && *str <= 'F');
       str++) {
    ret *= 16;
    if (*str >= 'a' && *str <= 'f')
      ret += *str - 'a' + 10;
    else if (*str >= 'A' && *str <= 'F')
      ret += *str - 'A' + 10;
    else
      ret += *str - '0';
  }
  
  return ret;
}

/*
 * basically a strdup_n for compatibility sake
 */
LOCAL	char	*string_copy(char * ptr)
{
  char	*ptrp;
  char		*ret, *retp;
  int		len;
  
  len = strlen(ptr);
  ret = (char *)malloc(len + 1);
  if (ret != NULL) {
    for (ptrp = ptr, retp = ret; *ptrp != NULLC;)
      *retp++ = *ptrp++;
    *retp = NULLC;
  }
  
  return ret;
}

/*
 * break STR and return an array of char * whose values are tokenized by TOK.
 * it passes back the number of tokens in TOKN.
 * NOTE: the return value should be freed later and the STR should stay around
 * until that time.
 */
LOCAL	char	**vectorize(char * str, char * tok, int * tokn)
{
  char	**vectp;
  char	*tmp, *tokp;
  int	tokc;
  
  /* count the tokens */
  tmp = string_copy(str);
  tokp = strtok(tmp, tok);
  for (tokc = 0; tokp != NULL; tokc++)
    tokp = strtok(NULL, tok);
  free(tmp);
  
  *tokn = tokc;
  
  if (tokc == 0)
    return NULL;
  
  /* allocate the pointer grid */
  vectp = (char **)malloc(sizeof(char *) * tokc);
  
  /* load the tokens into the list */
  vectp[0] = strtok(str, tok);
  
  for (tokc = 1; tokc < *tokn; tokc++)
    vectp[tokc] = strtok(NULL, tok);
  
  return vectp;
}

/*
 * display printable chars from BUF of SIZE, non-printables as \%03o
 */
LOCAL	char	*expand_buf(void * buf, int size)
{
  static char	out[DUMP_SPACE_BUF];
  int		sizec;
  void		*bufp;
  char	 	*outp;
  
  for (sizec = 0, outp = out, bufp = (void *)buf; sizec < size;
       sizec++, bufp = (char *)bufp + 1) {
    char	*specp;
    
    /* handle special chars */
    if (outp + 2 >= out + sizeof(out))
      break;
    
    /* search for special characters */
    for (specp = (char*)(SPECIAL_CHARS + 1); *(specp - 1) != NULLC; specp += 2)
      if (*specp == *(char *)bufp)
	break;
    
    /* did we find one? */
    if (*(specp - 1) != NULLC) {
      if (outp + 2 >= out + sizeof(out))
	break;
      (void)sprintf(outp, "\\%c", *(specp - 1));
      outp += 2;
      continue;
    }
    
    if (*(unsigned char *)bufp < 128 && isprint(*(char *)bufp)) {
      if (outp + 1 >= out + sizeof(out))
	break;
      *outp = *(char *)bufp;
      outp += 1;
    }
    else {
      if (outp + 4 >= out + sizeof(out))
	break;
      (void)sprintf(outp, "\\%03o", *(unsigned char *)bufp);
      outp += 4;
    }
  }
  
  *outp = NULLC;
  return out;
}

/****************************** usage routines *******************************/

/*
 * print a short-format usage message
 */
LOCAL	void	usage_short(argv_t * args, int flag)
{
  argv_t	*argp;
  int		len, colc = 0;
  char		mark = ARGV_FALSE, *prefix;
  
  if (argv_error_stream == NULL)
    return;
  
  /* print the usage message header */
  (void)fprintf(argv_error_stream, "%s%s", USAGE_LABEL, argv_program);
  colc += USAGE_LABEL_LENGTH + strlen(argv_program);
  
  /*
   * print all of the boolean arguments first.
   * NOTE: we assume they all fit on the line
   */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++) {
    
    /* skip or-specifiers */
    if (argp->ar_short_arg == ARGV_OR
	|| argp->ar_short_arg == ARGV_XOR)
      continue;
    
    /* skip non booleans */
    if (HAS_ARG(argp->ar_type))
      continue;
    
    /* skip args with no short component */
    if (argp->ar_short_arg == NULLC)
      continue;
    
    if (! mark) {
      len = 2 + SHORT_PREFIX_LENGTH;
      prefix = " [";
      
      /* we check for -2 here because we should have 1 arg and ] on line */
      if (colc + len > SCREEN_WIDTH - 2) {
	(void)fprintf(argv_error_stream, "\n%*.*s",
		      (int)USAGE_LABEL_LENGTH, (int)USAGE_LABEL_LENGTH, "");
	colc = USAGE_LABEL_LENGTH;
	
	/* if we are the start of a line, skip any starting spaces */
	if (*prefix == ' ') {
	  prefix++;
	  len--;
	}
      }
      
      (void)fprintf(argv_error_stream, "%s%s", prefix, SHORT_PREFIX);
      colc += len;
      mark = ARGV_TRUE;
    }
    
    len = 1;
    /* we check for -1 here because we should need ] */
    if (colc + len > SCREEN_WIDTH - 1) {
      (void)fprintf(argv_error_stream, "]\n%*.*s",
		    (int)USAGE_LABEL_LENGTH, (int)USAGE_LABEL_LENGTH, "");
      colc = USAGE_LABEL_LENGTH;
      
      /* restart the short option list */
      (void)fprintf(argv_error_stream, "[%s", SHORT_PREFIX);
      colc += 1 + SHORT_PREFIX_LENGTH;
    }
    
    (void)fprintf(argv_error_stream, "%c", argp->ar_short_arg);
    colc++;
  }
  
  if (mark) {
    (void)fprintf(argv_error_stream, "]");
    colc++;
  }
  
  /* print remaining (non-boolean) arguments */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++) {
    int		var_len;
    char	*var_str, *postfix;
    
    /* skip or-specifiers */
    if (argp->ar_short_arg == ARGV_OR
	|| argp->ar_short_arg == ARGV_XOR)
      continue;
    
    /* skip booleans types */
    if (! HAS_ARG(argp->ar_type))
      continue;
    
    if (argp->ar_var_label == NULL) {
      if (ARGV_TYPE(argp->ar_type) == ARGV_BOOL_ARG) {
	var_str = BOOL_ARG_LABEL;
	var_len = BOOL_ARG_LENGTH;
      }
      else {
	var_len = UNKNOWN_ARG_LENGTH;
	var_str = UNKNOWN_ARG;
      }
    }
    else {
      var_len = strlen(argp->ar_var_label);
      var_str = argp->ar_var_label;
    }
    
    if (argp->ar_short_arg == ARGV_MAND) {
      /* print the mandatory argument desc */
      len = 1 + var_len;
      prefix = " ";
      postfix = "";
    }
    else if (argp->ar_short_arg == ARGV_MAYBE) {
      /* print the maybe argument desc */      
      len = 2 + var_len + 1;
      prefix = " [";
      postfix = "]";
    }
    else {
      /* handle options with arguments */
      
      /* " [" + short_prefix + char */
      len = 2 + SHORT_PREFIX_LENGTH + 1;
      prefix = " [";
      
      /* do we need to wrap */
      if (colc + len > SCREEN_WIDTH) {
	(void)fprintf(argv_error_stream, "\n%*.*s",
		      (int)USAGE_LABEL_LENGTH, (int)USAGE_LABEL_LENGTH, "");
	colc = USAGE_LABEL_LENGTH;
	
	/* if we are the start of a line, skip any starting spaces */
	if (*prefix == ' ') {
	  prefix++;
	  len--;
	}
      }
      (void)fprintf(argv_error_stream, "%s%s%c",
		    prefix, SHORT_PREFIX, argp->ar_short_arg);
      colc += len;
      
      len = 1 + var_len + 1;
      prefix = " ";
      postfix = "]";
    }
    
    if (colc + len > SCREEN_WIDTH) {
      (void)fprintf(argv_error_stream, "\n%*.*s",
		    (int)USAGE_LABEL_LENGTH, (int)USAGE_LABEL_LENGTH, "");
      colc = USAGE_LABEL_LENGTH;
      
      /* if we are the start of a line, skip any starting spaces */
      if (*prefix == ' ') {
	prefix++;
	len--;
      }
    }
    
    (void)fprintf(argv_error_stream, "%s%s%s", prefix, var_str, postfix);
    colc += len;
  }
  
  (void)fprintf(argv_error_stream, "\n");
  
  if (flag == GLOBAL_USAGE_SHORTREM)
    (void)fprintf(argv_error_stream,
		  "%*.*sUse the '%s%s' argument for more assistance.\n",
		  (int)USAGE_LABEL_LENGTH, (int)USAGE_LABEL_LENGTH, "",
		  LONG_PREFIX, USAGE_ARG);
}


/*
 * print a short-format usage message
 */
LOCAL	void	display_arg(FILE * stream, argv_t * argp, int max,
			    int * colc)
{
  int	var_len, len;
  
  if (argp->ar_var_label == NULL)
    var_len = 0;
  else
    var_len = strlen(argp->ar_var_label);
  
  switch (ARGV_TYPE(argp->ar_type)) {
    
  case ARGV_BOOL:
  case ARGV_BOOL_NEG:
  case ARGV_INCR:
    break;
    
  case ARGV_BOOL_ARG:
    (void)fprintf(stream, "%s", BOOL_ARG_LABEL);
    (*colc) += BOOL_ARG_LENGTH;
    break;
    
  case ARGV_CHAR:
  case ARGV_CHARP:
  case ARGV_FLOAT:
  case ARGV_SHORT:
  case ARGV_INT:
  case ARGV_U_INT:
  case ARGV_LONG:
  case ARGV_U_LONG:
  case ARGV_BIN:
  case ARGV_OCT:
  case ARGV_HEX:
    if (argp->ar_var_label == NULL) {
      len = max - *colc;
      (void)fprintf(stream, "%-.*s", len, UNKNOWN_ARG);
      *colc += MIN(len, UNKNOWN_ARG_LENGTH);
    }
    else {
      len = max - *colc;
      (void)fprintf(stream, "%-.*s", len, argp->ar_var_label);
      *colc += MIN(len, var_len);
    }
    break;
  }
}

/*
 * display an option entry ARGP to STREAM while counting COLC
 */
LOCAL	void	display_option(FILE * stream, argv_t * argp, int * colc)
{
  if (stream == NULL)
    return;
  
  (void)fputc('[', stream);
  (*colc)++;
  
  /* arg maybe does not have a -? preface */
  if (argp->ar_short_arg != ARGV_MAYBE) {
    (void)fprintf(stream, "%s%c",
		  SHORT_PREFIX, argp->ar_short_arg);
    *colc += SHORT_PREFIX_LENGTH + 1;
    
    if (HAS_ARG(argp->ar_type)) {
      /* display optional argument */
      (void)fputc(' ', stream);
      (*colc)++;
    }
  }
  
  display_arg(stream, argp, LONG_COLUMN - 1, colc);
  (void)fputc(']', stream);
  (*colc)++;
}

/*
 * print a long-format usage message
 */
LOCAL	void	usage_long(argv_t * args)
{
  argv_t	*argp;
  int		colc, len;
  
  if (argv_error_stream == NULL)
    return;
  
  /* print the usage message header */
  (void)fprintf(argv_error_stream, "%s%s\n", USAGE_LABEL, argv_program);
  
  /* run through the argument structure */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++) {
    
    /* skip or specifiers */
    if (argp->ar_short_arg == ARGV_OR || argp->ar_short_arg == ARGV_XOR)
      continue;
    
    /* indent to the short-option colc */
    (void)fprintf(argv_error_stream, "%*.*s", SHORT_COLUMN, SHORT_COLUMN, "");
    
    /* start column counter */
    colc = SHORT_COLUMN;
    
    /* print the short-arg stuff if there */
    if (argp->ar_short_arg == NULLC) {
      (void)fputc('[', argv_error_stream);
      colc++;
    }
    else {
      if (argp->ar_short_arg == ARGV_MAND)
	display_arg(argv_error_stream, argp, COMMENT_COLUMN, &colc);
      else
	display_option(argv_error_stream, argp, &colc);
      
      /* put the long-option message on the correct column */
      if (colc < LONG_COLUMN) {
	(void)fprintf(argv_error_stream, "%*.*s",
		      LONG_COLUMN - colc, LONG_COLUMN - colc, "");
	colc = LONG_COLUMN;
      }
    }
    
    /* print the long-option message */
    if (argp->ar_long_arg != NULL) {
      len = COMMENT_COLUMN - colc - (LONG_PREFIX_LENGTH + 1);
      if (argp->ar_short_arg != NULLC) {
	(void)fprintf(argv_error_stream, "%s", LONG_LABEL);
	colc += LONG_LABEL_LENGTH;
	len -= LONG_LABEL_LENGTH;
      }
      (void)fprintf(argv_error_stream, "%s%-.*s",
		    LONG_PREFIX, len, argp->ar_long_arg);
      colc += LONG_PREFIX_LENGTH + MIN(len, strlen(argp->ar_long_arg));
    }
    
    /* add the optional argument if no short-arg */
    if (argp->ar_short_arg == NULLC) {
      if (HAS_ARG(argp->ar_type)) {
	(void)fputc(' ', argv_error_stream);
	colc++;
      }
      
      /* display any optional arguments */
      display_arg(argv_error_stream, argp, COMMENT_COLUMN - 1, &colc);
      (void)fputc(']', argv_error_stream);
      colc++;
    }
    
    /* print the comment */
    if (argp->ar_comment != NULL) {
      /* put the comment message on the correct column */
      if (colc < COMMENT_COLUMN) {
	(void)fprintf(argv_error_stream, "%*.*s",
		      COMMENT_COLUMN - colc,
		      COMMENT_COLUMN - colc, "");
	colc = COMMENT_COLUMN;
      }
      
      len = SCREEN_WIDTH - colc - COMMENT_LABEL_LENGTH;
      (void)fprintf(argv_error_stream, "%s%-.*s",
		    COMMENT_LABEL, len, argp->ar_comment);
    }
    
    (void)fprintf(argv_error_stream, "\n");
  }
}

/*
 * do the usage depending on FLAG
 */
LOCAL	void	do_usage(argv_t * args, int flag)
{
  if (argv_error_stream == NULL)
    return;
  
  if (flag == GLOBAL_USAGE_SEE)
    (void)fprintf(argv_error_stream,
		  "%*.*sUse the '%s%s' argument for assistance.\n",
		  (int)USAGE_LABEL_LENGTH, (int)USAGE_LABEL_LENGTH, "",
		  LONG_PREFIX, USAGE_ARG);
  else if (flag == GLOBAL_USAGE_SHORT || flag == GLOBAL_USAGE_SHORTREM)
    usage_short(args, flag);
  else if (flag == GLOBAL_USAGE_LONG || flag == GLOBAL_USAGE_ALL)
    usage_long(args);
  
  if (flag == GLOBAL_USAGE_ALL) {
    (void)fprintf(argv_error_stream, "\n");
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' for default usage information.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, USAGE_ARG);
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' for short usage information.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, USAGE_SHORT_ARG);
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' for long usage information.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, USAGE_LONG_ARG);
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' for all usage information.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, USAGE_ALL_ARG);
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' to display the help message.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, HELP_ARG);
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' to display the version message.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, VERSION_ARG);
    (void)fprintf(argv_error_stream,
		  "%*.*sUse '%s%s' to display the options and their values.\n",
		  SHORT_COLUMN, SHORT_COLUMN, "",
		  LONG_PREFIX, DISPLAY_ARG);
  }
}

/******************************* preprocessing *******************************/

/*
 * preprocess argument array ARGS of NUM_ARGS entries and set the MAND
 * and MAYBE boolean arrays.  returns [NO]ERROR
 */
LOCAL	int	preprocess_array(argv_t * args, int num_args)
{
  argv_t	*argp;
  char		mand_array = ARGV_FALSE, maybe_field = ARGV_FALSE;
  
  /* count the args and find the first mandatory */
  for (argp = args; argp < args + num_args; argp++) {
    
    /* clear internal flags */
    argp->ar_type &= ~ARGV_FLAG_USED;
    
    /* do we have a mandatory-array? */
    if (argp->ar_short_arg == ARGV_MAND) {
      if (mand_array) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, no ARGV_MAND's can follow a MAND or MAYBE array\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
      if (maybe_field) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, no ARGV_MAND's can follow a ARGV_MAYBE\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
      
#if 0
      if (argp->ar_long_arg != NULL) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, ARGV_MAND's should not have long-options\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
#endif
      
      if (argp->ar_type & ARGV_ARRAY)
	mand_array = ARGV_TRUE;
    }
    
    /* do we have a maybe field? */
    if (argp->ar_short_arg == ARGV_MAYBE) {
      if (mand_array) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, no ARGV_MAYBE's can follow a MAND or MAYBE array\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
      
      maybe_field = ARGV_TRUE;
      if (argp->ar_type & ARGV_ARRAY)
	mand_array = ARGV_TRUE;
    }
    
    /* handle initializing the argument array */
    if (argp->ar_type & ARGV_ARRAY) {
      if (! HAS_ARG(argp->ar_type)) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, cannot have an array of boolean values\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
      if (ARGV_TYPE(argp->ar_type) == ARGV_INCR) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, cannot have an array of incremental values\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
#if 0 /* dkb 10/22/97 */
      arrp->aa_entryn = 0;
#endif /* dkb 10/22/97 */
    }
    
#if 0
    /* must have a valid ar_short_arg */
    if (argp->ar_short_arg == NULLC) {
      if (argv_error_stream != NULL)
	(void)fprintf(argv_error_stream,
		      "%s: %s, short-option character is '\\0'\n",
		      argv_program, INTERNAL_ERROR_NAME);
      if (argv_interactive)
	(void)exit(EXIT_CODE);
      return ERROR;
    }
#endif
    
    /* verify variable pointer */
    if (argp->ar_variable == NULL
	&& argp->ar_short_arg != ARGV_OR
	&& argp->ar_short_arg != ARGV_XOR) {
      if (argv_error_stream != NULL)
	(void)fprintf(argv_error_stream,
		      "%s: %s, NULL variable specified in arg array\n",
		      argv_program, INTERNAL_ERROR_NAME);
      if (argv_interactive)
	(void)exit(EXIT_CODE);
      return ERROR;
    }
    
    /* verify [X]OR's */
    if (argp->ar_short_arg == ARGV_OR
	|| argp->ar_short_arg == ARGV_XOR) {
      
      /* that they are not at the start or end of list */
      if (argp == args || argp >= (args + num_args - 1)) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, ARGV_[X]OR entries cannot be at start or end of array\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
      
      /* that two aren't next to each other */
      if ((argp - 1)->ar_short_arg == ARGV_OR
	  || (argp - 1)->ar_short_arg == ARGV_XOR) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, two ARGV_[X]OR entries cannot be next to each other\n",
			argv_program, INTERNAL_ERROR_NAME);
	if (argv_interactive)
	  (void)exit(EXIT_CODE);
	return ERROR;
      }
    }
  }
  
  return NOERROR;
}

/*
 * translate string argument ARG into VAR depending on its TYPE
 */
LOCAL	void	translate_value(char * arg, ARGV_PNT var,
				short type)
{
  argv_array_t	*arrp;
  argv_type_t	*typep;
  int		val_type = ARGV_TYPE(type), size = 0;
  
  /* find the type and the size for array */
  for (typep = argv_types; typep->at_value != 0; typep++)
    if (typep->at_value == val_type) {
      size = typep->at_size;
      break;
    }
  
  if (typep->at_value == 0) {
    if (argv_error_stream != NULL)
      (void)fprintf(argv_error_stream, "%s: illegal variable type %d\n",
		    __FILE__, val_type);
    return;
  }
  
  if (type & ARGV_ARRAY) {
    arrp = (argv_array_t *)var;
    
    if (arrp->aa_entryn == 0)
      arrp->aa_entries = (char *)malloc(ARRAY_INCR * size);
    else if (arrp->aa_entryn % ARRAY_INCR == 0)
      arrp->aa_entries =
	(char *)realloc(arrp->aa_entries, (arrp->aa_entryn + ARRAY_INCR) *
			size);
    
    var = (char *)(arrp->aa_entries) + arrp->aa_entryn * size;
    arrp->aa_entryn++;
  }
  
  /* translate depending on type */
  switch (val_type) {
    
  case ARGV_BOOL:
    /* if no close argument, set to true */
    if (arg == NULL)
      *(char *)var = ARGV_TRUE;
    else if (*(char *)arg == 't' || *(char *)arg == 'T'
	     || *(char *)arg == 'y' || *(char *)arg == 'Y'
	     || *(char *)arg == '1')
      *(char *)var = ARGV_TRUE;
    else
      *(char *)var = ARGV_FALSE;
    break;
    
  case ARGV_CHAR:
    *(char *)var = *(char *)arg;
    break;
    
  case ARGV_CHARP:
    *(char **)var = string_copy((char *)arg);
    break;
    
  case ARGV_FLOAT:
    *(float *)var = (float)atof(arg);
    break;
    
  case ARGV_SHORT:
    *(short *)var = (short)atoi(arg);
    break;
    
  case ARGV_INT:
    *(int *)var = atoi(arg);
    break;
    
  case ARGV_U_INT:
    *(unsigned int *)var = atoi(arg);
    break;
    
  case ARGV_LONG:
    *(long *)var = atol(arg);
    break;
    
  case ARGV_U_LONG:
    *(unsigned long *)var = atol(arg);
    break;
    
  case ARGV_BIN:
    *(int *)var = btoi(arg);
    break;
    
  case ARGV_OCT:
    *(int *)var = otoi(arg);
    break;
    
  case ARGV_HEX:
    *(int *)var = htoi(arg);
    break;
    
  case ARGV_BOOL_NEG:
    /* if no close argument, set to false */
    if (arg == NULL)
      *(char *)var = ARGV_FALSE;
    else if (*(char *)arg == 't' || *(char *)arg == 'T'
	     || *(char *)arg == 'y' || *(char *)arg == 'Y'
	     || *(char *)arg == '1')
      *(char *)var = ARGV_TRUE;
    else
      *(char *)var = ARGV_FALSE;
    break;
    
  case ARGV_INCR:
    /* if no close argument then increment else set the value */
    if (arg == NULL)
      (*(int *)var)++;
    else
      *(int *)var = atoi(arg);
    break;
    
  case ARGV_BOOL_ARG:
    if (*(char *)arg == 't' || *(char *)arg == 'T'
	|| *(char *)arg == 'y' || *(char *)arg == 'Y'
	|| *(char *)arg == '1')
      *(char *)var = ARGV_TRUE;
    else
      *(char *)var = ARGV_FALSE;
    break;
  }
}

/*
 * translate value from VAR into string STR depending on its TYPE
 */
LOCAL	void	display_value(ARGV_PNT var, short type)
{
  int		val_type = ARGV_TYPE(type);
  
  /* translate depending on type */
  switch (val_type) {
    
  case ARGV_BOOL:
  case ARGV_BOOL_NEG:
  case ARGV_BOOL_ARG:
    if (*(char *)var)
      (void)printf("ARGV_TRUE");
    else
      (void)printf("ARGV_FALSE");
    break;
    
  case ARGV_CHAR:
    (void)printf("'%s'", expand_buf((char *)var, 1));
    break;
    
  case ARGV_CHARP:
    {
      int	len;
      if (*(char **)var == NULL)
	(void)printf("(null)");
      else {
	len = strlen(*(char **)var);
	(void)printf("\"%s\"", expand_buf(*(char **)var, len));
      }
    }
    break;
    
  case ARGV_FLOAT:
    (void)printf("%f", *(float *)var);
    break;
    
  case ARGV_SHORT:
    (void)printf("%d", *(short *)var);
    break;
    
  case ARGV_INT:
  case ARGV_INCR:
    (void)printf("%d", *(int *)var);
    break;
    
  case ARGV_U_INT:
    (void)printf("%u", *(unsigned int *)var);
    break;
    
  case ARGV_LONG:
    (void)printf("%ld", *(long *)var);
    break;
    
  case ARGV_U_LONG:
    (void)printf("%lu", *(unsigned long *)var);
    break;
    
    /* this should be a routine */
  case ARGV_BIN:
    {
      int	bitc;
      char	first = ARGV_FALSE;
      
      (void)fputc('0', argv_error_stream);
      
      if (*(int *)var != 0) {
	(void)fputc('b', argv_error_stream);
	
	for (bitc = sizeof(int) * BITS_IN_BYTE - 1; bitc >= 0; bitc--) {
	  int	bit = *(int *)var & (1 << bitc);
	  
	  if (bit == 0) {
	    if (first)
	      (void)fputc('0', argv_error_stream);
	  }
	  else {
	    (void)fputc('1', argv_error_stream);
	    first = ARGV_TRUE;
	  }
	}
      }
    }
    break;
    
  case ARGV_OCT:
    (void)printf("%#o", *(int *)var);
    break;
    
  case ARGV_HEX:
    (void)printf("%#x", *(int *)var);
    break;
  }
}

/*
 * translate value from VAR into string STR depending on its TYPE
 */
LOCAL	void	display_variables(argv_t * args)
{
  argv_t	*argp;
  argv_type_t	*typep;
  
  /* run through the argument structure */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++) {
    int		colc, val_type = ARGV_TYPE(argp->ar_type);
    
    /* skip or specifiers */
    if (argp->ar_short_arg == ARGV_OR || argp->ar_short_arg == ARGV_XOR)
      continue;
    
    colc = 0;
    if (argp->ar_short_arg == ARGV_MAND)
      display_arg(argv_error_stream, argp, COMMENT_COLUMN, &colc);
    else
      display_option(argv_error_stream, argp, &colc);
    
    /* put the type in the correct column */
    if (colc < LONG_COLUMN) {
      (void)printf("%*.*s", LONG_COLUMN - colc, LONG_COLUMN - colc, "");
      colc = LONG_COLUMN;
    }
    
    /* find the type */
    typep = NULL;
    for (typep = argv_types; typep->at_value != 0; typep++) {
      if (typep->at_value == ARGV_TYPE(argp->ar_type)) {
	int	len, tlen;
	
	len = COMMENT_COLUMN - colc - 1;
	tlen = strlen(typep->at_name);
	(void)printf(" %-.*s", len, typep->at_name);
	colc += MIN(len, tlen);
	if (argp->ar_type & ARGV_ARRAY) {
	  (void)printf("%s", ARRAY_LABEL);
	  colc += sizeof(ARRAY_LABEL) - 1;
	}
	break;
      }
    }
    
    if (colc < COMMENT_COLUMN) {
      (void)printf("%*.*s", COMMENT_COLUMN - colc, COMMENT_COLUMN - colc, "");
      colc = COMMENT_COLUMN;
    }
    
    if (argp->ar_type & ARGV_ARRAY) {
      argv_array_t	*arrp;
      int		entryc, size = 0;
      
      /* find the type and the size for array */
      if (typep == NULL) {
	(void)printf("%s: illegal variable type %d\n", __FILE__, val_type);
	continue;
      }
      size = typep->at_size;
      arrp = (argv_array_t *)argp->ar_variable;
      
      if (arrp->aa_entryn == 0)
	(void)printf("no entries");
      else {
	for (entryc = 0; entryc < arrp->aa_entryn; entryc++) {
	  ARGV_PNT	var;
	  if (entryc > 0)
	    (void)fputc(',', argv_error_stream);
	  var = (char *)(arrp->aa_entries) + entryc * size;
	  display_value(var, val_type);
	}
      }
    }
    else
      display_value(argp->ar_variable, val_type);
    (void)fputc('\n', argv_error_stream);
  }
}

/************************** checking used arguments **************************/

/*
 * check out if WHICH argument from ARGS has an *or* specified
 * attached to it.  returns [NO]ERROR
 */
LOCAL	int	check_or(argv_t * args, argv_t * which)
{
  argv_t	*argp, *matchp = NULL;
  
  /* check ORs below */
  for (argp = which - 2; argp >= args; argp -= 2) {
    if ((argp + 1)->ar_short_arg != ARGV_OR
	&& (argp + 1)->ar_short_arg != ARGV_XOR)
      break;
    if (argp->ar_type & ARGV_FLAG_USED) {
      matchp = argp;
      break;
    }
  }
  
  /* check ORs above */
  if (matchp == NULL) {
    /* NOTE: we assume that which is not pointing now to ARGV_LAST */
    for (argp = which + 2;
	 argp->ar_short_arg != ARGV_LAST
	 && (argp - 1)->ar_short_arg != ARGV_LAST;
	 argp += 2) {
      if ((argp - 1)->ar_short_arg != ARGV_OR
	  && (argp - 1)->ar_short_arg != ARGV_XOR)
	break;
      if (argp->ar_type & ARGV_FLAG_USED) {
	matchp = argp;
	break;
      }
    }
  }
  
  /* did we not find a problem? */
  if (matchp == NULL)
    return NOERROR;
  
  if (argv_error_stream == NULL)
    return ERROR;
  
  (void)fprintf(argv_error_stream,
		"%s: %s, specify only one of the following:\n",
		argv_program, USAGE_ERROR_NAME);
  
  /* little hack to print the one that matched and the one we were checking */
  for (;;) {
    if (matchp->ar_long_arg == NULL)
      (void)fprintf(argv_error_stream, "     %s%c\n",
		    SHORT_PREFIX, matchp->ar_short_arg);
    else
      (void)fprintf(argv_error_stream, "     %s%c (%s%s)\n",
		    SHORT_PREFIX, matchp->ar_short_arg,
		    LONG_PREFIX, matchp->ar_long_arg);
    
    if (matchp == which)
      break;
    matchp = which;
  }
  
  return ERROR;
}

/*
 * find all the XOR arguments and make sure each group has at least
 * one specified in it.  returns [NO]ERROR
 */
LOCAL	int	check_xor(argv_t * args)
{
  argv_t	*startp = NULL, *argp;
  
  /* run through the list of arguments */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++) {
    
    /* only check the XORs */
    if (argp->ar_short_arg != ARGV_XOR)
      continue;
    
    startp = argp;
    
    /*
     * NOTE: we are guaranteed that we are on a XOR so there is
     * something below and above...
     */
    if ((argp - 1)->ar_type & ARGV_FLAG_USED)
      startp = NULL;
    
    /* run through all XORs */
    for (;;) {
      argp++;
      if (argp->ar_type & ARGV_FLAG_USED)
	startp = NULL;
      if ((argp + 1)->ar_short_arg != ARGV_XOR)
	break;
      argp++;
    }
    
    /* were none of the xor's filled? */
    if (startp != NULL)
      break;
  }
  
  /* did we not find a problem? */
  if (startp == NULL)
    return NOERROR;
  
  /* argp points to the first XOR which failed */
  if (argv_error_stream == NULL)
    return ERROR;
  
  (void)fprintf(argv_error_stream, "%s: %s, must specify one of:\n",
		argv_program, USAGE_ERROR_NAME);
  
  for (argp = startp;; argp += 2) {
    /*
     * NOTE: we are guaranteed that we are on a XOR so there is
     * something below and above...
     */
    (void)fprintf(argv_error_stream, "     %s%c",
		  SHORT_PREFIX, (argp - 1)->ar_short_arg);
    if ((argp - 1)->ar_long_arg != NULL)
      (void)fprintf(argv_error_stream, " (%s%s)",
		    LONG_PREFIX, (argp - 1)->ar_long_arg);
    (void)fprintf(argv_error_stream, "\n");
    
    if (argp->ar_short_arg != ARGV_XOR)
      break;
  }
  
  return ERROR;
}

/*
 * check to see if we have any mandatory arguments left.  returns [NO]ERROR
 */
LOCAL	int	check_mand(argv_t * args)
{
  argv_t	*argp;
  int		slotc = 0;
  
  /* see if there are any mandatory args left */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++)
    if (argp->ar_short_arg == ARGV_MAND
	&& (! (argp->ar_type & ARGV_FLAG_USED))
	&& ! (argp->ar_type & ARGV_ARRAY))
      slotc++;
  
  if (slotc > 0) {
    if (argv_error_stream != NULL)
      (void)fprintf(argv_error_stream,
		    "%s: %s, %d more mandatory argument%s must be specified\n",
		    argv_program, USAGE_ERROR_NAME,
		    slotc, (slotc == 1 ? "" : "s"));
    return ERROR;
  }
  
  return NOERROR;
}

/*
 * check for any missing argument options.  returns [NO]ERROR
 */
LOCAL	int	check_opt(void)
{
  int	queuec;
  
  queuec = QUEUE_COUNT();
  if (queuec > 0) {
    if (argv_error_stream != NULL)
      (void)fprintf(argv_error_stream,
		    "%s: %s, %d more option-argument%s must be specified\n",
		    argv_program, USAGE_ERROR_NAME,
		    queuec, (queuec == 1 ? "" : "s"));
    return ERROR;
  }
  
  return NOERROR;
}

/**************************** argument processing ****************************/

/*
 * process an argument in MATCHP which looking at GRID. sets okayp to
 * FALSE if the argument was not okay.
 */
LOCAL	void	do_arg(argv_t * grid, argv_t * matchp, char * closep,
		       char * okayp)
{
  if (global_multi == GLOBAL_MULTI_REJECT) {
    /*
     * have we used this one before?
     * NOTE: should this be a warning or a non-error altogether?
     */
    if (matchp->ar_type & ARGV_FLAG_USED
	&& (! (matchp->ar_type & ARGV_ARRAY))
	&& ARGV_TYPE(matchp->ar_type) != ARGV_INCR) {
      if (argv_error_stream != NULL)
	(void)fprintf(argv_error_stream,
		      "%s: %s, you've already specified the '%c' argument\n",
		      argv_program, USAGE_ERROR_NAME,
		      matchp->ar_short_arg);
      *okayp = ARGV_FALSE;
    }
  }
  
  /* we used this argument */
  matchp->ar_type |= ARGV_FLAG_USED;
  
  /* check arguments that must be OR'd */
  if (check_or(grid, matchp) != NOERROR) {
    /*
     * don't return here else we might generate an XOR error
     * because the argument wasn't specified
     */
    *okayp = ARGV_FALSE;
  }
  
  /*
   * If we have a close argument, pass to translate.  If it is a
   * boolean or increment variable, then pass in a value of null
   * else queue it for needing a value argument.
   */
  if (global_close == GLOBAL_CLOSE_ENABLE && closep != NULL)
    translate_value(closep, matchp->ar_variable, matchp->ar_type);
  else if (! HAS_ARG(matchp->ar_type))
    translate_value(NULL, matchp->ar_variable, matchp->ar_type);
  else if (global_close == GLOBAL_CLOSE_ENABLE && closep != NULL)
    translate_value(closep, matchp->ar_variable, matchp->ar_type);
  else
    QUEUE_ENQUEUE(matchp);
}

/*
 * process a list of arguments ARGV and ARGV as it applies to ARGS.
 * on NUM_ARGS members.  OKAYP is a pointer to the boolean error
 * marker.
 */
LOCAL	void	do_list(argv_t * grid, int argc, char ** argv,
			char * okayp)
{
  argv_t	*gridp, *matchp;
  int		len, charc, unwantc = 0;
  char		last_arg = ARGV_FALSE;
  char		*closep = NULL, **argp;
  
  /* run throught rest of arguments */
  for (argp = argv; argp < argv + argc; argp++) {
    
    /* have we reached the LAST_ARG marker? */
    if (! last_arg && strcmp(LAST_ARG, *argp) == 0) {
      last_arg = ARGV_TRUE;
      continue;
    }
    
    /* check for close equals marker */
    if (! last_arg && global_close == GLOBAL_CLOSE_ENABLE) {
      closep = strchr(*argp, ARG_EQUALS);
      /* if we find the special char then punch the null and set pointer */
      if (closep != NULL) {
	*closep = NULLC;
	closep++;
      }
    }
    
    /* are we processing a long option? */
    if (! last_arg && strncmp(LONG_PREFIX, *argp, LONG_PREFIX_LENGTH) == 0) {
      
      /* get length of rest of argument */
      len = strlen(*argp) - LONG_PREFIX_LENGTH;
      
      /* we need more than the prefix */
      if (len <= 0) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, empty long-option prefix '%s'\n",
			argv_program, USAGE_ERROR_NAME, *argp);
	
	*okayp = ARGV_FALSE;
	continue;
      }
      
      matchp = NULL;
      
      /* run though long options looking for a match */
      for (gridp = grid; gridp->ar_short_arg != ARGV_LAST; gridp++) {
	if (gridp->ar_long_arg == NULL)
	  continue;
	
	if (strncmp(*argp + LONG_PREFIX_LENGTH,
		    gridp->ar_long_arg, len) == 0) {
	  if (matchp != NULL) {
	    if (argv_error_stream != NULL)
	      (void)fprintf(argv_error_stream,
			    "%s: %s, '%s' might be '%s' or '%s'\n",
			    argv_program, USAGE_ERROR_NAME, *argp,
			    gridp->ar_long_arg, matchp->ar_long_arg);
	    
	    *okayp = ARGV_FALSE;
	    break;
	  }
	  
	  /* record a possible match */
	  matchp = gridp;
	  
	  /* don't break, need to see if another one matches */
	}
      }
      
      /* if we found a match but quit then we must have found two matches */
      if (matchp != NULL && gridp->ar_short_arg != ARGV_LAST)
	continue;
      
      /* did we not find long-option match? */
      if (matchp == NULL) {
	
	/* check for special usage value */
	if (strncmp(USAGE_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0
	    || strncmp(HELP_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    do_usage(grid, global_usage);
	    (void)exit(0);
	  }
	  continue;
	}
	
	/* check for special short-usage value */
	if (strncmp(USAGE_SHORT_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    do_usage(grid, GLOBAL_USAGE_SHORT);
	    (void)exit(0);
	  }
	  continue;
	}
	
	/* check for special long-usage value */
	if (strncmp(USAGE_LONG_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    do_usage(grid, GLOBAL_USAGE_LONG);
	    (void)exit(0);
	  }
	  continue;
	}
	
	/* check for special long-usage value */
	if (strncmp(USAGE_ALL_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    do_usage(grid, GLOBAL_USAGE_ALL);
	    (void)exit(0);
	  }
	  continue;
	}
	
	/* check for special help value */
	if (strncmp(HELP_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    if (argv_error_stream != NULL) {
	      if (argv_help_string == NULL)
		(void)fprintf(argv_error_stream,
			      "%s: I'm sorry, no help is available.\n",
			      argv_program);
	      else
		(void)fprintf(argv_error_stream, "%s: %s\n",
			      argv_program, argv_help_string);
	    }
	    (void)exit(0);
	  }
	  continue;
	}
	
	/* check for special version value */
	if (strncmp(VERSION_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    if (argv_error_stream != NULL) {
	      if (argv_version_string == NULL)
		(void)fprintf(argv_error_stream,
			      "%s: no version information is available.\n",
			      argv_program);
	      else
		(void)fprintf(argv_error_stream, "%s: %s%s\n",
			      argv_program, VERSION_LABEL,
			      argv_version_string);
	    }
	    (void)exit(0);
	  }
	  continue;
	}
	
	/* check for display arguments value */
	if (strncmp(DISPLAY_ARG, *argp + LONG_PREFIX_LENGTH, len) == 0) {
	  if (argv_interactive) {
	    if (argv_error_stream != NULL)
	      display_variables(grid);
	    (void)exit(0);
	  }
	  continue;
	}
	
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, unknown long option '%s'.\n",
			argv_program, USAGE_ERROR_NAME, *argp);
	
	*okayp = ARGV_FALSE;
	continue;
      }
      
      do_arg(grid, matchp, closep, okayp);
      continue;
    }
    
    /* are we processing a short option? */
    if (! last_arg && strncmp(SHORT_PREFIX, *argp, SHORT_PREFIX_LENGTH) == 0) {
      
      /* get length of rest of argument */
      len = strlen(*argp) - SHORT_PREFIX_LENGTH;
      
      /* we need more than the prefix */
      if (len <= 0) {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: %s, empty short-option prefix '%s'\n",
			argv_program, USAGE_ERROR_NAME, *argp);
	*okayp = ARGV_FALSE;
	continue;
      }
      
      /* run through the chars in this option */
      for (charc = 0; charc < len; charc++) {
	
	/* run through the arg list looking for a match */
	for (matchp = grid; matchp->ar_short_arg != ARGV_LAST; matchp++)
	  if (matchp->ar_short_arg == (*argp)[SHORT_PREFIX_LENGTH + charc])
	    break;
	
	/* did we not find argument? */
	if (matchp->ar_short_arg == ARGV_LAST) {
	  
	  /* check for special usage value */
	  if ((*argp)[SHORT_PREFIX_LENGTH + charc] == USAGE_CHAR_ARG) {
	    if (argv_interactive) {
	      do_usage(grid, global_usage);
	      (void)exit(0);
	    }
	    continue;
	  }
	  
	  /* create an error string */
	  if (argv_error_stream != NULL)
	    (void)fprintf(argv_error_stream,
			  "%s: %s, unknown short option '%s%c'.\n",
			  argv_program, USAGE_ERROR_NAME, SHORT_PREFIX,
			  (*argp)[SHORT_PREFIX_LENGTH + charc]);
	  *okayp = ARGV_FALSE;
	  continue;
	}
	
	do_arg(grid, matchp, closep, okayp);
      }
      
      continue;
    }
    
    /* could this be a value? */
    if (grid->ar_short_arg != ARGV_LAST && QUEUE_COUNT() > 0) {
      QUEUE_DEQUEUE(matchp);
      translate_value(*argp, matchp->ar_variable, matchp->ar_type);
      continue;
    }
    
    /* process mandatory args if some left to process */
    for (gridp = grid; gridp->ar_short_arg != ARGV_LAST; gridp++)
      if (gridp->ar_short_arg == ARGV_MAND
	  && (! (gridp->ar_type & ARGV_FLAG_USED)))
	break;
    if  (gridp->ar_short_arg != ARGV_LAST) {
      /* absorb another mand. arg */
      translate_value(*argp, gridp->ar_variable, gridp->ar_type);
      if (! (gridp->ar_type & ARGV_ARRAY))
	gridp->ar_type |= ARGV_FLAG_USED;
      continue;
    }
    
    /* process maybe args if some left to process */
    for (gridp = grid; gridp->ar_short_arg != ARGV_LAST; gridp++)
      if (gridp->ar_short_arg == ARGV_MAYBE
	  && (! (gridp->ar_type & ARGV_FLAG_USED)))
	break;
    if  (gridp->ar_short_arg != ARGV_LAST) {
      /* absorb another maybe arg */
      translate_value(*argp, gridp->ar_variable, gridp->ar_type);
      if (! (gridp->ar_type & ARGV_ARRAY))
	gridp->ar_type |= ARGV_FLAG_USED;
      continue;
    }
  
    /* default is an error */
    unwantc++;
    *okayp = ARGV_FALSE;
  }
  
  if (unwantc > 0 && argv_error_stream != NULL)
    (void)fprintf(argv_error_stream,
		  "%s: %s, %d unwanted additional argument%s\n",
		  argv_program, USAGE_ERROR_NAME,
		  unwantc, (unwantc == 1 ? "" : "s"));
}

/****************************** env processing *******************************/

/*
 * handle the args from the ENV variable
 */
LOCAL	void	do_env_args(argv_t * args, char * okayp)
{
  int	envc, envn;
  char	**vectp, env_name[256], *environp;
  
  /* create the env variable */
  (void)sprintf(env_name, ENVIRON_FORMAT, argv_program);
  
  /* NOTE: by default the env name is all uppercase */
  for (environp = env_name; *environp != NULLC; environp++)
    if (islower(*environp))
      *environp = toupper(*environp);
  
  environp = getenv(env_name);
  if (environp == NULL)
    return;
  
  /* break the list into tokens and do the list */
  environp = string_copy(environp);
  vectp = vectorize(environp, " \t", &envn);
  if (vectp != NULL) {
    do_list(args, envn, vectp, okayp);
    
    /* free token list */
    for (envc = 0; envc < envn; envc++)
      free(vectp[envc]);
    free(vectp);
  }
  free(environp);
}

/*
 * process the global env variable
 */
LOCAL	void	process_env(void)
{
  static char	done = ARGV_FALSE;
  char		*environp, *envp, *arg;
  int		len;
  
  /* make sure we only do this once */
  if (done)
    return;
  else
    done = ARGV_TRUE;
  
  /* get the argv information */
  environp = getenv(GLOBAL_NAME);
  if (environp == NULL)
    return;
  
  /* save a copy of it */
  environp = string_copy(environp);
  arg = environp;
  
  for (;;) {
    envp = strtok(arg, " \t,:");
    if (envp == NULL)
      break;
    arg = NULL;
    
    len = strlen(GLOBAL_CLOSE);
    if (strncmp(GLOBAL_CLOSE, envp, len) == 0) {
      envp += len;
      if (strcmp(envp, "disable") == 0)
	global_close = GLOBAL_CLOSE_DISABLE;
      else if (strcmp(envp, "enable") == 0)
	global_close = GLOBAL_CLOSE_ENABLE;
      else {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: illegal env variable '%s' '%s' argument '%s'\n",
			__FILE__, GLOBAL_NAME, GLOBAL_CLOSE, envp);
      }
      continue;
    }
    
    len = strlen(GLOBAL_ENV);
    if (strncmp(GLOBAL_ENV, envp, len) == 0) {
      envp += len;
      if (strcmp(envp, "none") == 0)
	global_env = GLOBAL_ENV_NONE;
      else if (strcmp(envp, "before") == 0)
	global_env = GLOBAL_ENV_BEFORE;
      else if (strcmp(envp, "after") == 0)
	global_env = GLOBAL_ENV_AFTER;
      else {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: illegal env variable '%s' '%s' argument '%s'\n",
			__FILE__, GLOBAL_NAME, GLOBAL_ENV, envp);
      }
      continue;
    }
    
    len = strlen(GLOBAL_ERROR);
    if (strncmp(GLOBAL_ERROR, envp, len) == 0) {
      envp += len;
      if (strcmp(envp, "none") == 0)
	global_error = GLOBAL_ERROR_NONE;
      else if (strcmp(envp, "see") == 0)
	global_error = GLOBAL_ERROR_SEE;
      else if (strcmp(envp, "short") == 0)
	global_error = GLOBAL_ERROR_SHORT;
      else if (strcmp(envp, "shortrem") == 0)
	global_error = GLOBAL_ERROR_SHORTREM;
      else if (strcmp(envp, "long") == 0)
	global_error = GLOBAL_ERROR_LONG;
      else if (strcmp(envp, "all") == 0)
	global_error = GLOBAL_ERROR_ALL;
      else {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: illegal env variable '%s' '%s' argument '%s'\n",
			__FILE__, GLOBAL_NAME, GLOBAL_ERROR, envp);
      }
      continue;
    }
    
    len = strlen(GLOBAL_MULTI);
    if (strncmp(GLOBAL_MULTI, envp, len) == 0) {
      envp += len;
      if (strcmp(envp, "reject") == 0)
	global_multi = GLOBAL_MULTI_REJECT;
      else if (strcmp(envp, "accept") == 0)
	global_multi = GLOBAL_MULTI_ACCEPT;
      else {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: illegal env variable '%s' '%s' argument '%s'\n",
			__FILE__, GLOBAL_NAME, GLOBAL_MULTI, envp);
      }
      continue;
    }
    
    len = strlen(GLOBAL_USAGE);
    if (strncmp(GLOBAL_USAGE, envp, len) == 0) {
      envp += len;
      if (strcmp(envp, "short") == 0)
	global_usage = GLOBAL_USAGE_SHORT;
      else if (strcmp(envp, "shortrem") == 0)
	global_usage = GLOBAL_USAGE_SHORTREM;
      else if (strcmp(envp, "long") == 0)
	global_usage = GLOBAL_USAGE_LONG;
      else if (strcmp(envp, "all") == 0)
	global_usage = GLOBAL_USAGE_ALL;
      else {
	if (argv_error_stream != NULL)
	  (void)fprintf(argv_error_stream,
			"%s: illegal env variable '%s' '%s' argument '%s'\n",
			__FILE__, GLOBAL_NAME, GLOBAL_USAGE, envp);
      }
      continue;
    }
    
    if (argv_error_stream != NULL)
      (void)fprintf(argv_error_stream,
		    "%s: illegal env variable '%s' setting '%s'\n",
		    __FILE__, GLOBAL_NAME, envp);
  }
  
  free(environp);
}

/*
 * processes ARGS from ARGC and ARGV.
 * returns 0 if no error else -1.
 */
LOCAL	int	process_args(argv_t * args, int argc, char ** argv)
{
  int		num_args;
  char	*progp;
  char		okay = ARGV_TRUE;
  argv_t	*argp;
  
  process_env();
  
  if (args == NULL)
    args = empty;
  
  if (argc < 0) {
    if (argv_error_stream != NULL)
      (void)fprintf(argv_error_stream,
		    "%s: %s, argc argument to argv_process is %d\n",
		    __FILE__, INTERNAL_ERROR_NAME, argc);
    if (argv_interactive)
      (void)exit(EXIT_CODE);
    return ERROR;
  }
  
  if (argv == NULL) {
    if (argv_error_stream != NULL)
      (void)fprintf(argv_error_stream,
		    "%s: %s, argv argument to argv_process is NULL\n",
		    __FILE__, INTERNAL_ERROR_NAME);
    if (argv_interactive)
      (void)exit(EXIT_CODE);
    return ERROR;
  }
  
  /* set global variables */
  argv_argv = argv;
  argv_argc = argc;
  
  /* build the program name from the argv[0] path */
  {
    char	*tmpp;
    
    progp = *argv;
    for (tmpp = *argv; *tmpp != NULLC; tmpp++)
      if (*tmpp == '/')
	progp = tmpp + 1;
  }
  
  /* so we can step on the environmental space */
  (void)strncpy(argv_program, progp, PROGRAM_NAME);
  
  /* count the args */
  num_args = 0;
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++)
    num_args++;
  
  /* verify the argument array */
  if (preprocess_array(args, num_args) != NOERROR)
    return ERROR;
  
  /* allocate our value queue */
  if (num_args > 0)
    QUEUE_ALLOC(argv_t *, num_args);
  
  /* do the env args before? */
  if (global_env == GLOBAL_ENV_BEFORE)
    do_env_args(args, &okay);
  
  /* do the external args */
  do_list(args, argc - 1, argv + 1, &okay);
  
  /* do the env args after? */
  if (global_env == GLOBAL_ENV_AFTER)
    do_env_args(args, &okay);
  
  /* make sure the XOR and MAND args and argument-options are okay */
  if (check_mand(args) != NOERROR)
    okay = ARGV_FALSE;
  if (check_opt() != NOERROR)
    okay = ARGV_FALSE;
  if (check_xor(args) != NOERROR)
    okay = ARGV_FALSE;
  
  /* if we allocated the space then free it */
  if (num_args > 0)
    QUEUE_FREE();
  
  /* was there an error? */
#ifndef MPI_NEGF  /* some cluster versions of mpi append other command line options
		     that are not recognized. */
  if (! okay) {
    if (argv_error_stream != NULL)
      do_usage(args, global_error);
    if (argv_interactive)
      (void)exit(EXIT_CODE);
    return ERROR;
  }
#endif

  return NOERROR;
}

/***************************** exported routines *****************************/

/*
 * processes ARGC number of arguments from ARGV depending on argument
 * info array ARGS (if null then an empty array is used).  this
 * routine will not modify the argv array in any way.  NOTE: it will
 * modify the args array by setting various flags in the type field.
 * returns 0 if no error else -1.
 */
EXPORT	int	argv_process(argv_t * args, int argc, char ** argv)
{
  if (process_args(args, argc, argv) == NOERROR)
    return NOERROR;
  else
    return ERROR;
}

/*
 * processes arguments sent in via the STRING that a web-server might
 * send to program in ARG0.  returns 0 on noerror else -1.
 */
EXPORT	int	argv_web_process_string(argv_t * args, char * arg0,
					char * string)
{
  char *strp;
  char *copy, *copyp, **argv;
  int  argc, ret, alloced;
  
  /* copy incoming string so we can punch nulls */
  copy = (char*)malloc(strlen(string) + 1);
  
  /* create argv array */
  alloced = ARG_MALLOC_INCR;
  argv = (char **)malloc(sizeof(char *) * alloced);
  
  argc = 0;
  argv[argc++] = (char *)arg0;
  if (*string != NULLC) {
    if (argc >= alloced) {
      alloced += ARG_MALLOC_INCR;
      argv = (char **)realloc(argv, sizeof(char *) * alloced);
    }
    argv[argc++] = copy;
  }
  
  copyp = copy;
  for (strp = string;; strp++) {
    int		val;
    
    /* are we done? */
    if (*strp == NULLC) {
      *copyp = NULLC;
      break;
    }
    
    /* next arg marker */
    if (*strp == '&' || *strp == '=') {
      *copyp++ = NULLC;
      if (argc >= alloced) {
	alloced += ARG_MALLOC_INCR;
	argv = (char **)realloc(argv, sizeof(char *) * alloced);
      }
      argv[argc++] = copyp;
      continue;
    }
    
    /* a space */
    if (*strp == '+') {
      *copyp++ = ' ';
      continue;
    }
    
    /* no binary character, than it is normal */
    if (*strp != '%') {
      *copyp++ = *strp;
      continue;
    }      
    
    strp++;
    
    if (*strp >= 'a' && *strp <= 'f')
      val = 10 + *strp - 'a';
    else if (*strp >= 'A' && *strp <= 'F')
      val = 10 + *strp - 'A';
    else if (*strp >= '0' && *strp <= '9')
      val = *strp - '0';
    else
      continue;
    
    strp++;
    
    if (*strp >= 'a' && *strp <= 'f')
      val = val * 16 + (10 + *strp - 'a');
    else if (*strp >= 'A' && *strp <= 'F')
      val = val * 16 + (10 + *strp - 'A');
    else if (*strp >= '0' && *strp <= '9')
      val = val * 16 + (*strp - '0');
    else
      strp--;
    
    *copyp++ = (char)val;
  }
  
  ret = process_args(args, argc, argv);
  
  free(copy);
  free(argv);
  
  if (ret == NOERROR)
    return NOERROR;
  else
    return ERROR;
}

/*
 * processes arguments sent in via the QUERY_STRING environmental
 * variable that a web-server might send to program in ARG0.  returns
 * 0 on noerror else -1.
 */
EXPORT	int	argv_web_process(argv_t * args, char * arg0)
{
  char	*env, *work = NULL;
  int	ret, len;
  
  env = getenv("REQUEST_METHOD");
  if (env != NULL && strcmp(env, "POST") == 0) {
    env = getenv("CONTENT_LENGTH");
    if (env != NULL) {
      len = atoi(env);
      if (len > 0) {
	work = (char *) malloc(len + 1);
	/*(void)read(STDIN, work, len);*/
	work[len] = NULLC;
      }
    }
  }
  
  if (work == NULL) {
    env = getenv("QUERY_STRING");
    
    /* if it is not set or empty, then nothing to do */
    if (env == NULL || *env == NULLC) {
      work = (char *) malloc(1);
      work[0] = NULLC;
    }
    else
      work = string_copy(env);
  }
  
  ret = argv_web_process_string(args, arg0, work);
  free(work);
  
  if (ret == NOERROR)
    return NOERROR;
  else
    return ERROR;
}

/*
 * print the standard usage messages for argument array ARGS (if null
 * then an empty array is used).  WHICH chooses between long or short
 * messages (see argv.h).
 * NOTE: if this is called before argv_process then the program name
 * may be messed up.
 */
EXPORT	void	argv_usage(argv_t * args, int which)
{
  process_env();
  
  if (args == NULL)
    args = empty;
  
  if (which == ARGV_USAGE_SHORT)
    usage_short(args, 0);
  else if (which == ARGV_USAGE_LONG)
    usage_long(args);
  else
    /* default/env settings */
    do_usage(args, global_usage);
}

/*
 * see if ARG argument was used in a previous call to argv_process on
 * ARGS, returns 1 if yes else 0
 */
EXPORT	int	argv_was_used(argv_t * args, char arg)
{
  argv_t	*argp;
  
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++)
    if (argp->ar_short_arg == arg) {
      if (argp->ar_type & ARGV_FLAG_USED)
	return 1;
      else
	return 0;
    }
  
  return 0;
}

/*
 * frees up any allocations in ARGS that may have been done by
 * argv_process.  This should be done at the end of the program or
 * after all the arguments have been referenced.
 */
EXPORT	void	argv_cleanup(argv_t * args)
{
  argv_t	*argp;
  int		entryc;
  
  if (args == NULL)
    return;
  
  /* run through the argument structure */
  for (argp = args; argp->ar_short_arg != ARGV_LAST; argp++) {
    /* handle any arrays */
    if (argp->ar_type & ARGV_ARRAY) {
      argv_array_t	*arrp = (argv_array_t *)argp->ar_variable;
      
      /* free any entries */
      if (arrp->aa_entryn > 0) {
	if (ARGV_TYPE(argp->ar_type) == ARGV_CHARP) {
	  for (entryc = 0; entryc < arrp->aa_entryn; entryc++)
	    free(ARGV_ARRAY_ENTRY(*arrp, char *, entryc));
	}
	free(arrp->aa_entries);
      }
      arrp->aa_entries = NULL;
      arrp->aa_entryn = 0;
      continue;
    }
    
    /* handle individual charps */
    if (argp->ar_type & ARGV_FLAG_USED
	&& ARGV_TYPE(argp->ar_type) == ARGV_CHARP) {
      free(*(char **)argp->ar_variable);
      continue;
    }
  }
}

/*
 * copy all the args (after the 0th), one after the other, into BUF of
 * MAX_SIZE.  NOTE: you can get the 0th argument from argv_argv[0].
 */
EXPORT	void	argv_copy_args(char * buf, int max_size)
{
  char	**argvp, *bufp = buf, *argp;
  int	argc, sizec = max_size;
  
  process_env();
  
  if (max_size == 0)
    return;
  
  *bufp = NULLC;
  
  if (argv_argv == NULL || max_size == 1)
    return;
  
  for (argvp = argv_argv + 1, argc = 1; argc < argv_argc; argvp++, argc++) {
    
    if (sizec < 2)
      break;
    
    if (argvp > argv_argv + 1) {
      *bufp++ = ' ';
      sizec--;
    }
    
    for (argp = *argvp; *argp != NULLC && sizec >= 2; sizec--)
      *bufp++ = *argp++;
  }
  
  *bufp = NULLC;
}
