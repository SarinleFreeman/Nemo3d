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
$Header: /repo/nemo3d/src/util/util_stringops.c,v 1.5 2004/10/27 19:43:10 gekco
Exp $
*****************************************************************************/

#include "util_stringops.h"

/*
 * Allocates new memory to copy the character string into.
 */
char *n3d_strdup_n(const char *str) {

  if (str) {
    return strcpy((char *)str_malloc(strlen(str) + 1), str);
  } else {
    return (char *)NULL;
  }
}

/*
 * Allocates new memory space, removes trailing white spaces from the string,
 * and then copies the string into the memory space.
 */
char *n3d_copy_str(const char *str) {

  size_t length;
  char *s, *s2;

  if (str != NULL) {
    length = strlen(str);
    if (length <= 0) {
      return NULL;
    }
    s = (char *)str_calloc((length + 1), sizeof(char));
    strcpy(s, str);

    /* Remove trailing whitespaces */
    for (s2 = s + length - 1; (isspace(*s2) && s2 != s); s2--) {
      *s2 = '\0';
    }
    return s;

  } else {
    return NULL;
  }
}

/*
 * Contactenates two or more character strings. It allocates space for the
 * resulting string and NULL terminates it. Since it uses str_calloc() to
 * allocate space, the calling routing should call str_free() to free up the
 * string.
 */
char *n3d_strcat(char *firstToken, ...) {
  int charLen;
  char *fullStr, *token;
  va_list v;

  if (firstToken) {
    charLen = strlen(firstToken);

    /* Get total length of firstToken and the rest of the tokens */
    va_start(v, firstToken);
    while ((token = va_arg(v, char *))) {
      charLen += strlen(token);
    }

    /* Allocate new memory to hold all the tokens and set the initial value to
     * be firstToken */
    fullStr = (char *)str_calloc(charLen + 1, sizeof(char));
    strcat(fullStr, firstToken);

    /* Concatenate the rest of the tokens into fullStr */
    va_start(v, firstToken);
    while ((token = va_arg(v, char *))) {
      strcat(fullStr, token);
    }

    va_end(v);

    return fullStr;
  }

  else {
    return NULL;
  }
}

/*
 * Same as strcmp() except that a check is done for the two strings before it
 * calls strcmp(). If both strings are null, a 0 is returned.
 */
int n3d_strcmp_null_check(const char *str1, const char *str2) {

  if (!str1 && !str2) {
    return 0;
  } else if (!str1 && str2) {
    return 1;
  } else if (!str2) {
    return -1;
  }

  return strcmp(str1, str2);
}

/* Function retunrs allocated memory that needs to be free'ed with str_free */
char *n3d_FileNameGet(char *filename) {
  char *periodPosition;
  char *filenameDup;

  filenameDup = n3d_strdup_n(filename);
  if ((periodPosition = strrchr(filenameDup, '.'))) {
    *periodPosition = (char)NULL;
  }
  return filenameDup;
}

/*
 * Sets the file type extension for the nemo3d output.
 */
/* the allocated memory needs to be free'ed by str_free */
void n3d_FileTypeSet(char **ptrFile, char *fileType, int overWrite) {

  char *filename;

  if (!*ptrFile || !fileType)
    return;

  if (strrchr(*ptrFile, '.') && !overWrite)
    return;

  filename = n3d_FileNameGet(*ptrFile);
  str_free(*ptrFile);
  *ptrFile = n3d_strcat(filename, ".", fileType, NULL);
  str_free(filename);
}

/*
 * Returns the number of string names
 */
int n3d_StringListNum(char **strList) {
  int numLabels;
  char **str;

  if (strList) {
    /* Get number of labels */
    for (numLabels = 0, str = strList; (str && str[numLabels]); numLabels++)
      ;
    return numLabels;
  } else {
    return 0;
  }
}

/*
 * Appends a string to an existing list of strings.
 * If no list exists, it will create one.
 */
char **n3d_StringListAdd(char ***str_list_ptr, const char *new_string) {

  char **str_ptr, **str_new;
  char **str_list = *str_list_ptr;
  int numStrings, i;

  if (str_list) {
    /* Get number of strings */
    for (str_ptr = str_list, numStrings = 0; *str_ptr; str_ptr++, numStrings++)
      ;

    str_new = (char **)str_calloc(numStrings + 2, sizeof(char *));

    for (i = 0; i < numStrings; i++)
      str_new[i] = str_list[i];

    str_new[numStrings] =
        (new_string) ? n3d_strdup_n(new_string) : n3d_strdup_n(" ");
    str_free(str_list);
  } else {
    str_new = (char **)str_calloc(2, sizeof(char *));
    str_new[0] = n3d_strdup_n(new_string);
  }

  *str_list_ptr = str_new;
  return str_new;
}
