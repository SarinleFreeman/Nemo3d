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
$Header: /repo/nemo3d/src/io/fileformat_old.c,v 1.1 2007/03/25 16:24:11 hoonryu
Exp $
*****************************************************************************/

#ifdef MACOSX
#include <sys/types.h>
#endif
#include <netinet/in.h>

#include "fileformat.h"

#ifdef AIX_SPECIFIC
#define bzero(a, b) memset(a, 0, b)
#endif

void parseRecordFormat(const char *label, int *Nint, int *Nreal, int *Nfloat,
                       int *Ncol_var, int *sz_var, int *Nvar, char *fieldType) {
  const char strR[] = "real ";
  const char strF[] = "float ";
  const char strI[] = "int ";
  int markR = 0, markI = 0, markF = 0;

  *Nint = *Nreal = *Nvar = *Nfloat = 0;

  for (const char *c = label; (*c) != '\0'; c++) {
    if (*c == '{') {
      Ncol_var[*Nvar] = 1;
      sz_var[*Nvar] = 0;
      (*Nvar)++;
    } else if (*c == ',') {
      Ncol_var[*Nvar - 1]++;
    }
    /* check for match with real */
    if (*c == ' ' && markR == 4) {
      markR = 0;
      fieldType[*Nint + *Nreal + *Nfloat] = 'R';
      sz_var[*Nvar - 1] += sizeof(double);
      (*Nreal)++;
    } else if (*c == strR[markR]) {
      markR++;
    } else
      markR = 0;
    /* check for match with int */
    if (*c == ' ' && markI == 3) {
      markI = 0;
      fieldType[*Nint + *Nreal + *Nfloat] = 'Z';
      sz_var[*Nvar - 1] += sizeof(int);
      (*Nint)++;
    } else if (*c == strI[markI]) {
      markI++;
    } else
      markI = 0;
    /* check for match with float */
    if (*c == ' ' && markF == 5) {
      markF = 0;
      fieldType[*Nint + *Nreal + *Nfloat] = 'S';
      sz_var[*Nvar - 1] += sizeof(nml_float);
      (*Nfloat)++;
    } else if (*c == strF[markF]) {
      markF++;
    } else
      markF = 0;
  }

  if (*Nvar > MAX_VARS) {
    die("Must specify less than %d variable arguments to 'writeSiteInfo'\n",
        MAX_VARS);
  }
}

FILE *readHeader(const char *filename, const char *recFmt, int *Ndim, int *dim,
                 char *dim_lbl[], int *Nuser, char *usr_comment[]) {
  const char strB[] = "<HDR>";
  const char strE[] = "<\\HDR>";
  const char strRD[] = "record dimension";
  const char strUD[] = "user data";
  char c;
  int markB = 0, markE = 0, markRD = 0, markUD = 0;

  FILE *fp = fopen(filename, "rb");
  if (!fp)
    die("Could not open file %s for reading\n", filename);

  *Ndim = *Nuser = 0;

  /* first get to beginning of header info */
  do {
    c = getc(fp);

    if (c == '>' && markB == 4) {
      break;
    } else if (c == strB[markB]) {
      markB++;
    } else
      markB = 0;
  } while (c != EOF);

  /* now look for '=' and copy what follows into recFmt */
  do {
    c = getc(fp);
  } while (c != EOF && c != '=');
  fscanf(fp, "%[^\n]s", recFmt);

  do {
    c = getc(fp);
    /* match record dimension */
    if (c == 'n' && markRD == 15) {
      do {
        c = getc(fp);
      } while (c != EOF && c != '=');
      fscanf(fp, "%d %[^\n]s", dim + (*Ndim), dim_lbl[*Ndim]);
      (*Ndim)++;
    } else if (c == strRD[markRD]) {
      markRD++;
    } else
      markRD = 0;
    /* match user data */
    if (c == 'a' && markUD == 8) {
      do {
        c = getc(fp);
      } while (c != EOF && c != '=');
      fscanf(fp, "%[^\n]s", usr_comment[*Nuser]);
      (*Nuser)++;
    } else if (c == strUD[markUD]) {
      markUD++;
    } else
      markUD = 0;
    /* match end tag */
    if (c == '>' && markE == 5) {
      break;
    } else if (c == strE[markE]) {
      markE++;
    } else
      markE = 0;

  } while (c != EOF);

  dim_lbl[*Ndim] = NULL;
  dim[*Ndim] = 0;
  usr_comment[*Nuser] = NULL;

  return fp;
}

void writeHeader(const char *filename, const char *recFmt, int *dim,
                 char *dim_lbl[], char *usr_comment[]) {
  char dum[30], dum2[30];

  FILE *fp = (!strcmp(filename, "STDOUT") ? stdout : fopen(filename, "w"));
  fprintf(fp, "<HDR>\n");
  fprintf(fp, "%-20s = %s\n", "record format", recFmt);
  for (int i = 1; *dim > 0; dim++, dim_lbl++, i++) {
    sprintf(dum2, "%d", i);
    strcat(strcpy(dum, "record dimension "), dum2);
    fprintf(fp, "%-20s = %-10d %s", dum, *dim, *dim_lbl);
    fprintf(fp, "\n");
  }
  for (int i = 1; *usr_comment; usr_comment++, i++) {
    sprintf(dum2, "%d", i);
    strcat(strcpy(dum, "user data "), dum2);
    fprintf(fp, "%-20s = %s", dum, *usr_comment);
    fprintf(fp, "\n");
  }
  fprintf(fp, "<\\HDR>");
  fclose(fp);
}

/* writes header to stdout */
void displayHeader(const char *filename) {
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar, Nfloat;
  char fieldType[21];

  bzero(fieldType, sizeof(char) * 21);
  for (int i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  FILE *fp;
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);
  writeHeader("STDOUT", label, dim, dim_lbl, usr_comment);
  printf("\n\n");
  fclose(fp);

  for (int i = 0; i < MAX_VARS; i++) {
    delete dim_lbl[i];
    delete usr_comment[i];
  }
}

void fmtdat(int argc, char *argv[]) {
  /* multi-line string literals are deprecated */
  string info = "USAGE:  fmtdat [-a?  | -b | -N | -H | -Wf2Pdb | -?] file";
  info.append("        -h:   display header");
  info.append(
      "        -b:   convert from ascii format to double precision binary");
  info.append("        -a?:  convert from ? precision binary format to ascii");
  info.append("        -N?:  convert binary data to network byte ordering");
  info.append("        -H?:  convert binary data to hardware byte ordering");
  info.append(
      "        -dx:  convert binary data to dx data for 3D visualization");
  info.append("        -Wf2Pdb wf_file pdb_file:  merge binary wavefunction "
              "file with pdb file and output a pdb format file");
  info.append("        -?:   convert binary data to ? precision\n");
  info.append("        Note:  For single precision, ?='1' or ''");
  info.append("        Note:  For double precision, ?='2'");

  if (argc < 3) {
    cout << info << endl;
    return;
  }

  char *filename = argv[argc - 1];

  for (int j = 1; j < argc; j++) {
    if (!strcmp(argv[j], "-h")) {
      displayHeader(filename);
    }
    if (!strcmp(argv[j], "-b")) {
      convertAsciiToBin(filename);
    }
#define ConvertWfToPdb
#ifdef ConvertWfToPdb
    /*add a function convertWfToPdb to convert Wf format to Pdb format*/
    else if (!strcmp(argv[j], "-Wf2Pdb")) {
      char *filename2 = argv[argc - 2];
      convertWfToPdb(filename2, filename);
    }
#endif
    else if (!strcmp(argv[j], "-a1")) {
      convertBinToAscii(filename, 1, 1);
    } else if (!strcmp(argv[j], "-a2")) {
      if (strstr(filename, "_dx_") != NULL) {
        printf("Cannot use -a2 for dx file. Use -dx option.!\n");
        exit(1);
      }
      convertBinToAscii(filename, 2, 1);
    } else if (!strcmp(argv[j], "-N1")) {
      reformatBinaryData(filename, "N", 1);
    } else if (!strcmp(argv[j], "-N2")) {
      reformatBinaryData(filename, "N", 2);
    } else if (!strcmp(argv[j], "-H1")) {
      reformatBinaryData(filename, "H", 1);
    } else if (!strcmp(argv[j], "-H2")) {
      reformatBinaryData(filename, "H", 2);
    } else if (!strcmp(argv[j], "-1")) {
      reformatBinaryData(filename, "1", 2);
    } else if (!strcmp(argv[j], "-2")) {
      reformatBinaryData(filename, "2", 2);
    }
    // added by hoon
    else if (!strcmp(argv[j], "-dx")) {

      if (strstr(filename, "_dx_") == NULL) {
        printf("Use -dx option for dx file only!\n");
        exit(1);
      }
      convertBinToAscii(filename, 2, 1);
      writeDXfile(filename);
    }
  }
}

void writeDXfile(char *filename) {

  //  printf("\n%s\n", filename);

  FILE *fp;

  char c, tension[5], shape[7], ev[10], lattice[10], xmax[10], ymax[10],
      zmax[10];
  char evalue[20], xbuffer[19], ybuffer[19], zbuffer[19], databuffer[20],
      newfilename[100];
  int count = 0;
  int RNum = 0;
  int evN, xi, yi, zi, index, cellxmax, cellymax, cellzmax, index1d;
  double value, eLvalue, lattice_x, lattice_y, lattice_z;
  long initpos;
  long datainitpos;

  strcpy(xbuffer, "");
  strcpy(ybuffer, "");
  strcpy(zbuffer, "");
  strcpy(databuffer, "");
  strcpy(ev, "");
  strcpy(lattice, "");
  strcpy(xmax, "");
  strcpy(ymax, "");
  strcpy(zmax, "");
  strcpy(evalue, "");
  strcpy(newfilename, "");
  strcpy(tension, "");
  strcpy(shape, "");

  evN = 0;
  xi = 0;
  yi = 0;
  zi = 0;
  index = 0;
  cellxmax = 0;
  cellymax = 0;
  cellzmax = 0;
  index1d = 0;
  value = 0;
  eLvalue = 0;
  lattice_x = 0;
  lattice_y = 0;
  lattice_z = 0;
  initpos = 0;
  datainitpos = 0;

  fp = NULL;

  fp = fopen(filename, "rt");

  if (fp == NULL) {
    printf("\n Ascii File Open Error!\n");
    exit(1);
  }

  initpos = ftell(fp);

  // find # of row
  while (1) {
    c = fgetc(fp);
    if (c == '\n')
      RNum++;
    else if (c == EOF)
      break;
  }

  RNum = RNum - 9;
  // printf("\n%d", RNum);
  fseek(fp, initpos, SEEK_SET);

  // get necessary info.
  // printf("\n%s\n", "point1");
  while (1) {
    c = fgetc(fp);
    if (c == '\n')
      count++;
    if (count == 5 && c == 'm') {
      count = 0;
      break;
    }
  }

  fgetc(fp);
  fgetc(fp);
  fread(xmax, 9, 1, fp);
  xmax[9] = '\0';
  // printf("\n%s\n", xmax);
  fgetc(fp);
  fgetc(fp);
  fread(ymax, 9, 1, fp);
  ymax[9] = '\0';
  fgetc(fp);
  fgetc(fp);
  fread(zmax, 9, 1, fp);
  zmax[9] = '\0';

  // printf("\n%s\n", "point2-3");

  fseek(fp, initpos, SEEK_SET);
  count = 0;

  // printf("\n%s\n", "point2");

  if (strstr(filename, ".nd_dx_e") != NULL) {
    while (1) {
      c = fgetc(fp);
      if (c == '\n')
        count++;
      if (count == 6 && c == '=') {
        count = 0;
        break;
      }
    }

    fgetc(fp);
    fgetc(fp);
    fread(tension, 4, 1, fp);
    tension[4] = '\0';
  }

  else if (strstr(filename, ".nd_dx_shape") != NULL) {
    while (1) {
      c = fgetc(fp);
      if (c == '\n')
        count++;
      if (count == 6 && c == '=') {
        count = 0;
        break;
      }
    }

    fgetc(fp);
    fgetc(fp);
    fread(shape, 6, 1, fp);
    shape[6] = '\0';
  }

  else {
    while (1) {
      c = fgetc(fp);
      if (c == '\n')
        count++;
      if (count == 6 && c == 'v') {
        count = 0;
        break;
      }
    }
    fgetc(fp);
    fgetc(fp);
    fread(ev, 9, 1, fp);
    ev[9] = '\0';
    fgetc(fp);
    fgetc(fp);
    fread(evalue, 19, 1, fp);
    evalue[19] = '\0';
  }

  fseek(fp, initpos, SEEK_SET);
  // printf("\n%s\n", "point3");

  while (1) {
    c = fgetc(fp);
    if (c == '\n')
      count++;
    if (count == 7 && c == 'T') {
      count = 0;
      break;
    }
  }

  fgetc(fp);
  fgetc(fp);
  fread(lattice, 9, 1, fp);
  lattice[9] = '\0';
  lattice_x = atof(lattice);
  fgetc(fp);
  fgetc(fp);
  fread(lattice, 9, 1, fp);
  lattice[9] = '\0';
  lattice_y = atof(lattice);
  fgetc(fp);
  fgetc(fp);
  fread(lattice, 9, 1, fp);
  lattice[9] = '\0';
  lattice_z = atof(lattice);

  cellxmax = atof(xmax);
  cellymax = atof(ymax);
  cellzmax = atof(zmax);

  // printf("\n%s\n", "point4");

  if (strstr(filename, ".nd_dx_e") == NULL) {
    evN = atof(ev);
    eLvalue = atof(evalue);
  }

  // printf("\n%s\n", "point5");
  double Lvalue[cellxmax * cellymax * cellzmax];
  // printf("\n%s\n", "point6");

  fseek(fp, initpos, SEEK_SET);
  count = 0;
  // printf("\n%s\n", "point5");
  // skip first 9 lines
  //
  //
  while (1) {
    if (fgetc(fp) == '\n')
      count++;
    if (count == 9) {
      count = 0;
      break;
    }
  }
  // printf("\n%s\n", "point6");
  //  move to the initial position
  datainitpos = ftell(fp);
  fgetc(fp);
  fgetc(fp);
  fgetc(fp);
  for (index = 1; index <= RNum; index++) {
    read_oneLine(fp, xbuffer, ybuffer, zbuffer, databuffer, &xi, &yi, &zi,
                 &value);
    //  printf("\n%s %s %s %s %d %d %d %f\n", xbuffer, ybuffer, zbuffer,
    //  databuffer, xi, yi, zi, value);

    if (index == 1) {
      index1d = (xi - 1) * cellymax * cellzmax + (yi - 1) * cellzmax + (zi - 1);
      for (count = 0; count < index1d; count++) {
        Lvalue[count] = 0.0;
      }
    }

    else {
      for (count = index1d + 1; count < (xi - 1) * cellymax * cellzmax +
                                            (yi - 1) * cellzmax + (zi - 1);
           count++) {
        Lvalue[count] = 0.0;
      }
    }

    index1d = (xi - 1) * cellymax * cellzmax + (yi - 1) * cellzmax + (zi - 1);
    Lvalue[count] = value;

    if (index == RNum) {
      for (count = index1d + 1; count <= cellxmax * cellymax * cellzmax - 1;
           count++) {
        Lvalue[count] = 0.0;
      }
    }
  }

  fclose(fp); // we don't need temporarily-generated file any more.
  fp = NULL;

  // printf("\n%s\n", "point7");
  fp = fopen(filename, "wt");

  if (strstr(filename, ".nd_dx_shape") != NULL) {
    writeDXHeader_For_Tension_Profile(fp, cellxmax, cellymax, cellzmax,
                                      lattice_x, lattice_y, lattice_z, shape);
  }

  else if (strstr(filename, ".nd_dx_e") != NULL) {
    writeDXHeader_For_Tension_Profile(fp, cellxmax, cellymax, cellzmax,
                                      lattice_x, lattice_y, lattice_z, tension);

  } else
    writeDXHeader_For_WF_Profile(fp, cellxmax, cellymax, cellzmax, evN, eLvalue,
                                 lattice_x, lattice_y, lattice_z);

  for (count = 0; count < cellxmax * cellymax * cellzmax; count++) {
    fprintf(fp, "%g\n", Lvalue[count]);
  }
  writeDXTail(fp);

  // printf("\n%d %d %d \n", cellxmax, cellymax, cellzmax);

  fclose(fp);
}

void read_oneLine(FILE *fp, char *xbuffer, char *ybuffer, char *zbuffer,
                  char *databuffer, int *xi, int *yi, int *zi, double *value) {

  fread(xbuffer, 18, 1, fp);
  xbuffer[18] = '\0';
  fgetc(fp);
  fgetc(fp);
  fgetc(fp);
  fread(ybuffer, 18, 1, fp);
  ybuffer[18] = '\0';
  fgetc(fp);
  fgetc(fp);
  fgetc(fp);
  fread(zbuffer, 18, 1, fp);
  zbuffer[18] = '\0';
  fgetc(fp);
  fgetc(fp);
  fread(databuffer, 19, 1, fp);
  xbuffer[19] = '\0';
  fgetc(fp);
  fgetc(fp);
  fgetc(fp);
  fgetc(fp);
  *xi = atof(xbuffer);
  *yi = atof(ybuffer);
  *zi = atof(zbuffer);
  *value = atof(databuffer);
}

void writeDXHeader_For_WF_Profile(FILE *fp, int cellxmax, int cellymax,
                                  int cellzmax, int evN, double eLvalue,
                                  double lattice_x, double lattice_y,
                                  double lattice_z) {

  fprintf(fp, "#### Generated directly from NEMO3D\n");
  fprintf(fp, "#### eigenvalue # %d\n", evN);
  fprintf(fp, "#### Eigenvalue %g\n", eLvalue);
  fprintf(fp, "object 1 class gridpositions counts %d %d %d\n", cellxmax,
          cellymax, cellzmax);
  fprintf(fp, "origin 0 0 0\n");
  fprintf(fp, "delta  %g 0 0\n", lattice_x);
  fprintf(fp, "delta  0 %g 0\n", lattice_y);
  fprintf(fp, "delta  0 0 %g\n", lattice_z);
  fprintf(fp, "object 2 class gridconnections counts %d %d %d\n", cellxmax,
          cellymax, cellzmax);
  fprintf(fp, "attribute \"element type\" string \"cubes\"\n");
  fprintf(fp, "attribute \"ref\" string \"positions\"\n");
  fprintf(fp, "object 3 class array type float rank 0 items %d data follows\n",
          cellxmax * cellymax * cellzmax);
}

void writeDXHeader_For_Tension_Profile(FILE *fp, int cellxmax, int cellymax,
                                       int cellzmax, double lattice_x,
                                       double lattice_y, double lattice_z,
                                       char *tension) {

  fprintf(fp, "#### Generated directly from NEMO3D\n");
  fprintf(fp, "#### Tension: %s\n", tension);
  fprintf(fp, "object 1 class gridpositions counts %d %d %d\n", cellxmax,
          cellymax, cellzmax);
  fprintf(fp, "origin 0 0 0\n");
  fprintf(fp, "delta  %g 0 0\n", lattice_x);
  fprintf(fp, "delta  0 %g 0\n", lattice_y);
  fprintf(fp, "delta  0 0 %g\n", lattice_z);
  fprintf(fp, "object 2 class gridconnections counts %d %d %d\n", cellxmax,
          cellymax, cellzmax);
  fprintf(fp, "attribute \"element type\" string \"cubes\"\n");
  fprintf(fp, "attribute \"ref\" string \"positions\"\n");
  fprintf(fp, "object 3 class array type float rank 0 items %d data follows\n",
          cellxmax * cellymax * cellzmax);
}

void writeDXTail(FILE *fp) {

  fprintf(fp, "attribute \"dep\" string \"positions\"\n");
  fprintf(fp, "object \"regular positions regular connections\" class field\n");
  fprintf(fp, "component \"positions\" value 1\n");
  fprintf(fp, "component \"connections\" value 2\n");
  fprintf(fp, "component \"data\" value 3\n");
  fprintf(fp, "end\n");
}

void reformatBinaryData(char *filename, const char *flags, int prec) {
  int sz_chunk = 1 << 20;
  int szReal = (prec == 1 ? sizeof(float) : sizeof(double));
  int szReal_long = szReal / sizeof(long);

  short _convertHW = 0, _convertNET = 0, _convertSNGL = 0, _convertDBL = 0;

  FILE *fp, *fp_new;
  int i, j, k, Nrow = 1;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar, Nfloat;
  int N_chunk, sz_rem;
  char fieldType[21], filename_new[100];
  int sz_datstr;

  for (const char *ptr = flags; *ptr; ptr++) {
    switch (*ptr) {
    case 'N':
      strcat(strcpy(filename_new, filename), "_net");
      _convertNET = 1;
      break;
    case 'H':
      strcat(strcpy(filename_new, filename), "_hw");
      _convertHW = 1;
      break;
    case '1':
      strcat(strcpy(filename_new, filename), "_sngl");
      _convertSNGL = 1;
      break;
    case '2':
      strcat(strcpy(filename_new, filename), "_dbl");
      _convertDBL = 1;
      break;
    default:
      die("Unrecognized flag.");
    }
  }

  for (i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  if (!fp)
    die("Could not find file %s\n", filename);

  bzero(fieldType, sizeof(char) * 21);
  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);

  if (Nint)
    die("Can't handle integers yet");
  if (Nfloat)
    die("Can't handle floats yet");

  writeHeader(filename_new, label, dim, dim_lbl, usr_comment);

  for (i = 0; i < Ndim; i++) {
    Nrow *= dim[i];
  }

  for (i = 0; i < MAX_VARS; i++) {
    delete dim_lbl[i];
    delete usr_comment[i];
  }

  fp_new = fopen(filename_new, "a");

  /* break up reads/writes into smaller chunks */
  N_chunk = Nrow / (sz_chunk);
  sz_rem = Nrow % (sz_chunk);
  sz_datstr = Nreal * szReal + Nint * sizeof(int) + Nfloat * sizeof(nml_float);

  if ((_convertNET || _convertHW) && !_convertSNGL) {
    unsigned long *dat1_ptr, *dat1, *dat2, *dat2_ptr;
    dat1 = (unsigned long int *)malloc(sz_datstr * sz_chunk);
    dat2 = (unsigned long *)malloc(sizeof(unsigned long) * szReal_long *
                                   sz_chunk * Nreal);

    /* do full chunks */
    for (k = 0; k < N_chunk; k++) {
      dat2_ptr = dat2;

      fread(dat1, sz_datstr, sz_chunk, fp);
      dat1_ptr = dat1;

      if (_convertNET) {
        for (i = 0; i < sz_chunk * Nreal; i++) {
          for (j = 0; j < szReal_long; j++) {
            dat2_ptr[szReal_long - 1 - j] = htonl(*dat1_ptr);
            dat1_ptr++;
          }
          dat2_ptr += szReal_long;
        }
      } else if (_convertHW) {
        for (i = 0; i < sz_chunk * Nreal; i++) {
          for (j = 0; j < szReal_long; j++) {
            dat2_ptr[szReal_long - 1 - j] = ntohl(*dat1_ptr);
            dat1_ptr++;
          }
          dat2_ptr += szReal_long;
        }
      }
      fwrite(dat2, sizeof(long), szReal_long * sz_chunk * Nreal, fp_new);
    }

    /* do what's left over */
    fread(dat1, sz_datstr, sz_rem, fp);
    dat1_ptr = dat1;
    dat2_ptr = dat2;

    if (_convertNET) {
      for (i = 0; i < sz_rem * Nreal; i++) {
        for (j = 0; j < szReal_long; j++) {
          dat2_ptr[szReal_long - 1 - j] = htonl(*dat1_ptr);
          dat1_ptr++;
        }
        dat2_ptr += szReal_long;
      }
    } else if (_convertHW) {
      for (i = 0; i < sz_rem * Nreal; i++) {
        for (j = 0; j < szReal_long; j++) {
          dat2_ptr[szReal_long - 1 - j] = ntohl(*dat1_ptr);
          dat1_ptr++;
        }
        dat2_ptr += szReal_long;
      }
    }
    fwrite(dat2, sizeof(long), szReal_long * sz_rem * Nreal, fp_new);

    free(dat1);
    free(dat2);
  } else if (_convertSNGL) {
    double *dat1_ptr, *dat1;
    float *dat2, *dat2_ptr;
    dat1 = new double[sz_datstr / sizeof(double) * sz_chunk];
    dat2 = new float[sz_datstr / sizeof(double) * sz_chunk];

    /* do full chunks */
    for (k = 0; k < N_chunk; k++) {
      dat2_ptr = dat2;

      fread(dat1, sz_datstr, sz_chunk, fp);
      dat1_ptr = dat1;

      for (i = 0; i < sz_chunk * Nreal; i++) {
        *dat2_ptr = (float)(*dat1_ptr);
        dat1_ptr++;
        dat2_ptr++;
      }
      fwrite(dat2, sizeof(float), sz_chunk * Nreal, fp_new);
    }

    /* do what's left over */
    fread(dat1, sz_datstr, sz_rem, fp);
    dat1_ptr = dat1;
    dat2_ptr = dat2;

    for (i = 0; i < sz_rem * Nreal; i++) {
      *dat2_ptr = (float)(*dat1_ptr);
      dat1_ptr++;
      dat2_ptr++;
    }
    fwrite(dat2, sizeof(float), sz_rem * Nreal, fp_new);

    delete dat1;
    delete dat2;
  } else if (_convertDBL) {
    float *dat1_ptr, *dat1;
    double *dat2, *dat2_ptr;
    dat1 = new float[sz_datstr / sizeof(double) * sz_chunk];
    dat2 = new double[sz_datstr / sizeof(double) * sz_chunk];

    /* do full chunks */
    for (k = 0; k < N_chunk; k++) {
      dat2_ptr = dat2;

      fread(dat1, sz_datstr / 2, sz_chunk, fp);
      dat1_ptr = dat1;

      for (i = 0; i < sz_chunk * Nreal; i++) {
        *dat2_ptr = (double)(*dat1_ptr);
        dat1_ptr++;
        dat2_ptr++;
      }
      fwrite(dat2, sizeof(double), sz_chunk * Nreal, fp_new);
    }

    /* do what's left over */
    fread(dat1, sz_datstr / 2, sz_rem, fp);
    dat1_ptr = dat1;
    dat2_ptr = dat2;

    for (i = 0; i < sz_rem * Nreal; i++) {
      *dat2_ptr = (double)(*dat1_ptr);
      dat1_ptr++;
      dat2_ptr++;
    }
    fwrite(dat2, sizeof(double), sz_rem * Nreal, fp_new);

    delete dat1;
    delete dat2;
  }

  fclose(fp);
  fclose(fp_new);
  strcpy(filename, filename_new);
}

void convertBinToAscii(char *filename, int prec, int _inclHdr) {
  int sz_chunk = 1 << 20;
  int szReal = (prec == 1 ? sizeof(float) : sizeof(double));

  FILE *fp, *fp_new;
  int i, j, k, Nrow = 1;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar, Nfloat;
  int N_chunk, sz_rem;
  char fieldType[21], filename_new[100];
  int sz_datstr, rel_sz_real, rel_sz_int, rel_sz_float;
  char *dat, *dat0;

  rel_sz_real = szReal / sizeof(char);
  rel_sz_int = sizeof(int) / sizeof(char);
  rel_sz_float = sizeof(nml_float) / sizeof(char);

  bzero(fieldType, sizeof(char) * 21);

  strcat(strcpy(filename_new, filename), "_ascii");
  if (!_inclHdr)
    strcat(filename_new, "_stripped");

  for (i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);

  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);

  writeHeader(filename_new, label, dim, dim_lbl, usr_comment);

  for (i = 0; i < Ndim; i++) {
    Nrow *= dim[i];
  }

  fp_new = fopen(filename_new, (_inclHdr ? "a" : "w"));
  fprintf(fp_new, "\n");

  /* break up reads/writes into smaller chunks */
  N_chunk = Nrow / (sz_chunk);
  sz_rem = Nrow % (sz_chunk);
  sz_datstr = Nreal * szReal + Nint * sizeof(int) + Nfloat * sizeof(nml_float);
  int rel_sz_datstr = sz_datstr / sizeof(char);
  dat0 = new char[rel_sz_datstr * sz_chunk];

  /* do full chunks */
  for (k = 0; k < N_chunk; k++) {
    fread(dat0, sz_datstr, sz_chunk, fp);
    dat = dat0;
    for (i = 0; i < sz_chunk; i++) {
      for (j = 0; j < Nreal + Nint + Nfloat; j++) {
        if (fieldType[j] == 'R') {
          if (prec == 2)
            fprintf(fp_new, "%21.12e", *((double *)dat));
          else
            fprintf(fp_new, "%16.7e", *((float *)dat));
          dat += rel_sz_real;
        }
        if (fieldType[j] == 'S') {
          fprintf(fp_new, "%16.7e", *((float *)dat));
          dat += rel_sz_float;
        }
        if (fieldType[j] == 'Z') {
          fprintf(fp_new, "%10d", *((int *)dat));
          dat += rel_sz_int;
        }
      }
      fprintf(fp_new, "\n");
    }
  }

  /* do what's left over */
  fread(dat0, sz_datstr, sz_rem, fp);
  dat = dat0;
  for (i = 0; i < sz_rem; i++) {
    for (j = 0; j < Nreal + Nint + Nfloat; j++) {
      if (fieldType[j] == 'R') {
        if (prec == 2)
          fprintf(fp_new, "%21.12e", *((double *)dat));
        else
          fprintf(fp_new, "%16.7e", *((float *)dat));
        dat += rel_sz_real;
      }
      if (fieldType[j] == 'S') {
        fprintf(fp_new, "%16.7e", *((float *)dat));
        dat += rel_sz_float;
      }
      if (fieldType[j] == 'Z') {
        fprintf(fp_new, "%10d", *((int *)dat));
        dat += rel_sz_int;
      }
    }
    fprintf(fp_new, "\n");
  }

  fclose(fp);
  fclose(fp_new);

  for (i = 0; i < MAX_VARS; i++) {
    delete dim_lbl[i];
    delete usr_comment[i];
  }

  delete dat0;
  strcpy(filename, filename_new);
}
#ifdef ConvertWfToPdb
/* this function has two arguments, the first one should be a wavefunction file
 * name, the second should be a pdb file name
 * the wavefunction should be in the binary format
 * the pdb file is created by running nemo3d with input option StructurePDB set
 * as true
 * the output of this function is a pdb format file with wavefunction values
 */

void convertWfToPdb(char *filename, char *filename2) {
  int sz_chunk = 1 << 20;
  int szReal = (sizeof(real));

  FILE *fp, *fp_new, *fp_pdb;
  int i, j, k, Nrow = 1;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar, Nfloat;
  int N_chunk, sz_rem;
  char fieldType[21], filename_new[200];
  int sz_datstr, rel_sz_real, rel_sz_int, rel_sz_float;
  char *dat, *dat0;
  float r1, r2, r3, r4, x, y, z;
  unsigned int na;
  int num;
  char str[50], str2[10];
  rel_sz_real = szReal / sizeof(char);
  rel_sz_int = sizeof(int) / sizeof(char);
  rel_sz_float = sizeof(nml_float) / sizeof(char);

  bzero(fieldType, sizeof(char) * 21);

  strcat(strcpy(filename_new, filename), "_pdb_ascii");

  for (i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);

  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);
  fp_pdb = fopen(filename2, "r");
  fp_new = fopen(filename_new, ("a"));
  for (i = 0; i < 6; i++) {
    fscanf(fp_pdb, "%s%f%f%f%f\n", str, &r1, &r2, &r3, &r4);
    fprintf(fp_new, "%s      %1.6f  %1.6f  %1.6f        %1.5f\n", str, r1, r2,
            r3, r4);
  }

  for (i = 0; i < Ndim; i++) {
    Nrow *= dim[i];
  }

  fprintf(fp_new, "\n");
  /* break up reads/writes into smaller chunks */
  N_chunk = Nrow / (sz_chunk);
  sz_rem = Nrow % (sz_chunk);
  sz_datstr = Nreal * szReal + Nint * sizeof(int) + Nfloat * sizeof(nml_float);
  int rel_sz_datstr = sz_datstr / sizeof(char);
  dat0 = new char[rel_sz_datstr * sz_chunk];

  /* do full chunks */
  for (k = 0; k < N_chunk; k++) {
    fread(dat0, sz_datstr, sz_chunk, fp);
    dat = dat0;
    for (i = 0; i < sz_chunk; i++) {
      for (j = 0; j < Nreal + Nint + Nfloat; j++) {
        {
          /* fscanf(fp_pdb, "%s  "FORMATOUT"%s%d
           * %f%f%f"COORDINATEREAD"\n",str2,&na,str,&num,&x,&y,&z,&r1,&r2);*/
          fscanf(fp_pdb, "%s  " FORMATOUT "%s%d    %f%f%f" COORDINATEREAD "\n",
                 str2, &na, str, &num, &x, &y, &z);
          fprintf(fp_new, "ATOM " FORMATOUT " %s%12d    %8.3f%8.3f%8.3f%14.7e",
                  na, str, num, x, y, z, *((float *)dat));
          fflush(fp_new);
          dat += rel_sz_real;
        }
      }
      fprintf(fp_new, "\n");
    }
  }
  /* do what's left over */
  fread(dat0, sz_datstr, sz_rem, fp);
  dat = dat0;
  for (i = 0; i < sz_rem; i++) {
    for (j = 0; j < Nreal + Nint + Nfloat; j++) {
      {
        /* fscanf(fp_pdb, "%s
         * "FORMATOUT"%s%d%f%f%f"COORDINATEREAD"\n",str2,&na,str,&num,&x,&y,&z,&r1,&r2);
         */
        fscanf(fp_pdb, "%s  " FORMATOUT "%s%d%f%f%f" COORDINATEREAD "\n", str2,
               &na, str, &num, &x, &y, &z);
        fprintf(fp_new, "ATOM " FORMATOUT " %s%12d    %8.3f%8.3f%8.3f%14.7e",
                na, str, num, x, y, z, *((double *)dat));
        fflush(fp_new);
        dat += rel_sz_real;
      }
    }
    fprintf(fp_new, "\n");
  }

  fclose(fp);
  fscanf(fp_pdb, "%s%d\n", str, &na);
  fprintf(fp_new, "%s  " FORMATOUT " \n", str, na);
  while (fscanf(fp_pdb, "%s%d%d\n", str, &na, &num) != EOF) {
    if (strcmp(str, "END") != 0)
      fprintf(fp_new, "%s " FORMATOUT " " FORMATOUT "\n", str, na, num);
    else
      fprintf(fp_new, "END");
  }
  fclose(fp_new);
  fclose(fp_pdb);

  for (i = 0; i < MAX_VARS; i++) {
    delete dim_lbl[i];
    delete usr_comment[i];
  }

  delete dat0;
  // strcpy(filename, filename_new);
}
#endif
void convertAsciiToBin(char *filename) {
  int sz_chunk = 1 << 20;

  FILE *fp, *fp_new;
  int i, j, k, Nrow = 1;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar, Nfloat;
  int N_chunk, sz_rem;
  char fieldType[21], filename_new[100];
  int sz_datstr, rel_sz_real, rel_sz_int, rel_sz_float;
  char *dat, *dat0;

  rel_sz_real = sizeof(real) / sizeof(char);
  rel_sz_int = sizeof(int) / sizeof(char);
  rel_sz_float = sizeof(nml_float) / sizeof(char);

  bzero(fieldType, sizeof(char) * 21);

  strcat(strcpy(filename_new, filename), "_bin");

  for (i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);

  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);

  writeHeader(filename_new, label, dim, dim_lbl, usr_comment);

  for (i = 0; i < Ndim; i++) {
    Nrow *= dim[i];
  }

  fp_new = fopen(filename_new, "a");
  /* break up reads/writes into smaller chunks */
  N_chunk = Nrow / (sz_chunk);
  sz_rem = Nrow % (sz_chunk);
  sz_datstr =
      Nreal * sizeof(double) + Nint * sizeof(int) + Nfloat * sizeof(nml_float);
  int rel_sz_datstr = sz_datstr / sizeof(char);
  dat0 = new char[rel_sz_datstr * sz_chunk];

  /* do full chunks */
  for (k = 0; k < N_chunk; k++) {
    dat = dat0;
    for (i = 0; i < sz_chunk; i++) {
      for (j = 0; j < Nreal + Nint + Nfloat; j++) {
        if (fieldType[j] == 'R') {
          fscanf(fp, "%le", ((real *)dat));
          dat += rel_sz_real;
        }
        if (fieldType[j] == 'S') {
          fscanf(fp, "%g", ((nml_float *)dat));
          dat += rel_sz_float;
        }
        if (fieldType[j] == 'Z') {
          fscanf(fp, "%d", ((int *)dat));
          dat += rel_sz_int;
        }
      }
    }
    fwrite(dat0, sz_datstr, sz_chunk, fp_new);
  }

  /* do what's left over */
  dat = dat0;
  for (i = 0; i < sz_rem; i++) {
    for (j = 0; j < Nreal + Nint + Nfloat; j++) {
      if (fieldType[j] == 'R') {
        fscanf(fp, "%le", ((real *)dat));
        dat += rel_sz_real;
      }
      if (fieldType[j] == 'S') {
        fscanf(fp, "%g", ((nml_float *)dat));
        dat += rel_sz_float;
      }
      if (fieldType[j] == 'Z') {
        fscanf(fp, "%d", ((int *)dat));
        dat += rel_sz_int;
      }
    }
  }
  fwrite(dat0, sz_datstr, sz_rem, fp_new);

  fclose(fp);
  fclose(fp_new);

  for (i = 0; i < MAX_VARS; i++) {
    delete dim_lbl[i];
    delete usr_comment[i];
  }

  delete dat0;
  strcpy(filename, filename_new);
}

// function is similar to convertBinToAscii however it extracts
// a matrix of data from the file named filename and reurns that matrix
double *getBinFileData(const char *filename, int prec, int *rdim, int *cdim) {
  int sz_chunk = 1 << 20;
  int szReal = (prec == 1 ? sizeof(float) : sizeof(double));

  FILE *fp;
  int i, j, k, Nrow = 1;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar, Nfloat;
  int N_chunk, sz_rem;
  char fieldType[21];
  int sz_datstr, rel_sz_real, rel_sz_int, rel_sz_float;
  char *dat, *dat0;
  double *result = NULL, *result_ptr;

  rel_sz_real = szReal / sizeof(char);
  rel_sz_int = sizeof(int) / sizeof(char);
  rel_sz_float = sizeof(nml_float) / sizeof(char);

  bzero(fieldType, sizeof(char) * 21);

  for (i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);

  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);

  for (i = 0; i < Ndim; i++) {
    Nrow *= dim[i];
  }

  result = (double *)calloc((Nrow * (Nreal + Nint + Nfloat)), sizeof(double));
  result_ptr = result;
  *rdim = Nrow;
  *cdim = Nreal + Nint + Nfloat;

  /* break up reads/writes into smaller chunks */
  N_chunk = Nrow / (sz_chunk);
  sz_rem = Nrow % (sz_chunk);
  sz_datstr = Nreal * szReal + Nint * sizeof(int) + Nfloat * sizeof(nml_float);
  int rel_sz_datstr = sz_datstr / sizeof(char);
  dat0 = new char[rel_sz_datstr * sz_chunk];

  /* do full chunks */
  for (k = 0; k < N_chunk; k++) {
    fread(dat0, sz_datstr, sz_chunk, fp);
    dat = dat0;
    for (i = 0; i < sz_chunk; i++) {
      for (j = 0; j < Nreal + Nint + Nfloat; j++) {
        if (fieldType[j] == 'R') {
          if (prec == 2)
            *result_ptr = *((double *)dat);
          else
            *result_ptr = *((float *)dat);
          dat += rel_sz_real;
          result_ptr++;
        }
        if (fieldType[j] == 'S') {
          *result_ptr = *((float *)dat);
          dat += rel_sz_float;
          result_ptr++;
        }
        if (fieldType[j] == 'Z') {
          *result_ptr = *((int *)dat);
          dat += rel_sz_int;
          result_ptr++;
        }
      }
    }
  }

  /* do what's left over */
  fread(dat0, sz_datstr, sz_rem, fp);
  dat = dat0;
  for (i = 0; i < sz_rem; i++) {
    for (j = 0; j < Nreal + Nint + Nfloat; j++) {
      if (fieldType[j] == 'R') {
        if (prec == 2)
          *result_ptr = *((double *)dat);
        else
          *result_ptr = *((float *)dat);
        dat += rel_sz_real;
        result_ptr++;
      }
      if (fieldType[j] == 'S') {
        *result_ptr = *((float *)dat);
        dat += rel_sz_float;
        result_ptr++;
      }
      if (fieldType[j] == 'Z') {
        *result_ptr = *((int *)dat);
        dat += rel_sz_int;
        result_ptr++;
      }
    }
  }

  fclose(fp);

  for (i = 0; i < MAX_VARS; i++) {
    delete dim_lbl[i];
    delete usr_comment[i];
  }

  delete dat0;

  return result;
}

// The definition below establishes an unmangled C function compiled with C++
// such that this function can be called from a C-linked executable
#ifdef __cplusplus
extern "C" {
#endif
double *getBinFileData_c(const char *filename, int prec, int *rdim, int *cdim) {
  double *result_vec = getBinFileData(filename, prec, rdim, cdim);
#if 0
  {
    int i,j;
    printf("cdim=%d rdim=%d\n",*cdim,*rdim);
    for (i=0;i< *rdim;i++){
      for (j=0;j< *cdim;j++){
	printf("%g ",result_vec[j+i*  *cdim]);
      }
      printf("\n");
    }
    printf("\n");
  }
#endif

  return result_vec;
}
#ifdef __cplusplus
}
#endif
