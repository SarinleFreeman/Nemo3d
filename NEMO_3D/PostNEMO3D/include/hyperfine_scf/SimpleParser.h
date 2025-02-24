/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D
PostProcessing package.
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
  Seungwon Lee (Seungwon.Lee@jpl.nasa.gov)

Written by:  Seungwon Lee
*****************************************************************************/


#ifndef _SIMPLEPARSER_H
#define _SIMPLEPARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
//#include <strstream.h> /* This header file is no longer available after         */
                         /* Lear update and it is not used in the postprocess     */
                         /* code. So commenting it out doesn't cause any problem. */
                         /* If uncommented it gives an error */
#include <sstream>
using namespace std;

char* readLine(char *s, int max, istream &fp);
int getwords(vector<string>& slist, istream &fp);
unsigned parsewords(char *inbuf, vector<string>& slist);
unsigned parsewords(char *inbuf, list<string>& slist);

#endif

