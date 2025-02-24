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
$Header: /repo/nemo3d/src/util/util_command_BX.c,v 1.3 2003/10/08 16:18:53 hook Exp $
*****************************************************************************/

#include "util_command_BX.h"


CMD_OPTIONS CmdOptions = { 
   NULL, 	/* input 	*/
   ARGV_FALSE, 	/* silent 	*/
   ARGV_FALSE, 	/* verbose 	*/
   NULL 	/* debuglaunch 	*/
};
   

#define input_arg_set ARGV_MAND


/* The command-line argument array.
 * Note: If in X Windows, DO NOT use the following arguments: -i, -b, -t
 */
argv_t args[] = {
  {input_arg_set, "input",    	 ARGV_CHARP, 	&CmdOptions.input,       "filename", 	"input filename" },
  { 'S',	  "silent",   	 ARGV_BOOL,  	&CmdOptions.silent,   	 NULL,    	"no status line output" },
  { 'v',          "verbose",     ARGV_BOOL,  	&CmdOptions.verbose,   	 NULL,    	"print out active command line options" },
  { 'D', 	  "debuglaunch", ARGV_CHARP, 	&CmdOptions.debuglaunch, "displayname",	"launch a debugger on display name" },

  { ARGV_LAST},
};


Boolean argv_check(void) {

   if (CmdOptions.verbose && (mpi_n3d_id==0)) {

#define print_bool(x) ( x ? printf("%s = TRUE\n", #x) :  printf("%s = FALSE\n", #x))
   
      printf("Input = %s\n", CmdOptions.input);
      printf("Debug Launch = %s\n", CmdOptions.debuglaunch);
   }
 
   return TRUE;
}


