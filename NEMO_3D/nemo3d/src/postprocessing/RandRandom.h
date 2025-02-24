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
$Header: /repo/nemo3d/src/postprocessing/RandRandom.h,v 1.1 2004/08/26 21:28:07 swlee Exp $
*****************************************************************************/

#ifndef RANDRANDOM_H
#define RANDRANDOM_H
#include <math.h>
#include <stdlib.h>

/*! \class RandRandom
 *  A wrapper class to generate a random number using C srand and rand.
 */
class RandRandom {

public:
  typedef double Return_t;

  RandRandom(bool gen_seed = false):thisStreamID(0), nStreams(1)
  {
    init(gen_seed);
  }

  RandRandom(int i, int nstr, bool gen_seed=false):thisStreamID(i), nStreams(nstr)
  {
    init(gen_seed);
  }
 
  ~RandRandom() { }

  void init(bool gen_seed=false);

  inline Return_t getRandom() { return Return_t(rand())*rand_max_inv;}

  inline Return_t operator()() { return getRandom(); }

  inline int irand() { return rand();}

  inline void bivariate(Return_t& g1, Return_t& g2) {
    Return_t v1, v2, r;
    do {
    v1 = 2.0e0*getRandom() - 1.0e0;
    v2 = 2.0e0*getRandom() - 1.0e0;
    r = v1*v1+v2*v2;
    } while(r > 1.0e0);
    Return_t fac = sqrt(-2.0e0*log(r)/r);
    g1 = v1*fac;
    g2 = v2*fac;
  }

private:
  static const Return_t rand_max_inv;
  int thisStreamID;
  int nStreams;
  int thisSeed;
};
#endif
