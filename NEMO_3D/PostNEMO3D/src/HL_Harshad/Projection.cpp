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


#include "Projection.h"
#include "SimpleParser.h"

void Projection::get_data(Complex_t* coef_e, Complex_t* coef_h,
                      int* atomid, double* lat)
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _atom = atomid;
   _lattice = lat;
}

void Projection::compute_projection(void)
{
  double one_over_sqrt_three = 1.0/sqrt(3.0);
  double sqrt_two_over_sqrt_three = sqrt(2.0)/sqrt(3.0);
  double one_over_sqrt_six = 1.0/sqrt(6.0);

  Complex_t I(0.0,1.0);

  for(int n=0;n<e_max+h_max;n++){

    double Tmp0=0.0, Tmp1=0.0, Tmp2=0.0, Tmp3=0.0, Tmp4=0.0, Tmp5=0.0, Tmp6=0.0, Tmp7=0.0;

    for(int g=g_min;g<g_max;g++){

      for(int s=0; s<2; s++){
        for(int orbital=0; orbital<o_max; orbital++) {
          if(n<e_max) coef(s,orbital)=cf_e(n,g,s,orbital);
          else        coef(s,orbital)=cf_h(n-e_max,g,s,orbital);
        }
      }


      Tmp0+=( coef(0,0)*conj(coef(0,0)) ).real();
      Tmp0+=( coef(1,0)*conj(coef(1,0)) ).real(); // s-band

      Tmp1+=( coef(0,1)*conj(coef(0,1)) ).real();
      Tmp1+=( coef(1,1)*conj(coef(1,1)) ).real(); // s*-band
      Complex_t hh_u = coef(0,2)-I*coef(0,3); //heavy-hole J=3/2, Jz=3/2
      Complex_t hh_d = coef(1,2)+I*coef(1,3); //heavy-hole J=3/2, Jz=-3/2
      Complex_t lh_u = sqrt_two_over_sqrt_three*coef(0,4) 
                     - one_over_sqrt_six*( coef(1,2)-I*coef(1,3) ); //light-hole J=3/2, Jz=1/2
      Complex_t lh_d = sqrt_two_over_sqrt_three*coef(1,4) 
                     + one_over_sqrt_six*( coef(0,2)+I*coef(0,3) ); //light-hole J=3/2, Jz=-1/2
      Complex_t ss_u = one_over_sqrt_three*( coef(0,4)+coef(1,2)-I*coef(1,3) ); //split-off J=1/2, Jz=1/2
      Complex_t ss_d = one_over_sqrt_three*( coef(0,2)+I*coef(0,3)-coef(1,4) ); //split-off J=1/2, Jz=-1/2   
      
      Tmp2+=( hh_u*conj(hh_u) + hh_d*conj(hh_d) ).real(); 
      Tmp3+=( lh_u*conj(lh_u) + lh_d*conj(lh_d) ).real(); 
      Tmp4+=( ss_u*conj(ss_u) + ss_d*conj(ss_d) ).real(); 

      if(o_max==10){
        for(int o1=5; o1<10; o1++){
          Tmp5+=( coef(0,o1)*conj(coef(0,o1)) ).real();
          Tmp5+=( coef(1,o1)*conj(coef(1,o1)) ).real(); //d-band
        }
      }

      for(int orb=0; orb<o_max; orb++){
        Tmp6+= ( coef(0,orb)*conj(coef(0,orb)) ).real(); // spin up
        Tmp7+= ( coef(1,orb)*conj(coef(1,orb)) ).real(); // spin down
      }       

    }

    projection(n,0)=Tmp0;
    projection(n,1)=Tmp1;
    projection(n,2)=Tmp2*0.5;
    projection(n,3)=Tmp3;
    projection(n,4)=Tmp4;
    projection(n,5)=Tmp5;
    projection(n,6)=Tmp6;
    projection(n,7)=Tmp7;
  }
}
