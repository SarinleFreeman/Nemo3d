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
$Header: /repo/nemo3d/src/postprocessing/AngularMomentum.cpp,v 1.1 2004/08/26 21:28:07 swlee Exp $
*****************************************************************************/

#include "AngularMomentum.h"

void AngularMomentum::get_data(complex* coef_e, complex* coef_h,
                      int* atom_id, int* shape_rank, int* nn, double* lat)
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _atom = atom_id;
   _shape = shape_rank;
   _lattice = lat;
   _nn = nn;

}

void AngularMomentum::set_origin(double X, double Y, double Z)
{
  X0 = X; Y0 = Y; Z0 = Z;
  cout<<"origin for angular momentum operator: "<<X0<<" "<<Y0<<" "<<Z0<<endl;
}

void AngularMomentum::compute_angular_momentum_spin(const char& direction)
{
   cout<<endl;
   cout<<"=========================================================="<<endl;
   for(int ie=0; ie<e_max; ie++){
     complex Le=compute_global_angular_momentum(ie, 'e', direction);
     complex le=compute_local_angular_momentum(ie, ie, 'e', direction);
     complex Se=compute_spin(ie, ie, 'e', direction);
     if(my_id==0) {
       cout<<"Electron Level "<<ie+1<<endl;
       cout<<"Global Angular Momentum: "<< Le.r <<endl;
       cout<<"Local Angular Momentum: "<< le.r <<endl;
       cout<<"Spin: "<< Se.r <<endl;
       cout<<endl;
     }
   }
   for(int ih=0; ih<h_max; ih++){
     complex Lh=compute_global_angular_momentum(ih, 'h', direction);
     complex lh=compute_local_angular_momentum(ih, ih, 'h', direction);
     complex Sh=compute_spin(ih, ih, 'h', direction);
     if(my_id==0) {
       cout<<"Hole Level "<<ih+1<<endl;
       cout<<"Global Angular Momentum: "<< Lh.r<<endl;
       cout<<"Local Angular Momentum: "<< lh.r<<endl;
       cout<<"Spin: "<< Sh.r<<endl;
       cout<<endl;
     }
   }
}

complex AngularMomentum::compute_global_angular_momentum_square(int n1,
		         const char& carrier, const char& direction)
{
  complex global_L2=Complex(0.0,0.0);
  if(carrier!='e' && carrier!='h') {
    cout<<"Choose carrier-type either e for electron or h for hole"<<endl;
    return global_L2;
  }
  else if(direction!='x' && direction!='y' && direction!='z'){
    cout<<"Choose direction x, y, or z for angular momentum"<<endl;
    return global_L2;
  }
  else {
   cvectr X_WF = Cvectr(d->seg_ln[my_id]);
   cvectr Y_WF = Cvectr(d->seg_ln[my_id]);
   cvectr WF = Cvectr(d->seg_ln[my_id]);
    complex cf;
    int n=0;
    for(int iat=g_min; iat<g_max; iat++){
      double X, Y;
      X=lattice(iat,0)-X0; 
      Y=lattice(iat,1)-Y0; 
      for(int spin=0; spin<2; spin++){
        for(int iorb=0; iorb<o_max; iorb++){
          if(carrier=='e') cf = cf_e(n1,iat,spin,iorb);
          else             cf = cf_h(n1,iat,spin,iorb);
          WF[n].r = cf.r;
          WF[n].i = cf.i;
          X_WF[n].r = cf.r*X; 
          X_WF[n].i = cf.i*X; 
          Y_WF[n].r = cf.r*Y; 
          Y_WF[n].i = cf.i*Y; 
          n++;
        }
      }
    }
    cvectr H_X_WF = Cvectr(d->seg_ln[my_id]); 
    cvectr H_Y_WF = Cvectr(d->seg_ln[my_id]); 
    Hmult_spds_complete(H_X_WF, d, 0.0, X_WF);
    Hmult_spds_complete(H_Y_WF, d, 0.0, Y_WF);
    cvectr Lz_WF = Cvectr(d->seg_ln[my_id]);
    cvectr X_Lz_WF = Cvectr(d->seg_ln[my_id]);
    cvectr Y_Lz_WF = Cvectr(d->seg_ln[my_id]);
    n=0;
    for(int iat=g_min; iat<g_max; iat++){
      double X, Y;
      X=lattice(iat,0)-X0; 
      Y=lattice(iat,1)-Y0; 
      for(int spin=0; spin<2; spin++){
        for(int iorb=0; iorb<o_max; iorb++){
          Lz_WF[n].r = X*H_Y_WF[n].r - Y*H_X_WF[n].r;
          Lz_WF[n].i = X*H_Y_WF[n].i - Y*H_X_WF[n].i;
          X_Lz_WF[n].r = X*Lz_WF[n].r;
          X_Lz_WF[n].i = X*Lz_WF[n].i;
          Y_Lz_WF[n].r = Y*Lz_WF[n].r;
          Y_Lz_WF[n].i = Y*Lz_WF[n].i;
          n++;
        }
      }
    } 
    cvectr H_X_Lz_WF = Cvectr(d->seg_ln[my_id]);
    cvectr H_Y_Lz_WF = Cvectr(d->seg_ln[my_id]);
    Hmult_spds_complete(H_X_Lz_WF, d, 0.0, X_Lz_WF);
    Hmult_spds_complete(H_Y_Lz_WF, d, 0.0, Y_Lz_WF);
    complex XHYLz = vect_dot_vect(X_WF, H_Y_Lz_WF, d);
    complex YHXLz = vect_dot_vect(Y_WF, H_X_Lz_WF, d);
    global_L2.r = -174.77*(XHYLz.r-YHXLz.r);
    global_L2.i = -174.77*(XHYLz.i-YHXLz.i);

    rm_cvectr(&WF);
    rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&H_X_WF); rm_cvectr(&H_Y_WF);
    rm_cvectr(&Lz_WF); rm_cvectr(&X_Lz_WF); rm_cvectr(&Y_Lz_WF); 
    rm_cvectr(&H_X_Lz_WF); rm_cvectr(&H_Y_Lz_WF);
    return global_L2;
  }
}

complex AngularMomentum::compute_global_angular_momentum(int n1, 
                         const char& carrier, const char& direction)
{
  complex global_L=Complex(0.0,0.0);
  if(carrier!='e' && carrier!='h') {
    cout<<"Choose carrier-type either e for electron or h for hole"<<endl;
    return global_L;
  }
  else if(direction!='x' && direction!='y' && direction!='z'){
    cout<<"Choose direction x, y, or z for angular momentum"<<endl;
    return global_L;
  }
  else {
    cvectr X_WF = Cvectr(d->seg_ln[my_id]);
    cvectr Y_WF = Cvectr(d->seg_ln[my_id]);
    cvectr Z_WF = Cvectr(d->seg_ln[my_id]);
    cvectr WF = Cvectr(d->seg_ln[my_id]);
    complex cf;
    int n=0;
    for(int iat=g_min; iat<g_max; iat++){
      double X, Y, Z;
      X=lattice(iat,0)-X0; 
      Y=lattice(iat,1)-Y0; 
      Z=lattice(iat,2)-Z0; 
      for(int spin=0; spin<2; spin++){
        for(int iorb=0; iorb<o_max; iorb++){
          if(carrier=='e') cf = cf_e(n1,iat,spin,iorb);
          else             cf = cf_h(n1,iat,spin,iorb);
          X_WF[n].r = cf.r*X; 
          X_WF[n].i = cf.i*X; 
          Y_WF[n].r = cf.r*Y; 
          Y_WF[n].i = cf.i*Y; 
          Z_WF[n].r = cf.r*Z; 
          Z_WF[n].i = cf.i*Z; 
          WF[n].r = cf.r;
          WF[n].i = cf.i;
          n++;
        }
      }
    }
    cvectr H_WF = Cvectr(d->seg_ln[my_id]);
    Hmult_spds_complete(H_WF, d, 0.0, WF);
/*
    complex WF_H_WF = vect_dot_vect(WF, H_WF, d);
    complex WF_WF = vect_dot_vect(WF, WF, d);
    masterPrint("Energy of %i level: %e",n1, WF_H_WF.r);  
    masterPrint("Wavefunction Norm of %i level: %e",n1, WF_WF.r);  
*/
    if(direction=='x'){
      cvectr H_Y_WF = Cvectr(d->seg_ln[my_id]);
      cvectr H_Z_WF = Cvectr(d->seg_ln[my_id]);
      Hmult_spds_complete(H_Y_WF, d, 0.0, Y_WF);
      Hmult_spds_complete(H_Z_WF, d, 0.0, Z_WF);
      complex Y_H_Z = vect_dot_vect(Y_WF, H_Z_WF, d);
      complex Z_H_Y = vect_dot_vect(Z_WF, H_Y_WF, d);
      complex global_L;
      global_L.r=13.22*(Z_H_Y.i-Y_H_Z.i);
      global_L.i=13.22*(Y_H_Z.r-Z_H_Y.r);
      rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&Z_WF); rm_cvectr(&WF);
      rm_cvectr(&H_Y_WF); rm_cvectr(&H_Z_WF);
      return global_L;
    }
    else if(direction=='y'){
      cvectr H_X_WF = Cvectr(d->seg_ln[my_id]);
      cvectr H_Z_WF = Cvectr(d->seg_ln[my_id]);
      Hmult_spds_complete(H_X_WF, d, 0.0, X_WF);
      Hmult_spds_complete(H_Z_WF, d, 0.0, Z_WF);
      complex Z_H_X = vect_dot_vect(Z_WF, H_X_WF, d);
      complex X_H_Z = vect_dot_vect(X_WF, H_Z_WF, d);
      complex global_L;
      global_L.r=13.22*(X_H_Z.i-Z_H_X.i);
      global_L.i=13.22*(Z_H_X.r-X_H_Z.r);
      rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&Z_WF); rm_cvectr(&WF);
      rm_cvectr(&H_X_WF); rm_cvectr(&H_Z_WF);
      return global_L;
    }
    else {
      cvectr H_X_WF = Cvectr(d->seg_ln[my_id]);
      cvectr H_Y_WF = Cvectr(d->seg_ln[my_id]);
      Hmult_spds_complete(H_X_WF, d, 0.0, X_WF);
      Hmult_spds_complete(H_Y_WF, d, 0.0, Y_WF);
      complex X_H_Y = vect_dot_vect(X_WF, H_Y_WF, d);
      complex Y_H_X = vect_dot_vect(Y_WF, H_X_WF, d);
      complex global_L;
      global_L.r=13.22*(Y_H_X.i-X_H_Y.i);
      global_L.i=13.22*(X_H_Y.r-Y_H_X.r);
      rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&Z_WF); rm_cvectr(&WF);
      rm_cvectr(&H_X_WF); rm_cvectr(&H_Y_WF);
      return global_L;
    }
  }
}

complex AngularMomentum::compute_local_angular_momentum(int n1, int n2, 
                         const char& carrier, const char& direction)
{
  complex local_l=Complex(0.0,0.0);
  if(carrier!='e' && carrier!='h') {
    cout<<"Choose carrier-type either e for electron or h for hole"<<endl;
    return local_l;
  }
  else if(direction!='x' && direction!='y' && direction!='z'){
    cout<<"Choose direction x, y, or z for angular momentum"<<endl;
    return local_l;
  }
  else {
    complex c1, c2, l;
    for(int at=g_min; at<g_max; at++){
      for(int spin=0; spin<2; spin++){
        for(int iorb=0; iorb<o_max; iorb++){
          if(carrier=='e') c1=cf_e(n1,at,spin,iorb);
          else c1=cf_h(n1,at,spin,iorb);
          for(int jorb=0; jorb<o_max; jorb++){ 
            if(carrier=='e') c2=cf_e(n2,at,spin,jorb);
            else c2=cf_h(n2,at,spin,jorb);
            if(direction=='x')      l = MatrixLx[iorb*10+jorb];
            else if(direction=='y') l = MatrixLy[iorb*10+jorb];
            else                    l = MatrixLz[iorb*10+jorb];
            if(l.i!=0.0) {
              local_l.r = local_l.r - (c1.r*c2.i-c1.i*c2.r)*l.i;
              local_l.i = local_l.i + (c1.r*c2.r+c1.i*c2.i)*l.i;
            }
          } 
        }  
      }
    }
    complex total_l;
    MPI_Allreduce(&local_l, &total_l, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return total_l;
  }
}

complex AngularMomentum::compute_spin(int n1, int n2, 
                         const char& carrier, const char& direction)
{
  complex local_s = Complex(0.0,0.0);
  if(carrier!='e' && carrier!='h') {
    cout<<"Choose carrier-type either e for electron or h for hole"<<endl;
    return local_s;
  }
  else if(direction!='x' && direction!='y' && direction!='z'){
    cout<<"Choose direction x, y, or z for angular momentum"<<endl;
    return local_s;
  }
  else{
    for(int at=g_min; at<g_max; at++){
      for(int iorb=0; iorb<o_max; iorb++){
        complex c1_up, c1_down, c2_up, c2_down; 
        if(carrier=='e') {
          c1_up=cf_e(n1,at,0,iorb);
          c2_up=cf_e(n2,at,0,iorb);
          c1_down=cf_e(n1,at,1,iorb);
          c2_down=cf_e(n2,at,1,iorb);
        }
        else {
          c1_up=cf_h(n1,at,0,iorb);
          c2_up=cf_h(n2,at,0,iorb);
          c1_down=cf_h(n1,at,1,iorb);
          c2_down=cf_h(n2,at,1,iorb);
        } 
        if(direction=='x'){
          local_s.r +=c1_up.r*c2_down.r+c1_up.i*c2_down.i
                   +c1_down.r*c2_up.r+c1_down.i*c2_up.i;
          local_s.i +=c1_up.r*c2_down.i-c1_up.i*c2_down.r
                   +c1_down.r*c2_up.i-c1_down.i*c2_up.r;
        }
        else if(direction=='y'){
          local_s.r +=c1_up.r*c2_down.i-c1_up.i*c2_down.r
                   +c1_down.i*c2_up.r-c1_down.r*c2_up.i;
          local_s.i +=c1_down.r*c2_up.r+c1_down.i*c2_up.i
                   -c1_up.r*c2_down.r-c1_up.i*c2_down.i;
        }
        else {
          local_s.r +=c1_up.r*c2_up.r+c1_up.i*c2_up.i
                   -c1_down.r*c2_down.r-c1_down.i*c2_down.i;
          local_s.i +=c1_up.r*c2_up.i-c1_up.i*c2_up.r
                    -c1_down.r*c2_down.i+c1_down.i*c2_down.r; 
        }
      }
    }
    local_s.r *=0.5; local_s.i *=0.5;
    complex total_s;
    MPI_Allreduce(&local_s, &total_s, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return total_s;  
  }
}

void AngularMomentum::prepare_MatrixL(void)
{
  MatrixLx = Cvectr(100);
  MatrixLy = Cvectr(100);
  MatrixLz = Cvectr(100);
  for(int i=0; i<10; i++)
    for(int j=0; j<10; j++) {
      MatrixLx[i*10+j].r = 0.0;
      MatrixLx[i*10+j].i = 0.0;
      MatrixLy[i*10+j].r = 0.0;
      MatrixLy[i*10+j].i = 0.0;
      MatrixLz[i*10+j].r = 0.0;
      MatrixLz[i*10+j].i = 0.0;
  }
  MatrixLz[23].i=-1.0;
  MatrixLz[32].i=1.0;
  MatrixLz[58].i=2.0;
  MatrixLz[67].i=1.0;
  MatrixLz[76].i=-1.0;
  MatrixLz[85].i=-2.0;

  MatrixLx[34].i=-1.0;
  MatrixLx[43].i=1.0;
  MatrixLx[57].i=-1.0;
  MatrixLx[68].i=-1.0;
  MatrixLx[69].i=-sqrt(3.0);
  MatrixLx[75].i=1.0;
  MatrixLx[86].i=1.0;
  MatrixLx[96].i=sqrt(3.0);

  MatrixLy[24].i=1.0;
  MatrixLy[42].i=-1.0;
  MatrixLy[56].i=1.0;
  MatrixLy[65].i=-1.0;
  MatrixLy[78].i=-1.0;
  MatrixLy[79].i=sqrt(3.0);
  MatrixLy[87].i=1.0;
  MatrixLy[97].i=-sqrt(3.0);
}
