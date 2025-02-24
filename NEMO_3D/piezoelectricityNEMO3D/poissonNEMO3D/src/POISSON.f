CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Copyright (C) 2005 Shaikh Shahid Ahmed
C This application is free software, which you can redistribute and/or modify
C under the terms of the GNU Lesser General Public License as published by the
C Free Software Foundation; either version 2.1 of the License, or (at your
C option) any later version.
C
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C
C You should have received a copy of the GNU Lesser General Public License
C along with this library; see the file COPYING. If not, write to the
C Free Software Foundation, Inc.,
C 59 Temple Place, Suite 330,
C Boston, MA  02111-1307  USA
C
C For additional information, please contact
C Shaikh S. Ahmed 
C Department of Electrical and Computer Engineering
C Purdue University
C West Lafayette, IN, USA.
C Email: ssahmed@purdue.edu
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      MAIN PROGRAM FOR 3D POISSON SOLVER ON A ZINCBLENDE LATTICE    C
C      SHAIKH SHAHID AHMED/PURDUE UNIVERSITY/Jan. 2006               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit real*8(A-H, O-Z)
     
      include 'params.h'
      include 'common_par.inc'
      
C     Grids:

      real*8 X_VEC(-1:i_max)
      real*8 Y_VEC(-1:j_max)
      real*8 Z_VEC(k_neg:k_max)

C     Doping and others:

      real*8 DOP(0:i_max,0:j_max,k_neg:k_max)
      integer FLAG_DOM(0:i_max,0:j_max,k_neg:k_max)
      real*8 RO(0:i_max,0:j_max,k_neg:k_max)
      real*8 FI(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 N(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 P(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 eps(-1:i_max+1,0:j_max,k_neg:k_max+1)

C     Coefficients from the finite difference scheme:
C     (Poisson's equation when converted to algebraic equations)

      real*8 B_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 C_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 D_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 E_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 F_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 G_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 H_vec(0:i_max,0:j_max,k_neg:k_max)

      real*8 pot_err(1:5000)
 
      character*100 device_name
      integer k_iter

      NAMELIST /PARAM_LIST/iamdummy



C.....Program begins here
C.....Data:
     
C.....'input_poisson3d' file must be copied to the run directory      

      open(unit=1,file='input_poisson3d',status='old')

      rewind(1)
      read(1,NML=PARAM_LIST)
      read(1,*)device_name
      read(1,*)V_max
      read(1,*)eps_o
      read(1,*)eps_si
      read(1,*)eps_ox
      read(1,*)eps_p
      read(1,*)eps_inas
      read(1,*)eps_ingaas
      read(1,*)eps_gaas
      read(1,*)eps_sige
      read(1,*)eps_si_strained


C.....Used format
     
      eps_si = eps_si*eps_o
      eps_ox = eps_ox*eps_o
      eps_p = eps_p*eps_o
      eps_inas = eps_inas*eps_o
      eps_ingaas = eps_ingaas*eps_o
      eps_gaas = eps_gaas*eps_o
      eps_sige = eps_sige*eps_o
      eps_si_strained = eps_si_strained*eps_o

      q = 1.61e-19
      tolf  = V_max
      debye = 1.0D0
      Rni = 1.0
C     Normalization purpose     
      T = 300.0 
      pi  = 4.D0*DATAN(1.D0)
      k_b = 1.38066D-23
      h   = 6.62617D-34
      hb  = h/2.D0/pi
      Vt = k_b*T/q
C     Normalization parameters end
    
      if(device_name.eq.'SiP_DEF'.or.device_name.eq.'SiP_HOON')then
        Rni = 1.45D16
        tolf  = V_max/Vt
        debye = dsqrt(2.D0*eps_o*Vt/(q*RNi)) 
C        Vt = 1.0
C        Rni = 1.0
C        tolf  = V_max/Vt
C        debye = 1.0
      endif

      print*,'Potential covergence (tolerance)  = ',tolf

C.....CALL THE SUBROUTINE THAT SETS UP THE STRUCTURE:

      CALL STRUCTURE (device_name,X_VEC,Y_VEC,Z_VEC,
     &                DOP,FLAG_DOM,RO,FI,N,P,eps)
      print*,'Done with initialization of the structure'
      print*,'Overall No. of grid points (nx,ny,nz): ',nx,ny,nz


C.....Poisson coefficent calculation

      CALL COEFF_P(X_VEC, Y_VEC, Z_VEC, FLAG_DOM,
     &             B_VEC, C_VEC, D_VEC, E_VEC, F_VEC, G_VEC, H_VEC
     &            ,eps)
      print*,'Calculated the Poisson solver coefficients'
      

C.....POISSON EQUATION SOLUTION: (Solving linearized Equation)   

      k_iter = 0
      delta_max = 5.D0  
C.....(just to enter the loop a large value is set first)

      DO WHILE((DELTA_MAX.GT.TOLF))

       k_iter = k_iter+1  
       CALL ILU_POISSON(Y_VEC, B_VEC, C_VEC,
     1    D_VEC, E_VEC, F_VEC, G_VEC, H_VEC, FLAG_DOM,
     2    RO, FI, N, P)

C.....Calculate charge based on the new potential:      
       do i = 0, nx
        do j = 0, ny
         do k = -k_ox, nz
           if(flag_dom(i,j,k).eq.4)then
            p(i,j,k) = 0.D-5
            n(i,j,k) = 0.D-5
           if(device_name.eq.'SiP_DEF'.or.device_name.eq.'SiP_HOON')then
            n(i,j,k) = exp(fi(i,j,k))
            p(i,j,k) = exp(-fi(i,j,k))
           endif
            ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
           endif
         enddo
        enddo
       enddo

       print *,'Iteration #: ', k_iter
       print *,'Grid point with max. error: ', i_err, j_err, k_err
       print *,'Max. error', delta_max

      ENDDO  

C.....Creating data files:
      print*,'Writing down results in data files'
      print*,'  '
      CALL write_data_eff(device_name,x_vec,y_vec,z_vec,n,p,dop,
     &                    fi,eps,flag_dom)


C.....Final comments:
      print *,' '
      print *,'Number of iterations in the potential loop: ',k_iter
      print *,'End of program...'
      print *,'------------------------------------------'

      END

       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SUBROUTINE STRUCTURE.F IN WHICH THE LAYOUT OF THE STRUCTURE       C
C         IS DEFINED AND APPROPRIATE BOUNDARY CONDITIONS ARE SET          C
C          Definition of various domains in the device:                   C
C                flag_dom = 4 (bulk region)                               C
C                flag_dom = 5 (electrodes)                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE STRUCTURE(device_name,X_VEC,Y_VEC,Z_VEC,
     &                     DOP,FLAG_DOM,RO,FI,N,P,eps)
      implicit real*8(A-H,O-Z)

      include 'params.h'
      include 'common_par.inc'

C.....Grids:
      real*8 x_vec(-1:i_max)
      real*8 y_vec(-1:j_max)
      real*8 z_vec(k_neg:k_max)

C.....Others:
      real*8 DOP(0:i_max,0:j_max,k_neg:k_max)
      integer FLAG_DOM(0:i_max,0:j_max,k_neg:k_max)
      real*8 RO(0:i_max,0:j_max,k_neg:k_max)
      real*8 FI(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 N(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 P(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 eps(-1:i_max+1,0:j_max,k_neg:k_max+1)

      real*8 flag_reg(-1:i_max+1,0:j_max,k_neg:k_max+1)

      real*8 xion(number_appx),yion(number_appx)
      real*8 zion(number_appx),qion(number_appx)
      real*8 flag(number_appx),atom(number_appx)
      real*8 x_dop(number_appx),y_dop(number_appx)
      real*8 z_dop(number_appx),charge_dop(number_appx)

      real*8 px(number_appx),py(number_appx),pz(number_appx)
      real*8 polx(-1:i_max+1,0:j_max,k_neg:k_max+1)    
      real*8 poly(-1:i_max+1,0:j_max,k_neg:k_max+1)    
      real*8 polz(-1:i_max+1,0:j_max,k_neg:k_max+1)    

      real*8 rion_number(-1:i_max+1,0:j_max,k_neg:k_max+1)    
      integer iflag_atom(-1:i_max+1,0:j_max,k_neg:k_max+1)    
      character*100 device_name
      character(100) extension,FN
      real*8 norm_doping,atom_number,norm_factor 
      
      character*80 strain_file_cation,atom_type_file
      character*80 ratom_file
      integer argc

      NAMELIST /STRUC_LIST_DOT/iamdummy
      NAMELIST /STRUC_LIST_CQD/iamdummy
      NAMELIST /STRUC_LIST_SiP/iamdummy
      NAMELIST /STRUC_LIST_SiP_DEF/iamdummy
      NAMELIST /STRUC_LIST_SiP_HOON/iamdummy
      


   5  format(2x,9(e12.6,2x))
   8  format(2x,8(e12.6,2x))
   6  format(9(e12.6,1x))

   51 format(2x,12(e12.6,2x))
   81 format(2x,11(e12.6,2x))

C.....Random number generation constant
      iso = 1345

C.....Read data for a particular device 
      select case (device_name)
      case('DOT')
        FN = 'DOT'
        rewind(1)
        read(1,NML=STRUC_LIST_DOT)
        read(1,*)a0  
        read(1,*)discret_factor_x 
        read(1,*)discret_factor_y 
        read(1,*)discret_factor_z 
        read(1,*)piezo_model 
        read(1,*)interpolation 
        read(1,*)z_contour 
      case('CQD')
        FN = 'CQD'
        rewind(1)
        read(1,NML=STRUC_LIST_CQD)
        read(1,*)a0  
        read(1,*)discret_factor_x 
        read(1,*)discret_factor_y 
        read(1,*)discret_factor_z 
        read(1,*)z_ox 
        read(1,*)k_ox 
        read(1,*)piezo_model 
        read(1,*)interpolation 
        read(1,*)z_contour 
      case('SiP')
        FN = 'SiP'
        rewind(1)
        read(1,NML=STRUC_LIST_SiP)
        read(1,*)a0  
        read(1,*)discret_factor_x 
        read(1,*)discret_factor_y 
        read(1,*)discret_factor_z 
        read(1,*)n_dopants 
        read(1,*)x_dop(1) 
        read(1,*)y_dop(1) 
        read(1,*)z_dop(1) 
        read(1,*)charge_dop(1) 
        read(1,*)x_dop(2) 
        read(1,*)y_dop(2) 
        read(1,*)z_dop(2) 
        read(1,*)charge_dop(2) 
        read(1,*)z_ox 
        read(1,*)k_ox 
        read(1,*)z_contour 
        read(1,*)gate_voltage 
      case('SiP_DEF')
        FN = 'SiP_DEF'
        rewind(1)
        read(1,NML=STRUC_LIST_SiP_DEF)
        read(1,*)x_length 
        read(1,*)y_width 
        read(1,*)z_depth 
        read(1,*)a0  
        read(1,*)discret_factor_x 
        read(1,*)discret_factor_y 
        read(1,*)discret_factor_z 
        read(1,*)n_dopants 
        read(1,*)x_dop(1) 
        read(1,*)y_dop(1) 
        read(1,*)z_dop(1) 
        read(1,*)charge_dop(1) 
        read(1,*)x_dop(2) 
        read(1,*)y_dop(2) 
        read(1,*)z_dop(2) 
        read(1,*)charge_dop(2) 
        read(1,*)z_ox 
        read(1,*)k_ox 
        read(1,*)z_contour 
        read(1,*)gate_voltage 
      case('SiP_HOON')
        FN = 'SiP_HOON'
        rewind(1)
        read(1,NML=STRUC_LIST_SiP_HOON)
        read(1,*)a0  
        read(1,*)discret_factor_x 
        read(1,*)discret_factor_y 
        read(1,*)discret_factor_z 
        read(1,*)z_contour 
        read(1,*)y_contour 
        read(1,*)x_contour 
      endselect

      close(1)
C....(So here we close the input file and end of reading)

      extension = trim(FN)
      strain_file = trim(FN)



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....InAs/GaAs Single/Coupled Dots for Piezoelectric calculation
C.....No gate oxide and applied potential is considered
C.....Steps: (1) readin the strain information and domain size and grids 
C.....       (2) calculate the polarization
C.....       (3) calculate the charge density

      select case(device_name)

      case('DOT')
         
      z_ox = 0.0
      k_ox = 0
      dx1 = a0/discret_factor_x
      dy1 = a0/discret_factor_y
      dz1 = a0/discret_factor_z

      istrainfile = 1       

      if(istrainfile.eq.1)then       
         argc = iargc()
         if(argc.ne.3) then
          write(*,*)'Usage: ../bin/poisson3d.ex <strain_file> 
     &               <strain_file_cation> <atom_file>'
          write(*,*)' '
          call exit(1)
         endif

C........Reading the strain3d file 
         call getarg(1,strain_file)
C.......(need to know total atom number and domain size)
C         call domain_grid(z_ox,nn,x_vec,y_vec,z_vec) 

         xmax = -1.d38
         xmin = 1.d38
         ymax = -1.d38
         ymin = 1.d38
         zmax = -1.d38
         zmin = 1.d38

         nn = 0
C.......(nn is the atom counter)
         open(unit=1,file=strain_file,status='old')
   3     read(1,*,err=2,end=2)x,y,z,exx,eyy,ezz,exy,exz,eyz
         nn=nn+1
         x=x*1.d-9
         xmax = dmax1(x,xmax)
         xmin = dmin1(x,xmin)
         y=y*1.d-9
         ymax = max(y,ymax)
         ymin = min(y,ymin)
         z=z*1.d-9
         zmax = max(z,zmax)
         zmin = min(z,zmin)

         go to 3      
   2     close(1)         
      
C.....Calculate the grid and vectors
 
         domain_length = xmax-xmin 
         domain_width = ymax-ymin 
         domain_depth = zmax-zmin
         print*,' '
         print*,'Domain length (nm)=',domain_length*1.0e9
         print*,'Domain width (nm)=',domain_width*1.0e9
         print*,'Domain depth (nm)=',domain_depth*1.0e9
         print*,' '

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int(domain_depth/dz1 + 0.5)

         do i = 0, nx-1
            x_vec(i) = dx1
         enddo
         print *,'x-grid calculated, nx =',nx
     
         do j = 0, ny-1
            y_vec(j) = dy1
         enddo
         print *,'y-grid calculated, ny =',ny

         do k = 0, nz-1
            z_vec(k) = dz1
         enddo
         print *,'z-grid calculated, nz =',nz
         print*,' '
      
C........For ploting strain profiles for ALL atoms (cations and anions)
         open(unit=7,file="strain1Dz_txt"// trim(strain_file),
     1        status='unknown')
         open(unit=8,file="strain2Dxz_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=4,file="strain1Dxy_txt"//trim(strain_file),
     1        status='unknown')
         open(unit=9,file="strain1Dxmy_txt"// trim(strain_file),
     1       status='unknown')

         open(unit=1,file=strain_file,status='old')
         do ii=1,nn

          read(1,*)x,y,z,exx,eyy,ezz,exy,exz,eyz
           
          x=x*1e-9
          y=y*1e-9
          z=z*1e-9 
          if(int(x/dx1+0.5).eq.nx/2.and.
     &      int(y/dy1+0.5).eq.ny/2) 
     &     write(7,8)z,exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

          if(int(y/dy1+0.5).eq.ny/2)
     &     write(8,5)x,z,exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

          if((y/x.ge.0.98.and.y/x.le.1.02)
     &      .and.(int(z/dz1+0.5).eq.nz/2)) 
     &     write(4,5)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

          y=ymax-y       
          if((y/x.ge.0.98.and.y/x.le.1.02)
     &     .and.(int(z/dz1+0.5).eq.nz/2))          
     &    write(9,5)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

         enddo
         close(4)
         close(7)
         close(8)
         close(1)         
         close(9)

         print*,'reading strain3d file complete'

C........Reading the atom types        
         call getarg(3,atom_type_file)
         open(unit=1,file=atom_type_file,status='old')
         nn = 0
C.......(nn counts the cations only.
C....... atom() contains the cation atoms only.
C....... Currently polarization is calculated for Cations ONLY 
C....... Following the recipe of Zunger's paper/PRB 71, 045318, 2005)
   21    read(1,*,err=22,end=22)iatom
         if(iatom.eq.10.or.iatom.eq.14)then
          nn = nn+1 
          atom(nn) = iatom
         endif 
         go to 21
   22    close(1)
          
         print *,'reading atomtype file complete'
      
C........Here the polarization and charge calculation begin.... 
C........Reading the strain3d_cation file        
         call getarg(2, strain_file_cation)
         pol_const = -0.045D0   
C.......(deafult linear polz constant for InAs)

         nn = 0
         open(unit=1,file=strain_file_cation,status='old')
         open(unit=4,file="strain1Dxyc_txt"//trim(strain_file),
     1        status='unknown')
         open(unit=7,file="strain1Dzc_txt"// trim(strain_file),
     1        status='unknown')
         open(unit=8,file="strain2Dxzc_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=9,file="strain1Dxmyc_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=10,file="AtomGallium_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=11,file="AtomIndiam_txt"// trim(strain_file),
     1       status='unknown')
   20    read(1,*,err=9,end=9)x,y,z,exx,eyy,ezz,exy,exz,eyz
          nn=nn+1
          xion(nn)=x*1.d-9
          yion(nn)=y*1.d-9
          zion(nn)=z*1.d-9
          if(atom(nn).eq.14)then       
C         (InAs)
           if(piezo_model.eq.1)then
            pol_const_e14=-0.045
           elseif(piezo_model.eq.2.or.piezo_model.eq.3
     &           .or.piezo_model.eq.4)then
            pol_const_e14=-0.115

            pol_const_b114=-0.531  
            pol_const_b124=-4.076 
            pol_const_b156=-0.120
 
            pol_const_b141 = pol_const_b114 
            pol_const_b142 = pol_const_b124 
            pol_const_b165 = pol_const_b156 

            pol_const_b225 = pol_const_b114 
            pol_const_b336 = pol_const_b114 
            pol_const_b252 = pol_const_b225 
            pol_const_b363 = pol_const_b336 

            pol_const_b215 = pol_const_b124 
            pol_const_b316 = pol_const_b124 
            pol_const_b326 = pol_const_b124 
            pol_const_b235 = pol_const_b124 
            pol_const_b134 = pol_const_b124 
            pol_const_b251 = pol_const_b215 
            pol_const_b361 = pol_const_b316 
            pol_const_b362 = pol_const_b326 
            pol_const_b253 = pol_const_b235 
            pol_const_b143 = pol_const_b134 

            pol_const_b345 = pol_const_b156 
            pol_const_b354 = pol_const_b156 
            pol_const_b246 = pol_const_b156 
            pol_const_b264 = pol_const_b156 
           endif 
           flag(nn)=1          
           write(11,*)x,y,z         
          elseif(atom(nn).eq.10)then   
C         (GaAs)
           if(piezo_model.eq.1)then
            pol_const_e14=-0.16
           elseif(piezo_model.eq.2.or.piezo_model.eq.3
     &           .or.piezo_model.eq.4)then
            pol_const_e14=-0.230

            pol_const_b114=-0.439  
            pol_const_b124=-3.765  
            pol_const_b156=-0.492

            pol_const_b141 = pol_const_b114 
            pol_const_b142 = pol_const_b124 
            pol_const_b165 = pol_const_b156 

            pol_const_b225 = pol_const_b114 
            pol_const_b336 = pol_const_b114 
            pol_const_b252 = pol_const_b225 
            pol_const_b363 = pol_const_b336 

            pol_const_b215 = pol_const_b124 
            pol_const_b316 = pol_const_b124 
            pol_const_b326 = pol_const_b124 
            pol_const_b235 = pol_const_b124 
            pol_const_b134 = pol_const_b124 
            pol_const_b251 = pol_const_b215 
            pol_const_b361 = pol_const_b316 
            pol_const_b362 = pol_const_b326 
            pol_const_b253 = pol_const_b235 
            pol_const_b143 = pol_const_b134 

            pol_const_b345 = pol_const_b156 
            pol_const_b354 = pol_const_b156 
            pol_const_b246 = pol_const_b156 
            pol_const_b264 = pol_const_b156 
           endif 
           flag(nn)=0
           write(10,*)x,y,z         
          endif 
      
          if(piezo_model.eq.1)then
           px(nn)=pol_const_e14*2*eyz
           py(nn)=pol_const_e14*2*exz
           pz(nn)=pol_const_e14*2*exy
          elseif(piezo_model.eq.2)then
           px(nn)=pol_const_e14*2*eyz + 0.0*( 
     &       pol_const_b114*exx*eyz + pol_const_b141*eyz*exx + 
     &       pol_const_b124*eyy*eyz + pol_const_b142*eyz*eyy + 
     &       pol_const_b156*exz*exy + pol_const_b165*exy*exz + 
     &       pol_const_b134*ezz*eyz + pol_const_b143*eyz*ezz ) 

           py(nn)=pol_const_e14*2*exz + 0.0*( 
     &       pol_const_b215*exx*exz + pol_const_b251*exz*exx + 
     &       pol_const_b225*eyy*exz + pol_const_b252*exz*eyy + 
     &       pol_const_b235*ezz*exz + pol_const_b253*exz*ezz + 
     &       pol_const_b246*eyz*exy + pol_const_b264*exy*eyz ) 

           pz(nn)=pol_const_e14*2*exy + 0.0*( 
     &       pol_const_b316*exx*exy + pol_const_b361*exy*exx + 
     &       pol_const_b326*eyy*exy + pol_const_b362*exy*eyy + 
     &       pol_const_b336*ezz*exy + pol_const_b363*exy*ezz + 
     &       pol_const_b345*eyz*exz + pol_const_b354*exz*eyz ) 
          elseif(piezo_model.eq.3)then
           px(nn)=pol_const_e14*0*eyz + 0.5*( 
     &       pol_const_b114*exx*eyz + pol_const_b141*eyz*exx + 
     &       pol_const_b124*eyy*eyz + pol_const_b142*eyz*eyy + 
     &       pol_const_b156*exz*exy + pol_const_b165*exy*exz + 
     &       pol_const_b134*ezz*eyz + pol_const_b143*eyz*ezz ) 

           py(nn)=pol_const_e14*0*exz + 0.5*( 
     &       pol_const_b215*exx*exz + pol_const_b251*exz*exx + 
     &       pol_const_b225*eyy*exz + pol_const_b252*exz*eyy + 
     &       pol_const_b235*ezz*exz + pol_const_b253*exz*ezz + 
     &       pol_const_b246*eyz*exy + pol_const_b264*exy*eyz ) 

           pz(nn)=pol_const_e14*0*exy + 0.5*( 
     &       pol_const_b316*exx*exy + pol_const_b361*exy*exx + 
     &       pol_const_b326*eyy*exy + pol_const_b362*exy*eyy + 
     &       pol_const_b336*ezz*exy + pol_const_b363*exy*ezz + 
     &       pol_const_b345*eyz*exz + pol_const_b354*exz*eyz ) 
          elseif(piezo_model.eq.4)then
           px(nn)=pol_const_e14*2*eyz + 0.5*( 
     &       pol_const_b114*exx*eyz + pol_const_b141*eyz*exx + 
     &       pol_const_b124*eyy*eyz + pol_const_b142*eyz*eyy + 
     &       pol_const_b156*exz*exy + pol_const_b165*exy*exz + 
     &       pol_const_b134*ezz*eyz + pol_const_b143*eyz*ezz ) 

           py(nn)=pol_const_e14*2*exz + 0.5*( 
     &       pol_const_b215*exx*exz + pol_const_b251*exz*exx + 
     &       pol_const_b225*eyy*exz + pol_const_b252*exz*eyy + 
     &       pol_const_b235*ezz*exz + pol_const_b253*exz*ezz + 
     &       pol_const_b246*eyz*exy + pol_const_b264*exy*eyz ) 

           pz(nn)=pol_const_e14*2*exy + 0.5*( 
     &       pol_const_b316*exx*exy + pol_const_b361*exy*exx + 
     &       pol_const_b326*eyy*exy + pol_const_b362*exy*eyy + 
     &       pol_const_b336*ezz*exy + pol_const_b363*exy*ezz + 
     &       pol_const_b345*eyz*exz + pol_const_b354*exz*eyz ) 
          endif

          if((y/x.ge.0.98.and.y/x.le.1.02)
     &   .and.(int(zion(nn)/dz1+0.5).eq.nz/2)) 
     &    write(4,51)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz
     &              ,px(nn),py(nn),pz(nn)
          if(int(xion(nn)/dx1+0.5).eq.nx/2.and.
     &     int(yion(nn)/dy1+0.5).eq.ny/2) 
     &     write(7,81)z,exx,eyy,ezz,exy,exz,eyz,(exx+eyy+ezz)
     &              ,px(nn),py(nn),pz(nn)
          if(int(yion(nn)/dy1+0.5).eq.ny/2)
     &    write(8,51)x,z,exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz
     &              ,px(nn),py(nn),pz(nn)
          y=ymax*1.e9 - y       
          if((y/x.ge.0.98.and.y/x.le.1.02)
     &     .and.(int(z*1.e-9/dz1+0.5).eq.nz/2))          
     &    write(9,51)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz
     &              ,px(nn),py(nn),pz(nn)

         go to 20      
   9     close(1)         
         close(4)
         close(7)
         close(8)
         close(9)
         close(10)
         close(11)

         print *,'reading strain3d_cation file complete'

C.......Finding polarization for each GRID point
C.......Independent of polarization model (first/higher order)
C......(A) Initialize
        do i = 0,nx
        do j = 0,ny
        do k = 0,nz
          polx(i,j,k) = 0.         
          poly(i,j,k) = 0.   
          polz(i,j,k) = 0.   
        enddo
        enddo
        enddo
        
C......(B) Interpolate (cations only)
        do ii=1,nn
         i=int(xion(ii)/dx1)
         j=int(yion(ii)/dy1)
         k=int(zion(ii)/dz1)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1

         flag_reg(i,j,k) = flag(ii)
C       (flag_reg 0 for GaAs and 1 for InAs)

        if(interpolation.eq.0)then 
         polx(i,j,k) = polx(i,j,k) + px(ii)
         poly(i,j,k) = poly(i,j,k) + py(ii)
         polz(i,j,k) = polz(i,j,k) + pz(ii)
        endif
  
        if(interpolation.eq.1)then 
         polx(i,j,k)       = polx(i,j,k) + px(ii)*0.125
         polx(i,j+1,k)     = polx(i,j+1,k) + px(ii)*0.125
         polx(i+1,j,k)     = polx(i+1,j,k) + px(ii)*0.125
         polx(i+1,j+1,k)   = polx(i+1,j+1,k) + px(ii)*0.125
         polx(i,j,k+1)     = polx(i,j,k+1) + px(ii)*0.125
         polx(i,j+1,k+1)   = polx(i,j+1,k+1) + px(ii)*0.125
         polx(i+1,j,k+1)   = polx(i+1,j,k+1) + px(ii)*0.125
         polx(i+1,j+1,k+1) = polx(i+1,j+1,k+1) + px(ii)*0.125
         
         poly(i,j,k)       = poly(i,j,k) + py(ii)*0.125
         poly(i,j+1,k)     = poly(i,j+1,k) + py(ii)*0.125
         poly(i+1,j,k)     = poly(i+1,j,k) + py(ii)*0.125
         poly(i+1,j+1,k)   = poly(i+1,j+1,k) + py(ii)*0.125
         poly(i,j,k+1)     = poly(i,j,k+1) + py(ii)*0.125
         poly(i,j+1,k+1)   = poly(i,j+1,k+1) + py(ii)*0.125
         poly(i+1,j,k+1)   = poly(i+1,j,k+1) + py(ii)*0.125
         poly(i+1,j+1,k+1) = poly(i+1,j+1,k+1) + py(ii)*0.125
      
         polz(i,j,k)       = polz(i,j,k) + pz(ii)*0.125
         polz(i,j+1,k)     = polz(i,j+1,k) + pz(ii)*0.125
         polz(i+1,j,k)     = polz(i+1,j,k) + pz(ii)*0.125
         polz(i+1,j+1,k)   = polz(i+1,j+1,k) + pz(ii)*0.125
         polz(i,j,k+1)     = polz(i,j,k+1) + pz(ii)*0.125
         polz(i,j+1,k+1)   = polz(i,j+1,k+1) + pz(ii)*0.125
         polz(i+1,j,k+1)   = polz(i+1,j,k+1) + pz(ii)*0.125
         polz(i+1,j+1,k+1) = polz(i+1,j+1,k+1) + pz(ii)*0.125
        endif

        if(interpolation.eq.2)then
         x0 = (xion(ii) - i*dx1)/dx1
         x1 = 1.0 - x0
         y0 = (yion(ii) - j*dy1)/dy1
         y1 = 1.0 - y0
         z0 = (zion(ii) - k*dz1)/dz1
         z1 = 1.0 - z0
 
         polx(i,j,k)       = polx(i,j,k) + px(ii)*x0*y0*z0
         polx(i,j+1,k)     = polx(i,j+1,k) + px(ii)*x0*y1*z0
         polx(i+1,j,k)     = polx(i+1,j,k) + px(ii)*x1*y0*z0
         polx(i+1,j+1,k)   = polx(i+1,j+1,k) + px(ii)*x1*y1*z0
         polx(i,j,k+1)     = polx(i,j,k+1) + px(ii)*x0*y0*z1
         polx(i,j+1,k+1)   = polx(i,j+1,k+1) + px(ii)*x0*y1*z1
         polx(i+1,j,k+1)   = polx(i+1,j,k+1) + px(ii)*x1*y0*z1
         polx(i+1,j+1,k+1) = polx(i+1,j+1,k+1) + px(ii)*x1*y1*z1
         
         poly(i,j,k)       = poly(i,j,k) + py(ii)*x0*y0*z0
         poly(i,j+1,k)     = poly(i,j+1,k) + py(ii)*x0*y1*z0
         poly(i+1,j,k)     = poly(i+1,j,k) + py(ii)*x1*y0*z0
         poly(i+1,j+1,k)   = poly(i+1,j+1,k) + py(ii)*x1*y1*z0
         poly(i,j,k+1)     = poly(i,j,k+1) + py(ii)*x0*y0*z1
         poly(i,j+1,k+1)   = poly(i,j+1,k+1) + py(ii)*x0*y1*z1
         poly(i+1,j,k+1)   = poly(i+1,j,k+1) + py(ii)*x1*y0*z1
         poly(i+1,j+1,k+1) = poly(i+1,j+1,k+1) + py(ii)*x1*y1*z1
      
         polz(i,j,k)       = polz(i,j,k) + pz(ii)*x0*y0*z0
         polz(i,j+1,k)     = polz(i,j+1,k) + pz(ii)*x0*y1*z0
         polz(i+1,j,k)     = polz(i+1,j,k) + pz(ii)*x1*y0*z0
         polz(i+1,j+1,k)   = polz(i+1,j+1,k) + pz(ii)*x1*y1*z0
         polz(i,j,k+1)     = polz(i,j,k+1) + pz(ii)*x0*y0*z1
         polz(i,j+1,k+1)   = polz(i,j+1,k+1) + pz(ii)*x0*y1*z1
         polz(i+1,j,k+1)   = polz(i+1,j,k+1) + pz(ii)*x1*y0*z1
         polz(i+1,j+1,k+1) = polz(i+1,j+1,k+1) + pz(ii)*x1*y1*z1
        endif

        enddo
        close(1)         
        close(9)

      endif

C.....Calculate charge density for Poisson equation
 
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      open(unit=133,file='ChargeDensity2Dxz_txt'//trim(strain_file),
     &     status='unknown')
      open(unit=134,file='ChargeDensity2Dxy_txt'//trim(strain_file),
     &     status='unknown')
      
      do i = 0,nx
      do j = 0,ny
      do k = 0,nz
       if(istrainfile.eq.1)then       
C         dpx = (polx(i+1,j,k) - polx(i,j,k))/dx1 !/(4./discret_factor_x) 
C         dpy = (poly(i,j+1,k) - poly(i,j,k))/dy1 !/(4./discret_factor_y)
C         dpz = (polz(i,j,k+1) - polz(i,j,k))/dz1 !/(4./discret_factor_z) 
C       (Differentiation on grids (numerical appx formula used))
         dpx = (1./6./dx1) * (11*polx(i,j,k) - 18*polx(i-1,j,k) 
     &                        + 9*polx(i-2,j,k) - 2*polx(i-3,j,k))  
         dpy = (1./6./dy1) * (11*poly(i,j,k) - 18*poly(i,j-1,k) 
     &                        + 9*poly(i,j-2,k) - 2*poly(i,j-3,k))  
         dpz = (1./6./dz1) * (11*polz(i,j,k) - 18*polz(i,j,k-1) 
     &                        + 9*polz(i,j,k-2) - 2*polz(i,j,k-3))  

         rion_number(i,j,k) = -(dpx+dpy+dpz) 
         denn = rion_number(i,j,k)
       endif

C      'dop()' is specifically used for Poisson solver
       if(denn.gt.0.0)then
         dop(i,j,k) = -denn - 0.0e0   
       else
         dop(i,j,k) = -denn + 0.0e0
       endif
      
       if(j.eq.ny/2) write(133,*)i,k,-dop(i,j,k)
       if(k.eq.k_fix) write(134,*)i,j,-dop(i,j,k)

       flag_dom(i,j,k) = 4
       n(i,j,k) = 0.D-5
       p(i,j,k) = 0.D-5
       ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
       if(flag_reg(i,j,k).eq.1)eps(i,j,k) = eps_ingaas
       if(flag_reg(i,j,k).eq.0)eps(i,j,k) = eps_gaas
      enddo
      enddo
      enddo
      close(133)
      close(134)
     



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....CQD with Gate Oxide
C.....InAs/GaAs Coupled Quantum Dot for Piezoelectric calculation
C.....There is gate oxide layer for applying external potential
      case('CQD')
         
      dx1 = a0/discret_factor_x
      dy1 = a0/discret_factor_y
      dz1 = a0/discret_factor_z

      istrainfile = 1       

      if(istrainfile.eq.1)then       
         argc = iargc()
         if(argc.ne.3) then
          write(*,*)'Usage: ../bin/poisson3d.ex <strain3d_file> 
     &               <strain3d_file_cation> <atom_file>'
          write(*,*)' '
          call exit(1)
         endif

C........Reading the strain3d file for all atoms
C........We want to find the domain size and number of total atoms
         call getarg(1, strain_file)

         xmax = -1.d38
         xmin = 1.d38
         ymax = -1.d38
         ymin = 1.d38
         zmax = -1.d38
         zmin = 1.d38

         nn = 0
         
         open(unit=1,file=strain_file,status='old')
   7     read(1,*,err=18,end=18)x,y,z,exx,eyy,ezz,exy,exz,eyz
         nn=nn+1
         x=x*1.d-9
         xmax = dmax1(x,xmax)
         xmin = dmin1(x,xmin)
         y=y*1.d-9
         ymax = max(y,ymax)
         ymin = min(y,ymin)
         z=z*1.d-9
         zmax = max(z,zmax)
         zmin = min(z,zmin)

         go to 7      
   18    close(1)         
      
C.....Calculate the grid and vectors
 
         domain_length = xmax-xmin 
         domain_width = ymax-ymin 
         domain_depth = zmax-zmin + z_ox
         print*,' '
         print*,'Domain length (nm) =',domain_length*1.0e9
         print*,'Domain width (nm) =',domain_width*1.0e9
         print*,'Domain depth (nm) =',domain_depth*1.0e9

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int((zmax-zmin)/dz1 + 0.5)

         do i = 0, nx-1
            x_vec(i) = dx1
         enddo
         print *,'x-grid calculated, nx =',nx
     
         do j = 0, ny-1
            y_vec(j) = dy1
         enddo
         print *,'y-grid calculated, ny =',ny

         if(z_ox.ne.0)then
          zo = 0.1e-9             
C        (initial grid step)
          r_in = 1.1D0
          r = r_coeff(r_in,zo,k_ox,z_ox)
          z_vec(-1) = zo
          do k = -2, -k_ox, -1
            z_vec(k) = r*z_vec(k+1)
          enddo
         endif

         do k = 0, nz-1
            z_vec(k) = dz1
         enddo
         print *,'z-grid calculated, nz =',nz
      
C........For ploting strain profiles for ALL atoms (cations and anions)
         open(unit=7,file="strain1Dz_txt"// trim(strain_file),
     1        status='unknown')
         open(unit=8,file="strain2Dxz_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=4,file="strain1Dxy_txt"//trim(strain_file),
     1        status='unknown')
         open(unit=9,file="strain1Dxmy_txt"// trim(strain_file),
     1       status='unknown')

         open(unit=1,file=strain_file,status='old')
         do ii=1,nn

          read(1,*)x,y,z,exx,eyy,ezz,exy,exz,eyz
           
          x=x*1e-9
          y=y*1e-9
          z=z*1e-9 
          if(int(x/dx1+0.5).eq.nx/2.and.
     &      int(y/dy1+0.5).eq.ny/2) 
     &     write(7,8)z,exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

          if(int(y/dy1+0.5).eq.ny/2)
     &     write(8,5)x,z,exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

          if((y/x.ge.0.98.and.y/x.le.1.02)
     &      .and.(int(z/dz1+0.5).eq.nz/2)) 
     &     write(4,5)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

          y=ymax-y       
          if((y/x.ge.0.98.and.y/x.le.1.02)
     &     .and.(int(z/dz1+0.5).eq.nz/2))          
     &    write(9,5)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz

         enddo
         close(4)
         close(7)
         close(8)
         close(1)         
         close(9)

         print *,'reading strain3d file complete'

C........Reading the atom types        
         call getarg(3,atom_type_file)
         open(unit=1,file=atom_type_file,status='old')
         nn = 0
C.......(nn counts the cations only.
C....... atom() contains the cation atoms only.
C....... Currently polarization is calculated for Cations ONLY 
C....... Following the recipe of Zunger's paper/PRB 71, 045318, 2005)
   31    read(1,*,err=32,end=32)iatom
         if(iatom.eq.10.or.iatom.eq.14)then  !!cation
          nn = nn+1 
          atom(nn) = iatom
         endif 
         go to 31
   32    close(1)
          
         print *,'reading atomtype file complete'
       
 
C........Here the polarization and charge calculation begin.... 
C........Reading the strain3d_cation file        
         call getarg(2, strain_file_cation)
         pol_const = -0.045D0  
C........Default linear value for InAs

         nn = 0
         open(unit=1,file=strain_file_cation,status='old')
         open(unit=4,file="strain1Dxyc_txt"//trim(strain_file),
     1        status='unknown')
         open(unit=7,file="strain1Dzc_txt"// trim(strain_file),
     1        status='unknown')
         open(unit=8,file="strain2Dxzc_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=9,file="strain1Dxmyc_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=10,file="AtomGallium_txt"// trim(strain_file),
     1       status='unknown')
         open(unit=11,file="AtomIndiam_txt"// trim(strain_file),
     1       status='unknown')
   30    read(1,*,err=10,end=10)x,y,z,exx,eyy,ezz,exy,exz,eyz
          nn=nn+1
          xion(nn)=x*1.d-9
          yion(nn)=y*1.d-9
          zion(nn)=z*1.d-9
          if(atom(nn).eq.14)then
C          (InAs)
            pol_const=-0.045D0  
            flag(nn)=1
            write(11,*)x,y,z         
          elseif(atom(nn).eq.10)then
C          (GaAs)
            pol_const=-0.16D0   
            flag(nn)=0
            write(10,*)x,y,z         
          endif 
      
          px(nn)=pol_const*2.*eyz
          py(nn)=pol_const*2.*exz
          pz(nn)=pol_const*2.*exy
          if((y/x.ge.0.98.and.y/x.le.1.02)
     &   .and.(int(zion(nn)/dz1+0.5).eq.nz/2)) 
     &    write(4,51)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz
     &              ,px(nn),py(nn),pz(nn)
          if(int(xion(nn)/dx1+0.5).eq.nx/2.and.
     &     int(yion(nn)/dy1+0.5).eq.ny/2) 
     &     write(7,81)z,exx,eyy,ezz,exy,exz,eyz,(exx+eyy+ezz)
     &              ,px(nn),py(nn),pz(nn)
          if(int(yion(nn)/dy1+0.5).eq.ny/2)
     &    write(8,51)x,z,exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz
     &              ,px(nn),py(nn),pz(nn)
          y=ymax*1.e9 - y       
          if((y/x.ge.0.98.and.y/x.le.1.02)
     &     .and.(int(z*1.e-9/dz1+0.5).eq.nz/2))          
     &    write(9,51)x,sqrt(x*x+y*y),exx,eyy,ezz,exy,exz,eyz,exx+eyy+ezz
     &              ,px(nn),py(nn),pz(nn)

         go to 30      
   10    close(1)         
         close(4)
         close(7)
         close(8)
         close(9)
         close(10)
         close(11)

         print *,'reading strain3d_cation file complete'

C.......Finding polarization for each GRID point
C.......Independent of polarization model (first/higher order)
C......(A) Initialize
        do i = 0,nx
        do j = 0,ny
        do k = 0,nz
          polx(i,j,k) = 0.         
          poly(i,j,k) = 0.   
          polz(i,j,k) = 0.   
        enddo
        enddo
        enddo
        
C......(B) Interpolate (cations only)
        do ii=1,nn
         i=int(xion(ii)/dx1)
         j=int(yion(ii)/dy1)
         k=int(zion(ii)/dz1)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
C        (flag_reg 0 for GaAs and 1 for InAs)
         flag_reg(i,j,k) = flag(ii)

        if(interpolation.eq.0)then 
         polx(i,j,k) = polx(i,j,k) + px(ii)
         poly(i,j,k) = poly(i,j,k) + py(ii)
         polz(i,j,k) = polz(i,j,k) + pz(ii)
        endif
  
        if(interpolation.eq.1)then 
         polx(i,j,k)       = polx(i,j,k) + px(ii)*0.125
         polx(i,j+1,k)     = polx(i,j+1,k) + px(ii)*0.125
         polx(i+1,j,k)     = polx(i+1,j,k) + px(ii)*0.125
         polx(i+1,j+1,k)   = polx(i+1,j+1,k) + px(ii)*0.125
         polx(i,j,k+1)     = polx(i,j,k+1) + px(ii)*0.125
         polx(i,j+1,k+1)   = polx(i,j+1,k+1) + px(ii)*0.125
         polx(i+1,j,k+1)   = polx(i+1,j,k+1) + px(ii)*0.125
         polx(i+1,j+1,k+1) = polx(i+1,j+1,k+1) + px(ii)*0.125
         
         poly(i,j,k)       = poly(i,j,k) + py(ii)*0.125
         poly(i,j+1,k)     = poly(i,j+1,k) + py(ii)*0.125
         poly(i+1,j,k)     = poly(i+1,j,k) + py(ii)*0.125
         poly(i+1,j+1,k)   = poly(i+1,j+1,k) + py(ii)*0.125
         poly(i,j,k+1)     = poly(i,j,k+1) + py(ii)*0.125
         poly(i,j+1,k+1)   = poly(i,j+1,k+1) + py(ii)*0.125
         poly(i+1,j,k+1)   = poly(i+1,j,k+1) + py(ii)*0.125
         poly(i+1,j+1,k+1) = poly(i+1,j+1,k+1) + py(ii)*0.125
      
         polz(i,j,k)       = polz(i,j,k) + pz(ii)*0.125
         polz(i,j+1,k)     = polz(i,j+1,k) + pz(ii)*0.125
         polz(i+1,j,k)     = polz(i+1,j,k) + pz(ii)*0.125
         polz(i+1,j+1,k)   = polz(i+1,j+1,k) + pz(ii)*0.125
         polz(i,j,k+1)     = polz(i,j,k+1) + pz(ii)*0.125
         polz(i,j+1,k+1)   = polz(i,j+1,k+1) + pz(ii)*0.125
         polz(i+1,j,k+1)   = polz(i+1,j,k+1) + pz(ii)*0.125
         polz(i+1,j+1,k+1) = polz(i+1,j+1,k+1) + pz(ii)*0.125
        endif

        if(interpolation.eq.2)then
         x0 = (xion(ii) - i*dx1)/dx1
         x1 = 1.0 - x0
         y0 = (yion(ii) - j*dy1)/dy1
         y1 = 1.0 - y0
         z0 = (zion(ii) - k*dz1)/dz1
         z1 = 1.0 - z0
 
         polx(i,j,k)       = polx(i,j,k) + px(ii)*x0*y0*z0
         polx(i,j+1,k)     = polx(i,j+1,k) + px(ii)*x0*y1*z0
         polx(i+1,j,k)     = polx(i+1,j,k) + px(ii)*x1*y0*z0
         polx(i+1,j+1,k)   = polx(i+1,j+1,k) + px(ii)*x1*y1*z0
         polx(i,j,k+1)     = polx(i,j,k+1) + px(ii)*x0*y0*z1
         polx(i,j+1,k+1)   = polx(i,j+1,k+1) + px(ii)*x0*y1*z1
         polx(i+1,j,k+1)   = polx(i+1,j,k+1) + px(ii)*x1*y0*z1
         polx(i+1,j+1,k+1) = polx(i+1,j+1,k+1) + px(ii)*x1*y1*z1
         
         poly(i,j,k)       = poly(i,j,k) + py(ii)*x0*y0*z0
         poly(i,j+1,k)     = poly(i,j+1,k) + py(ii)*x0*y1*z0
         poly(i+1,j,k)     = poly(i+1,j,k) + py(ii)*x1*y0*z0
         poly(i+1,j+1,k)   = poly(i+1,j+1,k) + py(ii)*x1*y1*z0
         poly(i,j,k+1)     = poly(i,j,k+1) + py(ii)*x0*y0*z1
         poly(i,j+1,k+1)   = poly(i,j+1,k+1) + py(ii)*x0*y1*z1
         poly(i+1,j,k+1)   = poly(i+1,j,k+1) + py(ii)*x1*y0*z1
         poly(i+1,j+1,k+1) = poly(i+1,j+1,k+1) + py(ii)*x1*y1*z1
      
         polz(i,j,k)       = polz(i,j,k) + pz(ii)*x0*y0*z0
         polz(i,j+1,k)     = polz(i,j+1,k) + pz(ii)*x0*y1*z0
         polz(i+1,j,k)     = polz(i+1,j,k) + pz(ii)*x1*y0*z0
         polz(i+1,j+1,k)   = polz(i+1,j+1,k) + pz(ii)*x1*y1*z0
         polz(i,j,k+1)     = polz(i,j,k+1) + pz(ii)*x0*y0*z1
         polz(i,j+1,k+1)   = polz(i,j+1,k+1) + pz(ii)*x0*y1*z1
         polz(i+1,j,k+1)   = polz(i+1,j,k+1) + pz(ii)*x1*y0*z1
         polz(i+1,j+1,k+1) = polz(i+1,j+1,k+1) + pz(ii)*x1*y1*z1
        endif

        enddo
        close(1)         
        close(9)

      endif

      if(istrainfile.ne.1)then       
         argc = iargc()
         if(argc.ne.2) then
          write(*,*)'Use: ../bin/poisson3d.ex <rAtom_file> <aType_file>'
          write(*,*)' '
          call exit(1)
         endif


         call getarg(1, ratom_file)

         xmax = -1.d38
         xmin = 1.d38
         ymax = -1.d38
         ymin = 1.d38
         zmax = -1.d38
         zmin = 1.d38

         nn = 0
         
         open(unit=1,file=ratom_file,status='old')
   11    read(1,*,err=12,end=12)x,y,z
         nn=nn+1
         x=x*1.d-9
         xion(nn) = x
         xmax = dmax1(x,xmax)
         xmin = dmin1(x,xmin)
         y=y*1.d-9
         yion(nn) = y
         ymax = max(y,ymax)
         ymin = min(y,ymin)
         z=z*1.d-9
         zion(nn) = z
         zmax = max(z,zmax)
         zmin = min(z,zmin)

         go to 11      
   12    close(1)         
         print *,'reading ratom (position) file complete'
      
C.....Calculate the grid and vectors
    
         domain_length = xmax-xmin 
         domain_width = ymax-ymin 
         domain_depth = zmax-zmin + z_ox
         print*,' '
         print*,'Domain length =',domain_length
         print*,'Domain width =',domain_width
         print*,'Domain depth =',domain_depth

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int((zmax-zmin)/dz1 + 0.5)

         do i = 0, nx-1
            x_vec(i) = dx1
         enddo
         print *,'x-grid calculated, nx =',nx
     
         do j = 0, ny-1
            y_vec(j) = dy1
         enddo
         print *,'y-grid calculated, ny =',ny

         if(z_ox.ne.0)then
         zo = 0.1e-9             
         r_in = 1.1D0
         r = r_coeff(r_in,zo,k_ox,z_ox)
         z_vec(-1) = zo
         do k = -2, -k_ox, -1
            z_vec(k) = r*z_vec(k+1)
         enddo
         endif

         do k = 0, nz-1
            z_vec(k) = dz1
         enddo
         print *,'z-grid calculated, nz =',nz
 
     
C........Reading the atom types        
         call getarg(2,atom_type_file)
         open(unit=1,file=atom_type_file,status='old')
         nn = 0
   14    read(1,*,err=13,end=13)iatom
         if(iatom.eq.10.or.iatom.eq.14)then  !!cation
          nn = nn+1 
          atom(nn) = iatom
         endif 
         go to 14
   13    close(1)
          
         print *,'reading atomtype file complete'

      endif
C     end the non-strain-file condition)


 
C.....Calculate charge density
 
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      open(unit=133,file='ChargeDensity2Dxz_txt'//trim(strain_file),
     &     status='unknown')
      open(unit=134,file='ChargeDensity2Dxy_txt'//trim(strain_file),
     &     status='unknown')
      
      do i = 0,nx
      do j = 0,ny
      do k = 0,nz
       if(istrainfile.eq.1)then       
         dpx = (1./6./dx1) * (11*polx(i,j,k) - 18*polx(i-1,j,k) 
     &                        + 9*polx(i-2,j,k) - 2*polx(i-3,j,k))  
         dpy = (1./6./dy1) * (11*poly(i,j,k) - 18*poly(i,j-1,k) 
     &                        + 9*poly(i,j-2,k) - 2*poly(i,j-3,k))  
         dpz = (1./6./dz1) * (11*polz(i,j,k) - 18*polz(i,j,k-1) 
     &                        + 9*polz(i,j,k-2) - 2*polz(i,j,k-3))  
         rion_number(i,j,k) = - (dpx+dpy+dpz) 
         denn = rion_number(i,j,k)
       endif
       if(istrainfile.ne.1)then       
          rion_number(i,j,k) = 0.0
          denn = rion_number(i,j,k)
       endif

       if(denn.gt.0.0)then
            dop(i,j,k) = -denn - 0.0e0   
       else
            dop(i,j,k) = -denn + 0.0e0
       endif
       if(j.eq.ny/2) write(133,*)i,k,-dop(i,j,k)
       if(k.eq.k_fix) write(134,*)i,j,-dop(i,j,k)

       flag_dom(i,j,k) = 4
       n(i,j,k) = 0.D-5
       p(i,j,k) = 0.D-5
       ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
       if(flag_reg(i,j,k).eq.1)eps(i,j,k) = eps_ingaas
       if(flag_reg(i,j,k).eq.0)eps(i,j,k) = eps_gaas
      enddo
      enddo
      enddo
      close(133)
      close(134)

C     Gate oxide reion 
      if(z_ox.ne.0)then

       do i = 0,nx
       do j = 0,ny
       do k = -k_ox,0
        dop(i,j,k) = 0.0e0   
        flag_dom(i,j,k) = 1
        n(i,j,k) = 0.D0
        p(i,j,k) = 0.D0
        eps(i,j,k) = eps_ox
        ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
       enddo
       enddo
       enddo
     
C......boundary condition 
       k = -k_ox
       DO i = 0, nx
       DO j = 0, ny
         flag_dom(i,j,k) = 5
         fi(i,j,k) = gate_voltage
         ro(i,j,k) = fi(i,j,k)
       ENDDO
       ENDDO

      endif




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....P/Si system
C.....Phosphorous impurity and Silicon based quantum computing 

      case('SiP')
         
         argc = iargc()
         if(argc.ne.1) then
          write(*,*)'Use: ../bin/poisson3d.ex <rAtom_file>'
          write(*,*)' '
          call exit(1)
         endif

         call getarg(1, ratom_file)

         xmax = -1.d38
         xmin = 1.d38
         ymax = -1.d38
         ymin = 1.d38
         zmax = -1.d38
         zmin = 1.d38

         nn = 0
         
         open(unit=1,file=ratom_file,status='old')
   41    read(1,*,err=42,end=42)x,y,z
         nn=nn+1
         x=x*1.d-9
         xion(nn) = x
         xmax = dmax1(x,xmax)
         xmin = dmin1(x,xmin)
         y=y*1.d-9
         yion(nn) = y
         ymax = max(y,ymax)
         ymin = min(y,ymin)
         z=z*1.d-9
         zion(nn) = z
         zmax = max(z,zmax)
         zmin = min(z,zmin)
      
C       (for Silicon atoms...)
         flag(nn) = 0
         qion(nn)= 0     

         go to 41      
   42    close(1)         
         print *,'reading ratom (position) file complete'
      
C.....Calculate the grid and vectors
    
         domain_length = xmax-xmin 
         domain_width = ymax-ymin 
         domain_depth = zmax-zmin + z_ox
         print*,' '
         print*,'Domain length (nm) =',domain_length*1.0e9
         print*,'Domain width (nm) =',domain_width*1.0e9
         print*,'Domain depth (nm) =',domain_depth*1.0e9

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int((zmax-zmin)/dz1 + 0.5)

         do i = 0, nx-1
            x_vec(i) = dx1
         enddo
         print *,'x-grid calculated, nx =',nx
     
         do j = 0, ny-1
            y_vec(j) = dy1
         enddo
         print *,'y-grid calculated, ny =',ny

         if(z_ox.ne.0)then
         zo = 0.1e-9             
C        (initial grid step)
         r_in = 1.1D0
         r = r_coeff(r_in,zo,k_ox,z_ox)
         z_vec(-1) = zo
         do k = -2, -k_ox, -1
            z_vec(k) = r*z_vec(k+1)
         enddo
         endif

         do k = 0, nz-1
            z_vec(k) = dz1
         enddo
         print *,'z-grid calculated, nz =',nz
 
        norm_doping = 1.0
        cell_volume = dx1*dy1*dz1
        norm_factor = cell_volume*norm_doping
        norm_factor = 1./norm_factor
        open(unit=52,file='Phosphorous_dopants_txt'//trim(extension),
     1                status='unknown')

        do i = 1,n_dopants
           nn=nn+1
           xion(nn)=x_dop(i)
           yion(nn)=y_dop(i)
           zion(nn)=z_dop(i)
           qion(nn)=charge_dop(i)*norm_factor*q      
           write(52,91)x_dop(i),y_dop(i),z_dop(i),
     &             charge_dop(i)*norm_factor*q
C         (for Phosphorous atoms)
           flag(nn) = 1    
        enddo            
  91    format(4(1x,e12.6))
        close(52)

        print*,'Total atom (+ Phosphorous) number = ',nn
      
        do i = 0,nx
        do j = 0,ny
        do k = 0,nz
           rion_number(i,j,k) = 0.   
C         (charge density)
        enddo
        enddo
        enddo
        
      do ii=1,nn
         i=int(xion(ii)/dx1 + 0.5)
         j=int(yion(ii)/dy1 + 0.5)
         k=int(zion(ii)/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
         rion_number(i,j,k)=rion_number(i,j,k) + qion(ii)  
         flag_reg(i,j,k) = flag(ii) 
      enddo
   
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      open(unit=233,file='ChargeDensity2Dxz_txt'//trim(extension),
     1           status='unknown')
      open(unit=232,file='ChargeDensity2Dxy_txt'//trim(extension),
     1           status='unknown')
      do i = 0,nx
      do j = 0,ny
      do k = 0,nz
         denn = rion_number(i,j,k)
         if(denn.gt.0.0)then
            dop(i,j,k) = -denn - 0.0e0   
         else
            dop(i,j,k) = -denn + 0.0e0
         endif
        if(k.eq.k_fix)write(232,*)i,j,-dop(i,j,k)
        if(j.eq.ny/2)write(233,*)i,k,-dop(i,j,k)
           flag_dom(i,j,k) = 4
           n(i,j,k) = 0.D-5
           p(i,j,k) = 0.D-5
           ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
           if(flag_reg(i,j,k).eq.1)eps(i,j,k) = eps_p
           if(flag_reg(i,j,k).eq.0)eps(i,j,k) = eps_si
      enddo
      enddo
      enddo
      close(232)
      close(233)

C     Oxide region
      do i = 0,nx
      do j = 0,ny
      do k = -k_ox,0
        dop(i,j,k) = 0.0e0   
        flag_dom(i,j,k) = 1
        n(i,j,k) = 0.D0
        p(i,j,k) = 0.D0
        eps(i,j,k) = eps_ox
        ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
      enddo
      enddo
      enddo
     
C.....boundary condition 
      k = -k_ox
      DO i = 0, nx
      DO j = 0, ny
         flag_dom(i,j,k) = 5
         fi(i,j,k) = gate_voltage
         ro(i,j,k) = fi(i,j,k)
      ENDDO
      ENDDO


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....SiP_DEF
C.....Here we are creating the Zincblende crystal ourselves
C.....Also, we are assigning arbitrary charges
      case('SiP_DEF')
         
         argc = iargc()
         if(argc.ne.0) then
          write(*,*)'Usage: ../bin/poisson3d.ex'
          write(*,*)' '
          call exit(1)
         endif
         
         dx1 = a0/discret_factor_x
         dy1 = a0/discret_factor_y
         dz1 = a0/discret_factor_z

         domain_length = x_length
         domain_width = y_width
         domain_depth = z_ox+z_depth

         print*,' '
         print*,'domain length (nm) =',domain_length*1.0e9
         print*,'domain width (nm) =',domain_width*1.0e9
         print*,'domain depth (nm) =',domain_depth*1.0e9

         if(z_ox.gt.0)then      
           zo = 0.1e-9             
C         (initial grid step)
           r_in = 1.1D0
           r = r_coeff(r_in,zo,k_ox,z_ox)
           z_vec(-1) = zo
           do k = -2, -k_ox, -1
              z_vec(k) = r*z_vec(k+1)/debye
           enddo
         endif

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int(z_depth/dz1 + 0.5) 

         do i = 0, nx-1
          x_vec(i) = dx1/debye
         enddo
         print *,'x-grid calculated, nx =',nx
     
         do j = 0, ny-1
          y_vec(j) = dy1/debye
         enddo
         print *,'y-grid calculated, ny =',ny

         do k = 0, nz-1
          z_vec(k) = dz1/debye
         enddo
         print *,'z-grid calculated, nz = ',nz

         norm_doping = Rni
         cell_volume = dx1*dy1*dz1
         norm_factor = cell_volume*norm_doping
         norm_factor = 1./norm_factor
 
         limx = int(domain_length*2/dx1/4 + 0.5)
         limy = int(domain_width*2/dy1/4 + 0.5)
         limz = int(z_depth*2/dz1/4 + 0.5)
        
         open(unit=50,file='charge_txt'//trim(extension)
     1                ,status='unknown')
         open(unit=51,file='charge_next_txt'//trim(extension),
     1                status='unknown')
         open(unit=52,file='charge_dopants_txt'//trim(extension),
     1                status='unknown')

        nn = 0
        do n1 = -2*limx,2*limx,1
        do n2 = -2*limy,2*limy,1       
        do n3 = -2*limz,2*limz,1
          x = a0/2*(n1+n3) 
          y = a0/2*(n1+n2)
          z = a0/2*(n2+n3)
C        (assigning arbitrary charges)                  
 23       rr = ran(iso)
          if(rr.lt.1.D-5)go to 23 
          if(rr.gt.0.5)charge = ran(iso)
          if(rr.le.0.5)charge = -ran(iso)
          charge = 0
          if(x.ge.0.and.x.le.domain_length.and.y.ge.0.and.y.le.
     &      domain_width.and.z.ge.0.and.z.le.z_depth)then
            nn=nn+1
            xion(nn)=x
            yion(nn)=y
            zion(nn)=z
            qion(nn)=charge*norm_factor/10000     
            write(50,80)x,y,z,charge*norm_factor/10000
            flag(nn) = 0
          endif
           
          x = x + a0/4
          y = y + a0/4
          z = z + a0/4
 24       rr = ran(iso)
          if(rr.lt.1.D-5)go to 24 
          if(rr.gt.0.5)charge = ran(iso)
          if(rr.le.0.5)charge = -ran(iso)
          charge = 0
          if(x.ge.0.and.x.le.domain_length.and.y.ge.0.and.y.le.
     &      domain_width.and.z.ge.0.and.z.le.z_depth)then
            nn=nn+1
            xion(nn)=x
            yion(nn)=y
            zion(nn)=z
            qion(nn)=charge*norm_factor/10000      
            write(51,80)x,y,z,charge*norm_factor/10000
            flag(nn) = 0
          endif
        enddo
        enddo
        enddo

        do i = 1,n_dopants
           nn=nn+1
           xion(nn)=x_dop(i)
           yion(nn)=y_dop(i)
           zion(nn)=z_dop(i)
           qion(nn)=charge_dop(i)*norm_factor      
           write(52,80)x_dop(i),y_dop(i),z_dop(i),
     &             charge_dop(i)*norm_factor
           flag(nn) = 1
        enddo            
  80    format(4(1x,e12.6))
        close(50)
        close(51)
        close(52)

        print*,'Total number of atoms (+dopants) = ',nn
      
        do i = 0,nx
        do j = 0,ny
        do k = 0,nz
           rion_number(i,j,k) = 0.   
C         (charge density)
        enddo
        enddo
        enddo
        
        do ii=1,nn
         i=int(xion(ii)/dx1 + 0.5)
         j=int(yion(ii)/dy1 + 0.5)
         k=int(zion(ii)/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
         rion_number(i,j,k)=rion_number(i,j,k) + qion(ii)  
         flag_reg(i,j,k) = flag(ii) 
        enddo
   
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      open(unit=233,file='ChargeDensity2Dxz_txt'//trim(extension),
     1           status='unknown')
      open(unit=232,file='ChargeDensity2Dxy_txt'//trim(extension),
     1           status='unknown')
      do i = 0,nx
      do j = 0,ny
      do k = 0,nz
         denn = rion_number(i,j,k)
         if(denn.gt.0.0)then
            dop(i,j,k) = -denn  - Rni/norm_doping   
            fi(i,j,k) = log(-dop(i,j,k))
         else
            dop(i,j,k) = -denn  + Rni/norm_doping
            fi(i,j,k) = -log(dop(i,j,k))
         endif
        if(k.eq.k_fix)write(232,*)i,j,-dop(i,j,k)
        if(j.eq.ny/2)write(233,*)i,k,-dop(i,j,k)
           flag_dom(i,j,k) = 4
           n(i,j,k) = exp(fi(i,j,k))
           p(i,j,k) = exp(-fi(i,j,k))
C           fi(i,j,k) = 0.0
C           n(i,j,k) = 0.D-5
C           p(i,j,k) = 0.D-5
           ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
           if(flag_reg(i,j,k).eq.1)eps(i,j,k) = eps_p
           if(flag_reg(i,j,k).eq.0)eps(i,j,k) = eps_si
      enddo
      enddo
      enddo
      close(232)
      close(233)

      if(z_ox.gt.0)then      
C     Gate Oxide region
      do i = 0,nx
      do j = 0,ny
      do k = -k_ox,0
        dop(i,j,k) = 0.0   
        flag_dom(i,j,k) = 1
        fi(i,j,k) = 0.0
        n(i,j,k) = 0.0
        p(i,j,k) = 0.0
        eps(i,j,k) = eps_ox
        ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
      enddo
      enddo
      enddo
     
C.....boundary condition 
      k = -k_ox
      DO i = 0, nx
      DO j = 0, ny
         flag_dom(i,j,k) = 5
         fi(i,j,k) = gate_voltage/Vt
         ro(i,j,k) = fi(i,j,k)
      ENDDO
      ENDDO
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....P/Si system: HOON: Multiple P atoms

      case('SiP_HOON')
         
         argc = iargc()
         if(argc.ne.1) then
          write(*,*)'Use: ../bin/poisson3d.ex <file with x,y,z,q,atom>'
          write(*,*)' '
          call exit(1)
         endif

         call getarg(1, ratom_file)

         xmax = -1.d38
         xmin = 1.d38
         ymax = -1.d38
         ymin = 1.d38
         zmax = -1.d38
         zmin = 1.d38

         dx1 = a0/discret_factor_x
         dy1 = a0/discret_factor_y
         dz1 = a0/discret_factor_z

         norm_doping = Rni
         cell_volume = dx1*dy1*dz1
         norm_factor = cell_volume*norm_doping
         norm_factor = 1./norm_factor

         nn = 0
         open(unit=1,file=ratom_file,status='old')
   16    read(1,*,err=17,end=17)x,y,z,charge,atom_number
         nn=nn+1
         x=x*1.d-9
         xion(nn) = x
         xmax = dmax1(x,xmax)
         xmin = dmin1(x,xmin)
         y=y*1.d-9
         yion(nn) = y
         ymax = max(y,ymax)
         ymin = min(y,ymin)
         z=z*1.d-9
         zion(nn) = z
         zmax = max(z,zmax)
         zmin = min(z,zmin)
         qion(nn) = -charge*norm_factor
         if(atom_number.eq.7.D0)then
           flag(nn)=0
         elseif(atom_number.eq.8.D0)then
           flag(nn)=1
         endif
         go to 16      
   17    close(1)         
         print *,'reading input file complete'
         print*,'Total atom (including Phosphorous) number = ',nn

C.....Calculate the grid and vectors
         
         z_ox=0  
         k_ox=0  
         domain_length = xmax-xmin 
         domain_width = ymax-ymin 
         domain_depth = zmax-zmin + z_ox
         print*,' '
         print*,'Domain length (nm) =',domain_length*1.0e9
         print*,'Domain width (nm) =',domain_width*1.0e9
         print*,'Domain depth (nm) =',domain_depth*1.0e9

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int((zmax-zmin)/dz1 + 0.5)

         do i = 0, nx-1
            x_vec(i) = dx1/debye
         enddo
         print *,'x-grid calculated, nx =',nx
     
         do j = 0, ny-1
            y_vec(j) = dy1/debye
         enddo
         print *,'y-grid calculated, ny =',ny

         do k = 0, nz-1
            z_vec(k) = dz1/debye
         enddo
         print *,'z-grid calculated, nz =',nz
 
      
        do i = 0,nx
        do j = 0,ny
        do k = 0,nz
           rion_number(i,j,k) = 0.   
C         (charge density assigned to grids)
        enddo
        enddo
        enddo
        
      do ii=1,nn
         i=int(xion(ii)/dx1 + 0.5)
         j=int(yion(ii)/dy1 + 0.5)
         k=int(zion(ii)/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
         rion_number(i,j,k)=rion_number(i,j,k) + qion(ii)  
         flag_reg(i,j,k) = flag(ii) 
      enddo
   
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      j_fix = int(y_contour*1.e-9/dy1 + 0.5)
      i_fix = int(x_contour*1.e-9/dx1 + 0.5)
      open(unit=232,file='ChargeDensity2Dxy_txt'//trim(extension),
     1           status='unknown')
      open(unit=233,file='ChargeDensity2Dxz_txt'//trim(extension),
     1           status='unknown')
      open(unit=234,file='ChargeDensity2Dyz_txt'//trim(extension),
     1           status='unknown')
      do i = 0,nx
      do j = 0,ny
      do k = 0,nz
         denn = rion_number(i,j,k)
         if(denn.gt.0.0)then
            dop(i,j,k) = -denn - Rni/norm_doping
            if(flag_reg(i,j,k).eq.0)then
              fi(i,j,k) = log(Rni/norm_doping)
            elseif(flag_reg(i,j,k).eq.1)then
              fi(i,j,k) = log(-dop(i,j,k))
            endif
         else
            dop(i,j,k) = -denn + Rni/norm_doping
            if(flag_reg(i,j,k).eq.0)then
              fi(i,j,k) = -log(Rni/norm_doping)
            elseif(flag_reg(i,j,k).eq.1)then
              fi(i,j,k) = -log(dop(i,j,k))
            endif
         endif
        if(k.eq.k_fix)write(232,*)i,j,-dop(i,j,k)
        if(j.eq.j_fix)write(233,*)i,k,-dop(i,j,k)
        if(i.eq.i_fix)write(234,*)j,k,-dop(i,j,k)
           flag_dom(i,j,k) = 4
           n(i,j,k) = exp(fi(i,j,k))
           p(i,j,k) = exp(-fi(i,j,k))
           fi(i,j,k) = 0.0
           n(i,j,k) = 0.D-5
           p(i,j,k) = 0.D-5
           ro(i,j,k) = n(i,j,k) - p(i,j,k) + dop(i,j,k)
           if(flag_reg(i,j,k).eq.1)eps(i,j,k) = eps_p
           if(flag_reg(i,j,k).eq.0)eps(i,j,k) = eps_si
      enddo
      enddo
      enddo
      close(232)
      close(233)
      close(234)


      end select
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....Common for all cases/devices
C.....For all cases plot the charge for ions

      open(unit=202,file='charge_txt'//trim(extension),
     1               status='unknown')
      open(unit=203,file='charge_negative_txt'//trim(extension),
     &     status='unknown')

      do ii = 1,nn
       if(FN.eq.'DOT'.or.FN.eq.'CQD')then      
         i=int(xion(ii)/dx1 + 0.5)
         j=int(yion(ii)/dy1 + 0.5)
         k=int(zion(ii)/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
         qion(ii) = rion_number(i,j,k) 
       endif
      if(device_name.eq.'SiP_DEF')then
         qion(ii) = qion(ii)*q*norm_doping
      elseif(device_name.eq.'SiP_HOON')then
         qion(ii) = qion(ii)*q*norm_doping
      endif 
       write(202,43)xion(ii),yion(ii),zion(ii),qion(ii)
       if(qion(ii).lt.0)write(203,43)xion(ii),yion(ii),zion(ii),qion(ii)
      enddo
  43  format(5(1x,e12.5))
      close(203)
      close(202)

       
C.....OPEN FILES FOR WRITING THE AXES:
      open(unit =  1, file = 'x_axis'//trim(strain_file),
     &     status = 'unknown')
      open(unit =  2, file = 'y_axis'//trim(strain_file),
     &     status = 'unknown')
      open(unit =  3, file = 'z_axis'//trim(strain_file),
     &     status = 'unknown')

      x = 0.D0
      do i = 0, nx
        write(1,*)x*debye
        x = x + x_vec(i)
      enddo
      close(1)
      y = 0.D0
      do j = 0, ny
        write(2,*)y*debye
        y = y + y_vec(j)
      enddo
      close(2)       
      z = 0.D0
      do k = -k_ox, nz
        write(3,*)z*debye
        z = z + z_vec(k)
      enddo
      close(3)
       

      return 
      end 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Subroutines/Functions used in the Structure code

      subroutine domain_grid(z_ox,x_vec,y_vec,z_vec)

      implicit real*8(A-H, O-Z)

      include 'params.h'
      include 'common_par.inc'

      real*8 x_vec(-1:i_max)
      real*8 y_vec(-1:j_max)
      real*8 z_vec(k_neg:k_max)


C.....Find the domain size

         xmax = -1.e38
         xmin = 1.e38
         ymax = -1.e38
         ymin = 1.e38
         zmax = -1.e38
         zmin = 1.e38

         nn = 0
         
         open(unit=1,file=strain_file,status='old')
   7     read(1,*,err=18,end=18)x,y,z,exx,eyy,ezz,exy,exz,eyz
         nn = nn + 1
         x=x*1.e-9
         xmax = max(x,xmax)
         xmin = min(x,xmin)
         y=y*1.e-9
         ymax = max(y,ymax)
         ymin = min(y,ymin)
         z=z*1.e-9
         zmax = max(z,zmax)
         zmin = min(z,zmin)

         go to 7      
   18    close(1)         
      
C.....Calculate the grid and vectors

         domain_length = xmax-xmin 
         domain_width = ymax-ymin 
         domain_depth = zmax-zmin + z_ox
         print*,' '
         print*,'Domain length (nm) =',domain_length*1.0e9
         print*,'Domain width (nm) =',domain_width*1.0e9
         print*,'Domain depth (nm) =',domain_depth*1.0e9

         nx = int(domain_length/dx1 + 0.5)
         ny = int(domain_width/dy1 + 0.5)
         nz = int((zmax-zmin)/dz1 + 0.5)
       
         do i = 0, nx-1
            x_vec(i) = dx1
         enddo
         print*,'x-grid calculated, nx =', nx
     
         do j = 0, ny-1
            y_vec(j) = dy1
         enddo
         print *,'y-grid calculated, ny =', ny

         if(z_ox.ne.0)then
          zo = 0.1e-9             
C        (initial grid step)
          r_in = 1.1e0
          r = r_coeff(r_in,zo,k_ox,z_ox)
          z_vec(-1) = zo
          do k = -2, -k_ox, -1
            z_vec(k) = r*z_vec(k+1)
          enddo
         endif

         do k = 0, nz-1
            z_vec(k) = dz1
         enddo
         print *,'z-grid calculated, nz =', nz

      return
      end



      FUNCTION R_COEFF(R_IN, XO, N, XT)
      implicit real*8(A-H, O-Z)
      LOGICAL FLAGCONV
      flagconv = .false.
      IF(N.GT.1)THEN
        r = r_in
        do while(.not.flagconv)
          if(r.eq.1.D0)then
            rnum = dfloat(N)*xo - xt
            denom = 0.5D0*xo*dfloat(N)*(dfloat(N)-1.D0)
          else
            term1 = 1.D0/(r-1.D0)
            term2 = r**dfloat(N)-1.D0
            rnum  = xo*term1*term2 - xt
            term3 = dfloat(N)*r**dfloat(N-1)
            denom = xo*(term1*term3-term1*term1*term2)
          endif
          r_new = r - rnum/denom
          error = dabs(r-r_new)/dabs(r_new)
          if(error.le.1.D-15)then
            flagconv = .true.
          else
            r = r_new
          endif
          if(r.lt.0.D0)then
            r = 1.D0
            flagconv = .true.
          endif
        enddo
      ELSE
        r = 1.D0
      ENDIF
      r_coeff = r
      
      return
      end

      FUNCTION rand(iso)
      INTEGER iso,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!Long period (> 2×1018) random number generator of L’Ecuyer with Bays-Durham shu.e
!and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
!of the endpoint values). Call with iso a negative integer to initialize; thereafter, do not
!alter iso between successive deviates in a sequence. RNMX should approximate the largest
!.oating value that is less than 1.
      INTEGER iso2,j,k,iv(NTAB),iy
      SAVE iv,iy,iso2
      DATA iso2/123456789/, iv/NTAB*0/, iy/0/

        if (iso.le.0) then 
        iso=max(-iso,1) 
        iso2=iso
        do j=NTAB+8,1,-1 
        k=iso/IQ1

!Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
!Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machinereadable
!files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books or CDROMs, visit website
!http://www.nr.com or call 1-800-872-7423 (North America only), or send email to directcustserv@cambridge.org (outside North America).
        iso=IA1*(iso-k*IQ1)-k*IR1
        if (iso.lt.0) iso=iso+IM1
        if (j.le.NTAB) iv(j)=iso
        enddo 
        iy=iv(1)
        endif
        k=iso/IQ1 
        iso=IA1*(iso-k*IQ1)-k*IR1 
        if (iso.lt.0) iso=iso+IM1
        k=iso2/IQ2
        iso2=IA2*(iso2-k*IQ2)-k*IR2 
        if (iso2.lt.0) iso2=iso2+IM2
        j=1+iy/NDIV 
        iy=iv(j)-iso2 
        iv(j)=iso
        if(iy.lt.1)iy=iy+IMM1
        rand=min(AM*iy,RNMX) 
        return
        END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                   IMPLEMENTATION OF ILU METHOD
C              FOR THE SOLUTION OF THE POISSON EQUATION
C          It does not use an internal loop for convergence.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE ILU_POISSON(Y_VEC, B_VEC, C_VEC,
     1    D_VEC, E_VEC, F_VEC, G_VEC, H_VEC, FLAG_DOM, RO, FI,
     2    N, P)

      implicit real*8(A-H, O-Z)       
      include 'params.h'
      INCLUDE 'common_par.inc'

      real*8 Y_VEC(-1:j_max)

C     Coefficients from the finite difference scheme:
      real*8 B_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 C_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 D_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 E_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 F_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 G_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 H_vec(0:i_max,0:j_max,k_neg:k_max)

C     Elements of the L- and U- matrices:
      real*8 d1(-1:i_max+1,-1:j_max,k_neg:k_max+1)

C     Elements used in the SiP method:
      real*8 R_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 V_vec(-1:i_max+1,0:j_max,k_neg:k_max+1)

C     Others:
      integer FLAG_DOM(0:i_max,0:j_max,k_neg:k_max)
      real*8 RO(0:i_max,0:j_max,k_neg:k_max)
      real*8 FI(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 DELTA(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 N(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 P(-1:i_max+1,0:j_max,k_neg:k_max+1)



C......CALCULATION OF THE MATRIX ELEMENT d1(i,j,k):

      DO I = 0, nx
      DO J = 0, ny
      DO K = -K_OX, nz
      IF(FLAG_DOM(I,J,K).NE.0)THEN
        IF(FLAG_DOM(I,J,K).EQ.4)THEN
          E_vec_new = E_vec(i,j,k)-(n(i,j,k) + p(i,j,k))
        ELSE
          E_vec_new = E_vec(i,j,k)
        ENDIF       
        IF(B_VEC(I,J,K).NE.0)THEN
          sum1 = B_vec(i,j,k)*H_vec(i,j,k-1)/d1(i,j,k-1)
        ELSE
          sum1 = 0.
        ENDIF
        IF(C_VEC(I,J,K).NE.0)THEN
          sum2 = C_vec(i,j,k)*G_vec(i,j-1,k)/d1(i,j-1,k)
        ELSE
          sum2 = 0.
        ENDIF
        IF(D_VEC(I,J,K).NE.0)THEN
          sum3 = D_vec(i,j,k)*F_vec(i-1,j,k)/d1(i-1,j,k)
        ELSE
          sum3 = 0.
        ENDIF
        d1(i,j,k) =  E_vec_new - (sum1 + sum2 + sum3)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

C......CALCULATE THE IMPROVED FI-VECTOR USING INCOMPLETE LU-DECOMPOSITION:
C......Calculate the vector R_VEC:

      DO I = 0, nx
      DO J = 0, ny
      DO K = -K_OX, nz
      IF(FLAG_DOM(I,J,K).NE.0)THEN
        IF(B_VEC(I,J,K).NE.0)THEN
          term1 = F_vec(i,j,k-1)*fi(i+1,j,k-1)
          term2 = G_vec(i,j,k-1)*fi(i,j+1,k-1)
          sum1  = B_vec(i,j,k)*(term1 + term2)/d1(i,j,k-1)
        ELSE
          sum1 = 0.
        ENDIF
        IF(C_VEC(I,J,K).NE.0)THEN
          term1 = F_vec(i,j-1,k)*fi(i+1,j-1,k)
          term2 = H_vec(i,j-1,k)*fi(i,j-1,k+1)
          sum2  = C_vec(i,j,k)*(term1 + term2)/d1(i,j-1,k)
        ELSE
          sum2 = 0.
        ENDIF
        IF(D_VEC(I,J,K).NE.0)THEN
          term1 = G_vec(i-1,j,k)*fi(i-1,j+1,k)
          term2 = H_vec(i-1,j,k)*fi(i-1,j,k+1)
          sum3  = D_vec(i,j,k)*(term1 + term2)/d1(i-1,j,k)
        ELSE
          sum3 = 0.
        ENDIF
        IF(FLAG_DOM(I,J,K).EQ.4)THEN
          sum4 = (n(i,j,k) + p(i,j,k))*fi(i,j,k)
        ELSE
          sum4 = 0.
        ENDIF
        R_vec(i,j,k) = RO(i,j,k) + sum1 + sum2 + sum3 - sum4
      ENDIF
      ENDDO
      ENDDO
      ENDDO

C......Calculate the vector V_VEC:

      DO I = 0, nx
      DO J = 0, ny
      DO K = -K_OX, nz
      IF(FLAG_DOM(I,J,K).NE.0)THEN
        IF(B_VEC(I,J,K).NE.0)THEN
          sum1 = B_vec(i,j,k)*V_vec(i,j,k-1)/d1(i,j,k-1)
        ELSE
          sum1 = 0.
        ENDIF
        IF(C_VEC(I,J,K).NE.0)THEN
          sum2 = C_vec(i,j,k)*V_vec(i,j-1,k)/d1(i,j-1,k)
        ELSE
          sum2 = 0.
        ENDIF
        IF(D_VEC(I,J,K).NE.0)THEN
          sum3 = D_vec(i,j,k)*V_vec(i-1,j,k)/d1(i-1,j,k)
        ELSE
          sum3 = 0.
        ENDIF
        V_vec(i,j,k) = R_vec(i,j,k) - (sum1 + sum2 + sum3)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

C......Calculate the new FI vector:

      DO I = nx, 0, -1
      DO J = ny, 0, -1       
      DO K = nz, -K_OX, -1
      IF(FLAG_DOM(I,J,K).NE.0)THEN
        IF(F_VEC(I,J,K).NE.0)THEN
          sum1 = F_vec(i,j,k)*fi(i+1,j,k)
        ELSE
          sum1 = 0.
        ENDIF
        IF(G_VEC(I,J,K).NE.0)THEN
          sum2 = G_vec(i,j,k)*fi(i,j+1,k)
        ELSE
          sum2 = 0.
        ENDIF
        IF(H_VEC(I,J,K).NE.0)THEN
          sum3 = H_vec(i,j,k)*fi(i,j,k+1)
        ELSE
          sum3 = 0.
        ENDIF
        sum = sum1 + sum2 + sum3
        term = (V_vec(i,j,k) - sum)/d1(i,j,k)
        delta(i,j,k) = term - fi(i,j,k)
        fi(i,j,k) = term
      ENDIF
      ENDDO
      ENDDO
      ENDDO

C......CALCULATE THE ERROR :

      delta_max = 0.
      DO I = 0, nx
      DO J = 0, ny
      DO K = -K_OX, nz
      IF(FLAG_DOM(I,J,K).NE.0)THEN
        error = dabs(delta(i,j,k))
        IF(error.gt.delta_max)THEN
          i_err = i
          j_err = j
          k_err = k
          delta_max = error
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO              


      RETURN
      END              


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        SUBROUTINE THAT CALLS VARIOUS SUBROUTINES FOR THE CALCULATION
C               OF THE COEFFICIENTS B, C, D, E, F, G, H
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE COEFF_P(X_VEC, Y_VEC, Z_VEC, FLAG_DOM,
     1            B_VEC, C_VEC, D_VEC, E_VEC, F_VEC, G_VEC, H_VEC
     2            , eps)
      implicit real*8(A-H, O-Z)
       
      include 'params.h'
      INCLUDE 'common_par.inc'

      real*8 X_VEC(-1:i_max)
      real*8 Y_VEC(-1:j_max)
      real*8 Z_VEC(k_neg:k_max)

      integer FLAG_DOM(0:i_max,0:j_max,k_neg:k_max)

C     Coefficients from the finite difference scheme:
      real*8 B_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 C_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 D_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 E_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 F_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 G_vec(0:i_max,0:j_max,k_neg:k_max)
      real*8 H_vec(0:i_max,0:j_max,k_neg:k_max)

      real*8 eps(-1:i_max+1,0:j_max,k_neg:k_max+1)


      DO I = 0, nx
      DO J = 0, ny       
      DO K = -K_OX, nz


      select case(flag_dom(i,j,k))

        case (1,4)
          CALL BULK(X_VEC, Y_VEC, Z_VEC, coeff_B, coeff_C, coeff_D,
     1              coeff_E, coeff_F, coeff_G, coeff_H, I, J, K
     2              , eps)

        case (5)
          coeff_B = 0.
          coeff_C = 0.
          coeff_D = 0.
          coeff_E = 1.
          coeff_F = 0.
          coeff_G = 0.
          coeff_H = 0.

      end select

        B_vec(i,j,k) = coeff_B
        C_vec(i,j,k) = coeff_C
        D_vec(i,j,k) = coeff_D
        E_vec(i,j,k) = coeff_E
        F_vec(i,j,k) = coeff_F
        G_vec(i,j,k) = coeff_G
        H_vec(i,j,k) = coeff_H

      ENDDO
      ENDDO
      ENDDO 

      RETURN
      END         


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Bulk Coefficient 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE BULK(X_VEC, Y_VEC, Z_VEC, coeff_B, coeff_C, coeff_D,
     1                coeff_E, coeff_F, coeff_G, coeff_H, I, J, K
     2                , eps)
      implicit real*8(A-H, O-Z)
       
      include 'params.h'    
      INCLUDE 'common_par.inc'

      real*8 X_VEC(-1:i_max)
      real*8 Y_VEC(-1:j_max)
      real*8 Z_VEC(k_neg:k_max)

      real*8 eps(-1:i_max+1,0:j_max,k_neg:k_max+1)

     
      method = 2
      
      IF (method.EQ.1)THEN

       if(i.eq.0)then
          coeff_d = 0.
          coeff_f = eps(i,j,k)/(x_vec(i)*x_vec(i)) + 
     1              (eps(i+1,j,k)-eps(i,j,k))/(x_vec(i)*x_vec(i))

       elseif(i.eq.nx)then
          coeff_f = 0.
          coeff_d = eps(i,j,k)/(x_vec(i-1)*x_vec(i-1))

       else
          coeff_d = eps(i,j,k)/(x_vec(i-1)*(x_vec(i)+x_vec(i-1)))
          coeff_f = eps(i,j,k)/(x_vec(i)*(x_vec(i)+x_vec(i-1)))+
     1            (eps(i+1,j,k)-eps(i,j,k))/(2.*x_vec(i)*x_vec(i))
       endif

       if(j.eq.0)then
          coeff_c = 0.
          coeff_g = eps(i,j,k)/(y_vec(j)*y_vec(j)) +
     1              (eps(i,j+1,k)-eps(i,j,k))/(y_vec(j)*y_vec(j))

       elseif(j.eq.ny)then
          coeff_g = 0.
          coeff_c = eps(i,j,k)/(y_vec(j-1)*y_vec(j-1))

       else
          coeff_c = eps(i,j,k)/(y_vec(j-1)*(y_vec(j)+y_vec(j-1)))
          coeff_g = eps(i,j,k)/(y_vec(j)*(y_vec(j)+y_vec(j-1))) +
     1              (eps(i,j+1,k)-eps(i,j,k))/(2.*y_vec(j)*y_vec(j))
       endif 

       if(k.eq.-k_ox)then
          coeff_b = 0.
          coeff_h = eps(i,j,k)/(z_vec(k)*z_vec(k)) +
     1              (eps(i,j,k+1)-eps(i,j,k))/(z_vec(k)*z_vec(k))

       elseif(k.eq.nz)then
          coeff_h = 0.
          coeff_b = eps(i,j,k)/(z_vec(k-1)*z_vec(k-1))

       else
          coeff_b = eps(i,j,k)/(z_vec(k-1)*(z_vec(k)+z_vec(k-1)))
          coeff_h = eps(i,j,k)/(z_vec(k)*(z_vec(k)+z_vec(k-1))) +
     1              (eps(i,j,k+1)-eps(i,j,k))/(2.*z_vec(k)*z_vec(k))
      endif

      coeff_e = - (coeff_b+coeff_c+coeff_d+coeff_f+
     1             coeff_g+coeff_h)

      ENDIF


      IF (method.EQ.2)THEN

       if(i.eq.0)then
          coeff_d = 0.
          coeff_f = 4.*eps(i,j,k)/(x_vec(i)*(x_vec(i)+x_vec(i-1)))
       elseif(i.eq.nx)then
          coeff_f = 0.
          coeff_d = 4.*eps(i-1,j,k)/(x_vec(i-1)*(x_vec(i)+x_vec(i-1)))
       else
          coeff_d = 2*eps(i-1,j,k)/(x_vec(i-1)*(x_vec(i)+x_vec(i-1)))
          coeff_f = 2*eps(i,j,k)/(x_vec(i)*(x_vec(i)+x_vec(i-1)))
       endif

       if(j.eq.0)then
          coeff_c = 0.
          coeff_g = 4.*eps(i,j,k)/(y_vec(j)*(y_vec(j)+y_vec(j-1)))
       elseif(j.eq.ny)then
          coeff_g = 0.
          coeff_c = 4.*eps(i,j-1,k)/(y_vec(j-1)*(y_vec(j)+y_vec(j-1)))
       else
          coeff_c = 2.*eps(i,j-1,k)/(y_vec(j-1)*(y_vec(j)+y_vec(j-1)))
          coeff_g = 2.*eps(i,j,k)/(y_vec(j)*(y_vec(j)+y_vec(j-1))) 
       endif

       if(k.eq.-k_ox)then
          coeff_b = 0.
          coeff_h = 4.*eps(i,j,k)/(z_vec(k)*(z_vec(k)+z_vec(k-1)))
       elseif(k.eq.nz)then
          coeff_h = 0.
          coeff_b = 4.*eps(i,j,k-1)/(z_vec(k-1)*(z_vec(k)+z_vec(k-1)))
       else
          coeff_b = 2.*eps(i,j,k-1)/(z_vec(k-1)*(z_vec(k)+z_vec(k-1)))
          coeff_h = 2.*eps(i,j,k)/(z_vec(k)*(z_vec(k)+z_vec(k-1))) 
       endif

      coeff_e = - (coeff_b+coeff_c+coeff_d+coeff_f+
     1             coeff_g+coeff_h)

      ENDIF 

     
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Subroutine for writing results:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE write_data_eff(device_name,x_vec,y_vec,z_vec,n,p,dop,
     &    fi,eps,flag_dom)
      implicit real*8(a-h,o-z)

      include 'params.h'
      include 'common_par.inc'


      real*8 X_VEC(-1:i_max)
      real*8 Y_VEC(-1:j_max)
      real*8 Z_VEC(k_neg:k_max)
      real*8 Ec_band(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 N(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 P(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 RO(0:i_max,0:j_max,k_neg:k_max)
      real*8 DOP(0:i_max,0:j_max,k_neg:k_max)
      real*8 fi(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 el_field(0:i_max,k_neg:k_max)
      real*8 eps(-1:i_max+1,0:j_max,k_neg:k_max+1)
      real*8 sumx(-1:i_max+1)
      real*8 sumy(0:j_max)
      real*8 sum2(0:j_max,k_neg:k_max+1)
      real*8 sum3(-1:i_max+1,k_neg:k_max+1)
      integer FLAG_DOM(0:i_max,0:j_max,k_neg:k_max)
      
      character (100) extension,FN
      character (100) device_name

                
      i_fix = nx/2
      j_fix = ny/2        
      k_fix = nz/2

      select case (device_name)

      case('DOT')

      FN = 'DOT'
      extension = trim(FN)
      extension = trim(strain_file)

      open(unit=15,file='potential3D_txt'//trim(extension))
      open(unit=16,file='potential3Dneg_txt'//trim(extension))
      open(unit=1,file='potential')
C    ('potential' file is required for NEMO 3-D)
      open(unit=2,file='potential1dxy_txt'//trim(extension))
      open(unit=6,file='potential1dxmy_txt'//trim(extension))
      open(unit=3,file='potential1dz_txt'//trim(extension))
      open(unit=4,file='potential2dxy_txt'//trim(extension))

      open(unit=5,file=strain_file,status='old')
  73  read(5,*,err=55,end=55)x,y,z,exx,eyy,ezz,exy,exz,eyz
         i=int(x*1.e-9/dx1 + 0.5)
         j=int(y*1.e-9/dy1 + 0.5)
         k=int(z*1.e-9/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
       
         pot1 = fi(i,j,k)

         write(15,71)x,y,z,pot1
         if(pot1.lt.0)write(16,71)x,y,z,pot1

         write(1,4)pot1

         if((y/x.ge.0.99.and.y/x.le.1.01)
     &      .and.(int(z/dz1+0.5).eq.nz/2)) 
     &      write(2,7)x*1e9,sqrt(x*x+y*y)*1e9,pot1
         if(((ymax-y)/x.ge.0.98.and.(ymax-y)/x.le.1.02)
     &      .and.(int(z/dz1+0.5).eq.nz/2)) 
     &      write(6,7)x*1e9,sqrt(x*x+y*y)*1e9,pot1
         if(int(x/dx1+0.5).eq.nx/2.and.int(y/dy1+0.5).eq.ny/2)
     &      write(3,8)z*1e9,pot1
         if(int(z/dz1+0.5).eq.nz/2)
     &      write(4,7)x*1e9,y*1e9,pot1

        goto 73
  55  close(5)
  
  70  format(9(1x,e12.6))
  71  format(4(1x,e12.6))
  4   format(e12.6)
  7   format(3(2x,e12.6))
  8   format(2(2x,e12.6))
      close(1)       
      close(2)       
      close(3)       
      close(4)    
      close(6)   
      close(15)   
      close(16)   

 
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      open(unit=1,file='potential_contour_xy_'//trim(extension))
      do j = 0, ny
      do i = 0, nx
         write(1,*)fi(i,j,k_fix)
      enddo
      enddo
      close(1)

      j_fix = ny/2
      open(unit=1,file='potential_contour_xz_'//trim(extension))
      do i = 0, nx
      do k = 0, nz
         write(1,*)fi(i,j_fix,k)
      enddo
      enddo
      close(1)

      i_fix = nx/2
      j_fix = ny/2
      open(unit=1,file='potential_faiZ'//trim(extension))
      do k = 0, nz
         write(1,*)fi(i_fix,j_fix,k)
      enddo
      close(1)

      open(unit=1,file='potential_faiXY'//trim(extension))
      do i = 0, nx
      do j = 0, ny
         x=i*dx1*1e9
         y=j*dy1*1e9
         if(i.eq.j)write(1,*)x,sqrt(x*x+y*y),fi(i,j,k_fix)
      enddo
      enddo
      close(1)

      open(unit=1,file='potential_faiXmY'//trim(extension))
      do i = 0, nx
      do j = 0, ny
         x=i*dx1*1e9
         y=j*dy1*1e9
         if(i.eq.ny-j)write(1,*)x,sqrt(x*x+y*y),fi(i,j,k_fix)
      enddo
      enddo
      close(1)



C.....SiP...............
      case('SiP_DEF')

      FN = 'SiP_DEF'
      extension = trim(FN)

      open(unit=15,file='potential3D_txt'//trim(extension))
      open(unit=16,file='potential3Dneg_txt'//trim(extension))
      open(unit=1,file='potential')
C    ('potential_txt_' file is required for NEMO 3-D.
C      Please converter the name to simply 'potential')

      open(unit=5,file='charge_txt'//trim(extension),status='old')
  74  read(5,*,err=56,end=56)x,y,z,q1
         i=int(x/dx1 + 0.5)
         j=int(y/dy1 + 0.5)
         k=int(z/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
       
         pot1 = fi(i,j,k)*Vt*1000

         write(1,4)pot1
         write(15,71)x,y,z,pot1
         if(pot1.lt.0)write(16,71)x,y,z,pot1

        goto 74
  56  close(5)
  
      close(1)       
      close(2)       
      close(3)       
      close(4)    
      close(6)   
      close(15)   
      close(16)   

 
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)
      open(unit=1,file='potential_contour_xy_'//trim(extension))
      do j = 0, ny
      do i = 0, nx
         write(1,*)fi(i,j,k_fix)*Vt*1000
      enddo
      enddo
      close(1)

      j_fix = ny/2
      open(unit=1,file='potential_contour_xz_'//trim(extension))
      do i = 0, nx
      do k = -k_ox, nz
         write(1,*)fi(i,j_fix,k)*Vt*1000
      enddo
      enddo
      close(1)


C.....SiP_HOON...............
      case('SiP_HOON')

      FN = 'SiP_HOON'
      extension = trim(FN)

      i_fix = int(x_contour*1.e-9/dx1 + 0.5)
      j_fix = int(y_contour*1.e-9/dy1 + 0.5)
      k_fix = int(z_contour*1.e-9/dz1 + 0.5)

      open(unit=1,file='potential')
C    ('potential' file is required for NEMO 3-D input)
      open(unit=2,file='potential1dx_txt'//trim(extension))
      open(unit=3,file='potential1dy_txt'//trim(extension))
      open(unit=4,file='potential2dz_txt'//trim(extension))
      open(unit=15,file='potential3D_txt'//trim(extension))
      open(unit=16,file='potential3Dneg_txt'//trim(extension))

      open(unit=5,file='charge_txt'//trim(extension),status='old')
  75  read(5,*,err=76,end=76)x,y,z,q1
         i=int(x/dx1 + 0.5)
         j=int(y/dy1 + 0.5)
         k=int(z/dz1 + 0.5)
         if(i.lt.0)i=0
         if(i.ge.nx)i=nx-1
         if(j.lt.0)j=0
         if(j.ge.ny)j=ny-1
         if(k.lt.0)k=0
         if(k.ge.nz)k=nz-1
       
         pot1 = fi(i,j,k)
         write(1,4)pot1

         if(j.eq.j_fix.and.k.eq.k_fix)write(2,7)x*1e9,pot1
         if(i.eq.i_fix.and.k.eq.k_fix)write(3,7)y*1e9,pot1
         if(i.eq.i_fix.and.j.eq.j_fix)write(4,7)z*1e9,pot1

         write(15,71)x,y,z,pot1
         if(pot1.lt.0)write(16,71)x,y,z,pot1

        goto 75
  76  close(5)
  
      close(1)       
      close(2)       
      close(3)       
      close(4)    
      close(15)   
      close(16)   

 
      open(unit=1,file='potential_contour_xy_'//trim(extension))
      do j = 0, ny
      do i = 0, nx
         write(1,*)fi(i,j,k_fix)
      enddo
      enddo
      close(1)

      open(unit=1,file='potential_contour_xz_'//trim(extension))
      do i = 0, nx
      do k = 0, nz
         write(1,*)fi(i,j_fix,k)
      enddo
      enddo
      close(1)

      open(unit=1,file='potential_contour_yz_'//trim(extension))
      do j = 0, ny
      do k = 0, nz
         write(1,*)fi(i_fix,j,k)
      enddo
      enddo
      close(1)

      open(unit=1,file='potential_faiZ'//trim(extension))
      do k = 0, nz
         write(1,*)fi(i_fix,j_fix,k)
      enddo
      close(1)
      open(unit=1,file='potential_faiX'//trim(extension))
      do i = 0, nx
         write(1,*)fi(i,j_fix,k_fix)
      enddo
      close(1)
      open(unit=1,file='potential_faiY'//trim(extension))
      do j = 0, ny
         write(1,*)fi(i_fix,j,k_fix)
      enddo
      close(1)


      endselect



      
      return
      END

