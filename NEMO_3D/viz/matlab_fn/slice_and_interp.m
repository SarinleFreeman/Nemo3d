function [r1, r2, f12] = slice_and_interp(g, usr)

%  slice_and_interp:  [r1, r2, f12] <-- (g, usr)
%
%==============================================================================
%
%  DESCRIPTION:  Matlab, by default, has difficulty with non-cartesian grids.
%                Plotting functions defined only at atomic sites is further
%                complicated by strain.  This function takes an nx4 matrix
%                of the form g = [x y z f(x,y,z)], obtains a slice along
%                a plane parallel to xy, yz, or xz, and returns a new f that 
%                is interpolated on a (typically much finer) rectangular 
%                2D grid.
%
%  INPUT:
%
%    real(r,4) g:    columns of g are x, y, z, f(x,y,z)
%
%    string usr:     this string is evaluated within matlab and is used
%                    to set the following variables:
%
%                    limits = a 4-vector of the form 
%                             [r1_min, r1_max, r2_min, r2_max]
%                             used to sepcify the domain of the interpolating 
%                             grid.  By default, it is defined to be the 
%                             smallest domain that fits all points in the 
%                             slice of f
%
%                    dmsh = a real number that equals the grid spacing
%
%                    The restriction of g to a plane is accomplished in
%                    the same manner as for restrict().
%                    Type 'help restrict' for more details
%
%  OUTPUT:
%
%    real(1,m) r1:   interpolation points along direction 1
%
%    real(n,1) r2:   interpolation points along direction 2
%
%    real(n,m) f12:  interpolation of g on 2D grid defined by r1 and r2
%
%  EXAMPLES:
%
%    Restrict an nx4 matrix g to the plane z=3.5 and interpolate:
%
%       [x,y,g_interp] = slice_and_interp(g,'z=3.5;');
%
%    Do as above, except restrict g to the plane closest to z=3.5 that
%    contains an atom, and only interpolate on the rectangle formed by
%    the points (0,0), (0,1), (2,0), (2,1), and refine the interpolated
%    mesh spacing to 0.01:
%
%       [x,y,g_interp] = slice_and_interp(g,'z=3.5; limits=[0 1 0 2]; ...
%                                            dmesh=0.01; snap=1;');
%
%  DEPENDENCIES:  restrict
%
%==============================================================================

if ( exist('usr') );  eval(usr);  end

existX = exist('x'); existY = exist('y'); existZ = exist('z');

if (~existX & ~existY & ~existZ)
   disp('ERROR:  A plane MUST be specified (e.g. usr=''z=3.5;'')');
   r1=[]; r2=[]; f12=[];
   return
end

if ( ~exist('dmsh') ); dmsh=0.06; end;

f = restrict(g, usr);

if ( existX )
   if ( existY | existZ); disp('Must specify only ONE plane'); return; end
   r1_0 = f(:,2);
   r2_0 = f(:,3);
end

if ( existY )
   if ( existX | existZ); disp('Must specify only ONE plane'); return; end
   r1_0 = f(:,1);
   r2_0 = f(:,3);
end

if ( existZ )
   if ( existX | existY); disp('Must specify only ONE plane'); return; end
   r1_0 = f(:,1);
   r2_0 = f(:,2);
end

if (~exist('limits')); limits = [min(r1_0) max(r1_0) min(r2_0) max(r2_0)]; end

r1=[limits(1):dmsh:limits(2)]; r2=[limits(3):dmsh:limits(4)]';

f12=griddata(r1_0,r2_0,f(:,4),r1,r2,'cubic');

