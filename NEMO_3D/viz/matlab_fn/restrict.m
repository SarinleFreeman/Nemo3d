function [gr] = restrict(g, usr)

%  restrict:  [real(n,4) gr] <-- (real(n,4) g, string usr)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    returns a subset of rows of the matrix g=[x y z f(x,y,z)]
%    determined by the contraints specified by usr.
%
%  INPUT:
%
%    nx4 matrix g:   columns of g are x, y, z, f(x,y,z)
%
%    string usr:     this string is evaluated within matlab and is used
%                    to set the following variables:
%
%                         {x, y, z, f, dx, dy, dz, df, dr, fLim, snap}
%                    
%                    The coordinate variables {x, y, z} specify a "thick"
%                    plane to which f is restricted.  The thicknesses of 
%                    these planes is first determined by 'dr'.  
%                    Alternatively, one can specify different thicknesses 
%                    for each dimension by setting {dx, dy, dz}.
%                    Although restricting f(x,y,z) in the same way as the
%                    coordinates is allowed, a more convenient form
%                    of restriction to specify minimum and maximum values
%                    of f though the 2-vector fLim.  Finally, if 'snap' 
%                    is defined, the midpoint of the thick planes are 
%                    snapped to the closest point.  For example, if
%
%                                      usr = 'x=5.0; snap=1;'
%
%                    the plane to which f is restricted will be reset
%                    to the closest value in the first column of g.
%
%  OUTPUT:
%
%    real(n,4) gr:   return data in the form of an nx4 matrix
%
%  EXAMPLES:
%
%    Restriction of g to a thin slab at z=3.0 with thickness 0.25:
%
%               g_restricted = restrict(g,'z=3.0; dz=0.25;');
%
%    Alternatively if dz is not specified, a default value of 0.15 is used.
%
%    Suppose we want to restrict g to a line approximately determined 
%    by z=3.0 and x=4.2, except that due to strain, we don't know exactly
%    where the lattice planes reside.  We can then snap the line to
%    the closest one that contains an actual atom as follows:
%
%               g_restricted = restrict(g,'z=3.0; x=4.2; snap=1;');
%
%    Suppose we want to find all g such that f(x,y,z) lies in the range
%    [0.6 1]:
%
%               g_restricted = restrict(g,'fLim=[0.6 1];');
%
%    Suppose we want to find all g such that f(x,y,z) is less than or 
%    equal to 0.2:
%
%               g_restricted = restrict(g,'fLim=[-inf 0.2];');
%
%  LIMITATIONS:  at the moment restricting planes have to be parallel
%                to either x, y, or z axes.
%
%  DEPENDENCIES:  none
%
%==============================================================================

%sanity check
if (size(g,2)~=4)
   disp(['ERROR: g should have 4 cols, not ' num2str(size(g,2))]);
   return;
end

if ( exist('usr') );  eval(usr); end
if ( ~exist('dr') ); dr=0.15; end
gr = g;

if ( exist('x') )
   gx = gr(:,1);
   if ( ~exist('dx') ); dx = dr; end
   if ( exist('snap') )
      [jnk indx] = min(abs(gx-x));
      x = gr(indx,1);
   end
   xMin = x - dx;
   xMax = x + dx;
   gr = gr(find(gx<xMax & gx>xMin),:);
end

if ( exist('y') )
   gy = gr(:,2);
   if ( ~exist('dy') ); dy = dr; end
   if ( exist('snap') )
      [jnk indx] = min(abs(gy-y));
      y = gr(indx,2);
   end
   yMin = y - dy;
   yMax = y + dy;
   gr = gr(find(gy<yMax & gy>yMin),:);
end

if ( exist('z') )
   gz = gr(:,3);
   if ( ~exist('dz') ); dz = dr; end
   if ( exist('snap') )
      [jnk indx] = min(abs(gz-z));
      z = gr(indx,3);
   end
   zMin = z - dz;
   zMax = z + dz;
   gr = gr(find(gz<zMax & gz>zMin),:);
end

if ( exist('fLim') )
   gf = gr(:,4);
   if ( fLim(1)==-inf )
      gr = gr(find(gf<fLim(2)),:);
   elseif ( fLim(2)==inf )
      gr = gr(find(gf>fLim(1)),:);
   else
      gr = gr(find(gf<fLim(2) & gf>fLim(1)),:);
   end
elseif ( exist('f') )
   gf = gr(:,4);
   if ( ~exist('df') ); df = dr; end
   if ( exist('fLim') )
      [jnk indx] = min(abs(gf-f));
      f = gr(indx,3);
   fMin = f - df;
   fMax = f + df;
   end
   gr = gr(find(gf<fMax & gf>fMin),:);
end

