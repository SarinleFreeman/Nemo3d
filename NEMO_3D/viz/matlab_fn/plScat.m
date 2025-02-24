function [H, Mv] = plScat(g, ampl, usr)

%  plScat:  [plot* H, movie Mv] <-- (real(n,4) g, vector ampl, string usr)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Generates an array of plot handles to scatter plots generated
%    by restricting the values of col 4 of g.  Optionally generate a movie.
%
%  INPUT:
%
%    nx4 matrix g:  columns of g are x, y, z, f(x,y,z)
%
%    vector ampl:   represents a sequence of contiguous disjoint intervals
%                   to which f(x,y,z) is restricted.  For example,
%                    
%                           ampl = [.3 .5 1];
%
%                   will produce one scatter plots such that 0.3 < f < 0.5 
%                   and another such that 0.5 < f < 1.0.  These two plots
%                   are colored differently.
%    
%    string usr:    this string is evaluated within matlab and is used
%                   to set the following optional variables:
%
%                   nrm:  if not defined, f is not normalized
%                         if equal to 2, f is normalized by the 2-norm
%                         if equal to inf, f is normalized by max(f)
%                    
%                   mov:  if defined, a new plot is created, making a movie.
%                         if not, the plots are superimposed.  
%
%                   big:  if defined, the plotting domain is enlarged to
%                         accommodate rotations
%             
%                   notr: if defined, no translation to centroid of figure
%                         is made.
%
%  OUTPUT:
%
%    plot(m) H:     array of handles to plots
%
%    movie  Mv:     if a movie is generated, this is non-null.
%
%  EXAMPLES:
%
%    Typical simple plot:  H = plScat(g, [1 .6 .5 .4])
%
%    Typical plot using normalized psi:  H = plScat(g, [1 .02 .005], 'nrm=2;')
%
%    Typical movie:  [H, Mv] = plScat(g,[1:-.01:.1],'mov=1;')
%                    movie(Mv);
%
%  COMMON ERRORS:  
%  
%    an error frequently made is that none of the amplitude ranges contains
%    any points in g(:,4).
%
%  DEPENDENCIES:  none
%
%==============================================================================

if ( exist('usr') );  eval(usr);  end

if (size(ampl,2)==1); ampl=[ampl 1.0]; end
N = size(ampl,2)-1;

colourMtrx=[0 0 1; 0 1 0; 1 0 0];

Req = g(:,1:3);
if ( exist('nrm') )
  if (nrm==2); psi2 = g(:,4) / norm(g(:,4)); end
  if (nrm==inf); psi2 = g(:,4) / max(abs(g(:,4))); end
else
  psi2 = g(:,4);
end

R=[];
iCntr=[1];
count=0;
for j=1:N
   r = Req(psi2 > ampl(j+1) & psi2 <= ampl(j),:);
   dum = iCntr(max(size(iCntr)))+size(r,1);
   if (dum ~= iCntr)
      iCntr = [iCntr iCntr(max(size(iCntr)))+size(r,1)];
      R = [R; r];
      count = count + 1;
   end
end
N = count;

if ( size(R,1)==0 )
   disp('Could not find any points in the range specified');
   disp(['ampl = ' num2str(ampl)]);
   return;
end

if ( ~exist('notr') ) 
   R = R - ones(size(R,1),1)*mean(R);
end

if ( exist('big') ) 
   rmaxmax = max(max(R));
   rminmin = min(min(R));
   xmin = 1.1*min(rminmin,-rmaxmax);
   xmax = 1.1*max(rmaxmax,-rminmin);
   ymin = 1.1*min(rminmin,-rmaxmax);
   ymax = 1.1*max(rmaxmax,-rminmin);
   zmin = 1.1*min(rminmin,-rmaxmax);
   zmax = 1.1*max(rmaxmax,-rminmin);
else
   rmin = min(R);
   rmax = max(R);
   xmin = rmin(1);
   xmax = rmax(1);
   ymin = rmin(2);
   ymax = rmax(2);
   zmin = rmin(3);
   zmax = rmax(3);
end

Mv = moviein(N);
H=[];
for j=1:N
   r = R(iCntr(j):iCntr(j+1)-1,:);
   h=plot3(r(:,1),r(:,2),r(:,3),'.b'); hold on
   H=[H h];
   set(h,'MarkerSize',18);
   if ( exist('mov') );  
      axis([xmin xmax ymin ymax zmin zmax])
      xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)'); 
      grid on;
      set(gca,'Box', 'on')
      Mv(:,j) = getframe;  
   else
      set(h,'Color',colourMtrx(1+rem(j-1,3),:));
   end
end

axis([xmin xmax ymin ymax zmin zmax])
xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)'); 
grid on;
set(gca,'Box', 'on')
hold off;

