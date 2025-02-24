function plPsi2_slice(wf, pltType, usr)

%  plPsi2_slice:  [] <-- (wf, pltType, usr)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Plots a slice of a wave function in various possible formats
%
%  INPUT:
%
%    nx4 matrix wf:  columns of wf are x, y, z, f(x,y,z)
%
%    int pltType:    set of integers indicating what type of figure to plot.
%                    Possible set members include:
%
%                       1:  b/w mesh plot
%                       2:  color 2D plot
%                       3:  greyscale 2D plot
%                       4:  contour plot
%
%    string usr:    this string is evaluated within matlab and is used
%                   to set the following optional variables:
%
%                   fil:  if defined, print to file, 
%                         if not empty, treat as a base name
%                    
%                   limits = a 4-vector of the form 
%                            [r1_min, r1_max, r2_min, r2_max]
%                            used to sepcify the domain of the interpolating 
%                            grid.  By default, it is defined to be the 
%                            smallest domain that fits all points in the 
%                            slice of f
%
%
%  EXAMPLES:
%
%    color 2D plot:  plPsi2_slice(wf, 2, 'z=3.5; fil=''wf'';');
%
%
%  DEPENDENCIES:  prettyPlot, slice_and_interp(restrict)
%
%==============================================================================

if ( exist('usr') );  eval(usr); else; usr=''; end

if ( ~exist('pltType') );  pltType=[1 2 3 4];  end

if ( max(pltType) > 1)
   usr_loc = [usr ' dmesh=0.06;'];
   [r1, r2, f12] = slice_and_interp(wf, usr_loc);
end

if (ismember(1, pltType))
   usr_loc = [usr ' dmesh=0.4;'];
   [r1_msh, r2_msh, f12_msh] = slice_and_interp(wf, usr_loc);

   figure; set(gcf,'Position',[5 730 400 400]);
   mesh(r1_msh,r2_msh,f12_msh);
   set(gcf,'Colormap',~white)
   view(-23,23); 
   if ( exist('limits') ); axis(limits); end
   xlabel('x (nm)'); ylabel('y (nm)'); zlabel('|\Psi({\bf{r}})|^2 (au)');
   if (exist('fil'))
      usr_loc = [usr 'bx=1; gr=1; sh=''square''; ps=[''' fil '_msh''];'];
      prettyPlot(usr_loc);
      eval(['print -djpeg90 ' fil '_msh.jpg']);
   else
      prettyPlot('bx=1; gr=1;');
   end
end
   
if (ismember(2, pltType))
   figure; set(gcf,'Position',[5 730 400 400]);
   mesh(r1,r2,f12); shading interp
   set(gcf,'Colormap',jet);
   view(0,90); 
   xlabel('x (nm)'); ylabel('y (nm)');
   if ( exist('limits') ); axis(limits);  end
   if (exist('fil'))
      usr_loc = [usr 'bx=1; gr=1; sh=''square''; ps=[''' fil '_clr''];'];
      prettyPlot(usr_loc);
      eval(['print -djpeg90 ' fil '_clr.jpg']);
   else
      prettyPlot('bx=1; gr=1;');
   end
end
   
if (ismember(3, pltType))
   figure; set(gcf,'Position',[5 730 400 400]);
   surf(r1,r2,f12);; shading interp
   set(gcf,'Colormap',1-gray);
   view(0,90); 
   if ( exist('limits') ); axis(limits); end
   xlabel('x (nm)'); ylabel('y (nm)'); 
   if (exist('fil'))
      usr_loc = [usr 'bx=1; gr=1; sh=''square''; ps=[''' fil '_grey''];'];
      prettyPlot(usr_loc);
      eval(['print -djpeg90 ' fil '_grey.jpg']);
   else
      prettyPlot('bx=1; gr=1;');
   end
end
   
if (ismember(4, pltType))
   figure; set(gcf,'Position',[5 730 400 400]);
   contour(r1,r2,f12,15, 'k');
   
   xlabel('x (nm)'); ylabel('y (nm)');
   if ( exist('limits') ); axis(limits);  end
   if (exist('fil'))
      usr_loc = [usr 'bx=1; gr=0; sh=''square''; ps=[''' fil '_cntr''];'];
      prettyPlot(usr_loc);
      eval(['print -djpeg90 ' fil '_cntr.jpg']);
   else
      prettyPlot('bx=1; gr=0;');
   end
end
   
