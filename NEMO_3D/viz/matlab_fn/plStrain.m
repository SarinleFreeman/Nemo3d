function plStrain(R, ep, whatToPlot, drct, rn)

%  plStrain:  [] <-- (R, ep, whatToPlot, drct, rn)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Plots strain along a line specified by R.
%    Caveat:  This function has not been well-tested.  Also, the interface
%             to this function is likely to change in the future.
%             Use at your own risk!!!
%
%  INPUT:
%
%    real(n) R:  position along an axis.
%
%    real(3,n) ep:  diagonal components of strain tensor as a function
%                   of position determined by R_line.  For example,
%                   ep_line(:,n) is the diagonal part of the strain 
%                   tensor at R(n).
%
%    int(m) whatToPlot:  a vector that specifies the type of plot to make.
%                        The vector may contain an integer from 1 to 4, where
%
%                             1 plots ep_x in blue
%                             2 plots ep_y in green
%                             3 plots ep_z in red
%                             4 plots tr{ep} in black
%
%    string drct: Axis of line.  Possible values include {'x', 'y', 'z'}
%
%    string rn:  Name of the run.  If not specified, no eps file is generated.
%
%  EXAMPLES:  here is a basic use of this function --
%
%    plStrain(R, ep, [1:3], 'x', 'test');
%
%  DEPENDENCIES: prettyPlot
%
%==============================================================================

%pl_types={'a_{xx}','a_{yy}','a_{zz}','Tr{a}'};
pl_types={'\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}','Tr\{\epsilon\}'};

indx=1;
if (ismember(1,whatToPlot)); hl(indx)=plot(R,ep(1,:),'.b'); indx=indx+1; hold on; end
if (ismember(2,whatToPlot)); hl(indx)=plot(R,ep(2,:),'.g'); indx=indx+1; hold on; end
if (ismember(3,whatToPlot)); hl(indx)=plot(R,ep(3,:),'xr'); indx=indx+1; hold on; end
if (ismember(4,whatToPlot)); hl(indx)=plot(R,sum(ep),'.k-'); indx=indx+1; hold on; end
legend(hl,pl_types(whatToPlot));

set(hl,'MarkerSize',[0.2]);
set(gca,'XLim',[min(R) max(R)]);

if ( exist('rn') )
   disp('rn exists');
   dum = ['ms=6; ps=[''' rn '''];'];
%   dum = ['ms=6; ps=[''' rn '_str_' drct '''];'];
else
   dum='ms=6';
end

xlabel(['Position along ' drct ' (nm)']);
%ylabel('Lattice constant (nm)');
%ylabel('Hydrostatic strain');
ylabel('\epsilon');
prettyPlot(dum);
