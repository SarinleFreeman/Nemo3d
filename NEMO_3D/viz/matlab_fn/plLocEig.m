function [h] = plLocEig(dat, ltype, usr)

%  plLocEig:  [h] <-- (dat, ltype, usr)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Plots local bandedges along a thin tube.
%
%  INPUT:
%
%    string or mtrx dat:  if dat is a string:
%                           ascii input file consisting of 8 columns where
%
%                              col 1-3 are positions x, y, and z
%                              col 4   is  an eigenenergy
%                              col 5   is  an integer specifying atom species  
%
%                           This file can be constructed by concatenating 
%                           run3d output files of type
%                           <file>.nd_Ek_proc_*_kx0.000_ky_0.000_kz_0.000_r
%                           and deleting the first three columns
%
%                         if dat is a matrix, it can assume 2 possible forms:
%
%                           (1) dat may be a matrix as described above.
%
%                           (2) dat may be the output of locEig and have
%                               the form [R,E]
%
%    string ltype:  line type.  Possible values include '.', ',', 'o', etc.
%
%    string usr:  string evaluated in matlab; is used to define
%                 variables that are passed to locEig and to prettyPlot
%                 (see documentation for these functions for more details)
%                 Additional parameters that may be set include:
%
%                no_legend:  if this exists, no legend is produced.
%       
%
%  OUTPUT:
%
%    plot* h:    pointer to list of pointers to plots
%
%  EXAMPLES:  
%
%    here is a basic use of this function --
%
%      plLocEig('test.nd_Ek_r', '.', 'drct=''x''; Elim=[-1 1.65]; ...
%                Rlim=[-1 99 -1 99 -1 99]; atomType=''all'';')
%
%    equivalently, we may call plLocEig as follows:
%
%      [R,E]=locEig('test.nd_Ek_r','x',[-1 1.65], [-1 99 -1 99 -1 99], 'all');
%      h=plLocEig([R,E],'.');
%
%  DEPENDENCIES: prettyPlot
%
%==============================================================================

if ( exist('usr') );  
   eval(usr);  
   if (~exist('sz_cell')); sz_cell = 160; end;
   if (~exist('atomType')); atomType = 'all'; end;
   if (~exist('Rlim')); Rlim = [-1 999 -1 999 -1 999]; end;
   if (~exist('Elim')); Elim = [-1 1.65]; end;
else
   usr='';
end

if (exist('drct'))
    [R,E]=locEig(dat, drct, Elim, Rlim, atomType);
else
    R = dat(:,1);
    E = dat(:, 2:size(dat,2));
end

Neig = size(E,2);

colour = [0 0 0; .6 0 0; 0 .6 0; 0 0 .6];
keystr = {'SO','H1','H2','EL'};

for j=1:Neig
   h(j)=plot(R,E(:,j),[ltype 'r']); hold on
   set(h(j), 'Color', colour(mod(j, 4)+1, :))
end

xlabel(['position along x (nm)']);
ylabel('energy (eV)');
if (Neig==4 & ~exist('no_legend')); legend(h, keystr); end

prettyPlot(usr);
