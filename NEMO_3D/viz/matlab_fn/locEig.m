function [R, E] = locEig(dat, drct, Elim, Rlim, atomType)

%  locEig:  [R, E] <-- (dat, drct, Elim, Rlim, atomType)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Filters out local eigenvalue data, and returns eigenvalues
%    contained within a thin tube.
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
%                         if dat is a matrix:
%                           It must have the same form as the data
%                           described above.
%
%    string drct: Axis of tube.  Possible values include {'x', 'y', 'z'}
%
%    real(2) Elim:  Limits of energy range to examine.  This need to be
%                   chosen very carefully.  Ideally, there should be exactly 8
%                   eigenvalues (degeneracies included) in this range.  
%                   In practice, some fault-tolerance is built in, so
%                   that if there are occasionally 7 or 9 eigenvalues,
%                   the code will still work.
%                   
%    real(6) Rlim:  Dimensions of tube to which data are constrained.
%                   Rlim has the form [xmin xmax ymin ymax zmin zmax]
%
%    string atomType:  atom species to plot; can be 'Ga', 'In', 'As', 
%                      or 'all' (default)
%
%  OUTPUT:
%
%    real(n) R:    position vector
%
%    real(n,m) E:  Local energy eigenvalues [E_1(R) ... E_n(R)] 
%
%  EXAMPLES:  here is a basic use of this function --
%
%    [R, E] = locEig('test.nd_Ek_r','x',[-1 1.65], ...
%                    [-1 99 -1 99 -1 99], 'all')
%
%  DEPENDENCIES: none
%
%==============================================================================

if (~exist('sz_cell')); sz_cell = 160; end;  % default: assume 20-oribital case
if (~exist('atomType')); atomType = 'all'; end;
if (~exist('Rlim')); Rlim = [-1 999 -1 999 -1 999]; end;
if (~exist('Elim')); Elim = [-1 1.65]; end;

if (ischar(dat))
   e0 = load(dat);
else
   e0 = dat;
end

if (size(e0,2)==8)
   e0 = e0(:,:);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% set up column index drct_int

if ( strcmp(drct,'x') )
   drct_int = 1;
elseif ( strcmp(drct,'y') )
   drct_int = 2;
elseif ( strcmp(drct,'z') )
   drct_int = 3;
else
   disp(['Confused by direction ' drct '; Bailing out']); return;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% filter by Energy range
% First find all absolute indices within energy range.  Then, determine 
% the indices within the primitive cell that occur most frequently (at 
% least half the time) and store the absolute indices in a matrix.

i_inErange = find(e0(:,4)>Elim(1) & e0(:,4)<Elim(2)) - 1;
N=hist(mod(i_inErange,sz_cell), [0:sz_cell]);
i_inCell = find( N > .5 * size(e0,1) / sz_cell );
Neig = size(i_inCell,2);

if (Neig ~= 8); disp(['WARNING:  found ' num2str(Neig) ' eigenvalues!']); end

i1 = [0:sz_cell:size(e0,1)-1] * ones(1,Neig) + ...
     ones(size(e0,1)/sz_cell,1)*i_inCell;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% filter by atom type
% Each row in the matrix i1 corresponds to a single primitive cell
% and thus each index within a row points to the same atom.  We therefore
% only need to check on the first column of i1 to determine the atom species.

aType = e0(i1(:,1),5);

if (strcmp(atomType,'all'))
   r1 = [1:size(i1,1)];
elseif (strcmp(atomType,'Ga'))
   r1 = find(aType==10);
elseif (strcmp(atomType,'In'))
   r1 = find(aType==14);
elseif (strcmp(atomType,'As'))
   r1 = find(aType==12);
end

i2 = i1(r1,:);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% filter by position

pos = e0(i2(:,1),[1:3]);
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);

r2 = find( x >= Rlim(1) &  x < Rlim(2) & ...
           y >= Rlim(3) &  y < Rlim(4) & ...
           z >= Rlim(5) &  z < Rlim(6));

i3 = i2(r2, :);
R = pos(r2,drct_int);
R = [R;R];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

E = [];
for j = 2:2:Neig
   E = [E [e0(i3(:,j-1),4); e0(i3(:,j),4)]];
end

