function [mean_var, std_var, Nsample] = distrib(fil, plotType, usr)

%  distrib:  [mean_var std_var Nsample] <-- (fil, plotType, usr)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Computes statistical information and optionally plots histogram
%    of sample data.  Typically, this function is used to analyze
%    distribution of ground state eigenvalues.  Means and standard
%    deviations are returned as output.
%
%  INPUT:
%
%    string fil:   ascii file consisting of three columns:
%                      (1) tag (usually seed information or concentration)
%                      (2) electron ground state energy
%                      (3) hole ground state energy
%
%    string plotType:  'sngl' -- each plot is made on a separate figure
%                      'sub'  -- plots are made on a single figure
%                      else   -- nothing is plotted
%
%    string usr:     this string is evaluated within matlab and is used
%                    to set the following variables:
%
%                         {numBins, eBnds, hBnds}
%                    
%                    numBins -- specifies number of bins for the histogram
%
%                    eBnds -- a 2-vector specifying min/max electron eigvals
%
%                    hBnds -- a 2-vector specifying min/max electron eigvals
%
%  OUTPUT:
%
%    real(4) mean_var:  a vector of the form
%                       [mean(col1), mean(col2), mean(col3), mean(col2-col3)]
%
%    real(4) std_var:   a vector of the form
%                       [std(col1), std(col2), std(col3), std(col2-col3)]
%
%    int(4) Nsample:    Size of sample space for each column
%
%  EXAMPLES:
%
%    distrib('fil', 'sngl', 'numBins=40; eBnds=[1.31 1.32]; hBnds=[.16 .165];')
%
%  DEPENDENCIES:  none
%
%==============================================================================

if ( exist('usr')); eval(usr); end

if (~exist('plotType')); plotType = 'none'; end
if (~exist('numBins')); numBins = 40; end
if (~exist('eBnds')); eBnds=[0.6 2.0]; end
if (~exist('hBnds')); hBnds=[-1.0 0.5]; end

% read in file and filter out outlyers
dat = load(fil);
indx = find(abs(dat(:,2)) < 10 & abs(dat(:,3)) < 10); 
indx_compl = find(abs(dat(:,2)) > 10 | abs(dat(:,3)) > 10);
dat = dat(indx,:);

row_ok_el = find(dat(:,2) > eBnds(1) & dat(:,2) < eBnds(2));
row_ok_hl = find(dat(:,3) > hBnds(1) & dat(:,3) < hBnds(2));
row_ok_eg = intersect(row_ok_el, row_ok_hl);

x=dat(:,1);
el=dat(row_ok_el,2);
hl=dat(row_ok_hl,3);
eg=dat(row_ok_eg,2)-dat(row_ok_eg,3);

mean_var = [mean(x), mean(el), mean(hl), mean(eg)];
std_var = [std(x), std(el), std(hl), std(eg)];
Nsample = [size(x,1), size(el,1), size(hl,1), size(eg,1)];

if ( strncmp(plotType, 'sngl', min(size(plotType,2),size('sngl',2))) )
   figure; hist(x,numBins);
   xlabel('concentration (x)'); ylabel('samples');
   figure; hist(eg,numBins);
   xlabel('energy gap (eV)'); ylabel('samples');
   figure; hist(el,numBins);
   xlabel('conduction band edge (eV)'); ylabel('samples');
   figure; hist(hl,numBins);
   xlabel('valence band edge (eV)'); ylabel('samples');
elseif ( strncmp(plotType, 'sub', min(size(plotType,2),size('sub',2))) )
   subplot(2,2,1), hist(el,numBins);
   xlabel('conduction band edge (eV)'); ylabel('samples');
   subplot(2,2,2), hist(hl,numBins);
   xlabel('valence band edge (eV)'); ylabel('samples');
   subplot(2,2,3), hist(eg,numBins);
   xlabel('energy gap (eV)'); ylabel('samples');
   subplot(2,2,4), hist(x,numBins);
   xlabel('concentration (x)'); ylabel('samples');
end

%[nc, xc] = hist(dat(:,1),numBins);
%[ne, xe] = hist(dat(:,2),numBins);
%[nh, xh] = hist(dat(:,3),numBins);
