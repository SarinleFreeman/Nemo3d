function svPlot(fil)

% plLanczProg:  [ptr_fig] <-- (fil)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Saves plot of current figure to file
%
%  INPUT:
%
%    string fil:  name of file
%
%  OUTPUT:
%
%    none
%
%  DEPENDENCIES: none
%
%==============================================================================

% gold copy
eval(['print -depsc ' fil '.eps']);

% for importation into ppt
system(['convert -density 200 -size 25% ' fil '.eps -quality 90 -transparent white ' fil '.png']);
