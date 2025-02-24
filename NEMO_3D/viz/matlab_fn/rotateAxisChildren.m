function [Mv] = rotateAxisChildren(Nth, usr)

%  rotateAxisChildren:  [movie Mv] <-- (int Nth)
%
%==============================================================================
%
%  DESCRIPTION:
%
%    Rotates all children of the current axis one revolution about the
%    z-axis in Nth steps.  Returns a movie (and optionally save in
%    .bmp format).
%
%  INPUT:
%    
%    int Nth:   number of steps in which revolution should occur
%
%    string usr:    this string is evaluated within matlab and is used
%                   to set the following optional variables:
%
%                   bmp:  directory where generated .bmp files will be placed
%
%  OUTPUT:
%
%    movie  Mv:     movie of rotation
%
%  DEPENDENCIES: none
%
%  BUGS:  bmp output is not working on volcanoes.  I'm not sure if this
%         is a library problem or a problem with matlab 6.0.
%                                                          
%==============================================================================

if ( exist('usr') );  eval(usr);  end

ch = get(gca,'Children');
Mv = moviein(Nth);

for j=1:Nth
   theta = 2*pi*j/Nth;
   for k=ch
      rotate(k, [0 0 1], 360/Nth);
   end
   Mv(:,j) = getframe;
   if ( exist('bmp') )
      unix(['mkdir -p ' bmp]);
      fil = [bmp '/' num2str(10*theta*180/pi,'%.0f') '.bmp' ];
      disp(['printing to ' fil]);
      eval(['print -dbmp256 ' fil]);
   end
end
