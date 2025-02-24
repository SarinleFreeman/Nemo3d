% EQUISURF	equisurf(M, limits, constant(s))
%
%	Visualisation of the solution to the equation f(x,y,z) = constant,
%	where f is a real function that is given as a three-dimensional
%	matrix M. The solution for well-behaved functions f is a surface, 
%	which is calculated by eyta and returned as four lists of polygons:
%	p3, p4, p5 and p6 that hold the corner coordinates for all triangles,
%	rectangles, five sided polygons etc that together build up the 
%	surface. eyta does not present the surface on the screen; the 
%	function SHOW should be used for that.
%
%	Arguments:
%
%	M is a three dimensional matrix, i.e. if the total number of elements
%	in the y-direction is ny, the element in M at position (x,y,z) is
%	denoted M((x-1)*ny + y, z). Note that y is the 'innermost index',
%	and z is the outermost. This is consistent with the meshgrid function.
%
%	x-dimension is the total number of elements in the x-direction.
%
%	limits is the vector [xmin xmax ymin ymax zmin zmax] that is used to
%	specify which region of the matrix M that shall be visualised. If
%	limits is empty, the whole M is used.
%
%	The fourth argument is a vector which specifies the constants for which
%	the surface is calculated. If several is given, all the surfaces will
%	appear in the same plot.
%
%	Return values:
%
%	eyta creates four global variables p3, p4, p5 and p6. They must not be
%	used for other purposes in the main program. 
%
%	If you have questions, ask Lennart Bengtsson (tfylb). 
	
function eyta(m, limits, constant)
    global p3 p4 p5 p6

    [p3,p4,p5,p6] = polygons(m, limits, constant);


% The functions 'equisurf' and 'show' are written by Lennart Bengtsson, CTH.

% The C code for the mex-functions 'polygons' and 'interp' can be found in the
% directory ~tfylb/Matlab/Grafik

