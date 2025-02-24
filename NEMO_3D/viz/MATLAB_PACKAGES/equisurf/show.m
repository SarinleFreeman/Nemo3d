% SHOW	Show a equipotential surface on the screen.
%	
%	'show' presents the result of a call to equisurf as a shaded surface on
%	the screen. It takes zero to five arguments, and the most general call
%	is show(view_dir, source_dir, params, 'shading', 'line colour') where 
%
%	view_dir is a 3-D unit vector that specifies the viewing direction and
%	points away from the viewer.
%
%	source_dir is a 3-D unit vector that specifies the lightning direction
%	and points away from the light source.
%
%	params is the vector [ambient, diffuse, specular, spread] that sets
%	the relative contributions due to ambient light, diffuse reflection,
%	specular reflection and specular spread. (See help on diffuse and
%	specular.)
%
%	shading is the string 'interp' (default) or 'flat'. 'interp' gives
%	interpolated shading, that usually looks best on the screen. If you
%	intend to print the figure on a printer 'flat' shading is less
%	time consuming than 'interp'.
%
%	line colour is a string that specifies the colour of the edges of the
%	polygons, i.e. 'black' or 'white'. It is also possible to specify the
%	colour with a vector with three elements [R G B]. The default is black.
%	To get black lines on a printer, you should use white colour since the
%	colour of the lines is inverted before printing. 
%	Note: If you want to provide a line colour, you must specify 'shading'
%	as well.
%
%	The view and source directions may alternatively be given as the
%	angles [az, el]. See help on view for more details.
%
%	The default value for view_dir is taken from the current setting of the
%	view axis property. The default setting of the source direction is
%	the view direction rotated 45 degrees in positive direction around the
%	z-axis.
%
%	The default params are [0.55 0.6 0.4 10].
%
%	Note: The command 'view' should not be used to change the viewing
%	direction dircetly, since the colours would not be correct. Use show
%	instead. 
%
%	The best results are obtained with a single-coloured colourmap i.e.
%	'colormap bone'.
%
%	If you have questions, ask Lennart Bengtsson (tfylb).

function show(v, s, p, m, colour)
    global p3 p4 p5 p6

    error(nargchk(0,5,nargin));

% Parse arguments

    sh = 'interp';
    ec = 'black';
    nin = nargin;	% nin = how many of the arguments v, s and p that are
    if nargin == 1,	% provided
	if isstr(v),
	    sh = v;
	    nin = 0;
	end
    elseif nargin == 2,
	if isstr(v),
	    sh = v;
	    ec = s;
	    nin = 0;
	elseif isstr(s),
	    sh = s;
	    nin = 1;
	end
    elseif nargin == 3,
	if isstr(s),
	    sh = s;
	    ec = p;
	    nin = 1;
	elseif isstr(p),
	    sh = p;
	    nin = 2;
	end
    elseif nargin == 4,
	if isstr(p),
	    sh = p;
	    ec = m;
	    nin = 2;
	elseif isstr(m),
	    sh = m;
	    nin = 3;
	else
	    error('The shading must be a string');
	end
    elseif nargin == 5,
	if isstr(m),
	    sh = m;
	    ec = colour;
	    nin = 3;
	else
	    error('The shading must be a string');
	end
    end

    if ~strcmp(sh, 'flat') & ~strcmp(sh, 'interp'),
	error('The shading can only be ''flat'' or ''interp''');
    end

    if nin < 3,
	p = [0.55 0.6 0.4 10];	% Default parameters fetched from surfl.m
    end

% Compute view and source directions

    if nin < 1,
	[az,el] = view;
	az = az * pi / 180 + pi;
	el = -el * pi / 180;
    else
	if length(v) == 3,
	    v = v / norm(v);
	    vr = norm(v(1:2));
	    if vr == 0,
		el = pi/2 * v(3);
		az = 0;
	    else
		el = atan(v(3)/vr);
		az = atan2(v(1),v(2)) + pi/2;
	    end
	else
	    az = v(1) * pi / 180 + pi;
	    el = -v(2) * pi / 180;
	end
    end
    v = [cos(el)*sin(az), -cos(el)*cos(az), sin(el)];

    if nin < 2 
	s = [az*180/pi - 135, -el*180/pi];
    end
    if length(s) == 2,
	saz = s(1) * pi / 180 + pi;
	sel = -s(2) * pi / 180;
	s(1) = cos(sel)*sin(saz);
	s(2) = -cos(sel)*cos(saz);
	s(3) = sin(sel);
    end

    p(1:3) = p(1:3) / sum(p(1:3));		       % Normalize parameters
    [c3,c4,c5,c6] = interp(p3, p4, p5, p6, v, s, p);   % Calculate interpolated
						       % shadings
% Begin to plot

    cax = newplot;
	
    h = [];
    [j,i] = size(p3);
    if i > 0,
	if strcmp(sh,'flat'),
	    c3 = mean(c3);
	end
	i = i / 3; 
	h = patch(p3(:,1:i), p3(:,i+1:2*i), p3(:,2*i+1:3*i), c3); 
    end

    [j,i] = size(p4);
    if i > 0,
	if strcmp(sh,'flat'),
	    c4 = mean(c4);
	end
	i = i / 3;
	h = [h; patch(p4(:,1:i), p4(:,i+1:2*i), p4(:,2*i+1:3*i), c4)]; 
    end

    [j,i] = size(p5);
    if i > 0,
	if strcmp(sh,'flat'),
	    c5 = mean(c5);
	end
	i = i / 3;
	h = [h; patch(p5(:,1:i), p5(:,i+1:2*i), p5(:,2*i+1:3*i), c5)];
    end

    [j,i] = size(p6);
    if i > 0,
	if strcmp(sh,'flat'),
	    c6 = mean(c6);
	end
	i = i / 3;
	h = [h; patch(p6(:,1:i), p6(:,i+1:2*i), p6(:,2*i+1:3*i), c6)];
    end

    cmax = max([max(max(c3)) max(max(c4)) max(max(c5)) max(max(c6))]);
				% Allthough non-standard, this improves the
    caxis([0 cmax]);		% contrast
    view((az-pi)*180/pi, -el*180/pi);

    if ~strcmp(ec, 'black'),
	set (h, 'Edgecolor', ec);
    end

% Author: Lennart Bengtsson, CTH, Sweden
