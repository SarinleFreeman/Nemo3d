% TEST
%
%	Equisurf demo
%
[x,y,z]  = meshgrid([-1:0.2:1],[-1:0.2:1],[-1:0.2:1]);
u = x.*x + y.*y - z.*z;
equisurf(u, [], [-0.5, 0, 0.5]);
show([-37,30]); colormap bone
