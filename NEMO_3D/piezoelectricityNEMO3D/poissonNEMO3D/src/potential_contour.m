clc;
clear all;

x = ones(59);
y = ones(59);
z = ones(59,59);
%
load x_axisstrain3d;
load y_axisstrain3d;
load potential_contour_xy_strain3d;
x = x_axisstrain3d(:,1)*1.E9;
y = y_axisstrain3d(:,1)*1.E9;
z(:)=potential_contour_xy_strain3d(:,1);

figure('color','white');
[C,h]=contour(x,y,z,30);
h=surf(x,y,z);

set(h,'LineWidth',1);
set(gcf,'colormap',jet);
axis tight;
axis('off');
set(gca,'color',[1 1 1]*0.9);
set(gca,'fontsize',12);
set(gca,'layer','top');
set(gca,'linewidth',1);
view([0 0 90]);
%grid on
shading interp;
xlabel('DISTANCE, x [nm]', 'fontsize', 14);
ylabel('DISTANCE, y [nm]', 'fontsize', 14);
zlabel('POTENTIAL [meV]', 'fontsize', 14);
saveas(h,'potential_contour_xy.fig');
print potential_contour_xy.eps -deps;
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1 = ones(55);
z1 = ones(55,59);
load z_axisstrain3d;
y1 = z_axisstrain3d(:,1)*1.E9;
load potential_contour_xz_strain3d;
z1(:)=potential_contour_xz_strain3d(:,1);
figure('color','white');
[C,h]=contour(x,y1,z1,30);
h=surf(x,y1,z1);

set(h,'LineWidth',1);
set(gcf,'colormap',jet);
axis tight;
axis('off');
set(gca,'color',[1 1 1]*0.9);
set(gca,'fontsize',12);
set(gca,'layer','top');
set(gca,'linewidth',1);
view([0 0 90]);
%grid on
shading interp;
xlabel('DISTANCE, x [nm]', 'fontsize', 14);
ylabel('DISTANCE, z [nm]', 'fontsize', 14);
zlabel('POTENTIAL [meV]', 'fontsize', 14);
saveas(h,'potential_contour_xz.fig');
print potential_contour_xz.eps -deps;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load potential_faiXYstrain3d_1;
load potential_faiXmYstrain3d_1;
load potential_faiXYstrain3d;
load potential_faiXmYstrain3d;
x = potential_faiXYstrain3d_1(:,1);
y = potential_faiXYstrain3d_1(:,2);
xx = sqrt(x.*x + y.*y);

yy1 = potential_faiXYstrain3d_1(:,3)*1000;
yy2 = potential_faiXmYstrain3d_1(:,3)*1000;
yy3 = potential_faiXYstrain3d(:,3)*1000;
yy4 = potential_faiXmYstrain3d(:,3)*1000;

figure('color','white');
plot(xx,yy1,xx,yy2,xx,yy3,xx,yy4);
set(h,'LineWidth',1);
%set(gcf,'colormap',jet);
axis tight;
%axis('off');
set(gca,'color',[1 1 1]*0.9);
set(gca,'fontsize',12);
set(gca,'layer','top');
set(gca,'linewidth',1);
%view([0 0 90]);
%grid on
shading interp;
xlabel('DISTANCE [nm]', 'fontsize', 14);
ylabel('PIEZOELECTRIC POTENTIAL [mV]', 'fontsize', 14);
%zlabel('POTENTIAL [meV]', 'fontsize', 14);
saveas(h,'potential_XY_comparison.fig');
print potential_XY_comparison.eps -deps;