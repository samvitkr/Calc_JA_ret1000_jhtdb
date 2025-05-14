close all
clear
load('../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
xp=ret*(lx*[0:nx-1]/nx-lx/2);
zp=ret*(lz*[0:nz-1]/nz-lz/2);
yp=ret*(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;
jcond=130;
jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_%03d.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp1=xp(itarget-wxx:itarget+wxx);
zp1=zp(ktarget-wzz:ktarget+wzz);
l1=min(m1.lambda2(1:end,1:end,jc),[],'all');
l2=min(m2.lambda2(1:end,1:end,jc),[],'all');
val=0.1;
[X,Z,Y]=(meshgrid(xp1,zp1,yp));

% mt=matfile('velgrad_transfer_flp_0070000.mat')
% m=matfile('lambdaflp_0070000.mat');
%mt=matfile('transferfields_0040000.mat')
%m=matfile('lambda_0040000.mat');
% ml=(mean(mean(m.lambda2,1),2));
% l=m.lambda2;
% lrms=rms(rms(m.lambda2-ml,1),2);

x1=150;
y1=150;
x2=2*450;
y2=350;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,2,1)
m=m1;
%polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
l1=max(nl(1:end,1:end,jc),[],'all');

subplot(1,2,1)
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl./l1,[2 1 3]),val,permute(nl,[2 1 3]))
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl./l1,[2 1 3]),-val,permute(nl,[2 1 3]))

% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),val,permute(nl,[2 1 3]))
% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),-val,permute(nl,[2 1 3]))

scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,45)
%camlight('left')
clim([-0.05 0.05])
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)
%zlim([30 600])
% xlim([-100 100])
% ylim([-100 100])
%zlim([0 200])
subplot(1,2,2)

m=m2;

%polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
l2=min(nl(1:end,1:end,jc),[],'all');

% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl./l2,[2 1 3]),val,permute(nl,[2 1 3]))
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl./l2,[2 1 3]),-val,permute(nl,[2 1 3]))

% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),val,permute(nl,[2 1 3]))
% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),-val,permute(nl,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading flat
lightangle(45,45)
%camlight('left')
%clim([-0.05 0.05] )
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('z')
ylabel('x')
zlabel('y')
view(45,45)

   fdne=sprintf('eddy_j_cond_flux_%03d.fig',jcond)
   saveas(h1,fdne)
% xlim([-100 100])
% ylim([-100 100])
%zlim([0 200])
%zlim([30 600])