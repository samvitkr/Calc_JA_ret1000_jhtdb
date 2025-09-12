clear 
close all

x1=100;
y1=10;
width=600;
height=300;
ut=0.0499;
dnu=1.0006e-3;

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
jcond=71;
yc=yv(jcond)+1;
ys=yp(jcond)
fvgp=sprintf("../data/conditionalp_jcond_inst_%03d_04.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_inst_%03d_01.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp=xp(itarget-wxx:itarget+wxx);
zp=zp(ktarget-wzz:ktarget+wzz);
%l1=min(m1.lambda2(1:end,1:end,jcond),[],'all');
%l2=min(m2.lambda2(1:end,1:end,jcond),[],'all');
val=200;
[X,Z,Y]=(meshgrid(xp,zp,yp));

x=permute(X,[2 1 3]);
z=permute(Z,[2 1 3]);
y=permute(Y,[2 1 3]);

oxd=permute(m1.dwdy-m1.dvdz,[2 1 3]);
oyd=permute(m1.dudz-m1.dwdx,[2 1 3]);
ozd=permute(m1.dvdx-m1.dudy,[2 1 3]);
%ld=permute(m1.lambda2,[2 1 3]);

oxu=permute(m2.dwdy-m2.dvdz,[2 1 3]);
oyu=permute(m2.dudz-m2.dwdx,[2 1 3]);
ozu=permute(m2.dvdx-m2.dudy,[2 1 3]);
%lu=permute(m2.lambda2,[2 1 3]);

%%

ltd=-val*ut^2/yc;
ltu=-val*ut^2/yc;

%%
 [startZ,startX,startY]=meshgrid( 0 ,0,ys);
vertsvd = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
vertsv2d = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);
vertsvu = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
vertsv2u = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);

[startZ,startX,startY]=meshgrid( 0 ,[-0.01:0.0025:0.0025]./dnu,ys);
% [startZ,startX,startY]=meshgrid( 0 ,[-0.005:0.0025:0.0025]./dnu,ys);
% 

vertsd = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
verts2d = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);

%[startZ,startX,startY]=meshgrid( [-0.01:0.01:0.01]./dnu,0 ,ys);
% [startZ,startX,startY]=meshgrid( 0 ,0,ys);

[startZ,startX,startY]=meshgrid( 0 ,[0.0025:0.0025:0.01]./dnu,ys);


vertsu = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
verts2u = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);

fd=figure;
fd.Position=[x1 y1 width height];
%%
subplot(1,2,1)
hold on;

lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')


l=streamline(vertsd);
set(l, 'Color', 'k'); 
set(l,'LineWidth',1.5);
l2=streamline(verts2d);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',1.5);
l=streamline(vertsvd);
set(l, 'Color', 'm'); 
set(l,'LineWidth',2);
l2=streamline(vertsv2d);
set(l2, 'Color', 'm'); 
set(l2,'LineWidth',2);

colormap redblue
axis equal
view(45,45);

ylim([-100 50])
xlim([-100 100])
zlim([0 200])

% ylim([-200 100])
% xlim([-150 150])
% zlim([0 400])

% c=colorbar;
% clim([-100 100])
%ylabel(c,"u'")
set(gca,'FontSize',11)
%ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(a)')
grid on
%%
subplot(1,2,2)
hold on;

camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')


l=streamline(vertsu);
set(l, 'Color', 'k'); 
set(l,'LineWidth',1.5);
l2=streamline(verts2u);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',1.5);

l=streamline(vertsvu);
set(l, 'Color', 'm'); 
set(l,'LineWidth',2);
l2=streamline(vertsv2u);
set(l2, 'Color', 'm'); 
set(l2,'LineWidth',2);
%colormap redblue
axis equal

view(45,45);
ylim([-50 100])
xlim([-100 100])
zlim([0 200])

% ylim([-100 200])
% xlim([-200 200])
% zlim([0 400])

set(gca,'FontSize',11)
% ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
%zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(b)')
grid on
%%
fdne=sprintf('eddy_j_instlines_magenta_%03d.fig',jcond)
%fdne=sprintf('eddy_j_cond_nl_%03d.fig',jcond)   
saveas(fd,fdne)