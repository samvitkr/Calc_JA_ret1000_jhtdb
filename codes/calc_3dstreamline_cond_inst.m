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
fvgp=sprintf("../data/conditionalp_jcond_inst_%03d_01.mat",jcond);
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
% % % % % load(fn)
% % % % % z=Z./dnu;
% % % % % y=Y./dnu;
% % % % % x=X./dnu;
% % % % % %ys=y(1,1,jcond-110);
% % % % % ys=Y(1,1,jcond)
% % % % % ltdm=min(ld(:,:,jcond),[],'all');
% % % % % ltum=min(lu(:,:,jcond),[],'all');
% % % % % % alpha=0.05;
% % % % % % alpha=0.02;
% % % % % % ltd=alpha*ltdm
% % % % % % ltu=alpha*ltum
ltd=-val*ut^2/yc;
ltu=-val*ut^2/yc;
% % % % % % ltu=-0.02;
% % % % % alpha=-0.0008;
% % % % % ltd=alpha/(ys)^2;
% % % % % ltu=alpha/(ys)^2
% % % % % idx=52;
%load('eddyset_d_j_156.mat')
%%
 [startZ,startX,startY]=meshgrid( 0 ,0,ys);

%[startZ,startX,startY]=meshgrid( 0 ,[-0.01:0.0025:0.0]./dnu,ys);

vertsvd = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
vertsv2d = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);
vertsd = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
verts2d = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);

%[startZ,startX,startY]=meshgrid( [-0.01:0.01:0.01]./dnu,0 ,ys);
% [startZ,startX,startY]=meshgrid( 0 ,0,ys);

%[startZ,startX,startY]=meshgrid( 0 ,[0.0:0.0025:0.01]./dnu,ys);

vertsvu = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
vertsv2u = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);
vertsu = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
verts2u = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);

fd=figure;
fd.Position=[x1 y1 width height];
%%
subplot(1,2,1)
hold on;
%
% isosurf=isosurface(z,x,y,ld,ltd*100);
% interpColors = interp3(z, x, y, oxd./ut,...wd.*0, ...
%     isosurf.vertices(:, 1), ...
%     isosurf.vertices(:, 2), ...
%     isosurf.vertices(:, 3));
% p = patch(isosurf);
% p.FaceColor = 'interp';        % Interpolated color
% p.EdgeColor = 'none';          % Remove edges
% p.FaceVertexCData = interpColors;    % Assign interpolated colors
% p.FaceAlpha = 0.7;  
% camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')
% l=streamline(vertsvd);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',2);
% l2=streamline(vertsv2d);
% set(l2, 'Color', 'k'); 
% set(l2,'LineWidth',2);

l=streamline(vertsd);
set(l, 'Color', 'k'); 
set(l,'LineWidth',2);
l2=streamline(verts2d);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',2);

colormap redblue
axis equal
view(45,45);
% % ylim([-100 50])
% % xlim([-100 100])
% % zlim([0 200])
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
%
% isosurf=isosurface(z,x,y,lu,ltu);
% interpColors = interp3(z, x, y, oxu./ut,...wd.*0, ...
%     isosurf.vertices(:, 1), ...
%     isosurf.vertices(:, 2), ...
%     isosurf.vertices(:, 3));
% p = patch(isosurf);
% p.FaceColor = 'interp';        % Interpolated color
% p.EdgeColor = 'none';          % Remove edges
% p.FaceVertexCData = interpColors;    % Assign interpolated colors
% p.FaceAlpha = 0.7;  
camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')


l=streamline(vertsu);
set(l, 'Color', 'k'); 
set(l,'LineWidth',2);
l2=streamline(verts2u);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',2);
colormap redblue
axis equal

% l=streamline(vertsd);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',2);
% l2=streamline(verts2d);
% set(l2, 'Color', 'k'); 
% set(l2,'LineWidth',2);
view(45,45);
% % ylim([-50 100])
% % xlim([-100 100])
% % zlim([0 200])
% ylim([-50 50])
% xlim([-50 50])
% zlim([0 150])
% c=colorbar;
% clim([-100 100])
%ylabel(c,"u'")
set(gca,'FontSize',11)
% ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
%zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(b)')
grid on
%%
fdne=sprintf('eddy_j_instlines_%03d.fig',jcond)
%fdne=sprintf('eddy_j_cond_nl_%03d.fig',jcond)   
saveas(fd,fdne)