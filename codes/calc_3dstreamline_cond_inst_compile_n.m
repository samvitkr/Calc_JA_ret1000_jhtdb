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
xg=ret*(lx*[0:nx-1]/nx-lx/2);
zg=ret*(lz*[0:nz-1]/nz-lz/2);
yp=ret*(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;
jcset=[47 54 71 105 130];
x1=100;
y1=-500;
width=600;
%height=600;
height=800;
fd=figure;

fd.Position=[x1 y1 width height ];
%t=tiledlayout(3,5);
t=tiledlayout(5,3);

t.TileSpacing="compact";
t.Padding="compact";
%jcond=130;
labels=arrayfun(@(k) ['(' char('a'+k-1) ')'],1:15,'UniformOutput',false);
for j=1:5
jcond=jcset(j);
yc=yv(jcond)+1;
ys=yp(jcond)
fvgp=sprintf("../data/conditionaln_jcond_inst_%03d_01.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_inst_%03d_02.mat",jcond);
fvgq=sprintf("../data/conditionaln_jcond_inst_%03d_03.mat",jcond);

m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
m3=matfile(fvgq,'Writable',true);

[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp=xg(itarget-wxx:itarget+wxx);
zp=zg(ktarget-wzz:ktarget+wzz);

val=200;
[X,Z,Y]=(meshgrid(xp,zp,yp));

x=permute(X,[2 1 3]);
z=permute(Z,[2 1 3]);
y=permute(Y,[2 1 3]);

oxd=permute(m1.dwdy-m1.dvdz,[2 1 3]);
oyd=permute(m1.dudz-m1.dwdx,[2 1 3]);
ozd=permute(m1.dvdx-m1.dudy,[2 1 3]);

oxu=permute(m2.dwdy-m2.dvdz,[2 1 3]);
oyu=permute(m2.dudz-m2.dwdx,[2 1 3]);
ozu=permute(m2.dvdx-m2.dudy,[2 1 3]);

oxq=permute(m3.dwdy-m3.dvdz,[2 1 3]);
oyq=permute(m3.dudz-m3.dwdx,[2 1 3]);
ozq=permute(m3.dvdx-m3.dudy,[2 1 3]);

%%
 [startZ,startX,startY]=meshgrid( 0 ,0,ys);

vertsvd = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
vertsv2d = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);
vertsd = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
verts2d = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);

vertsvu = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
vertsv2u = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);
vertsu = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
verts2u = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);

vertsvq = stream3(z,x,y,  ozq, oxq, oyq,startZ,startX,startY);
vertsv2q = stream3(z,x,y,-ozq,-oxq,-oyq,startZ,startX,startY);
vertsq = stream3(z,x,y,   ozq, oxq, oyq,startZ,startX,startY);
verts2q = stream3(z,x,y, -ozq,-oxq,-oyq,startZ,startX,startY);

% fd=figure;
% fd.Position=[x1 y1 width height];
%%
nexttile
hold on;

lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')

l=streamline(vertsd);
set(l, 'Color', 'k'); 
set(l,'LineWidth',1);
l2=streamline(verts2d);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',1);

colormap redblue
axis equal
view(0,0);

set(gca,'FontSize',9)
xlabel('$z^+$','interpreter','latex','FontSize',11)
%ylabel('$x^+$','interpreter','latex','FontSize',11)
zlabel('$y^+$','interpreter','latex','FontSize',11)
grid on

switch j
    case 1
        xlim([-150 150])
        zlim([0 150])
    case 2
        xlim([-150 150])
        zlim([0 200])
    case 3
        xlim([-150 150])
        zlim([0 400])
    case 4
        xlim([-200 200])
        zlim([0 500])
    case 5
        xlim([-400 400])
        zlim([0 710])
end
title( labels(3*(j-1)+1 ))


%%
nexttile
hold on;

camlight;                      % Add lighting
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')


l=streamline(vertsu);
set(l, 'Color', 'k'); 
set(l,'LineWidth',1);
l2=streamline(verts2u);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',1);
colormap redblue
axis equal

view(0,0);

set(gca,'FontSize',9)
xlabel('$z^+$','interpreter','latex','FontSize',11)
%ylabel('$x^+$','interpreter','latex','FontSize',11)
zlabel('$y^+$','interpreter','latex','FontSize',11)
grid on

switch j
    case 1
        xlim([-150 150])
        zlim([0 150])
    case 2
        xlim([-150 150])
        zlim([0 200])
    case 3
        xlim([-150 150])
        zlim([0 400])
    case 4
        xlim([-200 200])
        zlim([0 500])
    case 5
        xlim([-400 400])
        zlim([0 710])
end
title( labels(3*(j-1)+2 ))
%%
nexttile
hold on;

camlight;                      % Add lighting
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,50,'green','filled')


l=streamline(vertsq);
set(l, 'Color', 'k'); 
set(l,'LineWidth',1);
l2=streamline(verts2q);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',1);
colormap redblue
axis equal

view(0,0);

set(gca,'FontSize',9)
xlabel('$z^+$','interpreter','latex','FontSize',11)
%ylabel('$x^+$','interpreter','latex','FontSize',11)
zlabel('$y^+$','interpreter','latex','FontSize',11)
grid on

switch j
    case 1
        xlim([-150 150])
        zlim([0 150])
    case 2
        xlim([-150 150])
        zlim([0 200])
    case 3
        xlim([-150 150])
        zlim([0 400])
    case 4
        xlim([-200 200])
        zlim([0 500])
    case 5
        xlim([-400 400])
        zlim([0 710])
end
title( labels(3*(j) ) )
end
saveas(fd,'instantlines_n.fig')
exportgraphics(fd,'instantlines_n.eps','BackgroundColor','white')