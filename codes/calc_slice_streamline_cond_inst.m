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
height=1000;
fd=figure;

fd.Position=[x1 y1 width height ];
% % t=tiledlayout(3,5);
  t=tiledlayout(5,3);
% % % t=tiledlayout(1,3);
% t.TileSpacing="compact";
% t.Padding="compact";
%jcond=130;
labels=arrayfun(@(k) ['(' char('a'+k-1) ')'],1:15,'UniformOutput',false);
for j=2:2
jcond=jcset(j);
yc=yv(jcond)+1;
ys=yp(jcond)
% fvgp=sprintf("../data/conditionalp_jcond_inst_%03d_01.mat",jcond);
% fvgn=sprintf("../data/conditionalp_jcond_inst_%03d_02.mat",jcond);
% fvgq=sprintf("../data/conditionalp_jcond_inst_%03d_03.mat",jcond);

% fvgp=sprintf("../data/conditionalp_jcond_inst_%03d_01.mat",jcond);
% fvgn=sprintf("../data/conditionalp_jcond_inst_%03d_02.mat",jcond);
% fvgq=sprintf("../data/conditionalp_jcond_inst_%03d_03.mat",jcond);

fvgp=sprintf("../data/conditionaln_jcond_inst_%03d_01.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_inst_%03d_02.mat",jcond);
fvgq=sprintf("../data/conditionaln_jcond_inst_%03d_03.mat",jcond);

m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
m3=matfile(fvgq,'Writable',true);

[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;
kslice=wzz+1;

xp=xg(itarget-wxx:itarget+wxx);
zp=zg(ktarget-wzz:ktarget+wzz);

val=200;
[X,Z,Y]=(meshgrid(xp,zp,yp));
[Xz,Yz]=meshgrid(xp,yp);
[Xy,Zy]=meshgrid(xp,zp);

x=permute(X,[2 1 3]);
z=permute(Z,[2 1 3]);
y=permute(Y,[2 1 3]);
ry=y-ys;
r=sqrt(x.^2+z.^2+ry.^2);

oxd=permute(m1.dwdy-m1.dvdz,[2 1 3]);
oyd=permute(m1.dudz-m1.dwdx,[2 1 3]);
ozd=permute(m1.dvdx-m1.dudy,[2 1 3]);
ld = permute(m1.lambda2,[2 1 3]);

oxu=permute(m2.dwdy-m2.dvdz,[2 1 3]);
oyu=permute(m2.dudz-m2.dwdx,[2 1 3]);
ozu=permute(m2.dvdx-m2.dudy,[2 1 3]);
lu=permute(m2.lambda2,[2 1 3]);

oxq=permute(m3.dwdy-m3.dvdz,[2 1 3]);
oyq=permute(m3.dudz-m3.dwdx,[2 1 3]);
ozq=permute(m3.dvdx-m3.dudy,[2 1 3]);
lq=permute(m3.lambda2,[2 1 3]);

omd=sqrt(oxd.^2+oyd.^2+ozd.^2);
omu=sqrt(oxu.^2+oyu.^2+ozu.^2); 
omq=sqrt(oxq.^2+oyq.^2+ozq.^2); 


sinmp=0.9567;
cosmp=-0.2909;
thetadeg = atan(sinmp/cosmp)/pi*180+90;
rdotmp = x.*cosmp + ry.*(sinmp);
rho = sqrt(r.^2-rdotmp.^2);


rcodx=ry.*ozd-z.*oyd;
rcody=z.*oxd-x.*ozd;
rcodz=x.*oyd-ry.*oxd;

rcodx=rcodx./(rho.*omd);
rcody=rcody./(rho.*omd);
rcodz=rcodz./(rho.*omd);

rcoux=ry.*ozu-z.*oyu;
rcouy=z.*oxu-x.*ozu;
rcouz=x.*oyu-ry.*oxu;

rcoux=rcoux./(rho.*omu);
rcouy=rcouy./(rho.*omu);
rcouz=rcouz./(rho.*omu);

rcoqx=ry.*ozq-z.*oyq;
rcoqy=z.*oxq-x.*ozq;
rcoqz=x.*oyu-ry.*oxu;

rcoqx=rcoqx./(rho.*omq);
rcoqy=rcoqy./(rho.*omq);
rcoqz=rcoqz./(rho.*omq);

omtd=rcodx.*cosmp+rcody.*sinmp;
omtu=rcoux.*cosmp+rcouy.*sinmp;
omtq=rcoqx.*cosmp+rcoqy.*sinmp;

hslice=surf(linspace(-100,100,100),linspace(-100,100,100),zeros(100));
rotate(hslice,[-1,0,0],-thetadeg);
xd=get(hslice,'XData');
yd=get(hslice,'YData');
zd=get(hslice,'Zdata')+ys-5;
delete(hslice)

% rcodx=ry.*ozd-z.*oyd;
% rcody=z.*oxd-x.*ozd;
% rcodz=x.*oyd-ry.*oxd;
% 
% rcodx=rcodx./(r.*omd);
% rcody=rcody./(r.*omd);
% rcodz=rcodz./(r.*omd);
% 
% rcoux=ry.*ozu-z.*oyu;
% rcouy=z.*oxu-x.*ozu;
% rcouz=x.*oyu-ry.*oxu;
% 
% rcoux=rcoux./(r.*omu);
% rcouy=rcouy./(r.*omu);
% rcouz=rcouz./(r.*omu);
% 
% rcoqx=ry.*ozq-z.*oyq;
% rcoqy=z.*oxq-x.*ozq;
% rcoqz=x.*oyq-ry.*oxq;
% 
% rcoqx=rcoqx./r;
% rcoqy=rcoqy./r;
% rcoqz=rcoqz./r;


	ozdslicez= squeeze(m1.dvdx(kslice,:,:)-m1.dudy(kslice,:,:))';
	ozuslicez= squeeze(m2.dvdx(kslice,:,:)-m2.dudy(kslice,:,:))';
	ozqslicez= squeeze(m3.dvdx(kslice,:,:)-m3.dudy(kslice,:,:))';
ozdslicey= squeeze(m1.dvdx(:,:,jcond)-m1.dudy(:,:,jcond));
ozuslicey= squeeze(m2.dvdx(:,:,jcond)-m2.dudy(:,:,jcond));
ozqslicey= squeeze(m3.dvdx(:,:,jcond)-m3.dudy(:,:,jcond));


%%
 % [startZ,startX,startY]=meshgrid( 0 ,0,ys);
  [startZ,startX,startY]=meshgrid( 0 ,[-5:10:5],[ys-5:10:ys+5]);

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
% ff=figure;
% t=tiledlayout(1,3);
% t.TileSpacing="compact";
% t.Padding="compact";
%%
nexttile
hold on;


scatter3(0,0,ys+1,50,'green','filled')

% l=streamline(vertsd);
%slice(z,x,y,dnu*ozd./ut,[0],[],[])
% slice(z,x,y,omtd,[],[],[ys])
slice(z,x,y,omtd,xd,yd,zd)
shading flat
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1);
% l2=streamline(verts2d);
% set(l2, 'Color', 'k'); 
% set(l2,'LineWidth',1);

colormap redblue
axis equal
view(0,90);
 clim([-0.5 0.5])
 % clim([-0.01 0.01])
lighting gouraud
lightangle(45,45)
set(gca,'FontSize',8)
xlabel('$z^+$','interpreter','latex','FontSize',8)
ylabel('$x^+$','interpreter','latex','FontSize',11)
% zlabel('$y^+$','interpreter','latex','FontSize',8)
grid on

switch j
    case 1
        xlim([-100 100])
        zlim([ys-40 ys+40])
        ylim([-100 100])

    case 2
        xlim([-100 100])
zlim([ys-40 ys+40])
          ylim([-100 100])
  case 3
        xlim([-100 100])
zlim([ys-40 ys+40])       
ylim([-100 100])

    case 4
        xlim([-100 100])
zlim([ys-40 ys+40])               
ylim([-100 100])
 
    case 5
        xlim([-100 100])
zlim([ys-40 ys+40])
ylim([-100 100])

end
title( labels(3*(j-1)+1 ))


%%
nexttile
hold on;

camlight;                      % Add lighting

scatter3(0,0,ys+1,50,'green','filled')


% l=streamline(vertsu);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1);
% l2=streamline(verts2u);
% set(l2, 'Color', 'k'); 
% set(l2,'LineWidth',1);

% slice(z,x,y,dnu*ozu./ut,[0],[],[])
% slice(z,x,y,omtu,[],[],[ys])
slice(z,x,y,omtu,xd,yd,zd)

shading flat
colormap redblue
axis equal
 clim([-0.5 0.5])
% clim([-0.01 0.01])

view(45,45)
view(0,90);
lighting gouraud; 
lightangle(45,45)
set(gca,'FontSize',8)
xlabel('$z^+$','interpreter','latex','FontSize',8)
ylabel('$x^+$','interpreter','latex','FontSize',11)
% zlabel('$y^+$','interpreter','latex','FontSize',8)
grid on

switch j
    case 1
        xlim([-100 100])
        zlim([0 500])
        ylim([-100 100])

    case 2
        xlim([-100 100])
        zlim([0 500])

          ylim([-100 100])
  case 3
        xlim([-100 100])
        zlim([0 500])
        ylim([-100 100])

    case 4
        xlim([-100 100])
        zlim([0 500])
               ylim([-100 100])
 
    case 5
        xlim([-100 100])
        zlim([0 500])
                ylim([-100 100])

end
title( labels(3*(j-1)+2 ))
%%
nexttile
hold on;

camlight;                      % Add lighting

scatter3(0,0,ys+1,50,'green','filled')


% l=streamline(vertsq);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1);
% l2=streamline(verts2q);
% set(l2, 'Color', 'k'); 
% set(l2,'LineWidth',1);

% slice(z,x,y,dnu*ozq./ut,[0],[],[])
%slice(z,x,y,omtq,[],[],[ys])
slice(z,x,y,omtq,xd,yd,zd)

shading flat
colormap redblue
 clim([-0.5 0.5])
% clim([-0.01 0.01])

axis equal
view(45,45)
lighting gouraud; 
lightangle(45,45)
 view(0,90);

set(gca,'FontSize',8)
xlabel('$z^+$','interpreter','latex','FontSize',8)
ylabel('$x^+$','interpreter','latex','FontSize',11)
% zlabel('$y^+$','interpreter','latex','FontSize',8)
grid on

switch j
    case 1
        xlim([-100 100])
        zlim([0 500])
        ylim([-100 100])

    case 2
        xlim([-100 100])
        zlim([0 500])

          ylim([-100 100])
  case 3
        xlim([-100 100])
        zlim([0 500])
        ylim([-100 100])

    case 4
        xlim([-100 100])
        zlim([0 500])
               ylim([-100 100])
 
    case 5
        xlim([-100 100])
        zlim([0 500])
                ylim([-100 100])

end
title( labels(3*(j) ) )
end
% saveas(fd,'instantlines_ring_p.fig')
% exportgraphics(fd,'instantlines_p.eps','BackgroundColor','white')
 % saveas(fd,'instant_oz_slice_n.fig')
 % saveas(fd,'instant_omt_slice_p.fig')
  saveas(fd,'instant_omt_slice_n.fig')
