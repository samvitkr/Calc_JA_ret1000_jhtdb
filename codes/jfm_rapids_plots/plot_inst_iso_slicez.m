close all
clear all
jcond=105;
counter=5;
 fvgp=sprintf('../../data/conditionalp_dudyinst_jcond_%03d_counter_%03d.mat',jcond,counter)



 %fvgp=sprintf("../../data/conditionalp_jcond_inst_%03d_%02d.mat",jcond,counter);
% fvgn=sprintf("../../data/conditionaln_jcond_inst_%03d_%02d.mat",jcond,counter);
 mp=matfile(fvgp)
% mn=matfile(fvgn)



load('../../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
ut=0.0499;
dnu=1.0006e-3;
jcond=105;
xp=ret*(lx*[0:nx-1]/nx-lx/2);
zp=ret*(lz*[0:nz-1]/nz-lz/2);
yp=ret*(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;

yc=yv(jcond)+1;
ycp=yp(jcond);

[nzz, nxx, nyy]=size(mp.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp=xp(itarget-wxx:itarget+wxx);
zp=zp(ktarget-wzz:ktarget+wzz);
% l1=min(m1.lambda2(1:end,1:end,jcond),[],'all');
% l2=min(m2.lambda2(1:end,1:end,jcond),[],'all');
% val=0.025;
l1=-ut^2/yc^2;
l2=-ut^2/yc^2;
val=-30;
val=-20;
[X,Z,Y]=(meshgrid(xp,zp,yp));

x=permute(X,[2 1 3]);
z=permute(Z,[2 1 3]);
y=permute(Y,[2 1 3]);

 ozd=permute(mp.dvdx-mp.dudy,[2 1 3]).*dnu./ut;
 ld=permute(mp.lambda2,[2 1 3]);

 % ozu=permute(mn.dvdx-mn.dudy,[2 1 3]).*dnu./ut;
 % lu=permute(mn.lambda2,[2 1 3]);

xl1=-200;
xl2=200;
yl1=-300;
yl2=300;
zl1=50;
zl2=300;
cl1=-0.2;
cl2=0.2;
x1=100;
y1=100;


% xl1=-100;
% xl2=100;
% yl1=-100;
% yl2=100;
% zl1=0;
% zl2=150;
% cl1=-0.8;
% cl2=0.8;
% x1=100;
% y1=100;

width=500;
height=550;
f=figure('OuterPosition',[x1,y1,width,height]);
t=tiledlayout(2,1);
t.TileSpacing="tight";
t.Padding="tight";

nexttile
hold on
isosurface(z,x,y,ld,val,ozd)
 % scatter3(0,0,ycp,50,'green','filled')
 hold off
% 
 colormap redblue
 % clim([-1 1])

axis equal
view(45,45)
view(75,10)
zlim([zl1 zl2])
xlim([xl1 xl2])
ylim([yl1 yl2])
clim([cl1 cl2])

xlabel('z^+')
ylabel('x^+')
zlabel('y^+')
text(-250,-410,320,'(a)','FontSize',12)



fn=sprintf('../../data/conditional_zslice_inst_j_%03d_%02d.mat',jcond,counter)
load(fn)
Y=Y./dnu;
yc=Y(jcond,2);
X=X./dnu;
ss=size(X);
mdz=(ss(2)+1)/2;

Lcid=lcid.*(sign(ozd));
[vald ilc]=max(ozd,[],"all")

nexttile
hold on
pcolor(X,Y,ozd./(ut/dnu))
contour(X,Y,Lcid.*(Lcid>0),[1:10],':y','LineWidth',1.5)
contour(X,Y,Lcid.*(Lcid<0),[-10:-1],'-y','LineWidth',1)
% annotation('rectangle',[0.25,0.25,0.1,0.1],'Color','k')
annotation('rectangle',[0.25,0.33,0.135,0.063],'Color','k')
annotation('rectangle',[0.55,0.12,0.13,0.14],'Color','k')
text(-120,238,'(I)')
text(80,160,'(II)')
hold off
shading interp
axis equal
vald=round(vald);
 clim([-0.2 0.2])
ylabel('y^+')
xlabel('x^+')
set(gca, 'TickDir', 'out')%, 'XTickLabel', [])
xticks([-200:100:200])
xlim([-200 200])
ylim([50 300])
% yline(yc)
% xline(0)
c=colorbar;
colormap redblue
ylabel(c,'\omega_z^+')
text(-300,270,'(b)','FontSize',12)
% set(gca,'FontSize',11)


 figname=sprintf("isosurf_inst_slice_jcond_%03d.fig",jcond)
saveas(f,figname)
