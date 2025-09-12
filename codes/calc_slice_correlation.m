close all
clear
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
xp = [0:Nx-1]*Lx/(Nx)-Lx/2;
zp=  [0:1:Nz-1]*Lz/(Nz)-Lz/2;
dnu=1.0006e-3;
ut=0.0499;

xp=xp./dnu;
zp=zp./dnu;
load('../data/bsplinedata.mat')
yp = yv;
y=yp(1:Ny/2)';
% nbx=164;
% nbz=180;
% x=xp(Nx/2-nbx:Nx/2+nbx);
% z=zp(Nz/2-nbz:Nz/2+nbz);
[X,Z]=meshgrid(xp,zp);

load('../data/JHTDB_RET1000.mat')
U0=JHTDB_RET1000(end,2);
yp=JHTDB_RET1000(:,1);
y=yp./yp(end);
U=JHTDB_RET1000(:,2)./U0;
uu=JHTDB_RET1000(:,4)./U0.^2;
vv=JHTDB_RET1000(:,5)./U0.^2;
ww=JHTDB_RET1000(:,6)./U0.^2;
load('../data/vort_rms.mat')

oxrms=0.5*(oxrms+flipud(oxrms));
oxrms=oxrms(1:256);
oyrms=0.5*(oyrms+flipud(oyrms));
oyrms=oyrms(1:256);
ozrms=0.5*(ozrms+flipud(ozrms));
ozrms=ozrms(1:256);
%ozrms=sqrt(ozrms.^2+dUdy.^2);
%%
jcond=105;
yc=yp(jcond);
oxrms=oxrms(jcond)
oyrms=oyrms(jcond)
ozrms=ozrms(jcond)
% vrms=sqrt(vv(jcond))
  fr=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
 mr=matfile(fr);

 vvscale=mr.Rvv(1,1,jcond);
vrms=sqrt(vvscale);
 % ylvox=fftshift(mr.Rvdwdy(:,:,jcond)-mr.Rvdvdz(:,:,jcond))./vvscale;
 % ylvoy=fftshift(mr.Rvdudz(:,:,jcond)-mr.Rvdwdx(:,:,jcond))./vvscale;
 ylvoz=fftshift(mr.Rvdvdx(:,:,jcond)-mr.Rvdudy(:,:,jcond))./vvscale;

ozp=max(ylvoz,[],'all');
ozn=min(ylvoz,[],'all');
 ozav=-dUdy(jcond);

 vout=-ozav/ozp
 vin=ozav/-ozn

 oztotp = ozav+ylvoz.*vout;
 oztotn = ozav+ylvoz.*vin;
 %%
%  figure
%  subplot(1,2,1)
%  pcolor(Z,X,oztotp)
%  shading flat
%  axis equal
%  xlim([-300 300])
% ylim([-300 300])
% clim([ozav -ozav])
% subplot(1,2,2)
% pcolor(Z,X,oztotn)
% axis equal
% xlim([-300 300])
% ylim([-300 300])
% shading flat
% clim([ozav -ozav])
% [Xz,Y]=meshgrid(xp,yp);
%  zlvox=fftshift(squeeze(mr.Rvdwdy(1,:,:)-mr.Rvdvdz(1,:,:))./(vvscale),1)';
%  zlvoy=fftshift(squeeze(mr.Rvdudz(1,:,:)-mr.Rvdwdx(1,:,:))./(vvscale),1)';
%  zlvoz=fftshift(squeeze(mr.Rvdvdx(1,:,:)-mr.Rvdudy(1,:,:))./(vvscale),1)';


%  f=sprintf('../data/slicey_cor_j_%03d.mat',jcond);
%  m=matfile(f,'Writable',true);
%  m.rvox=yrvox;
%  m.rvoy=yrvoy;
%  m.rvoz=yrvoz;
% m.X=X;
% m.Z=Z;
% m.dUdy=dUdy;
% 
%   f=sprintf('../data/slicez_cor_j_%03d.mat',jcond);
%  m=matfile(f,'Writable',true);
%  m.rvox=yrvox;
%  m.rvoy=yrvoy;
%  m.rvoz=yrvoz;
%  m.dUdy=dUdy;
%  m.X=Xz;
%  m.Y=Y;

%% Q = -[ 1/2{(du/dx)^2+(dv/dy)^2+(dw/dz)^2} + (du/dz)(dw/dx) + (dv/dz)(dw/dy) + (du/dy)(dvdx)   ]

lvdudx = fftshift(squeeze(mr.Rvdudx(1,:,jcond)))./vvscale;
lvdvdx = fftshift(squeeze(mr.Rvdvdx(1,:,jcond)))./vvscale;
lvdwdx = fftshift(squeeze(mr.Rvdwdx(1,:,jcond)))./vvscale;

lvdudy = fftshift(squeeze(mr.Rvdudy(1,:,jcond)))./vvscale;
lvdvdy = fftshift(squeeze(mr.Rvdvdy(1,:,jcond)))./vvscale;
lvdwdy = fftshift(squeeze(mr.Rvdwdy(1,:,jcond)))./vvscale;

lvdudz = fftshift(squeeze(mr.Rvdudz(1,:,jcond)))./vvscale;
lvdvdz = fftshift(squeeze(mr.Rvdvdz(1,:,jcond)))./vvscale;
lvdwdz = fftshift(squeeze(mr.Rvdwdz(1,:,jcond)))./vvscale;
istart=1000;
iend=1050;
xg=xp(istart:iend);
%%
vel=ut.*[-40:0.5:40]';


[Xg,Vg] = meshgrid(xg,vel);

dudx=vel*lvdudx(istart:iend);
dvdx=vel*lvdvdx(istart:iend);
dwdx=vel*lvdwdx(istart:iend);

dudy=vel*lvdudy(istart:iend)+dUdy(jcond);
dvdy=vel*lvdvdy(istart:iend);
dwdy=vel*lvdwdy(istart:iend);

dudz=vel*lvdudz(istart:iend);
dvdz=vel*lvdvdz(istart:iend);
dwdz=vel*lvdwdz(istart:iend);

oz=dvdx-dudy;
% Q = -( 0.5*(dudx.^2+dvdy.^2+dwdz.^2) + dudz.*dwdx + dvdz.*dwdy + dudy.*dvdx);
Q = -( 0.5*(dudx.^2+dvdy.^2) + dudy.*dvdx);


S=zeros(3,3);
O=zeros(3,3);
D=zeros(3,3)

lambda2=dudx.*0;
lsw=dudx;
[n1, n2]=size(dudx);
for i =1:n1
    for j=1:n2
D(1,1)=dudx(i,j);
D(1,2)=dudy(i,j);
D(1,3)=dudz(i,j);

D(2,1)=dvdx(i,j);
D(2,2)=dvdy(i,j);
D(2,3)=dvdz(i,j);

D(3,1)=dwdx(i,j);
D(3,2)=dwdy(i,j);
D(3,3)=dwdz(i,j);

S=0.5*(D+D');
O=0.5*(D-D');
A = S*S + O*O;
		ll = sort(eig(A));
		lambda2(i,j) = ll(2);

        lsw(i,j) = max(abs(imag(eig(D))));

end
end


 checkQoz = (Q>0).*(oz>0);
 checkloz = (lambda2<0).*(oz>0);
checklswoz = (lsw>0).*(oz>0);
 % 
% max(checkQoz)
%% 
close all

f=figure

subplot(1,3,1)
pcolor(Xg,Vg./ut,checkQoz)
yline(vrms./ut)
yline(-vrms./ut)
xline(0)
clim([0 1])
xlabel('r_x^+')
ylabel('v^+')
shading flat
title('Q>0 & \omega_z>0')

subplot(1,3,2)
pcolor(Xg,Vg./ut,checkloz)
yline(vrms./ut)
yline(-vrms./ut)
xline(0)
clim([0 1])
xlabel('r_x^+')
ylabel('v^+')

shading flat
title('\lambda_2<0 & \omega_z>0')

subplot(1,3,3)
pcolor(Xg,Vg./ut,checklswoz)
yline(vrms./ut)
yline(-vrms./ut)
xline(0)
clim([0 1])
xlabel('r_x^+')
ylabel('v^+')
shading flat
title('\lambda_{ci}>0 & \omega_z>0')
colormap cool

figname=sprintf('slicey_coherent_oz_j_%03d.fig',jcond)

saveas(f,figname)
%%


Lvdudx = fftshift(squeeze(mr.Rvdudx(:,:,jcond)))./vvscale;
Lvdvdx = fftshift(squeeze(mr.Rvdvdx(:,:,jcond)))./vvscale;
Lvdwdx = fftshift(squeeze(mr.Rvdwdx(:,:,jcond)))./vvscale;

Lvdudy = fftshift(squeeze(mr.Rvdudy(:,:,jcond)))./vvscale;
Lvdvdy = fftshift(squeeze(mr.Rvdvdy(:,:,jcond)))./vvscale;
Lvdwdy = fftshift(squeeze(mr.Rvdwdy(:,:,jcond)))./vvscale;

Lvdudz = fftshift(squeeze(mr.Rvdudz(:,:,jcond)))./vvscale;
Lvdvdz = fftshift(squeeze(mr.Rvdvdz(:,:,jcond)))./vvscale;
Lvdwdz = fftshift(squeeze(mr.Rvdwdz(:,:,jcond)))./vvscale;

kstart=720;
kend=818;
%%
vcheck=10*ut;

dudxp=vcheck*Lvdudx(kstart:kend,istart:iend);
dvdxp=vcheck*Lvdvdx(kstart:kend,istart:iend);
dwdxp=vcheck*Lvdwdx(kstart:kend,istart:iend);

dudyp=vcheck*Lvdudy(kstart:kend,istart:iend)+dUdy(jcond);
dvdyp=vcheck*Lvdvdy(kstart:kend,istart:iend);
dwdyp=vcheck*Lvdwdy(kstart:kend,istart:iend);

dudzp=vcheck*Lvdudz(kstart:kend,istart:iend);
dvdzp=vcheck*Lvdvdz(kstart:kend,istart:iend);
dwdzp=vcheck*Lvdwdz(kstart:kend,istart:iend);
ozp=dvdxp-dudyp;

Qp = -( 0.5*(dudxp.^2+dvdyp.^2+dwdzp.^2) + dudzp.*dwdxp + dvdzp.*dwdyp + dudyp.*dvdxp);


dudxn=-vcheck*Lvdudx(kstart:kend,istart:iend);
dvdxn=-vcheck*Lvdvdx(kstart:kend,istart:iend);
dwdxn=-vcheck*Lvdwdx(kstart:kend,istart:iend);

dudyn=-vcheck*Lvdudy(kstart:kend,istart:iend)+dUdy(jcond);
dvdyn=-vcheck*Lvdvdy(kstart:kend,istart:iend);
dwdyn=-vcheck*Lvdwdy(kstart:kend,istart:iend);

dudzn=-vcheck*Lvdudz(kstart:kend,istart:iend);
dvdzn=-vcheck*Lvdvdz(kstart:kend,istart:iend);
dwdzn=-vcheck*Lvdwdz(kstart:kend,istart:iend);
ozn=dvdxn-dudyn;
Qn = -( 0.5*(dudxn.^2+dvdyn.^2+dwdzn.^2) + dudzn.*dwdxn + dvdzn.*dwdyn + dudyn.*dvdxn);

%%
close all
[X,Z]=meshgrid(xp(istart:iend),zp(kstart:kend));
figure
subplot(1,2,1)
pcolor(Z,X,ozp.*(Qp.*0+1).*(Qp>0.1))
shading flat
clim([-10 10])
yline(0)
xline(0)
xlabel('r_z^+')
ylabel('r_x^+')
axis equal
xlim([-200 200])
ylim([-200 200])
title('(a) outflow ')
colorbar
subplot(1,2,2)
pcolor(Z,X,ozn.*(Qn.*0+1).*(Qn>0.1))
shading flat
clim([-10 10])
xline(0)
yline(0)
xlabel('r_z^+')
ylabel('r_x^+')
axis equal
xlim([-200 200])
ylim([-200 200])
colormap redblue
title('(b) inflow ')
colorbar
sgtitle('\omega_z')