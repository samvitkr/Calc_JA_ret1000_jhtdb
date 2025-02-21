close all
clear
%jcset=[116, 135, 187, 198, 205];
jcond=41;
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
%kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
load('../data/bsplinedata.mat')
xp = [0:Nx-1]*Lx/(Nx)-Lx/2;
%kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz)-Lz/2;

nbx=60;
nbz=60;

x=xp(Nx/2-nbx:Nx/2+nbx);
z=zp(Nz/2-nbz:Nz/2+nbz);
y=yv(1:256)'+1;


[X,Z,Y]=meshgrid(x,z,y);

fnp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
mp=matfile(fnp);

fnn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
mn=matfile(fnn);

[nz1,nx1,ny1]=size(mp.u);
ncz=floor(nz1/2);
ncx=floor(nx1/2);


fn=sprintf('../data/lse_eddyset_j_%03d.mat',jcond)
m=matfile(fn,'Writable',true)
m.X=permute(X,[2 1 3]);
m.Y=permute(Y,[2 1 3]);
m.Z=permute(Z,[2 1 3]);

m.ud=	permute(mp.u(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.vd=	permute(mp.v(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.wd=	permute(mp.w(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.oxd=	permute(mp.dwdy(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:)-mp.dvdz(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.oyd=	permute(mp.dudz(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:)-mp.dwdx(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.ozd=	permute(mp.dvdx(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:)-mp.dudy(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.ld=	permute(mp.lambda2(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);


m.uu=	permute(mn.u(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.vu=	permute(mn.v(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.wu=	permute(mn.w(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.oxu=	permute(mn.dwdy(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:)-mn.dvdz(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.oyu=	permute(mn.dudz(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:)-mn.dwdx(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.ozu=	permute(mn.dvdx(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:)-mn.dudy(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);
m.lu=	permute(mn.lambda2(ncz-nbz:ncz+nbz,ncx-nbx:ncx+nbx,:),[2 1 3]);





