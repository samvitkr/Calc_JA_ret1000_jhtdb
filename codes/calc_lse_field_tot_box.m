close all
clear
%jcset=[116, 135, 187, 198, 205];
jcond=105;
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
%kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx)-Lx/2;
%kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz)-Lz/2;

nbx=164;
nbz=180;

x=xp(Nx/2-nbx:Nx/2+nbx);
z=zp(Nz/2-nbz:Nz/2+nbz);

load('../data/mean_profiles.mat')

Um=reshape(U,[1 1 Ny/2]);
dUdym=reshape(dUdy,[1 1 Ny/2]);
vozm=reshape(vozm,[1 1 Ny/2]);
woym=reshape(woym,[1 1 Ny/2]);

fn=sprintf('../data/lsevp_field_j_%03d.mat',jcond)
m=matfile(fn,'Writable',true)

fnp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
mp=matfile(fnp,'Writable',true);
mp.x=x;
mp.z=z;

fnn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
mn=matfile(fnn,'Writable',true);
mn.x=x;
mn.z=z;


mp.u=Um+	m.ulse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.v=		m.vlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.w=		m.wlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mp.dudx=	m.dudxlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.dvdx=	m.dvdxlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.dwdx=	m.dwdxlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mp.dudy=dUdym+	m.dudylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.dvdy=	m.dvdylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.dwdy=	m.dwdylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mp.dudz=	m.dudzlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.dvdz=	m.dvdzlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.dwdz=	m.dwdzlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mp.voz=vozm+	m.vozlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mp.woy=woym+	m.woylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);



mn.u=Um		-m.ulse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.v=		-m.vlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.w=		-m.wlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mn.dudx=	-m.dudxlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.dvdx=	-m.dvdxlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.dwdx=	-m.dwdxlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mn.dudy=dUdym	-m.dudylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.dvdy=	-m.dvdylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.dwdy=	-m.dwdylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mn.dudz=	-m.dudzlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.dvdz=	-m.dvdzlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.dwdz=	-m.dwdzlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);

mn.voz=vozm	-m.vozlse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
mn.woy=woym	-m.woylse(Nz/2-nbz:Nz/2+nbz,Nx/2-nbx:Nx/2+nbx,1:Ny/2);
