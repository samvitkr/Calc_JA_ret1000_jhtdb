close all
clear

Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
%kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
%kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz);
load('../data/bsplinedata.mat')
%C0= colmat0 ;
%C1= colmat1 ;
%C2= colmat2 ;
yp = yv;
clear colmat0 colmat1 colmat2  yv kk knots

%load('../data/bsplinedata.mat')

%Nx=512;
%Nz=384;
%Ny=220;
%jcond=71
%jc=jcond-Ny/2;

nf=1;
%phivv=		zeros(Nz,Nx,Ny/2);
%phivu=		zeros(Nz,Nx,Ny/2);
%phivw=		zeros(Nz,Nx,Ny/2);
%
%phivdudx=	zeros(Nz,Nx,Ny/2);
%phivdvdx=	zeros(Nz,Nx,Ny/2);
%phivdwdx=	zeros(Nz,Nx,Ny/2);
%
%phivdudy=	zeros(Nz,Nx,Ny/2);
%phivdvdy=	zeros(Nz,Nx,Ny/2);
%phivdwdy=	zeros(Nz,Nx,Ny/2);
%
%phivdudz=	zeros(Nz,Nx,Ny/2);
%phivdvdz=	zeros(Nz,Nx,Ny/2);
%phivdwdz=	zeros(Nz,Nx,Ny/2);
%
%%phivfx=		zeros(Nz,Nx,Ny/2);
%phivvoz=	zeros(Nz,Nx,Ny/2);
%phivwoy=	zeros(Nz,Nx,Ny/2);

tstart=2;
tend=5;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfieldpar_%02d.mat",time);
        m=matfile(fvel,'Writable',true)

	fvelgx=sprintf("../data/velgradx_%03d.mat",time);
        mgx=matfile(fvelgx,'Writable',true)

	fvelgy=sprintf("../data/velgrady_%03d.mat",time);
        mgy=matfile(fvelgy,'Writable',true)

	fvelgz=sprintf("../data/velgradz_%03d.mat",time);
        mgz=matfile(fvelgz,'Writable',true)

	ft=sprintf("../data/Transfer_%03d.mat",time);
        mt=matfile(ft,'Writable',true)
	%viscF=fft2(mt.visc(:,:,Ny/2+1:end))./(Nz*Nx);
	%vozF=fft2(mt.voz(:,:,Ny/2+1:end))./(Nz*Nx);
	%woyF=fft2(mt.woy(:,:,Ny/2+1:end))./(Nz*Nx);
	
	m.uF=	fft2(permute(m.ufield(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	m.vF=	fft2(permute(m.vfield(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	m.wF=	fft2(permute(m.wfield(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);

	%vfj=vF(:,:,jcond);

	mgx.dudxF=	fft2(permute(mgx.dudx(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	mgx.dvdxF=	fft2(permute(mgx.dvdx(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	mgx.dwdxF=	fft2(permute(mgx.dwdx(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	
	mgy.dudyF=	fft2(permute(mgy.dudy(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        mgy.dvdyF=	fft2(permute(mgy.dvdy(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        mgy.dwdyF=	fft2(permute(mgy.dwdy(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	
	mgz.dudzF=	fft2(permute(mgz.dudz(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        mgz.dvdzF=	fft2(permute(mgz.dvdz(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        mgz.dwdzF=	fft2(permute(mgz.dwdz(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	
	mt.vozF=  	fft2(permute(  mt.voz(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        mt.woyF=  	fft2(permute(  mt.woy(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        mt.viscF=   fft2(permute(  mt.visc(1:Ny,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);	


end

