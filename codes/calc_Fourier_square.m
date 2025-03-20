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

tstart=1;
tend=10;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time

	ff=sprintf("../data/vel_vort_square_F_%03d.mat",time);
        mf=matfile(ff,'Writable',true)

	fvel=sprintf("../data/velfieldpar_%02d.mat",time);
        m=matfile(fvel,'Writable',true)
	mf.u2F=	fft2(permute(m.ufield(1:Ny,1:Nx,1:Nz).^2,[3 2 1]))./(Nx*Nz);
	mf.v2F=	fft2(permute(m.vfield(1:Ny,1:Nx,1:Nz).^2,[3 2 1]))./(Nx*Nz);
	mf.w2F=	fft2(permute(m.wfield(1:Ny,1:Nx,1:Nz).^2,[3 2 1]))./(Nx*Nz);
	clear m

	fvo=sprintf("../data/vort_%03d.mat",time);
        mvo=matfile(fvo,'Writable',true)	
	mf.ox2F= fft2(permute(mvo.omegax(1:Ny,1:Nx,1:Nz).^2,[3 2 1]))./(Nx*Nz);
        mf.oy2F= fft2(permute(mvo.omegay(1:Ny,1:Nx,1:Nz).^2,[3 2 1]))./(Nx*Nz);
        mf.oz2F= fft2(permute(mvo.omegaz(1:Ny,1:Nx,1:Nz).^2,[3 2 1]))./(Nx*Nz);
	clear mvo
	clear mf
end

