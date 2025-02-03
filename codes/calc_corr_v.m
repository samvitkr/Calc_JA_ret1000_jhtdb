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
jcond=71
%jc=jcond-Ny/2;

nf=1;
phivv=		zeros(Nz,Nx,Ny/2);
phivu=		zeros(Nz,Nx,Ny/2);
phivw=		zeros(Nz,Nx,Ny/2);

phivdudx=	zeros(Nz,Nx,Ny/2);
phivdvdx=	zeros(Nz,Nx,Ny/2);
phivdwdx=	zeros(Nz,Nx,Ny/2);

phivdudy=	zeros(Nz,Nx,Ny/2);
phivdvdy=	zeros(Nz,Nx,Ny/2);
phivdwdy=	zeros(Nz,Nx,Ny/2);

phivdudz=	zeros(Nz,Nx,Ny/2);
phivdvdz=	zeros(Nz,Nx,Ny/2);
phivdwdz=	zeros(Nz,Nx,Ny/2);

%phivfx=		zeros(Nz,Nx,Ny/2);
phivvoz=	zeros(Nz,Nx,Ny/2);
phivwoy=	zeros(Nz,Nx,Ny/2);

tstart=1;
tend=5;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfieldpar_%02d.mat",time);
        m=matfile(fvel)

	fvelgx=sprintf("../data/velgradx_%03d.mat",time);
        mgx=matfile(fvelgx)

	fvelgy=sprintf("../data/velgrady_%03d.mat",time);
        mgy=matfile(fvelgy)

	fvelgz=sprintf("../data/velgradz_%03d.mat",time);
        mgz=matfile(fvelgz)

	ft=sprintf("../data/transferfields_%07d.mat",time);
        mt=matfile(ft)
	%viscF=fft2(mt.visc(:,:,Ny/2+1:end))./(Nz*Nx);
	%vozF=fft2(mt.voz(:,:,Ny/2+1:end))./(Nz*Nx);
	%woyF=fft2(mt.woy(:,:,Ny/2+1:end))./(Nz*Nx);
	
	uF=	fft2(permute(m.ufield(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	vF=	fft2(permute(m.vfield(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	wF=	fft2(permute(m.wfield(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);

	vfj=vF(:,:,jcond);

	dudxF=	fft2(permute(mgx.dudx(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	dvdxF=	fft2(permute(mgx.dvdx(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	dwdxF=	fft2(permute(mgx.dwdx(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	
	dudyF=	fft2(permute(mgy.dudy(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        dvdyF=	fft2(permute(mgy.dvdy(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        dwdyF=	fft2(permute(mgy.dwdy(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	
	dudzF=	fft2(permute(mgz.dudz(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        dvdzF=	fft2(permute(mgz.dvdz(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        dwdzF=	fft2(permute(mgz.dwdz(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
	
	vozF=  	fft2(permute(  mt.voz(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        woyF=  	fft2(permute(  mt.woy(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);
        %viscF=   fft2(permute(  mt.woy(1:Ny/2,1:Nx,1:Nz),[3 2 1]))./(Nx*Nz);	

	%vfourier=fft(vj,[],2)./(Nz);
	%vfj=m.vFourier(:,:,jcond);
	%ufj(1,1)=0;
	%vfj(1,1)=0;
	%for i=1:Nx
		phivv=phivv+conj(vfj).*vF;
		phivu=phivu+conj(vfj).*uF;
		phivw=phivw+conj(vfj).*wF;

        	phivdudx=phivdudx+conj(vfj).*dudxF;
        	phivdvdx=phivdvdx+conj(vfj).*dvdxF;
        	phivdwdx=phivdwdx+conj(vfj).*dwdxF;
        	phivdudy=phivdudy+conj(vfj).*dudyF;
        	phivdvdy=phivdvdy+conj(vfj).*dvdyF;
        	phivdwdy=phivdwdy+conj(vfj).*dwdyF;
        	phivdudz=phivdudz+conj(vfj).*dudzF;
        	phivdvdz=phivdvdz+conj(vfj).*dvdzF;
        	phivdwdz=phivdwdz+conj(vfj).*dwdzF;
		%phivfx=phivfx+conj(vfj).*viscF;
		phivvoz=phivvoz+conj(vfj).*vozF;
		phivwoy=phivwoy+conj(vfj).*woyF;
	%end

end

phivv=phivv./nf;
phivu=phivu./nf;
phivw=phivw./nf;

phivdudx=phivdudx./nf;
phivdvdx=phivdvdx./nf;
phivdwdx=phivdwdx./nf;
phivdudy=phivdudy./nf;
phivdvdy=phivdvdy./nf;
phivdwdy=phivdwdy./nf;
phivdudz=phivdudz./nf;
phivdvdz=phivdvdz./nf;
phivdwdz=phivdwdz./nf;
%phivfx=phivfx./nf;
phivvoz=phivvoz./nf;
phivwoy=phivwoy./nf;


fn=sprintf('../data/velgrad_corr_v_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);

mf.Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
mf.Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
mf.Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
mf.Rvdudx=ifft2(phivdudx*(Nz*Nx),'symmetric');
mf.Rvdvdx=ifft2(phivdvdx*(Nz*Nx),'symmetric');
mf.Rvdwdx=ifft2(phivdwdx*(Nz*Nx),'symmetric');
mf.Rvdudy=ifft2(phivdudy*(Nz*Nx),'symmetric');
mf.Rvdvdy=ifft2(phivdvdy*(Nz*Nx),'symmetric');
mf.Rvdwdy=ifft2(phivdwdy*(Nz*Nx),'symmetric');
mf.Rvdudz=ifft2(phivdudz*(Nz*Nx),'symmetric');
mf.Rvdvdz=ifft2(phivdvdz*(Nz*Nx),'symmetric');
mf.Rvdwdz=ifft2(phivdwdz*(Nz*Nx),'symmetric');
mf.Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
mf.Rvvoz=ifft2(phivvoz*(Nz*Nx),'symmetric');
mf.Rvwoy=ifft2(phivwoy*(Nz*Nx),'symmetric');

%fn=sprintf('../data/velgrad_corr_v_j_%03d.mat',jcond);
%mf=matfile(fn,"Writable",true);
%mf.Rvv=Rvv;
%mf.Rvu=Rvu;
%mf.Rvw=Rvw;
%mf.yCheb=yCheb(Ny/2+1:end);
mf.j=jcond;
%
%
%mf.Rvdudx=Rvdudx;
%mf.Rvdvdx=Rvdvdx;
%mf.Rvdwdx=Rvdwdx;
%mf.Rvdudy=Rvdudy;
%mf.Rvdvdy=Rvdvdy;
%mf.Rvdwdy=Rvdwdy;
%mf.Rvdudz=Rvdudz;
%mf.Rvdvdz=Rvdvdz;
%mf.Rvdwdz=Rvdwdz;
%
%mf.Rvfx=Rvfx;
%mf.Rvvoz=Rvvoz;
%mf.Rvwoy=Rvwoy;
