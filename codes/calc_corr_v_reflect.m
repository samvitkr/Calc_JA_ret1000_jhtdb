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
jct=Ny-jcond+1
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

	ft=sprintf("../data/Transfer_%03d.mat",time);
        mt=matfile(ft)
	vfj=m.vF(:,:,jcond);
	vfj(1,1)=0;
	vfjt=m.vF(:,:,jct);
        vfjt(1,1)=0;

		phivv=phivv+conj(vfj).*m.vF(1:Nz,1:Nx,1:Ny/2);
		phivu=phivu+conj(vfj).*m.uF(1:Nz,1:Nx,1:Ny/2);
		phivw=phivw+conj(vfj).*m.wF(1:Nz,1:Nx,1:Ny/2);

		phivv=phivv+flip(conj(vfjt).*m.vF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivu=phivu-flip(conj(vfjt).*m.uF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivw=phivw-flip(conj(vfjt).*m.wF(1:Nz,1:Nx,Ny/2+1:Ny),3);

        	phivdudx=phivdudx+conj(vfj).*mgx.dudxF(1:Nz,1:Nx,1:Ny/2);
        	phivdvdx=phivdvdx+conj(vfj).*mgx.dvdxF(1:Nz,1:Nx,1:Ny/2);
        	phivdwdx=phivdwdx+conj(vfj).*mgx.dwdxF(1:Nz,1:Nx,1:Ny/2);
        	phivdudy=phivdudy+conj(vfj).*mgy.dudyF(1:Nz,1:Nx,1:Ny/2);
        	phivdvdy=phivdvdy+conj(vfj).*mgy.dvdyF(1:Nz,1:Nx,1:Ny/2);
        	phivdwdy=phivdwdy+conj(vfj).*mgy.dwdyF(1:Nz,1:Nx,1:Ny/2);
        	phivdudz=phivdudz+conj(vfj).*mgz.dudzF(1:Nz,1:Nx,1:Ny/2);
        	phivdvdz=phivdvdz+conj(vfj).*mgz.dvdzF(1:Nz,1:Nx,1:Ny/2);
        	phivdwdz=phivdwdz+conj(vfj).*mgz.dwdzF(1:Nz,1:Nx,1:Ny/2);

		phivdudx=phivdudx-flip(conj(vfjt).*mgx.dudxF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivdvdx=phivdvdx+flip(conj(vfjt).*mgx.dvdxF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivdwdx=phivdwdx-flip(conj(vfjt).*mgx.dwdxF(1:Nz,1:Nx,Ny/2+1:Ny),3);

                phivdudy=phivdudy+flip(conj(vfjt).*mgy.dudyF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivdvdy=phivdvdy-flip(conj(vfjt).*mgy.dvdyF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivdwdy=phivdwdy+flip(conj(vfjt).*mgy.dwdyF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                
		phivdudz=phivdudz-flip(conj(vfjt).*mgz.dudzF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivdvdz=phivdvdz+flip(conj(vfjt).*mgz.dvdzF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivdwdz=phivdwdz-flip(conj(vfjt).*mgz.dwdzF(1:Nz,1:Nx,Ny/2+1:Ny),3);

		phivvoz=phivvoz+conj(vfj).*mt.vozF(1:Nz,1:Nx,1:Ny/2);
		phivwoy=phivwoy+conj(vfj).*mt.woyF(1:Nz,1:Nx,1:Ny/2);

		phivvoz=phivvoz-flip(conj(vfjt).*mt.vozF(1:Nz,1:Nx,Ny/2+1:Ny),3);
                phivwoy=phivwoy-flip(conj(vfjt).*mt.woyF(1:Nz,1:Nx,Ny/2+1:Ny),3);
	%end

end

phivv=0.5*phivv./nf;
phivu=0.5*phivu./nf;
phivw=0.5*phivw./nf;

phivdudx=0.5*phivdudx./nf;
phivdvdx=0.5*phivdvdx./nf;
phivdwdx=0.5*phivdwdx./nf;
phivdudy=0.5*phivdudy./nf;
phivdvdy=0.5*phivdvdy./nf;
phivdwdy=0.5*phivdwdy./nf;
phivdudz=0.5*phivdudz./nf;
phivdvdz=0.5*phivdvdz./nf;
phivdwdz=0.5*phivdwdz./nf;
phivvoz=0.5*phivvoz./nf;
phivwoy=0.5*phivwoy./nf;


fn=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
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
%mf.Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
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
