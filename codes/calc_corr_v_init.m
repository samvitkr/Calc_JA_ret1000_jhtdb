close all
clear
load('../data/ygrid.mat')

Nx=2048;
Nz=1536;
Ny=512;
jcond=53;
%jc=jcond-Ny/2;
tstart=1;
tstep=1;
tend=50;

fn=sprintf("phi_lse_v_j_%03d.mat",jcond);
m=matfile(fn,'Writable',true);
m.nt=0;
m.tstart=tstart;
m.tstep=tstep;
m.tend=tend;
m.jcond=jcond;
m.Nx=Nx;
m.Ny=Ny;
m.Nz=Nz;

m.phivv=zeros(Nz,Nx,Ny/2);
m.phivu=zeros(Nz,Nx,Ny/2);
m.phivw=zeros(Nz,Nx,Ny/2);
m.phivdudx=zeros(Nz,Nx,Ny/2);
m.phivdvdx=zeros(Nz,Nx,Ny/2);
m.phivdwdx=zeros(Nz,Nx,Ny/2);
m.phivdudy=zeros(Nz,Nx,Ny/2);
m.phivdvdy=zeros(Nz,Nx,Ny/2);
m.phivdwdy=zeros(Nz,Nx,Ny/2);
m.phivdudz=zeros(Nz,Nx,Ny/2);
m.phivdvdz=zeros(Nz,Nx,Ny/2);
m.phivdwdz=zeros(Nz,Nx,Ny/2);
m.phivvoz=zeros(Nz,Nx,Ny/2);
m.phivwoy=zeros(Nz,Nx,Ny/2);
%tstart=0000;
%tend=100000;
%tstep=1000;
%nf=(tend-tstart)/tstep+1;
%%load('lambda_stats.mat')
%for time=tstart:tstep:tend
%        time
%        fvel=sprintf("../data/velfields_%07d.mat",time);
%        m=matfile(fvel);
%	fvelg=sprintf("../data/velgrad_%07d.mat",time);
%        mg=matfile(fvelg);
%	ft=sprintf("../data/transferfields_%07d.mat",time);
%        mt=matfile(ft);
%	viscF=fft2(mt.visc(:,:,Ny/2+1:end))./(Nz*Nx);
%	vozF=fft2(mt.voz(:,:,Ny/2+1:end))./(Nz*Nx);
%	woyF=fft2(mt.woy(:,:,Ny/2+1:end))./(Nz*Nx);
%	vfj=m.vFourier(:,:,jcond);
%	ufj(1,1)=0;
%	vfj(1,1)=0;
%
%	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
%	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
%	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
%
%        phivdudx=phivdudx+conj(vfj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdvdx=phivdvdx+conj(vfj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdwdx=phivdwdx+conj(vfj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdudy=phivdudy+conj(vfj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdvdy=phivdvdy+conj(vfj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdwdy=phivdwdy+conj(vfj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdudz=phivdudz+conj(vfj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdvdz=phivdvdz+conj(vfj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
%        phivdwdz=phivdwdz+conj(vfj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);
%	phivfx=phivfx+conj(vfj).*viscF;
%	phivvoz=phivvoz+conj(vfj).*vozF;
%	phivwoy=phivwoy+conj(vfj).*woyF;
%
%end
%
%phivv=phivv./nf;
%phivu=phivu./nf;
%phivw=phivw./nf;
%
%phivdudx=phivdudx./nf;
%phivdvdx=phivdvdx./nf;
%phivdwdx=phivdwdx./nf;
%phivdudy=phivdudy./nf;
%phivdvdy=phivdvdy./nf;
%phivdwdy=phivdwdy./nf;
%phivdudz=phivdudz./nf;
%phivdvdz=phivdvdz./nf;
%phivdwdz=phivdwdz./nf;
%phivfx=phivfx./nf;
%phivvoz=phivvoz./nf;
%phivwoy=phivwoy./nf;
%
%Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
%Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
%Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
%
%
%Rvdudx=ifft2(phivdudx*(Nz*Nx),'symmetric');
%Rvdvdx=ifft2(phivdvdx*(Nz*Nx),'symmetric');
%Rvdwdx=ifft2(phivdwdx*(Nz*Nx),'symmetric');
%Rvdudy=ifft2(phivdudy*(Nz*Nx),'symmetric');
%Rvdvdy=ifft2(phivdvdy*(Nz*Nx),'symmetric');
%Rvdwdy=ifft2(phivdwdy*(Nz*Nx),'symmetric');
%Rvdudz=ifft2(phivdudz*(Nz*Nx),'symmetric');
%Rvdvdz=ifft2(phivdvdz*(Nz*Nx),'symmetric');
%Rvdwdz=ifft2(phivdwdz*(Nz*Nx),'symmetric');
%
%Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
%Rvvoz=ifft2(phivvoz*(Nz*Nx),'symmetric');
%Rvwoy=ifft2(phivwoy*(Nz*Nx),'symmetric');
%
%fn=sprintf('../data/velgrad_corr_v_j_%03d.mat',jcond);
%mf=matfile(fn,"Writable",true);
%mf.Rvv=Rvv;
%mf.Rvu=Rvu;
%mf.Rvw=Rvw;
%mf.yCheb=yCheb(Ny/2+1:end);
%mf.j=jc;
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
