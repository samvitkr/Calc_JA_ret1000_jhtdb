close all
clear

Nx=2048;
Ny=512;
Nz=1536;
jcond=71
jct=Ny-jcond+1

nf=1;
phivv2=		single(zeros(Nz,Nx,Ny/2));
phivu2=		single(zeros(Nz,Nx,Ny/2));
phivw2=		single(zeros(Nz,Nx,Ny/2));

phivox2=	single(zeros(Nz,Nx,Ny/2));
phivoy2=	single(zeros(Nz,Nx,Ny/2));
phivoz2=	single(zeros(Nz,Nx,Ny/2));


tstart=1;
tend=10;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')

j=1:Ny/2;
jcr=j+Ny/2;
for time=tstart:tstep:tend
%for time=1:1
        time
        fvel=sprintf("../data/velfieldpar_%02d.mat",time);
        m=matfile(fvel)
	vfj=single(m.vF(:,:,jcond));
        vfj(1,1)=0;
        vfjt=single(m.vF(:,:,jct));
        vfjt(1,1)=0;
	tic
	clear m

	ff=sprintf("../data/vel_vort_square_F_%03d.mat",time);
        mf=matfile(ff,'Writable',true)

	phivv2(:,:,j)=phivv2(:,:,j)+conj(vfj).*single(mf.v2F(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(mf.v2F(1:Nz,1:Nx,jcr)),3);
        phivu2(:,:,j)=phivu2(:,:,j)+conj(vfj).*single(mf.u2F(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(mf.u2F(1:Nz,1:Nx,jcr)),3);
        phivw2(:,:,j)=phivw2(:,:,j)+conj(vfj).*single(mf.w2F(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(mf.w2F(1:Nz,1:Nx,jcr)),3);


	phivox2(:,:,j)=phivox2(:,:,j)+conj(vfj).*single(mf.ox2F(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(mf.ox2F(1:Nz,1:Nx,jcr)),3);
        phivoy2(:,:,j)=phivoy2(:,:,j)+conj(vfj).*single(mf.oy2F(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(mf.oy2F(1:Nz,1:Nx,jcr)),3);
        phivoz2(:,:,j)=phivoz2(:,:,j)+conj(vfj).*single(mf.oz2F(1:Nz,1:Nx,j))-flip( conj(vfjt).*single(mf.oz2F(1:Nz,1:Nx,jcr)),3);
	clear mf
	toc
end


N=single((Nx*Nz)/(2*nf));

fn=sprintf('../data/corr_v_rms_reflect_j_%03d.mat',jcond);
%fn="test.mat"
mf=matfile(fn,"Writable",true)
%mf=matfile('testmem.mat','Writable',true)
%for 
jw =1:Ny/2;
%jw
	mf.Rvv2(1:Nz,1:Nx,jw)=ifft2(single((phivv2(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivv2;
	mf.Rvu2(1:Nz,1:Nx,jw)=ifft2(single((phivu2(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivu2;
	mf.Rvw2(1:Nz,1:Nx,jw)=ifft2(single((phivw2(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivw2;

        mf.Rvox2(1:Nz,1:Nx,jw)= ifft2(single((phivox2(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivox2;
	mf.Rvoy2(1:Nz,1:Nx,jw)=	ifft2(single((phivoy2(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivoy2;
	mf.Rvoz2(1:Nz,1:Nx,jw)=	ifft2(single((phivoz2(1:Nz,1:Nx,jw)*N)),'symmetric');clear phivoz2;
%end
mf.j=jcond;
toc
