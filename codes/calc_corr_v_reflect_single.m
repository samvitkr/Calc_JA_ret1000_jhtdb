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
xp =single([0:Nx-1]*Lx/(Nx));
%kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  single([0:1:Nz-1]*Lz/(Nz));
load('../data/bsplinedata.mat')
%C0= colmat0 ;
%C1= colmat1 ;
%C2= colmat2 ;
yp = single(yv);
clear colmat0 colmat1 colmat2  yv kk knots

%load('../data/bsplinedata.mat')

%Nx=512;
%Nz=384;
%Ny=220;
jcond=130
jct=Ny-jcond+1
%jc=jcond-Ny/2;

nf=1;
phivv=		single(zeros(Nz,Nx,Ny/2));
phivu=		single(zeros(Nz,Nx,Ny/2));
phivw=		single(zeros(Nz,Nx,Ny/2));

phivdudx=	single(zeros(Nz,Nx,Ny/2));
phivdvdx=	single(zeros(Nz,Nx,Ny/2));
phivdwdx=	single(zeros(Nz,Nx,Ny/2));

phivdudy=	single(zeros(Nz,Nx,Ny/2));
phivdvdy=	single(zeros(Nz,Nx,Ny/2));
phivdwdy=	single(zeros(Nz,Nx,Ny/2));

phivdudz=	single(zeros(Nz,Nx,Ny/2));
phivdvdz=	single(zeros(Nz,Nx,Ny/2));
phivdwdz=	single(zeros(Nz,Nx,Ny/2));

%phivfx=	single(	zeros(Nz,Nx,Ny/2));
phivvoz=	single(zeros(Nz,Nx,Ny/2));
phivwoy=	single(zeros(Nz,Nx,Ny/2));

tstart=1;
tend=10;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
   %     fvel=sprintf("../data/velfieldpar_%02d.mat",time);
   %     m=matfile(fvel)

   %     fvelgx=sprintf("../data/velgradx_%03d.mat",time);
   %     mgx=matfile(fvelgx)

   %     fvelgy=sprintf("../data/velgrady_%03d.mat",time);
   %     mgy=matfile(fvelgy)

   %     fvelgz=sprintf("../data/velgradz_%03d.mat",time);
   %     mgz=matfile(fvelgz)

   %     ft=sprintf("../data/Transfer_%03d.mat",time);
   %     mt=matfile(ft)


	fvel = sprintf("../data/velfieldpar_%02d.mat", time);
    m = memmapfile(fvel, 'Format', 'single', 'Writable', false);
    
    fvelgx = sprintf("../data/velgradx_%03d.mat", time);
    mgx = memmapfile(fvelgx, 'Format', 'single', 'Writable', false);
    
    fvelgy = sprintf("../data/velgrady_%03d.mat", time);
    mgy = memmapfile(fvelgy, 'Format', 'single', 'Writable', false);
    
    fvelgz = sprintf("../data/velgradz_%03d.mat", time);
    mgz = memmapfile(fvelgz, 'Format', 'single', 'Writable', false);
    
    ft = sprintf("../data/Transfer_%03d.mat", time);
    mt = memmapfile(ft, 'Format', 'single', 'Writable', false);
    



	vfj=single(m.Data.vF(:,:,jcond));
	vfj(1,1)=0;
	vfjt=single(m.Data.vF(:,:,jct));
        vfjt(1,1)=0;

		phivv(:)=phivv(:)+conj(vfj).*(m.Data.vF(1:Nz,1:Nx,1:Ny/2));
		phivu(:)=phivu(:)+conj(vfj).*(m.Data.uF(1:Nz,1:Nx,1:Ny/2));
		phivw(:)=phivw(:)+conj(vfj).*(m.Data.wF(1:Nz,1:Nx,1:Ny/2));

		phivv(:)=phivv(:)+flip(conj(vfjt).*(m.Data.vF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivu(:)=phivu(:)-flip(conj(vfjt).*(m.Data.uF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivw(:)=phivw(:)-flip(conj(vfjt).*(m.Data.wF(1:Nz,1:Nx,Ny/2+1:Ny),3));

        	phivdudx(:)=phivdudx(:)+conj(vfj).*(mgx.Data.dudxF(1:Nz,1:Nx,1:Ny/2));
        	phivdvdx(:)=phivdvdx(:)+conj(vfj).*(mgx.Data.dvdxF(1:Nz,1:Nx,1:Ny/2));
        	phivdwdx(:)=phivdwdx(:)+conj(vfj).*(mgx.Data.dwdxF(1:Nz,1:Nx,1:Ny/2));
        	phivdudy(:)=phivdudy(:)+conj(vfj).*(mgy.Data.dudyF(1:Nz,1:Nx,1:Ny/2));
        	phivdvdy(:)=phivdvdy(:)+conj(vfj).*(mgy.Data.dvdyF(1:Nz,1:Nx,1:Ny/2));
        	phivdwdy(:)=phivdwdy(:)+conj(vfj).*(mgy.Data.dwdyF(1:Nz,1:Nx,1:Ny/2));
        	phivdudz(:)=phivdudz(:)+conj(vfj).*(mgz.Data.dudzF(1:Nz,1:Nx,1:Ny/2));
        	phivdvdz(:)=phivdvdz(:)+conj(vfj).*(mgz.Data.dvdzF(1:Nz,1:Nx,1:Ny/2));
        	phivdwdz(:)=phivdwdz(:)+conj(vfj).*(mgz.Data.dwdzF(1:Nz,1:Nx,1:Ny/2));

		phivdudx(:)=phivdudx(:)-flip(conj(vfjt).*(mgx.Data.dudxF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivdvdx(:)=phivdvdx(:)+flip(conj(vfjt).*(mgx.Data.dvdxF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivdwdx(:)=phivdwdx(:)-flip(conj(vfjt).*(mgx.Data.dwdxF(1:Nz,1:Nx,Ny/2+1:Ny),3));

                phivdudy(:)=phivdudy(:)+flip(conj(vfjt).*(mgy.Data.dudyF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivdvdy(:)=phivdvdy(:)-flip(conj(vfjt).*(mgy.Data.dvdyF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivdwdy(:)=phivdwdy(:)+flip(conj(vfjt).*(mgy.Data.dwdyF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                
		phivdudz(:)=phivdudz(:)-flip(conj(vfjt).*(mgz.Data.dudzF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivdvdz(:)=phivdvdz(:)+flip(conj(vfjt).*(mgz.Data.dvdzF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivdwdz(:)=phivdwdz(:)-flip(conj(vfjt).*(mgz.Data.dwdzF(1:Nz,1:Nx,Ny/2+1:Ny),3));

		phivvoz(:)=phivvoz(:)+conj(vfj).*(mt.Data.vozF(1:Nz,1:Nx,1:Ny/2));
		phivwoy(:)=phivwoy(:)+conj(vfj).*(mt.Data.woyF(1:Nz,1:Nx,1:Ny/2));

		phivvoz(:)=phivvoz(:)-flip(conj(vfjt).*(mt.Data.vozF(1:Nz,1:Nx,Ny/2+1:Ny),3));
                phivwoy(:)=phivwoy(:)-flip(conj(vfjt).*(mt.Data.woyF(1:Nz,1:Nx,Ny/2+1:Ny),3));
	%end
clear m mgx mgy mgz mt vfj vfjt
end


N=single((Nx*Nz)/(2*nf));

fn=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true)

for j =1:Ny/2
j
mf.Rvv(:,:,j	)=	ifft2(single((phivv(:,:,j)*N)),'symmetric')	;
mf.Rvu(:,:,j	)=	ifft2(single((phivu(:,:,j)*N)),'symmetric')	;
mf.Rvw(:,:,j	)=	ifft2(single((phivw(:,:,j)*N)),'symmetric')	;

mf.Rvdudx(:,:,j	)=	ifft2(single((phivdudx(:,:,j)*N)),'symmetric');
mf.Rvdvdx(:,:,j	)=	ifft2(single((phivdvdx(:,:,j)*N)),'symmetric');
mf.Rvdwdx(:,:,j	)=	ifft2(single((phivdwdx(:,:,j)*N)),'symmetric');

mf.Rvdudy(:,:,j	)=	ifft2(single((phivdudy(:,:,j)*N)),'symmetric');
mf.Rvdvdy(:,:,j	)=	ifft2(single((phivdvdy(:,:,j)*N)),'symmetric');
mf.Rvdwdy(:,:,j	)=	ifft2(single((phivdwdy(:,:,j)*N)),'symmetric');

mf.Rvdudz(:,:,j	)=	ifft2(single((phivdudz(:,:,j)*N)),'symmetric');
mf.Rvdvdz(:,:,j	)=	ifft2(single((phivdvdz(:,:,j)*N)),'symmetric');
mf.Rvdwdz(:,:,j	)=	ifft2(single((phivdwdz(:,:,j)*N)),'symmetric');

mf.Rvvoz(:,:,j	)=	ifft2(single((phivvoz(:,:,j)*N)),'symmetric');
mf.Rvwoy(:,:,j	)=	ifft2(single((phivwoy(:,:,j)*N)),'symmetric');
end
mf.j=jcond;
