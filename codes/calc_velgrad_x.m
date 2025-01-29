Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz);
load('../data/bsplinedata.mat')
C0= colmat0 ;
C1= colmat1 ;
C2= colmat2 ;
yp = yv;
%nproc=6;
%nzproc=Nz/nproc;
nt=4;
tstart=1;
tend=4;
%kstart=zeros(nproc,1);
%for proc=1:nproc
%   kstart(proc)=(proc-1)*nzproc;
%end

dudx=single(zeros(Ny,Nx,Nz));
dvdx=single(zeros(Ny,Nx,Nz));
dwdx=single(zeros(Ny,Nx,Nz));
d2udx2=single(zeros(Ny,Nx,Nz));




%dudy=single(zeros(Ny,Nx,Nz));
%dvdy=single(zeros(Ny,Nx,Nz));
%dwdy=single(zeros(Ny,Nx,Nz));
%
%dudz=single(zeros(Ny,Nx,Nz));
%dvdz=single(zeros(Ny,Nx,Nz));
%dwdz=single(zeros(Ny,Nx,Nz));

for time=tstart:tend
%	fvel=sprintf("vel_%03d.mat",time)
	fvel=sprintf("../data/velfieldpar_%02d.mat",time)
	mv=matfile(fvel)

	for k =1:Nz
	k	
		ufieldslice=mv.ufield(:,:,k);
		vfieldslice=mv.vfield(:,:,k);
		wfieldslice=mv.wfield(:,:,k);

		fu(:,:)=fft(ufieldslice(:,:).').';
        	fv(:,:)=fft(vfieldslice(:,:).').';
        	fw(:,:)=fft(wfieldslice(:,:).').';


        	dfu(:,:)=(fu(:,:)).*(1i*kx);
        	dfv(:,:)=(fv(:,:)).*(1i*kx);
        	dfw(:,:)=(fw(:,:)).*(1i*kx);
        	d2fu(:,:)=(fu(:,:)).*(-kx.^2);


        	dudx(:,:,k) = ifft(dfu(:,:).').';
        	dvdx(:,:,k) = ifft(dfv(:,:).').';
        	dwdx(:,:,k) = ifft(dfw(:,:).').';

        	d2udx2(:,:,k) = ifft(d2fu(:,:).').';
	end
	fvelg=sprintf("../data/velgradx_%03d.mat",time)
	mvg=matfile(fvelg,'Writable',true);
	mvg.dudx=single(dudx);
	mvg.dvdx=single(dvdx);
	mvg.dwdx=single(dwdx);
	mvg.d2udx2=single(d2udx2);
	
end

