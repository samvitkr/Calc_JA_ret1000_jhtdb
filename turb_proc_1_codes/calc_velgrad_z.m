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
load('bsplinedata.mat')
C0= colmat0 ;
C1= colmat1 ;
C2= colmat2 ;
tstart=36;
tend=38;
yp = yv;
nproc=6;
nzproc=Nz/nproc;
nt=3;
kstart=zeros(nproc,1);
for proc=1:nproc
   kstart(proc)=(proc-1)*nzproc;
end

dudz=single(zeros(Ny,Nx,Nz));
dvdz=single(zeros(Ny,Nx,Nz));
dwdz=single(zeros(Ny,Nx,Nz));
d2udz2=single(zeros(Ny,Nx,Nz));




%dudy=single(zeros(Ny,Nx,Nz));
%dvdy=single(zeros(Ny,Nx,Nz));
%dwdy=single(zeros(Ny,Nx,Nz));
%
%dudz=single(zeros(Ny,Nx,Nz));
%dvdz=single(zeros(Ny,Nx,Nz));
%dwdz=single(zeros(Ny,Nx,Nz));

for time=tstart:tend
	fvel=sprintf("vel_%03d.mat",time)
	mv=matfile(fvel);

	for i =1:Nx
	i	
		ufieldslice=squeeze(mv.u(:,i,:));
		vfieldslice=squeeze(mv.v(:,i,:));
		wfieldslice=squeeze(mv.w(:,i,:));

		fu(:,:)=fft(ufieldslice(:,:).').';
        	fv(:,:)=fft(vfieldslice(:,:).').';
        	fw(:,:)=fft(wfieldslice(:,:).').';


        	dfu(:,:)=(fu(:,:)).*(1i*kz);
        	dfv(:,:)=(fv(:,:)).*(1i*kz);
        	dfw(:,:)=(fw(:,:)).*(1i*kz);
        	d2fu(:,:)=(fu(:,:)).*(-kz.^2);


        	dudz(:,i,:) = ifft(dfu(:,:).').';
        	dvdz(:,i,:) = ifft(dfv(:,:).').';
        	dwdz(:,i,:) = ifft(dfw(:,:).').';

        	d2udz2(:,i,:) = ifft(d2fu(:,:).').';
	end
	fvelg=sprintf("velgradz_%03d.mat",time)
	mvg=matfile(fvelg,'Writable',true);
	mvg.dudz=single(dudz);
	mvg.dvdz=single(dvdz);
	mvg.dwdz=single(dwdz);
	mvg.d2udz2=single(d2udz2);
	
end
