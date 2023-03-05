Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nt=3;
tstart=1;
tend=35;

jstart=1;
jend=180;
Nj=jend-jstart+1;
Delta=15;
Deltaz=1.25;
G=zeros(Nj,Nx);
Gz=zeros(Nj,Nz);

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

omegaxf=zeros(Nj,Nx,Nz);
omegayf=zeros(Nj,Nx,Nz);
omegazf=zeros(Nj,Nx,Nz);

load('bsplinedata.mat');
yl=yv(jstart:jend)+1;


for time=tstart:tend
%fgx=sprintf("velgradx_%03d.mat",time);
%fgy=sprintf("velgrady_%03d.mat",time);
%fgz=sprintf("velgradz_%03d.mat",time);

%mgx=matfile(fgx);
%mgy=matfile(fgy);
%mgz=matfile(fgz);

fo=sprintf("vort_%03d.mat",time)
mo=matfile(fo);
ff=sprintf("vortfiltered_%03d.mat",time)
mf=matfile(ff);

for j =1:Nj

        G(j,:)  = exp(-( yl(j).^2*(kx.^2)./Delta^2));
        Gz(j,:)= exp(-( yl(j).^2*(kz.^2)./Deltaz^2));

	omegaxf(j,:,:)=ifft(  fft(squeeze( mo.omegax(j,:,:))).*(G(j,:)'));
	omegayf(j,:,:)=ifft(  fft(squeeze( mo.omegay(j,:,:))).*(G(j,:)'));
	omegazf(j,:,:)=ifft(  fft(squeeze( mo.omegaz(j,:,:))).*(G(j,:)'));

	omegaxf(j,:,:)=ifft(  fft(squeeze( omegaxf(j,:,:))').*(Gz(j,:)'))';
	omegayf(j,:,:)=ifft(  fft(squeeze( omegayf(j,:,:))').*(Gz(j,:)'))';
	omegazf(j,:,:)=ifft(  fft(squeeze( omegazf(j,:,:))').*(Gz(j,:)'))';

end

%mo.omegax=single(zeros(Ny,Nx,Nz));
%mo.omegay=single(zeros(Ny,Nx,Nz));
%mo.omegaz=single(zeros(Ny,Nx,Nz));

mf.omegaxf=omegaxf;
mf.omegayf=omegayf;
mf.omegazf=omegazf;

end
