Nx=2048;
Ny=512;
Nz=1536;

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
load('bsplinedata.mat');
nu=5e-5;
jstart=1;
jend=180;
Nj=jend-jstart+1;
tstart=1;
tend=35;
%time=2

yl=yv(jstart:jend)+1;
Delta=15;
G=zeros(Nj,Nx);
Deltaz=1.25;
Gz=zeros(Nj,Nz);

for time=tstart:tend
	v = single(zeros(Nj,Nx,Nz));
	w = single(zeros(Nj,Nx,Nz));
	omegaz = single(zeros(Nj,Nx,Nz));
	omegay = single(zeros(Nj,Nx,Nz));
	viscous_x = single(zeros(Nj,Nx,Nz));
	convective_x = single(zeros(Nj,Nx,Nz));
	total = single(zeros(Nj,Nx,Nz));


	fgx=sprintf("velgradx_%03d.mat",time);
	fgy=sprintf("velgrady_%03d.mat",time);
	fgz=sprintf("velgradz_%03d.mat",time);
	
	mgx=matfile(fgx);
	mgy=matfile(fgy);
	mgz=matfile(fgz);
	
	fo=sprintf("vort_%03d.mat",time);
	mo=matfile(fo);
	
	fv=sprintf("vel_%03d.mat",time);
	mv=matfile(fv);
	
%	ft=sprintf("Transfer_%03d.mat",time)
%	mt=matfile(ft,'Writable',true);
	
	v=mv.v(jstart:jend,:,:);
	w=mv.w(jstart:jend,:,:);	
	omegaz=mo.omegaz(jstart:jend,:,:);
	omegay=mo.omegay(jstart:jend,:,:);

%	convective_x=(mv.v).*(mo.omegaz)-(mv.w).*(mo.omegay);
	viscous_x=nu*((mgx.d2udx2(jstart:jend,:,:))+(mgy.d2udy2(jstart:jend,:,:))+(mgz.d2udz2(jstart:jend,:,:)));
%	total=mt.convective_x + mt.viscous_x;

	for j =1:Nj

	G(j,:)	= exp(-( yl(j).^2*(kx.^2)./Delta^2));
	Gz(j,:) = exp(-( yl(j).^2*(kz.^2)./Deltaz^2));
	
	omegaz(j,:,:)=omegaz(j,:,:)-mean(mean(  squeeze(omegaz(j,:,:))  ));

	v(j,:,:)=ifft(  fft(squeeze( v(j,:,:))).*(G(j,:)'));
	w(j,:,:)=ifft(  fft(squeeze( w(j,:,:))).*(G(j,:)'));
	omegaz(j,:,:)=ifft(  fft(squeeze(omegaz(j,:,:))).*(G(j,:)'));
	omegay(j,:,:)=ifft(  fft(squeeze(omegay(j,:,:))).*(G(j,:)'));
	viscous_x(j,:,:)=ifft(  fft(squeeze( viscous_x(j,:,:))).*(G(j,:)'));
		
	v(j,:,:)=ifft(  fft(squeeze( v(j,:,:))').*(Gz(j,:)'))';
	w(j,:,:)=ifft(  fft(squeeze( w(j,:,:))').*(Gz(j,:)'))';
	omegaz(j,:,:)=ifft(  fft(squeeze( omegaz(j,:,:))').*(Gz(j,:)'))';
	omegay(j,:,:)=ifft(  fft(squeeze( omegay(j,:,:))').*(Gz(j,:)'))';
	viscous_x(j,:,:)=ifft(  fft(squeeze( viscous_x(j,:,:))').*(Gz(j,:)'))';
	end
	
	convective_x=v.*omegaz-w.*omegay;
	total=viscous_x+convective_x;




	fl=sprintf("Transfer_filtered_fluc_x15z1p25_%03d",time);
	ml=matfile(fl,'Writable',true);
	ml.convective_x=convective_x;
	ml.viscous_x=viscous_x;
	ml.total=total;
end

