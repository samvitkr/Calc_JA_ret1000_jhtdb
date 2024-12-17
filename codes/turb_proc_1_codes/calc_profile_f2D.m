Nx=2048;
Ny=512;
Nz=1536;

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

[Kx,Kz]=meshgrid(kx,kz);
Kx=Kx';
Kz=Kz';


%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
%load('bsplinedata.mat');
%load('kz_vals.mat');
load('filter_2D_flower.mat')
nu=5e-5;
%jstart=jloc(3);
%jend=jloc(21);
Nj=jend-jstart+1;
tstart=1;
tend=38;
%time=2

%yl=yv(jstart:jend)+1;
%Delta=15;
%G=zeros(Nj,Nx);
%Deltaz=1.25;
%Gz=zeros(Nj,Nz);
%
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
%	viscous_x=nu*((mgx.d2udx2(jstart:jend,:,:))+(mgy.d2udy2(jstart:jend,:,:))+(mgz.d2udz2(jstart:jend,:,:)));
%	total=mt.convective_x + mt.viscous_x;

	for j =1:Nj

%	jl=j-1+jstart;
%	hcx=abs(Kx)<kztaylor(jl);
%        hcz=abs(Kz)<kztaylor(jl);
%        f=exp( - ((  abs(Kx)+mf(jl).*(abs(Kz) ))).^2/axf(jl)^2+((abs(Kz) )-mf(jl).*abs(Kx)).^2/azf(jl)^2);

 %       f=f.*hcx.*hcz;

%	G(j,:)	= exp(-( yl(j).^2*(kx.^2)./Delta^2));
%	Gz(j,:) = exp(-( yl(j).^2*(kz.^2)./Deltaz^2));
	
%	omegaz(j,:,:)=omegaz(j,:,:);%-mean(mean(  squeeze(omegaz(j,:,:))  ));

%	v(j,:,:)=ifft2( 	fft2(squeeze( v(j,:,:))).*fbp(:,:,j)	,'symmetric');
%	w(j,:,:)=ifft2( 	fft2(squeeze( w(j,:,:))).*fbp(:,:,j)	,'symmetric');
%	omegaz(j,:,:)=ifft2(fft2(squeeze(omegaz(j,:,:))).*fbp(:,:,j)	,'symmetric');
%	omegay(j,:,:)=ifft2(fft2(squeeze(omegay(j,:,:))).*fbp(:,:,j)	,'symmetric');

	
	v(j,:,:)=ifft2(         fft2(squeeze( v(j,:,:))).*fbp(:,:,j)    ,'symmetric');
        w(j,:,:)=ifft2(         fft2(squeeze( w(j,:,:))).*fbp(:,:,j)    ,'symmetric');
        omegaz(j,:,:)=ifft2(fft2(squeeze(omegaz(j,:,:))).*fbp(:,:,j)    ,'symmetric');
        omegay(j,:,:)=ifft2(fft2(squeeze(omegay(j,:,:))).*fbp(:,:,j)    ,'symmetric');
	
	
	
	%	viscous_x(j,:,:)=ifft2(fft2(squeeze( viscous_x(j,:,:))).*f(:,:,j),'symmetric');
		
	end
	
%	convective_x=v.*omegaz-w.*omegay;
%	total=viscous_x+convective_x;




	fl=sprintf("Transfer_fhp_2Dinertial_%03d",time)
	ml=matfile(fl,'Writable',true);
	ml.voz=v.*omegaz;
	ml.woy=w.*omegay;
%	ml.convective_x=convective_x;
%	ml.viscous_x=viscous_x;
%	ml.total=total;
	ml.jstart=jstart;
	mk.jend=jend;
end

