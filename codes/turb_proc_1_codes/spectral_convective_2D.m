Nx=2048;
Ny=512;
Nz=1536;

ntstart=1;
ntend=38;
nt=ntend-ntstart+1;
%load('indy.mat')
%jloc=[ 38;53;75;92;106;119;172 ];
%jloc=[40;45;50;55;60;65;70;75];
indy=18;
jloc1=indy;
jloc2=512-flip(indy)+1;
jloc=[jloc1,jloc2];
Nj=length(jloc);
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
phi_v_oz=zeros(Nx,Nz,Nj);
phi_oy_w=zeros(Nx,Nz,Nj);
v_oz=phi_v_oz;
oy_w=phi_oy_w;
%oy=single(zeros(Ny,Nx,Nz));
%oz=oy;
%v=oy;
%w=oy;
for time=ntstart:ntend
	time
fvel=sprintf("vel_%03d.mat",time);
fvort=sprintf("vort_%03d.mat",time);
mvel=matfile(fvel);
mvort=matfile(fvort);
v=mvel.v;
w=mvel.w;
oy=mvort.omegay;
oz=mvort.omegaz;

	for jl =1:Nj
		j=jloc(jl);
		vslice=squeeze(v(j,:,:));
		wslice=squeeze(w(j,:,:));
		oyslice=squeeze(oy(j,:,:));
		ozslice=squeeze(oz(j,:,:));
		fv=fft2(vslice);
		fw=fft2(wslice);
		foy=fft2(oyslice);
		foz=fft2(ozslice);
		
%		size( (mean(fv.*(conj(foz)),2)./Nx) )
%		size( (mean(foy.*(conj(fw)),2)./Nx) )
		phi_v_oz(:,:,jl)=phi_v_oz(:,:,jl) + fv.*(conj(foz));
		phi_oy_w(:,:,jl)=phi_oy_w(:,:,jl)+ foy.*(conj(fw));
		v_oz(:,:,jl)=v_oz(:,:,jl) + vslice.*ozslice;
		oy_w(:,:,jl)=oy_w(:,:,jl) + oyslice.*wslice;
	

	end
end
phi_v_oz=phi_v_oz./(nt*Nx*Nz);
phi_oy_w=phi_oy_w./(nt*Nx*Nz);
v_oz=v_oz./nt;
oy_w=oy_w./nt;



m=matfile('spec_conv_2D_5.mat','Writable',true);

%m.conv(:,:,21:25)=phi_v_oz(:,:,21:25)-phi_oy_w(:,:,21:25);
m.phi_v_oz=phi_v_oz./(Nx*Nz);
m.phi_oy_w=phi_oy_w./(Nx*Nz);
m.v_oz=v_oz;
m.oy_w=oy_w;
m.jloc=jloc;
