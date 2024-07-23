Nx=2048;
Ny=512;
Nz=1536;
Nj=7;
jloc=[ 38;53;75;92;106;119;172 ];

load('bsplinedata.mat');
yl=yv(jloc)+1;
Delta=10;
G=zeros(Nj,Nx);

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nz-1]*Lz/(Nz);
phi_v_oz=zeros(Nj,Nz,28);
phi_oy_w=zeros(Nj,Nz,28);
v_oz=phi_v_oz;
oy_w=phi_oy_w;
oy=single(zeros(Ny,Nx,Nz));
oz=oy;
v=oy;
w=oy;
for time=1:28
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

		G(jl,:)  = exp(-( yl(jl).^2*(kx.^2)./Delta^2));
		j=jloc(jl);

		vslice=real(ifft(  fft(squeeze( v(j,:,:))).*(G(jl,:)')));
		vslice=vslice';

		wslice=real(ifft(  fft(squeeze( w(j,:,:))).*(G(jl,:)')));
                wslice=wslice';

		oyslice=real(ifft( fft(squeeze( oy(j,:,:))).*(G(jl,:)')));
                oyslice=oyslice';

		ozslice=real(ifft(  fft(squeeze( oz(j,:,:))).*(G(jl,:)')));
                ozslice=ozslice';
		%vslice=squeeze(v(j,:,:))';
		%wslice=squeeze(w(j,:,:))';
		%oyslice=squeeze(oy(j,:,:))';
		%ozslice=squeeze(oz(j,:,:))';
		fv=fft(vslice);
		fw=fft(wslice);
		foy=fft(oyslice);
		foz=fft(ozslice);
%		size( (mean(fv.*(conj(foz)),2)./Nx) )
%		size( (mean(foy.*(conj(fw)),2)./Nx) )
		phi_v_oz(jl,:,time)=(mean(fv.*(conj(foz)),2)./Nz).';
		phi_oy_w(jl,:,time)=(mean(foy.*(conj(fw)),2)./Nz).';
		v_oz(jl,:,time)=mean( vslice.*ozslice,2)';
		oy_w(jl,:,time)=mean( oyslice.*wslice,2)';
	

	end
end
%phi_v_oz=phi_v_oz./20.0;
%phi_oy_w=phi_oy_w./20.0;
m=matfile('spec_conv_xfil_z.mat','Writable',true);

m.conv=phi_v_oz-phi_oy_w;
m.phi_v_oz=phi_v_oz;
m.phi_oy_w=phi_oy_w;
m.v_oz=v_oz;
m.oy_w=oy_w;
