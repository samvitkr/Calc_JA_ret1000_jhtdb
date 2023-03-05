Nx=2048;
Ny=512;
Nz=1536;
Nj=7;
jloc=[ 38;53;75;92;106;119;172 ];
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
phi_v_oz=zeros(Nj,Nx,25);
phi_oy_w=zeros(Nj,Nx,25);
v_oz=phi_v_oz;
oy_w=phi_oy_w;
oy=single(zeros(Ny,Nx,Nz));
oz=oy;
v=oy;
w=oy;
for time=21:25
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
		fv=fft(vslice);
		fw=fft(wslice);
		foy=fft(oyslice);
		foz=fft(ozslice);
%		size( (mean(fv.*(conj(foz)),2)./Nx) )
%		size( (mean(foy.*(conj(fw)),2)./Nx) )
		phi_v_oz(jl,:,time)=(mean(fv.*(conj(foz)),2)./Nx).';
		phi_oy_w(jl,:,time)=(mean(foy.*(conj(fw)),2)./Nx).';
		v_oz(jl,:,time)=mean( vslice.*ozslice,2)';
		oy_w(jl,:,time)=mean( oyslice.*wslice,2)';
	

	end
end
%phi_v_oz=phi_v_oz./20.0;
%phi_oy_w=phi_oy_w./20.0;
m=matfile('spec_conv_field.mat','Writable',true);

m.conv(:,:,21:25)=phi_v_oz(:,:,21:25)-phi_oy_w(:,:,21:25);
m.phi_v_oz(:,:,21:25)=phi_v_oz(:,:,21:25);
m.phi_oy_w(:,:,21:25)=phi_oy_w(:,:,21:25);
m.v_oz(:,:,21:25)=v_oz(:,:,21:25);
m.oy_w(:,:,21:25)=oy_w(:,:,21:25);
