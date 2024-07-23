Nx=2048;
Ny=512;
Nz=1536;
jstart=1;
jend=512;
Nj=jend-jstart+1;



%jloc=[ 38;53;75;92;106;119;172 ];
jloc=[jstart:jend];


load('bsplinedata.mat');
yl=yv(jloc)+1;
Delta=15;

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
phi_v_oz=zeros(Nj,Nz);
phi_oy_w=zeros(Nj,Nz);
v_oz=phi_v_oz;
oy_w=phi_oy_w;
oy=single(zeros(Ny,Nx,Nz));
oz=oy;
v=oy;
w=oy;

G= single( zeros( Nj,Nx ));

nt=38;
for time=1:nt
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

%		G(jl,:)  = exp(-( yl(jl).^2*(kx.^2)./Delta^2));
%
%                vslice=real(ifft(  fft(squeeze( v(j,:,:))).*(G(jl,:)')));
%                vslice=vslice';
%
%                wslice=real(ifft(  fft(squeeze( w(j,:,:))).*(G(jl,:)')));
%                wslice=wslice';
%
%                oyslice=real(ifft( fft(squeeze( oy(j,:,:))).*(G(jl,:)')));
%                oyslice=oyslice';
%
%                ozslice=real(ifft(  fft(squeeze( oz(j,:,:))).*(G(jl,:)')));
%                ozslice=ozslice';




		vslice=squeeze(v(j,:,:))';
		wslice=squeeze(w(j,:,:))';
		oyslice=squeeze(oy(j,:,:))';
		ozslice=squeeze(oz(j,:,:))';
		fv=fft(vslice);
		fw=fft(wslice);
		foy=fft(oyslice);
		foz=fft(ozslice);
%		size( (mean(fv.*(conj(foz)),2)./Nx) )
%		size( (mean(foy.*(conj(fw)),2)./Nx) )
		phi_v_oz(jl,:)= phi_v_oz(jl,:)+(mean(fv.*(conj(foz)),2)./Nz).';
		phi_oy_w(jl,:)= phi_oy_w(jl,:)+(mean(foy.*(conj(fw)),2)./Nz).';
		v_oz(jl,:)=v_oz(jl,:)+mean( vslice.*ozslice,2)';
		oy_w(jl,:)=oy_w(jl,:)+mean( oyslice.*wslice,2)';
	

	end
end
%phi_v_oz=phi_v_oz./20.0;
%phi_oy_w=phi_oy_w./20.0;
m=matfile('spec_conv_avz_full.mat','Writable',true);

m.conv= ( phi_v_oz-phi_oy_w )./nt;
m.phi_v_oz= phi_v_oz./nt;
m.phi_oy_w= phi_oy_w./nt;
m.v_oz=v_oz./nt;
m.oy_w=oy_w./nt;
m.jloc=jloc;
