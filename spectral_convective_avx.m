Nx=2048;
Ny=512;
Nz=1536;
jstart=1;
jend=256;
Nj=jend-jstart+1;
%jloc=[ 38;53;75;92;106;119;172 ];
jloc=[jstart:jend];
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
phi_v_oz=zeros(Nj,Nx);
phi_oy_w=zeros(Nj,Nx);
v_oz=phi_v_oz;
oy_w=phi_oy_w;
oy=single(zeros(Ny,Nx,Nz));
oz=oy;
v=oy;
w=oy;
nt=35;
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
		phi_v_oz(jl,:)= phi_v_oz(jl,:)+(mean(fv.*(conj(foz)),2)./Nx).';
		phi_oy_w(jl,:)= phi_oy_w(jl,:)+(mean(foy.*(conj(fw)),2)./Nx).';
		v_oz(jl,:)=v_oz(jl,:)+mean( vslice.*ozslice,2)';
		oy_w(jl,:)=oy_w(jl,:)+mean( oyslice.*wslice,2)';
	

	end
end
%phi_v_oz=phi_v_oz./20.0;
%phi_oy_w=phi_oy_w./20.0;
m=matfile('spec_conv_avx.mat','Writable',true);

m.conv= ( phi_v_oz-phi_oy_w )./nt;
m.phi_v_oz= phi_v_oz./nt;
m.phi_oy_w= phi_oy_w./nt;
m.v_oz=v_oz./nt;
m.oy_w=oy_w./nt;
m.jloc=jloc;
