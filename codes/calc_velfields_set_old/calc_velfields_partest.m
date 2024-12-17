%%
clear
numWorkers = 32; % Use 24 processors on the current node

% Start the parallel pool with the specified number of workers
parpool('local', numWorkers);

%load('velfield.mat')
nu=5e-5;
%n=size(ufield);
%Ny=n(1);
%Nx=n(2);
%Nz=n(3);

Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
%kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
%kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz);
load('bsplinedata.mat')
%C0= colmat0 ;
%C1= colmat1 ;
%C2= colmat2 ;
yp = yv;
clear colmat0 colmat1 colmat2  yv kk knots
%%
authkey = 'edu.jhu.skumar67-bc933816';
dataset = 'channel';
NoTInt   = 'None' ; % No temporal interpolation
PCHIPInt = 'PCHIP'; % Piecewise cubic Hermit interpolation in time

% ---- Spatial Interpolation Flags for getVelocity & getVelocityAndPressure ----
NoSInt = 'None'; % No spatial interpolation
Lag4   = 'Lag4'; % 4th order Lagrangian interpolation in space
Lag6   = 'Lag6'; % 6th order Lagrangian interpolation in space
Lag8   = 'Lag8'; % 8th order Lagrangian interpolation in space

%nproc=6;
nproc=8;
nzproc=Nz/nproc;
nt=6;
tstart=21;
tend=21;
proc=1

ufield=single(zeros(Ny,Nx,Nz/nproc));
vfield=ufield;
wfield=ufield;

ufieldslice = zeros(Ny,Nx);
vfieldslice = zeros(Ny,Nx);
wfieldslice = zeros(Ny,Nx);
%%
p1 = xp;
p3 = xp*0;
%
pointset=[];
for j =1:Ny
    p2 = xp*0+yp(j);
    points = [p1;p2;p3];
    pointset = [ pointset,points];
end
npoints = Nx*Ny;
%tic



for time=tstart:tend

time

	parfor k =1:nzproc
	ufieldslice = zeros(Ny,Nx);
vfieldslice = zeros(Ny,Nx);
wfieldslice = zeros(Ny,Nx);
	    k
	    pointset=[];
	for j =1:Ny
    	p2 = xp*0+yp(j);
    	points = [p1;p2;p3];
    	pointset = [ pointset,points];
	end
	npoints = Nx*Ny;
	pointset(3,:)=zp((proc-1)*nzproc+k);
	tic
	vel =  getVelocity (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
	toc
	%toc
	    vel1 =  vel(1,:);
	    vel2 =  vel(2,:);
	    vel3 =  vel(3,:);
	    for j =1:Ny
	        ufieldrow = vel1( (j-1)*Nx+1:j*Nx );
	        vfieldrow = vel2( (j-1)*Nx+1:j*Nx );
	        wfieldrow = vel3( (j-1)*Nx+1:j*Nx );
	        ufieldslice(j,:)=ufieldrow;
	        vfieldslice(j,:)=vfieldrow;
	        wfieldslice(j,:)=wfieldrow;
	    end
		ufield(:,:,k)=ufieldslice(:,:);
		vfield(:,:,k)=vfieldslice(:,:);
		wfield(:,:,k)=wfieldslice(:,:);
	end
 fn=sprintf("velfield_%02d_%03d.mat",proc,time);
 mn=matfile(fn,'Writable',true);
 mn.ufield=single(ufield);
 mn.vfield=single(vfield);
 mn.wfield=single(wfield);
end
delete(gcp('nocreate'));


%    kefieldslice=0.5*(ufieldslice.^2 + vfieldslice.^2 +wfieldslice.^2);
%    
%    %%
%        fke(:,:)=fft(kefieldslice(:,:).').';
%        fu(:,:)=fft(ufieldslice(:,:).').';
%        fv(:,:)=fft(vfieldslice(:,:).').';
%        fw(:,:)=fft(wfieldslice(:,:).').';
%        dfke(:,:)=(fke(:,:)).*(1i*kx);
%        dfu(:,:)=(fu(:,:)).*(1i*kx);
%        dfv(:,:)=(fv(:,:)).*(1i*kx);
%        dfw(:,:)=(fw(:,:)).*(1i*kx);
%        d2fu(:,:)=(fu(:,:)).*(-kx.^2);
%        dke(:,:) = ifft(dfke(:,:).').';
%        dudxslice(:,:) = ifft(dfu(:,:).').';
%        dvdxslice(:,:) = ifft(dfv(:,:).').';
%        dwdxslice(:,:) = ifft(dfw(:,:).').';
%        
%        d2udx2(:,:) = ifft(d2fu(:,:).').';
%        coeffu=C0\ufieldslice(:,:);
%        coeffv=C0\vfieldslice(:,:);
%        coeffw=C0\wfieldslice(:,:);
%        dudyslice(:,:)=C1*coeffu;
%        dvdyslice(:,:)=C1*coeffv;
%        dwdyslice(:,:)=C1*coeffw;
%        d2udy2(:,:)=C2*coeffu;
%        convective_slice(:,:)=0.5*dke-(ufieldslice.*dudxslice+vfieldslice.*dudyslice);
%        viscous_slice(:,:)= nu*( d2udx2+d2udy2 );    
%%%
%toc
%mvel.ufield(:,:,k)=single(ufieldslice);
%mvel.vfield(:,:,k)=single(vfieldslice);
%mvel.wfield(:,:,k)=single(wfieldslice);
%mt.convective(:,:,k)=mt.convective(:,:,k)+single(convective_slice);
%mt.viscous(:,:,k)=mt.viscous(:,:,k)+single(viscous_slice);
%mvelg.dudx(:,:,k)=single(dudxslice);
%mvelg.dvdx(:,:,k)=single(dvdxslice);
%mvelg.dwdx(:,:,k)=single(dwdxslice);
%mvelg.dudy(:,:,k)=single(dudyslice);
%mvelg.dvdy(:,:,k)=single(dvdyslice);
%mvelg.dwdy(:,:,k)=single(dwdyslice);
%toc
%

 
