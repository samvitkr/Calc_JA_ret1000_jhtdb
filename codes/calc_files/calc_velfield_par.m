%%
clear
numWorkers = 4; % Use 48 processors on the current node

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
load('../data/bsplinedata.mat')
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
nproc=Nz/numWorkers;
nzproc=numWorkers;
%nt=6;
tstart=11;
tend=11;
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


ufield=single(zeros(Ny,Nx,Nz/nproc));
vfield=ufield;
wfield=ufield;
ufieldslice = zeros(Ny,Nx);
vfieldslice = zeros(Ny,Nx);
wfieldslice = zeros(Ny,Nx);
for time=tstart:tend
	time
	procstart=1;
	fn=sprintf("../data/velfieldpar_%02d.mat",time);
        mn=matfile(fn,'Writable',true);
        mn.ufield=single(zeros(Ny,Nx,Nz));
        mn.vfield=single(zeros(Ny,Nx,Nz));
        mn.wfield=single(zeros(Ny,Nx,Nz));
	mn.procstart=procstart;

	for proc=procstart:nproc
		kstart=(proc-1)*nzproc+1;
		kend=proc*nzproc;
		proc
		parfor k =1:nzproc
%		for k =1:nzproc
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
		%	tic
			success=false;
			while ~success
                        	try
                        	vel =  getVelocity (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
                        	success = true;
                        	fprintf('successful for k %d \n',k);
                        	catch ME
                        	fprintf('Retrying for k %d \n',k);
                        	pause(5)
                        	end
                	end
		%	toc
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
	%		pause(300)
		end
	%	fn=sprintf("../data/velfieldpar_%02d_%03d.mat",proc,time);
	% 	mn=matfile(fn,'Writable',true);
	 	mn.ufield(:,:,kstart:kend)=single(ufield);
	 	mn.vfield(:,:,kstart:kend)=single(vfield);
	 	mn.wfield(:,:,kstart:kend)=single(wfield);
		procstart=mod(proc,nproc)+1
		mn.procstart=procstart;
	end
	tstart=tstart+1;
	mn.tstart=tstart;
end

delete(gcp('nocreate'));

