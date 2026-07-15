clear
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
tend=25;
proc=1

timeset=[ 12.9675,13.0000, 13.0325 ];

ufield=single(zeros(Ny,Nz,3));
%vfield=ufield;
%wfield=ufield;

ufieldslice = zeros(Ny,Nz);
vfieldslice = zeros(Ny,Nz);
wfieldslice = zeros(Ny,Nz);
pin=zeros(Ny,Nz);
pout=zeros(Ny,Nz);
%%
p1 = zp*0+xp(Nx/2);
p3 = zp*0;
%
pointset=[];
for j =1:Ny
    p2 = zp*0+yp(j);
    points = [p1;p2;p3];
    pointset = [ pointset,points];
end
npoints = Nz*Ny;

for t=1:3
time=timeset(t)
     %   for k =1:nzproc
     %       k
        tic
        vel =  getVelocity (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
        toc
        %toc
            vel1 =  vel(1,:);
            vel2 =  vel(2,:);
            vel3 =  vel(3,:);
            for j =1:Ny
                ufieldrow = vel1( (j-1)*Nz+1:j*Nz );
                vfieldrow = vel2( (j-1)*Nz+1:j*Nz );
                wfieldrow = vel3( (j-1)*Nz+1:j*Nz );
                ufieldslice(j,:)=ufieldrow;
                vfieldslice(j,:)=vfieldrow;
                wfieldslice(j,:)=wfieldrow;
            end
                ufield(:,:,t)=ufieldslice(:,:);
                %vfield(:,:,k)=vfieldslice(:,:);
                %wfield(:,:,k)=wfieldslice(:,:);
     % end
end

 fn=sprintf("uinflow_013.mat");
 mn=matfile(fn,'Writable',true);
 mn.ufield_mid=single(ufield);
%
%
%time=timeset(2);
%tic
%pressure =  getPressure (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
%toc
%
%for j =1:Ny
%               prow  = pressure( (j-1)*Nz+1:j*Nz );
%	       pin(j,:)=prow;
%end
%p1=zp*0+Lx;
%pointset=[];
%for j =1:Ny
%    p2 = zp*0+yp(j);
%    points = [p1;p2;p3];
%    pointset = [ pointset,points];
%end
%tic
%pressure =  getPressure (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
%toc
%
%for j =1:Ny
%               prow  = pressure( (j-1)*Nz+1:j*Nz );
%               pout(j,:)=prow;
%end
%mn.pin=pin;
%mn.pout=pout;
