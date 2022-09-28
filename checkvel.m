clear 
ucheck=zeros(1,1536);
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;

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
zp=  [0:1:Nz-1]*Lz/(Nz);
xp= [0:Nx-1]*Lx/(Nx);
p1=zp*0+xp(20);
p2=zp*0+yp(20);
p3=zp;
pointset = [p1;p2;p3];
time=22;
npoints=Nz;
vel =  getVelocity (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
%%
m=matfile('vel_022.mat');
mu=squeeze(m.u(20,20,:));
hold on
plot(vel(1,:),'-')
plot(mu,'.r')
hold off

