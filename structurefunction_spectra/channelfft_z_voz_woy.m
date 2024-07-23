clear all;
close all;

authkey = 'edu.jhu.skumar67-bc933816';
dataset = 'channel5200';
NoTInt   = 'None' ; % No temporal interpolation
PCHIPInt = 'PCHIP'; % Piecewise cubic Hermit interpolation in time

% ---- Spatial Interpolation Flags for getVelocity & getVelocityAndPressure ----
NoSInt = 'None'; % No spatial interpolation
Lag4   = 'Lag4'; % 4th order Lagrangian interpolation in space
Lag6   = 'Lag6'; % 6th order Lagrangian interpolation in space
Lag8   = 'Lag8'; % 8th order Lagrangian interpolation in space
Nx = 10240;
Nz = 7680;
Nzproc=96;
nproc=Nz/Nzproc;
%Nt = 11;
Nt=1;

Lx = 8*pi;
Lz = 3*pi;

kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
xp = [0:Nx-1]*Lx/(Nx);
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz);

timeseries =[0:1:10];
ret = 5185.897;
yplus = 1504.44;
yloc = yplus/ret;
yp = -1+yloc;


p1 = zp*0;
p3 = zp;
p2=zp*0+yp;

%pointset=[];
%for k =1:Nzproc
%    p3 = xp*0+zp(k);
%    points = [p1;p2;p3];
%    pointset = [ pointset,points];
%end
npoints = Nz;
%

vozF=zeros(1,Nz);
woyF=zeros(1,Nz);




mspec = matfile('spectra_voz_woyx_z_1504_5200.mat','Writable',true);
mspec.vozF=zeros(1,Nz);
mspec.woyF=zeros(1,Nz);
mspec.ns=0;
istart=10;
iend=10240;
iskip=10;
ni=(iend-istart)/iskip+1;
ns=0;
for t =0:(Nt-1)
    time=t
    
   for i=istart:iskip:iend
       i 
        pointset=[];
       % for i =1:Nx
            p1 = zp*0+xp(i);
            points = [p1;p2;p3];
            pointset = [ points];
       % end
        
        
        
        tic
        vel =  getVelocity (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
        toc
        
        
        v =  vel(2,:);
        w =  vel(3,:);
        fv=fft(v);
        fw=fft(w);
      
       	tic	
	velgrad =  getVelocityGradient (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
	toc	
        
	oy=velgrad(3,:)-velgrad(7,:);
 	oz=velgrad(4,:)-velgrad(2,:);	
      
	foy=fft(oy);
        foz=fft(oz);

	
	vozF=vozF+fv.*conj(foz)./Nz;
	woyF=woyF+fw.*conj(foy)./Nz;
	
        
        
   %end
   ns=ns+1;
    mspec.vozF=mspec.vozF+vozF;
    mspec.woyF=mspec.woyF+woyF;
    mspec.ns=ns;
end
end
    %mspec.vozF=mspec.vozF./(ni*Nt);
    %mspec.woyF=mspec.woyF./(ni*Nt);
