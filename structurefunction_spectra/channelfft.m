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
Nt = 11;
Lx = 8*pi;
Lz = 3*pi;

kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
xp = [0:Nx-1]*Lx/(Nx);
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz);

timeseries =[0:1:10];
ret = 5186;
yplus = 2000;
yloc = yplus/ret;
yp = -1+yloc;


p1 = xp;
p3 = xp*0;
p2=xp*0+yp;

%pointset=[];
%for k =1:Nzproc
%    p3 = xp*0+zp(k);
%    points = [p1;p2;p3];
%    pointset = [ pointset,points];
%end
npoints = Nx*Nzproc;
%
uslice=zeros(Nx,Nzproc);
vslice=zeros(Nx,Nzproc);
wslice=zeros(Nx,Nzproc);

fuuxslice=zeros(Nx,Nzproc);
fvvxslice=zeros(Nx,Nzproc);
fwwxslice=zeros(Nx,Nzproc);

fuvxslice=zeros(Nx,Nzproc);
fvwxslice=zeros(Nx,Nzproc);
fuwxslice=zeros(Nx,Nzproc);



mspec = matfile('Full_spectra_x_2000.mat','Writable',true);
mspec.fuux=zeros(Nx,1);
mspec.fvvx=zeros(Nx,1);
mspec.fwwx=zeros(Nx,1);

mspec.fuvx=zeros(Nx,1);
mspec.fvwx=zeros(Nx,1);
mspec.fuwx=zeros(Nx,1);


for t =0:(Nt-1)
    time=t
    
    for proc=1:nproc
        proc
        pointset=[];
        for k =1:Nzproc
            p3 = xp*0+zp( (proc-1)*Nzproc+k);
            points = [p1;p2;p3];
            pointset = [ pointset,points];
        end
        
        
        
        tic
        vel =  getVelocity (authkey, dataset, time, Lag8, PCHIPInt, npoints, pointset);
        toc
        
        
        vel1 =  vel(1,:);
        vel2 =  vel(2,:);
        vel3 =  vel(3,:);
        for k =1:Nzproc
            ufieldrow = vel1( (k-1)*Nx+1:k*Nx );
            vfieldrow = vel2( (k-1)*Nx+1:k*Nx );
            wfieldrow = vel3( (k-1)*Nx+1:k*Nx );
            uslice(:,k)=ufieldrow;
            vslice(:,k)=vfieldrow;
            wslice(:,k)=wfieldrow;
        end
        fu(:,:)=fft(uslice);
        fv(:,:)=fft(vslice);
        fw(:,:)=fft(wslice);
        
        fus=conj(fu);
        fvs=conj(fv);
        fws=conj(fw);
        
        fuuxslice=fuuxslice+fu.*fus./Nx;
        fvvxslice=fvvxslice+fv.*fvs./Nx;
        fwwxslice=fwwxslice+fw.*fws./Nx;
        
        fuvxslice=fuvxslice+fu.*fvs./Nx;
        fvwxslice=fvwxslice+fv.*fws./Nx;
        fuwxslice=fuwxslice+fu.*fws./Nx;
        
        
    end
    
    mspec.fuux=mspec.fuux+sum(fuuxslice,2);
    mspec.fvvx=mspec.fvvx+sum(fvvxslice,2);
    mspec.fwwx=mspec.fwwx+sum(fwwxslice,2);
    
    mspec.fuvx=mspec.fuvx+sum(fuvxslice,2);
    mspec.fvwx=mspec.fvwx+sum(fvwxslice,2);
    mspec.fuwx=mspec.fuwx+sum(fuwxslice,2);
    
end

    mspec.fuux=mspec.fuux./(Nz*Nt);
    mspec.fvvx=mspec.fvvx./(Nz*Nt);
    mspec.fwwx=mspec.fwwx./(Nz*Nt);
    
    mspec.fuvx=mspec.fuvx./(Nz*Nt);
    mspec.fvwx=mspec.fvwx./(Nz*Nt);
    mspec.fuwx=mspec.fuwx./(Nz*Nt);
