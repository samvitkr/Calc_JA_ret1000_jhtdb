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
Nt=11;

Lx = 8*pi;
Lz = 3*pi;

kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
xp = [0:Nx-1]*Lx/(Nx);
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
zp=  [0:1:Nz-1]*Lz/(Nz);

timeseries =[0:1:10];
ret = 5185.897;
yplus = 2500;
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




mspec = matfile('spectra_voz_woyx_z_2500_5200.mat','Writable',true);
%mspec.vozF=zeros(1,Nz);
%mspec.woyF=zeros(1,Nz);
istart=210-150;
iend=10210-150;
iskip=500;
ni=(iend-istart)/iskip+1

%for t =0:(Nt-1)
%    time=t
  pointset=[];
npoints=Nz*ni;

uslice=zeros(ni,Nz);
vslice=zeros(ni,Nz);
wslice=zeros(ni,Nz);
dvdxslice=zeros(ni,Nz);
dwdxslice=zeros(ni,Nz);
dudyslice=zeros(ni,Nz);
dudzslice=zeros(ni,Nz);
	for i=istart:iskip:iend
       i 
 %       pointset=[];
       % for i =1:Nx
            	p1 = zp*0+xp(i);
            	points = [p1;p2;p3];
	    	pointsip= [ zp*0+xp(i+1);p2;p3 ];
	    	pointsim=[ zp*0+xp(i-1);p2;p3 ];
	    	pointsjp= [ p1;p2+0.5/ret;p3 ];
		pointsjm=[ p1;p2-0.5/ret;p3 ];	    
            pointset = [pointset,points,pointsip,pointsim,pointsjp,pointsjm];
	    size(pointset)
       % end
	end
for t =0:(Nt-1)
    time=t        
        
        
        tic
        veli =  getVelocity (authkey, dataset, time, Lag8, NoTInt, 5*npoints, pointset);
       % veli=zeros(3,5*npoints);
	toc
	for i=1:ni
		ilstart=(i-1)*Nz*5;
		uslice(i,:)=veli(1,ilstart+1:ilstart+Nz);
		vslice(i,:)=veli(2,ilstart+1:ilstart+Nz);
		wslice(i,:)=veli(3,ilstart+1:ilstart+Nz);
		dvdxslice(i,:)=(veli(2,ilstart+Nz+1:ilstart+2*Nz)-veli(2,ilstart+2*Nz+1:ilstart+3*Nz))./(xp(4)-xp(2));
		dwdxslice(i,:)=(veli(3,ilstart+Nz+1:ilstart+2*Nz)-veli(3,ilstart+2*Nz+1:ilstart+3*Nz))./(xp(4)-xp(2));
                dudyslice(i,:)=(veli(1,ilstart+3*Nz+1:ilstart+4*Nz)-veli(1,ilstart+4*Nz+1:ilstart+5*Nz))./(1/ret);

	end
	fu(:,:)=fft(uslice(:,:).').';
        fv(:,:)=fft(vslice(:,:).').';
	fw(:,:)=fft(wslice(:,:).').';	
	size(fu)
	size(kz)
	dfu(:,:)=(fu(:,:)).*(1i*kz);
	size(dfu)
        dudzslice(:,:) = ifft(dfu(:,:).','symmetric').';
	oyslice=dudzslice-dwdxslice;
        ozslice=dvdxslice-dudyslice;	
	foy(:,:)=fft(oyslice(:,:).').';
        foz(:,:)=fft(ozslice(:,:).').';
	vozFslice=fv.*( conj(foz) )./Nz;
	woyFslice=fw.*( conj(foy) )./Nz;
	vozF=mean(vozFslice,1);
	woyF=mean(woyFslice,1);

	%       	vel=veli(:,1:Nz);
%       	velip=veli(:,Nz+1:2*Nz);
%	velim=veli(:,2*Nz+1:3*Nz);
% 	veljp=veli(:,3*Nz+1:4*Nz);
%	veljm=veli(:,4*Nz+1:5*Nz);	

%        ul =  vel(1,:);
%        vl =  vel(2,:);
%        wl =  vel(3,:);

	
% 	fu=fft(u);
% 	fv=fft(v);
%       fw=fft(w);
      
%        tic
%        velip =  getVelocity (authkey, dataset, time, NoSInt, NoTInt, npoints, [ zp*0+xp(i+1);p2;p3 ]);
%        toc
	
%	tic
%        velim =  getVelocity (authkey, dataset, time, NoSInt, NoTInt, npoints, [ zp*0+xp(i-1);p2;p3 ]);
%        toc

%	tic
%        veljp =  getVelocity (authkey, dataset, time, Lag8, NoTInt, npoints, [ p1;p2+0.5/ret;p3 ]);
%        toc

%       tic
%        veljm =  getVelocity (authkey, dataset, time, Lag8, NoTInt, npoints, [ p1;p2-0.5/ret;p3 ]);
%        toc

%	dvdx=(velip(2,:)-velim(2,:))./(xp(4)-xp(2));
%	dwdx=(velip(3,:)-velim(3,:))./(xp(4)-xp(2));
%	dudy=(veljp(1,:)-veljm(1,:))./(1/ret);
%	dfu=fu.*(1i*kz);
%	dudz=ifft(dfu,'symmetric');
%	oy=dudz-dwdx;
%	oz=dvdx-dudy;
	%       	tic	
%	velgrad =  getVelocityGradient (authkey, dataset, time, NoSInt, NoTInt, npoints, points);
%	toc	
%
% %                  getVelocityGradient (authkey, dataset, time,  FD4Lag4, NoTInt, npoints, points);
%
%	oy=velgrad(3,:)-velgrad(7,:);
% 	oz=velgrad(4,:)-velgrad(2,:);	
%      
%	foy=fft(oy);
%        foz=fft(oz);

	
%	vozF=vozF+fv.*conj(foz)./Nz;
%	woyF=woyF+fw.*conj(foy)./Nz;
	
        
  toc      
   %end
    
    mspec.vozF=mspec.vozF+vozF;
    mspec.woyF=mspec.woyF+woyF;

end
    %mspec.vozF=mspec.vozF./(Nt);
    %mspec.woyF=mspec.woyF./(Nt);
