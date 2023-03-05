clear
close all
winmax=10;
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
xp = [0:Nx-1]*Lx/(Nx);
zp=  [0:1:Nz-1]*Lz/(Nz);
dnu=1.0006e-3;
ut = 0.0499;
nu=dnu*ut;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
%kxint=floor(kx(end));
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

% kxw=zeros(1,Nx/2+2*window);
% kzw=zeros(1,Nz/2+2*window);
% kxw(window+1:window+Nx/2)=kx(1:Nx/2);
% kzw(window+1:window+Nz/2)=kx(1:Nz/2);



dkx=kx(2)-kx(1);
dkz=kz(2)-kz(1);
load('bsplinedata.mat');
% load('JHTDB_RET1000.mat');
% load('spec_conv_2D.mat');
m=matfile('spec_conv_2D_5.mat');
%m.jloc=[ 38;53;75;92;106;119;172 ];
[Kx,Kz]=meshgrid(kx,kz);
Kx=Kx';
Kz=Kz';

% [Kxw,Kzw]=meshgrid(kxw,kzw);
% Kxw=Kxw';
% Kzw=Kzw';


y=1+yv;
yp=y./dnu;
yl=y(m.jloc);
njl=length(yl);
conv=m.phi_oy_w;
%conv=(conv./(Nx*Nz));
conv=conv./(dkx*dkz);
%conv=2*conv(1:Nx/2,1:Nz/2);
%psc=conv(:,:,end);
%pscm=psc;
psc1=conv;
psc2=flip(psc1);
psc3=flip(psc1,2);
psc4=flip(psc2,2);
pscm=psc1+psc2+psc3+psc4;
pscm=real(pscm(1:Nx/2,1:Nz/2,:));

(dkx*dkz)*sum(sum(pscm),2)./ut^2

pscfull=zeros(Nx,Nz);

mw=matfile('smooth_2d_spectra_5.mat','Writable',true);

winmax=3;
window=winmax;
%mw.winmax=winmax;
wconv=zeros(winmax+Nx/2,winmax+Nz/2,njl);
for jl=1:njl
	jl
	pscfull=zeros(Nx,Nz);
	pscfull(Nx/2+1:end,Nz/2+1:end)=pscm(:,:,jl);
	pscfull=pscfull+flip(pscfull,1)+flip(pscfull,2)+ flip(flip(pscfull,1),2);
	window=winmax;
	phishift=zeros(Nx+2*window,Nz+2*window);
	phishift(window+1:window+Nx,window+1:window+Nz)=pscfull;
	%wconv(1:Nx/2,1:Nz/2,1)=pscm(:,:,jl);
	
	sphi=phishift;
	count=1;
	for i =1:window
	    
	     sphi=sphi+circshift(phishift,i);
	     sphi=sphi+circshift(phishift,-i);
	     sphi=sphi+circshift(phishift,i,2);
	     sphi=sphi+circshift(phishift,-i,2);
	 
	count=count+4;
	    for j =1:window
	    sphi=sphi+circshift(circshift(phishift,i),j,2);
	    sphi=sphi+circshift(circshift(phishift,i),-j,2);
	    sphi=sphi+circshift(circshift(phishift,-i),j,2);
	    sphi=sphi+circshift(circshift(phishift,-i),j,2);
	
	count=count+4;
	
	    end
	end
	count
	% sphi=sphi./(8*window+1);
	sphi=sphi./count;
	%[nx, nz]=size(sphi)
	% winmax-window+1
	%wconv(winmax-window+1:winmax-window+nx, winmax-window+1:winmax-window+nz ,window+1)=sphi;
	sphi=0.25*(sphi+flip(sphi,1)+flip(sphi,2)+ flip(flip(sphi,1),2));
	sphiq=sphi(window+Nx/2+1:end,window+Nz/2+1:end);
	[nx, nz]=size(sphiq)
	wconv(1:nx,1:nz,jl)=sphiq;
	%phiq=phishift(window+Nx/2+1:end,window+Nz/2+1:end);
	
	% 
	(dkx*dkz)*sum(sum(sphiq ),2)/ut^2
end
%(dkx*dkz)*sum(sum(phiq ),2)/ut^2

% 
% (dkx*dkz)*0.25*sum(sum(sphi ),2)/ut^2
%(dkx*dkz)*0.25*sum(sum(phishift ),2)/ut^2
% 
% (dkx*dkz)*sum(sum(sphi(1:window+Nx/2,1:window+Nz/2) ),2)/ut^2
% (dkx*dkz)*sum(sum(sphi(1:window+Nx/2,window+Nz/2+1:end) ),2)/ut^2
% (dkx*dkz)*sum(sum(sphi(window+Nx/2+1:end,1:window+Nz/2 ) ),2)/ut^2
 %(dkx*dkz)*sum(sum(sphi(window+Nx/2+1:end,window+Nz/2+1:end) ),2)/ut^2
% 
% 
% 
% (dkx*dkz)*sum(sum(phishift ),2)/ut^2





%%

%dw=wconv(:,:,2:end).*0;
%distw=zeros(1,winmax-1);
% for w=1:winmax
%dw(:,:,w)=wconv(:,:,w+1)-wconv(:,:,w);
%r1=rms(dw(:,:,w),1);
%distw(w)=rms(r1,2);
%
% subplot(2,3,w)
% pcolor(-wconv(:,:,w)')
% shading flat
% caxis([-1e-6 1e-6])
% end
% % % % 
 mw.woy=wconv;
 %mw.distw=distw;
 mw.jl=m.jloc;
% % % %  
% % % %  
% % % % figure
% % % % %pcolor(Kx(1:Nx/2,1:Nz/2),Kz(1:Nx/2,1:Nz/2),-sphi(window+1:window+Nx/2,window+1:window+Nz/2)./ut^2)
% % % % %shading flat
% % % % %caxis([-1e-6 1e-6])
% % % %  plot(distw,'o-')
% % % %  
% % % % % figure
% % % % % pcolor(Kx(1:Nx/2,1:Nz/2),Kz(1:Nx/2,1:Nz/2),-phishift(window+1:window+Nx/2,window+1:window+Nz/2)./ut^2)
% % % % % %pcolor(phishift./ut^2)
% % % % % shading flat
% % % % % caxis([-1e-6 1e-6])
% % % % 
% % % % % (dkx*dkz)*sum(sum(psc ),2)/ut^2
% % % % % (dkx*dkz)*sum(sum(psc(1:Nx/2,1:Nz/2) ),2)/ut^2
% % % % % (dkx*dkz)*sum(sum(psc(1:Nx/2,Nz/2+1:end)),2)/ut^2
% % % % % (dkx*dkz)*sum(sum(psc(Nx/2+1:end,Nz/2+1:end)),2)/ut^2
% % % % % (dkx*dkz)*sum(sum(psc(Nx/2+1:end,1:Nz/2)),2)/ut^2
