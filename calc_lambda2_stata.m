Nx=2048;
Ny=512;
Nz=1536;

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
load('bsplinedata.mat');

nproc=8;
jstart=1;
jend=256;
Nj=jend-jstart+1;
tstart=1;
tend=35;
nt=tend-tstart+1;

yl=yv(jstart:jend)+1;
Delta=15;
Deltaz=1.25;
G=zeros(Nj,Nx);
Gz=zeros(Nj,Nz);
lrms=zeros(Nj,1);
lm=zeros(Nj,1);

mlp=matfile('lambda_stats.mat','Writable',true);
for time=tstart:tend
	fl=sprintf("lambda_%03d",time)
	ml=matfile(fl);
	
	lm=lm+mean(mean(ml.lambda2(jstart:jend,:,:),3),2);
%	l=ml.lambda2-lm;
%	rt=rms( l,3);
%	rt2=rms(rt,2);
%	l2rms_filtered(:,time)=rt2;

end
lm=lm./nt;

size(lm)
size(ml.lambda2(jstart:jend,:,:))
for time=tstart:tend
	fl=sprintf("lambda_%03d",time)
        ml=matfile(fl);
	lfluc=ml.lambda2(jstart:jend,:,:)-lm;
	lrms=lrms+mean( mean ( lfluc.^2,3 ),2);
end

lrms=lrms./nt;
lrms=sqrt(lrms);
mlp.lm=lm;
mlp.lrms=lrms;
