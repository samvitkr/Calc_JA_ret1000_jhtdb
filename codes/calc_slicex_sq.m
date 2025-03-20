close all
clear
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
xp = [0:Nx-1]*Lx/(Nx)-Lx/2;
zp=  [0:1:Nz-1]*Lz/(Nz)-Lz/2;
load('../data/bsplinedata.mat')
yp = yv;
y=yp(1:Ny/2)';
nbx=164;
nbz=180;
x=xp(Nx/2-nbx:Nx/2+nbx);
z=zp(Nz/2-nbz:Nz/2+nbz);
[Z,Y]=meshgrid(z,y);
jcond=71;

fnp=sprintf('../data/lsevp_sq_field_j_%03d.mat',jcond)
m=matfile(fnp);

%fn=sprintf('../data/lsevp_sq_field_j_%03d.mat',jcond)
%m=matfile(fn,'Writable',true);

%fnn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
%mu=matfile(fnn);

islice=nbx+1;
%islice=nbx-25:nbx+25;
	fx=sprintf('../data/lse_xslice_sq_i_%03d_j_%03d.mat',islice,jcond)
%	fx=sprintf('../data/lse_xsliceset_j_%03d.mat',jcond)
	mx=matfile(fx,'Writable',true)
	
	mx.u2lse=  squeeze(m.u2lse(:,islice,:))';	
	mx.v2lse=  squeeze(m.v2lse(:,islice,:))';
	mx.w2lse=  squeeze(m.w2lse(:,islice,:))';
	mx.ox2lse= squeeze(m.ox2lse(:,islice,:))';
	mx.oz2lse= squeeze(m.oy2lse(:,islice,:))';
	mx.oy2lse= squeeze(m.oz2lse(:,islice,:))';
	mx.Z=Z;
	mx.Y=Y;

%%
