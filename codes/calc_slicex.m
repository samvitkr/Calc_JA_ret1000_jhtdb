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
jcond=41;

fnp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
m=matfile(fnp);

fnn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
mu=matfile(fnn);

islice=nbx+1;
%islice=nbx-25:nbx+25;
	fx=sprintf('../data/lse_xslice_i_%03d_j_%03d.mat',islice,jcond)
%	fx=sprintf('../data/lse_xsliceset_j_%03d.mat',jcond)
	mx=matfile(fx,'Writable',true)
	
	ud=  squeeze(m.u(:,islice,:))';	
	vd=  squeeze(m.v(:,islice,:))';
	wd=  squeeze(m.w(:,islice,:))';
	oxd= squeeze(m.dwdy(:,islice,:)-m.dvdz(:,islice,:))';
	ozd= squeeze(m.dvdx(:,islice,:)-m.dudy(:,islice,:))';
	oyd= squeeze(m.dudz(:,islice,:)-m.dwdx(:,islice,:))';
	ld=  squeeze(m.lambda2(:,islice,:))';
	vozd=squeeze(m.v(:,islice,:).*(m.dvdx(:,islice,:)-m.dudy(:,islice,:)))';
        woyd=squeeze(m.w(:,islice,:).*(m.dudz(:,islice,:)-m.dwdx(:,islice,:)))';
        vozdc=squeeze(m.voz(:,islice,:))';
        woydc=squeeze(m.woy(:,islice,:))';

	uu= squeeze(mu.u(:,islice,:))';	
	vu=  squeeze(mu.v(:,islice,:))';
	wu=  squeeze(mu.w(:,islice,:))';
	oxu= squeeze(mu.dwdy(:,islice,:)-mu.dvdz(:,islice,:))';
	ozu= squeeze(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:))';
	oyu= squeeze(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:))';
	lu=  squeeze(mu.lambda2(:,islice,:))';
	vozu=squeeze(mu.v(:,islice,:).*(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:)))';
        woyu=squeeze(mu.w(:,islice,:).*(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:)))';
        vozuc=squeeze(mu.voz(:,islice,:))';
        woyuc=squeeze(mu.woy(:,islice,:))';


	mx.ud=ud;
	mx.uu=uu;	
	mx.vd=vd;
	mx.vu=vu;
	mx.oxd=oxd;
	mx.oxu=oxu;
	mx.ozd=ozd;
	mx.ozu=ozu;
	mx.wd=wd;
	mx.wu=wu;
	mx.oyd=oyd;
	mx.oyu=oyu;
	mx.vozd=vozd;
	mx.vozu=vozu;
	mx.woyd=woyd;
	mx.woyu=woyu;
	mx.vozdc=vozdc;
        mx.vozuc=vozuc;
        mx.woydc=woydc;
        mx.woyuc=woyuc;
	mx.ld=ld;
	mx.lu=lu;
	mx.Z=Z;
	mx.Y=Y;

%%
