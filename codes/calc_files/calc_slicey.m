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
[Z,X]=meshgrid(z,x);
jcond=130;

fnp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
m=matfile(fnp);

fnn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
mu=matfile(fnn);

jslice=jcond;

%islice=nbx-25:nbx+25;
	fx=sprintf('../data/lse_yslice_j_%03d_jc_%03d.mat',jslice,jcond)
%	fx=sprintf('../data/lse_xsliceset_j_%03d.mat',jcond)
	mx=matfile(fx,'Writable',true)
	
	ud=  squeeze(m.u(:,:,jslice))';	
	vd=  squeeze(m.v(:,:,jslice))';
	wd=  squeeze(m.w(:,:,jslice))';
	oxd= squeeze(m.dwdy(:,:,jslice)-m.dvdz(:,:,jslice))';
	ozd= squeeze(m.dvdx(:,:,jslice)-m.dudy(:,:,jslice))';
	oyd= squeeze(m.dudz(:,:,jslice)-m.dwdx(:,:,jslice))';
	ld=  squeeze(m.lambda2(:,:,jslice))';
	vozd=squeeze(m.v(:,:,jslice).*(m.dvdx(:,:,jslice)-m.dudy(:,:,jslice)))';
        woyd=squeeze(m.w(:,:,jslice).*(m.dudz(:,:,jslice)-m.dwdx(:,:,jslice)))';
        vozdc=squeeze(m.voz(:,:,jslice))';
        woydc=squeeze(m.woy(:,:,jslice))';

	uu= squeeze(mu.u(:,:,jslice))';	
	vu=  squeeze(mu.v(:,:,jslice))';
	wu=  squeeze(mu.w(:,:,jslice))';
	oxu= squeeze(mu.dwdy(:,:,jslice)-mu.dvdz(:,:,jslice))';
	ozu= squeeze(mu.dvdx(:,:,jslice)-mu.dudy(:,:,jslice))';
	oyu= squeeze(mu.dudz(:,:,jslice)-mu.dwdx(:,:,jslice))';
	lu=  squeeze(mu.lambda2(:,:,jslice))';
	vozu=squeeze(mu.v(:,:,jslice).*(mu.dvdx(:,:,jslice)-mu.dudy(:,:,jslice)))';
        woyu=squeeze(mu.w(:,:,jslice).*(mu.dudz(:,:,jslice)-mu.dwdx(:,:,jslice)))';
        vozuc=squeeze(mu.voz(:,:,jslice))';
        woyuc=squeeze(mu.woy(:,:,jslice))';


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
	mx.X=X;

%%
