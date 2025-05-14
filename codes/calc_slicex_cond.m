
close all
clear
load('../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
xp=(lx*[0:nx-1]/nx-lx/2);
zp=(lz*[0:nz-1]/nz-lz/2);
yp=(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;
jcond=47;
jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_%03d.mat",jcond);
m=matfile(fvgp,'Writable',true);
mu=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m.lambda2);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp=xp(itarget-wxx:itarget+wxx);
zp=zp(ktarget-wzz:ktarget+wzz);




[Z,Y]=meshgrid(zp,yp);

islice=wxx+1;
xp(islice)
fx=sprintf('../data/lse_xslice_cond_i_%03d_j_%03d.mat',islice,jcond)
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
