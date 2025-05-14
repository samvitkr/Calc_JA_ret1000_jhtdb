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
jcond=130;
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




[X,Z]=meshgrid(xp,zp);

islice=wxx+1;
xp(islice)
%fx=sprintf('../data/lse_xslice_cond_i_%03d_j_%03d.mat',islice,jcond)
%	fx=sprintf('../data/lse_xsliceset_j_%03d.mat',jcond)
	%mx=matfile(fx,'Writable',true)
	
	ud=  squeeze(m.u(:,:,jcond));	
	vd=  squeeze(m.v(:,:,jcond));
	wd=  squeeze(m.w(:,:,jcond));
	oxd= squeeze(m.dwdy(:,:,jcond)-m.dvdz(:,:,jcond));
	ozd= squeeze(m.dvdx(:,:,jcond)-m.dudy(:,:,jcond));
	oyd= squeeze(m.dudz(:,:,jcond)-m.dwdx(:,:,jcond));
	ld=  squeeze(m.lambda2(:,:,jcond));
	vozd=squeeze(m.v(:,:,jcond).*(m.dvdx(:,:,jcond)-m.dudy(:,:,jcond)));
        woyd=squeeze(m.w(:,:,jcond).*(m.dudz(:,:,jcond)-m.dwdx(:,:,jcond)));
        vozdc=squeeze(m.voz(:,:,jcond));
        woydc=squeeze(m.woy(:,:,jcond));

	uu= squeeze(mu.u(:,:,jcond));	
	vu=  squeeze(mu.v(:,:,jcond));
	wu=  squeeze(mu.w(:,:,jcond));
	oxu= squeeze(mu.dwdy(:,:,jcond)-mu.dvdz(:,:,jcond));
	ozu= squeeze(mu.dvdx(:,:,jcond)-mu.dudy(:,:,jcond));
	oyu= squeeze(mu.dudz(:,:,jcond)-mu.dwdx(:,:,jcond));
	lu=  squeeze(mu.lambda2(:,:,jcond));
	vozu=squeeze(mu.v(:,:,jcond).*(mu.dvdx(:,:,jcond)-mu.dudy(:,:,jcond)));
        woyu=squeeze(mu.w(:,:,jcond).*(mu.dudz(:,:,jcond)-mu.dwdx(:,:,jcond)));
        vozuc=squeeze(mu.voz(:,:,jcond));
        woyuc=squeeze(mu.woy(:,:,jcond));


% 	mx.ud=ud;
% 	mx.uu=uu;	
% 	mx.vd=vd;
% 	mx.vu=vu;
% 	mx.oxd=oxd;
% 	mx.oxu=oxu;
% 	mx.ozd=ozd;
% 	mx.ozu=ozu;
% 	mx.wd=wd;
% 	mx.wu=wu;
% 	mx.oyd=oyd;
% 	mx.oyu=oyu;
% 	mx.vozd=vozd;
% 	mx.vozu=vozu;
% 	mx.woyd=woyd;
% 	mx.woyu=woyu;
% 	mx.vozdc=vozdc;
%         mx.vozuc=vozuc;
%         mx.woydc=woydc;
%         mx.woyuc=woyuc;
% 	mx.ld=ld;
% 	mx.lu=lu;
% 	mx.Z=Z;
% 	mx.Y=Y;

figure
subplot(1,2,1)
pcolor(X,Z,woyd-vozd)
shading interp
clim([-0.1 0.1])
axis equal
xlim([-0.1 0.1])
ylim([-0.1 0.1])

subplot(1,2,2)
pcolor(X,Z,woyu-vozu)
shading interp
clim([-0.1 0.1])
axis equal
axis tight
xlim([-0.1 0.1])
ylim([-0.1 0.1])