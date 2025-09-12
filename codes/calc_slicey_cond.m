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
jcond=71;
jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_1_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
m=matfile(fvgp,'Writable',true);
mu=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m.lambda2);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp=xp(itarget-wxx:itarget+wxx);
zp=zp(ktarget-wzz:ktarget+wzz);




% [X,Z]=meshgrid(xp,zp);
[Z,X]=meshgrid(zp,xp);
islice=wxx+1;
xp(islice)
fx=sprintf('../data/lse_yslice_cond_j_%03d.mat',jcond)
%	fx=sprintf('../data/lse_xsliceset_j_%03d.mat',jcond)
	mx=matfile(fx,'Writable',true)
	
	ud=  squeeze(m.u(:,:,jcond))';	
	vd=  squeeze(m.v(:,:,jcond))';
	wd=  squeeze(m.w(:,:,jcond))';



	oxd= squeeze(m.dwdy(:,:,jcond)-m.dvdz(:,:,jcond))';
	ozd= squeeze(m.dvdx(:,:,jcond)-m.dudy(:,:,jcond))';
	oyd= squeeze(m.dudz(:,:,jcond)-m.dwdx(:,:,jcond))';
	ld=  squeeze(m.lambda2(:,:,jcond))';
	% vozd=squeeze(m.v(:,:,jcond).*(m.dvdx(:,:,jcond)-m.dudy(:,:,jcond))');
    %     woyd=squeeze(m.w(:,:,jcond).*(m.dudz(:,:,jcond)-m.dwdx(:,:,jcond))');
        vozdc=squeeze(m.voz(:,:,jcond))';
        woydc=squeeze(m.woy(:,:,jcond))';
dvdxd=squeeze(m.dvdx(:,:,jcond))';
dudyd=squeeze(m.dudy(:,:,jcond))';
dudxd=squeeze(m.dudx(:,:,jcond))';
dvdyd=squeeze(m.dvdy(:,:,jcond))';
dwdzd=squeeze(m.dwdz(:,:,jcond))';

	uu= squeeze(mu.u(:,:,jcond))';	
	vu=  squeeze(mu.v(:,:,jcond))';
	wu=  squeeze(mu.w(:,:,jcond))';

	oxu= squeeze(mu.dwdy(:,:,jcond)-mu.dvdz(:,:,jcond))';
	ozu= squeeze(mu.dvdx(:,:,jcond)-mu.dudy(:,:,jcond))';
	oyu= squeeze(mu.dudz(:,:,jcond)-mu.dwdx(:,:,jcond))';

	lu=  squeeze(mu.lambda2(:,:,jcond))';
	% vozu=squeeze(mu.v(:,:,jcond).*(mu.dvdx(:,:,jcond)-mu.dudy(:,:,jcond))');
    %     woyu=squeeze(mu.w(:,:,jcond).*(mu.dudz(:,:,jcond)-mu.dwdx(:,:,jcond))');

        vozuc=squeeze(mu.voz(:,:,jcond))';
        woyuc=squeeze(mu.woy(:,:,jcond))';

dvdxu=squeeze(mu.dvdx(:,:,jcond))';
dudyu=squeeze(mu.dudy(:,:,jcond))';

dudxu=squeeze(mu.dudx(:,:,jcond))';
dvdyu=squeeze(mu.dvdy(:,:,jcond))';
dwdzu=squeeze(mu.dwdz(:,:,jcond))';

v2d=squeeze(m.v2(:,:,jcond))';
w2d=squeeze(m.w2(:,:,jcond))';
oy2d=squeeze(m.oy2(:,:,jcond))';
oz2d=squeeze(m.oz2(:,:,jcond))';

v2u=squeeze(mu.v2(:,:,jcond))';
w2u=squeeze(mu.w2(:,:,jcond))';
oy2u=squeeze(mu.oy2(:,:,jcond))';
oz2u=squeeze(mu.oz2(:,:,jcond))';

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
	%mx.vozd=vozd;
	%mx.vozu=vozu;
	%mx.woyd=woyd;
	%mx.woyu=woyu;
	mx.vozdc=vozdc;
        mx.vozuc=vozuc;
        mx.woydc=woydc;
        mx.woyuc=woyuc;
	mx.ld=ld;
	mx.lu=lu;
	mx.X=X;
	mx.Z=Z;
    
    mx.dvdxd=dvdxd;
    mx.dvdxu=dvdxu;

    mx.dudyd=dudyd;
    mx.dudyu=dudyu;

mx.dudxd=dudxd;
mx.dvdyd=dvdyd;
mx.dwdzd=dwdzd;

mx.dudxu=dudxu;
mx.dvdyu=dvdyu;
mx.dwdzu=dwdzu;

mx.v2d=v2d;
mx.w2d=w2d;
mx.oy2d=oy2d;
mx.oz2d=oz2d;

mx.v2u=v2u;
mx.w2u=w2u;
mx.oy2u=oy2u;
mx.oz2u=oz2u;


% 
% figure
% subplot(1,2,1)
% pcolor(X,Z,woyd-vozd)
% shading interp
% clim([-0.1 0.1])
% axis equal
% xlim([-0.1 0.1])
% ylim([-0.1 0.1])
% 
% subplot(1,2,2)
% pcolor(X,Z,woyu-vozu)
% shading interp
% clim([-0.1 0.1])
% axis equal
% axis tight
% xlim([-0.1 0.1])
% ylim([-0.1 0.1])