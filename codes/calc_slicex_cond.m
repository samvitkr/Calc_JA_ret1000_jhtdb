
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




[Z,Y]=meshgrid(zp,yp);

islice=wxx+1;
xp(islice)
fx=sprintf('../data/lse_xslice_cond_j_%03d.mat',jcond)
%fx=sprintf('../data/lse_xslice_cond_i_%03d_j_%03d.mat',islice,jcond)
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

%u2d=  squeeze(m.u2(:,islice,:))';
%v2d=  squeeze(m.v2(:,islice,:))';
%w2d=  squeeze(m.w2(:,islice,:))';
%ox2d=  squeeze(m.ox2(:,islice,:))';
%oy2d=  squeeze(m.oy2(:,islice,:))';
%oz2d=  squeeze(m.oz2(:,islice,:))';

	uu= squeeze(mu.u(:,islice,:))';	
	vu=  squeeze(mu.v(:,islice,:))';
	wu=  squeeze(mu.w(:,islice,:))';
	oxu= squeeze(mu.dwdy(:,islice,:)-mu.dvdz(:,islice,:))';
	ozu= squeeze(mu.dvdx(:,islice,:)-mu.dudy(:,islice,:))';
	oyu= squeeze(mu.dudz(:,islice,:)-mu.dwdx(:,islice,:))';

%u2u=  squeeze(mu.u2(:,islice,:))';
%v2u=  squeeze(mu.v2(:,islice,:))';
%w2u=  squeeze(mu.w2(:,islice,:))';
%ox2u=  squeeze(mu.ox2(:,islice,:))';
%oy2u=  squeeze(mu.oy2(:,islice,:))';
%oz2u=  squeeze(mu.oz2(:,islice,:))';

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

    %mx.u2d=u2d;
    %mx.v2d=v2d;
    %mx.w2d=w2d;
    %mx.ox2d=ox2d;
    %mx.oy2d=oy2d;
    %mx.oz2d=oz2d;

    %mx.u2u=u2u;
    %mx.v2u=v2u;
    %mx.w2u=w2u;
    %mx.ox2u=ox2u;
    %mx.oy2u=oy2u;
    %mx.oz2u=oz2u;



%%
