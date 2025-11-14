close all
clear
load('../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
jcset=[47 54 71 105 130]

for jj=1:5
jcond=jcset(jj);
jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
% for counter = 1:1
fvgp=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);

for nn=1:2
	switch nn
		case 1
		m=matfile(fvgp);
		fx=sprintf('../data/conditionalp_dudysplit_zslice_j_%03d.mat',jcond)
		case 2
		m=matfile(fvgn);
	        fx=sprintf('../data/conditionaln_dudysplit_zslice_j_%03d.mat',jcond)
	end
 % mu=matfile(fvgn,'Writable',true);
	[nzz, nxx, nyy]=size(m.lambda2);
	wzz=(nzz-1)/2;
	wxx=(nxx-1)/2;
	xp=(lx*[0:nx-1]/nx-lx/2);
	zp=(lz*[0:nz-1]/nz-lz/2);
	yp=(yv(1:Ny)'+1);
	itarget=nx/2+1;
	ktarget=nz/2+1;
	xp=xp(itarget-wxx:itarget+wxx);
	zp=zp(ktarget-wzz:ktarget+wzz);	
	
	[X,Y]=meshgrid(xp,yp);
	
	islice=wxx+1;
	kslice=wzz+1;
	xp(islice)
	mx=matfile(fx,'Writable',true)
	
	ud=  squeeze(m.up(kslice,:,:))';	
	vd=  squeeze(m.vp(kslice,:,:))';
	wd=  squeeze(m.wp(kslice,:,:))';

	oxd= squeeze(m.dwdyp(kslice,:,:)-m.dvdzp(kslice,:,:))';
	ozd= squeeze(m.dvdxp(kslice,:,:)-m.dudyp(kslice,:,:))';
	oyd= squeeze(m.dudzp(kslice,:,:)-m.dwdxp(kslice,:,:))';
	
	dudxd=squeeze(m.dudxp(kslice,:,:))'; 
    	dvdxd=squeeze(m.dvdxp(kslice,:,:))';
    	dwdxd=squeeze(m.dwdxp(kslice,:,:))';

	dudyd=squeeze(m.dudyp(kslice,:,:))';
	dvdyd=squeeze(m.dvdyp(kslice,:,:))';
	dwdyd=squeeze(m.dwdyp(kslice,:,:))';

	dudzd=squeeze(m.dudzp(kslice,:,:))';
        dvdzd=squeeze(m.dvdzp(kslice,:,:))';
        dwdzd=squeeze(m.dwdzp(kslice,:,:))';

	ld=  squeeze(m.lambda2p(kslice,:,:))';
    	qd=  squeeze(m.Qp(kslice,:,:))';

        vozd=squeeze(m.vozp(kslice,:,:))';
        woyd=squeeze(m.woyp(kslice,:,:))';

	
	uu=  squeeze(m.un(kslice,:,:))';
        vu=  squeeze(m.vn(kslice,:,:))';
        wu=  squeeze(m.wn(kslice,:,:))';

        oxu= squeeze(m.dwdyn(kslice,:,:)-m.dvdzn(kslice,:,:))';
        ozu= squeeze(m.dvdxn(kslice,:,:)-m.dudyn(kslice,:,:))';
        oyu= squeeze(m.dudzn(kslice,:,:)-m.dwdxn(kslice,:,:))';

        dudxu=squeeze(m.dudxn(kslice,:,:))';
        dvdxu=squeeze(m.dvdxn(kslice,:,:))';
        dwdxu=squeeze(m.dwdxn(kslice,:,:))';

        dudyu=squeeze(m.dudyn(kslice,:,:))';
        dvdyu=squeeze(m.dvdyn(kslice,:,:))';
        dwdyu=squeeze(m.dwdyn(kslice,:,:))';

        dudzu=squeeze(m.dudzn(kslice,:,:))';
        dvdzu=squeeze(m.dvdzn(kslice,:,:))';
        dwdzu=squeeze(m.dwdzn(kslice,:,:))';

        lu=  squeeze(m.lambda2n(kslice,:,:))';
        qu=  squeeze(m.Qn(kslice,:,:))';

        vozu=squeeze(m.vozn(kslice,:,:))';
        woyu=squeeze(m.woyn(kslice,:,:))';



	mx.ud=ud;
	mx.vd=vd;
	mx.wd=wd;

	mx.oxd=oxd;
	mx.oyd=oyd;
	mx.ozd=ozd;

	mx.dudxd=dudxd;
	mx.dvdxd=dvdxd;
	mx.dwdxd=dwdxd;

	mx.dudyd=dudyd;
        mx.dvdyd=dvdyd;
        mx.dwdyd=dwdyd;

	mx.dudzd=dudzd;
        mx.dvdzd=dvdzd;
        mx.dwdzd=dwdzd;

	mx.vozd=vozd;
	mx.woyd=woyd;
    	mx.ld=ld;
    	mx.qd=qd;


	mx.uu=uu;
        mx.vu=vu;
        mx.wu=wu;

        mx.oxu=oxu;
        mx.oyu=oyu;
        mx.ozu=ozu;

        mx.dudxu=dudxu;
        mx.dvdxu=dvdxu;
        mx.dwdxu=dwdxu;

        mx.dudyu=dudyu;
        mx.dvdyu=dvdyu;
        mx.dwdyu=dwdyu;

        mx.dudzu=dudzu;
        mx.dvdzu=dvdzu;
        mx.dwdzu=dwdzu;

        mx.vozu=vozu;
        mx.woyu=woyu;
        mx.lu=lu;
        mx.qu=qu;

	mx.X=X;
	mx.Y=Y;

end
end
    % clear mx
    % clear n
% end
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
