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
jj=3;
% for counter=1:5
% jcond=jcset(jj);
% jc=jcond;
% %fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
% %fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
% % for counter = 1:1
% fvgp=sprintf("../data/conditionalp_jcond_inst_%03d_%02d.mat",jcond,counter);

jcond=105;
counter=5;
fvgp=sprintf('../data/conditionalp_dudyinst_jcond_%03d_counter_%03d.mat',jcond,counter)


% fvgn=sprintf("../data/conditionaln_jcond_inst_%03d_%02d.mat",jcond,counter);
fx=sprintf('../data/conditional_zslice_inst_j_%03d_%02d.mat',jcond,counter);
mx=matfile(fx,'Writable',true)

for nn=1:1
	switch nn
		case 1
		m=matfile(fvgp);
		case 2
		m=matfile(fvgn);
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
	switch nn
	case 1	
		mx.ud=  squeeze(m.u(kslice,:,:))';	
		mx.vd=  squeeze(m.v(kslice,:,:))';
		mx.wd=  squeeze(m.w(kslice,:,:))';

		mx.oxd= squeeze(m.dwdy(kslice,:,:)-m.dvdz(kslice,:,:))';
		mx.ozd= squeeze(m.dvdx(kslice,:,:)-m.dudy(kslice,:,:))';
		mx.oyd= squeeze(m.dudz(kslice,:,:)-m.dwdx(kslice,:,:))';
		
		mx.dudxd=squeeze(m.dudx(kslice,:,:))'; 
    		mx.dvdxd=squeeze(m.dvdx(kslice,:,:))';
    		mx.dwdxd=squeeze(m.dwdx(kslice,:,:))';

		mx.dudyd=squeeze(m.dudy(kslice,:,:))';
		mx.dvdyd=squeeze(m.dvdy(kslice,:,:))';
		mx.dwdyd=squeeze(m.dwdy(kslice,:,:))';

		mx.dudzd=squeeze(m.dudz(kslice,:,:))';
        	mx.dvdzd=squeeze(m.dvdz(kslice,:,:))';
        	mx.dwdzd=squeeze(m.dwdz(kslice,:,:))';

		mx.ld=  squeeze(m.lambda2(kslice,:,:))';
    		mx.qd=  squeeze(m.Q(kslice,:,:))';

        	mx.vozd=squeeze(m.voz(kslice,:,:))';
        	mx.woyd=squeeze(m.woy(kslice,:,:))';

	case 2	
		mx.uu=  squeeze(m.u(kslice,:,:))';
        	mx.vu=  squeeze(m.v(kslice,:,:))';
        	mx.wu=  squeeze(m.w(kslice,:,:))';

        	mx.oxu= squeeze(m.dwdy(kslice,:,:)-m.dvdz(kslice,:,:))';
        	mx.ozu= squeeze(m.dvdx(kslice,:,:)-m.dudy(kslice,:,:))';
        	mx.oyu= squeeze(m.dudz(kslice,:,:)-m.dwdx(kslice,:,:))';

        	mx.dudxu=squeeze(m.dudx(kslice,:,:))';
        	mx.dvdxu=squeeze(m.dvdx(kslice,:,:))';
        	mx.dwdxu=squeeze(m.dwdx(kslice,:,:))';

        	mx.dudyu=squeeze(m.dudy(kslice,:,:))';
        	mx.dvdyu=squeeze(m.dvdy(kslice,:,:))';
        	mx.dwdyu=squeeze(m.dwdy(kslice,:,:))';

        	mx.dudzu=squeeze(m.dudz(kslice,:,:))';
        	mx.dvdzu=squeeze(m.dvdz(kslice,:,:))';
        	mx.dwdzu=squeeze(m.dwdz(kslice,:,:))';

        	mx.lu=  squeeze(m.lambda2(kslice,:,:))';
        	mx.qu=  squeeze(m.Q(kslice,:,:))';

        	mx.vozu=squeeze(m.voz(kslice,:,:))';
        	mx.woyu=squeeze(m.woy(kslice,:,:))';
	end
	mx.X=X;
	mx.Y=Y;

end
% end
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
