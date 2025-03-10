jcond=171;
load('../data/ygrid.mat')
nx=512;
nz=384;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
fn=sprintf('../data/velgrad_corr_v_j_%03d.mat',jcond);
m=matfile(fn);
rvv=m.Rvv;
rvv=fftshift(fftshift(rvv,1),2);
[mxv,id]=max(rvv,[],'all');
[im,jm,km]=ind2sub(size(rvv),id);
rvv=rvv./mxv;
rvv5=((100.*abs(rvv)>5));

for i = 1:nx
    sl=squeeze(rvv5(:,i,:));
if(max(sl,[],"all")>0)
    i1=i-1;
    break
end
end

for k = 1:nz
    sl=squeeze(rvv5(k,:,:));
if(max(sl,[],"all")>0)
    k1=k-1;
    break
end
end

for i = nx:-1:1
    sl=squeeze(rvv5(:,i,:));
if(max(sl,[],"all")>0)
    i2=i+1;
    break
end
end

for k = nz:-1:1
    sl=squeeze(rvv5(k,:,:));
if(max(sl,[],"all")>0)
    k2=k+1;
    break
end
end

fnf=sprintf('../data/velgrad_v_lse_j_%03d.mat',jcond);
	ml=matfile(fnf,'Writable',true);
    ml.i1bb=i1;
    ml.i2bb=i2;
    ml.k1bb=k1;
    ml.k2bb=k2;
    %%
    [X,Z]=meshgrid(xp,zp);
    pcolor(X,Z,rvv(:,:,jcond-110))
    shading flat
    xline(xp(i1))
    xline(xp(i2))
    yline(zp(k1))
    yline(zp(k2))

    clim([-0.1 0.1])
    colormap redblue
    xlabel('x')
    ylabel('z')
    colorbar
    title('\rho_{vv}')

    axis equal
    xlim([- 2*pi 2*pi])
    ylim([-pi pi])