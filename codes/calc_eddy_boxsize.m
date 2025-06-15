close all
clear
load('../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
xp=ret*(lx*[0:nx-1]/nx-lx/2);
zp=ret*(lz*[0:nz-1]/nz-lz/2);
yp=ret*(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;
jcond=47;
yc=yv(jcond)+1;
ut=0.0499;
dnu=1.0006e-3;
jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_1_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;
zmid=wzz+1;
xmid=wxx+1;
xp1=xp(itarget-wxx:itarget+wxx);
zp1=zp(ktarget-wzz:ktarget+wzz);
l=-ut^2/yc^2;
%l2=-ut^2/yc^2;
val=6*l;
cutoff=0.01;
[X,Z,Y]=(meshgrid(xp1,zp1,yp));
lj1p=squeeze(sum(m1.lambda2(:,:,:)<val,[1 2]));
lj2p=squeeze(sum(m2.lambda2(:,:,:)<val,[1 2]));

lj1pm=max(lj1p);
lj2pm=max(lj2p);
for j=jcond:-1:1
	if (lj1p(j)<cutoff*lj1pm)
		break
	end
end
jbot1=j;
for j=jcond:Ny
        if (lj1p(j)<cutoff*lj1pm)
                break
        end
end
jtop1=j;

for j=jcond:-1:1
	if (lj2p(j)<cutoff*lj1pm)
		break
	end
end
jbot2=j;
for j=jcond:Ny
        if (lj2p(j)<cutoff*lj1pm)
                break
        end
end
jtop2=j;


lk1p=squeeze(sum(m1.lambda2(:,:,jbot1:jtop1)<val,[2 3]));
lk2p=squeeze(sum(m2.lambda2(:,:,jbot2:jtop2)<val,[2 3]));
lk1pm=max(lk1p);
lk2pm=max(lk2p);

li1p=squeeze(sum(m1.lambda2(:,:,jbot1:jtop1)<val,[1 3]));
li2p=squeeze(sum(m2.lambda2(:,:,jbot2:jtop2)<val,[1 3]));
li1pm=max(li1p);
li2pm=max(li2p);

for k=zmid:-1:1
        if(lk1p(k)<cutoff*lk1pm)
                break
        end
end
kstart1=k;

for k=zmid:nzz
	if(lk1p(k)<cutoff*lk1pm)
		break
	end
end
kend1=k;

for k=zmid:-1:1
        if(lk2p(k)<cutoff*lk2pm)
                break
        end
end
kstart2=k;

for k=zmid:nzz
        if(lk2p(k)<cutoff*lk2pm)
                break
        end
end
kend2=k;

for i=xmid:-1:1
        if(li1p(i)<cutoff*li1pm)
                break
        end
end
istart1=i;

for i=xmid:nxx
        if(li1p(i)<cutoff*li1pm)
                break
        end
end
iend1=i;

for i=xmid:-1:1
        if(li2p(i)<cutoff*li2pm)
                break
        end
end
istart2=i;

for i=xmid:nxx
        if(li2p(i)<cutoff*li2pm)
                break
        end
end
iend2=i;

m1.xstart=xp1(istart1);
m2.xstart=xp1(istart2);
m1.ystart=yp(jbot1);
m2.ystart=yp(jbot2);
m1.zstart=zp1(kstart1);
m2.zstart=zp1(kstart2);

m1.xend=xp1(iend1);
m2.xend=xp1(iend2);
m1.yend=yp(jtop1);
m2.yend=yp(jtop2);
m1.zend=zp1(kend1);
m2.zend=zp1(kend2);

m1.xsize=m1.xend-m1.xstart;
m1.ysize=m1.yend-m1.ystart;
m1.zsize=m1.zend-m1.zstart;

m2.xsize=m2.xend-m2.xstart;
m2.ysize=m2.yend-m2.ystart;
m2.zsize=m2.zend-m2.zstart;
