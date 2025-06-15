close all
clear
jset=[47 54 71 105 130]
xsizep=zeros(5,1);
xsizen=zeros(5,1);
ysizep=zeros(5,1);
ysizen=zeros(5,1);
zsizep=zeros(5,1);
zsizen=zeros(5,1);
load('../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
yc=yv(jset)+1;
dnu=1.0006e-3;
yp=yc./dnu;
ut=0.0499;
ar=lx*lz;
da=ar/(nx*nz);
arp=dnu^2;

nump=zeros(5,1);
numn=zeros(5,1);
nlp=zeros(5,1);
nln=zeros(5,1);
areap=zeros(5,1);
arean=zeros(5,1);
nldprofiles=zeros(Ny,5);
nluprofiles=zeros(Ny,5);

val=6;
for j=1:5
    jcond=jset(j);
    l=-val*ut^2/yc(j)^2;

    fvgp=sprintf("../data/conditionalp_jcond_1_%03d.mat",jcond);
    fvgn=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
    m1=matfile(fvgp)
    m2=matfile(fvgn)
    [nzz, nxx, nyy]=size(m1.voz);
    nzc=0.5*(nzz+1)
    nxc=0.5*(nxx+1)

%     nld=da*squeeze(sum( (m1.voz-m1.woy).*(m1.lambda2./l>1),[1 2]));
    vold=da*squeeze(sum( (m1.lambda2./l>1),[1 2] ));
%     nlu=da*squeeze(sum( (m2.voz-m2.woy).*(m2.lambda2./l>1),[1 2]));
    volu=da*squeeze(sum( (m2.lambda2./l>1),[1 2] ));

    nld=m1.voz(nzc,nxc,jcond)-m1.woy(nzc,nxc,jcond);
    nlu=m2.voz(nzc,nxc,jcond)-m2.woy(nzc,nxc,jcond);

    areap(j)=vold(jcond);
    arean(j)=volu(jcond);

    nlp(j)=nld;
    nln(j)=nlu;
%     nlp(j)=nld(jcond);
%     nln(j)=nlu(jcond);

    xsizep(j)=m1.xsize;
    xsizen(j)=m2.xsize;

    ysizep(j)=m1.ysize;
    ysizen(j)=m2.ysize;

    zsizep(j)=m1.zsize;
    zsizen(j)=m2.zsize;

    nump(j)=length(m1.event);
    numn(j)=length(m2.event);

       nldprofiles(:,j)=m1.nlav;
       nluprofiles(:,j)=m2.nlav;
    
    
end
numt=20;


m=matfile('eddy_boxsize_ar_nl.mat','Writable',true)
m.xsizep=xsizep;
m.ysizep=ysizep;
m.zsizep=zsizep;

m.xsizen=xsizen;
m.ysizen=ysizen;
m.zsizen=zsizen;

m.nump=nump;
m.numn=numn;

m.areap=areap;
m.arean=arean;

m.nlp=nlp;
m.nln=nln;
m.nldprofiles=nldprofiles;
m.nluprofiles=nluprofiles;


%figure
%subplot(1,2,1)
%hold on
%plot(xsizep./yp,yp,'-r')
%plot(ysizep./yp,yp,'-b')
%plot(zsizep./yp,yp,'-k')
%hold off
%
%
%subplot(1,2,2)
%hold on
%plot(xsizen./yp,yp,'-r')
%plot(ysizen./yp,yp,'-b')
%plot(zsizen./yp,yp,'-k')
%hold off
