close all
clear
%load('ygrid.mat')

Nx=2048;
Ny=512;
Nz=1536;

Lx=8*pi;
Lz=3*pi;

%jcond=130;
%xbox=0.8;
%zbox=0.6;

%jcond=105;
%xbox=0.7;
%zbox=0.4;

jcond=71;
xbox=0.6;
zbox=0.3;

% jcond=54;
% xbox=0.5;
% zbox=0.25;

%jcond=47;
%xbox=0.5;
%zbox=0.2;

jc=Ny-jcond+1
dx=Lx/Nx;
dz=Lz/Nz;
nf=1;

tstart=1;
%tend=00000;
tend=10;
%tend=0;
tstep=1;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')

Kcoord = 1:Nz;
Icoord = 1:Nx;

itarget=Nx/2+1;
ktarget=Nz/2+1;
wini=round(xbox/dx);
wink=round(zbox/dz);
winiav=round(0.5*wini);
winkav=round(0.5*wink);
nzav=2*winkav+1;
nxav=2*winiav+1;
event_location=[];
counter=0;
vmul=10;
load('../../data/JHTDB_RET1000.mat')
uvav=abs((JHTDB_RET1000(:,3))./JHTDB_RET1000(end,2)^2);
%vrms=sqrt(0.5*(mp.vv(jcond,1)+mp.vv(jc,1)));
vthreshold=vmul*uvav(jcond)
clear JHTDB_RET1000  

un=	single(zeros(nzav,nxav,Ny/2));
vn=	single(zeros(nzav,nxav,Ny/2));
wn=	single(zeros(nzav,nxav,Ny/2));

dudxn=	single(zeros(nzav,nxav,Ny/2));
dvdxn=	single(zeros(nzav,nxav,Ny/2));
dwdxn=	single(zeros(nzav,nxav,Ny/2));

dudyn=	single(zeros(nzav,nxav,Ny/2));
dvdyn=	single(zeros(nzav,nxav,Ny/2));
dwdyn=	single(zeros(nzav,nxav,Ny/2));

dudzn=	single(zeros(nzav,nxav,Ny/2));
dvdzn=	single(zeros(nzav,nxav,Ny/2));
dwdzn=	single(zeros(nzav,nxav,Ny/2));

vozn=	single(zeros(nzav,nxav,Ny/2));
woyn=	single(zeros(nzav,nxav,Ny/2));

s=[Nz Nx];
vjav=zeros(Nz,Nx);

for time=tstart:tstep:tend
    	time

	if time == 3 || time == 4
        fprintf('Skipping corrupted time step %d...\n', time);
        continue;
    	end

%    fvel=sprintf("../data/velfields_%07d.mat",time);
	fvel=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat",time)
	m=matfile(fvel)
%
	size(m.vfield)
     	vb=permute(single(m.vfield(jcond,:,:)),[3 2 1]);
    	
	vt=permute(single(m.vfield(jc,:,:))   ,[3 2 1]);
	vt=-flip(vt,1);

	ufieldb= single(		permute(m.ufield(1:Ny/2,:,:)	,[3 2 1]));
	
	ufieldt= single(	flip(	permute(m.ufield(Ny/2+1:end,:,:),[3 2 1]),3));
	ufieldt= flip(ufieldt,1);			

 	vfieldb= single(		permute(m.vfield(1:Ny/2,:,:)	,[3 2 1]));
 	
	vfieldt= single(	flip(	permute(m.vfield(Ny/2+1:end,:,:),[3 2 1]),3));
	vfieldt= -flip(vfieldt,1);

	wfieldb=single(		permute(m.wfield(1:Ny/2,:,:)	,[3 2 1]));
	
	wfieldt=single(	flip(	permute(m.wfield(Ny/2+1:end,:,:),[3 2 1]),3));
    	wfieldt= -flip(wfieldt,1);

    	clear m

	upb = ufieldb(:,:,jcond)-mean(ufieldb(:,:,jcond),'all');
	
	upt = flip(ufieldt(:,:,jcond)-mean(ufieldt(:,:,jcond),'all'),1);
	upt = flip(upt,1);

	uvb = upb.*vb.*(vb>0);
	uvt = upt.*vt.*(vt>0);
	
	vj = uvb.*(uvb<-vthreshold);
	vjt= uvt.*(uvt<-vthreshold);

	fvelgx=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradx_%03d.mat",time);
	mgx=matfile(fvelgx)

	dudxb=single(             permute(mgx.dudx(1:Ny/2,:,:)    ,[3 2 1]));
    	dudxt=single(     flip(   permute(mgx.dudx(Ny/2+1:end,:,:),[3 2 1]),3));
    	dvdxb=single(             permute(mgx.dvdx(1:Ny/2,:,:)    ,[3 2 1]));
    	dvdxt=single(     flip(   permute(mgx.dvdx(Ny/2+1:end,:,:),[3 2 1]),3));
    	dwdxb=single(             permute(mgx.dwdx(1:Ny/2,:,:)    ,[3 2 1]));
    	dwdxt=single(     flip(   permute(mgx.dwdx(Ny/2+1:end,:,:),[3 2 1]),3));

	dudxt =  flip(dudxt,1);
	dvdxt = -flip(dvdxt,1);
	dwdxt = -flip(dwdxt,1);

	clear mgx

	fvelgy=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgrady_%03d.mat",time);
	mgy=matfile(fvelgy)
	dudyb=single(             permute(mgy.dudy(1:Ny/2,:,:)    ,[3 2 1]));
        dudyt=single(     flip(   permute(mgy.dudy(Ny/2+1:end,:,:),[3 2 1]),3));
        dvdyb=single(             permute(mgy.dvdy(1:Ny/2,:,:)    ,[3 2 1]));
        dvdyt=single(     flip(   permute(mgy.dvdy(Ny/2+1:end,:,:),[3 2 1]),3));
        dwdyb=single(             permute(mgy.dwdy(1:Ny/2,:,:)    ,[3 2 1]));
        dwdyt=single(     flip(   permute(mgy.dwdy(Ny/2+1:end,:,:),[3 2 1]),3));

	dudyt = -flip(dudyt,1);
	dvdyt =  flip(dvdyt,1);
	dwdyt =  flip(dwdyt,1);

	clear mgy

        fvelgz=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradz_%03d.mat",time);
        mgz=matfile(fvelgz)
	dudzb=single(             permute(mgz.dudz(1:Ny/2,:,:)    ,[3 2 1]));
        dudzt=single(     flip(   permute(mgz.dudz(Ny/2+1:end,:,:),[3 2 1]),3));
        dvdzb=single(             permute(mgz.dvdz(1:Ny/2,:,:)    ,[3 2 1]));
        dvdzt=single(     flip(   permute(mgz.dvdz(Ny/2+1:end,:,:),[3 2 1]),3));
        dwdzb=single(             permute(mgz.dwdz(1:Ny/2,:,:)    ,[3 2 1]));
        dwdzt=single(     flip(   permute(mgz.dwdz(Ny/2+1:end,:,:),[3 2 1]),3));
        
	dudzt=-flip(dudzt,1);
	dvdzt= flip(dvdzt,1);
	dwdzt= flip(dwdzt,1);

	clear mgz

	ft=sprintf("/vast/geyink1/skumar67/Ret_1000_data/Transfer_%03d.mat",time);
	mt=matfile(ft);
	vozb=single(             permute(mt.voz(1:Ny/2,:,:)    ,[3 2 1]));
        vozt=single(     flip(   permute(mt.voz(Ny/2+1:end,:,:),[3 2 1]),3));
        woyb=single(             permute(mt.woy(1:Ny/2,:,:)    ,[3 2 1]));
        woyt=single(     flip(   permute(mt.woy(Ny/2+1:end,:,:),[3 2 1]),3));

	vozt = flip(vozt,1);
	woyt = flip(woyt,1);

	clear mt

%
    %% towards the wall
    %bottom half
    disp('bot')

    [M,I] = min(vj(:));

    [kloc, iloc] = ind2sub(s,I);
    
    vjc=vj;
%	for ii=1:2
   while(abs(M)>abs(vthreshold))
	event_location=[event_location;kloc iloc jcond time];
        counter=counter+1
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;

        vjc=circshift(vjc,[kdelta idelta]);
        vjc(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;

        vjc=circshift(vjc,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

	Is = circshift(Icoord,idelta);
        Ks = circshift(Kcoord,kdelta);
        Ic = Is(itarget-winiav:itarget+winiav);
        Kc = Ks(ktarget-winkav:ktarget+winkav);

        ufieldb=circshift( ufieldb ,[kdelta idelta]);
        vfieldb=circshift( vfieldb ,[kdelta idelta]);
        wfieldb=circshift( wfieldb ,[kdelta idelta]);

        un=	un      +ufieldb(Kc,Ic,:);
        vn=	vn      +vfieldb(Kc,Ic,:);
        wn=	wn      +wfieldb(Kc,Ic,:);

    	dudxn=dudxn     +  dudxb(Kc,Ic,:);
        dvdxn=dvdxn     +  dvdxb(Kc,Ic,:);
        dwdxn=dwdxn     +  dwdxb(Kc,Ic,:);

    	dudyn=dudyn     +  dudyb(Kc,Ic,:);
        dvdyn=dvdyn     +  dvdyb(Kc,Ic,:);
        dwdyn=dwdyn     +  dwdyb(Kc,Ic,:);

    	dudzn=dudzn     +  dudzb(Kc,Ic,:);
        dvdzn=dvdzn     +  dvdzb(Kc,Ic,:);
        dwdzn=dwdzn     +  dwdzb(Kc,Ic,:);

        vozn=	vozn	+   vozb(Kc,Ic,:);
        woyn=	woyn	+   woyb(Kc,Ic,:);

    end
    clear ufieldb vfieldb wfieldb
    clear dudxb dvdxb dwdxb
    clear dudyb dvdyb dwdyb
    clear dudzb dvdzb dwdzb
    clear vozb woyb

%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    top half
    disp('top')
    [M,I] = min(vjt(:));
    [kloc, iloc] = ind2sub(s,I);
    vjc=vjt;

    while(abs(M)>abs(vthreshold))
	    event_location=[event_location;kloc iloc jc time];
        counter=counter+1
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;
        temp=circshift(vjc,[kdelta idelta]);
        temp(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;
        vjc=circshift(temp,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

	Is = circshift(Icoord,idelta);
        Ks = circshift(Kcoord,kdelta);
        Ic = Is(itarget-winiav:itarget+winiav);
        Kc = Ks(ktarget-winkav:ktarget+winkav);

        un=	un	+ufieldt(Kc,Ic,:);
        vn=	vn	+vfieldt(Kc,Ic,:);
        wn=	wn	+wfieldt(Kc,Ic,:);

    	dudxn=dudxn	+  dudxt(Kc,Ic,:);
        dvdxn=dvdxn     +  dvdxt(Kc,Ic,:);
        dwdxn=dwdxn     +  dwdxt(Kc,Ic,:);

        dudyn=dudyn     +  dudyt(Kc,Ic,:);
        dvdyn=dvdyn     +  dvdyt(Kc,Ic,:);
        dwdyn=dwdyn     +  dwdyt(Kc,Ic,:);

        dudzn=dudzn     +  dudzt(Kc,Ic,:);
        dvdzn=dvdzn     +  dvdzt(Kc,Ic,:);
        dwdzn=dwdzn     +  dwdzt(Kc,Ic,:);

        vozn=	vozn	+   vozt(Kc,Ic,:);
        woyn=	woyn	+   woyt(Kc,Ic,:);

    end
	clear ufieldt vfieldt wfieldt
    clear dudxt dvdxt dwdxt
    clear dudyt dvdyt dwdyt
    clear dudzt dvdzt dwdzt
    clear vozt woyt
end
%counter

fc=sprintf("../../data/conditionalq2_ref_jcond_%03d.mat",jcond);
%fc=sprintf("../data/test.mat")
mc=matfile(fc,'Writable',true);
mc.event=event_location;
mc.u=un./counter;
mc.v=vn./counter;
mc.w=wn./counter;

mc.dudx=dudxn./counter;
mc.dvdx=dvdxn./counter;
mc.dwdx=dwdxn./counter;

mc.dudy=dudyn./counter;
mc.dvdy=dvdyn./counter;
mc.dwdy=dwdyn./counter;

mc.dudz=dudzn./counter;
mc.dvdz=dvdzn./counter;
mc.dwdz=dwdzn./counter; 

mc.voz=vozn./counter;
mc.woy=woyn./counter;
