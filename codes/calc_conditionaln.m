close all
clear
%load('ygrid.mat')

Nx=2048;
Ny=512;
Nz=1536;

Lx=8*pi;
Lz=3*pi;

jcond=130;
xbox=0.8;
zbox=0.6;

%jcond=105;
%xbox=0.7;
%zbox=0.4;

%jcond=71
%xbox=0.6;
%zbox=0.3;

%jcond=54;
%xbox=0.5;
%zbox=0.25;

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
vmul=1;
load('../data/JHTDB_RET1000.mat')
vrms=sqrt(JHTDB_RET1000(:,5))./JHTDB_RET1000(:,2);
%vrms=sqrt(0.5*(mp.vv(jcond,1)+mp.vv(jc,1)));
vthreshold=vmul*vrms(jcond);
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
%    fvel=sprintf("../data/velfields_%07d.mat",time);
	fvel=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat",time)
	m=matfile(fvel)
%
size(m.vfield)
     vj=permute(single(m.vfield(jcond,:,:)),[3 2 1]);
    vjt=permute(single(m.vfield(jc,:,:))   ,[3 2 1]);

    ufieldb=single(		permute(m.ufield(1:Ny/2,:,:)	,[3 2 1]));
    ufieldt=single(	flip(	permute(m.ufield(Ny/2+1:end,:,:),[3 2 1]),3));
    vfieldb=single(		permute(m.vfield(1:Ny/2,:,:)	,[3 2 1]));
    vfieldt=single(	flip(	permute(m.vfield(Ny/2+1:end,:,:),[3 2 1]),3));
    wfieldb=single(		permute(m.wfield(1:Ny/2,:,:)	,[3 2 1]));
    wfieldt=single(	flip(	permute(m.wfield(Ny/2+1:end,:,:),[3 2 1]),3));
    clear m


	fvelgx=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradx_%03d.mat",time);
	mgx=matfile(fvelgx)

	dudxb=single(             permute(mgx.dudx(1:Ny/2,:,:)    ,[3 2 1]));
    	dudxt=single(     flip(   permute(mgx.dudx(Ny/2+1:end,:,:),[3 2 1]),3));
    	dvdxb=single(             permute(mgx.dvdx(1:Ny/2,:,:)    ,[3 2 1]));
    	dvdxt=single(     flip(   permute(mgx.dvdx(Ny/2+1:end,:,:),[3 2 1]),3));
    	dwdxb=single(             permute(mgx.dwdx(1:Ny/2,:,:)    ,[3 2 1]));
    	dwdxt=single(     flip(   permute(mgx.dwdx(Ny/2+1:end,:,:),[3 2 1]),3));


	clear mgx

	fvelgy=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgrady_%03d.mat",time);
	mgy=matfile(fvelgy)
	dudyb=single(             permute(mgy.dudy(1:Ny/2,:,:)    ,[3 2 1]));
        dudyt=single(     flip(   permute(mgy.dudy(Ny/2+1:end,:,:),[3 2 1]),3));
        dvdyb=single(             permute(mgy.dvdy(1:Ny/2,:,:)    ,[3 2 1]));
        dvdyt=single(     flip(   permute(mgy.dvdy(Ny/2+1:end,:,:),[3 2 1]),3));
        dwdyb=single(             permute(mgy.dwdy(1:Ny/2,:,:)    ,[3 2 1]));
        dwdyt=single(     flip(   permute(mgy.dwdy(Ny/2+1:end,:,:),[3 2 1]),3));
	clear mgy

        fvelgz=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradz_%03d.mat",time);
        mgz=matfile(fvelgz)
	dudzb=single(             permute(mgz.dudz(1:Ny/2,:,:)    ,[3 2 1]));
        dudzt=single(     flip(   permute(mgz.dudz(Ny/2+1:end,:,:),[3 2 1]),3));
        dvdzb=single(             permute(mgz.dvdz(1:Ny/2,:,:)    ,[3 2 1]));
        dvdzt=single(     flip(   permute(mgz.dvdz(Ny/2+1:end,:,:),[3 2 1]),3));
        dwdzb=single(             permute(mgz.dwdz(1:Ny/2,:,:)    ,[3 2 1]));
        dwdzt=single(     flip(   permute(mgz.dwdz(Ny/2+1:end,:,:),[3 2 1]),3));
        clear mgz

	ft=sprintf("/vast/geyink1/skumar67/Ret_1000_data/Transfer_%03d.mat",time);
	mt=matfile(ft);
	vozb=single(             permute(mt.voz(1:Ny/2,:,:)    ,[3 2 1]));
        vozt=single(     flip(   permute(mt.voz(Ny/2+1:end,:,:),[3 2 1]),3));
        woyb=single(             permute(mt.woy(1:Ny/2,:,:)    ,[3 2 1]));
        woyt=single(     flip(   permute(mt.woy(Ny/2+1:end,:,:),[3 2 1]),3));
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

        ufieldb=circshift( ufieldb ,[kdelta idelta]);
        vfieldb=circshift( vfieldb ,[kdelta idelta]);
        wfieldb=circshift( wfieldb ,[kdelta idelta]);

        dudxb=circshift( dudxb ,[kdelta idelta]);
        dvdxb=circshift( dvdxb ,[kdelta idelta]);
        dwdxb=circshift( dwdxb ,[kdelta idelta]);

        dudyb=circshift( dudyb ,[kdelta idelta]);
        dvdyb=circshift( dvdyb ,[kdelta idelta]);
        dwdyb=circshift( dwdyb ,[kdelta idelta]);

        dudzb=circshift( dudzb ,[kdelta idelta]);
        dvdzb=circshift( dvdzb ,[kdelta idelta]);
        dwdzb=circshift( dwdzb ,[kdelta idelta]);

        vozb=circshift( vozb ,[kdelta idelta]);
        woyb=circshift( woyb ,[kdelta idelta]);

        un=	un      +ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        vn=	vn      +vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        wn=	wn      +wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudxn=dudxn     +dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdxn=dvdxn     +dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdxn=dwdxn     +dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudyn=dudyn     +dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdyn=dvdyn     +dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdyn=dwdyn     +dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudzn=dudzn     +dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdzn=dvdzn     +dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdzn=dwdzn     +dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        vozn=	vozn	+vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        woyn=	woyn	+woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        ufieldb=circshift( ufieldb ,-[kdelta idelta]);
        vfieldb=circshift( vfieldb ,-[kdelta idelta]);
        wfieldb=circshift( wfieldb ,-[kdelta idelta]);
        dudxb=circshift( dudxb ,-[kdelta idelta]);
        dvdxb=circshift( dvdxb ,-[kdelta idelta]);
        dwdxb=circshift( dwdxb ,-[kdelta idelta]);
        dudyb=circshift( dudyb ,-[kdelta idelta]);
        dvdyb=circshift( dvdyb ,-[kdelta idelta]);
        dwdyb=circshift( dwdyb ,-[kdelta idelta]);
        dudzb=circshift( dudzb ,-[kdelta idelta]);
        dvdzb=circshift( dvdzb ,-[kdelta idelta]);
        dwdzb=circshift( dwdzb ,-[kdelta idelta]);
        vozb=circshift( vozb ,-[kdelta idelta]);
        woyb=circshift( woyb ,-[kdelta idelta]);

    end
    clear ufieldb vfieldb wfieldb
    clear dudxb dvdxb dwdxb
    clear dudyb dvdyb dwdyb
    clear dudzb dvdzb dwdzb
    clear vozb woyb

%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    top half
    disp('top')
    [M,I] = max(vjt(:));
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
        [M,I] = max(vjc(:));
        [kloc, iloc] = ind2sub(s,I);
%
        ufieldt	=circshift( ufieldt	,[kdelta idelta]);
        vfieldt	=circshift( vfieldt	,[kdelta idelta]);
        wfieldt	=circshift( wfieldt	,[kdelta idelta]);

    	dudxt=circshift( dudxt ,[kdelta idelta]);
        dvdxt=circshift( dvdxt ,[kdelta idelta]);
        dwdxt=circshift( dwdxt ,[kdelta idelta]);

        dudyt=circshift( dudyt ,[kdelta idelta]);
        dvdyt=circshift( dvdyt ,[kdelta idelta]);
        dwdyt=circshift( dwdyt ,[kdelta idelta]);

        dudzt=circshift( dudzt ,[kdelta idelta]);
        dvdzt=circshift( dvdzt ,[kdelta idelta]);
        dwdzt=circshift( dwdzt ,[kdelta idelta]);

        vozt	=circshift( vozt 	,[kdelta idelta]);
        woyt	=circshift( woyt 	,[kdelta idelta]);

        un=	un	+ufieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        vn=	vn	-vfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        wn=	wn	+wfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudxn=dudxn        +dudxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdxn=dvdxn        -dvdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdxn=dwdxn        +dwdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        dudyn=dudyn        -dudyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdyn=dvdyn        +dvdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdyn=dwdyn        -dwdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        dudzn=dudzn        +dudzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdzn=dvdzn        -dvdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdzn=dwdzn        +dwdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        vozn=	vozn	+vozt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        woyn=	woyn	+woyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        ufieldt	=circshift( ufieldt,-[kdelta idelta]);
        vfieldt	=circshift( vfieldt,-[kdelta idelta]);
        wfieldt	=circshift( wfieldt,-[kdelta idelta]);

    	dudxt=circshift( dudxt,-[kdelta idelta]);
        dvdxt=circshift( dvdxt,-[kdelta idelta]);
        dwdxt=circshift( dwdxt,-[kdelta idelta]);

        dudyt=circshift( dudyt,-[kdelta idelta]);
        dvdyt=circshift( dvdyt,-[kdelta idelta]);
        dwdyt=circshift( dwdyt,-[kdelta idelta]);

        dudzt=circshift( dudzt,-[kdelta idelta]);
        dvdzt=circshift( dvdzt,-[kdelta idelta]);
        dwdzt=circshift( dwdzt,-[kdelta idelta]);

        vozt	=circshift( vozt 	,-[kdelta idelta]);
        woyt	=circshift( woyt 	,-[kdelta idelta]);
    end
	clear ufieldt vfieldt wfieldt
    clear dudxt dvdxt dwdxt
    clear dudyt dvdyt dwdyt
    clear dudzt dvdzt dwdzt
    clear vozt woyt
end
%counter

fc=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
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
