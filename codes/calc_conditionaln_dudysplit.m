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

%jcond=71;
%xbox=0.6;
%zbox=0.3;

%jcond=54;
%xbox=0.5;
%zbox=0.25;

jcond=47;
xbox=0.5;
zbox=0.2;

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
eventn_location=[];
eventp_location=[];
counter=0;
countern=0;
counterp=0;
vmul=1;
load('../data/JHTDB_RET1000.mat')
vrms=sqrt(JHTDB_RET1000(:,5))./JHTDB_RET1000(end,2);
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


up=     single(zeros(nzav,nxav,Ny/2));
vp=     single(zeros(nzav,nxav,Ny/2));
wp=     single(zeros(nzav,nxav,Ny/2));

dudxp=  single(zeros(nzav,nxav,Ny/2));
dvdxp=  single(zeros(nzav,nxav,Ny/2));
dwdxp=  single(zeros(nzav,nxav,Ny/2));

dudyp=  single(zeros(nzav,nxav,Ny/2));
dvdyp=  single(zeros(nzav,nxav,Ny/2));
dwdyp=  single(zeros(nzav,nxav,Ny/2));

dudzp=  single(zeros(nzav,nxav,Ny/2));
dvdzp=  single(zeros(nzav,nxav,Ny/2));
dwdzp=  single(zeros(nzav,nxav,Ny/2));

vozp=   single(zeros(nzav,nxav,Ny/2));
woyp=   single(zeros(nzav,nxav,Ny/2));

u=     single(zeros(nzav,nxav,Ny/2));
v=     single(zeros(nzav,nxav,Ny/2));
w=     single(zeros(nzav,nxav,Ny/2));

dudx=  single(zeros(nzav,nxav,Ny/2));
dvdx=  single(zeros(nzav,nxav,Ny/2));
dwdx=  single(zeros(nzav,nxav,Ny/2));

dudy=  single(zeros(nzav,nxav,Ny/2));
dvdy=  single(zeros(nzav,nxav,Ny/2));
dwdy=  single(zeros(nzav,nxav,Ny/2));

dudz=  single(zeros(nzav,nxav,Ny/2));
dvdz=  single(zeros(nzav,nxav,Ny/2));
dwdz=  single(zeros(nzav,nxav,Ny/2));

voz=   single(zeros(nzav,nxav,Ny/2));
woy=   single(zeros(nzav,nxav,Ny/2));


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


  %% towards the wall
    %bottom half
    disp('bot')

    [M,I] = min(vj(:));

    [kloc, iloc] = ind2sub(s,I);
    
    dudybflag=(dudyb(kloc,iloc,jcond))


    vjc=vj;
%	for ii=1:2
   while(abs(M)>abs(vthreshold))
	event_location=[event_location;kloc iloc jcond time];
        counter=counter+1
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;

	dudybflag=(dudyb(kloc,iloc,jcond))

        vjc=circshift(vjc,[kdelta idelta]);
        vjc(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;


        ufieldb=circshift( ufieldb ,[kdelta idelta]);
        vfieldb=circshift( vfieldb ,[kdelta idelta]);
        wfieldb=circshift( wfieldb ,[kdelta idelta]);

        dudxb=circshift( dudxb ,[kdelta idelta]);
        dvdxb=circshift( dvdxb ,[kdelta idelta]);
        dwdxb=circshift( dwdxb ,[kdelta idelta]);

        dudyb=circshift( dudyb ,[kdelta idelta]);
        dvdyb=circshift( dvdyb ,[kdelta idelta]);
        dwdyb=circshift( dwdyb ,[kdelta idelta]);

%	dudybflagcheck=dudyb(ktarget,itarget,jcond)

        dudzb=circshift( dudzb ,[kdelta idelta]);
        dvdzb=circshift( dvdzb ,[kdelta idelta]);
        dwdzb=circshift( dwdzb ,[kdelta idelta]);

        vozb=circshift( vozb ,[kdelta idelta]);
        woyb=circshift( woyb ,[kdelta idelta]);

        u=	u      +ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        v=	v      +vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        w=	w      +wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudx=dudx     +dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdx=dvdx     +dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdx=dwdx     +dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudy=dudy     +dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdy=dvdy     +dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdy=dwdy     +dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudz=dudz     +dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdz=dvdz     +dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdz=dwdz     +dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        voz=voz	+vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        woy=woy	+woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

	%% retrograde
	if(dudybflag<0)
		eventn_location=[eventn_location;kloc iloc jcond time];
		countern=countern+1;
		un=      un      +ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	vn=      vn      +vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	wn=      wn      +wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudxn=dudxn     +dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdxn=dvdxn     +dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdxn=dwdxn     +dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudyn=dudyn     +dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
		%dudyncheck = dudyn(winkav+1,winiav+1,jcond)

        	dvdyn=dvdyn     +dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdyn=dwdyn     +dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudzn=dudzn     +dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdzn=dvdzn     +dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdzn=dwdzn     +dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	vozn=vozn +vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	woyn=woyn +woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	else
		eventp_location=[eventp_location;kloc iloc jcond time];
		counterp = counterp+1;
		up=      up      +ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	vp=      vp      +vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	wp=      wp      +wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudxp=dudxp     +dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdxp=dvdxp     +dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdxp=dwdxp     +dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudyp=dudyp     +dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdyp=dvdyp     +dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdyp=dwdyp     +dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudzp=dudzp     +dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdzp=dvdzp     +dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdzp=dwdzp     +dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	vozp=vozp 	+vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	woyp=woyp 	+woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	end

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
	vjc=circshift(vjc,[-kdelta -idelta]);
        [M,I] = min(vjc(:));
        [kloc, iloc] = ind2sub(s,I);

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
    dudytflag = (dudyt(kloc,iloc,jcond));
    vjc=vjt;

    while(abs(M)>abs(vthreshold))
	event_location=[event_location;kloc iloc jc time];
        counter=counter+1
        kdelta=ktarget-kloc;
        idelta=itarget-iloc;
        temp=circshift(vjc,[kdelta idelta]);
        temp(ktarget-wink:ktarget+wink,itarget-wini:itarget+wini)=NaN;
	dudytflag = (dudyt(kloc,iloc,jcond))
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

%	dudytflagcheck=dudyt(ktarget,itarget,jcond)

        dudzt=circshift( dudzt ,[kdelta idelta]);
        dvdzt=circshift( dvdzt ,[kdelta idelta]);
        dwdzt=circshift( dwdzt ,[kdelta idelta]);

        vozt	=circshift( vozt 	,[kdelta idelta]);
        woyt	=circshift( woyt 	,[kdelta idelta]);

        u=	u	+ufieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        v=	v	-vfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        w=	w	+wfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

    	dudx=dudx        +dudxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdx=dvdx        -dvdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdx=dwdx        +dwdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        dudy=dudy        -dudyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdy=dvdy        +dvdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdy=dwdy        -dwdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        dudz=dudz        +dudzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dvdz=dvdz        -dvdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        dwdz=dwdz        +dwdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        voz=voz	+vozt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        woy=woy	+woyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);


	%% retrograde
	if (dudytflag>0)
		eventn_location=[eventn_location;kloc iloc jc time];
		countern=countern+1;
		un=      un       +ufieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	vn=      vn       -vfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	wn=      wn       +wfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudxn=dudxn        +dudxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdxn=dvdxn        -dvdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdxn=dwdxn        +dwdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudyn=dudyn        -dudyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdyn=dvdyn        +dvdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdyn=dwdyn        -dwdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

		%dudyncheck = dudyn(winkav+1,winiav+1,jcond)

        	dudzn=dudzn        +dudzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdzn=dvdzn        -dvdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdzn=dwdzn        +dwdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	vozn=vozn +vozt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	woyn=woyn +woyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	else
		eventp_location=[eventp_location;kloc iloc jc time];
		counterp=counterp+1;
		up=      up       +ufieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	vp=      vp       -vfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	wp=      wp       +wfieldt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudxp=dudxp        +dudxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdxp=dvdxp        -dvdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdxp=dwdxp        +dwdxt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudyp=dudyp        -dudyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdyp=dvdyp        +dvdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdyp=dwdyp        -dwdyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	dudzp=dudzp        +dudzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dvdzp=dvdzp        -dvdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	dwdzp=dwdzp        +dwdzt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);

        	vozp=vozp +vozt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
        	woyp=woyp +woyt(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
	end


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
	vjc=circshift(temp,[-kdelta -idelta]);

        [M,I] = max(vjc(:));
        [kloc, iloc] = ind2sub(s,I);
    	end
	clear ufieldt vfieldt wfieldt
    clear dudxt dvdxt dwdxt
    clear dudyt dvdyt dwdyt
    clear dudzt dvdzt dwdzt
    clear vozt woyt
end
%counter

fc=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);
%fc=sprintf("../data/test.mat")
mc=matfile(fc,'Writable',true);
mc.event=event_location;
mc.u=u./counter;
mc.v=v./counter;
mc.w=w./counter;

mc.dudx=dudx./counter;
mc.dvdx=dvdx./counter;
mc.dwdx=dwdx./counter;

mc.dudy=dudy./counter;
mc.dvdy=dvdy./counter;
mc.dwdy=dwdy./counter;

mc.dudz=dudz./counter;
mc.dvdz=dvdz./counter;
mc.dwdz=dwdz./counter; 

mc.voz=voz./counter;
mc.woy=woy./counter;



mc.eventn=eventn_location;
mc.un=un./countern;
mc.vn=vn./countern;
mc.wn=wn./countern;

mc.dudxn=dudxn./countern;
mc.dvdxn=dvdxn./countern;
mc.dwdxn=dwdxn./countern;

mc.dudyn=dudyn./countern;
mc.dvdyn=dvdyn./countern;
mc.dwdyn=dwdyn./countern;

mc.dudzn=dudzn./countern;
mc.dvdzn=dvdzn./countern;
mc.dwdzn=dwdzn./countern;

mc.vozn=vozn./countern;
mc.woyn=woyn./countern;



mc.eventp=eventp_location;
mc.up=up./counterp;
mc.vp=vp./counterp;
mc.wp=wp./counterp;

mc.dudxp=dudxp./counterp;
mc.dvdxp=dvdxp./counterp;
mc.dwdxp=dwdxp./counterp;

mc.dudyp=dudyp./counterp;
mc.dvdyp=dvdyp./counterp;
mc.dwdyp=dwdyp./counterp;

mc.dudzp=dudzp./counterp;
mc.dvdzp=dvdzp./counterp;
mc.dwdzp=dwdzp./counterp;

mc.vozp=vozp./counterp;
mc.woyp=woyp./counterp;
