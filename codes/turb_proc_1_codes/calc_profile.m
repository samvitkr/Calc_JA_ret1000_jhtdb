Ny=512;
Nt=38;
load('bsplinedata.mat')
%conv_mean=zeros(Ny,3);
%visc_mean=zeros(Ny,3);
%total_mean=zeros(Ny,3);
mvoz=zeros(Ny,Nt);
mwoy=zeros(Ny,Nt);
mc=zeros(Ny,Nt);
mv=zeros(Ny,Nt);

m=matfile('transfer_profile_38.mat','Writable',true);
%mc=m.mc;
%mv=m.mv;
for time=1:38
	ft=sprintf("Transfer_%03d.mat",time)
	mt=matfile(ft);
	mcp=mean(mean(mt.convective_x,3),2);
	mvp=mean(mean(mt.viscous_x,3),2);
	mvozp=mean(mean(mt.voz,3),2);
	mwoyp=mean(mean(mt.woy,3),2);

	%mtot=mean(mean(mt.total_x,3),2);
   	fvel=sprintf("vel_%03d.mat",time)
	mvel=matfile(fvel);
    	mu =mean(mean( mvel.u,3),2);
%    ubar=trapz(yv,mu)/2;
    
	mc(:,time)=mcp;
	mv(:,time)=mvp;
	mvoz(:,time)=mvozp;
	mwoy(:,time)=mwoyp;

	%total_mean=mtot;
end

m.mvoz=mvoz;
m.mwoy=mwoy;
m.mconv=mc;
m.mvisc=mv;
m.mtot=(mc+mv);
