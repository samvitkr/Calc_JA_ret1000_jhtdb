Ny=512;
Nt=20;
load('bsplinedata.mat')
%conv_mean=zeros(Ny,3);
%visc_mean=zeros(Ny,3);
%total_mean=zeros(Ny,3);
mc=zeros(Ny,Nt);
mv=zeros(Ny,Nt);
m=matfile('transfer_profile_25.mat','Writable',true);
%mc=m.mc;
%mv=m.mv;
for time=1:25
	ft=sprintf("Transfer_%03d.mat",time)
	mt=matfile(ft);
	mcp=mean(mean(mt.convective_x,3),2);
	mvp=mean(mean(mt.viscous_x,3),2);
	%mtot=mean(mean(mt.total_x,3),2);
   fvel=sprintf("vel_%03d.mat",time)
	mvel=matfile(fvel);
    mu =mean(mean( mvel.u,3),2);
    ubar=trapz(yv,mu)/2;
    
	mc(:,time)=mcp.*ubar;
	mv(:,time)=mvp.*ubar;
	%total_mean=mtot;
end
m.mconv=mc;
m.mvisc=mv;
m.mtot=(mc+mv);
