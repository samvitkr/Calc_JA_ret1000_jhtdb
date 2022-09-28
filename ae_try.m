Ny=512;
ut = 0.0499;
dnu=1.0006e-3;
njstart=1;
njend=256;
nj=njend-njstart+1;
%load('transfer_profile.mat');
%edges=[-3:0.1:0];
%edges=flip(-logspace(-1,0.3010,20),2);
%edges=[-2:0.05:-0.05];
%edges=[0.05:0.05:2];
%load('smoothedges.mat');
%edges=[-0.1:0.01:0.1];
%edges=sme*0.1;
%edges=[-inf,-0.5:0.1:-0.2,sme*0.1,0.2:0.1:1,inf]
%edges=[-inf,-50:10:-10,-9:1:9,10:10:50,inf];
%edges=[-inf,-20:1:20,inf];
%edges=[-inf,-0.05:0.01:0.05,inf];
%edges=[-inf,-10:1:10,inf];
edges=[-inf,-1:0.5:1,inf];
ne=length(edges);
%p=zeros(Ny/2,ne-1);
p=zeros(nj,ne-1);
pconv=p;
pvisc=p;
ptot=p;
pu=p;
pv=p;
pconv_u=p;
pconv_v=p;
pvisc_u=p;
pvisc_v=p;

mb=matfile('bsplinedata.mat');
yp = mb.yv;
yl=yp+1;
nstart=1;
nend=35;
load('JHTDB_RET1000.mat');
%load('lambda_stats.mat')
for time=nstart:nend
	time
%	fl=sprintf("lambda_filtered_15x1p25z_%03d",time);
%	fl=sprintf("lambda_%03d",time);
%  	fl=sprintf("uv_%03d",time);
	fvel=sprintf("vel_%03d.mat",time);
%	ml=matfile(fl);
	mvel=matfile(fvel);
%	l=mvel.v(njstart:njend,:,:)./ut;
	
	ur=sqrt( JHTDB_RET1000(njstart:njend,4) ).*ut;
	vr=sqrt( JHTDB_RET1000(njstart:njend,5) ).*ut;
	
	vp= mvel.v(njstart:njend,:,:)./vr;
%	up=mvel.v;
	up=(mvel.u(njstart:njend,:,:) -mean( mean( mvel.u(njstart:njend,:,:),3),2))./ur;
%	l=up(njstart:njend,:,:)./ut;
%	l=l./ur;
	%	uv=mvel.v(njstart:njend,:,:).*up(njstart:njend,:,:);
%	urvr=sqrt(JHTDB_RET1000(njstart:njend,4).*JHTDB_RET1000(njstart:njend,5));
%	uv=uv(:,:,:)./urvr;
%	l=uv./ut^2;
%	[ccm,~,binsm]=histcounts(up,[-inf 0]);
%	[ccp,~,binsp]=histcounts(up,[0 inf]);
%	l=ml.lambda2(njstart:njend,:,:);
	%size(l)
	%size(lrms)
%	l=l./lrms(njstart:njend);
%	l=l(:,:,:).*(y(njstart:njend).^2);
%	yl=yl./ut^2;
%	l=ml.Q*dnu/ut;
%    l=(ml.uv)./ut^2;
%	l=uv./ut^2;
	ft=sprintf("Transfer_%03d.mat",time)

%	ft=sprintf("Transfer_filtered_fluc_x15z1p25_%03d.mat",time)
%uprofile=mean( mean( ml.u,3),2);
    	%l=( ml.u - uprofile )./ut;
	mcv=matfile(ft);
	conv=mcv.convective_x(njstart:njend,:,:);
    	visc=mcv.viscous_x(njstart:njend,:,:);
%	tot=conv+visc;
%	l=tot./ut^2;
	%%
	
	for i =1:ne-1
		i

		l=up;
		[counts, ~, bins] = histcounts(l(:,:,:), [edges(i) edges(i+1)]);
		size(l)
		size(p)
		size(mean(mean(bins,3),2))
		pu(:,i)=pu(:,i)+mean(mean(bins,3),2);
		pconv_u(:,i)=pconv_u(:,i)+mean(mean(bins.*conv(:,:,:),3),2);%./p(:,i);
        	pvisc_u(:,i)=pvisc_u(:,i)+mean(mean(bins.*visc(:,:,:),3),2);%./p(:,i);
		
		l=vp;
                [counts, ~, bins] = histcounts(l(:,:,:), [edges(i) edges(i+1)]);
                size(l)
                size(p)
                size(mean(mean(bins,3),2))
                pv(:,i)=pv(:,i)+mean(mean(bins,3),2);
                pconv_v(:,i)=pconv_v(:,i)+mean(mean(bins.*conv(:,:,:),3),2);%./p(:,i);
                pvisc_v(:,i)=pvisc_v(:,i)+mean(mean(bins.*visc(:,:,:),3),2);%./p(:,i);


%		pconvm(:,i)=pconvm(:,i)+mean(mean( (binsm.*bins).*conv,3),2);%./p(:,i);
%               pviscm(:,i)=pviscm(:,i)+mean(mean( (binsm.*bins).*visc,3),2);%./p(:,i);
		
%		pconvp(:,i)=pconvp(:,i)+mean(mean( (binsp.*bins).*conv,3),2);%./p(:,i);
%               pviscp(:,i)=pviscp(:,i)+mean(mean( (binsp.*bins).*visc,3),2);%./p(:,i);
		

	end
end



nt=nend-nstart+1;

pu=pu./nt;
pconv_u=pconv_u/nt;
pvisc_u=pvisc_u/nt;

pv=pv./nt;
pconv_v=pconv_v/nt;
pvisc_v=pvisc_v/nt;

%pconvm=pconvm/20.0;
%pviscm=pviscm/20.0;
%pconvp=pconvp/20.0;
%pviscp=pviscp/20.0;

%mae=matfile('conditioned_ae_upv_v','Writable',true);
mae=matfile('conditioned_velrms','Writable',true);
mae.edges=edges;
mae.pu=pu;
mae.pconv_u=pconv_u;
mae.pvisc_u=pvisc_u;

mae.pv=pv;
mae.pconv_v=pconv_v;
mae.pvisc_v=pvisc_v;

%mae.pconvp=pconvp;
%mae.pviscp=pviscp;

%mae.pconvm=pconvm;
%mae.pviscm=pviscm;


