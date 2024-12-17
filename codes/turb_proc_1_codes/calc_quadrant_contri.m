Ny=512;
ut = 0.0499;
dnu=1.0006e-3;
njstart=1;
njend=256;
nj=njend-njstart+1;

%edges=[-inf,0,inf];
%ne=length(edges);

p=zeros(nj,4);
vozq=p;
woyq=p;
viscq=p;
mb=matfile('bsplinedata.mat');
yp = mb.yv;
yl=yp+1;
nstart=1;
nend=38;
load('JHTDB_RET1000.mat');
%load('lambda_stats.mat')
for time=nstart:nend

	time
	fvel=sprintf("vel_%03d.mat",time);
	mvel=matfile(fvel);
	
	ur=sqrt( JHTDB_RET1000(njstart:njend,4) ).*ut;
	vr=sqrt( JHTDB_RET1000(njstart:njend,5) ).*ut;
	
	vp= mvel.v(njstart:njend,:,:)./vr;
	up=(mvel.u(njstart:njend,:,:) -mean( mean( mvel.u(njstart:njend,:,:),3),2))./ur;
	
	q1=(vp>0).*(up>0);
	q2=(vp>0).*(up<0);
	q3=(vp<0).*(up<0);
	q4=(vp<0).*(up>0);

	ft=sprintf("Transfer_%03d.mat",time);
	mcv=matfile(ft);

	vozq1=mcv.voz(njstart:njend,:,:).*q1;%mcv.convective_x(njstart:njend,:,:);
    	vozq2=mcv.voz(njstart:njend,:,:).*q2;%mcv.viscous_x(njstart:njend,:,:);
	vozq3=mcv.voz(njstart:njend,:,:).*q3;
	vozq4=mcv.voz(njstart:njend,:,:).*q4;

	woyq1=mcv.woy(njstart:njend,:,:).*q1;%mcv.convective_x(njstart:njend,:,:);
        woyq2=mcv.woy(njstart:njend,:,:).*q2;%mcv.viscous_x(njstart:njend,:,:);
        woyq3=mcv.woy(njstart:njend,:,:).*q3;
        woyq4=mcv.woy(njstart:njend,:,:).*q4;	

	vozq(:,1)=vozq(:,1)+mean( mean(  vozq1,3  ) ,2 );
	vozq(:,2)=vozq(:,2)+mean( mean(  vozq2,3  ) ,2 );
	vozq(:,3)=vozq(:,3)+mean( mean(  vozq3,3  ) ,2 );
	vozq(:,4)=vozq(:,4)+mean( mean(  vozq4,3  ) ,2 );	

	woyq(:,1)=woyq(:,1)+mean( mean(  woyq1,3  ) ,2 );
        woyq(:,2)=woyq(:,2)+mean( mean(  woyq2,3  ) ,2 );
        woyq(:,3)=woyq(:,3)+mean( mean(  woyq3,3  ) ,2 );
        woyq(:,4)=woyq(:,4)+mean( mean(  woyq4,3  ) ,2 );	

end



nt=nend-nstart+1;

mae=matfile('conditioned_parts_quad','Writable',true);
mae.vozq=vozq./nt;
mae.woyq=woyq./nt;

%mae.pconvp=pconvp;
%mae.pviscp=pviscp;

%mae.pconvm=pconvm;
%mae.pviscm=pviscm;


