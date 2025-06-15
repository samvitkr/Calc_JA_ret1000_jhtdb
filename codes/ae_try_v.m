Ny=512;
ut = 0.0499;
dnu=1.0006e-3;
njstart=1;
njend=256;
nj=njend-njstart+1;
edges=[-inf,-2:0.25:2,inf];
ne=length(edges);
%p=zeros(Ny/2,ne-1);
p=zeros(nj,ne-1);
pvoz=p;
pwoy=p;
pv=p;

mb=matfile('../data/bsplinedata.mat');
yp = mb.yv;
yl=yp+1;
nstart=1;
nend=10;
load('../data/JHTDB_RET1000.mat');
%load('lambda_stats.mat')
for time=nstart:nend
	time
	fvel=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat",time)
	mvel=matfile(fvel);
	
	vr=sqrt( JHTDB_RET1000(njstart:njend,5) ).*ut;
	
	vp= mvel.vfield(njstart:njend,:,:)./vr;
	ft=sprintf("/vast/geyink1/skumar67/Ret_1000_data/Transfer_%03d.mat",time)
	mcv=matfile(ft);

	voz=mcv.voz(njstart:njend,:,:);
	woy=mcv.woy(njstart:njend,:,:);
	%%
	
	for i =1:ne-1
		i

		
		l=vp;
                [counts, ~, bins] = histcounts(l(:,:,:), [edges(i) edges(i+1)]);
                size(l)
                size(p)
                size(mean(mean(bins,3),2))
                pv(:,i)=pv(:,i)+mean(mean(bins,3),2);
                pvoz(:,i)=pvoz(:,i)+mean(mean(bins.*voz(:,:,:),3),2);%./p(:,i);
                pwoy(:,i)=pwoy(:,i)+mean(mean(bins.*woy(:,:,:),3),2);%./p(:,i);


		

	end
end



nt=nend-nstart+1;


pv=pv./nt;
pvoz=pvoz/nt;
pwoy=pwoy/nt;

mae=matfile('conditioned_velrms_profile','Writable',true);
mae.edges=edges;
mae.pv=pv;
mae.pvoz=pvoz;
mae.pwoy=pwoy;

