Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nu=5e-5;
tstart=1;
tend=38;

load('filter_2D_flower.mat')
nu=5e-5;
Nj=jend-jstart+1;
nt=tend-tstart+1;

vozlp=(zeros(Nj,1));
woylp=(zeros(Nj,1));
vozhp=(zeros(Nj,1));
woyhp=(zeros(Nj,1));
voz=(zeros(Ny,1));
woy=(zeros(Ny,1));

lmhp=(zeros(Nj,1));
lmlp=(zeros(Nj,1));

for time =tstart:tend

	ft=sprintf("Transfer_%03d",time)
        mt=matfile(ft);
        voz=voz+mean(mean(mt.voz,3),2);
        woy=woy+mean(mean(mt.woy,3),2);	

	flp=sprintf("Transfer_f_2Dinertial_%03d",time)
	ml=matfile(flp);
        vozlp=vozlp+mean(mean(ml.voz,3),2);
        woylp=woylp+mean(mean(ml.woy,3),2);

	fhp=sprintf("Transfer_fhp_2Dinertial_%03d",time)
        mh=matfile(fhp);
        vozhp=vozhp+mean(mean(mh.voz,3),2);
        woyhp=woyhp+mean(mean(mh.woy,3),2);
	
	flhp=sprintf("lambda_fhp_2Dinertial_%03d",time)
	mlhp=matfile(flhp);	
	lmhp  =lmhp+mean(mean( mlhp.lambda2,3),2);

	fl=sprintf("lambda_f_2Dinertial_%03d",time)
        mllp=matfile(fl);
        lmlp  =lmlp+mean(mean( mllp.lambda2,3),2);

end

fm=sprintf("Transfer_lambda_mean_f2D.mat",time);
m=matfile(fm,'Writable',true)
m.voz=voz./nt;
m.vozlp=vozlp./nt;
m.vozhp=vozhp./nt;
m.woy=woy./nt;
m.woylp=woylp./nt;
m.woyhp=woyhp./nt;
m.lmlp=lmlp./nt;
m.lmhp=lmhp./nt;



