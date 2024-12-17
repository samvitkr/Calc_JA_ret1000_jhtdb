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
lrmslp=(zeros(Nj,1));
lrmshp=(zeros(Nj,1));
llpsq=(zeros(Nj,1));
lhpsq=(zeros(Nj,1));
for time =tstart:tend

	m=matfile('Transfer_lambda_mean_f2D.mat','Writable',true)
	lmhp=m.lmhp;
	lmlp=m.lmlp;


	flhp=sprintf("lambda_fhp_2Dinertial_%03d",time)
        mlhp=matfile(flhp);
	lhpsq=lhpsq+(rms( rms( mlhp.lambda2-lmhp ,3) ,2)).^2;

        fl=sprintf("lambda_f_2Dinertial_%03d",time)
        mllp=matfile(fl);
	llpsq=llpsq+(rms( rms( mllp.lambda2-lmlp ,3) ,2)).^2;

end
lhpsq=lhpsq./nt;
llpsq=llpsq./nt;
m.lrmslp=sqrt(llpsq);
m.lrmshp=sqrt(lhpsq);
