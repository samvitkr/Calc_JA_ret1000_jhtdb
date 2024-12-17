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
lp_lrms=(zeros(Nj,1));
hp_lrms=(zeros(Nj,1));
lp_0=(zeros(Nj,1));
hp_0=(zeros(Nj,1));

volfrac_lp_lrms=(zeros(Nj,1));
volfrac_hp_lrms=(zeros(Nj,1));
volfrac_lp_0=(zeros(Nj,1));
volfrac_hp_0=(zeros(Nj,1));


m=matfile('Transfer_lambda_mean_f2D.mat','Writable',true)
        lrmshp=m.lrmshp;
        lrmslp=m.lrmslp;
for time =tstart:tend
        
       	flhp=sprintf("lambda_fhp_2Dinertial_%03d",time)
       	mlhp=matfile(flhp);
	lhp=mlhp.lambda2./lrmshp;

       	fllp=sprintf("lambda_f_2Dinertial_%03d",time)
       	mllp=matfile(fllp);
	llp=mllp.lambda2./lrmslp;

	flt=sprintf("Transfer_f_2Dinertial_%03d",time)
	mtlp=matfile(flt);
%	lp_0= lp_0 + mean(squeeze(mean( ( mtlp.voz-mtlp.woy  ).* ( llp<0 )   ,3 )),2);
%	lp_lrms= lp_lrms + mean(squeeze(mean( ( mtlp.voz-mtlp.woy  ).* ( llp<-1 )   ,3 )),2);
	
	volfrac_lp_0=volfrac_lp_0+ mean(squeeze(mean(  ( llp<0 )   ,3 )),2);
	volfrac_lp_lrms=volfrac_lp_lrms + mean(squeeze(mean(  ( llp<-1 )   ,3 )),2);



	fht=sprintf("Transfer_fhp_2Dinertial_%03d",time)
        mthp=matfile(fht);
%	hp_0= hp_0 + mean(squeeze(mean( ( mthp.voz-mthp.woy  ).* ( lhp<0 )   ,3 )),2);
%        hp_lrms= hp_lrms + mean(squeeze(mean( ( mthp.voz-mthp.woy  ).* ( lhp<-1 ),3 )),2);

	volfrac_hp_0= volfrac_hp_0 + mean(squeeze(mean(  ( lhp<0 )   ,3 )),2);
        volfrac_hp_lrms= volfrac_hp_lrms + mean(squeeze(mean( ( lhp<-1 ),3 )),2);

end

%m.lp_0=lp_0./nt;
%m.lp_lrms=lp_lrms./nt;
%m.hp_0=hp_0./nt;
%m.hp_lrms=hp_lrms./nt;


m.volfrac_lp_0   =volfrac_lp_0./nt;
m.volfrac_lp_lrms=volfrac_lp_lrms./nt;
m.volfrac_hp_0   =volfrac_hp_0./nt;
m.volfrac_hp_lrms=volfrac_hp_lrms./nt;
