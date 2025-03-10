nt=10;
tstart=1;
tend=10;
ny=512;

oxm=zeros(ny,1);
oym=zeros(ny,1);
ozm=zeros(ny,1);
lambda2m=zeros(ny,1);
oxrms=zeros(ny,1);
oyrms=zeros(ny,1);
ozrms=zeros(ny,1);
lambda2rms=zeros(ny,1);

for time=tstart:tend
	fo=sprintf("../data/vort_%03d.mat",time)
	mo=matfile(fo);
	oxm=oxm		+squeeze(mean(mo.omegax,	[2 3]));
	oym=oym		+squeeze(mean(mo.omegay,	[2 3]));
	ozm=ozm		+squeeze(mean(mo.omegaz,	[2 3]));
	oxrms=oxrms	+squeeze(mean(mo.omegax.^2,	[2 3]));
	oyrms=oyrms	+squeeze(mean(mo.omegay.^2,	[2 3]));
	ozrms=ozrms	+squeeze(mean(mo.omegaz.^2,	[2 3]));
	clear mo
%	fl=sprintf("../data/lambda_%03d",time)
%	ml=matfile(fl);
%	lambda2m=lambda2m	+squeeze(mean(ml.lambda2,   [2 3]));
%	lambda2rms=lambda2rms	+squeeze(mean(ml.lambda2.^2,[2 3]));
%	clear ml

end

oxm=oxm./nt;
oym=oym./nt;
ozm=ozm./nt;
oxrms=oxrms./nt;
oyrms=oyrms./nt;
ozrms=ozrms./nt;
%lambda2m=lambda2m./nt;
%lambda2rms=lambda2rms./nt;

oxrms		=sqrt(abs(oxrms-oxm.^2	));
oyrms		=sqrt(abs(oyrms-oym.^2	));
ozrms		=sqrt(abs(ozrms-ozm.^2	));
%lambda2rms	=sqrt(abs(lambda2rms-lambda2m.^2));

m=matfile("../data/vort_rms.mat",'Writable',true)
m.oxrms=oxrms;
m.oyrms=oyrms;
m.ozrms=ozrms;
%m.lambda2rms=lambda2rms;

