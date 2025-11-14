jcond=130;
Ny=512;
jc=Ny-jcond+1
mp=matfile('../data/mean_profiles.mat');
mvr=matfile('../data/vort_rms.mat');

dUdy=mp.dUdy;

ozrms=mvr.ozrms;
ozrms = 0.5*(ozrms+flipud(ozrms));
ozrms=ozrms(1:256);
ne = 201;
% edges = linspace(-dUdy(jcond)-2*ozrms(jcond),-dUdy(jcond)+2*ozrms(jcond),ne);
e1=-dUdy(jcond)-2*ozrms(jcond);
e2=-dUdy(jcond)+2*ozrms(jcond);
% e1r = round(e1*100)/100;
% e2r = round(e2*100)/100;

edges=[e1r:1e-2:e2r];

nbins = numel(edges)-1;
totalcounts=zeros(1,nbins);
binwidth= edges(4)-edges(3);
tstart=1;
%tend=00000;
tend=10;
%tend=0;
tstep=1;
nf=(tend-tstart)/tstep+1;


for time=tstart:tstep:tend
    time

    fvo=sprintf("/vast/geyink1/skumar67/Ret_1000_data/vort_%03d.mat",time)
	m=matfile(fvo)
    ozb=m.omegaz(jcond,:,:);
    ozt=-m.omegaz(jc,:,:);
    counts=histcounts(ozb,edges);
    totalcounts=totalcounts+counts;
    counts=histcounts(ozt,edges);
    totalcounts=totalcounts+counts;

end
pdfoz = totalcounts./(sum(totalcounts)*binwidth);
binCenters = edges(1:ne-1)+binwidth/2;

 fvpdf=sprintf("../data/pdfoz_j_%03d.mat",jcond)
	m=matfile(fvpdf,'Writable',true);
    m.pdfoz=pdfoz;
    m.binCenters=binCenters;
    m.edges=edges;
    m.ozrms=ozrms(jcond);
    m.dUdy=dUdy(jcond);
