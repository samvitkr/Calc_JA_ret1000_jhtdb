clear 
close all
Nx=512;
Nz=384;
Ny=220;
jcond=171
tstart=0000;
tend=100000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
nbin=256;
umin=-0.25;
umax=0.25;
vmin=-0.25;
vmax=0.25;
density=zeros(nbin,nbin);

%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfields_%07d.mat",time);
         m=matfile(fvel);
         ul=reshape( m.ufield(:,:,jcond)-mean(m.ufield(:,:,jcond),"all") ,Nx*Nz,1 );
         vl=reshape( m.vfield(:,:,jcond) ,Nx*Nz,1 );
        uv=[ul,vl];
        [bandwidth,densityinst,ubin,vbin]=kde2d(uv,nbin,[umin vmin],[umax vmax]);
        density=density+densityinst;
end
density=density./nf;

fn=sprintf('../data/jointpdf_uv_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.density=density;
mf.bandwidth=bandwidth;
mf.ubin=ubin;
mf.vbin=vbin;
