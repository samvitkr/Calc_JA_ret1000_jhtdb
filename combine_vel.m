Nx=2048;
Ny=512;
Nz=1536;
nproc=8;
nzproc=Nz/nproc;
nt=3;
kstart=zeros(nproc,1);
for proc=1:nproc
    kstart(proc)=(proc-1)*nzproc;
end

for time=33:35
    u=single(zeros(Ny,Nx,Nz));
    v=single(zeros(Ny,Nx,Nz));
    w=single(zeros(Ny,Nx,Nz));
    
    fn=strings(nproc,1);
    for proc=1:nproc
        
        fn(proc)=sprintf("velfield_%02d_%03d.mat",proc,time)
        m=matfile(fn(proc));
        u(:,:,kstart(proc)+1:kstart(proc)+nzproc)=m.ufield(:,:,1:nzproc);
        v(:,:,kstart(proc)+1:kstart(proc)+nzproc)=m.vfield(:,:,1:nzproc);
        w(:,:,kstart(proc)+1:kstart(proc)+nzproc)=m.wfield(:,:,1:nzproc);
    end
    fvel=sprintf("vel_%03d.mat",time)
    
    mv=matfile(fvel,'Writable',true);
    mv.u=u;
    mv.v=v;
    mv.w=w;


	fuvn=sprintf("uv_%03d",time)
        mfuv=matfile(fuvn,'Writable',true);
	mfuv.uv=single( u.*v );
end

clear all
%calc_uv;clear all
calc_velgrad_x;clear all;
calc_velgrad_y;clear all;
calc_velgrad_z;clear all;
calc_vort;clear all;
calc_conv_visc;clear all;
calc_lambda2;clear all

