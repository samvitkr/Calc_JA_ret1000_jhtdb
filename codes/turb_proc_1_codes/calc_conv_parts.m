Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nt=3;
nu=5e-5;
tstart=38;
tend=38;
for time=tstart:tend
    
fgx=sprintf("velgradx_%03d.mat",time);
fgy=sprintf("velgrady_%03d.mat",time);
fgz=sprintf("velgradz_%03d.mat",time);

mgx=matfile(fgx);
mgy=matfile(fgy);
mgz=matfile(fgz);   

fo=sprintf("vort_%03d.mat",time);
mo=matfile(fo);
   
fv=sprintf("vel_%03d.mat",time);
mv=matfile(fv);

ft=sprintf("Transfer_%03d.mat",time)
mt=matfile(ft,'Writable',true);

mt.voz=single(zeros(Ny,Nx,Nz));
mt.woy=single(zeros(Ny,Nx,Nz));


mt.voz=(mv.v).*(mo.omegaz);
mt.woy=(mv.w).*(mo.omegay);


end
