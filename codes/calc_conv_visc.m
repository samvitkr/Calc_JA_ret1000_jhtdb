Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nt=3;
nu=5e-5;
tstart=1;
tend=5;
for time=tstart:tend
    
fgx=sprintf("../data/velgradx_%03d.mat",time);
fgy=sprintf("../data/velgrady_%03d.mat",time);
fgz=sprintf("../data/velgradz_%03d.mat",time);

mgx=matfile(fgx);
mgy=matfile(fgy);
mgz=matfile(fgz);   

fo=sprintf("../data/vort_%03d.mat",time);
mo=matfile(fo);
   
fv=sprintf("../data/velfieldpar_%02d.mat",time);
mv=matfile(fv);

ft=sprintf("../data/Transfer_%03d.mat",time)
mt=matfile(ft,'Writable',true);


voz=single(zeros(Ny,Nx,Nz));
woy=single(zeros(Ny,Nx,Nz));
visc=single(zeros(Ny,Nx,Nz));


mt.voz=(mv.vfield).*(mo.omegaz);
mt.woy=(mv.wfield).*(mo.omegay);
mt.visc=nu*((mgx.d2udx2)+(mgy.d2udy2)+(mgz.d2udz2));


end
