Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nt=5;
tstart=6;
tend=10;
for time=tstart:tend
fgx=sprintf("../data/velgradx_%03d.mat",time);
fgy=sprintf("../data/velgrady_%03d.mat",time);
fgz=sprintf("../data/velgradz_%03d.mat",time);

mgx=matfile(fgx);
mgy=matfile(fgy);
mgz=matfile(fgz);

fo=sprintf("../data/vort_%03d.mat",time)
mo=matfile(fo,'Writable',true);
mo.omegax=single(zeros(Ny,Nx,Nz));
mo.omegay=single(zeros(Ny,Nx,Nz));
mo.omegaz=single(zeros(Ny,Nx,Nz));

mo.omegax=mgy.dwdy-mgz.dvdz;
mo.omegay=mgz.dudz-mgx.dwdx;
mo.omegaz=mgx.dvdx-mgy.dudy;
end
