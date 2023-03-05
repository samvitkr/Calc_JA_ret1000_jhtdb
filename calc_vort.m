Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nt=3;
tstart=35;
tend=35;
for time=tstart:tend
fgx=sprintf("velgradx_%03d.mat",time);
fgy=sprintf("velgrady_%03d.mat",time);
fgz=sprintf("velgradz_%03d.mat",time);

mgx=matfile(fgx);
mgy=matfile(fgy);
mgz=matfile(fgz);

fo=sprintf("vort_%03d.mat",time)
mo=matfile(fo,'Writable',true);
mo.omegax=single(zeros(Ny,Nx,Nz));
mo.omegay=single(zeros(Ny,Nx,Nz));
mo.omegaz=single(zeros(Ny,Nx,Nz));

mo.omegax=mgy.dwdy-mgz.dvdz;
mo.omegay=mgz.dudz-mgx.dwdx;
mo.omegaz=mgx.dvdx-mgy.dudy;
end
