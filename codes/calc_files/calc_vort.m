Nx=2048;
Ny=512;
Nz=1536;
nproc=6;
nzproc=Nz/nproc;
nt=5;
tstart=3;
tend=4;
baseDir = "/vast/geyink1/skumar67/Ret_1000_data";
for time=tstart:tend
fgx=sprintf("velgradx_%03d.mat",time);
fgy=sprintf("velgrady_%03d.mat",time);
fgz=sprintf("velgradz_%03d.mat",time);

mgx=matfile(fullfile(baseDir,fgx));
mgy=matfile(fullfile(baseDir,fgy));
mgz=matfile(fullfile(baseDir,fgz));

fo=sprintf("vort_%03d.mat",time)
mo=matfile(fullfile(baseDir,fo),'Writable',true);
mo.omegax=single(zeros(Ny,Nx,Nz));
mo.omegay=single(zeros(Ny,Nx,Nz));
mo.omegaz=single(zeros(Ny,Nx,Nz));

mo.omegax=mgy.dwdy-mgz.dvdz;
mo.omegay=mgz.dudz-mgx.dwdx;
mo.omegaz=mgx.dvdx-mgy.dudy;

oz = mo.omegaz(105,:,:);
zero_fraction = sum(oz == 0, 'all') / numel(oz);
        if zero_fraction > 0.01 % If >1% of data is exact zeros
            fprintf('  -> WARNING: Skipping j=%d at time %d. Data is corrupted (%.1f%% missing).\n', j, time, zero_fraction * 100);
            continue;
        end

end
