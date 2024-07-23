tstart=14;
tend=14;
for time=tstart:tend
        fvel=sprintf("vel_%03d.mat",time)
        mv=matfile(fvel);
	umslice=mean(mv.u,3);
	fvm=sprintf("umslice_%03d.mat",time)
        mvg=matfile(fvm,'Writable',true);
	mvg.um=umslice;
end
