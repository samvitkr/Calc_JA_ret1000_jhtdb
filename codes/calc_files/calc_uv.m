Nx=2048;
Ny=512;
Nz=1536;
tstart=2;
tend=3;
uv=single(zeros(Ny,Nx,Nz));
for time=tstart:tend
	
	fn=sprintf("vel_%03d",time);
     fuvn=sprintf("uv_%03d",time)
	mfn=matfile(fn);
	mfuv=matfile(fuvn,'Writable',true);
	mfuv.uv=single( (mfn.u).*(mfn.v) );


end
