clear
close all
Nx=2048;
Ny=512;
Nz=1536;
uin=zeros(Ny,Nz);

for time=2:4
	fvel=sprintf("vel_%03d.mat",time)
	mv=matfile(fvel);	
	uin(:,:)=squeeze(mv.u(:,1,:));
	fuin=sprintf("uinflow_%03d.mat",time);
	muin=matfile(fuin,'Writable',true);
	muin.uin=uin;
end

