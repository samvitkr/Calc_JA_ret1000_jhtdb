close all
clear
%jcset=[116, 135, 187, 198, 205];
jcset=171;
for jc=1:1
	jcond=jcset(jc)
	
	nx=512;
	nz=384;
	ny=220;
	ulse=zeros(nz,nx,ny/2);
	vlse=zeros(nz,nx,ny/2);
	wlse=zeros(nz,nx,ny/2);
	fxlse=zeros(nz,nx,ny/2);
	
	lx=4*pi;
	lz=2*pi;
	xp=lx*[0:nx-1]/nx-lx/2;
	zp=lz*[0:nz-1]/nz-lz/2;
	
	load('../data/mean_profiles.mat')
	
	%U=reshape(Um,[1 1 ny]);
	%dUdy=reshape(dUdy,[1 1 ny]);
	%viscm=reshape(viscm,[ 1 1 ny]);
	%U=U(:,:,111:end);
	%dUdy=dUdy(:,:,111:end);
	%viscm=viscm(:,:,111:end);
	
	fp=sprintf('../data/jointpdf_uv_j_%03d.mat',jcond);
	mp=matfile(fp);
	
	fn=sprintf('../data/velgrad_v_lse_j_%03d.mat',jcond);
	ml=matfile(fn);
	%%
	
	%u=mp.u2;
	v=mp.vmax;
	
	   ulse=fftshift(fftshift(v*ml.L11,1),2);...+w*ml.L13;
	   vlse=fftshift(fftshift(v*ml.L21,1),2);...+w*ml.L23;
	   wlse=fftshift(fftshift(v*ml.L31,1),2);...+w*ml.L33;
	
	dudxlse=fftshift(fftshift(v*ml.L41,1),2);
	dvdxlse=fftshift(fftshift(v*ml.L51,1),2);
	dwdxlse=fftshift(fftshift(v*ml.L61,1),2);
	
	dudylse=fftshift(fftshift(v*ml.L71,1),2);
	dvdylse=fftshift(fftshift(v*ml.L81,1),2);
	dwdylse=fftshift(fftshift(v*ml.L91,1),2);
	
	dudzlse=fftshift(fftshift(v*ml.L101,1),2);
	dvdzlse=fftshift(fftshift(v*ml.L111,1),2);
	dwdzlse=fftshift(fftshift(v*ml.L121,1),2);
	
	fxlse =fftshift(fftshift(v*ml.L131,1),2);
	vozlse=fftshift(fftshift(v*ml.L141,1),2);
	woylse=fftshift(fftshift(v*ml.L151,1),2);
	
	yp=yCheb+1;
	[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
	fn=sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
	m=matfile(fn,'Writable',true);
	%m.ucond=u;
	m.vcond=v;
	m.u=ulse;
	m.v=vlse;
	m.w=wlse;
	m.dudx=dudxlse;
	m.dvdx=dvdxlse;
	m.dwdx=dwdxlse;
	m.dudy=dudylse;
	m.dvdy=dvdylse;
	m.dwdy=dwdylse;
	m.dudz=dudzlse;
	m.dvdz=dvdzlse;
	m.dwdz=dwdzlse;
	m.X=X;
	m.Y=Y;
	m.Z=Z;
	m.fx=fxlse;
	m.voz=vozlse;
	m.woy=woylse;
	
	%%
	
	%u=mp.u4;	
	v=mp.vmin;
	
	ulse   =fftshift(fftshift(v*ml.L11,1),2);...+w*ml.L13;
	vlse   =fftshift(fftshift(v*ml.L21,1),2);...+w*ml.L23;
	wlse   =fftshift(fftshift(v*ml.L31,1),2);...+w*ml.L33;
	
	dudxlse=fftshift(fftshift(v*ml.L41,1),2);
	dvdxlse=fftshift(fftshift(v*ml.L51,1),2);
	dwdxlse=fftshift(fftshift(v*ml.L61,1),2);
	
	dudylse=fftshift(fftshift(v*ml.L71,1),2);
	dvdylse=fftshift(fftshift(v*ml.L81,1),2);
	dwdylse=fftshift(fftshift(v*ml.L91,1),2);
	
	dudzlse=fftshift(fftshift(v*ml.L101,1),2);
	dvdzlse=fftshift(fftshift(v*ml.L111,1),2);
	dwdzlse=fftshift(fftshift(v*ml.L121,1),2);
	
	fxlse  =fftshift(fftshift(v*ml.L131,1),2);
	vozlse  =fftshift(fftshift(v*ml.L141,1),2);
	woylse  =fftshift(fftshift(v*ml.L151,1),2);
	
	%%
	yp=yCheb+1;
	[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
	fn=sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
	m=matfile(fn,'Writable',true);
	%m.ucond=u;
	m.vcond=v;
	m.u=ulse;
	m.v=vlse;
	m.w=wlse;
	m.dudx=dudxlse;
	m.dvdx=dvdxlse;
	m.dwdx=dwdxlse;
	m.dudy=dudylse;
	m.dvdy=dvdylse;
	m.dwdy=dwdylse;
	m.dudz=dudzlse;
	m.dvdz=dvdzlse;
	m.dwdz=dwdzlse;
	m.X=X;
	m.Y=Y;
	m.Z=Z;
	m.fx=fxlse;
	m.voz=vozlse;
	m.woy=woylse;
end
