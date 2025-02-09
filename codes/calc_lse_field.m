close all
clear
%jcset=[116, 135, 187, 198, 205];
jcond=71;

load('../data/JHTDB_RET1000.mat')

vrmsprofile=sqrt(JHTDB_RET1000(:,5))./(JHTDB_RET1000(end,2));
vrms=vrmsprofile(jcond)%

%nx=2048;
%nz=1536;
%ny=512;

fn=sprintf('../data/lse_coeff_j_%03d.mat',jcond);
ml=matfile(fn)

fn=sprintf('../data/lsevp_field_j_%03d.mat',jcond)
m=matfile(fn,'Writable',true);

v=vrms;

m.ulse=fftshift(fftshift(v*ml.L11,1),2);...+w*ml.L13;
m.vlse=fftshift(fftshift(v*ml.L21,1),2);...+w*ml.L23;
m.wlse=fftshift(fftshift(v*ml.L31,1),2);...+w*ml.L33;

m.dudxlse=fftshift(fftshift(v*ml.L41,1),2);
m.dvdxlse=fftshift(fftshift(v*ml.L51,1),2);
m.dwdxlse=fftshift(fftshift(v*ml.L61,1),2);

m.dudylse=fftshift(fftshift(v*ml.L71,1),2);
m.dvdylse=fftshift(fftshift(v*ml.L81,1),2);
m.dwdylse=fftshift(fftshift(v*ml.L91,1),2);

m.dudzlse=fftshift(fftshift(v*ml.L101,1),2);
m.dvdzlse=fftshift(fftshift(v*ml.L111,1),2);
m.dwdzlse=fftshift(fftshift(v*ml.L121,1),2);

m.vozlse=fftshift(fftshift(v*ml.L131,1),2);
m.woylse=fftshift(fftshift(v*ml.L141,1),2);

%yp=yCheb+1;
%[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
%fn=sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
%m=matfile(fn,'Writable',true);
%%m.ucond=u;
%m.vcond=v;
%m.u=ulse;
%m.v=vlse;
%m.w=wlse;
%m.dudx=dudxlse;
%m.dvdx=dvdxlse;
%m.dwdx=dwdxlse;
%m.dudy=dudylse;
%m.dvdy=dvdylse;
%m.dwdy=dwdylse;
%m.dudz=dudzlse;
%m.dvdz=dvdzlse;
%m.dwdz=dwdzlse;
%m.X=X;
%m.Y=Y;
%m.Z=Z;
%m.fx=fxlse;
%m.voz=vozlse;
%m.woy=woylse;
%
%%%
%
%%u=mp.u4;	
%v=mp.vmin;
%
%ulse   =fftshift(fftshift(v*ml.L11,1),2);...+w*ml.L13;
%vlse   =fftshift(fftshift(v*ml.L21,1),2);...+w*ml.L23;
%wlse   =fftshift(fftshift(v*ml.L31,1),2);...+w*ml.L33;
%
%dudxlse=fftshift(fftshift(v*ml.L41,1),2);
%dvdxlse=fftshift(fftshift(v*ml.L51,1),2);
%dwdxlse=fftshift(fftshift(v*ml.L61,1),2);
%
%dudylse=fftshift(fftshift(v*ml.L71,1),2);
%dvdylse=fftshift(fftshift(v*ml.L81,1),2);
%dwdylse=fftshift(fftshift(v*ml.L91,1),2);
%
%dudzlse=fftshift(fftshift(v*ml.L101,1),2);
%dvdzlse=fftshift(fftshift(v*ml.L111,1),2);
%dwdzlse=fftshift(fftshift(v*ml.L121,1),2);
%
%fxlse  =fftshift(fftshift(v*ml.L131,1),2);
%vozlse  =fftshift(fftshift(v*ml.L141,1),2);
%woylse  =fftshift(fftshift(v*ml.L151,1),2);
%
%%%
%yp=yCheb+1;
%[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
%fn=sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
%m=matfile(fn,'Writable',true);
%%m.ucond=u;
%m.vcond=v;
%m.u=ulse;
%m.v=vlse;
%m.w=wlse;
%m.dudx=dudxlse;
%m.dvdx=dvdxlse;
%m.dwdx=dwdxlse;
%m.dudy=dudylse;
%m.dvdy=dvdylse;
%m.dwdy=dwdylse;
%m.dudz=dudzlse;
%m.dvdz=dvdzlse;
%m.dwdz=dwdzlse;
%m.X=X;
%m.Y=Y;
%m.Z=Z;
%m.fx=fxlse;
%m.voz=vozlse;
%m.woy=woylse;
