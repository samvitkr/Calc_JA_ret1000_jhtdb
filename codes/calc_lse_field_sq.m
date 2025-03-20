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

fn=sprintf('../data/lse_coeff_sq_j_%03d.mat',jcond);
ml=matfile(fn)

fn=sprintf('../data/lsevp_sq_field_j_%03d.mat',jcond)
m=matfile(fn,'Writable',true);

v=vrms;

m.u2lse=fftshift(fftshift(v*ml.L11,1),2);...+w*ml.L13;
m.v2lse=fftshift(fftshift(v*ml.L21,1),2);...+w*ml.L23;
m.w2lse=fftshift(fftshift(v*ml.L31,1),2);...+w*ml.L33;

m.ox2lse=fftshift(fftshift(v*ml.L41,1),2);
m.oy2lse=fftshift(fftshift(v*ml.L51,1),2);
m.oz2lse=fftshift(fftshift(v*ml.L61,1),2);

