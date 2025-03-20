clear all
clc
clear
%numWorkers = 16; % Use 48 processors on the current node

% Start the parallel pool with the specified number of workers
%parpool('local', numWorkers);

Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
jcond=71;

fn=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
mn=matfile(fn)
Rvv=mn.Rvv
uij=[mn.Rvv(1,1,jcond)];...,Rvw(1,1,jc)];...
    %Rwu(1,1,jc),Rwv(1,1,jc),Rww(1,1,jc)];
uij=uij.';
clear mn

fnrms=sprintf('../data/corr_v_rms_reflect_j_%03d.mat',jcond);
load(fnrms)

fn=sprintf('../data/lse_coeff_sq_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);

mf.L11=            single(zeros(Nz,Nx, 2));
mf.L21=            single(zeros(Nz,Nx, 2));
mf.L31=            single(zeros(Nz,Nx, 2));
mf.L41=            single(zeros(Nz,Nx, 2));
mf.L51=            single(zeros(Nz,Nx, 2));
mf.L61=            single(zeros(Nz,Nx, 2));

L11=		single(zeros(Nz,Nx));
L21=		single(zeros(Nz,Nx));
L31=		single(zeros(Nz,Nx));
L41=		single(zeros(Nz,Nx));
L51=		single(zeros(Nz,Nx));
L61=		single(zeros(Nz,Nx));


for j=1:Ny/2
    j
    for k=1:Nz
        for i=1:Nx
            Rij=[Rvu2(k,i,j),Rvv2(k,i,j),Rvw2(k,i,j),...
		Rvox2(k,i,j),Rvoy2(k,i,j),Rvoz2(k,i,j)];
            Rij=Rij.';
            L=Rij/(uij);
            L11(k,i)=	single(real(L(1,1)));
            L21(k,i)=	single(real(L(2,1)));
            L31(k,i)=	single(real(L(3,1)));
            L41(k,i)=	single(real(L(4,1)));
	    L51(k,i)=	single(real(L(5,1)));
	    L61(k,i)=	single(real(L(6,1)));

        end
    end
	 mf.L11(1:Nz,1:Nx,j)=L11; 
         mf.L21(1:Nz,1:Nx,j)=L21; 
         mf.L31(1:Nz,1:Nx,j)=L31; 
         mf.L41(1:Nz,1:Nx,j)=L41; 
         mf.L51(1:Nz,1:Nx,j)=L51; 
         mf.L61(1:Nz,1:Nx,j)=L61; 

		
end


