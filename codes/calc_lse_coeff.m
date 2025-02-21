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
jcond=41;

fn=sprintf('../data/lse_coeff_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);



mf.L11=            single(zeros(Nz,Nx, 2));
mf.L21=            single(zeros(Nz,Nx, 2));
mf.L31=            single(zeros(Nz,Nx, 2));
mf.L41=            single(zeros(Nz,Nx, 2));
mf.L51=            single(zeros(Nz,Nx, 2));
mf.L61=            single(zeros(Nz,Nx, 2));
mf.L71=            single(zeros(Nz,Nx, 2));
mf.L81=            single(zeros(Nz,Nx, 2));
mf.L91=            single(zeros(Nz,Nx, 2));
mf.L101=           single(zeros(Nz,Nx, 2));
mf.L111=           single(zeros(Nz,Nx, 2));
mf.L121=           single(zeros(Nz,Nx, 2));
mf.L131=           single(zeros(Nz,Nx, 2));
mf.L141=           single(zeros(Nz,Nx, 2));



L11=		single(zeros(Nz,Nx));
L21=		single(zeros(Nz,Nx));
L31=		single(zeros(Nz,Nx));
L41=		single(zeros(Nz,Nx));
L51=		single(zeros(Nz,Nx));
L61=		single(zeros(Nz,Nx));
L71=		single(zeros(Nz,Nx));
L81=		single(zeros(Nz,Nx));
L91=		single(zeros(Nz,Nx));
L101=		single(zeros(Nz,Nx));
L111=		single(zeros(Nz,Nx));
L121=		single(zeros(Nz,Nx));
L131=		single(zeros(Nz,Nx));
L141=		single(zeros(Nz,Nx));



%m=matfile(fn);
fn=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
%m=matfile(fn);
load(fn);
%jc=j;
uij=[Rvv(1,1,jcond)];...,Rvw(1,1,jc)];...
    %Rwu(1,1,jc),Rwv(1,1,jc),Rww(1,1,jc)];
uij=uij.';

for j=1:Ny/2
    j
%L11=            single(zeros(Nz,Nx);
%L21=            single(zeros(Nz,Nx);
%L31=            single(zeros(Nz,Nx);
%L41=            single(zeros(Nz,Nx);
%L51=            single(zeros(Nz,Nx);
%L61=            single(zeros(Nz,Nx);
%L71=            single(zeros(Nz,Nx);
%L81=            single(zeros(Nz,Nx);
%L91=            single(zeros(Nz,Nx);
%L101=           single(zeros(Nz,Nx);
%L111=           single(zeros(Nz,Nx);
%L121=           single(zeros(Nz,Nx);
%L131=           single(zeros(Nz,Nx);
%L141=           single(zeros(Nz,Nx);
    for k=1:Nz
%	    k
        for i=1:Nx
            Rij=[Rvu(k,i,j),Rvv(k,i,j),Rvw(k,i,j),...
		Rvdudx(k,i,j),Rvdvdx(k,i,j),Rvdwdx(k,i,j),...
                Rvdudy(k,i,j),Rvdvdy(k,i,j),Rvdwdy(k,i,j),...
                Rvdudz(k,i,j),Rvdvdz(k,i,j),Rvdwdz(k,i,j),...
		Rvvoz(k,i,j),Rvwoy(k,i,j)];
               % Rwu(k,i,j),Rwv(k,i,j),Rww(k,i,j)];
            Rij=Rij.';
            L=Rij/(uij);
            L11(k,i)=	single(real(L(1,1)));
            L21(k,i)=	single(real(L(2,1)));
            L31(k,i)=	single(real(L(3,1)));
            L41(k,i)=	single(real(L(4,1)));
	    L51(k,i)=	single(real(L(5,1)));
	    L61(k,i)=	single(real(L(6,1)));
	    L71(k,i)=	single(real(L(7,1)));
	    L81(k,i)=	single(real(L(8,1)));
	    L91(k,i)=	single(real(L(9,1)));
	    L101(k,i)=	single(real(L(10,1)));
            L111(k,i)=	single(real(L(11,1)));
            L121(k,i)=	single(real(L(12,1)));
	    L131(k,i)=	single(real(L(13,1)));
	    L141(k,i)=	single(real(L(14,1)));

        end
    end
	 mf.L11(1:Nz,1:Nx,j)=L11; 
         mf.L21(1:Nz,1:Nx,j)=L21; 
         mf.L31(1:Nz,1:Nx,j)=L31; 
         mf.L41(1:Nz,1:Nx,j)=L41; 
         mf.L51(1:Nz,1:Nx,j)=L51; 
         mf.L61(1:Nz,1:Nx,j)=L61; 
         mf.L71(1:Nz,1:Nx,j)=L71; 
         mf.L81(1:Nz,1:Nx,j)=L81; 
         mf.L91(1:Nz,1:Nx,j)=L91; 
        mf.L101(1:Nz,1:Nx,j)=L101; 
        mf.L111(1:Nz,1:Nx,j)=L111; 
        mf.L121(1:Nz,1:Nx,j)=L121; 
        mf.L131(1:Nz,1:Nx,j)=L131; 
        mf.L141(1:Nz,1:Nx,j)=L141;

		
end


%fn=sprintf('../data/lse_coeff_j_%03d.mat',jcond);
%mf=matfile(fn,"Writable",true);
%mf.L11=real(L11);
%mf.L21=real(L21);
%mf.L31=real(L31);
%mf.L41=real(L41);
%mf.L51=real(L51);
%mf.L61=real(L61);
%mf.L71=real(L71);
%mf.L81=real(L81);
%mf.L91=real(L91);
%mf.L101=real(L101);
%mf.L111=real(L111);
%mf.L121=real(L121);
%mf.L131=real(L131);
%mf.L141=real(L141);
