jcond=71;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)

%m1=matfile(fvgp,'Writable',true)
%m2=matfile(fvgn,'Writable',true)

m1=matfile('../data/test.mat','Writable',true)
%fvgq=[fvgp fvgn];
for nn=1:2

	switch nn
	case 1
		mvg=m1;
	%case 2
	%	mvg=m2;
	end
%fvg=fvgoz;
%fvg=fvgq(nn);	
%mvg=matfile(fvg,'Writable',true)


[Nz Nx Ny]=size(mvg.u);

S_11	=single(zeros(Nz,Nx,Ny));
S_12	=single(zeros(Nz,Nx,Ny));
S_13	=single(zeros(Nz,Nx,Ny));
S_22	=single(zeros(Nz,Nx,Ny));
S_23	=single(zeros(Nz,Nx,Ny));
S_33	=single(zeros(Nz,Nx,Ny));
O_21	=single(zeros(Nz,Nx,Ny));
O_13	=single(zeros(Nz,Nx,Ny));
O_32	=single(zeros(Nz,Nx,Ny));
lambda2	=single(zeros(Nz,Nx,Ny));
Q	=single(zeros(Nz,Nx,Ny));

S_11=mvg.dudx;
S_12=0.5*(mvg.dudy+mvg.dvdx);
S_13=0.5*(mvg.dudz+mvg.dwdx);
S_22=mvg.dvdy;
S_23=0.5*(mvg.dwdy+mvg.dvdz);
S_33=mvg.dwdz;

omegaZ=mvg.dvdx-mvg.dudy;
omegaY=mvg.dudz-mvg.dwdx;
omegaX=mvg.dwdy-mvg.dvdz;
O_21 = 0.5*omegaZ(:,:,:);
O_13 = 0.5*omegaY(:,:,:);
O_32 = 0.5*omegaX(:,:,:);

O = zeros(3,3);
S = zeros(3,3);

for j =1:Nz
	j
        for i =1:Nx
        for k =1:Ny
        S(1,1) = S_11(j,i,k);%mvg.dudx(i,j,k);
        S(1,2) = S_12(j,i,k);%0.5*( mvg.dudy(i,j,k) +mvg.dvdx(i,j,k) );
        S(1,3) = S_13(j,i,k);%0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
        S(2,1) = S(1,2);
        S(2,2) = S_22(j,i,k);% mvg.dvdy(i,j,k);
        S(2,3) = S_23(j,i,k);%0.5*( mvelgz.dvdz(i,j,kstart+k)+mvg.dwdy(i,j,k));
        S(3,1) = S(1,3);%,ks 0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
        S(3,2) = S(2,3);
        S(3,3) = S_33(j,i,k);%mvelgz.dwdz(i,j,kstart+k);

        O(1,3) = O_13(j,i,k);%0.5*mo.omega_y(i,j,kstart+k);
        O(2,1) = O_21(j,i,k);%0.5*mo.omega_z(i,j,kstart+k);
        O(3,2) = O_32(j,i,k);% 0.5*mo.omega_x(i,j,kstart+k);
        O(1,2) =-O(2,1);
        O(2,3) =-O(3,2);
        O(3,1) =-O(1,3);

        A = S*S + O*O;
        B = O*O';
	C = S*S';

        ll = sort(eig(A));
        lambda2(j,i,k) = single(ll(2));
        Q(j,i,k)= 0.5*(trace(B)-trace(C));
        end
        end
end
mvg.lambda2=single(lambda2);
mvg.Q=single(Q);
end
