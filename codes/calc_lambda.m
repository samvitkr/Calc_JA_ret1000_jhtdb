clear
mvel1_a = matfile('../turbmat_master_par_1/Full_velfield.mat');
mvel1_b = matfile('../turbmat_master_par_1/Full_velfield_2.mat');
Nx = length(mvel1_a.xp);
Ny = length(mvel1_a.yp);
Nz = length(mvel1_a.zp);
%S_11 = single(zeros(Ny,Nx,Nz));
%S_12 = single(zeros(Ny,Nx,Nz));
%S_13 = single(zeros(Ny,Nx,Nz));
%S_22 = single(zeros(Ny,Nx,Nz));
%S_23 = single(zeros(Ny,Nx,Nz));
%S_33 = single(zeros(Ny,Nx,Nz));
lambda2 = single(zeros(Ny,Nx,Nz));

mvelgz=matfile('Full_velgrad_z.mat');
mo=matfile('Full_vorticity.mat');

mvelg1_a = matfile('../turbmat_master_par_1/Full_velgradfield.mat');
mvelg1_b = matfile('../turbmat_master_par_1/Full_velgradfield_2.mat');

mvelg2_a = matfile('../turbmat_master_par_2/Full_velgradfield.mat');
mvelg2_b = matfile('../turbmat_master_par_2/Full_velgradfield_2.mat');

mvelg3_a = matfile('../turbmat_master_par_3/Full_velgradfield.mat');
mvelg3_b = matfile('../turbmat_master_par_3/Full_velgradfield_2.mat');

mvelg4_a = matfile('../turbmat_master_par_4/Full_velgradfield_1.mat');
mvelg4_b = matfile('../turbmat_master_par_4/Full_velgradfield_2.mat');

O = zeros(3,3);
S = zeros(3,3);
%O_21 = 0.5*mo.omega_z;
%O_13 = 0.5*mo.omega_y;
%O_32 = 0.5*mo.omega_x;
%S_33 = mvelgz.dwdz;
%
for proc=1:4

	switch proc
	case 1
	ma=mvelg1_a;
	mb=mvelg1_b;
	case 2
	ma=mvelg2_a;
	mb=mvelg2_b;
	case 3
	ma=mvelg3_a;
	mb=mvelg3_b;
	case 4
	ma=mvelg4_a;
	mb=mvelg4_b;
	end


	ma
	mb
%	proc
	
	kstart=(proc-1)*384;

%	S_11(:,:,kstart+1:kstart+259) = ma.dudx(:,:,1:259);
%	S_12(:,:,kstart+1:kstart+259) = 0.5*( ma.dudy(:,:,1:259) +ma.dvdx(:,:,1:259) );
%	S_13(:,:,kstart+1:kstart+259) = 0.5*( mvelgz.dudz(:,:,kstart+1:kstart+259)+ma.dwdx(:,:,1:259));
%	S_22(:,:,kstart+1:kstart+259) = ma.dvdy(:,:,1:259)
%	S_23(:,:,kstart+1:kstart+259) = 0.5*( mvelgz.dvdz(:,:,kstart+1:kstart+259)+ma.dwdy(:,:,1:259));


%	S_11(:,:,kstart+260:kstart+384) = mb.dudx(:,:,260:384);
%        S_12(:,:,kstart+260:kstart+384) = 0.5*( mb.dudy(:,:,260:384) +mb.dvdx(:,:,260:384) );
%        S_13(:,:,kstart+260:kstart+384) = 0.5*( mvelgz.dudz(:,:,kstart+260:kstart+384)+mb.dwdx(:,:,260:384));
%        S_22(:,:,kstart+260:kstart+384) = mb.dvdy(:,:,260:384)
%        S_23(:,:,kstart+260:kstart+384) = 0.5*( mvelgz.dvdz(:,:,kstart+260:kstart+384)+mb.dwdy(:,:,260:384));
	

 
	tic

	for i =1:Nx
		for j = 1:Ny
		for k = 1:384

		if(k<260) mvg=ma;
		else mvg=mb;
		end

		S(1,1) = mvg.dudx(i,j,k);
		S(1,2) = 0.5*( mvg.dudy(i,j,k) +mvg.dvdx(i,j,k) );
		S(1,3) = 0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
		S(2,1) = S(1,2);
		S(2,2) =  mvg.dvdy(i,j,k);
		S(2,3) = 0.5*( mvelgz.dvdz(i,j,kstart+k)+mvg.dwdy(i,j,k));
		S(3,1) =  0.5*( mvelgz.dudz(i,j,kstart+k)+mvg.dwdx(i,j,k));
		S(3,2) = S(2,3);
		S(3,3) = mvelgz.dwdz(i,j,kstart+k);
		
		O(1,3)= 0.5*mo.omega_y(i,j,kstart+k);
		O(2,1) =0.5*mo.omega_z(i,j,kstart+k);
		O(3,2) = 0.5*mo.omega_x(i,j,kstart+k);
		O(1,2) = -O(2,1);
		O(2,3) = -O(3,2);
                O(3,1) = -O(1,3);

		A = S*S + O*O;
		ll = sort(eig(A));
		lambda2(i,j,kstart+k) = ll(2);

		end
		end
	end
toc
proc
	

end
tic
mw=matfile('Full_lambda2.mat','Writable',true);
mw.lambda2=single(lambda2);
toc


