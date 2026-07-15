jcset=[47 54 71 105 130]
for jj=4:4
    jcond=jcset(jj)

    %fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
    %fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)

    fvgn=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);
    fvgp=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
    % fvgn=sprintf("../data/conditionaln_jcond_inst_%03d_%02d.mat",jcond,counter);

fc = sprintf("../data/conditionalp_bottom_split_allvars_jcond_%03d.mat", jcond);
mc = matfile(fc, 'Writable', true);


    % m1=matfile(fvgp,'Writable',true)
    % m2=matfile(fvgn,'Writable',true)
    %mvg=matfile(fvgp,"Writable",true)
mvg=mc;

    %m1=matfile('../data/test.mat','Writable',true)
    %fvgq=[fvgp fvgn];
    for nn=1:2
        % fvgp=sprintf("../data/conditionalp_dudyinst_jcond_%03d_counter_%03d.mat",jcond,nn);
    	switch nn
        	case 1
        	dudx=mvg.dudx_neg;
                dvdx=mvg.dvdx_neg;
                dwdx=mvg.dwdx_neg;

                dudy=mvg.dudy_neg;
                dvdy=mvg.dvdy_neg;
                dwdy=mvg.dwdy_neg;

                dudz=mvg.dudz_neg;
                dvdz=mvg.dvdz_neg;
                dwdz=mvg.dwdz_neg;
        	case 2
        	dudx=mvg.dudx_pos;
                dvdx=mvg.dvdx_pos;
                dwdx=mvg.dwdx_pos;

                dudy=mvg.dudy_pos;
                dvdy=mvg.dvdy_pos;
                dwdy=mvg.dwdy_pos;

                dudz=mvg.dudz_pos;
                dvdz=mvg.dvdz_pos;
                dwdz=mvg.dwdz_pos;
            case 3
                dudx=mvg.dudx;
                dvdx=mvg.dvdx;
                dwdx=mvg.dwdx;

                dudy=mvg.dudy;
                dvdy=mvg.dvdy;
                dwdy=mvg.dwdy;

                dudz=mvg.dudz;
                dvdz=mvg.dvdz;
                dwdz=mvg.dwdz;
    	end
        %fvg=fvgoz;
        %fvg=fvgq(nn);
        %mvg=matfile(fvg,'Writable',true)


        [Nz Nx Ny]=size(mvg.u_pos);

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

        S_11=dudx;
        S_12=0.5*(dudy+dvdx);
        S_13=0.5*(dudz+dwdx);
        S_22=dvdy;
        S_23=0.5*(dwdy+dvdz);
        S_33=dwdz;

        omegaZ=dvdx-dudy;
        omegaY=dudz-dwdx;
        omegaX=dwdy-dvdz;
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
        switch nn
            case 1
                mvg.lambda2n=single(lambda2);
                mvg.Qn=single(Q);
            case 2
                mvg.lambda2p=single(lambda2);
                mvg.Qp=single(Q);
            case 3
                mvg.lambda2=single(lambda2);
                mvg.Q=single(Q);
        end
    end

end

