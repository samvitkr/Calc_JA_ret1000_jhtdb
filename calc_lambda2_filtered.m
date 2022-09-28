Nx=2048;
Ny=512;
Nz=1536;

Lx=  8*pi;
Lz = 3*pi;
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];

%kxint=floor(kx(end));
xp = [0:Nx-1]*Lx/(Nx);
load('bsplinedata.mat');

nproc=8;
jstart=1;
jend=180;
Nj=jend-jstart+1;
tstart=23;
tend=35;


yl=yv(jstart:jend)+1;
Delta=15;
Deltaz=1.25;
G=zeros(Nj,Nx);
Gz=zeros(Nj,Nz);

for time=tstart:tend
	S_11 = single(zeros(Nj,Nx,Nz));
	S_12 = single(zeros(Nj,Nx,Nz));
	S_13 = single(zeros(Nj,Nx,Nz));
	S_22 = single(zeros(Nj,Nx,Nz));
	S_23 = single(zeros(Nj,Nx,Nz));
	S_33 = single(zeros(Nj,Nx,Nz));
	
	O_21 = single(zeros(Nj,Nx,Nz));
	O_13 = single(zeros(Nj,Nx,Nz));
	O_32 = single(zeros(Nj,Nx,Nz));

	lambda2 = single(zeros(Nj,Nx,Nz));
	Q = single(zeros(Nj,Nx,Nz));
	
	fvgx=sprintf('velgradx_%03d.mat',time);
	fvgy=sprintf('velgrady_%03d.mat',time);
	fvgz=sprintf('velgradz_%03d.mat',time);
	fo=sprintf('vort_%03d',time);
	

	mvgx=matfile(fvgx);
	mvgy=matfile(fvgy);
	mvgz=matfile(fvgz);
	mo=matfile(fo);

	time
	tic
	S_11=      mvgx.dudx(jstart:jend,:,:);
	S_12=0.5*( mvgy.dudy(jstart:jend,:,:) + mvgx.dvdx(jstart:jend,:,:) );
	S_13=0.5*( mvgz.dudz(jstart:jend,:,:) + mvgx.dwdx(jstart:jend,:,:) );
	S_22=      mvgy.dvdy(jstart:jend,:,:);
	S_23=0.5*( mvgy.dwdy(jstart:jend,:,:) + mvgz.dvdz(jstart:jend,:,:) );
	S_33=      mvgz.dwdz(jstart:jend,:,:);
 		
	O_21 = 0.5*mo.omegaz(jstart:jend,:,:);
	O_13 = 0.5*mo.omegay(jstart:jend,:,:);
	O_32 = 0.5*mo.omegax(jstart:jend,:,:);	

	toc

	O = zeros(3,3);
	S = zeros(3,3);

	for j =1:Nj

	G(j,:)	= exp(-( yl(j).^2*(kx.^2)./Delta^2));
	Gz(j,:)= exp(-( yl(j).^2*(kz.^2)./Deltaz^2));



	S_11(j,:,:)=ifft(  fft(squeeze( S_11(j,:,:))).*(G(j,:)'));
	S_12(j,:,:)=ifft(  fft(squeeze( S_12(j,:,:))).*(G(j,:)'));
	S_13(j,:,:)=ifft(  fft(squeeze( S_13(j,:,:))).*(G(j,:)'));
	S_22(j,:,:)=ifft(  fft(squeeze( S_22(j,:,:))).*(G(j,:)'));
	S_23(j,:,:)=ifft(  fft(squeeze( S_23(j,:,:))).*(G(j,:)'));
	S_33(j,:,:)=ifft(  fft(squeeze( S_33(j,:,:))).*(G(j,:)'));

	O_21(j,:,:)=ifft(  fft(squeeze( O_21(j,:,:))).*(G(j,:)'));
	O_13(j,:,:)=ifft(  fft(squeeze( O_13(j,:,:))).*(G(j,:)'));
	O_32(j,:,:)=ifft(  fft(squeeze( O_32(j,:,:))).*(G(j,:)'));


	S_11(j,:,:)=ifft(  fft(squeeze( S_11(j,:,:))').*(Gz(j,:)'))';
	S_12(j,:,:)=ifft(  fft(squeeze( S_12(j,:,:))').*(Gz(j,:)'))';
	S_13(j,:,:)=ifft(  fft(squeeze( S_13(j,:,:))').*(Gz(j,:)'))';
	S_22(j,:,:)=ifft(  fft(squeeze( S_22(j,:,:))').*(Gz(j,:)'))';
	S_23(j,:,:)=ifft(  fft(squeeze( S_23(j,:,:))').*(Gz(j,:)'))';
	S_33(j,:,:)=ifft(  fft(squeeze( S_33(j,:,:))').*(Gz(j,:)'))';
	
	O_21(j,:,:)=ifft(  fft(squeeze( O_21(j,:,:))').*(Gz(j,:)'))';
	O_13(j,:,:)=ifft(  fft(squeeze( O_13(j,:,:))').*(Gz(j,:)'))';
	O_32(j,:,:)=ifft(  fft(squeeze( O_32(j,:,:))').*(Gz(j,:)'))';


		for k =1:Nz
		for i =1:Nx
		
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
	fl=sprintf("lambda_filtered_15x1p25z_%03d",time)
	ml=matfile(fl,'Writable',true);
	ml.lambda2=single(lambda2);
	ml.Q=single(Q);

	toc
end

