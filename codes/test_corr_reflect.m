close all
clear
m1=matfile('../data/corr_v_j_071.mat')

m2=matfile('../data/corr_v_reflect_j_071.mat')

jtest=100;

% figure
% subplot(2,1,1)
% pcolor(m1.Rvu(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvu(:,:,jtest))
% shading flat
% colorbar
% title('Rvu')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvv(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvv(:,:,jtest))
% shading flat
% colorbar
% title('Rvv')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvw(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvw(:,:,jtest))
% shading flat
% colorbar
% title('Rvw')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdudx(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdudx(:,:,jtest))
% shading flat
% colorbar
% title('Rvdudx')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdvdx(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdvdx(:,:,jtest))
% shading flat
% colorbar
% title('Rvdvdx')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdwdx(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdwdx(:,:,jtest))
% shading flat
% colorbar
% title('Rvdwdx')
% 
% 
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdudy(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdudy(:,:,jtest))
% shading flat
% colorbar
% title('Rvdudy')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdvdy(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdvdy(:,:,jtest))
% shading flat
% colorbar
% title('Rvdvdy')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdwdy(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdwdy(:,:,jtest))
% shading flat
% colorbar
% title('Rvdwdy')
% 
% 
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdudz(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdudz(:,:,jtest))
% shading flat
% colorbar
% title('Rvdudz')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdvdz(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdvdz(:,:,jtest))
% shading flat
% colorbar
% title('Rvdvdz')
% 
% figure
% subplot(2,1,1)
% pcolor(m1.Rvdwdz(:,:,jtest))
% shading flat
% colorbar
% subplot(2,1,2)
% pcolor(m2.Rvdwdz(:,:,jtest))
% shading flat
% colorbar
% title('Rvdwdz')


figure
subplot(2,1,1)
pcolor(m1.Rvvoz(:,:,jtest))
shading flat
colorbar
subplot(2,1,2)
pcolor(m2.Rvvoz(:,:,jtest))
shading flat
colorbar
title('Rvvoz')

figure
subplot(2,1,1)
pcolor(m1.Rvwoy(:,:,jtest))
shading flat
colorbar
subplot(2,1,2)
pcolor(m2.Rvwoy(:,:,jtest))
shading flat
colorbar
title('Rvwoy')