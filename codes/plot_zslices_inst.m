close all
clear
dnu=1.0006e-3;
jcond=105;
jc=jcond;
%figure
for counter = 6:6
%fvgp=sprintf('../data/conditionalp_dudy_zslice_j_%03d_counter_%03d.mat',jcond,counter)
fvgp=sprintf('../data/conditionalp_dudyinst_jcond_%03d_counter_%03d.mat',jcond,counter)
m=matfile(fvgp);
%ysp=m.Y(jcond,10)./dnu;
%subplot(3,4,counter)
figure
hold on
% pcolor(m.X./dnu,m.Y./dnu,m.ozd)
% contour(m.X./dnu,m.Y./dnu,m.qd,[0.1 0.5],'k')
% scatter(0,ysp,'black','filled','o','LineWidth',1)
pcolor(m.Q(:,:,105).*(m.Q(:,:,105)>0))
shading flat

hold off

%axis equal
%ylim([ysp-100 ysp+100])
%xlim([-100 100])
clim([-20 20])
%xlabel('x^+')
%ylabel('y^+')

shading flat
colormap redblue
end
c=colorbar;
ylabel(c,'\omega_z')
sgtitle('Outflow AND dudy<0')


%%
figure
subplot(1,2,1)
fvgp=sprintf('../data/conditionalp_dudy_zslice_j_%03d.mat',jcond)
m=matfile(fvgp);
ysp=m.Y(jcond,10)./dnu;
hold on
pcolor(m.X./dnu,m.Y./dnu,m.ozd)
contour(m.X./dnu,m.Y./dnu,m.qd,[0.1 0.5],'k')
scatter(0,ysp,'black','filled','o','LineWidth',1)
hold off

axis equal
ylim([ysp-100 ysp+100])
xlim([-100 100])
clim([-10 10])
xlabel('x^+')
ylabel('y^+')

shading flat
colormap redblue

c=colorbar;
ylabel(c,'\omega_z')
title('Outflow AND dudy<0')

fvgp=sprintf('../data/conditionaln_dudy_zslice_j_%03d.mat',jcond)
m=matfile(fvgp);
ysp=m.Y(jcond,10)./dnu;
subplot(1,2,2)
hold on
pcolor(m.X./dnu,m.Y./dnu,m.ozd)
contour(m.X./dnu,m.Y./dnu,m.qd,[0.1 0.5],'k')
scatter(0,ysp,'black','filled','o','LineWidth',1)
hold off

axis equal
ylim([ysp-100 ysp+100])
xlim([-100 100])
clim([-10 10])
xlabel('x^+')
ylabel('y^+')

shading flat
colormap redblue

c=colorbar;
ylabel(c,'\omega_z')
title('Inflow AND dudy<0')