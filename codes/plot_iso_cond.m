close all
clear
load('../data/bsplinedata.mat')
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
xp=ret*(lx*[0:nx-1]/nx-lx/2);
zp=ret*(lz*[0:nz-1]/nz-lz/2);
yp=ret*(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;
jcond=130;
yc=yv(jcond)+1;
ut=0.0499;
dnu=1.0006e-3;
jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_%03d.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;
xp1=xp(itarget-wxx:itarget+wxx);
zp1=zp(ktarget-wzz:ktarget+wzz);
% l1=min(m1.lambda2(1:end,1:end,jc),[],'all');
% l2=min(m2.lambda2(1:end,1:end,jc),[],'all');
l1=-ut^2/yc^2;
l2=-ut^2/yc^2;
val=6;
[X,Z,Y]=(meshgrid(xp1,zp1,yp));

% mt=matfile('velgrad_transfer_flp_0070000.mat')
% m=matfile('lambdaflp_0070000.mat');
%mt=matfile('transferfields_0040000.mat')
%m=matfile('lambda_0040000.mat');
% ml=(mean(mean(m.lambda2,1),2));
% l=m.lambda2;
% lrms=rms(rms(m.lambda2-ml,1),2);

x1=150;
y1=150;
x2=1200;
y2=400;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
m=m1;
%polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
nld=(nl).*(m.lambda2./l1>val);
vold=m.lambda2./l1>val;
nldtot=squeeze(mean((nl),[1 2]));
cl=floor(max([abs(min(nld,[],'all')./ut^2), abs(max(nld,[],'all')./ut^2)]));


subplot(1,5,[1 2])
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l1,[2 1 3]),val,permute(nl./ut^2,[2 1 3]))
% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),val,permute(nl,[2 1 3]))
% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),-val,permute(nl,[2 1 3]))

scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading interp
lightangle(45,45)
%camlight('left')
%clim([-1 1]*1e-1)
colorbar 
colormap redblue
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')

xlabel('$z^+$','interpreter','latex','FontSize',11)
ylabel('$x^+$','interpreter','latex','FontSize',11)
zlabel('$y^+$','interpreter','latex','FontSize',11)
c=colorbar;
set(gca,'FontSize',11)
%set(gca,'Position',[0.03 0.15 0.3 0.8])
%ylabel(c,"$(v\omega_z - w\omega_y)/(-u_{\tau}^2/H)$",'interpreter','latex','FontSize',11)

c.Ticks=[-cl, 0, cl];
view(45,45)
clim([-cl cl])
zlim([m.ystart-3 m.yend+3])
xlim([m.zstart-3 m.zend+3])
ylim([m.xstart-3 m.xend+3])

zticks([floor(m.ystart) floor(m.yend)+1])
xticks([floor(m.zstart) ,0,floor(m.zend)+1])
yticks([floor(m.xstart) ,0,floor(m.xend)+1])

% ylim([-100 100])
% zlim([0 200])
% zlim([yp(jc)/2 yp(jc)*1.5])
title('(a)')
subplot(1,5,[3 4])
m=m2;

%polywork=(m.fx).*(m.u)+(m.fy).*(m.v)+(m.fz).*(m.w);
nl=m.woy-m.voz;
ox=m.dwdy-m.dvdz;
nlu=(nl).*(m.lambda2./l1>val);
volu=(m.lambda2./l1>val);
nlutot=squeeze(mean((nl),[1 2]));
cl=floor(max([ abs(min(nlu,[],'all')./ut^2) , abs(max(nlu,[],'all')./ut^2)]));

% h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
hold on
isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(m.lambda2./l2,[2 1 3]),val,permute(nl./ut^2,[2 1 3]))
% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),val,permute(nl,[2 1 3]))
% isosurface(permute(Z,[2 1 3]),permute(X,[2 1 3]),permute(Y,[2 1 3]),permute(nl,[2 1 3]),-val,permute(nl,[2 1 3]))
scatter3(0,0,yp(jc),80,'green','filled')
hold off
%(mt.voz-mt.woy)./(-ut^2) )
%pbaspect([4*pi 2*pi 2])
axis equal
axis tight
shading interp
lightangle(45,45)
%camlight('left')
% clim([-cl cl])
colorbar 
colormap redblue
c=colorbar;
set(gca,'FontSize',11)
%set(gca,'Position',[0.4 0.15 0.3 0.8])

ylabel(c,"$(v\omega_z - w\omega_y)/(-u_{\tau}^2/H)$",'interpreter','latex','FontSize',11)
c.Ticks=[-cl, 0, cl];
view(45,45)
clim([-cl cl])
%print(h1,'isotryflp','-dpng');
%saveas(h1,'iso_lambda_flp_70000.fig')
xlabel('$z^+$','interpreter','latex','FontSize',11)
ylabel('$x^+$','interpreter','latex','FontSize',11)
%zlabel('$y^+$','interpreter','latex','FontSize',11)
view(45,45)
% xlim([-100 100])
zlim([m.ystart-3 m.yend+3])
xlim([m.zstart-3 m.zend+3])
ylim([m.xstart-3 m.xend+3])

zticks([floor(m.ystart) floor(m.yend)+1])
xticks([floor(m.zstart) ,0,floor(m.zend)+1])
yticks([floor(m.xstart) ,0,floor(m.xend)+1])
title('(b)')

% zlim([0 200])
%zlim([30 600])
%% flux contri
G=[0 0.5 0]
ar=lx*lz;
da=ar/(nx*nz);
arp=dnu^2;
load('../data/mean_profiles.mat')
nldav=da*squeeze(sum(nld ,[1 2] ));
nluav=da*squeeze(sum(nlu ,[1 2] ));
voldav=da*squeeze(sum(vold ,[1 2] ));
voluav=da*squeeze(sum(volu ,[1 2] ));
nltot=woym-vozm;
%figure
subplot(1,5,5)
hold on
plot(nldav./(ut^2),yp,'-k','LineWidth',1.5)
plot(nluav./(ut^2),yp,'-.','color',G,'LineWidth',1.5)
yline(yc*ret)
% plot(nldtot./ut^2-nltot./ut^2,yp,'-.r')
% plot(nlutot./ut^2-nltot./ut^2,yp,'-.b')
hold off
ylabel('$y^+$','interpreter','latex')
xlabel('$\langle v\omega_z - w\omega_y \rangle_{x,z}/(-u_{\tau}^2H)$','interpreter','latex')
trapz(yp,nldav)
trapz(yp,nluav)
ylim([ m1.ystart-10 m1.yend+10 ])
yticks( [ floor(m1.ystart) , round(yp(jc)) ,floor(m1.yend)+1] )
xline(0)
title('(c)')
set(gca,'YaxisLocation','right')
%set(gca,'Position',[0.75 0.15 0.2 0.8])

fdne=sprintf('eddy_flux_j_%03d.fig',jcond)
saveas(h1,fdne)


% % grid on
% % subplot(1,2,2)
% % hold on
% % plot(100*voldav./ar,yp,'-r')
% % plot(100*voluav./ar,yp,'-b')
% % yline(yc*ret)
% % % plot(nldtot./ut^2-nltot./ut^2,yp,'-.r')
% % % plot(nlutot./ut^2-nltot./ut^2,yp,'-.b')
% % hold off
% % ylabel('y^+')
% % xlabel('area fraction (%)')
% % trapz(yp,nldav)
% % trapz(yp,nluav)
% % ylim([yp(jc)/2 1.5*yc*ret])
% % grid on