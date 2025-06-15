close all
clear
load('../data/bsplinedata.mat')
mpr=matfile('eddy_boxsize_ar_nl.mat','Writable',true)
mpr.nldprofiles(1:51,4)=0;
mpr.nluprofiles(1:51,4)=0;
mpr.nldprofiles(1:81,5)=0;
mpr.nluprofiles(1:81,5)=0;
G=[0 0.5 0];
nx=2048;
nz=1536;
Ny=256;
lx=8*pi;
lz=3*pi;
ret=1000;
xg=ret*(lx*[0:nx-1]/nx-lx/2);
zg=ret*(lz*[0:nz-1]/nz-lz/2);
yp=ret*(yv(1:Ny)'+1);
itarget=nx/2+1;
ktarget=nz/2+1;
ut=0.0499;
dnu=1.0006e-3;
jcset=[47 54 71 105 130];

x1=100;
y1=10;
width=1150;
%height=600;
height=400;
fd=figure;
xb1=-125;
xb2=110;
zb1=-120;
zb2=120;
fd.Position=[x1 y1 width height];
%t=tiledlayout(3,5);
t=tiledlayout(2,5);

t.TileSpacing="compact";
t.Padding="compact";
%jcond=130;
labels=arrayfun(@(k) ['(' char('a'+k-1) ')'],1:15,'UniformOutput',false);

for j=1:5
    jcond=jcset(j);
yc=yv(jcond)+1;

jc=jcond;
%fvgp=sprintf('../data/lsevp_field_tot_j_%03d.mat',jcond)
%fvgn=sprintf('../data/lsevn_field_tot_j_%03d.mat',jcond)
fvgp=sprintf("../data/conditionalp_jcond_1_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m1.u);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;
xp=xg;
zp=zg;
xp1=xg(itarget-wxx:itarget+wxx);
zp1=zg(ktarget-wzz:ktarget+wzz);
yc=yv(jcond)+1;

ys=yp(jcond);

fvgp=sprintf("../data/conditionalp_jcond_1_%03d.mat",jcond);
fvgn=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
m1=matfile(fvgp,'Writable',true);
m2=matfile(fvgn,'Writable',true);
[nzz, nxx, nyy]=size(m1.lambda2);
wzz=(nzz-1)/2;
wxx=(nxx-1)/2;

xp=xp(itarget-wxx:itarget+wxx);
zp=zp(ktarget-wzz:ktarget+wzz);
% l1=min(m1.lambda2(1:end,1:end,jcond),[],'all');
% l2=min(m2.lambda2(1:end,1:end,jcond),[],'all');
% val=0.025;
l1=-ut^2/yc^2;
l2=-ut^2/yc^2;
val=6;
[X,Z,Y]=(meshgrid(xp,zp,yp));

x=permute(X,[2 1 3]);
z=permute(Z,[2 1 3]);
y=permute(Y,[2 1 3]);

oxd=permute(m1.dwdy-m1.dvdz,[2 1 3]);
oyd=permute(m1.dudz-m1.dwdx,[2 1 3]);
ozd=permute(m1.dvdx-m1.dudy,[2 1 3]);
nld=permute(m1.voz-m1.woy,[2 1 3]);
ld=permute(m1.lambda2,[2 1 3])./l1;

oxu=permute(m2.dwdy-m2.dvdz,[2 1 3]);
oyu=permute(m2.dudz-m2.dwdx,[2 1 3]);
ozu=permute(m2.dvdx-m2.dudy,[2 1 3]);
lu=permute(m2.lambda2,[2 1 3])./l2;
nlu=permute(m2.voz-m2.woy,[2 1 3]);
ltd=val;
ltu=val;
yb1=min([m1.ystart m2.ystart]);
yb2=max([m1.yend m2.yend]);
%%
%subplot(3,5,j)
nexttile
hold on;
%
isosurf=isosurface(z,x,y,ld,ltd);
interpColors = interp3(z, x, y, dnu*oxd./ut,...wd.*0, ...interpColors = interp3(z, x, y, -nld./ut^2,...wd.*0, ...
    isosurf.vertices(:, 1), ...
    isosurf.vertices(:, 2), ...
    isosurf.vertices(:, 3));
cld=floor(max(100*interpColors,[],'all'))/100

p = patch(isosurf);
p.FaceColor = 'interp';        % Interpolated color
p.EdgeColor = 'none';          % Remove edges
p.FaceVertexCData = interpColors;    % Assign interpolated colors
p.FaceAlpha = 0.8;  
camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,30,'green','filled')
view(45,45);
axis equal

zlim([yb1-10 yb2+10])
xlim([zb1-1 zb2+1])
ylim([xb1-1 xb2+1])

 zticks([round(m1.ystart) round(m1.yend)])
 xticks([zb1,0,zb2 ])
 yticks([xb1,0,xb2 ])

 % xticks([round(m1.zstart) ,0,round(m1.zend)])
% yticks([round(m1.xstart) ,0,round(m1.xend)])
clim([-0.1 0.1])
colormap saffrongreen

if j==5
 c=colorbar;
 ylabel(c,"$\langle \omega_x^+ \rangle_+$",'interpreter','latex','FontSize',11)
c.Location='eastoutside';
 end
xlabel('$z^+$','interpreter','latex','FontSize',11)
ylabel('$x^+$','interpreter','latex','FontSize',11)
if j==1
    zlabel('$y^+$','interpreter','latex','FontSize',11)
    
end
%title('(a)')
title( labels(j) )
%%


%subplot(3,5,j+5)
nexttile(j+5)

hold on;
%
isosurf=isosurface(z,x,y,lu,ltu);
interpColors = interp3(z, x, y, dnu*oxu./ut,...wd.*0, ...interpColors = interp3(z, x, y, -nld./ut^2,...wd.*0, ...
    isosurf.vertices(:, 1), ...
    isosurf.vertices(:, 2), ...
    isosurf.vertices(:, 3));
clu=floor(max(100*interpColors,[],'all'))/100

p = patch(isosurf);
p.FaceColor = 'interp';        % Interpolated color
p.EdgeColor = 'none';          % Remove edges
p.FaceVertexCData = interpColors;    % Assign interpolated colors
p.FaceAlpha = 0.8;  
camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
scatter3(0,0,ys,30,'green','filled')
view(135,45);
axis equal

% zlim([yb1-10 yb2+10])
% xlim([zb1 zb2])
% ylim([xb1 xb2])
% zticks([round(m2.ystart) round(m2.yend)])
% xticks([round(m2.zstart) ,0,round(m2.zend)])
% yticks([round(m2.xstart) ,0,round(m2.xend)])
zlim([yb1-10 yb2+10])
xlim([zb1-1 zb2+1])
ylim([xb1-1 xb2+1])

 zticks([round(m1.ystart) round(m1.yend)])
 xticks([zb1,0,zb2 ])
 yticks([xb1,0,xb2 ])

clim([-0.1 0.1])
colormap saffrongreen
 if j==5
 c=colorbar;
 ylabel(c,"$\langle\omega_x^+\rangle_-$",'interpreter','latex','FontSize',11)
c.Location='eastoutside';
 end
%set(gca,'FontSize',11)

xlabel('$z^+$','interpreter','latex','FontSize',11)
ylabel('$x^+$','interpreter','latex','FontSize',11)
if j==1 
    zlabel('$y^+$','interpreter','latex','FontSize',11)
end
%title('(a)')
title( labels(j+5) )
%%
%subplot(3,5,10+j)
% % nexttile(j+10)
% % hold on
% % plot(mpr.nldprofiles(:,j)./(-ut^2),yp,'-k','LineWidth',1.5)
% % plot(mpr.nluprofiles(:,j)./(-ut^2),yp,'-.','LineWidth',1.5,'Color',G)
% % hold off
% % xline(0)
% % yline(yp(jcond))
% % xlim([-0.2 0.4])
% % ylim([round(m1.ystart) round(m1.yend)])
% % yticks([ round(m1.ystart), round(yp(jcond)), round(m1.yend)] )
% % title(labels(j+10))
% % if j==1
% %     ylabel('$y^+$','interpreter','latex','FontSize',11)
% %     legend('v_+','v_-','Location','southeast')
% %     legend boxoff
% % end
% % xlabel('$\langle v\omega_z-w\omega_y\rangle/(-u_{\tau}^2/H)$','interpreter','latex')
end
%saveas(fd,'eddy_iso_flux_compile.fig')
saveas(fd,'eddy_iso_compile.fig')
