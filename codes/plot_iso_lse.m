clear 
close all
jcond=130;
fn=sprintf('../data/lse_eddyset_j_%03d.mat',jcond)
x1=100;
y1=100;
width=1150;
height=300;
ut=0.0499;
dnu=1.0006e-3;
%%
load(fn)
z=Z./dnu;
y=Y./dnu;
x=X./dnu;
%ys=y(1,1,jcond-110);
ys=Y(1,1,jcond)

z=Z./ys;
y=Y./ys;
x=X./ys;
ltdm=min(ld(:,:,jcond),[],'all');
ltum=min(lu(:,:,jcond),[],'all');
%alpha=0.05;
% ltd=alpha*ltdm
% ltu=alpha*ltum
alpha=-0.0004;
ltd=alpha/(ys)^2;
ltu=alpha/(ys)^2
idx=52;
%load('eddyset_d_j_156.mat')
xl=2;
yl=5;
zl=2;


%%
% [startZ,startX,startY]=meshgrid(0.3./dnu,[-0.2:0.025:0.2]./dnu,ys./dnu);
% vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);


fd=figure;
fd.Position=[x1 y1 width height];
subplot(1,2,1)
hold on;
%
isosurf=isosurface(z,x,y,ld,ltd);
interpColors = interp3(z, x, y, dnu*ys.*oxd./ut,...wd.*0, ...
    isosurf.vertices(:, 1), ...
    isosurf.vertices(:, 2), ...
    isosurf.vertices(:, 3));
p = patch(isosurf);
p.FaceColor = 'interp';        % Interpolated color
p.EdgeColor = 'none';          % Remove edges
p.FaceVertexCData = interpColors;    % Assign interpolated colors
p.FaceAlpha = 0.7;  
camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
%%
scatter3(0,0,y(2,2,jcond),100,'green','filled')
%%
%isonormals(z,x,y,ld,-0.005);
%l=streamline(vertsv);
%set(l, 'Color', 'k'); 
%set(l,'LineWidth',2);

%%
colormap redblue
axis equal
view(45,35);
             % Smooth shading
grid on;
shading flat

ylim([-yl yl])
xlim([-xl xl])
%zlim([0 200])
zlim([0 zl])
c=colorbar;
%clim([-20 20])
%ylabel(c,"u'")
set(gca,'FontSize',12)
ylabel(c,"$\bar{y}\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z/\bar{y}$','interpreter','latex','FontSize',14)
ylabel('$x/\bar{y}$','interpreter','latex','FontSize',14)
zlabel('$y/\bar{y}$','interpreter','latex','FontSize',14)
title('(a)')

%%
%[startZ,startX,startY]=meshgrid(0.3./dnu,[-0.2:0.025:0.2]./dnu,ys./dnu);
%vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);

% [startZp,startXp,startYp]=meshgrid(-0.15:0.2:0.5,0.8,0.05:0.05:1);
% vertsp = stream3(z,x,y,wds,uds,vds,startZp,startXp,startYp);


%%
subplot(1,2,2)
hold on;
%
%ax1=axes;
isosurf=isosurface(z,x,y,lu,ltu);
interpColors = interp3(z, x, y,dnu*ys.*oxu./ut,...wu.*0, ...
    isosurf.vertices(:, 1), ...
    isosurf.vertices(:, 2), ...
    isosurf.vertices(:, 3));
p = patch(isosurf);
p.FaceColor = 'interp';        % Interpolated color
p.EdgeColor = 'none';          % Remove edges
p.FaceVertexCData = interpColors;    % Assign interpolated colors
p.FaceAlpha = 0.7;  
camlight;                      % Add lighting
%light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
lighting gouraud; 
lightangle(-45,90)
%%
scatter3(0,0,y(2,2,jcond),100,'green','filled')

%%
% l=streamline(vertsv2);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',2);


%%
colormap redblue

axis equal
%view(135,35);
view(45,35)             % Smooth shading
grid on;
shading flat
ylim([-yl yl])
xlim([-xl xl])
%zlim([0 200])
zlim([0 zl])
%clim([-0.3 0.3])
%clim([-0.03 0.03])
%clim([-0.2 0.2])

c=colorbar;
%clim([-10 10])
%ylabel(c,"u'")
set(gca,'FontSize',12)
ylabel(c,"$\bar{y}\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z/\bar{y}$','interpreter','latex','FontSize',14)
ylabel('$x/\bar{y}$','interpreter','latex','FontSize',14)
zlabel('$y/\bar{y}$','interpreter','latex','FontSize',14)
title('(b)')
%%
    fdne=sprintf('eddy_scaled_j_%03d.fig',jcond)
   saveas(fd,fdne)

%% calculate bounding box