clear 
close all
jcond=71;
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
ys=Y(1,1,jcond);
ltdm=min(ld(:,:,jcond),[],'all');
ltum=min(lu(:,:,jcond),[],'all');
alpha=0.05;
alpha=0.02;
ltd=alpha*ltdm;
ltu=alpha*ltum;
ltd=-0.02;
ltu=-0.02;
idx=52;
%%
[startZ,startX,startY]=meshgrid(0.3./dnu,[-0.2:0.025:0.2]./dnu,ys./dnu);
vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);

fd=figure;
fd.Position=[x1 y1 width height];
subplot(1,2,1)
hold on;
%
isosurf1 = isosurface(z,x,y,ld,ltd);
interpColors1 = interp3(z, x, y, oxd./ut, ...
    isosurf1.vertices(:, 1), ...
    isosurf1.vertices(:, 2), ...
    isosurf1.vertices(:, 3));

p1 = patch(isosurf1);
p1.FaceColor = 'interp';        % Interpolated color
p1.EdgeColor = 'none';          % Remove edges
p1.FaceVertexCData = interpColors1;    
p1.FaceAlpha = 0.7;  
camlight;  
lighting gouraud; 
lightangle(-45,90);

%%
scatter3(0,0,y(2,2,jcond),100,'green','filled')
%%
l1 = streamline(vertsv);
set(l1, 'Color', 'k'); 
set(l1,'LineWidth',2);

%%
colormap redblue
axis equal
view(45,35);
grid on;
shading flat

ylim([-600 200])
xlim([-150 150])
zlim([0 250])
c=colorbar;
set(gca,'FontSize',12)
ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(a)')

%% Second Isosurface & Streamlines
vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);

subplot(1,2,2)
hold on;
%
isosurf2 = isosurface(z,x,y,lu,ltu);
interpColors2 = interp3(z, x, y,oxu./ut, ...
    isosurf2.vertices(:, 1), ...
    isosurf2.vertices(:, 2), ...
    isosurf2.vertices(:, 3));

p2 = patch(isosurf2);
p2.FaceColor = 'interp';
p2.EdgeColor = 'none';
p2.FaceVertexCData = interpColors2;
p2.FaceAlpha = 0.7;  
camlight;  
lighting gouraud; 
lightangle(-45,90);

%%
scatter3(0,0,y(2,2,jcond),100,'green','filled')

%%
l2 = streamline(vertsv2);
set(l2, 'Color', 'k'); 
set(l2,'LineWidth',2);

%%
colormap redblue
axis equal
view(45,35);
grid on;
shading flat
ylim([-300 300])
xlim([-150 150])
zlim([0 250])
c=colorbar;
set(gca,'FontSize',12)
ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(b)')

%%
fdne=sprintf('eddy_j_%03d.fig',jcond);
%saveas(fd,fdne);

%% Save Data for Export
% Extract streamline data
streamlineData1 = {};
for i = 1:length(vertsv)
    streamlineData1{i} = [vertsv{i}(:,1), vertsv{i}(:,2), vertsv{i}(:,3)];
end

streamlineData2 = {};
for i = 1:length(vertsv2)
    streamlineData2{i} = [vertsv2{i}(:,1), vertsv2{i}(:,2), vertsv2{i}(:,3)];
end

% Structure to save
exportData.isosurface1.vertices = isosurf1.vertices;
exportData.isosurface1.faces = isosurf1.faces;
exportData.isosurface1.faceColorData = interpColors1;

exportData.isosurface2.vertices = isosurf2.vertices;
exportData.isosurface2.faces = isosurf2.faces;
exportData.isosurface2.faceColorData = interpColors2;

exportData.streamlines1 = streamlineData1;
exportData.streamlines2 = streamlineData2;

save('exported_data.mat', 'exportData');
disp('Isosurface and streamline data saved.');

