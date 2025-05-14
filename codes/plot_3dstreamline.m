clear 
close all
jcond=71;
fn=sprintf('../data/lse_eddyset_j_%03d.mat',jcond)
fn2=sprintf('./eddy_isosurf_stream_j_%03d.mat',jcond)
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
ltdm=min(ld(:,:,jcond),[],'all');
ltum=min(lu(:,:,jcond),[],'all');
% alpha=0.05;
% alpha=0.02;
% ltd=alpha*ltdm
% ltu=alpha*ltum
% ltd=-0.02;
% ltu=-0.02;
alpha=-0.0008;
ltd=alpha/(ys)^2;
ltu=alpha/(ys)^2
idx=52;
%load('eddyset_d_j_156.mat')
%%
[startZ,startX,startY]=meshgrid(0.3./dnu,[-0.2:0.025:0.2]./dnu,ys./dnu);
vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);


fd=figure;
fd.Position=[x1 y1 width height];
subplot(1,2,1)
hold on;
%
isosurf=isosurface(z,x,y,ld,ltd);
interpColors = interp3(z, x, y, oxd./ut,...wd.*0, ...
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
l=streamline(vertsv);
set(l, 'Color', 'k'); 
set(l,'LineWidth',2);

%%
colormap redblue
axis equal
view(45,35);
             % Smooth shading
grid on;
shading flat

ylim([-200 150])
xlim([-100 100])
%zlim([0 200])
zlim([0 160])
c=colorbar;
%clim([-20 20])
%ylabel(c,"u'")
set(gca,'FontSize',12)
ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(a)')

vertices1=single(isosurf.vertices);
indices1=uint32(isosurf.faces-1);
faceColor1=single(interpColors);
maxLength = max(cellfun(@(s) size(s, 1), vertsv));
numStreams = numel(vertsv);
streamlineArray1 = NaN(maxLength,3,numStreams); % 3 for x, y, z coordinates
for i = 1:numStreams
    streamlineLength = size(vertsv{i}, 1);
    streamlineArray1(1:streamlineLength,:,i) = vertsv{i}; % Fill in streamlines
end
%%
%[startZ,startX,startY]=meshgrid(0.3./dnu,[-0.2:0.025:0.2]./dnu,ys./dnu);
vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
% [startZp,startXp,startYp]=meshgrid(-0.15:0.2:0.5,0.8,0.05:0.05:1);
% vertsp = stream3(z,x,y,wds,uds,vds,startZp,startXp,startYp);
%%
subplot(1,2,2)
hold on;
%
%ax1=axes;
isosurf=isosurface(z,x,y,lu,ltu);
interpColors = interp3(z, x, y,oxu./ut,...wu.*0, ...
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
l=streamline(vertsv2);
set(l, 'Color', 'k'); 
set(l,'LineWidth',2);
%%
colormap redblue
axis equal
%view(135,35);
view(45,35)             % Smooth shading
grid on;
shading flat
ylim([-200 200])
xlim([-100 100])
zlim([0 160])
%clim([-0.3 0.3])
%clim([-0.03 0.03])
%clim([-0.2 0.2])
c=colorbar;
%clim([-10 10])
%ylabel(c,"u'")
set(gca,'FontSize',12)
ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(b)')
vertices2=single(isosurf.vertices);
indices2=uint32(isosurf.faces-1);
faceColor2=single(interpColors);

maxLength = max(cellfun(@(s) size(s, 1), vertsv2));
numStreams = numel(vertsv2);
streamlineArray2 = NaN(maxLength,3,numStreams); % 3 for x, y, z coordinates
for i = 1:numStreams
    streamlineLength = size(vertsv2{i}, 1);
    streamlineArray2(1:streamlineLength,:,i) = vertsv2{i}; % Fill in streamlines
end


%%
%    fdne=sprintf('eddy_j_%03d.fig',jcond)
%    saveas(fd,fdne)
% %%
% save(fn2,'vertices1','indices1','faceColor1','streamlineArray1','vertices2','indices2','faceColor2','streamlineArray2')


%% calculate bounding box
