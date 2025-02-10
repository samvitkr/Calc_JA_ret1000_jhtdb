clear 
close all
jcond=71;
fn=sprintf('../data/lse_eddyset_j_%03d.mat',jcond)
x1=100;
y1=100;
width=1100;
height=350;
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
alpha=0.05;
ltd=alpha*ltdm;
ltu=alpha*ltum;
idx=52;
%load('eddyset_d_j_156.mat')
%%
[startZ,startX,startY]=meshgrid(0.15./dnu,[-0.1:0.05:0.1]./dnu,ys./dnu);
vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);

% [startZ,startX,startY]=meshgrid(0.15,0,ys-0.2:ys:ys+0.2);
% vertsv2 = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);

%streamline(vertsv)
% Plot the vector field and streamlines
fd=figure;
fd.Position=[x1 y1 width height];
subplot(1,2,1)
hold on;
%
isosurf=isosurface(z,x,y,ld,ltd);
interpColors = interp3(z, x, y, oxd,...wd.*0, ...
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
% l2=streamline(vertsv2);
% set(l2, 'Color', 'k'); 
% set(l2,'LineWidth',2);
%%
xs=0.8;
oxds=oxd.*0;
uds=ud.*0;
vds=ud.*0;
wds=ud.*0;
for i =1:101
   % oxds(i,:,:)=oxd(50,1:end,1:end);
   oxds(i,:,:)=ud(idx,1:end,1:end);
   uds(i,:,:)=ud(idx,1:end,1:end).*0;
   vds(i,:,:)=vd(idx,1:end,1:end);
   wds(i,:,:)=wd(idx,1:end,1:end);

end
% % % slicePlane = slice(z,x,y, oxds, [], xs, []); % Slice at x=0
% % % %slicePlane2 = slice(z,x,y, oxd, -0.5, [], []); % Slice at x=0
% % % 
% % % set(slicePlane, 'EdgeColor', 'none');  % Remove grid lines on the slice
% % % set(slicePlane, 'FaceColor', 'interp'); % Smooth interpolation
% % % set(slicePlane, 'FaceLighting', 'none'); % Removes shading for brighter appearance
% % % 
% % % %set(slicePlane2, 'EdgeColor', 'none');  % Remove grid lines on the slice
% % % %set(slicePlane2, 'FaceColor', 'interp'); % Smooth interpolation
% % % %set(slicePlane2, 'FaceLighting', 'none'); % Removes shading for brighter appearance
% % % quiverX = squeeze(x(idx,1:5:end, 1:5:end)).*0+xs;
% % % quiverY = squeeze(y(idx,1:5:end, 1:5:end));
% % % quiverZ = squeeze(z(idx,1:5:end, 1:5:end));
% % % quiverU = squeeze(ud(idx, 1:5:end, 1:5:end)).*0;
% % % quiverV = squeeze(vd(idx, 1:5:end, 1:5:end));
% % % quiverW = squeeze(wd(idx,1:5:end, 1:5:end));
% % % 
% % % % Reshape to 2D for quiver3
% % % % quiverX = squeeze(quiverX);
% % % % quiverY = squeeze(quiverY);
% % % % quiverZ = squeeze(quiverZ);
% % % % quiverU = squeeze(quiverU);
% % % % quiverV = squeeze(quiverV);
% % % % quiverW = squeeze(quiverW);
% % % scaleFactor=4;
% % % quiver3(quiverZ, quiverX, quiverY, quiverW, quiverU, quiverV,scaleFactor, 'k','LineWidth',1.5); % Black arrows
%%
colormap redblue
axis equal
view(45,25);
             % Smooth shading
grid on;
shading flat
ylim([-0.3 0.15]./dnu)
xlim([-0.12 0.12]./dnu)
%yticks([-1:0.2:0.8])
zlim([0 0.2]./dnu)
%clim([-0.3 0.3001])
%clim([-0.05 0.05])
clim([-0.2 0.2])
xlabel('z^+')
ylabel('x^+')
zlabel('y^+')
c=colorbar;
%ylabel(c,"u'")
ylabel(c,"\omega_x")
set(gca,'FontSize',11)
%%
[startZ,startX,startY]=meshgrid(0.15./dnu,[-0.1:0.05:0.1]./dnu,ys./dnu);
vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);

[startZp,startXp,startYp]=meshgrid(-0.15:0.2:0.5,0.8,0.05:0.05:1);
vertsp = stream3(z,x,y,wds,uds,vds,startZp,startXp,startYp);

%streamline(vertsv)
% Plot the vector field and streamlines
%fu=figure;
%%
subplot(1,2,2)
hold on;
%
%ax1=axes;
isosurf=isosurface(z,x,y,lu,ltu);
interpColors = interp3(z, x, y,oxu,...wu.*0, ...
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
l=streamline(vertsv2);
set(l, 'Color', 'k'); 
set(l,'LineWidth',2);

% l=streamline(vertsp,[]);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1);
%%
% % % %ax2=axes;
% % % oxus=oxd.*0;
% % % for i =1:101
% % %     %oxus(i,:,:)=oxu(50,1:end,1:end);
% % %     oxus(i,:,:)=uu(idx,1:end,1:end);
% % % end
% % % slicePlane = slice(z,x,y, oxus, [], 0.8, []); % Slice at x=0
% % % %slicePlane2 = slice(z,x,y, oxu, -0.5, [], []); % Slice at x=0
% % % 
% % % set(slicePlane, 'EdgeColor', 'none');  % Remove grid lines on the slice
% % % set(slicePlane, 'FaceColor', 'interp'); % Smooth interpolation
% % % set(slicePlane, 'FaceLighting', 'none'); % Removes shading for brighter appearance
% % % 
% % % %set(slicePlane2, 'EdgeColor', 'none');  % Remove grid lines on the slice
% % % %set(slicePlane2, 'FaceColor', 'interp'); % Smooth interpolation
% % % %set(slicePlane2, 'FaceLighting', 'none'); % Removes shading for brighter appearance
% % % 
% % % quiverX = squeeze(x(idx,1:5:end, 1:5:end)).*0+xs;
% % % quiverY = squeeze(y(idx,1:5:end, 1:5:end));
% % % quiverZ = squeeze(z(idx,1:5:end, 1:5:end));
% % % quiverU = squeeze(uu(idx, 1:5:end, 1:5:end)).*0;
% % % quiverV = squeeze(vu(idx, 1:5:end, 1:5:end));
% % % quiverW = squeeze(wu(idx,1:5:end, 1:5:end));
% % % 
% % % 
% % %  scaleFactor=4;
% % %  quiver3(quiverZ, quiverX, quiverY, quiverW, quiverU, quiverV,scaleFactor, 'k','LineWidth',1.5);
%%
colormap redblue
%colormap(ax1, jet);  % Use 'jet' colormap for the slice
%colormap(ax2, redblue);  % Use 'jet' colormap for the slice
axis equal
view(135,25);
             % Smooth shading
grid on;
shading flat
ylim([-0.3 0.15]./dnu)
xlim([-0.12 0.12]./dnu)
%yticks([-1:0.2:0.8])
zlim([0 0.2]./dnu)
%clim([-0.3 0.3])
%clim([-0.03 0.03])
clim([-0.2 0.2])
c=colorbar;
%ylabel(c,"u'")
ylabel(c,'\omega_x')
xlabel('z^+')
ylabel('x^+')
zlabel('y^+')
set(gca,'FontSize',11)

%%
%  fdn=sprintf('eddy_j_%03d.fig',jcond)
% % fun=sprintf('eddy_u_j_%03d.fig',jcond)
%  saveas(fd,fdn)
   fdne=sprintf('eddy_j_%03d.fig',jcond)
% % fun=sprintf('eddy_u_j_%03d.fig',jcond)
  saveas(fd,fdne)
% % saveas(fu,fun)
% % % %saveas(fd,'eddy_d_j_156.fig')
% % % %saveas(fu,'eddy_u_j_156.fig')
%% calculate bounding box
