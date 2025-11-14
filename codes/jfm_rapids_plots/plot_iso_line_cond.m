clear 
close all

x1=100;
y1=100;
width=720;
height=350;
ut=0.0499;
dnu=1.0006e-3;
% mU=matfile('../data/JHTDB_RET1000.mat');
% Um=mU.JHTDB_RET1000(:,2);
% Um=Um./Um(end);

clear mU

load('../../data/bsplinedata.mat')
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
jcond=105;
yc=yv(jcond)+1;
ycp=yp(jcond);
% fvgp=sprintf("../../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
% fvgn=sprintf("../../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);

fvgp=sprintf("../../data/conditionalp_jcond_%03d.mat",jcond);
fvgn=sprintf("../../data/conditionaln_jcond_%03d.mat",jcond);

figname=sprintf("isosurf_v_dudy_jcond_%03d.fig",jcond);

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
val=5*l1;
[X,Z,Y]=(meshgrid(xp,zp,yp));

x=permute(X,[2 1 3]);
z=permute(Z,[2 1 3]);
y=permute(Y,[2 1 3]);


oxd=permute(m1.dwdy-m1.dvdz,[2 1 3]);
oyd=permute(m1.dudz-m1.dwdx,[2 1 3]);
ozd=permute(m1.dvdx-m1.dudy,[2 1 3]);
ld=permute(m1.lambda2,[2 1 3]);
oxu=permute(m2.dwdy-m2.dvdz,[2 1 3]);
oyu=permute(m2.dudz-m2.dwdx,[2 1 3]);
ozu=permute(m2.dvdx-m2.dudy,[2 1 3]);
lu=permute(m2.lambda2,[2 1 3]);


[startZ,startX,startY]=meshgrid(0,[-150:25:150],ycp);
   % [startZ,startX,startY]=meshgrid(0,0,ycp);

% % z_stop=-m1.zend-10;

vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
vertsvn = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);
%[startZ,startX,startY]=meshgrid(0.1./dnu,[-0.12:0.012:0.12]./dnu,ys);
% [startZ,startX,startY]=meshgrid(m1.zend+10,[m2.xstart-20:10:m2.xend+80],1.2*ys);
vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
vertsv2n = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);


xl1=-150;
xl2=150;
yl1=-200;
yl2=200;
zl1=50;
zl2=300;
cl1=-0.05;
cl2=0.05;


f=figure('OuterPosition',[x1,y1,width,height]);


t=tiledlayout(1,2);
t.TileSpacing="tight";
t.Padding="tight";
nexttile
%subplot(1,2,1)
verts = vertsv;
hold on

for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end

    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);

    % Interpolate w along streamline
    ws = interp3(z,x,y,ozd,xs,ys,zs);  

    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',0.5);
end
sd = [ys,zs];
verts = vertsvn;
for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end

    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);

    % Interpolate w along streamline
    ws = interp3(z,x,y,ozd,xs,ys,zs);  

    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',0.5);
end

 isosurf=isosurface(z,x,y,ld,val);
  interpColors = interp3(z, x, y, ozd./(ut./dnu),...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
      isosurf.vertices(:, 1), ...
      isosurf.vertices(:, 2), ...
      isosurf.vertices(:, 3));
  p = patch(isosurf);
p.FaceColor = 'interp';        % Interpolated color
p.EdgeColor = 'none';          % Remove edges
p.FaceVertexCData = interpColors;    % Assign interpolated colors
p.FaceAlpha = 0.8; 
p.FaceLighting="gouraud";
p.SpecularStrength=0.2;
p.AmbientStrength=0.6;
scatter3(0,0,ycp,50,'green','filled')
hold off

colormap redblue
axis equal
%view(0,65);
view(45,45)
zlim([zl1 zl2])
xlim([xl1 xl2])
ylim([yl1 yl2])
clim([cl1 cl2])
title('(a)')
xlabel('z^+')
ylabel('x^+')
zlabel('y^+')
%%
%subplot(1,2,2)
% title('$(a) outflow \ and \ du/dy>0$','Interpreter','latex')
nexttile
verts = vertsv2;
hold on
for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end

    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);

    % Interpolate w along streamline
    ws = interp3(z,x,y,ozu,xs,ys,zs);  

    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',0.5);
end
sd = [ys,zs];
verts = vertsv2n;
for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end

    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);

    % Interpolate w along streamline
    ws = interp3(z,x,y,ozu,xs,ys,zs);  

    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',0.5);
end
 isosurf=isosurface(z,x,y,lu,val);
  interpColors = interp3(z, x, y, ozu./(ut/dnu),...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
      isosurf.vertices(:, 1), ...
      isosurf.vertices(:, 2), ...
      isosurf.vertices(:, 3));
  p = patch(isosurf);
p.FaceColor = 'interp';        % Interpolated color
p.EdgeColor = 'none';          % Remove edges
p.FaceVertexCData = interpColors;    % Assign interpolated colors
p.FaceAlpha = 0.8; 
p.FaceLighting="gouraud";
p.SpecularStrength=0.2;
p.AmbientStrength=0.6;
scatter3(0,0,ycp,50,'green','filled')
hold off
colormap redblue
axis equal
%view(0,65);
view(45,45)
zlim([zl1 zl2])
xlim([xl1 xl2])
ylim([yl1 yl2])
clim([cl1 cl2])
c=colorbar;
ylabel(c,'\omega_z^+')
title('(b)')
xlabel('z^+')
ylabel('x^+')
%%
% saveas(f,'iso_lines_cond_j_105.fig')