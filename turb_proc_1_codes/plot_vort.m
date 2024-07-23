Nxtot=2048;
Nytot=512;
Nztot=1536;
Lx=  8*pi;
Lz = 3*pi;
xp = [0:Nxtot-1]*Lx/(Nxtot);
zp=  [0:1:Nztot-1]*Lz/(Nztot);
mb=matfile('bsplinedata.mat');
ut = 0.0499;
dnu=1.0006e-3;
yp=(1+mb.yv)./dnu;
mcv=matfile('Transfer_018.mat');
ml=matfile('lambda_018.mat');
nx=200;
ny=120;
nz=200;
[x,z,y] = meshgrid(zp(1:nz),xp(1:nx),yp(1:ny));
l=permute(ml.lambda2(1:ny,1:nx,1:nz),[2 3 1]);
q=permute(ml.Q(1:ny,1:nx,1:nz),[2 3 1]);
subplot(2,1,1)
pcolor(squeeze(ml.lambda2(ny,1:nx,1:nz)))
shading interp
%pbaspect([xp(nx),zp(nz) 1])
caxis([-100 100])
colormap jet
title('\lambda_2')


subplot(2,1,2)
pcolor(squeeze(ml.Q(ny,1:nx,1:nz)))
shading interp
%pbaspect([xp(nx),zp(nz) 1])
caxis([-100 100])
colormap jet
title('Q')

%[faces,verts,colors] = isosurface(z,x,y,q.*dnu/ut,-2,q.*dnu/ut);

%p=patch('Vertices', verts, 'Faces', faces, ...
%      'FaceVertexCData', colors, ...
%      'FaceColor','interp', ...
%      'edgecolor', 'interp');

%p.FaceColor = 'interp';
%p.EdgeColor = 'none';

%camlight
%lighting gouraud
%xlabel('x/h')
%ylabel('z/h')
%zlabel('y^+')
%view(-45,45)
%pbaspect([xp(nx),zp(nz) 1])
%caxis([-100 100])
%colormap jet
