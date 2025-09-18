clear 
close all

x1=100;
y1=100;
width=1150;
height=300;
ut=0.0499;
dnu=1.0006e-3;

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

%jcond=130;
jcond=105
% jcond=71;
%jcond=54;
% jcond=47;
yc=yv(jcond)+1;

ys=yp(jcond);
fvgp=sprintf("../data/conditionalp_jcond_omegaz_%03d.mat",jcond);
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

% xb1=min([m1.xstart m2.xstart]);
% yb1=min([m1.ystart m2.ystart]);
% zb1=min([m1.zstart m2.zstart]);
% 
% xb2=max([m1.xend m2.xend]);
% yb2=max([m1.yend m2.yend]);
% zb2=max([m1.zend m2.zend]);


%%
% % % % % load(fn)
% % % % % z=Z./dnu;
% % % % % y=Y./dnu;
% % % % % x=X./dnu;
% % % % % %ys=y(1,1,jcond-110);
% % % % % ys=Y(1,1,jcond)
% % % % % ltdm=min(ld(:,:,jcond),[],'all');
% % % % % ltum=min(lu(:,:,jcond),[],'all');
% % % % % % alpha=0.05;
% % % % % % alpha=0.02;
% % % % % % ltd=alpha*ltdm
% % % % % % ltu=alpha*ltum
% ltd=-0.8;
% ltu=-0.8;
ltd=val;
ltu=val;
% % % % % % ltu=-0.02;
% % % % % alpha=-0.0008;
% % % % % ltd=alpha/(ys)^2;
% % % % % ltu=alpha/(ys)^2
% % % % % idx=52;
%load('eddyset_d_j_156.mat')
%%
%[startZ,startX,startY]=meshgrid(0.1./dnu,[-0.12:0.012:0.12]./dnu,ys/1.5);
% % % [startZ,startX,startY]=meshgrid(m1.zend+10,[m1.xstart-50:10:m1.xend+50],ys/1.5);
% % % z_stop=-m1.zend-10;
% % % 
% % % vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
% % % 
% % % 
% % % for i = 1:length(vertsv)
% % %     thisStream = vertsv{i};
% % %     if isempty(thisStream), continue; end
% % %     % Find index where z reaches or exceeds z_stop
% % %     idx = find(thisStream(:,1) <= z_stop, 1);
% % %     if ~isempty(idx)
% % %         vertsv{i} = thisStream(1:idx, :);  % truncate streamline
% % %     end
% % % end
% % % 
% % % %[startZ,startX,startY]=meshgrid(0.1./dnu,[-0.12:0.012:0.12]./dnu,ys);
% % % [startZ,startX,startY]=meshgrid(m1.zend+10,[m2.xstart-20:10:m2.xend+80],1.2*ys);
% % % 
% % % vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
% % % 
% % % for i = 1:length(vertsv2)
% % %     thisStream = vertsv2{i};
% % %     if isempty(thisStream), continue; end
% % %     % Find index where z reaches or exceeds z_stop
% % %     idx = find(thisStream(:,1) <= z_stop, 1);
% % %     if ~isempty(idx)
% % %         vertsv2{i} = thisStream(1:idx, :);  % truncate streamline
% % %     end
% % % end
% % % 
% % % fd=figure;
% % % fd.Position=[x1 y1 width height];
% % % %%

subplot(1,2,1)
% % % if (jcond==105)
% % %     ld(:,:,1:51)=0;
% % %     lu(:,:,1:51)=0;
% % % end
% % % if (jcond==130)
% % %     ld(:,:,1:81)=0;
% % %     lu(:,:,1:81)=0;
% % % end
hold on;
%
isosurf=isosurface(z,x,y,ld,ltd);
interpColors = interp3(z, x, y, -nld./ut^2,...dnu*oxd./ut,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
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
scatter3(0,0,ys,100,[0.5 0 0.5],'filled')
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1.5);
% % % % colormap redblue
% % % % axis equal
% % % % view(45,45);


% ylim([-150 150])
% xlim([-100 100])
% %zlim([0 200])
% zlim([0 150])
% yticks([-150:50:150])
% zlim([50 400])
% zlim([m1.ystart-10 m1.yend+10])
% xlim([m1.zstart-10 m1.zend+10])
% ylim([m1.xstart-20 m1.xend+20])

% % % zlim([yb1-10 yb2+10])
% % % xlim([zb1-10 zb2+10])
% % % ylim([xb1-20 xb2+20])
% % % 
% % % zticks([round(m1.ystart) round(m1.yend)])
% % % xticks([round(m1.zstart) ,0,round(m1.zend)])
% % % yticks([round(m1.xstart) ,0,round(m1.xend)])
%c1=colorbar;
clim([-cld cld])
%c1.Ticks=[-cld, 0, cld];

%clim([-1 1]*1e-1)
%clim([-150 150])

%ylabel(c,"u'")
set(gca,'FontSize',11)
%ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)

xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
title('(a)')
%%
% % % vertices1=single(isosurf.vertices);
% % % indices1=uint32(isosurf.faces-1);
% % % faceColor1=single(interpColors);
% % % maxLength = max(cellfun(@(s) size(s, 1), vertsv));
% % % numStreams = numel(vertsv);
% % % streamlineArray1 = NaN(maxLength,3,numStreams); % 3 for x, y, z coordinates
% % % for i = 1:numStreams
% % %     streamlineLength = size(vertsv{i}, 1);
% % %     streamlineArray1(1:streamlineLength,:,i) = vertsv{i}; % Fill in streamlines
% % % end

%%
subplot(1,2,2)
hold on;
%
isosurf=isosurface(z,x,y,lu,ltu);
interpColors = interp3(z, x, y, -nlu./ut^2,...dnu*oxu./ut,...wd.*0, ...interpColors = interp3(z, x, y, -nlu./ut^2,...wd.*0, ...
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
scatter3(0,0,ys,100,[0.5 0 0.5],'filled')
% l=streamline(vertsv2);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1.5);
colormap redblue
axis equal
%view(45,45);
view(135,45)
% ylim([-150 150])
% xlim([-100 100])
% %zlim([0 200])
% zlim([0 150])
% yticks([-150:50:150])

c=colorbar;
%clim([-1 1]*1e-1)
%clim([-150 150])
%ylabel(c,"u'")
set(gca,'FontSize',11)
%ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
%ylabel(c,"$(v\omega_z-w\omega_y)/(-u_{\tau}^2/H)$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)

% zlim([50 400])

% zlim([m2.ystart-10 m2.yend+10])
% xlim([m2.zstart-10 m2.zend+10])
% ylim([m2.xstart-10 m2.xend+40])

% zlim([m1.ystart-10 m1.yend+10])
% xlim([m1.zstart-10 m1.zend+10])
% ylim([m1.xstart-20 m1.xend+20])
% 
% zticks([round(m2.ystart) round(m2.yend)])
% xticks([round(m2.zstart) ,0,round(m2.zend)])
% yticks([round(m2.xstart) ,0,round(m2.xend)])

% zlim([yb1-10 yb2+10])
% xlim([zb1-10 zb2+10])
% ylim([xb1-20 xb2+20])
% 
% zticks([round(m2.ystart) round(m2.yend)])
% xticks([round(m2.zstart) ,0,round(m2.zend)])
% yticks([round(m2.xstart) ,0,round(m2.xend)])

% % % c=colorbar;
% % % clim([-cld cld])
% % %  c.Ticks=[-cld, 0, cld];
% % %  ylabel(c,"$\omega_x^+$",'interpreter','latex','FontSize',14)
% % % 
% % % %zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % title('(b)')
% % % %%
% % % vertices2=single(isosurf.vertices);
% % % indices2=uint32(isosurf.faces-1);
% % % faceColor2=single(interpColors);
% % % 
% % % maxLength = max(cellfun(@(s) size(s, 1), vertsv2));
% % % numStreams = numel(vertsv2);
% % % streamlineArray2 = NaN(maxLength,3,numStreams); % 3 for x, y, z coordinates
% % % for i = 1:numStreams
% % %     streamlineLength = size(vertsv2{i}, 1);
% % %     streamlineArray2(1:streamlineLength,:,i) = vertsv2{i}; % Fill in streamlines
% % % end
% % % %%
% % % fn2=sprintf('./eddy_isosurf_stream_j_%03d.mat',jcond)
% % % mfn2=matfile(fn2,'Writable',true)
% % % mfn2.faceflux1=faceColor1;
% % % mfn2.faceflux2=faceColor2;
% % % 
% % % %  save(fn2,'vertices1','indices1','faceColor1','streamlineArray1','vertices2','indices2','faceColor2','streamlineArray2','cld','clu')
% % % 
% % % fdne=sprintf('eddy_j_cond_%03d.fig',jcond)
% % % %    fdne=sprintf('eddy_j_cond_nl_%03d.fig',jcond)   
% % % %saveas(fd,fdne)
