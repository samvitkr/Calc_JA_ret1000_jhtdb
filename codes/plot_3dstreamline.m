clear 
close all

x1=100;
y1=100;
width=1150;
height=300;
ut=0.0499;
dnu=1.0006e-3;
mU=matfile('../data/JHTDB_RET1000.mat');
Um=mU.JHTDB_RET1000(:,2);
Um=Um./Um(end);

clear mU

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

% jcond=130;
 jcond=71;
% jcond=105
% jcond=71;
 % jcond=54;
  % jcond=47;
yc=yv(jcond)+1;

ycp=yp(jcond);
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
ry=y-ycp;
r=sqrt(x.^2+z.^2+ry.^2);
oxd=permute(m1.dwdy-m1.dvdz,[2 1 3]);
oyd=permute(m1.dudz-m1.dwdx,[2 1 3]);
ozd=permute(m1.dvdx-m1.dudy,[2 1 3]);
omd=sqrt(oxd.^2+oyd.^2+ozd.^2);



nld=permute(m1.voz-m1.woy,[2 1 3]);
ld=permute(m1.lambda2,[2 1 3])./l1;
vozd=permute(m1.voz,[2 1 3]);
woyd=permute(m1.woy,[2 1 3]);
ud=permute(m1.u,[2 1 3]);
vd=permute(m1.v,[2 1 3]);
wd=permute(m1.w,[2 1 3]);


vdozd=vd.*ozd;
wdoyd=wd.*oyd;

oxu=permute(m2.dwdy-m2.dvdz,[2 1 3]);
oyu=permute(m2.dudz-m2.dwdx,[2 1 3]);
ozu=permute(m2.dvdx-m2.dudy,[2 1 3]);
omu=sqrt(oxu.^2+oyu.^2+ozu.^2);

lu=permute(m2.lambda2,[2 1 3])./l2;
nlu=permute(m2.voz-m2.woy,[2 1 3]);
vozu=permute(m2.voz,[2 1 3]);
woyu=permute(m2.woy,[2 1 3]);

uu=permute(m2.u,[2 1 3]);
vu=permute(m2.v,[2 1 3]);
wu=permute(m2.w,[2 1 3]);

% vdc = vd(26,26,jcond);
% vuc = vu(26,26,jcond);
% udc = ud(26,26,jcond)-Um(jcond);
% uuc = uu(26,26,jcond)-Um(jcond);
% vmd=sqrt(vdc.^2+udc.^2);
% vmu=sqrt(vuc.^2+uuc.^2);
% vdcap=vdc/vmd;
% udcap=udc/vmd;
% vucap=vuc/vmu;
% uucap=uuc/vmu;

vuozu=vu.*ozu;
wuoyu=wu.*oyu;

sinmp=0.9567;
cosmp=-0.2909;
rdotmp = x.*cosmp + ry.*(sinmp);
rho = sqrt(r.^2-rdotmp.^2);


rcodx=ry.*ozd-z.*oyd;
rcody=z.*oxd-x.*ozd;
rcodz=x.*oyd-ry.*oxd;

rcodx=rcodx./(rho.*omd);
rcody=rcody./(rho.*omd);
rcodz=rcodz./(rho.*omd);

rcoux=ry.*ozu-z.*oyu;
rcouy=z.*oxu-x.*ozu;
rcouz=x.*oyu-ry.*oxu;

rcoux=rcoux./(rho.*omu);
rcouy=rcouy./(rho.*omu);
rcouz=rcouz./(rho.*omu);


omtd=rcodx.*cosmp+rcody.*sinmp;
omtu=rcoux.*cosmp+rcouy.*sinmp;

xb1=min([m1.xstart m2.xstart]);
yb1=min([m1.ystart m2.ystart]);
zb1=min([m1.zstart m2.zstart]);

xb2=max([m1.xend m2.xend]);
yb2=max([m1.yend m2.yend]);
zb2=max([m1.zend m2.zend]);

hslice=surf(linspace(-100,100,100),linspace(-100,100,100),zeros(100));
rotate(hslice,[-1,0,0],-16);
xd=get(hslice,'XData');
yd=get(hslice,'YData');
zd=get(hslice,'Zdata')+ycp-5;
delete(hslice)
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

  % [startZ,startX,startY]=meshgrid(0,[-80:20:80],ycp);
   [startZ,startX,startY]=meshgrid(0,0,ycp);

z_stop=-m1.zend-10;

vertsv = stream3(z,x,y,ozd,oxd,oyd,startZ,startX,startY);
vertsvn = stream3(z,x,y,-ozd,-oxd,-oyd,startZ,startX,startY);




%[startZ,startX,startY]=meshgrid(0.1./dnu,[-0.12:0.012:0.12]./dnu,ys);
% [startZ,startX,startY]=meshgrid(m1.zend+10,[m2.xstart-20:10:m2.xend+80],1.2*ys);

vertsv2 = stream3(z,x,y,ozu,oxu,oyu,startZ,startX,startY);
vertsv2n = stream3(z,x,y,-ozu,-oxu,-oyu,startZ,startX,startY);



% fd=figure;
% fd.Position=[x1 y1 width height];
%%
subplot(1,2,1)
% --- Colored streamlines by wd ---
% subplot(1,2,1)
% title('$(a) \langle v \rangle $ ','Interpreter','latex')
title('$(a) \omega_{\theta}$','Interpreter','latex')
verts = vertsv;
hold on
for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end

    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);

    % Interpolate w along streamline
    ws = interp3(z,x,y,omtd,xs,ys,zs);  

    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
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
    ws = interp3(z,x,y,omtd,xs,ys,zs);  

    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
end
quiver3(0,0,ycp,0,cosmp,sinmp,100)
% subplot(1,3,1)

hold on;

% isosurf=isosurface(z,x,y,ld,ltd);
% interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
%     isosurf.vertices(:, 1), ...
%     isosurf.vertices(:, 2), ...
%     isosurf.vertices(:, 3));
% cld=floor(max(100*interpColors,[],'all'))/100
lighting gouraud; 
% lightangle(-45,90)
 scatter3(0,0,ycp,50,'green','filled')
%sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% cld=max(abs(sl.CData(:)));
 %sl=slice(vertsv())
 % set(sl,'FaceColor','interp','FaceAlpha',0.8)
shading interp
% l=streamline(vertsv);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1);
% ln=streamline(vertsvn);
% set(ln, 'Color', 'k'); 
% set(ln,'LineWidth',1);
hold off

colormap redblue
axis equal
view(0,65);
zlim([yb1-10 yb2+10])
xlim([-100 100])
ylim([-100 100])
% c1=colorbar;
 % clim([-cld cld])
% c1.Ticks=[-cld, 0, cld];
%clim([-1 1]*1e-1)
%clim([-150 150])
% ylabel(c1," \omega_x^+ ")
% ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
set(gca,'FontSize',11)
%ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
clim([-1.5e-1 1.5e-1])
%%


% % % subplot(1,3,2)
% % % title('$(b) -\langle v \rangle \langle \omega_z \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vdozd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsvn;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vdozd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
% % % %%
% % % 
% % % subplot(1,3,3)
% % % title('$(b) -\langle v  \omega_z \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vozd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsvn;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vozd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
% % % 
% % % 
% % % 
% % % %%
% % % figure
% % % subplot(1,3,1)
% % % verts = vertsv;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsvn;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % title('(b) w')
% % % clim([-5e-3 5e-3])
% % % 
% % % %%
% % % 
% % % 
% % % subplot(1,3,2)
% % % title('$(b) \langle w \rangle \langle \omega_y \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wdoyd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsvn;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wdoyd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
% % % %%
% % % 
% % % subplot(1,3,3)
% % % title('$(b) \langle w  \omega_y \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,woyd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsvn;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,woyd,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])




%%

% figure

% --- Colored streamlines by wd ---
%subplot(1,3,1)
subplot(1,2,2)
title('$(b)  \omega_{\theta}  $ ','Interpreter','latex')

verts = vertsv2;
hold on
for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end
    
    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);
    
    % Interpolate w along streamline
    ws = interp3(z,x,y,omtu,xs,ys,zs);  
    
    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
end
su = [ys,zs];
verts = vertsv2n;
for k = 1:length(verts)
    thisStream = verts{k};
    if isempty(thisStream), continue; end
    
    xs = thisStream(:,1);
    ys = thisStream(:,2);
    zs = thisStream(:,3);
    
    % Interpolate w along streamline
    ws = interp3(z,x,y,omtu,xs,ys,zs);  
    
    % Plot colored line
    surface([xs xs],[ys ys],[zs zs],[ws ws], ...
        'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
end
 quiver3(0,0,ycp,0,cosmp,sinmp,100)

% subplot(1,3,1)

%hold on;
%
% isosurf=isosurface(z,x,y,ld,ltd);
% interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
%     isosurf.vertices(:, 1), ...
%     isosurf.vertices(:, 2), ...
%     isosurf.vertices(:, 3));
% cld=floor(max(100*interpColors,[],'all'))/100
lighting gouraud; 
% lightangle(-45,90)
 scatter3(0,0,ycp,50,'green','filled')
%sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
 % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% cld=max(abs(sl.CData(:)));
 %sl=slice(vertsv())
 % set(sl,'FaceColor','interp','FaceAlpha',0.8)
 % shading interp
% l=streamline(vertsv);
% set(l, 'Color', 'k'); 
% set(l,'LineWidth',1);
% ln=streamline(vertsvn);
% set(ln, 'Color', 'k'); 
% set(ln,'LineWidth',1);
hold off

colormap redblue
axis equal
view(0,65);
zlim([yb1-10 yb2+10])
xlim([-100 100])
ylim([-100 100])
% c1=colorbar;
 % clim([-cld cld])
% c1.Ticks=[-cld, 0, cld];
%clim([-1 1]*1e-1)
%clim([-150 150])
% ylabel(c1," \omega_x^+ ")
% ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
set(gca,'FontSize',11)
%ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
xlabel('$z^+$','interpreter','latex','FontSize',14)
ylabel('$x^+$','interpreter','latex','FontSize',14)
zlabel('$y^+$','interpreter','latex','FontSize',14)
clim([-5e-3 5e-3])
%%


% % % subplot(1,3,2)
% % % title('$(b) -\langle v \rangle \langle \omega_z \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv2;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vuozu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsv2n;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vuozu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
% % % %%
% % % 
% % % subplot(1,3,3)
% % % title('$(b) -\langle v  \omega_z \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv2;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vozu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsv2n;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,-vozu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
% % % 
% % % 
% % % 
% % % %%
% % % figure
% % % subplot(1,3,1)
% % % verts = vertsv2;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsv2n;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % title('(b) w')
% % % clim([-5e-3 5e-3])
% % % 
% % % %%
% % % 
% % % 
% % % subplot(1,3,2)
% % % title('$(b) \langle w \rangle \langle \omega_y \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv2;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wuoyu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsv2n;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,wuoyu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
% % % %%
% % % 
% % % subplot(1,3,3)
% % % title('$(b) \langle w  \omega_y \rangle$','Interpreter','latex')
% % % 
% % % verts = vertsv2;
% % % hold on
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,woyu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % verts = vertsv2n;
% % % for k = 1:length(verts)
% % %     thisStream = verts{k};
% % %     if isempty(thisStream), continue; end
% % % 
% % %     xs = thisStream(:,1);
% % %     ys = thisStream(:,2);
% % %     zs = thisStream(:,3);
% % % 
% % %     % Interpolate w along streamline
% % %     ws = interp3(z,x,y,woyu,xs,ys,zs);  
% % % 
% % %     % Plot colored line
% % %     surface([xs xs],[ys ys],[zs zs],[ws ws], ...
% % %         'FaceColor','none','EdgeColor','interp','LineWidth',1.5);
% % % end
% % % 
% % % % subplot(1,3,1)
% % % 
% % % %hold on;
% % % %
% % % % isosurf=isosurface(z,x,y,ld,ltd);
% % % % interpColors = interp3(z, x, y, -vozd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % % %     isosurf.vertices(:, 1), ...
% % % %     isosurf.vertices(:, 2), ...
% % % %     isosurf.vertices(:, 3));
% % % % cld=floor(max(100*interpColors,[],'all'))/100
% % % lighting gouraud; 
% % % % lightangle(-45,90)
% % %  scatter3(0,0,ycp,50,'green','filled')
% % % %sl=slice(z,x,y,-vozd./ut^2,[],[],[ys])
% % %  % sl=slice(z,x,y,woyd./ut^2,xd,yd,zd)
% % % % cld=max(abs(sl.CData(:)));
% % %  %sl=slice(vertsv())
% % %  % set(sl,'FaceColor','interp','FaceAlpha',0.8)
% % %  % shading interp
% % % % l=streamline(vertsv);
% % % % set(l, 'Color', 'k'); 
% % % % set(l,'LineWidth',1);
% % % % ln=streamline(vertsvn);
% % % % set(ln, 'Color', 'k'); 
% % % % set(ln,'LineWidth',1);
% % % hold off
% % % 
% % % colormap redblue
% % % axis equal
% % % view(0,65);
% % % zlim([yb1-10 yb2+10])
% % % xlim([-100 100])
% % % ylim([-100 100])
% % % % c1=colorbar;
% % %  % clim([-cld cld])
% % % % c1.Ticks=[-cld, 0, cld];
% % % %clim([-1 1]*1e-1)
% % % %clim([-150 150])
% % % % ylabel(c1," \omega_x^+ ")
% % % % ylabel(c1,"v \omega_z /(-u_{\tau}^2/H)")
% % % set(gca,'FontSize',11)
% % % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % % clim([-5e-3 5e-3])
%%
% maskd=(sd(:,1)>=-10)&(sd(:,1)<=10);
% sdf=sd(maskd,:);
% pd=polyfit(sdf(:,1),sdf(:,2),1);
% ydfit=polyval(pd,sdf(:,1));
% 
% masku=(su(:,1)>=-10)&(su(:,1)<=10);
% suf=su(masku,:);
% pu=polyfit(suf(:,1),suf(:,2),1);
% yufit=polyval(pu,suf(:,1));
% 
% md=pd(:,1)
% mu=pu(:,1)
% mudav=0.5*(md+mu);
% mdp=-1/md;
% mup=-1/mu;
% 
% atan(mdp)./pi*180
% atan(mup)./pi*180
% 
% mp=-1/mudav;
% sinmp=-mp./(sqrt(1+mp^2))
% cosmp=-1./(sqrt(1+mp^2))
% 
% sinmdp=-mdp./(sqrt(1+mdp^2))
% cosmdp=-1./(sqrt(1+mdp^2))
% 
% sinmup=-mup./(sqrt(1+mup^2))
% cosmup=-1./(sqrt(1+mup^2))
%%

% % subplot(1,3,2)
% % 
% % hold on;
% % %
% % % isosurf=isosurface(z,x,y,ld,ltd);
% % % interpColors = interp3(z, x, y, woyd./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % %     isosurf.vertices(:, 1), ...
% % %     isosurf.vertices(:, 2), ...
% % %     isosurf.vertices(:, 3));
% % % cld=floor(max(100*interpColors,[],'all'))/100
% % 
% % lighting gouraud; 
% % % lightangle(-45,90)
% % scatter3(0,0,ys,30,'green','filled')
% % 
% %  sl=slice(z,x,y,wd./ut,xd,yd,zd)
% % cld=max(abs(sl.CData(:)));
% % 
% % 
% %  set(sl,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.8)
% % shading interp
% % l=streamline(vertsv);
% % set(l, 'Color', 'k'); 
% % set(l,'LineWidth',1);
% % ln=streamline(vertsvn);
% % set(ln, 'Color', 'k'); 
% % set(ln,'LineWidth',1);
% % hold off
% % 
% % colormap redblue
% % axis equal
% % view(0,65);
% % zlim([yb1-10 yb2+10])
% % xlim([-100 100])
% % ylim([-100 100])
% % % c1=colorbar;
% %  clim([-cld cld])
% % % c1.Ticks=[-cld, 0, cld];
% % %clim([-1 1]*1e-1)
% % %clim([-150 150])
% % % ylabel(c1," \omega_y^+ ")
% % % ylabel(c1,"-w \omega_y /(-u_{\tau}^2/H)")
% % 
% % set(gca,'FontSize',11)
% % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % title('(b)')
% % %%
% % subplot(1,3,3)
% % 
% % hold on;
% % %
% % % isosurf=isosurface(z,x,y,ld,ltd);
% % % interpColors = interp3(z, x, y, -nld./ut^2,...wd.*0, ...interpColors = interp3(z, x, y, wd.*0, ...
% % %     isosurf.vertices(:, 1), ...
% % %     isosurf.vertices(:, 2), ...
% % %     isosurf.vertices(:, 3));
% % % cld=floor(max(100*interpColors,[],'all'))/100
% % lighting gouraud; 
% % % lightangle(-45,90)
% % scatter3(0,0,ys,30,'green','filled')
% % 
% %  sl=slice(z,x,y,wd.*oyd./ut^2,xd,yd,zd)
% %  cld=max(abs(sl.CData(:)));
% %  set(sl,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.8)
% %  %shading interp
% % l=streamline(vertsv);
% % set(l, 'Color', 'k'); 
% % set(l,'LineWidth',1);
% % ln=streamline(vertsvn);
% % set(ln, 'Color', 'k'); 
% % set(ln,'LineWidth',1);
% % hold off
% % 
% % colormap redblue
% % axis equal
% % view(0,65);
% % zlim([yb1-10 yb2+10])
% % xlim([-100 100])
% % ylim([-100 100])
% % % c1=colorbar;
% % clim([-cld cld])
% % % c1.Ticks=[-cld, 0, cld];
% % %clim([-1 1]*1e-1)
% % %clim([-150 150])
% % % ylabel(c1," (v \omega_z - w \omega_y)/(-u_{\tau}^2/H) ")
% % set(gca,'FontSize',11)
% % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % zlabel('$y^+$','interpreter','latex','FontSize',14)
% % title('(c)')
% % 
% % %%
% % figure
% % hold on;
% % %
% % isosurf=isosurface(z,x,y,lu,ltu);
% % % interpColors = interp3(z, x, y, -nlu./ut^2,...dnu*oxu./ut,...wd.*0, ...interpColors = interp3(z, x, y, -nlu./ut^2,...wd.*0, ...
% % %     isosurf.vertices(:, 1), ...
% % %     isosurf.vertices(:, 2), ...
% % %     isosurf.vertices(:, 3));
% % % clu=floor(max(100*interpColors,[],'all'))/100
% % 
% % % p = patch(isosurf);
% % % p.FaceColor = 'interp';        % Interpolated color
% % % p.EdgeColor = 'none';          % Remove edges
% % % p.FaceVertexCData = interpColors;    % Assign interpolated colors
% % % p.FaceAlpha = 0.8;  
% % camlight;                      % Add lighting
% % %light('Position', [0, 0, 0.5], 'Style', 'infinite');  % Light above and to the right
% % lighting gouraud; 
% % lightangle(-45,90)
% % scatter3(0,0,ys,100,'green','filled')
% % l=streamline(vertsv2);
% % 
% % % slice(z,x,y,-nlu./ut^2,[0],[],[])
% % % shading interp
% % 
% % 
% % set(l, 'Color', 'k'); 
% % set(l,'LineWidth',1);
% % ln=streamline(vertsv2n);
% % set(ln, 'Color', 'm'); 
% % set(ln,'LineWidth',1);
% % colormap redblue
% % axis equal
% % %view(45,45);
% % view(135,45)
% % % ylim([-150 150])
% % % xlim([-100 100])
% % % %zlim([0 200])
% % % zlim([0 150])
% % % yticks([-150:50:150])
% % 
% % c=colorbar;
% % %clim([-1 1]*1e-1)
% % %clim([-150 150])
% % %ylabel(c,"u'")
% % set(gca,'FontSize',11)
% % %ylabel(c,"$H\omega_x/u_{\tau}$",'interpreter','latex','FontSize',14)
% % %ylabel(c,"$(v\omega_z-w\omega_y)/(-u_{\tau}^2/H)$",'interpreter','latex','FontSize',14)
% % xlabel('$z^+$','interpreter','latex','FontSize',14)
% % ylabel('$x^+$','interpreter','latex','FontSize',14)
% % 
% % % zlim([50 400])
% % 
% % % zlim([m2.ystart-10 m2.yend+10])
% % % xlim([m2.zstart-10 m2.zend+10])
% % % ylim([m2.xstart-10 m2.xend+40])
% % 
% % % zlim([m1.ystart-10 m1.yend+10])
% % % xlim([m1.zstart-10 m1.zend+10])
% % % ylim([m1.xstart-20 m1.xend+20])
% % % 
% % % zticks([round(m2.ystart) round(m2.yend)])
% % % xticks([round(m2.zstart) ,0,round(m2.zend)])
% % % yticks([round(m2.xstart) ,0,round(m2.xend)])
% % 
% % zlim([yb1-10 yb2+10])
% % xlim([zb1-10 zb2+10])
% % ylim([xb1-20 xb2+20])
% % 
% % zticks([round(m2.ystart) round(m2.yend)])
% % xticks([round(m2.zstart) ,0,round(m2.zend)])
% % yticks([round(m2.xstart) ,0,round(m2.xend)])
% % 
% % c=colorbar;
% % % clim([-cld cld])
% %  % c.Ticks=[-cld, 0, cld];
% %  ylabel(c,"$(v\omega_z-w\omega_y)/(-u_{\tau}^2/H)$",'interpreter','latex','FontSize',14)
% % 
% % %zlabel('$y^+$','interpreter','latex','FontSize',14)
% % title('(b)')
% % %%
% % % fdn=sprintf("eddy_j_rings_slice_cond_%03d.fig",jcond)
% % % saveas(fd,fdn)
