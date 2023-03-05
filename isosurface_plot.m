%function[]=print_figure( outfile, file_format )  
%clear


close all
Nx=2048;
Ny=512;
Nz=1536;
Lx=  8*pi;
Lz = 3*pi;
xp = [0:Nx-1]*Lx/(Nx);
zp=  [0:1:Nz-1]*Lz/(Nz);
mb=matfile('bsplinedata.mat');
%load('isosurface_data_v.mat')
%load('transfer_profile.mat');
%mcv=matfile('Full_conv_visc.mat');
yp = mb.yv;
ut = 0.0499;
dnu=1.0006e-3;
nx=164;
%ny=130;
nz=165*3;
nystart=1;
nyend=180;
time=3;
ny=nyend-nystart+1;
[x,z,y] = meshgrid(zp(1:nz),xp(1:nx),yp(nystart:nyend)+1);
%time=3;
%for time=3:20
fvel=sprintf("vel_%03d.mat",time);
%fl=sprintf("uv_%03d.mat",time);
%fl=sprintf("lambda_%03d",time);
%fl=sprintf("lambda_filtered_15x1p25z_%03d.mat",time);
%fl="Transfer_filtered_002.mat";
ft=sprintf("Transfer_%03d.mat",time);
%ft=sprintf("Transfer_filtered_fluc_x15z1p25_%03d.mat",time);
%fvo=sprintf("vortfiltered_%03d.mat",time);
%ftitle=sprintf(" lambda_2=-10, time=%03d",time);
%figname=sprintf("l10_%03d.png",time);
%ml=matfile(fl);
mcv=matfile(ft);
%mvo=matfile(fvo);
mvel=matfile(fvel);
%mdata=matfile('isosurface_data_lrms_f15x1p25z.mat','Writable',true)
mdata=matfile('isosurface_data_vrms','Writable',true)
load('JHTDB_RET1000.mat');
%load('lambda_stats.mat')
%ml=matfile('conditioned_aver_field_l.mat');

 %conv=mcv.convective_x(1:120,1:164,1:327);
  %conv=conv(1:120,1:164,1:327);
%visc=mcv.viscous_x(1:120,1:164,1:327);
conv = mcv.convective_x(nystart:nyend,1:nx,1:nz) ;
visc = mcv.viscous_x(nystart:nyend,1:nx,1:nz) ;

%omegaz=mvo.omegazf(nystart:nyend,1:nx,1:nz) ;
%+ mcv.viscous_x(nystart:nyend,1:nx,1:nz);
%  ml=matfile('Full_lambda2_test.mat');
%uprofile=mean( mean( mvel.u,3),2);
%up=mvel.u(nystart:nyend,1:nx,1:nz) - uprofile(nystart:nyend);

%up=up./ut;
%ur=sqrt(JHTDB_RET1000(nystart:nyend,4));

%l=up./ur;
%[ccm,~,binsm]=histcounts(up(2:ny,1:nx,1:nz),[-inf 0]);
%[ccp,~,binsp]=histcounts(up(2:ny,1:nx,1:nz),[0 inf]);

%	l= mvel.v(nystart:nyend,1:nx,1:nz).*( mvel.u(nystart:nyend,1:nx,1:nz) - uprofile(nystart:nyend) );
%	l=l./ut^2;
%l=l./ut;
%	urvr=sqrt(JHTDB_RET1000(:,4).*JHTDB_RET1000(:,5));
%	urvr=urvr(nystart:nyend);
%	l=l./urvr;

%    l=(ml.uv(:,1:nx,1:nz))./ut^2;
%tot=mcv.total(nystart:nyend,1:nx,1:nz);
    %%
%l=ml.lambda2(nystart:nyend,1:nx,1:nz)./lrms(nystart:nyend);
%tot=mcv.total(:,1:nx,1:nz);
%tot=mcv.convective_x(:,1:nx,1:nz);
%ly=
  %  lcond=ml.lcond(2:ny,1:nx,1:nz);
  %  tot=ml.viscond(2:ny,1:nx,1:nz)+ml.convcond(2:ny,1:nx,1:nz);
  %  l=lcond;
  %  l(isnan(l))=0;
  
  %cp=permute(conv,[2 3 1]);
  %vp=permute(visc,[2 3 1]);

  %lp = permute(l,[2 3 1]);
 %ylp = permute(l,[2 3 1]).*(y.^2);

%	lpm=permute(binsm.*l,[2 3 1]);
%	lpp=permute(binsp.*l,[2 3 1]);

  %tp = permute(tot,[2 3 1]);

mdata.x=x;
mdata.y=y;
mdata.z=z;
mdata.l=l;
mdata.conv=conv;
mdata.visc=visc;
%mdata.omegaz=omegaz;
%cp2 = cp(1:200,1:300,1:100);
%data = sqrt(x.^2 + y.^2 + z.^2);
%cdata=Z;
%cdata = smooth3(rand(size(data)),'box',7);
%%

%ys=y./dnu;
%[faces,verts,colors] = isosurface(z,x,y./dnu, lp,-1,tp./ut^2);
%[facesm,vertsm,colorsm] = isosurface(z,x,y./dnu,lpm./ut^2,-1,tp./ut^2);
%[facesp,vertsp,colorsp] = isosurface(z,x,y./dnu,lpp./ut^2,-1,tp./ut^2);

%
%i=isosurface(X,Z,Y./dnu,conv,10);
%
%subplot(1,2,1)
%tiledlayout(1,2);

%nexttile
%p=patch('Vertices', verts, 'Faces', faces, ...
%      'FaceVertexCData', colors, ...
%      'FaceColor','interp', ...
%      'edgecolor', 'interp');

%p.FaceColor = 'interp';
%p.EdgeColor = 'none';
% % %view(150,30)
% % %daspect([1 1 1])
% % %axis tight
%camlight
%lighting gouraud
%xlabel('x/h')
%ylabel('z/h')
%zlabel('y^+')
%view(-30,30)
%%%
%pbaspect([2 1.5 0.5])
%view(3)
%lightangle(-15,25)
%ylim([0 1])
%xlim([0 2])

%
%ylim([0 1.5])
%zlim([1 500])
%title('$\bar{U}\cdot(\mathbf{u}\times\mathbf{\omega}  ) =100u_{\tau^2}\bar{U}/h$',...
%    'interpreter','latex')

% title('$-\nu \bar{U}\cdot \nabla \times \mathbf{\omega}=-100u_{\tau^2}\bar{U}/h$',...
%     'interpreter','latex')
%title('$Transfer = -u_{\tau}^2\bar{U}/h$','interpreter','latex')
%title('$\lambda_2^+=-0.5$','interpreter','latex')
%h=colorbar;
%ylabel(h,'$T/(u_{\tau}^2 \bar{U}/h)$','interpreter','latex')
%title('$ \bar{ \mathcal{I}  }  /(u_{\tau}^2 \bar{U}/h)$','interpreter','latex')
%caxis([-20 20])
%colormap jet
%set(gca,'FontSize',12);
%title("$\lambda_2<-\lambda_{rms}$",'interpreter','latex')
%title("$v<-v_{rms}$",'interpreter','latex')

    


%title('$\bar{U}\cdot(\mathbf{u}\times\mathbf{\omega}  ) =100u_{\tau^2}\bar{U}/h$',...
%    'interpreter','latex')

% title('$-\nu \bar{U}\cdot \nabla \times \mathbf{\omega}=-100u_{\tau^2}\bar{U}/h$',...
%     'interpreter','latex')
%title('$Transfer = -u_{\tau}^2\bar{U}/h$','interpreter','latex')
%title('$\lambda_2^+=-0.5$','interpreter','latex')
% h=colorbar;
% ylabel(h,'$T/(u_{\tau}^2 \bar{U}/h)$','interpreter','latex')
% 
% caxis([-200 200])
% colormap jet
% set(gca,'FontSize',12)


%%
%x0=100;
%y0=50;
%width=700;
%height=500;
%xlim([0 2])
%set(gcf,'position',[x0,y0,x0+width,y0+height])
%
%title("$ h^2 \bar{\lambda}_2/u_{\tau}^2=-70$",'interpreter','latex')
%%
%title("$ h^2 \bar{\lambda}_2/u_{\tau}^2=-70$",'interpreter','latex')

% subplot(1,2,2)
% %nexttile
% p=patch('Vertices', verts, 'Faces', faces, ...
%       'FaceVertexCData', colors, ...
%       'FaceColor','interp', ...
%       'edgecolor', 'interp');
% 
% p.FaceColor = 'interp';
% p.EdgeColor = 'none';
% % % %view(150,30)
% % % %daspect([1 1 1])
% % % %axis tight
% camlight
% lighting gouraud
% xlabel('x/h')
% ylabel('z/h')
% zlabel('y^+')
% view(30,30)
% pbaspect([2 1 1])
% %view(3)
% lightangle(-15,25)
% 
% xlim([0 2])
% ylim([0 1])
% zlim([1 1000])
% %title('$\bar{U}\cdot(\mathbf{u}\times\mathbf{\omega}  ) =100u_{\tau^2}\bar{U}/h$',...
% %    'interpreter','latex')
% 
% % title('$-\nu \bar{U}\cdot \nabla \times \mathbf{\omega}=-100u_{\tau^2}\bar{U}/h$',...
% %     'interpreter','latex')
% %title('$Transfer = -u_{\tau}^2\bar{U}/h$','interpreter','latex')
% %title('$\lambda_2^+=-0.5$','interpreter','latex')
% h=colorbar;
% ylabel(h,'$T/(u_{\tau}^2 \bar{U}/h)$','interpreter','latex')
% 
% caxis([-200 200])
% colormap jet
% set(gca,'FontSize',12)
% 
% 

%title("$u'v'<-1.75u_{rms}v_{rms}$",interpreter,'latex')

%title('$\bar{U}\cdot(\mathbf{u}\times\mathbf{\omega}  ) =100u_{\tau^2}\bar{U}/h$',...
%    'interpreter','latex')

% title('$-\nu \bar{U}\cdot \nabla \times \mathbf{\omega}=-100u_{\tau^2}\bar{U}/h$',...
%     'interpreter','latex')
%title('$Transfer = -u_{\tau}^2\bar{U}/h$','interpreter','latex')
%title('$\lambda_2^+=-0.5$','interpreter','latex')
% h=colorbar;
% ylabel(h,'y^+')
 %ylabel(h,'$ \bar{ \mathcal{I}  }  /(u_{\tau}^2 \bar{U}/h)$','interpreter','latex')

 %ylabel(h,"$ \overline{u_2\omega_3'-u_3\omega_2}/(u_{\tau}^2 \bar{U}/h) $",'interpreter','latex')
 
%print(outfile,file_format); 
 
 %end
 % 
% caxis([-200 200])
% colormap jet
% set(gca,'FontSize',12)


%%
% x0=100;
% y0=50;
% width=1400;
% height=500;
% xlim([0 2])
% set(gcf,'position',[x0,y0,x0+width,y0+height])
% %
% title("$\lambda=-1\bar{U}^2/H^2$",'interpreter','latex')
% title(ftitle);
%title('$\lambda_2^+=-1$','interpreter','latex')
%%
%print(gcf,'lambda_total_iso_1_back.png','-dpng','-r300');
% print(gcf,figname,'-dpng','-r300');
% close all
%end
%%
% mf=matfile('Flag_vort_quarter.mat');
% f=mf.flag(41:130,1:83,1:164);
% fp=permute(f,[2 3 1]);
% fl=double(reshape(fp,[],1));
% xl=reshape(x,[],1)./fl;
% yl=reshape(y./dnu,[],1)./fl;
% zl=reshape(z,[],1)./fl;
% tl=reshape(tp./ut^2,[],1)./fl;
% figure
% scatter3(zl,xl,yl,20,tl,'filled','s')
% pbaspect([1 1 0.25])
% %view(3)
% lightangle(-30,30)
% xlim([0 1])
% ylim([0 1])
% zlim([30 300])
% 
% 
% title('$\lambda_2^+=-0.25$','interpreter','latex')
% h=colorbar;
% caxis([-250 250])
% colormap jet
% ylabel(h,'total')
% view(-45,45)
