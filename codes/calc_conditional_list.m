close all
clear
Nx=2048;
Ny=512;
Nz=1536;

Lx=8*pi;
Lz=3*pi;

jcond=130;
jct=Ny-jcond+1;

x=(Lx*[0:Nx-1]/Nx-Lx/2);
z=(Lz*[0:Nz-1]/Nz-Lz/2);

fp=sprintf("../data/conditionalp_jcond_1_%03d.mat",jcond);
%fc=sprintf("../data/test.mat")
mp=matfile(fp,'Writable',true);
event_list_p=mp.event;
% clear fp mp

fn=sprintf("../data/conditionaln_jcond_1_%03d.mat",jcond);
%fc=sprintf("../data/test.mat")
mn=matfile(fn,'Writable',true);
event_list_n=mn.event;
% clear fn mn
dist_n=event_list_n(:,1).*0;
dist_p=event_list_p(:,1).*0;
theta_n=event_list_n(:,1).*0;
theta_p=event_list_p(:,1).*0;
tstart=1;
tend=10;
tstep=1;
 for time=tstart:tstep:tend
     time
     eptime=event_list_p(event_list_p(:,4)==time,:); 
     entime=event_list_n(event_list_n(:,4)==time,:); 
    epbot=eptime(eptime(:,3)==jcond,:);
    eptop=eptime(eptime(:,3)==jct,:);
    enbot=entime(entime(:,3)==jcond,:);
    entop=entime(entime(:,3)==jct,:) ;

    %%
    n1=length(enbot);
    n2=length(epbot);
    %%
    dbot=zeros(n1,1);
    thetabot=dbot;
    dx=zeros(n2,3);
    dz=zeros(n2,3);
    dist=zeros(n2,1);
    for i=1:n1
        z1=z(enbot(i,1));
        x1=x(enbot(i,2));

        zp=z(epbot(:,1));
        xp=x(epbot(:,2));

        dz(:,1)=abs(zp-z1);
        dz(:,2)=Lz-abs(zp-z1);
        dz(:,3)=min(dz(:,1:2),[],2);

        dx(:,1)=abs(xp-x1);
        dx(:,2)=Lx-abs(xp-x1);
        dx(:,3)= min(dx(:,1:2),[],2) ;

        dist=sqrt(dz(:,3).^2 + dx(:,3).^2);
        [dbot(i), im]=min( dist ,[],'all');
        thetabot(i) = atan(abs(dx(im,3)/dz(im,3)));
    end
    %%
    dbotp=zeros(n2,1);
    thetabotp=dbotp;
    dx=zeros(n1,3);
    dz=zeros(n1,3);
    distp=zeros(n1,1);
    for i=1:n2
        z1=z(epbot(i,1));
        x1=x(epbot(i,2));

        zp=z(enbot(:,1));
        xp=x(enbot(:,2));

        dz(:,1)=abs(zp-z1);
        dz(:,2)=Lz-abs(zp-z1);
        dz(:,3)=min(dz(:,1:2),[],2);

        dx(:,1)=abs(xp-x1);
        dx(:,2)=Lx-abs(xp-x1);
        dx(:,3)= min(dx(:,1:2),[],2) ;

        distp=sqrt(dz(:,3).^2 + dx(:,3).^2);
        [dbotp(i),imp]=min( distp ,[],'all');
        thetabotp(i)=atan(abs(dx(imp,3)/dz(imp,3)));
    end

 %%
    n1=length(entop);
    n2=length(eptop);
    %%
    dtop=zeros(n1,1);
    thetatop=dtop;
    dx=zeros(n2,3);
    dz=zeros(n2,3);
    dist=zeros(n2,1);
    for i=1:n1
        z1=z(entop(i,1));
        x1=x(entop(i,2));

        zp=z(eptop(:,1));
        xp=x(eptop(:,2));

        dz(:,1)=abs(zp-z1);
        dz(:,2)=Lz-abs(zp-z1);
        dz(:,3)=min(dz(:,1:2),[],2);

        dx(:,1)=abs(xp-x1);
        dx(:,2)=Lx-abs(xp-x1);
        dx(:,3)= min(dx(:,1:2),[],2) ;

        dist=sqrt(dz(:,3).^2 + dx(:,3).^2);
        [dtop(i),imt]=min( dist ,[],'all');
        thetatop(i)=atan(abs(dx(imt,3)/dz(imt,3)));

    end
    d=[dbot;dtop];
    theta=[thetabot;thetatop];
 dist_n(event_list_n(:,4)==time)=d; 
 theta_n(event_list_n(:,4)==time)=theta; 

 %%
 dtopp=zeros(n2,1);
 thetatopp=dtopp;
    dx=zeros(n1,3);
    dz=zeros(n1,3);
    distp=zeros(n1,1);
    for i=1:n2
        z1=z(eptop(i,1));
        x1=x(eptop(i,2));

        zp=z(entop(:,1));
        xp=x(entop(:,2));

        dz(:,1)=abs(zp-z1);
        dz(:,2)=Lz-abs(zp-z1);
        dz(:,3)=min(dz(:,1:2),[],2);

        dx(:,1)=abs(xp-x1);
        dx(:,2)=Lx-abs(xp-x1);
        dx(:,3)= min(dx(:,1:2),[],2) ;

        distp=sqrt(dz(:,3).^2 + dx(:,3).^2);
        [dtopp(i),imtp]=min( distp ,[],'all');
        thetatopp(i)=atan(abs(dx(imtp,3)/dz(imtp,3)));
     end
 dp=[dbotp;dtopp];
thetap=[thetabotp;thetatopp];

dist_p(event_list_p(:,4)==time)= dp;
theta_p(event_list_p(:,4)==time)=thetap;
  end

%   mp.dist_from_closest_n=dist_p;
%   mn.dist_from_closest_p=dist_n;
%  fr=sprintf('../data/corr_v_reflect_j_%03d.mat',jcond);
% mr=matfile(fr)
% rvv=fftshift(mr.Rvv(:,:,jcond));
% 
 fd = sprintf('./dist_corr_cond_j_%03d.mat',jcond)
 md=matfile(fd,'Writable',true)
% md.dist_p=dist_p;
% md.dist_n=dist_n;
md.theta_n=theta_n;
md.theta_p=theta_p;
% md.rvv=rvv;
% 
% [X,Z]=meshgrid(x,z);
% md.X=X;
% md.Z=Z;