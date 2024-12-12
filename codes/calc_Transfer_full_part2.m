%%
clear
mt = matfile('Full_Transfer.mat','Writable',true);
save('Full_Transfer.mat','-v7.3');
mvel = matfile('Full_velfield.mat','Writable',true);
save('Full_velfield.mat','-v7.3');
mvelg = matfile('Full_velgradfield.mat','Writable',true);
save('Full_velgradfield.mat','-v7.3');
[Ny,Nx,Nz]=size(mvel.ufield);
for i=1:Nx
uyz=squeeze(mvel.ufield(:,i,:));
vyz=squeeze(mvel.vfield(:,i,:));
wyz=squeeze(mvel.wfield(:,i,:));

fzu(:,:)=fft(uyz(:,:).').';
fzv(:,:)=fft(vyz(:,:).').';
fzw(:,:)=fft(wyz(:,:).').';

dfzu(:,:)=(fzu(:,:)).*(1i*kz);
d2fzu(:,:)=(fzu(:,:)).*(-kz.^2);
dfzv(:,:)=(fzv(:,:)).*(1i*kz);
dfzw(:,:)=(fzw(:,:)).*(1i*kz);

dudzslice(:,:) = ifft(dfzu(:,:).').';
dvdzslice(:,:) = ifft(dfzv(:,:).').';
dwdzslice(:,:) = ifft(dfzw(:,:).').';

d2udz2(:,:) = ifft(d2fzu(:,:).').';
convective_xy(:,1,:)=-wyz.*dudzslice;
viscous_xy(:,1,:)=nu*d2udz2;
mt.convective(:,i,:)=mt.convective(:,i,:)+convective_xy;
mt.viscous(:,i,:)=mt.viscous(:,i,:)+viscous_xy;
mvelg.dudz(:,:,k)=single(dudzslice);
mvelg.dvdz(:,:,k)=single(dvdzslice);
mvelg.dwdz(:,:,k)=single(dwdzslice);
end
uyav=trapz(zp',uyz,2);
ubar=trapz(yp,uyav)/(2*3*pi);
mvel.ubar=ubar;
% clear uyz vyz wyz ...
%       fzu fzv fzw...
%       dfzu d2fzu dudzu d2udz2 convective_xy viscous_xy

 
