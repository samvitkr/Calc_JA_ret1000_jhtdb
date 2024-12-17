Nx=2048;
Nz=1536;
m0=matfile('smooth_2d_spectra_indy2_w0.mat','Writable',true)
m1=matfile('smooth_2d_spectra_indy2_w1.mat','Writable',true)
m2=matfile('smooth_2d_spectra_indy2_w2.mat','Writable',true)
m3=matfile('smooth_2d_spectra_indy2.mat'   ,'Writable',true)
m4=matfile('smooth_2d_spectra_indy2_w4.mat','Writable',true)
m5=matfile('smooth_2d_spectra_indy2_w5.mat','Writable',true)
m6=matfile('smooth_2d_spectra_indy2_w7.mat','Writable',true)
m7=matfile('smooth_2d_spectra_indy2_w0.mat','Writable',true)


size(m7.wconv)

specl2=zeros(6,30);
d0=squeeze(rms(rms(abs(m1.wconv(1:Nx/2,1:Nz/2,:)-m0.wconv(1:Nx/2,1:Nz/2,:)),1),2));
d1=squeeze(rms(rms(abs(m2.wconv(1:Nx/2,1:Nz/2,:)-m1.wconv(1:Nx/2,1:Nz/2,:)),1),2));
d2=squeeze(rms(rms(abs(m3.wconv(1:Nx/2,1:Nz/2,:)-m2.wconv(1:Nx/2,1:Nz/2,:)),1),2));
d3=squeeze(rms(rms(abs(m4.wconv(1:Nx/2,1:Nz/2,:)-m3.wconv(1:Nx/2,1:Nz/2,:)),1),2));
d4=squeeze(rms(rms(abs(m5.wconv(1:Nx/2,1:Nz/2,:)-m4.wconv(1:Nx/2,1:Nz/2,:)),1),2));
d5=squeeze(rms(rms(abs(m6.wconv(1:Nx/2,1:Nz/2,:)-m5.wconv(1:Nx/2,1:Nz/2,:)),1),2));
d6=squeeze(rms(rms(abs(m7.wconv(1:Nx/2,1:Nz/2,:)-m6.wconv(1:Nx/2,1:Nz/2,:)),1),2));

m=matfile('spec_dist_l2.mat','Writable',true)
m.d0=d0;
m.d1=d1;
m.d2=d2;
m.d3=d3;
m.d4=d4;
m.d5=d5;
m.d6=d6;

