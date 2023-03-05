m=matfile('vel_001.mat');
u=squeeze(m.u(200,:,20));
v=squeeze(m.v(200,:,20));
fu=fft(u);
fv=fft(v);
fuvdirect=fft(u.*v-mean( u.*v ));
fuv=(fu.*conj(fv));


