m3=matfile('vel_003.mat');
m4=matfile('vel_004.mat');
m1=matfile('vel_001.mat');
u1=squeeze(m1.u(100,100,:));
u3=squeeze(m3.u(100,100,:));
u4=squeeze(m4.u(100,100,:));
close all
hold on
plot(u1,'-k')
% plot(u3,'-b')
% plot(u4,'-r')
for i=1:5
    xline(256*i)
    xline(256*i-64)
end
hold off