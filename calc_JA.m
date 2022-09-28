 fn=sprintf("uinflow_013.mat");
 mn=matfile(fn,'Writable',true);
 mn.timeset=[ 12.9675,13.0000, 13.0325 ];
 load('bsplinedata.mat')
 u=mn.ufield;
u1=mean(u(:,:,1),2);
u2=mean(u(:,:,2),2);
u3=mean(u(:,:,3),2);
Pin=mean(mn.pin(:,:),2);
Pout=mean(mn.pout,2);

U1=trapz(yv,u1)/2;
U2=trapz(yv,u2)/2;
U3=trapz(yv,u3)/2;
dUdt=(U3-U1)/( 13.0325-12.9675 );

mc=matfile("Transfer_013.mat");
tfield=mc.total;
T=mean(mean(tfield,3),2);
ret=9.9984509e+02;
m=matfile("JAcheck_013.mat",'Writable',true);
m.U1=u1;
m.U2=u2;
m.U3=u3;
m.T=T;
m.dUdt=dUdt;
m.ret=ret;
m.timeset=[ 12.9675,13.0000, 13.0325 ];
m.pin=Pin;
m.pout=Pout;





