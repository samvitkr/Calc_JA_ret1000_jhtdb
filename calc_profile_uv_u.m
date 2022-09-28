Ny=512;
Nt=3;
m=matfile('transfer_profile.mat','Writable',true);
uv_mean=m.muv;
v_mean=m.mv;
u_mean=m.mu;

for time=10:19
	ft=sprintf("vel_%03d.mat",time)
	mt=matfile(ft);
	muv=mean(mean( (mt.u).*(mt.v),3),2);
	mu =mean(mean( mt.u,3),2);
	mv =mean(mean( mt.v,3),2);
	uv_mean( :,time )=muv;
	v_mean( :,time ) =mv;
	u_mean( :,time ) =mu;
end
%m=matfile('transfer_profile.mat','Writable',true);
m.muv=uv_mean;
m.mv=v_mean;
m.mu=u_mean;
