%jcond=105;
ut = 0.0499;
dnu=1.0006e-3;
ret=1/dnu;

jcset=[47 54 71 105 130]
jj=4;
for counter=5:5
	jcond=jcset(jj)
	%fn=sprintf("lse_zslice_cond_j_%03d.mat",jcond)
	%fn=sprintf('../data/conditionaln_dudysplit_zslice_j_%03d.mat',jcond)
	fn=sprintf('../data/conditional_zslice_inst_j_%03d_%02d.mat',jcond,counter);
	m=matfile(fn,'Writable',true);
	[n1 n2]=size(m.ud);
	lcid = m.ud.*0;
	lciu=lcid;
	dd=zeros(2,2);
	du=zeros(2,2);
	for i =1:n1
	    for j=1:n2
	    dd(1,1)=m.dudxd(i,j);
	    dd(1,2)=m.dudyd(i,j);
	    dd(2,1)=m.dvdxd(i,j);
	    dd(2,2)=m.dvdyd(i,j);
	    l=max(abs(imag(eig(dd))));
	    lcid(i,j)=l;
	
	    % du(1,1)=m.dudxu(i,j);
	    % du(1,2)=m.dudyu(i,j);
	    % du(2,1)=m.dvdxu(i,j);
	    % du(2,2)=m.dvdyu(i,j);
	    % l=max(abs(imag(eig(du))));
	    % lciu(i,j)=l;
	
	
	    end
	end
	m.lcid=lcid;
	m.lciu=lciu;	
	%fn=sprintf('../data/conditionalp_dudysplit_zslice_j_%03d.mat',jcond)
	%m=matfile(fn,'Writable',true);
	%[n1 n2]=size(m.ud);
	%lcid = m.ud.*0;
	%lciu=lcid;
	%dd=zeros(2,2);
	%du=zeros(2,2);
	%for i =1:n1
	%    for j=1:n2
	%    dd(1,1)=m.dudxd(i,j);
	%    dd(1,2)=m.dudyd(i,j);
	%    dd(2,1)=m.dvdxd(i,j);
	%    dd(2,2)=m.dvdyd(i,j);
	%    l=max(abs(imag(eig(dd))));
	%    lcid(i,j)=l;
	%
	%    du(1,1)=m.dudxu(i,j);
	%    du(1,2)=m.dudyu(i,j);
	%    du(2,1)=m.dvdxu(i,j);
	%    du(2,2)=m.dvdyu(i,j);
	%    l=max(abs(imag(eig(du))));
	%    lciu(i,j)=l;
	%
	%
	%    end
	%end
	%
	%m.lcid=lcid;
	%m.lciu=lciu;
end
