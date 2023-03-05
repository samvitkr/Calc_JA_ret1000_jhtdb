


Nx=2048;
Nz=1536;
Ny=512;
ut = 0.0499;
dnu=1.0006e-3;
%fl=int8(zeros(Ny,Nx,Nz));

%edges=[-3:0.1:0];
%edges=flip(-logspace(-1,0.3010,20),2);
%edges=[-2:0.05:-0.05];
%edges=[0.05:0.05:2];
%load('smoothedges.mat');
%edges=[-0.1:0.01:0.1];
%edges=sme*0.1;
edges=[-inf,-2.4]
ne=length(edges);

lcond=single(zeros(Ny,Nx,Nz));
pcond=lcond;
viscond=lcond;
convcond=lcond;
mb=matfile('bsplinedata.mat');
yp = mb.yv;
for time=2:20
        time
        fl=sprintf("lambda_%03d",time);
        ml=matfile(fl);
       l=ml.lambda2;
%        l=ml.Q*dnu/ut;
        ft=sprintf("Transfer_%03d.mat",time)

        mcv=matfile(ft);
        conv=mcv.convective_x;
        visc=mcv.viscous_x;
        %%

%        for i =1:ne-1
                i=1;
                [counts, ~, bins] = histcounts(l, [edges(i) edges(i+1)]);
%                pcond = pcond + single(bins);
%                lcond = lcond + single(bins.*l);%./p(:,i);
		viscond  = viscond  + single(bins.*visc);
		convcond = convcond + single(bins.*conv);
%        end
end
%pcond=pcond./19.0;
%lcond=lcond./19.0;
viscond=viscond./19.0;
convcond=convcond./19.0;

mae=matfile('conditioned_aver_field_l','Writable',true);
%mae.edgesq=edges;
%mae.pcond=pcond;
%mae.lcond=lcond;
mae.viscond=viscond;
mae.convcond=convcond;
   
