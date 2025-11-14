jcond=105;
fc=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

fc=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

jcond=130;
fc=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

fc=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

jcond=71;
fc=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

fc=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

jcond=54;
fc=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)

fc=sprintf("../data/conditionaln_jcond_dudysplit_%03d.mat",jcond);
mc=matfile(fc)
size(mc.eventp)
size(mc.eventn)


%time=1;
%fvel = sprintf("../data/velfieldpar_%02d.mat", time);
%    m = memmapfile(fvel, 'Format', 'single', 'Writable', false);
%
%    disp(m);
%    disp(fieldnames(m.Data))
%    fvelgx = sprintf("../data/velgradx_%03d.mat", time);
%    mgx = memmapfile(fvelgx, 'Format', 'single', 'Writable', false);
%
%    fvelgy = sprintf("../data/velgrady_%03d.mat", time);
%    mgy = memmapfile(fvelgy, 'Format', 'single', 'Writable', false);
%
%    fvelgz = sprintf("../data/velgradz_%03d.mat", time);
%    mgz = memmapfile(fvelgz, 'Format', 'single', 'Writable', false);
%
%    ft = sprintf("../data/Transfer_%03d.mat", time);
%    mt = memmapfile(ft, 'Format', 'single', 'Writable', false);
