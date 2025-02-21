time=1;
fvel = sprintf("../data/velfieldpar_%02d.mat", time);
    m = memmapfile(fvel, 'Format', 'single', 'Writable', false);

    disp(m);
    disp(fieldnames(m.Data))
    fvelgx = sprintf("../data/velgradx_%03d.mat", time);
    mgx = memmapfile(fvelgx, 'Format', 'single', 'Writable', false);

    fvelgy = sprintf("../data/velgrady_%03d.mat", time);
    mgy = memmapfile(fvelgy, 'Format', 'single', 'Writable', false);

    fvelgz = sprintf("../data/velgradz_%03d.mat", time);
    mgz = memmapfile(fvelgz, 'Format', 'single', 'Writable', false);

    ft = sprintf("../data/Transfer_%03d.mat", time);
    mt = memmapfile(ft, 'Format', 'single', 'Writable', false);
