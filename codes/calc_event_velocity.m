% Load event lists
jcond=71;
fc=sprintf("../data/conditionalp_jcond_dudysplit_%03d.mat",jcond);
%fc=sprintf("../data/test.mat")
mc=matfile(fc,'Writable',true);

S = load(fc, 'eventn', 'eventp');
eventn = S.eventn;   % N1 x 4, [kloc, iloc, jloc, time]
eventp = S.eventp;   % N2 x 4, [kloc, iloc, jloc, time]

N1 = size(eventn, 1);
N2 = size(eventp, 1);

% Output: [u, v] for each event
uv_n = zeros(N1, 2);   % N1 x 2
uv_p = zeros(N2, 2);   % N2 x 2

% Loop over time slices (1 to 10)
for t = 1:10
	t
    % Find which events occur at this time
    idx_n = find(eventn(:,4) == t);
    idx_p = find(eventp(:,4) == t);

    % If no events at this time, skip loading
    if isempty(idx_n) && isempty(idx_p)
        continue;
    end

    % Load u, v for this time only
    fname = sprintf('/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat', t);
    data  = load(fname, 'ufield', 'vfield');
    u = data.ufield;
    v = data.vfield;

    % Size of u and v (assumed same)
    sz = size(u);   % [Nk, Ni, Nj]

    % ------- eventn at time t -------
    if ~isempty(idx_n)
        kn = eventn(idx_n, 1);
        in = eventn(idx_n, 2);
        jn = eventn(idx_n, 3);

        lin_idx_n = sub2ind(sz, kn, in, jn);

        uv_n(idx_n, 1) = u(lin_idx_n);
        uv_n(idx_n, 2) = v(lin_idx_n);
    end

    % ------- eventp at time t -------
    if ~isempty(idx_p)
        kp = eventp(idx_p, 1);
        ip = eventp(idx_p, 2);
        jp = eventp(idx_p, 3);

        lin_idx_p = sub2ind(sz, kp, ip, jp);

        uv_p(idx_p, 1) = u(lin_idx_p);
        uv_p(idx_p, 2) = v(lin_idx_p);
    end
end
mc.uv_n=uv_n;
mc.uv_p=uv_p;

% Now uv_n is N1x2 and uv_p is N2x2 with [u,v] at each (kloc,iloc,jloc,time).
% Optional:
% save('uv_event_vectors.mat', 'uv_n', 'uv_p');

