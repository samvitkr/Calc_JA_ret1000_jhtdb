
close all
clear

Nx=2048;
Ny=512;
Nz=1536;

Lx=8*pi;
Lz=3*pi;
dx=Lx/Nx;
dz=Lz/Nz;

% 1. Define all configurations to process in one batch
configs = [
    130, 0.8, 0.6;
    105, 0.7, 0.4;
    71,  0.6, 0.3;
    54,  0.5, 0.25;
    47,  0.5, 0.2
];
num_layers = size(configs, 1);

tstart=1;
tend=10;
tstep=1;
vmul=1;

% Precompute thresholds
load('../../data/JHTDB_RET1000.mat')
vrms=sqrt(JHTDB_RET1000(:,5))./JHTDB_RET1000(end,2);
clear JHTDB_RET1000

% 2. Preallocate 2D Cell Arrays [Row 1 = v positive, Row 2 = v negative]
un = cell(2, num_layers); vn = cell(2, num_layers); wn = cell(2, num_layers);
upn = cell(2, num_layers); dudypn = cell(2, num_layers);
dudxn = cell(2, num_layers); dvdxn = cell(2, num_layers); dwdxn = cell(2, num_layers);
dudyn = cell(2, num_layers); dvdyn = cell(2, num_layers); dwdyn = cell(2, num_layers);
dudzn = cell(2, num_layers); dvdzn = cell(2, num_layers); dwdzn = cell(2, num_layers);
vozn = cell(2, num_layers); woyn = cell(2, num_layers);

vthreshold = zeros(num_layers, 1);
event_location = cell(2, num_layers);
counter = zeros(2, num_layers);

for k = 1:num_layers
    jcond = configs(k, 1); xbox = configs(k, 2); zbox = configs(k, 3);
    vthreshold(k) = vmul * abs(vrms(jcond));

    wini=round(xbox/dx); wink=round(zbox/dz);
    winiav=round(0.5*wini); winkav=round(0.5*wink);
    nzav=2*winkav+1; nxav=2*winiav+1;

    for e = 1:2
        un{e,k}=single(zeros(nzav,nxav,Ny/2)); vn{e,k}=single(zeros(nzav,nxav,Ny/2)); wn{e,k}=single(zeros(nzav,nxav,Ny/2));
        upn{e,k}=single(zeros(nzav,nxav,Ny/2)); dudypn{e,k}=single(zeros(nzav,nxav,Ny/2));
        dudxn{e,k}=single(zeros(nzav,nxav,Ny/2)); dvdxn{e,k}=single(zeros(nzav,nxav,Ny/2)); dwdxn{e,k}=single(zeros(nzav,nxav,Ny/2));
        dudyn{e,k}=single(zeros(nzav,nxav,Ny/2)); dvdyn{e,k}=single(zeros(nzav,nxav,Ny/2)); dwdyn{e,k}=single(zeros(nzav,nxav,Ny/2));
        dudzn{e,k}=single(zeros(nzav,nxav,Ny/2)); dvdzn{e,k}=single(zeros(nzav,nxav,Ny/2)); dwdzn{e,k}=single(zeros(nzav,nxav,Ny/2));
        vozn{e,k}=single(zeros(nzav,nxav,Ny/2)); woyn{e,k}=single(zeros(nzav,nxav,Ny/2));
        event_location{e,k} = [];
    end
end

s=[Nz Nx];

% Helper function for 180-deg spatial rotation (Flips dim 3/y and dim 1/z)
rot_spat = @(x) flip(flip(x, 3), 1);

% Main Time Loop
for time=tstart:tstep:tend
    if time == 3 || time == 4
        fprintf('Skipping corrupted time step %d...\n', time);
        continue;
    end

    % Load Data Blocks with Top-Half Spatial Rotation & Sign Rules Applied Immediately
    fvel=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat",time);
    m=matfile(fvel);
    ufieldb=single(permute(m.ufield(1:Ny/2,:,:),[3 2 1])); ufieldt=single(rot_spat(permute(m.ufield(Ny/2+1:end,:,:),[3 2 1])));
    vfieldb=single(permute(m.vfield(1:Ny/2,:,:),[3 2 1])); vfieldt=-single(rot_spat(permute(m.vfield(Ny/2+1:end,:,:),[3 2 1])));
    wfieldb=single(permute(m.wfield(1:Ny/2,:,:),[3 2 1])); wfieldt=-single(rot_spat(permute(m.wfield(Ny/2+1:end,:,:),[3 2 1])));
    clear m

    mgx=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradx_%03d.mat",time));
    dudxb=single(permute(mgx.dudx(1:Ny/2,:,:),[3 2 1])); dudxt=single(rot_spat(permute(mgx.dudx(Ny/2+1:end,:,:),[3 2 1])));
    dvdxb=single(permute(mgx.dvdx(1:Ny/2,:,:),[3 2 1])); dvdxt=-single(rot_spat(permute(mgx.dvdx(Ny/2+1:end,:,:),[3 2 1])));
    dwdxb=single(permute(mgx.dwdx(1:Ny/2,:,:),[3 2 1])); dwdxt=-single(rot_spat(permute(mgx.dwdx(Ny/2+1:end,:,:),[3 2 1])));
    clear mgx

    mgy=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgrady_%03d.mat",time));
    dudyb=single(permute(mgy.dudy(1:Ny/2,:,:),[3 2 1])); dudyt=-single(rot_spat(permute(mgy.dudy(Ny/2+1:end,:,:),[3 2 1])));
    dvdyb=single(permute(mgy.dvdy(1:Ny/2,:,:),[3 2 1])); dvdyt=single(rot_spat(permute(mgy.dvdy(Ny/2+1:end,:,:),[3 2 1])));
    dwdyb=single(permute(mgy.dwdy(1:Ny/2,:,:),[3 2 1])); dwdyt=single(rot_spat(permute(mgy.dwdy(Ny/2+1:end,:,:),[3 2 1])));
    clear mgy

    mgz=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradz_%03d.mat",time));
    dudzb=single(permute(mgz.dudz(1:Ny/2,:,:),[3 2 1])); dudzt=-single(rot_spat(permute(mgz.dudz(Ny/2+1:end,:,:),[3 2 1])));
    dvdzb=single(permute(mgz.dvdz(1:Ny/2,:,:),[3 2 1])); dvdzt=single(rot_spat(permute(mgz.dvdz(Ny/2+1:end,:,:),[3 2 1])));
    dwdzb=single(permute(mgz.dwdz(1:Ny/2,:,:),[3 2 1])); dwdzt=single(rot_spat(permute(mgz.dwdz(Ny/2+1:end,:,:),[3 2 1])));
    clear mgz

    mt=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/Transfer_%03d.mat",time));
    vozb=single(permute(mt.voz(1:Ny/2,:,:),[3 2 1])); vozt=single(rot_spat(permute(mt.voz(Ny/2+1:end,:,:),[3 2 1])));
    woyb=single(permute(mt.woy(1:Ny/2,:,:),[3 2 1])); woyt=single(rot_spat(permute(mt.woy(Ny/2+1:end,:,:),[3 2 1])));
    clear mt

    U_mean_b = mean(ufieldb, [1 2]); dudy_mean_b = mean(dudyb, [1 2]);
    U_mean_t = mean(ufieldt, [1 2]); dudy_mean_t = mean(dudyt, [1 2]);

    for k = 1:num_layers
        jcond = configs(k, 1); xbox = configs(k, 2); zbox = configs(k, 3);
        jc = Ny - jcond + 1;

        wini=round(xbox/dx); wink=round(zbox/dz);
        winiav=round(0.5*wini); winkav=round(0.5*wink);

        vb = vfieldb(:,:,jcond); vt = vfieldt(:,:,jcond);

        % Define v-only masks
        vj_pos_b = vb.*(vb > vthreshold(k));
        vj_neg_b = vb.*(vb < -vthreshold(k));
        vj_pos_t = vt.*(vt > vthreshold(k));
        vj_neg_t = vt.*(vt < -vthreshold(k));

        for eType = 1:2
            if eType == 1
                vjc_bot = vj_pos_b; vjc_top = vj_pos_t;
            else
                vjc_bot = vj_neg_b; vjc_top = vj_neg_t; 
            end

            % --- Process Bottom Events ---
            [bot_k, bot_i] = extract_events(-abs(vjc_bot), vthreshold(k), wink, wini, Nz, Nx);

            for e = 1:length(bot_k)
                kloc = bot_k(e); iloc = bot_i(e);
                w_event = wfieldb(kloc, iloc, jcond); % Get w at event center
                
                event_location{eType,k} = [event_location{eType,k}; kloc iloc jcond time];
                counter(eType,k) = counter(eType,k) + 1;

                z_idx = mod((kloc - winkav : kloc + winkav) - 1, Nz) + 1;
                x_idx = mod((iloc - winiav : iloc + winiav) - 1, Nx) + 1;

                % Slice the raw 3D windows
                u_w=ufieldb(z_idx,x_idx,:); v_w=vfieldb(z_idx,x_idx,:); w_w=wfieldb(z_idx,x_idx,:);
                dudx_w=dudxb(z_idx,x_idx,:); dvdx_w=dvdxb(z_idx,x_idx,:); dwdx_w=dwdxb(z_idx,x_idx,:);
                dudy_w=dudyb(z_idx,x_idx,:); dvdy_w=dvdyb(z_idx,x_idx,:); dwdy_w=dwdyb(z_idx,x_idx,:);
                dudz_w=dudzb(z_idx,x_idx,:); dvdz_w=dvdzb(z_idx,x_idx,:); dwdz_w=dwdzb(z_idx,x_idx,:);
                voz_w=vozb(z_idx,x_idx,:); woy_w=woyb(z_idx,x_idx,:);
                up_w = u_w - U_mean_b; dudyp_w = dudy_w - dudy_mean_b;

                % Apply conditional spanwise reflection to force w > 0
                [u_w,v_w,w_w, dudx_w,dvdx_w,dwdx_w, dudy_w,dvdy_w,dwdy_w, dudz_w,dvdz_w,dwdz_w, voz_w,woy_w, up_w,dudyp_w] = ...
                    cond_span_flip(w_event, u_w,v_w,w_w, dudx_w,dvdx_w,dwdx_w, dudy_w,dvdy_w,dwdy_w, dudz_w,dvdz_w,dwdz_w, voz_w,woy_w, up_w,dudyp_w);

                un{eType,k}=un{eType,k}+u_w; vn{eType,k}=vn{eType,k}+v_w; wn{eType,k}=wn{eType,k}+w_w;
                dudxn{eType,k}=dudxn{eType,k}+dudx_w; dvdxn{eType,k}=dvdxn{eType,k}+dvdx_w; dwdxn{eType,k}=dwdxn{eType,k}+dwdx_w;
                dudyn{eType,k}=dudyn{eType,k}+dudy_w; dvdyn{eType,k}=dvdyn{eType,k}+dvdy_w; dwdyn{eType,k}=dwdyn{eType,k}+dwdy_w;
                dudzn{eType,k}=dudzn{eType,k}+dudz_w; dvdzn{eType,k}=dvdzn{eType,k}+dvdz_w; dwdzn{eType,k}=dwdzn{eType,k}+dwdz_w;
                vozn{eType,k}=vozn{eType,k}+voz_w; woyn{eType,k}=woyn{eType,k}+woy_w;
                upn{eType,k}=upn{eType,k}+up_w; dudypn{eType,k}=dudypn{eType,k}+dudyp_w;
            end

            % --- Process Top Events ---
            [top_k, top_i] = extract_events(-abs(vjc_top), vthreshold(k), wink, wini, Nz, Nx);

            for e = 1:length(top_k)
                kloc = top_k(e); iloc = top_i(e);
                kloc_global = Nz - kloc + 1; % Log the un-flipped global Z coordinate
                
                % Check w_event from the ALREADY rotated top-wall reference frame
                w_event = wfieldt(kloc, iloc, jcond); 

                event_location{eType,k} = [event_location{eType,k}; kloc_global iloc jc time];
                counter(eType,k) = counter(eType,k) + 1;

                z_idx = mod((kloc - winkav : kloc + winkav) - 1, Nz) + 1;
                x_idx = mod((iloc - winiav : iloc + winiav) - 1, Nx) + 1;

                % Slice the raw 3D windows
                u_w=ufieldt(z_idx,x_idx,:); v_w=vfieldt(z_idx,x_idx,:); w_w=wfieldt(z_idx,x_idx,:);
                dudx_w=dudxt(z_idx,x_idx,:); dvdx_w=dvdxt(z_idx,x_idx,:); dwdx_w=dwdxt(z_idx,x_idx,:);
                dudy_w=dudyt(z_idx,x_idx,:); dvdy_w=dvdyt(z_idx,x_idx,:); dwdy_w=dwdyt(z_idx,x_idx,:);
                dudz_w=dudzt(z_idx,x_idx,:); dvdz_w=dvdzt(z_idx,x_idx,:); dwdz_w=dwdzt(z_idx,x_idx,:);
                voz_w=vozt(z_idx,x_idx,:); woy_w=woyt(z_idx,x_idx,:);
                up_w = u_w - U_mean_t; dudyp_w = dudy_w - dudy_mean_t;

                % Apply conditional spanwise reflection to force w > 0
                [u_w,v_w,w_w, dudx_w,dvdx_w,dwdx_w, dudy_w,dvdy_w,dwdy_w, dudz_w,dvdz_w,dwdz_w, voz_w,woy_w, up_w,dudyp_w] = ...
                    cond_span_flip(w_event, u_w,v_w,w_w, dudx_w,dvdx_w,dwdx_w, dudy_w,dvdy_w,dwdy_w, dudz_w,dvdz_w,dwdz_w, voz_w,woy_w, up_w,dudyp_w);

                un{eType,k}=un{eType,k}+u_w; vn{eType,k}=vn{eType,k}+v_w; wn{eType,k}=wn{eType,k}+w_w;
                dudxn{eType,k}=dudxn{eType,k}+dudx_w; dvdxn{eType,k}=dvdxn{eType,k}+dvdx_w; dwdxn{eType,k}=dwdxn{eType,k}+dwdx_w;
                dudyn{eType,k}=dudyn{eType,k}+dudy_w; dvdyn{eType,k}=dvdyn{eType,k}+dvdy_w; dwdyn{eType,k}=dwdyn{eType,k}+dwdy_w;
                dudzn{eType,k}=dudzn{eType,k}+dudz_w; dvdzn{eType,k}=dvdzn{eType,k}+dvdz_w; dwdzn{eType,k}=dwdzn{eType,k}+dwdz_w;
                vozn{eType,k}=vozn{eType,k}+voz_w; woyn{eType,k}=woyn{eType,k}+woy_w;
                upn{eType,k}=upn{eType,k}+up_w; dudypn{eType,k}=dudypn{eType,k}+dudyp_w;
            end
        end
    end
    clear ufieldb ufieldt vfieldb vfieldt wfieldb wfieldt dudxb dudxt dvdxb dvdxt dwdxb dwdxt dudyb dudyt dvdyb dvdyt dwdyb dwdyt dudzb dudzt dvdzb dvdzt dwdzb dwdzt vozb vozt woyb woyt
end

% Save Results into safely named files
for k = 1:num_layers
    jcond = configs(k, 1);
    
    % Save Positive v Events (Updrafts)
    mc_pos = matfile(sprintf("../../data/conditional_vposw_jcond_%03d.mat", jcond), 'Writable', true);
    mc_pos.event = event_location{1,k}; mc_pos.u = un{1,k} ./ counter(1,k); mc_pos.v = vn{1,k} ./ counter(1,k); mc_pos.w = wn{1,k} ./ counter(1,k);
    mc_pos.up = upn{1,k} ./ counter(1,k); mc_pos.dudyp = dudypn{1,k} ./ counter(1,k);
    mc_pos.dudx = dudxn{1,k} ./ counter(1,k); mc_pos.dvdx = dvdxn{1,k} ./ counter(1,k); mc_pos.dwdx = dwdxn{1,k} ./ counter(1,k);
    mc_pos.dudy = dudyn{1,k} ./ counter(1,k); mc_pos.dvdy = dvdyn{1,k} ./ counter(1,k); mc_pos.dwdy = dwdyn{1,k} ./ counter(1,k);
    mc_pos.dudz = dudzn{1,k} ./ counter(1,k); mc_pos.dvdz = dvdzn{1,k} ./ counter(1,k); mc_pos.dwdz = dwdzn{1,k} ./ counter(1,k);
    mc_pos.voz = vozn{1,k} ./ counter(1,k); mc_pos.woy = woyn{1,k} ./ counter(1,k);

    % Save Negative v Events (Downdrafts)
    mc_neg = matfile(sprintf("../../data/conditional_vnegw_jcond_%03d.mat", jcond), 'Writable', true);
    mc_neg.event = event_location{2,k}; mc_neg.u = un{2,k} ./ counter(2,k); mc_neg.v = vn{2,k} ./ counter(2,k); mc_neg.w = wn{2,k} ./ counter(2,k);
    mc_neg.up = upn{2,k} ./ counter(2,k); mc_neg.dudyp = dudypn{2,k} ./ counter(2,k);
    mc_neg.dudx = dudxn{2,k} ./ counter(2,k); mc_neg.dvdx = dvdxn{2,k} ./ counter(2,k); mc_neg.dwdx = dwdxn{2,k} ./ counter(2,k);
    mc_neg.dudy = dudyn{2,k} ./ counter(2,k); mc_neg.dvdy = dvdyn{2,k} ./ counter(2,k); mc_neg.dwdy = dwdyn{2,k} ./ counter(2,k);
    mc_neg.dudz = dudzn{2,k} ./ counter(2,k); mc_neg.dvdz = dvdzn{2,k} ./ counter(2,k); mc_neg.dwdz = dwdzn{2,k} ./ counter(2,k);
    mc_neg.voz = vozn{2,k} ./ counter(2,k); mc_neg.woy = woyn{2,k} ./ counter(2,k);
end

% =========================================================================
% LOCAL FUNCTIONS (Must remain at the very bottom of the script)
% =========================================================================

function [k_locs, i_locs] = extract_events(vjc, threshold, wink, wini, Nz, Nx)
    % 1. Find all potential candidates exceeding the threshold
    candidate_idx = find(vjc <= -threshold);
    candidate_vals = vjc(candidate_idx);

    % 2. Sort candidates from strongest (most negative) to weakest
    [~, sort_order] = sort(candidate_vals, 'ascend');
    sorted_idx = candidate_idx(sort_order);

    % 3. Preallocate outputs and set up a tracking mask
    valid_mask = true(Nz, Nx);
    k_locs = zeros(length(sorted_idx), 1);
    i_locs = zeros(length(sorted_idx), 1);
    count = 0;

    % 4. Iterate down the sorted list
    for i = 1:length(sorted_idx)
        idx = sorted_idx(i);

        % If this pixel hasn't been blanked out by a stronger, earlier event
        if valid_mask(idx)
            count = count + 1;
            [k, i_col] = ind2sub([Nz, Nx], idx);

            k_locs(count) = k;
            i_locs(count) = i_col;

            % Blank out the surrounding spatial window
            z_mask = mod((k - wink : k + wink) - 1, Nz) + 1;
            x_mask = mod((i_col - wini : i_col + wini) - 1, Nx) + 1;
            valid_mask(z_mask, x_mask) = false;
        end
    end

    % Trim preallocated arrays to exact size
    k_locs = k_locs(1:count);
    i_locs = i_locs(1:count);
end

function [u, v, w, dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz, voz, woy, up, dudyp] = ...
    cond_span_flip(w_event, u, v, w, dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz, voz, woy, up, dudyp)
    % Applies a spatial flip (dim 1) and necessary tensor sign changes if w < 0
    if w_event < 0
        u = flip(u, 1);
        v = flip(v, 1);
        w = -flip(w, 1);  % Forces conditional average w to be positive
        
        dudx = flip(dudx, 1);
        dvdx = flip(dvdx, 1);
        dwdx = -flip(dwdx, 1);
        
        dudy = flip(dudy, 1);
        dvdy = flip(dvdy, 1);
        dwdy = -flip(dwdy, 1);
        
        dudz = -flip(dudz, 1);
        dvdz = -flip(dvdz, 1);
        dwdz = flip(dwdz, 1);
        
        voz = flip(voz, 1); % Non-linear term (sign unchanged)
        woy = flip(woy, 1); % Non-linear term (sign unchanged)
        
        up = flip(up, 1);
        dudyp = flip(dudyp, 1);
    end
end

