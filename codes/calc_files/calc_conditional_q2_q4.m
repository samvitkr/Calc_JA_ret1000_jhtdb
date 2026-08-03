

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

% Precompute thresholds
load('../../data/JHTDB_RET1000.mat')
uvav=abs((JHTDB_RET1000(:,3))./JHTDB_RET1000(end,2)^2);
clear JHTDB_RET1000

% 2. Preallocate 2D Cell Arrays [Row 1 = Q2, Row 2 = Q4]
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
    vthreshold(k) = abs(10 * uvav(jcond));

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
	time
    if time == 3 || time == 4
        fprintf('Skipping corrupted time step %d...\n', time);
        continue;
    end

    % Load Data Blocks
    fvel=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat",time);
    m=matfile(fvel);
    ufieldb=single(permute(m.ufield(1:Ny/2,:,:),[3 2 1]));% ufieldt=single(rot_spat(permute(m.ufield(Ny/2+1:end,:,:),[3 2 1])));
    vfieldb=single(permute(m.vfield(1:Ny/2,:,:),[3 2 1]));% vfieldt=-single(rot_spat(permute(m.vfield(Ny/2+1:end,:,:),[3 2 1])));
    wfieldb=single(permute(m.wfield(1:Ny/2,:,:),[3 2 1]));% wfieldt=-single(rot_spat(permute(m.wfield(Ny/2+1:end,:,:),[3 2 1])));
    clear m

    mgx=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradx_%03d.mat",time));
    dudxb=single(permute(mgx.dudx(1:Ny/2,:,:),[3 2 1]));% dudxt=single(rot_spat(permute(mgx.dudx(Ny/2+1:end,:,:),[3 2 1])));
    dvdxb=single(permute(mgx.dvdx(1:Ny/2,:,:),[3 2 1]));% dvdxt=-single(rot_spat(permute(mgx.dvdx(Ny/2+1:end,:,:),[3 2 1])));
    dwdxb=single(permute(mgx.dwdx(1:Ny/2,:,:),[3 2 1]));% dwdxt=-single(rot_spat(permute(mgx.dwdx(Ny/2+1:end,:,:),[3 2 1])));
    clear mgx

    mgy=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgrady_%03d.mat",time));
    dudyb=single(permute(mgy.dudy(1:Ny/2,:,:),[3 2 1]));% dudyt=-single(rot_spat(permute(mgy.dudy(Ny/2+1:end,:,:),[3 2 1])));
    dvdyb=single(permute(mgy.dvdy(1:Ny/2,:,:),[3 2 1]));% dvdyt=single(rot_spat(permute(mgy.dvdy(Ny/2+1:end,:,:),[3 2 1])));
    dwdyb=single(permute(mgy.dwdy(1:Ny/2,:,:),[3 2 1]));% dwdyt=single(rot_spat(permute(mgy.dwdy(Ny/2+1:end,:,:),[3 2 1])));
    clear mgy

    mgz=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradz_%03d.mat",time));
    dudzb=single(permute(mgz.dudz(1:Ny/2,:,:),[3 2 1]));% dudzt=-single(rot_spat(permute(mgz.dudz(Ny/2+1:end,:,:),[3 2 1])));
    dvdzb=single(permute(mgz.dvdz(1:Ny/2,:,:),[3 2 1]));% dvdzt=single(rot_spat(permute(mgz.dvdz(Ny/2+1:end,:,:),[3 2 1])));
    dwdzb=single(permute(mgz.dwdz(1:Ny/2,:,:),[3 2 1]));% dwdzt=single(rot_spat(permute(mgz.dwdz(Ny/2+1:end,:,:),[3 2 1])));
    clear mgz

    mt=matfile(sprintf("/vast/geyink1/skumar67/Ret_1000_data/Transfer_%03d.mat",time));
    vozb=single(permute(mt.voz(1:Ny/2,:,:),[3 2 1]));% vozt=single(rot_spat(permute(mt.voz(Ny/2+1:end,:,:),[3 2 1])));
    woyb=single(permute(mt.woy(1:Ny/2,:,:),[3 2 1]));% woyt=single(rot_spat(permute(mt.woy(Ny/2+1:end,:,:),[3 2 1])));
    clear mt

    U_mean_b = mean(ufieldb, [1 2]); dudy_mean_b = mean(dudyb, [1 2]);
%   U_mean_t = mean(ufieldt, [1 2]); dudy_mean_t = mean(dudyt, [1 2]);

    for k = 1:num_layers
        jcond = configs(k, 1); xbox = configs(k, 2); zbox = configs(k, 3);
        jc = Ny - jcond + 1;

        wini=round(xbox/dx); wink=round(zbox/dz);
        winiav=round(0.5*wini); winkav=round(0.5*wink);

        vb = vfieldb(:,:,jcond); %vt = vfieldt(:,:,jcond);
        upb = ufieldb(:,:,jcond) - U_mean_b(1,1,jcond);
       %upt = ufieldt(:,:,jcond) - U_mean_t(1,1,jcond);

        uvb_q2 = upb.*vb.*(vb>0); vj_q2 = uvb_q2.*(uvb_q2 < -vthreshold(k));
        uvb_q4 = upb.*vb.*(vb<0); vj_q4 = uvb_q4.*(uvb_q4 < -vthreshold(k));
       %uvt_q2 = upt.*vt.*(vt>0); vjt_q2 = uvt_q2.*(uvt_q2 < -vthreshold(k));
       %uvt_q4 = upt.*vt.*(vt<0); vjt_q4 = uvt_q4.*(uvt_q4 < -vthreshold(k));

        for eType = 1:2
            if eType == 1, vjc_bot = vj_q2; %vjc_top = vjt_q2;
            else, vjc_bot = vj_q4; %vjc_top = vjt_q4; 
	    end

            % --- Process Bottom Events ---
            [bot_k, bot_i] = extract_events(vjc_bot, vthreshold(k), wink, wini, Nz, Nx);
            
            for e = 1:length(bot_k)
                kloc = bot_k(e); iloc = bot_i(e);
                event_location{eType,k} = [event_location{eType,k}; kloc iloc jcond time];
                counter(eType,k) = counter(eType,k) + 1;
                
                z_idx_av = mod((kloc - winkav : kloc + winkav) - 1, Nz) + 1;
                x_idx_av = mod((iloc - winiav : iloc + winiav) - 1, Nx) + 1;
                
                un{eType,k}     = un{eType,k}     + ufieldb(z_idx_av, x_idx_av, :);
                vn{eType,k}     = vn{eType,k}     + vfieldb(z_idx_av, x_idx_av, :);
                wn{eType,k}     = wn{eType,k}     + wfieldb(z_idx_av, x_idx_av, :);
                dudxn{eType,k}  = dudxn{eType,k}  + dudxb(z_idx_av, x_idx_av, :);
                dvdxn{eType,k}  = dvdxn{eType,k}  + dvdxb(z_idx_av, x_idx_av, :);
                dwdxn{eType,k}  = dwdxn{eType,k}  + dwdxb(z_idx_av, x_idx_av, :);
                dudyn{eType,k}  = dudyn{eType,k}  + dudyb(z_idx_av, x_idx_av, :);
                dvdyn{eType,k}  = dvdyn{eType,k}  + dvdyb(z_idx_av, x_idx_av, :);
                dwdyn{eType,k}  = dwdyn{eType,k}  + dwdyb(z_idx_av, x_idx_av, :);
                dudzn{eType,k}  = dudzn{eType,k}  + dudzb(z_idx_av, x_idx_av, :);
                dvdzn{eType,k}  = dvdzn{eType,k}  + dvdzb(z_idx_av, x_idx_av, :);
                dwdzn{eType,k}  = dwdzn{eType,k}  + dwdzb(z_idx_av, x_idx_av, :);
                vozn{eType,k}   = vozn{eType,k}   + vozb(z_idx_av, x_idx_av, :);
                woyn{eType,k}   = woyn{eType,k}   + woyb(z_idx_av, x_idx_av, :);
                upn{eType,k}    = upn{eType,k}    + (ufieldb(z_idx_av, x_idx_av, :) - U_mean_b);
                dudypn{eType,k} = dudypn{eType,k} + (dudyb(z_idx_av, x_idx_av, :) - dudy_mean_b);
            end
            
           % % --- Process Top Events ---
           % [top_k, top_i] = extract_events(vjc_top, vthreshold(k), wink, wini, Nz, Nx);
           % 
           % for e = 1:length(top_k)
           %     kloc = top_k(e); iloc = top_i(e);
           %     kloc_global = Nz - kloc + 1; % Log the global un-flipped z-coordinate for reference
           %     
           %     event_location{eType,k} = [event_location{eType,k}; kloc_global iloc jc time];
           %     counter(eType,k) = counter(eType,k) + 1;
           %     
           %     z_idx_av = mod((kloc - winkav : kloc + winkav) - 1, Nz) + 1;
           %     x_idx_av = mod((iloc - winiav : iloc + winiav) - 1, Nx) + 1;
           %     
           %     un{eType,k}     = un{eType,k}     + ufieldt(z_idx_av, x_idx_av, :);
           %     vn{eType,k}     = vn{eType,k}     + vfieldt(z_idx_av, x_idx_av, :);
           %     wn{eType,k}     = wn{eType,k}     + wfieldt(z_idx_av, x_idx_av, :);
           %     dudxn{eType,k}  = dudxn{eType,k}  + dudxt(z_idx_av, x_idx_av, :);
           %     dvdxn{eType,k}  = dvdxn{eType,k}  + dvdxt(z_idx_av, x_idx_av, :);
           %     dwdxn{eType,k}  = dwdxn{eType,k}  + dwdxt(z_idx_av, x_idx_av, :);
           %     dudyn{eType,k}  = dudyn{eType,k}  + dudyt(z_idx_av, x_idx_av, :);
           %     dvdyn{eType,k}  = dvdyn{eType,k}  + dvdyt(z_idx_av, x_idx_av, :);
           %     dwdyn{eType,k}  = dwdyn{eType,k}  + dwdyt(z_idx_av, x_idx_av, :);
           %     dudzn{eType,k}  = dudzn{eType,k}  + dudzt(z_idx_av, x_idx_av, :);
           %     dvdzn{eType,k}  = dvdzn{eType,k}  + dvdzt(z_idx_av, x_idx_av, :);
           %     dwdzn{eType,k}  = dwdzn{eType,k}  + dwdzt(z_idx_av, x_idx_av, :);
           %     vozn{eType,k}   = vozn{eType,k}   + vozt(z_idx_av, x_idx_av, :);
           %     woyn{eType,k}   = woyn{eType,k}   + woyt(z_idx_av, x_idx_av, :);
           %     upn{eType,k}    = upn{eType,k}    + (ufieldt(z_idx_av, x_idx_av, :) - U_mean_t);
           %     dudypn{eType,k} = dudypn{eType,k} + (dudyt(z_idx_av, x_idx_av, :) - dudy_mean_t);
           % end
        end
    end
   %clear ufieldb ufieldt vfieldb vfieldt wfieldb wfieldt dudxb dudxt dvdxb dvdxt dwdxb dwdxt dudyb dudyt dvdyb dvdyt dwdyb dwdyt dudzb dudzt dvdzb dvdzt dwdzb dwdzt vozb vozt woyb woyt
end

for k = 1:num_layers
    jcond = configs(k, 1);
    mc2 = matfile(sprintf("../../data/conditionalq2_jcond_%03d.mat", jcond), 'Writable', true);
    mc2.event = event_location{1,k}; mc2.u = un{1,k} ./ counter(1,k); mc2.v = vn{1,k} ./ counter(1,k); mc2.w = wn{1,k} ./ counter(1,k);
    mc2.up = upn{1,k} ./ counter(1,k); mc2.dudyp = dudypn{1,k} ./ counter(1,k);
    mc2.dudx = dudxn{1,k} ./ counter(1,k); mc2.dvdx = dvdxn{1,k} ./ counter(1,k); mc2.dwdx = dwdxn{1,k} ./ counter(1,k);
    mc2.dudy = dudyn{1,k} ./ counter(1,k); mc2.dvdy = dvdyn{1,k} ./ counter(1,k); mc2.dwdy = dwdyn{1,k} ./ counter(1,k);
    mc2.dudz = dudzn{1,k} ./ counter(1,k); mc2.dvdz = dvdzn{1,k} ./ counter(1,k); mc2.dwdz = dwdzn{1,k} ./ counter(1,k);
    mc2.voz = vozn{1,k} ./ counter(1,k); mc2.woy = woyn{1,k} ./ counter(1,k);

    mc4 = matfile(sprintf("../../data/conditionalq4_jcond_%03d.mat", jcond), 'Writable', true);
    mc4.event = event_location{2,k}; mc4.u = un{2,k} ./ counter(2,k); mc4.v = vn{2,k} ./ counter(2,k); mc4.w = wn{2,k} ./ counter(2,k);
    mc4.up = upn{2,k} ./ counter(2,k); mc4.dudyp = dudypn{2,k} ./ counter(2,k);
    mc4.dudx = dudxn{2,k} ./ counter(2,k); mc4.dvdx = dvdxn{2,k} ./ counter(2,k); mc4.dwdx = dwdxn{2,k} ./ counter(2,k);
    mc4.dudy = dudyn{2,k} ./ counter(2,k); mc4.dvdy = dvdyn{2,k} ./ counter(2,k); mc4.dwdy = dwdyn{2,k} ./ counter(2,k);
    mc4.dudz = dudzn{2,k} ./ counter(2,k); mc4.dvdz = dvdzn{2,k} ./ counter(2,k); mc4.dwdz = dwdzn{2,k} ./ counter(2,k);
    mc4.voz = vozn{2,k} ./ counter(2,k); mc4.woy = woyn{2,k} ./ counter(2,k);
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



