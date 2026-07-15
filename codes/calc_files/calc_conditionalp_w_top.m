

%% average_bottom_split_allvars_p.m
% Calculates separate averages for all variables for bottom-half 
% 'positive' events (from calc_conditionalp.m) based on the sign of w.

close all;
clear;

%% 1. Parameters (Must match calc_conditionalp.m)
Nx=2048; Ny=512; Nz=1536;
Lx=8*pi; Lz=3*pi;
jcond=105;
jct = Ny-jcond+1;
xbox=0.7; zbox=0.4;

dx=Lx/Nx; dz=Lz/Nz;
itarget=Nx/2+1; ktarget=Nz/2+1;

wini=round(xbox/dx); wink=round(zbox/dz);
winiav=round(0.5*wini); winkav=round(0.5*wink);
nzav=2*winkav+1; nxav=2*winiav+1;

%% 2. Load the event list (UPDATED TO 'conditionalp')
event_file = sprintf("../data/conditionalp_jcond_1_%03d.mat", jcond);
load(event_file, 'event'); 

%% 3. Initialize Accumulators
% Helper function to create empty arrays
init_arr = @() single(zeros(nzav, nxav, Ny/2));

% Init pos accumulators (w > 0)
un_pos=init_arr(); vn_pos=init_arr(); wn_pos=init_arr();
dudxn_pos=init_arr(); dvdxn_pos=init_arr(); dwdxn_pos=init_arr();
dudyn_pos=init_arr(); dvdyn_pos=init_arr(); dwdyn_pos=init_arr();
dudzn_pos=init_arr(); dvdzn_pos=init_arr(); dwdzn_pos=init_arr();
vozn_pos=init_arr(); woyn_pos=init_arr();
count_pos = 0;

% Init neg accumulators (w < 0)
un_neg=init_arr(); vn_neg=init_arr(); wn_neg=init_arr();
dudxn_neg=init_arr(); dvdxn_neg=init_arr(); dwdxn_neg=init_arr();
dudyn_neg=init_arr(); dvdyn_neg=init_arr(); dwdyn_neg=init_arr();
dudzn_neg=init_arr(); dvdzn_neg=init_arr(); dwdzn_neg=init_arr();
vozn_neg=init_arr(); woyn_neg=init_arr();
count_neg = 0;

%% 4. Iterate and Process Events
unique_times = unique(event(:,4));

for t = 1:length(unique_times)
    time = unique_times(t);
    if(time==3||time==4)
	    fprintf('Skipping corrupted time step %d...\n', time);
	    continue;
    end
    fprintf('Processing Time Step: %d\n', time);
    
    idx_t = (event(:,4) == time);
    events_t = event(idx_t, :);
    
    % Filter for top half
    is_top = (events_t(:,3) == jct);
    events_top = events_t(is_top, :);
    
    if isempty(events_bot), continue; end
    
    % Load all fields perfectly matching original code structure
    fvel=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velfieldpar_%02d.mat",time)
    m=matfile(fvel);
    ufieldb=single(permute(m.ufield(Ny/2+1:end,:,:),[3 2 1]));
    vfieldb=single(permute(m.vfield(Ny/2+1:end,:,:),[3 2 1]));
    wfieldb=single(permute(m.wfield(Ny/2+1:end,:,:),[3 2 1]));
    clear m
    
    fvelgx=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradx_%03d.mat",time)
    mgx=matfile(fvelgx);
    dudxb=single(permute(mgx.dudx(1+Ny/2:end,:,:),[3 2 1]));
    dvdxb=single(permute(mgx.dvdx(1+Ny/2:end,:,:),[3 2 1]));
    dwdxb=single(permute(mgx.dwdx(1+Ny/2:end,:,:),[3 2 1]));
    clear mgx
    
    fvelgy=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgrady_%03d.mat",time)
    mgy=matfile(fvelgy);
    dudyb=single(permute(mgy.dudy(1+Ny/2:end,:,:),[3 2 1]));
    dvdyb=single(permute(mgy.dvdy(1+Ny/2:end,:,:),[3 2 1]));
    dwdyb=single(permute(mgy.dwdy(1+Ny/2:end,:,:),[3 2 1]));
    clear mgy
    
    fvelgz=sprintf("/vast/geyink1/skumar67/Ret_1000_data/velgradz_%03d.mat",time)
    mgz=matfile(fvelgz);
    dudzb=single(permute(mgz.dudz(1+Ny/2:end,:,:),[3 2 1]));
    dvdzb=single(permute(mgz.dvdz(1+Ny/2:end,:,:),[3 2 1]));
    dwdzb=single(permute(mgz.dwdz(1+Ny/2:end,:,:),[3 2 1]));
    clear mgz
    
    ft=sprintf("/vast/geyink1/skumar67/Ret_1000_data/Transfer_%03d.mat",time)
    mt=matfile(ft);
    vozb=single(permute(mt.voz(1+Ny/2:end,:,:),[3 2 1]));
    woyb=single(permute(mt.woy(1+Ny/2:end,:,:),[3 2 1]));
    clear mt
    
    for e = 1:size(events_bot, 1)
        kloc = events_bot(e, 1);
        iloc = events_bot(e, 2);
        
        % Check condition BEFORE shifting
        w_at_event = wfieldb(kloc, iloc, jcond);
        
        kdelta = ktarget - kloc;
        idelta = itarget - iloc;
        
        % Forward shift (In-place memory safe, identical to your script)
        ufieldb=circshift( ufieldb ,[kdelta idelta]);
        vfieldb=circshift( vfieldb ,[kdelta idelta]);
        wfieldb=circshift( wfieldb ,[kdelta idelta]);
        dudxb=circshift( dudxb ,[kdelta idelta]);
        dvdxb=circshift( dvdxb ,[kdelta idelta]);
        dwdxb=circshift( dwdxb ,[kdelta idelta]);
        dudyb=circshift( dudyb ,[kdelta idelta]);
        dvdyb=circshift( dvdyb ,[kdelta idelta]);
        dwdyb=circshift( dwdyb ,[kdelta idelta]);
        dudzb=circshift( dudzb ,[kdelta idelta]);
        dvdzb=circshift( dvdzb ,[kdelta idelta]);
        dwdzb=circshift( dwdzb ,[kdelta idelta]);
        vozb=circshift( vozb ,[kdelta idelta]);
        woyb=circshift( woyb ,[kdelta idelta]);
        
        % Accumulate based on w_at_event sign
        if w_at_event > 0
            un_pos    = un_pos    + ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            vn_pos    = vn_pos    + vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            wn_pos    = wn_pos    + wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dudxn_pos = dudxn_pos + dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dvdxn_pos = dvdxn_pos + dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dwdxn_pos = dwdxn_pos + dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dudyn_pos = dudyn_pos + dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dvdyn_pos = dvdyn_pos + dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dwdyn_pos = dwdyn_pos + dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dudzn_pos = dudzn_pos + dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dvdzn_pos = dvdzn_pos + dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dwdzn_pos = dwdzn_pos + dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            vozn_pos  = vozn_pos  + vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            woyn_pos  = woyn_pos  + woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            count_pos = count_pos + 1;
            
        elseif w_at_event < 0
            un_neg    = un_neg    + ufieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            vn_neg    = vn_neg    + vfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            wn_neg    = wn_neg    + wfieldb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dudxn_neg = dudxn_neg + dudxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dvdxn_neg = dvdxn_neg + dvdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dwdxn_neg = dwdxn_neg + dwdxb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dudyn_neg = dudyn_neg + dudyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dvdyn_neg = dvdyn_neg + dvdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dwdyn_neg = dwdyn_neg + dwdyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dudzn_neg = dudzn_neg + dudzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dvdzn_neg = dvdzn_neg + dvdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            dwdzn_neg = dwdzn_neg + dwdzb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            vozn_neg  = vozn_neg  + vozb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            woyn_neg  = woyn_neg  + woyb(ktarget-winkav:ktarget+winkav,itarget-winiav:itarget+winiav,:);
            count_neg = count_neg + 1;
        end
        
        % Backward shift (undo, safe for memory)
        ufieldb=circshift( ufieldb ,-[kdelta idelta]);
        vfieldb=circshift( vfieldb ,-[kdelta idelta]);
        wfieldb=circshift( wfieldb ,-[kdelta idelta]);
        dudxb=circshift( dudxb ,-[kdelta idelta]);
        dvdxb=circshift( dvdxb ,-[kdelta idelta]);
        dwdxb=circshift( dwdxb ,-[kdelta idelta]);
        dudyb=circshift( dudyb ,-[kdelta idelta]);
        dvdyb=circshift( dvdyb ,-[kdelta idelta]);
        dwdyb=circshift( dwdyb ,-[kdelta idelta]);
        dudzb=circshift( dudzb ,-[kdelta idelta]);
        dvdzb=circshift( dvdzb ,-[kdelta idelta]);
        dwdzb=circshift( dwdzb ,-[kdelta idelta]);
        vozb=circshift( vozb ,-[kdelta idelta]);
        woyb=circshift( woyb ,-[kdelta idelta]);
    end
    clear ufieldb vfieldb wfieldb dudxb dvdxb dwdxb dudyb dvdyb dwdyb dudzb dvdzb dwdzb vozb woyb
end

%% 5. Save the Results Incrementally (like original script)
fprintf('Done processing. Events with w > 0: %d | Events with w < 0: %d\n', count_pos, count_neg);

% UPDATED OUTPUT FILENAME for the positive events (conditionalp)
fc = sprintf("../data/conditionalp_top_split_allvars_jcond_%03d.mat", jcond);
mc = matfile(fc, 'Writable', true);

if count_pos > 0
    mc.u_pos = un_pos./count_pos; mc.v_pos = vn_pos./count_pos; mc.w_pos = wn_pos./count_pos;
    mc.dudx_pos = dudxn_pos./count_pos; mc.dvdx_pos = dvdxn_pos./count_pos; mc.dwdx_pos = dwdxn_pos./count_pos;
    mc.dudy_pos = dudyn_pos./count_pos; mc.dvdy_pos = dvdyn_pos./count_pos; mc.dwdy_pos = dwdyn_pos./count_pos;
    mc.dudz_pos = dudzn_pos./count_pos; mc.dvdz_pos = dvdzn_pos./count_pos; mc.dwdz_pos = dwdzn_pos./count_pos;
    mc.voz_pos = vozn_pos./count_pos; mc.woy_pos = woyn_pos./count_pos;
    mc.count_pos = count_pos;
end

if count_neg > 0
    mc.u_neg = un_neg./count_neg; mc.v_neg = vn_neg./count_neg; mc.w_neg = wn_neg./count_neg;
    mc.dudx_neg = dudxn_neg./count_neg; mc.dvdx_neg = dvdxn_neg./count_neg; mc.dwdx_neg = dwdxn_neg./count_neg;
    mc.dudy_neg = dudyn_neg./count_neg; mc.dvdy_neg = dvdyn_neg./count_neg; mc.dwdy_neg = dwdyn_neg./count_neg;
    mc.dudz_neg = dudzn_neg./count_neg; mc.dvdz_neg = dvdzn_neg./count_neg; mc.dwdz_neg = dwdzn_neg./count_neg;
    mc.voz_neg = vozn_neg./count_neg; mc.woy_neg = woyn_neg./count_neg;
    mc.count_neg = count_neg;
end

disp('Averages successfully saved.');

