
jconds = [47, 54, 71, 105, 130];
nj = length(jconds); % Number of layers to process
Ny = 512;

baseDir = '/data/geyink1/skumar67/Calc_JA_ret1000_jhtdb/data';

mp = matfile(fullfile(baseDir, 'mean_profiles.mat'));
mvr = matfile(fullfile(baseDir, 'vort_rms.mat'));

dUdy_full = mp.dUdy;
ozrms_full = mvr.ozrms;
ozrms_full = 0.5 * (ozrms_full + flipud(ozrms_full));
ozrms_full = ozrms_full(1:256);

% 1. Preallocate cell arrays for layer-specific data
edges_cell = cell(nj, 1);
totalcounts_cell = cell(nj, 1);
binwidth_cell = cell(nj, 1);


for k = 1:nj
    j = jconds(k);
    
    % 1. Find the maximum absolute bound and round it
    max_bound = 5 * ozrms_full(j);
    max_edge = ceil(max_bound * 100) / 100; 
    
    % 2. Create a perfectly symmetric edge array around absolute zero
    half_edges = 0 : 1e-3 : max_edge;
    edges_cell{k} = [-fliplr(half_edges(2:end)), half_edges]; 
    
    % Initialize counters as before
    totalcounts_cell{k} = zeros(1, numel(edges_cell{k}) - 1);
    binwidth_cell{k} = edges_cell{k}(4) - edges_cell{k}(3); % Will exactly equal 1e-2
end


% % 2. Precompute edges and initialize counters for EACH jcond
% for k = 1:nj
%     j = jconds(k);
%     
%     e1 = -5 * ozrms_full(j);
%     e2 = 5 * ozrms_full(j);
%     e1r = round(e1 * 100) / 100;
%     e2r = round(e2 * 100) / 100;
%     
%     edges_cell{k} = [e1r:1e-2:e2r];
%     totalcounts_cell{k} = zeros(1, numel(edges_cell{k}) - 1);
%     binwidth_cell{k} = edges_cell{k}(4) - edges_cell{k}(3);
% end

tstart = 1;
tend = 10;
tstep = 1;

tic
% 3. Time is the OUTER loop to minimize disk operations
for time = tstart:tstep:tend
    fprintf('Processing time step: %d\n', time);
    
    fvo = sprintf("/vast/geyink1/skumar67/Ret_1000_data/vort_%03d.mat", time);
    
    % Open the matfile ONCE per time step
    m = matfile(fvo);
    
    % 4. Loop over layers (jconds) on the INSIDE
    for k = 1:nj
        j = jconds(k);
        jc = Ny - j + 1;
        
        % Extract only the two slices needed for this specific layer
        ozb = squeeze(m.omegaz(j, :, :));
        ozt = -squeeze(m.omegaz(jc, :, :));
        oz = [ozb;ozt];

	% --- CORRUPTION CHECK: Skip missing data blocks ---
        zero_fraction = sum(oz == 0, 'all') / numel(oz);
        if zero_fraction > 0.01 % If >1% of data is exact zeros
            fprintf('  -> WARNING: Skipping j=%d at time %d. Data is corrupted (%.1f%% missing).\n', j, time, zero_fraction * 100);
            continue;
        end
        % --------------------------------------------------

        [Na ,Nb]=size(oz);
        oz = reshape(oz,[Na*Nb,1]);
        oz = oz-mean(oz);
        counts = histcounts(oz, edges_cell{k});
        totalcounts_cell{k} = totalcounts_cell{k} + counts;
    
    end
    toc % Time spent per whole time-step file
end
%%
Ozm = -dUdy_full;
%5. Calculate PDFs and save results for each jcond
%close all
%figure
%hold on
for k = 1:5
    j = jconds(k);

    totalcounts = totalcounts_cell{k};
    binwidth = binwidth_cell{k};
    edges = edges_cell{k};
    ne = numel(edges);

    pdfoz = totalcounts ./ (sum(totalcounts) * binwidth);
    binCenters = edges(1:ne-1) + binwidth/2;
%plot( binCenters./Ozm(j),pdfoz,'.')
mu=sum(binCenters.*pdfoz)*binwidth;
var = sum( ((binCenters-mu).^2).*pdfoz)*binwidth;
sig = sqrt(var);
sk = sum( ((binCenters-mu).^3).*pdfoz)*binwidth/sig^3;
    fvpdf = sprintf("pdfozfluc_j_%03d.mat", j);
    m_out = matfile(fullfile(baseDir, fvpdf), 'Writable', true);

    m_out.pdfoz = pdfoz;
    m_out.binCenters = binCenters;
    m_out.edges = edges;
    m_out.ozrms = ozrms_full(j);
    m_out.dUdy = dUdy_full(j);
    m_out.sk = sk;
    m_out.sig = sig;
end
% xline(-1)
% xline(0)
% xline(1)
% hold off
fprintf('Done processing all layers!\n');


%set(gca,'YScale','log')
