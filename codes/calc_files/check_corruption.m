% --- Universal Directory Corruption Scanner ---
baseDir = "/vast/geyink1/skumar67/Ret_1000_data";

% Find all .mat files in the target directory
matFiles = dir(fullfile(baseDir, '*.mat'));

fprintf('==================================================\n');
fprintf(' Starting Automated Corruption Scan\n');
fprintf('==================================================\n\n');

for i = 1:length(matFiles)
    filename = matFiles(i).name;
    filepath = fullfile(baseDir, filename);
    
    try
        % Open the matfile without loading the data into RAM
        m = matfile(filepath);
        
        % Get information about all variables inside the file
        file_vars = whos(m);
        file_corrupted = false;
        
        for v = 1:length(file_vars)
            var_name = file_vars(v).name;
            sz = file_vars(v).size;
            
            % Skip small metadata variables (like grid vectors or scalars)
            if prod(sz) < 10000
                continue; 
            end
            
            % Safely extract a test slice to avoid memory overload
            % We will test j=105, or the maximum first dimension if it's smaller
            idx = min(105, sz(1)); 
            
            if length(sz) == 3
                test_data = squeeze(m.(var_name)(idx, :, :));
            elseif length(sz) == 2
                test_data = squeeze(m.(var_name)(idx, :));
            else
                % Fallback for 1D arrays
                test_data = m.(var_name); 
            end
            
            % Calculate the fraction of exact zeros
            zero_fraction = sum(test_data == 0, 'all') / numel(test_data);
            
            % Flag as corrupted if more than 1% of the slice is exactly zero
            if zero_fraction > 0.01
                fprintf('[X] CORRUPT: %-25s | Var: %-10s | %5.1f%% zeros\n', filename, var_name, zero_fraction * 100);
                file_corrupted = true;
            end
        end
        
        if ~file_corrupted
            fprintf('[OK] Clean:  %-25s\n', filename);
        end
        
    catch ME
        fprintf('[!] ERROR:   Could not process %s. (File may be unreadable)\n', filename);
    end
end

fprintf('\nScan Complete.\n');




