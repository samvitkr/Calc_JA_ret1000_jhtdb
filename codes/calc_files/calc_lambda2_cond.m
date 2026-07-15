close all;
clear;

jcset = [47 54 71 105 130];
%quadrants = {'q2', 'q4'};
quadrants = {'vpos', 'vneg'};
for jj = 1:length(jcset)
    jcond = jcset(jj);
    
    for q = 1:2
        quad = quadrants{q};
        fprintf('Calculating %s Lambda2 and Q for jcond = %d...\n', quad, jcond);
        
        % Target the files we just generated
        file_path = sprintf("../../data/conditional_%s_jcond_%03d.mat", quad, jcond);
        m = matfile(file_path, 'Writable', true);
        
        % Load all gradients into RAM
        dudx = m.dudx; dvdx = m.dvdx; dwdx = m.dwdx;
        dudy = m.dudy; dvdy = m.dvdy; dwdy = m.dwdy;
        dudz = m.dudz; dvdz = m.dvdz; dwdz = m.dwdz;
        
        % Load the fluctuating dudy
        dudyp = m.dudyp;
        
        % 1. Calculate TOTAL Lambda2 and Q
        [lam2, Q] = calc_vortex_criteria(dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz);
        m.lambda2 = single(lam2);
        m.Q = single(Q);
        
        % 2. Calculate FLUCTUATING Lambda2 and Q
        % Substitute dudyp for dudy. All other fluctuating gradients 
        % equal the total gradients since mean flow is purely U(y).
        [lam2f, Qf] = calc_vortex_criteria(dudx, dvdx, dwdx, dudyp, dvdy, dwdy, dudz, dvdz, dwdz);
        m.lambda2f = single(lam2f);
        m.Qf = single(Qf);
    end
end

fprintf('All Lambda2 and Q calculations complete!\n');

% =========================================================================
% Vectorized Function to Calculate Q and Lambda_2
% =========================================================================
function [lam2, Q] = calc_vortex_criteria(dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz)
    
    % 1. Compute Strain (S) and Rotation (O) Tensors
    S11 = dudx;
    S22 = dvdy;
    S33 = dwdz;
    S12 = 0.5 .* (dudy + dvdx);
    S13 = 0.5 .* (dudz + dwdx);
    S23 = 0.5 .* (dwdy + dvdz);
    
    O12 = 0.5 .* (dudy - dvdx); % Note: O12 = -O21
    O13 = 0.5 .* (dudz - dwdx); % Note: O13 = -O31
    O23 = 0.5 .* (dwdy - dvdz); % Note: O23 = -O32
    
    % 2. Q Criterion: 0.5 * (||O||^2 - ||S||^2)
    norm_S_sq = S11.^2 + S22.^2 + S33.^2 + 2.*(S12.^2 + S13.^2 + S23.^2);
    norm_O_sq = 2.*(O12.^2 + O13.^2 + O23.^2);
    Q = 0.5 .* (norm_O_sq - norm_S_sq);
    
    % 3. Lambda 2 Criterion
    % Construct Symmetric Matrix A = S^2 + O^2
    A11 = S11.^2 + S12.^2 + S13.^2 - O12.^2 - O13.^2;
    A22 = S12.^2 + S22.^2 + S23.^2 - O12.^2 - O23.^2;
    A33 = S13.^2 + S23.^2 + S33.^2 - O13.^2 - O23.^2;
    
    A12 = S11.*S12 + S12.*S22 + S13.*S23 - O13.*O23;
    A13 = S11.*S13 + S12.*S23 + S13.*S33 + O12.*O23;
    A23 = S12.*S13 + S22.*S23 + S23.*S33 + O12.*O13;
    
    % 4. Vectorized Analytical Eigenvalues for Symmetric 3x3 Matrices
    p1 = A12.^2 + A13.^2 + A23.^2;
    trace_A_div3 = (A11 + A22 + A33) ./ 3;
    
    p2 = (A11 - trace_A_div3).^2 + (A22 - trace_A_div3).^2 + (A33 - trace_A_div3).^2 + 2.*p1;
    p = sqrt(p2 ./ 6);
    
    % Avoid division by zero where matrix is perfectly diagonal
    p(p == 0) = eps;
    
    % Construct intermediate matrix B
    B11 = (A11 - trace_A_div3) ./ p;
    B22 = (A22 - trace_A_div3) ./ p;
    B33 = (A33 - trace_A_div3) ./ p;
    B12 = A12 ./ p;
    B13 = A13 ./ p;
    B23 = A23 ./ p;
    
    % Determinant of B
    r = (B11.*B22.*B33 + 2.*B12.*B13.*B23 - B11.*B23.^2 - B22.*B13.^2 - B33.*B12.^2) ./ 2;
    
    % Clamp r to [-1, 1] to guard against floating-point inaccuracies
    r = max(-1, min(1, r));
    
    % Angle for trigonometric roots
    phi = acos(r) ./ 3;
    
    % The roots (eigenvalues)
    eig1 = trace_A_div3 + 2.*p.*cos(phi);
    eig3 = trace_A_div3 + 2.*p.*cos(phi + 2*pi/3);
    eig2 = 3.*trace_A_div3 - eig1 - eig3; % Intermediate eigenvalue (Lambda 2)
    
    lam2 = eig2;
end
