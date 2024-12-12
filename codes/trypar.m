% parallel_sum_of_squares_no_profile.m
% Script for parallel computation on a cluster node with 24 processors

% Specify the number of workers (processors) to use
numWorkers = 24; % Use 24 processors on the current node

% Start the parallel pool with the specified number of workers
parpool('local', numWorkers);

% Define the large array for computation
N = 1e6;                 % Number of elements
data = rand(1, N);       % Create a random array

% Preallocate result array
result = zeros(1, N);    % Preallocate memory for results

% Parallel loop for element-wise computation
parfor i = 1:N
    result(i) = data(i)^2;  % Compute square of each element in parallel
end

% Sum the results outside the parallel loop
totalSum = sum(result);     % Sum of all squared elements

% Display the result
fprintf('Total Sum of Squares: %.4f\n', totalSum);

% Close the parallel pool after execution
delete(gcp('nocreate'));
