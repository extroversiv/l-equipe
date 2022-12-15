function [post, cplg, row_length] = adjacencymatrix2postsynapticvector(A)

    %% Calculates from the given adjacency matrix A the postsynaptic 
    %% vector, the coupling strengths and the length vector for each row.
    %% The adjacency matrix has the form A(N,N) where the fist index 
    %% denotes the presynaptic neuron (rows) and the second the post-
    %% synaptic neurons (columns).
    
    N = size(A, 1);
    
    row_length = zeros(1, N);
    post = zeros(1, sum(row_length));
    cplg = zeros(1, sum(row_length));
    
    ind = 0;
    tmp = 1:N;

    for n = 1:N        
        connected = (A(n, :) ~= 0);
        connected(n) = 0;   %delete autapses 
        
        row_length(n) = sum(connected);
        post((ind+1):(ind+row_length(n))) = tmp(connected);
        cplg((ind+1):(ind+row_length(n))) = A(n, connected);
        
        ind = ind + row_length(n);
    end

    
    % If the coupling is the same for all neurons, than cplg can be reduced
    % to a single number (homogeneous network)
    
    testUniqueCplg = unique(cplg);
    
    if length(testUniqueCplg) == 1
        cplg = testUniqueCplg;
    end

end
