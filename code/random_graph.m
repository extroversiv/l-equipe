function [post, row_length] = random_graph(K, N)

    A = cell(1, N);
    Atmp = 1:N;
    for n = 1:N
        A{n} = Atmp(rand(1,N) < K/(N-1));
        A{n}(A{n} == n) = []; %delete autapses
    end

    row_length = zeros(1, N);
    for n = 1:N
        row_length(n) = length(A{n});
    end

    post = zeros(1, sum(row_length));
    ind = 0;
    for n = 1:N
        post((ind+1):(ind+row_length(n))) = A{n};
        ind = ind + row_length(n);
    end

end
