function [k] = multigeo(A, toll)
    % PRE: A is a square real or complex matrix, toll is a positive real scalar
    % POST: k is the rank of the matrix A

    if size(A, 1) ~= size(A, 2)
        error("Matrix must be square");
    end

    
    [L, U, P] = lu(A);

    rank = 0;

    for i = 1:size(U, 1)
        if abs(U(i, i)) > toll
            rank = rank + 1;
        end
    end

    k = rank;
end