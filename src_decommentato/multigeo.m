function [k] = multigeo(A, l, toll)
    [n, m] = size(A);
    if n ~= m
        error("La matrice A deve essere quadrata.");
    end

    [~, U, ~] = lu(A - l * eye(n)); 

    rank = 0;
    for i =1:n
        if abs(U(i,i)) > toll
            rank = rank + 1;
        end
    end

    k = n-rank;
end