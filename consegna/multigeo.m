function [k] = multigeo(A, l, toll)
    % PRE: A è una matrice quadrata reale o complessa
    %      l è uno scalare complesso (autovalore)
    %      toll è un valore di tolleranza positivo
    % POST: k è la dimensione dell'autospazio associato all'autovalore l
    
    %controllo non necessario nel nostro caso, ma per best pra lo facciamo
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