function [f, g] = myobjective(z, A)
    % PRE: 
    % - z è uno scalare.
    % - A è una matrice quadrata di dimensione n x n.
    % POST:
    % - f è il determinante della matrice B = A - z * eye(n).
    % - g è l'inverso della traccia dell'inversa di B.

    n = size(A, 1); 
    B = A - z * eye(n); 
    [L, U, P] = lu(B); 

    P_vec = 1:n; 

    for i = 1:n
        P_vec(i) = find(P(i, :) == 1); 
    end

    s = sum(P_vec ~= 1:n); 
    det_P = (-1)^s; 
    det_P = det(P);
    det_B = det_P * prod(diag(U)); 

    f = det_B; 

    B_inv = U \ (L \ P); 
    g = 1/trace(B_inv); 
end