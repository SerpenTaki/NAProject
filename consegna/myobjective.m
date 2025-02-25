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
    
    s = 0;
    for i = 1:n
     if (P(i,i) == 0)
         for j = i:n
             if (P(j,i) == 1)
                 P(i,i) = 1;
                 P(j,i) = 0; 
                 P(j,j) = 1; 
                 P(i,j) = 0; 
                 s = s+1;
             end
         end
     end
     end

    det_P = (-1)^s; 
    det_P = det(P);
    det_B = det_P * prod(diag(U)); 

    f = det_B; 

    B_inv = U \ (L \ P); 
    g = 1/trace(B_inv); 
end