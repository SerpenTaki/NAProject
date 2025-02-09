function [f, g] = myobjective(z, A)
 
 
    
    
    n = size(A, 1);
    B = A - z * eye(n);
    
    
    [L, U, P] = lu(B);
    
    
    P_vec = 1:n; 
    
    for i = 1:n
        P_vec(i) = find(P(i, :) == 1);
    end
    
    s = sum(P_vec ~= 1:n); 
    det_P = (-1)^s; 
    
    
    det_B = det_P * prod(diag(U));
    
    
    f = det_B;
    
    
    B_inv = U \ (L \ P); 
    
    
    g = -trace(B_inv); 
end