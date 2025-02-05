function [k] = multigeo(A, l, toll)
    %PRE: A matrice quadrata reale o complessa, l scalare complesso, toll
    %scalare reale positivo
    if size(A,1) ~= size(A,2)
        error("Matrice non quadrata");
    end
    if det(A) == 0 
        error("Matrice non invertibile");
    end

    [L,U,P] = lu(A);
    rank = 1;

    for i = 1:size(U,1)
        if abs(U(i,i)) < toll
            break;
        end
        if U(i,i) == 0
            rank = rank+1;
        end    
    end

    k = rank;
    
end
    %POST: la funzione restituisce k reale positivo dimensione dell'autospazio