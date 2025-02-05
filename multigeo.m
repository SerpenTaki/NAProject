function [k] = multigeo(A, l, toll)
    % PRE: A è una matrice quadrata reale o complessa
    %      l è uno scalare complesso (autovalore)
    %      toll è un valore di tolleranza positivo
    % POST: k è la dimensione dell'autospazio associato all'autovalore l

    % Controllo che A sia quadrata
    [n, m] = size(A);
    if n ~= m
        error("La matrice A deve essere quadrata.");
    end

    % Calcolo della LU di (A - l*I)
    [~, U, ~] = lu(A - l * eye(n)); %uso ~ per variabili inutilizzate

    % Conta gli elementi diagonali di U che sono in modulo < toll
    k = sum(abs(diag(U)) < toll);
end
