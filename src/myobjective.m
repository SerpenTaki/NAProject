function [f, g] = myobjective(z, A)
 %PRE: z numero complesso, A matrice reale o complessa
 %POST: Calcola fA(z) e f'A(z) (gA(z)) utilizzando la fattorizzazione LU con pivoting parziale
    
    % Costruisci B(z) = A - zI
    n = size(A, 1);
    B = A - z * eye(n);
    
    % Fattorizzazione LU con pivoting parziale: PA = LU
    [L, U, P] = lu(B);
    
    % Calcolo di det(P) basato sul numero di scambi di riga
    P_vec = 1:n; % Vettore degli indici
    
    for i = 1:n % Trova l'indice della colonna in cui P(i, :) è uguale a 1
        P_vec(i) = find(P(i, :) == 1);
    end
    
    s = sum(P_vec ~= 1:n); % Numero di scambi
    det_P = (-1)^s; %Determinante positivo se numero di scambi pari, numero dispari determinate negativo
    
    % Calcolo di det(B) = det(A - zI) = det(U) / det(P)
    det_B = det(U)/det_P;

    % Calcolo di fA(z)
    f = det_B;
    
    % Risoluzione B^{-1} utilizzando la fattorizzazione LU
    B_inv = U \ (L \ P); % Equivalente a inv(B)
    
    % Calcolo di gA(z)
    g = 1/trace(B_inv); %la funzione trace calcola la somma degli elementi diagonali di una matrice
    %g = 1/trace(A-z*eye(n));
end