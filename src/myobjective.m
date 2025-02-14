function [f, g] = myobjective(z, A)
    % PRE: 
    % - z è uno scalare.
    % - A è una matrice quadrata di dimensione n x n.
    % POST:
    % - f è il determinante della matrice B = A - z * eye(n).
    % - g è l'inverso della traccia dell'inversa di B.

    n = size(A, 1); % Determina la dimensione della matrice A.
    B = A - z * eye(n); % Calcola la matrice B sottraendo z dalla diagonale di A.
    [L, U, P] = lu(B); % Esegue la fattorizzazione LU della matrice B.

    P_vec = 1:n; % Inizializza il vettore di permutazione.

    for i = 1:n
        P_vec(i) = find(P(i, :) == 1); % Trova la posizione degli 1 nella matrice di permutazione P.
    end

    s = sum(P_vec ~= 1:n); % Conta il numero di scambi necessari nella permutazione.
    det_P = (-1)^s; % Calcola il determinante della matrice di permutazione.
    det_P = det(P);
    det_B = det_P * prod(diag(U)); % Calcola il determinante di B.

    f = det_B; % Assegna il determinante di B a f.

    B_inv = U \ (L \ P); % Calcola l'inversa di B usando la fattorizzazione LU.
    g = 1/trace(B_inv); % Calcola l'inverso della traccia dell'inversa di B.
end