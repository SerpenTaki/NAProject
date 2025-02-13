clc;
clear all;
close all;

% Definizione dei test
TestCases = {
    struct('lambda', [1,1,1,4,1,1,5,1,1,4,5,4], 'lO', 3.9, 'toll', 1e-5),%1
    struct('lambda', [1, 1.01, 1.02, 5, 5.01, 5.02, 10, 10.001, 10.002, 10.003], 'lO', 10.1, 'toll', 1e-9),%2
    struct('lambda', [1, 1.01, 1.02, 5, 5.01, 5.02, 10, 10.001, 10.002, 10.003], 'lO', 1.0125, 'toll', 1e-9),%3 % non coverge
    struct('lambda', [1,1,1,4,5,6], 'lO', 3.5, 'toll', 1e-5),%4 % non converge
    struct('lambda', [1,1,1,4,5,6], 'lO', 3.9, 'toll', 1e-5),%5
    struct('lambda', [1,1,1,4,5,6,4], 'lO', 3.9, 'toll', 1e-5),%6
    struct('lambda', [1,1,1,4,5,6,4,4], 'lO', 3.5, 'toll', 1e-9),%7
    struct('lambda', [1,1,1,4,5,6,4,4,4], 'lO', 3.7, 'toll', 1e-9),%8 % converge con molteplicità sbagliata
    struct('lambda', [1,1,1,4,5,6,4,4,4,4], 'lO', 3.3, 'toll', 1e-9),%9 % converge con molteplicità sbagliata
    struct('lambda', [1,1,1,4,5,6,4,4,4], 'lO', 3.7, 'toll', 1e-9),%10
    struct('lambda', [1,1,1,4,5,6,4,4,4,4], 'lO', 3.5, 'toll', 1e-9),%11
    struct('lambda', [1,1,1,4,5,6,4,4,4,4,4], 'lO', 3.4, 'toll', 1e-9), %12 % Converge con molteplicità sbagliata
    struct('lambda', [1,1,1,4,5,6,4,4,4,4,4], 'lO', 3.5, 'toll', 1e-5), % 13
    struct('lambda', [1,1,1,4,5,6,4,4,4,4,4,4], 'lO', 3.4, 'toll', 1e-5), %14
    struct('lambda', [1,1,1,4,5,6,4,4,4,4,4,4], 'lO', 3.9, 'toll', 1e-5), %15
};

for testID = 1:length(TestCases)
    fprintf('\nEseguendo Test %d...\n', testID);
    
    lambda = TestCases{testID}.lambda;
    fprintf('[%s]\n', sprintf('%g ', lambda));

    lO = TestCases{testID}.lO;
    fprintf('lO = %d\n', TestCases{testID}.lO);

    toll = TestCases{testID}.toll;
    fprintf('toll = %d\n', TestCases{testID}.toll);

    % Creazione della matrice Jacobiana
    J = diag(lambda);
    for i = 2:length(lambda)
        if J(i,i) == J(i-1,i-1)
            J(i-1,i) = 1;
        end
    end
    
    n = length(lambda);
    Q = orth(randn(n)); % genera una matrice ortogonale Q di dim nxn
    A = Q' * J * Q;
    disp(J);
    
    % Parametri per il metodo di Newton
    it = 8;
    maxit = 15;
    
    % Chiamata al metodo multialg che restituisce anche il vettore degli step
    [l, m, flag] = multialg(A, lO, toll, it, maxit);
    
    % Visualizzazione dei risultati
    if flag == 1
        fprintf('Newton convergente.\nAutovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
    else
        fprintf('Metodo non convergente.\nUltimo autovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
    end
    
    % Calcolo della molteplicità geometrica
    disp(l);
    l1 = round(l); %limite della geometrica dovuta all'arrotondamento
    k = multigeo(A, l1, toll);
    fprintf('Molteplicità geometrica di %f: %d\n', l1, k);

    linea = 50; % Lunghezza della linea
    fprintf('%s\n', repmat('-', 1, linea));
end
