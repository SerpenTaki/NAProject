clc;
clear all;
close all;
    
% Definizione dei test
TestCases = {
    % struct('lambda', [1,1,1,1,1,4,4,4,4,4,4,4], 'lO', 3.8, 'toll', 1e-4),%1
    % struct('lambda', [1, 1.01, 1.02, 5, 5.01, 5.02, 10, 10.001, 10.002, 10.003], 'lO', 10.1, 'toll', 1e-4),%2
    % struct('lambda', [1, 1.01, 1.02, 5, 5.01, 5.02, 10, 10.001, 10.002, 10.003], 'lO', 1.0125, 'toll', 1e-4),%3
    % struct('lambda', [1,1,1,4,5,6], 'lO', 3.5, 'toll', 1e-4),%4 
    % struct('lambda', [1,1,1,4,5,6], 'lO', 3.9, 'toll', 1e-4),%5
    % struct('lambda', [1,1,1,4,5,6,4], 'lO', 3.9, 'toll', 1e-4),%6
    % struct('lambda', [1,1,1,4,5,6,4,4], 'lO', 3.5, 'toll', 1e-4),%7
    % struct('lambda', [1,1,1,4,5,6,4,4,4], 'lO', 3.7, 'toll', 1e-4),%8 
    % struct('lambda', [1,1,1,4,5,6,4,4,4,4], 'lO', 3.3, 'toll', 1e-4),%9 
    % struct('lambda', [1,1,1,4,5,6,4,4,4], 'lO', 3.7, 'toll', 1e-4),%10
    % struct('lambda', [1,1,1,4,5,6,4,4,4,4], 'lO', 3.5, 'toll', 1e-4),%11
    % struct('lambda', [1,1,1,4,5,6,4,4,4,4,4], 'lO', 4.3, 'toll', 1e-4), %12
    % struct('lambda', [1,1,1,4,5,6,4,4,4,4,4], 'lO', 3.5, 'toll', 1e-4), % 13
     struct('lambda', [1,1,1,4,5,6,4,4,4,4,4,4], 'lO', 3.4, 'toll', 1e-2), %14
    % struct('lambda', [4,4,4,4,4,4,4,4], 'lO', 0., 'toll', 1e-4), %15
};

for testID = 1:length(TestCases)
    fprintf('\nExecuting Test %d...\n', testID);
    
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
    A1 = A;
    disp(J);
    
    % Parametri per il metodo di Newton
    it = 4;
    maxit = 10;
    
    % Chiamata al metodo multialg che restituisce anche il vettore degli step
    [l, m, flag] = multialg(A, lO, toll, it, maxit);
    
    % Visualizzazione dei risultati
    if flag == 1
        fprintf('Newton convergence.\nCalculated eigenvalue: %f\nEstimated algebraic multiplicity: %d\n', l, m);

        % Calcolo della molteplicità geometrica

        l1 = round(l);
        toll = 1e-4;
        k = multigeo(A1, l1, toll);
        fprintf('Geometric multiplicity of %f: %d\n', l1, k);

    else
        fprintf('Non-convergent Method.\n');
    end
    
    linea = 50; % Lunghezza della linea
    fprintf('%s\n', repmat('-', 1, linea));
end
