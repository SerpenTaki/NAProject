clc;
clear all;
close all;
    
TestCases = {
    struct('lambda', [1,1,1,1,1,4,4,4,4,4,4,4], 'lO', 3.8, 'toll', 1e-4),%1
    struct('lambda', [1, 1.01, 1.02, 5, 5.01, 5.02, 10, 10.001, 10.002, 10.003], 'lO', 10.1, 'toll', 1e-4),%2
    struct('lambda', [1, 1.01, 1.02, 5, 5.01, 5.02, 10, 10.001, 10.002, 10.003], 'lO', 1.0125, 'toll', 1e-4),%3
    struct('lambda', [1,1,1,4,5,6], 'lO', 3.5, 'toll', 1e-4),%4 
    struct('lambda', [1,1,1,4,5,6], 'lO', 3.9, 'toll', 1e-4),%5
    struct('lambda', [1,1,1,4,5,6,4], 'lO', 3.9, 'toll', 1e-4),%6
    struct('lambda', [1,1,1,4,5,6,4,4], 'lO', 3.5, 'toll', 1e-4),%7
    struct('lambda', [1,1,1,4,5,6,4,4,4], 'lO', 3.7, 'toll', 1e-4),%8 
    struct('lambda', [1,1,1,4,5,6,4,4,4,4], 'lO', 3.3, 'toll', 1e-4),%9 
    struct('lambda', [1,1,1,4,5,6,4,4,4], 'lO', 3.7, 'toll', 1e-4),%10
    struct('lambda', [1,1,1,4,5,6,4,4,4,4], 'lO', 3.5, 'toll', 1e-4),%11
    struct('lambda', [1,1,1,4,5,6,4,4,4,4,4], 'lO', 4.3, 'toll', 1e-4), %12
    struct('lambda', [1,1,1,4,5,6,4,4,4,4,4], 'lO', 3.5, 'toll', 1e-4), % 13
     struct('lambda', [1,1,1,4,5,6,4,4,4,4,4,4], 'lO', 3.4, 'toll', 1e-4), %14
    struct('lambda', [4,4,4,4,4,4,4,4], 'lO', 3, 'toll', 1e-4), %15
};

for testID = 1:length(TestCases)
    fprintf('\nEseguendo Test %d...\n', testID);
    
    lambda = TestCases{testID}.lambda;
    fprintf('[%s]\n', sprintf('%g ', lambda));

    lO = TestCases{testID}.lO;
    fprintf('lO = %d\n', TestCases{testID}.lO);

    toll = TestCases{testID}.toll;
    fprintf('toll = %d\n', TestCases{testID}.toll);

    J = diag(lambda);
    for i = 2:length(lambda)
        if J(i,i) == J(i-1,i-1)
            J(i-1,i) = 1;
        end
    end
    
    n = length(lambda);
    Q = orth(randn(n));
    A = Q' * J * Q;
    A1 = A;
    disp(J);
    it = 5;
    maxit = 15;
    
    [l, m, flag] = multialg(A, lO, toll, it, maxit);
    
    if flag == 1
        fprintf('Newton convergente.\nAutovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
        l1 = round(l);
        toll = 1e-6;
        k = multigeo(A1, l1, toll);
        fprintf('Molteplicità geometrica di %f: %d\n', l1, k);

    else
        fprintf('Metodo non convergente\n');
    end
    
    linea = 50;
    fprintf('%s\n', repmat('-', 1, linea));
end
close all;