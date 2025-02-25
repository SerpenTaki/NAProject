    clc;
clear all;
close all;
    
TestCases = {
    struct('lambda', [4,4], 'lO', 3.5, 'toll', 1e-2),%1
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
    it = 4;
    maxit = 8;
    
    [l, m, flag] = multialg(A, lO, toll, it, maxit);
    
    if flag == 1
        fprintf('Newton convergente.\nAutovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
        l1 = round(l);
        toll = 1e-6;
        k = multigeo(A1, l1, toll);
        fprintf('Molteplicità geometrica di %f: %d\n', l1, k);

    else
        fprintf('Metodo non convergente.\nUltimo autovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
    end
    
    linea = 50;
    fprintf('%s\n', repmat('-', 1, linea));

end