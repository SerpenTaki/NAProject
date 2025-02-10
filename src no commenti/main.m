clc;
clear all;
close all;

lambda = [1,1,2,2,3,1,3,1,1,6]; %matrice > 6 errore

J = creaJacob(lambda);
%lambda = [2, 2, 3, 3, 3];
n = length(lambda);
%J = blkdiag([2 1; 0 2], [3 1 0; 0 3 1; 0 0 3]);
Q = orth(randn(n)); 
A = Q' * J * Q;
    
lO = 2;

toll = 1e-6;
it = 4;
maxit = 50;

k = multigeo(A, lO, toll);
fprintf('Molteplicità geometrica di %f: %d\n', lO, k);

[f, g] = myobjective(lO, A);
fprintf('f(%f) = %f, g(%f) = %f\n', lO, f, lO, g);

[l, m, flag] = multialg(A, lO, toll, it, maxit);
if flag
    fprintf('Autovalore calcolato: %f con molteplicità A %d\n', l, m);
else
    fprintf('Errore nel calcolo dell''autovalore.\n');
end