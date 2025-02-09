clc;
clear all;
close all;


n = 5;
    

lambda = [2, 2, 3, 3, 3];
J = blkdiag([2 1; 0 2], [3 1 0; 0 3 1; 0 0 3]);
Q = orth(randn(n)); 
A = Q' * J * Q;
    

lO = 1;







    

toll = 1e-6;
it = 4;
maxit = 50;
    

k = multigeo(A, lO, toll);
fprintf('Molteplicità geometrica di 
    

[f, g] = myobjective(lO, A);
fprintf('f(
    

[l, m, flag] = multialg(A, lO, toll, it, maxit);
if flag
    fprintf('Autovalore calcolato: 
else
    fprintf('Errore nel calcolo dell''autovalore.\n');
end
