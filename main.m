clc;
clear all;
close all;

% Dimensione della matrice
n = 5;
    
% Creazione di una matrice test con autovalori noti
lambda = [2, 2, 3, 3, 3];
J = blkdiag([2 1; 0 2], [3 1 0; 0 3 1; 0 0 3]);
Q = orth(randn(n)); %genera una matrice ortogonale Q di dim nxn
A = Q * J * Q';
    
% Autovalore target
l0 = 2;
    
% Parametri
toll = 1e-6;
it = 5;
maxit = 50;
    
% Test multgeo
k = multigeo(A, l0, toll);
fprintf('Molteplicità geometrica di %f: %d\n', l0, k);
    
% Test myobjective
[f, g] = myobjective(l0, A);
fprintf('f(%f) = %f, g(%f) = %f\n', l0, f, l0, g);
    
% Test multalg
[l, m, flag] = multialg(A, l0, toll, it, maxit);
if flag
    fprintf('Autovalore calcolato: %f con molteplicità %d\n', l, m);
else
    fprintf('Errore nel calcolo dell''autovalore.\n');
end