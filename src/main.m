clc;
clear all;
close all;

% Dimensione della matrice
n = 3;
    
% Creazione di una matrice test con autovalori noti
lambda = [2, 2, 3, 3, 3];
J = blkdiag([2 1; 0 2], [3 1 0; 0 3 1; 0 0 3]);
Q = orth(randn(n)); %genera una matrice ortogonale Q di dim nxn
%A = Q' * J * Q;
    
% Autovalore target
lO = 3;

%test manuale

A = [5,4,2; 0,3,-1; 0,0,3]; %Diagonalizzabile auto= 5, 3 
%A = [4,1,0; 0,4,1; 0,0,4]; % Non diagonalizzabile auto 4

% geo <= alg -> diagonlizzabile 
    
% Parametri
toll = 1e-6;
it = 4;
maxit = 50;
    
% Test multgeo
k = multigeo(A, lO, toll);
fprintf('Molteplicità geometrica di %f: %d\n', lO, k);
    
% Test myobjective
[f, g] = myobjective(lO, A);
fprintf('f(%f) = %f, g(%f) = %f\n', lO, f, lO, g);
    
% Test multalg
[l, m, flag] = multialg(A, lO, toll, it, maxit);
if flag
    fprintf('Autovalore calcolato: %f con molteplicità A %d\n', l, m);
else
    fprintf('Errore nel calcolo dell''autovalore.\n');
end
