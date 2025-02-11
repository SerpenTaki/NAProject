clc;
clear all;
close all;

% Definizione di una matrice di test
% Ad esempio, consideriamo una matrice 3x3 con autovalore ripetuto (defect)
%A = [4 1 0; 0 4 1; 0 0 4];
lambda = [1,1,1,4,4,4,5,6];
J = creaJacob(lambda);
n = length(lambda);
disp(J);
Q = orth(randn(n)); %genera una matrice ortogonale Q di dim nxn
A = Q' * J * Q;
% Punto iniziale per il metodo di Newton
lO = 3;  % inizialmente vicino all'autovalore 4

% Parametri per il metodo
toll = 1e-9;
it = 2;
maxit = 50;

% Chiamata al metodo multialg che restituisce anche il vettore degli step
[l, m, flag, steps] = multialg(A, lO, toll, it, maxit);

% Visualizzazione dei risultati
if flag == 1
    fprintf('Newton convergente.\nAutovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
else
    fprintf('Metodo non convergente.\nUltimo autovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
end

% Stampa di tutti gli s_k (Newton steps)
fprintf('\nTutti i passi di Newton (s_k):\n');
disp(steps);

l1 = 4;
%multigo
