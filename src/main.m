clc;
clear all;
close all;

% Definizione di una matrice di test
% Ad esempio, consideriamo una matrice 3x3 con autovalore ripetuto (defect)
%A = [4 1 0; 0 4 1; 0 0 4];
lambda = input('Inserisci i valori di lambda come un vettore (es. [1,1,1,4,5,6,4,4]): ');
J = creaJacob(lambda);
n = length(lambda);
disp(J);
Q = orth(randn(n)); %genera una matrice ortogonale Q di dim nxn
A = Q' * J * Q;
% Punto iniziale per il metodo di Newton
lO = input('Inserisci il valore iniziale lO (es. 3.53): ');  % inizialmente vicino all'autovalore 4

% Parametri per il metodo
toll = input('Inserisci il valore di tolleranza toll (es. 1e-5): ');
it = 3;
maxit = 5;

% Chiamata al metodo multialg che restituisce anche il vettore degli step
[l, m, flag] = multialg(A, lO, toll, it, maxit)

% Visualizzazione dei risultati
if flag == 1
    fprintf('Newton convergente.\nAutovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
else
    fprintf('Metodo non convergente.\nUltimo autovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
end

% Stampa di tutti gli s_k (Newton steps)
fprintf('\nTutti i passi di Newton (s_k):\n');
%disp(steps);

l1 = round(l);
k = multigeo(A, l1, toll);
fprintf('Molteplicità geometrica di %f: %d\n', l1, k);