clc;
clear all;
close all;


lambda = input('Inserisci i valori di lambda come un vettore (es. [1,1,1,4,5,6,4,4]): ');
J = creaJacob(lambda);
n = length(lambda);
disp(J);
Q = orth(randn(n));
A = Q' * J * Q;

lO = input('Inserisci il valore iniziale lO (es. 3.53): '); 


toll = input('Inserisci il valore di tolleranza toll (es. 1e-5): ');
it = 2;
maxit = 15;

[l, m, flag] = multialg(A, lO, toll, it, maxit)

if flag == 1
    fprintf('Newton convergente.\nAutovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
else
    fprintf('Metodo non convergente.\nUltimo autovalore calcolato: %f\nMolteplicità algebrica stimata: %d\n', l, m);
end


fprintf('\nTutti i passi di Newton (s_k):\n');


l1 = 4;
k = multigeo(A, l1, toll);
fprintf('Molteplicità geometrica di %f: %d\n', l1, k);