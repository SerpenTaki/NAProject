clc;
close all;
clear all;

% Inizializzo matrici sia a livello di dimensione che per contenuto
A = eye(2) * 2;
B = eye(3) * 3;
dimA = size(A,1);
dimB = size(B,1);
n = dimA + dimB;

C = zeros(n, n);  
V = zeros(n, 1);  

for i = 1:dimA %riempie il vettore V
    V(i) = A(i,i);
end
for i = 1:dimB
    V(dimA + i) = B(i,i);
end

%Costruisco Jacobi
for i = 1:n
    C(i,i) = V(i);
end

disp(C);

for i = 2:n
    if C(i,i) == C(i-1,i-1)
        C(i-1,i) = 1;
    end 
end

disp(C);