% creaJacob - Crea una matrice Jacobiana da un vettore
% PRE:
%    v - Un vettore di valori numerici
% POST:
%    J - Una matrice Jacobiana creata dal vettore di input


function J = creaJacob(v)

    J = zeros(size(v,1));
    J = J + diag(v);

    for i = 2 : length(v)
        if(J(i,i) == J(i-1,i-1))
            J(i-1,i) = 1;
        end
    end

end