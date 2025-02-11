function [l, m, flag] = multialg(A, lO, toll, it, maxit)
% MULTIALG Calcola un autovalore di A e la sua molteplicità algebrica.
%
%   [l, m, flag] = multialg(A, lO, toll, it, maxit)
%
% PRE:
%   A    -> Matrice
%   lO   -> Punto di partenza per il metodo di Newton
%   toll -> Tolleranza per il criterio di arresto
%   it   -> Numero minimo di iterazioni per il Newton classico (>=2)
%   maxit-> Numero massimo di iterazioni per il Newton modificato
%
% POST:
%   l    -> Autovalore calcolato
%   m    -> Molteplicità algebrica stimata (naturale positivo)
%   flag -> Flag: 1 se il metodo converge, 0 altrimenti
%
% Il metodo lavora in due fasi:
% 1. Newton classico per ottenere un'approssimazione iniziale.
% 2. Se l'approssimazione non soddisfa il criterio, stima la molteplicità m 
%    mediante il rapporto |f(z - f/g)|/|f(z)| e prosegue con Newton modificato.

% Inizializzazione
z = lO;
iter_values = [];

% --- Fase 1: Newton classico ---
% Forza un numero minimo di iterazioni (almeno 10) prima di verificare il criterio
min_iter = max(10, it);
for i = 1:it
    [f, g] = myobjective(z, A);
    iter_values = [iter_values, z];
    if (i >= min_iter) && (abs(f) < toll)
        % Se il criterio è soddisfatto, consideriamo la radice semplice
        l = z;
        m = 1;
        flag = 1;
        testGrafico(iter_values);
        return;
    end
    z = z - f/g;
end

% --- Stima della molteplicità ---
% Calcoliamo f(z) e, dopo un passo di Newton classico, f(z - f/g)
[f1, g1] = myobjective(z, A);
z_temp = z - f1/g1;
[f2, ~] = myobjective(z_temp, A);

% Calcolo del rapporto osservato
ratio = abs(f2) / abs(f1);

% Se il rapporto è troppo elevato (vicino a 1), probabilmente la radice è semplice.
if ratio > 0.5
    l = z;
    m = 1;
    flag = 1;
    testGrafico(iter_values);
    return;
end

% Nuovo criterio di stima:
% Cerchiamo m in un intervallo più ampio, ad esempio da 2 a 50, 
% che minimizza |(1-1/m)^m - ratio|
m_candidates = 2:50;
differences = abs((1 - 1./m_candidates).^m_candidates - ratio);
[~, idx] = min(differences);
m_est = m_candidates(idx);
m = m_est;

% --- Fase 2: Newton modificato ---
% Puoi scegliere se ripartire da lO oppure usare l'approssimazione ottenuta.
% In questo esempio usiamo l'approssimazione ottenuta.
calls = 0;
for i = 1:maxit
    [f, g] = myobjective(z, A);
    calls = calls + 1;
    iter_values = [iter_values, z];
    if abs(f) < toll
        l = z;
        flag = 1;
        testGrafico(iter_values);
        return;
    end
    z = z - m*(f/g);
    if calls >= 10 * maxit
        break;
    end
end

% Se non si raggiunge la convergenza, incrementiamo m e ripartiamo dal punto iniziale
while calls < 10 * maxit
    m = m + 1;  % Proviamo a incrementare la molteplicità
    z = lO;
    for i = 1:maxit
        [f, g] = myobjective(z, A);
        calls = calls + 1;
        iter_values = [iter_values, z];
        if abs(f) < toll
            l = z;
            flag = 1;
            testGrafico(iter_values);
            return;
        end
        z = z - m*(f/g);
        if calls >= 10 * maxit
            break;
        end
    end
end

% Se non convergiamo, restituiamo l'ultima approssimazione e flag 0.
l = z;
flag = 0;
testGrafico(iter_values);
end
