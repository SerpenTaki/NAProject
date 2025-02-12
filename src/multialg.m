function [l, m, flag] = multialg(A, lO, toll, it, maxit)
% MULTIALG Calcola un autovalore di A e la sua molteplicità algebrica.
%
%   [l, m, flag, steps] = multialg(A, lO, toll, it, maxit)
%
% PRE:
%   A    -> Matrice
%   lO   -> Punto di partenza per il metodo di Newton
%   toll -> Tolleranza per il criterio di arresto
%   it   -> Numero minimo di iterazioni per il Newton classico (>=2)
%   maxit-> Numero massimo di iterazioni per il Newton modificato
%
% POST:
%   l     -> Autovalore calcolato
%   m     -> Molteplicità algebrica stimata (naturale positivo)
%   flag  -> Flag: 1 se il metodo converge, 0 altrimenti
%   steps -> Vettore contenente tutti gli s_k passi di Newton calcolati
%
% Il metodo lavora in due fasi:
% 1. Newton classico per ottenere un'approssimazione iniziale.
% 2. Se l'approssimazione non soddisfa il criterio, stima la molteplicità m 
%    mediante il rapporto |f(z - g)|/|f(z)| (con g = f/f') e prosegue con
%    Newton modificato.

% Inizializzazione
z = lO;
iter_values = [];
steps = []; % vettore per memorizzare tutti gli step di Newton

% --- Fase 1: Newton classico ---
% Forza un numero minimo di iterazioni (almeno 10) prima di verificare il criterio
min_iter = max(10, it);
for i = 1:it
    [f, g] = myobjective(z, A);  % qui g = f(z)/f'(z)
    iter_values = [iter_values, z];
    
    % Calcolo del passo di Newton; in questo esempio usiamo g direttamente
    s = g;                    
    %steps = [steps, s];        % Memorizza lo step corrente
    
    % Aggiornamento delle variabili per gli ultimi due step
    if i == 1
        last_step = s;
    else
        penultimate_step = last_step;
        last_step = s;
    end
    
    if (i >= min_iter) && (abs(f) < toll)
        % Se il criterio di arresto è soddisfatto, usciamo: la radice è considerata semplice.
        l = z;
        m = 1;
        flag = 1;
        testGrafico(iter_values);
        return;
    end
    % Aggiornamento di Newton: z = z - s
    z = z - s;
end

% Al termine del ciclo, le variabili 'penultimate_step' e 'last_step'
% contengono rispettivamente il penultimo e l'ultimo passo calcolato.
fprintf('Penultimo step: %e\n', penultimate_step);
fprintf('Ultimo step: %e\n', last_step);
maxit = 50;
m = abs(penultimate_step / (penultimate_step - last_step));
m = round(m);
l=z;
flag=0;

% Stima di m tramite Newton (risolviamo (1-1/m)^m = ratio)

% --- Punto 3: Newton modificato ---
totalCalls = 0;
%converged = false;
% Iniziamo con il valore m stimato e ripartiamo da lO
m_modified = m;
converged = false;
while totalCalls < 10 * maxit
    z = lO;
  % riparto da lO per ogni nuovo tentativo con m_modified
    for j = 1:maxit
        [f, g] = myobjective(z, A);
        totalCalls = totalCalls + 1;
        iter_values = [iter_values, z];
        s = m_modified * g;  % passo modificato
        steps = [steps, s];
        if j>=2 && j<maxit-1
            a = steps(j-1);
            b = steps(j);
            if (a-b) < toll
                converged = true;
                flag=1;
                l = z;
                m = m_modified;
                return;
            end
        end
        
        z = z - s;
        
        if totalCalls >= 10 * maxit
            break;
        end
        if s == 0;
            return;
        end
    end
    if converged == false
        % Se non converge entro maxit passi con l'attuale m, incrementa m di 1 e ripeti.
        m_modified = m_modified + 1;
    end
end
  % restituisco il valore finale di m
testGrafico(iter_values);
end