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
    flag = 0;
    iter_values = [];
    steps = []; % vettore per memorizzare tutti gli step di Newton
    
    % --- Fase 1: Newton classico ---
    % Forza un numero minimo di iterazioni (almeno 10) prima di verificare il criterio
    for i = 1:it
        [f, g] = myobjective(z, A);  % qui g = f(z)/f'(z)
        iter_values = [iter_values, z];

        if abs(f) < toll
            % Se il criterio di arresto è soddisfatto, usciamo: la radice è considerata semplice.
            l = z;
            m = 1;
            flag = 1;
            testGrafico(iter_values,A);
            return;
        end

        if i == 1
            last_step = g;
        else
            penultimate_step = last_step;
            last_step = g;
        end

        % Aggiornamento di Newton: z = z - s
        z = z  + g;
    end
    
    % Al termine del ciclo, le variabili 'penultimate_step' e 'last_step'
    % contengono rispettivamente il penultimo e l'ultimo passo calcolato.
    %fprintf('Penultimo step: %e\n', penultimate_step);
    %fprintf('Ultimo step: %e\n', last_step);
    m = abs(penultimate_step / (penultimate_step - last_step));
    m = round(m);
    l=z;
    flag=0;
    
    % --- Punto 3: Newton modificato ---
    totalCalls = 0;
    z = lO;
    % Iniziamo con il valore m stimato e ripartiamo da lO
    newMax = 10 * maxit;
    last_step = inf;
    penultimate_step = inf;
    while totalCalls < newMax
        %fprintf("\nTOTAL CALLS ->%f\n",totalCalls);
        for j = 1:maxit

            if(totalCalls > newMax-1)
                break;
            end
            
            [f, g] = myobjective(z, A);
            z = z + g;
            s = g;
            totalCalls = totalCalls + 1;
            iter_values = [iter_values, z];
            steps = [steps, s];

            penultimate_step = last_step;
            last_step = g;
            diffe = norm(last_step - penultimate_step);
            
            if diffe < toll
                flag=1;
                l = z;
                testGrafico(iter_values,A);
                %fprintf("Ultimo->%f\npenultimo->%f\ndiff->%f\n",last_step,penultimate_step,diffe);
                return;
            end

        end   
            
        if s == 0
            testGrafico(iter_values,A);
            return;
        end
        % Se non converge entro maxit passi con l'attuale m, incrementa m di 1 e ripeti.

        g = m * g;
        m = m+1;
        totalCalls = totalCalls + 1;
    end
    flag = 0;
    l = z;
    testGrafico(iter_values,A);
    end
    
function [] = testGrafico(values, A)
   
    figure;
    semilogy(1:length(values), values, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('Number of Iterations');
    ylabel('Value of z');
    title('Convergence of Newton Method');
    grid on;
   
    figure;

    p = poly(A);

    z_min = min(values) - 1;
    z_max = max(values) + 1;

    fplot(@(z) polyval(p, z), [z_min, z_max], 'r', 'LineWidth', 2);
    xlabel('z');
    ylabel('Value of the polynomial');
    title('Characteristic Polynomial of A');
    grid on;
    hold on;

    zeri = roots(p);

    zeri_real = zeri(imag(zeri) == 0);

    plot(zeri_real, zeros(size(zeri_real)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    legend('Polynomial', 'Zeros');
    hold off;
end


