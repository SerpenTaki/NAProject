function [l,m,flag] = multialg(A,lO,toll,it,maxit)
    %PRE: 
    %   A -> Matrice
    %   lO -> Punto di partenza del metodo di Newton
    %   toll -> Tolleranza per arresto 
    %   it -> intero maggiore o uguale a 2
    %   maxit -> intero positivo maggiore di it
    %POST: 
    %   l -> Autovalore calcolato
    %   m ->  Naturale positivo, molteplicità
    %   flag -> Flag 1 per successo, 0 per errore
    %Calcola la molteplicità algebrica di un autovalore di A, con metodo di
    %Newtoon calcolando la radice del termine del polinomio caratteristico
    %dove compare l'autovale calcolato

    % Passi iniziali del metodo di Newton
    z = lO;
    iter_values = [];
    for i = 1:it
        [f, g] = myobjective(z, A);
        iter_values = [iter_values, z];
        if abs(f) < toll
            l = z;
            m = 1;
            flag = 1;
            testGrafico(iter_values);
            return;
        end
        z = z - f / g;
    end
    
    %Stima iniziale della molteplicità
    %velocita di convergenza di Newtoon radici multiple =>
    % |en+1| ~= C|en|^m-1/m C constante dipende dalla matrice

    [f1, g1] = myobjective(z, A);
    [f2, g2] = myobjective(z - f1 / g1, A);

    m = round(log(abs(f2 / f1)) / log(abs(g2 / g1))); 
    
    % Metodo di Newton modificato

    z = lO;
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
        z = z - m * (f / g); 
        if calls >= 10 * maxit
            break;
        end
    end
    
    % Tentativo di incremento di m se il criterio non è soddisfatto
    while calls < 10 * maxit
        m = m + 1; % (x- auto)^m
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
            z = z - m * (f / g);
            if calls >= 10 * maxit
                break;
            end
        end
    end
    
   
    % Se non si raggiunge la convergenza
    l = z;
    flag = 0;
    testGrafico(iter_values);
end