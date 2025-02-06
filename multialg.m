function [l,m,flag] = multialg(A,lO,toll,it,maxit)
    % Calcola la molteplicità algebrica di un autovalore di A
    
    % Passi iniziali del metodo di Newton
    z = lO;
    for i = 1:it
        [f, g] = myobjective(z, A);
        if abs(f) < toll
            l = z;
            m = 1;
            flag = 1;
            return;
        end
        z = z - f / g;
    end
    
    % Stima iniziale della molteplicità
    [f1, g1] = myobjective(z, A);
    [f2, g2] = myobjective(z - f1 / g1, A);
    m = round(log(abs(f2 / f1)) / log(abs(g2 / g1)));
    
    % Metodo di Newton modificato
    z = lO;
    calls = 0;
    for i = 1:maxit
        [f, g] = myobjective(z, A);
        calls = calls + 1;
        if abs(f) < toll
            l = z;
            flag = 1;
            return;
        end
        z = z - m * (f / g);
        if calls >= 10 * maxit
            break;
        end
    end
    
    % Tentativo di incremento di m se il criterio non è soddisfatto
    while calls < 10 * maxit
        m = m + 1;
        z = l0;
        for i = 1:maxit
            [f, g] = myobjective(z, A);
            calls = calls + 1;
            if abs(f) < toll
                l = z;
                flag = 1;
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
end