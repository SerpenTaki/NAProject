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

    z = lO;
    flag = 0;
    iter_values = [];
    steps = []; 
    
    for i = 1:it
        [f, g] = myobjective(z, A); 
        iter_values = [iter_values, z];
        if abs(f) < toll
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
        z = z  + g;
    end
    
    m = abs(penultimate_step / (penultimate_step - last_step));
    m = round(m);
    l=z;
    flag=0;
    
    totalCalls = 0;
    z = lO;
    newMax = 10 * maxit;
    last_step = inf;
    penultimate_step = inf;
    while totalCalls < newMax
        for j = 1:maxit
            if(totalCalls > newMax-1)
                break;
            end
            g = m * g;
            s = g;
            z = z + s;
            [f, g] = myobjective(z, A);
  
            
            totalCalls = totalCalls + 1;
            iter_values = [iter_values, z];
            steps = [steps, s];
            penultimate_step = last_step;
            last_step = g;
            diffe = abs(last_step - penultimate_step);
            if abs(g) < toll
                flag=1;
                l = z;
                testGrafico(iter_values,A);
                return;
            end

        end   
            
        if s == 0
            testGrafico(iter_values,A);
            return;
        end
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
    xlabel('Numero di Iterazioni');
    ylabel('Valore di z');
    title('Convergenza del Metodo di Newton');
    grid on;
   
    figure;

    p = poly(A);

    z_min = min(values) - 1;
    z_max = max(values) + 1;

    fplot(@(z) polyval(p, z), [z_min, z_max], 'r', 'LineWidth', 2);
    xlabel('z');
    ylabel('Valore del polinomio');
    title('Polinomio Caratteristico di A');
    grid on;
    hold on;

    zeri = roots(p);

    zeri_real = zeri(imag(zeri) == 0);

    plot(zeri_real, zeros(size(zeri_real)), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    legend('Polinomio', 'Zeri');
    hold off;
end


