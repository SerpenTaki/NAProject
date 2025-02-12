function [l, m, flag] = multialg(A, lO, toll, it, maxit)
    z = lO;
    iter_values = [];
    steps = [];
    
    min_iter = max(10, it);
    for i = 1:it
        [f, g] = myobjective(z, A); 
        iter_values = [iter_values, z];

        s = g;                    
        if i == 1
            last_step = s;
        else
            penultimate_step = last_step;
            last_step = s;
        end
        
        if (i >= min_iter) && (abs(f) < toll)
            l = z;
            m = 1;
            flag = 1;
            testGrafico(iter_values);
            return;
        end
        z = z - s;
    end
    
    fprintf('Penultimo step: %e\n', penultimate_step);
    fprintf('Ultimo step: %e\n', last_step);
    maxit = 50;
    m = abs(penultimate_step / (penultimate_step - last_step));
    m = round(m);
    l=z;
    flag=0;
    
    totalCalls = 0;

    m_modified = m;
    converged = false;
    while totalCalls < 10 * maxit
        z = lO;
        for j = 1:maxit
            [f, g] = myobjective(z, A);
            totalCalls = totalCalls + 1;
            iter_values = [iter_values, z];
            s = m_modified * g;  
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
            m_modified = m_modified + 1;
        end
    end

    testGrafico(iter_values);
    end