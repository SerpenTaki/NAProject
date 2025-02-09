function [l,m,flag] = multialg(A,lO,toll,it,maxit)
    
    
    
    
    
    
    
    
    
    
    
    
    

    
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
    
    
    
    

    [f1, g1] = myobjective(z, A);
    [f2, g2] = myobjective(z - f1 / g1, A);

    m = round(log(abs(f2 / f1)) / log(abs(g2 / g1))); 
    
    

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
    
    
    while calls < 10 * maxit
        m = m + 1; 
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
    
   
    
    l = z;
    flag = 0;
    testGrafico(iter_values);
end