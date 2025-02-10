function molteALG = calcola_molteplicita(p, lambda)

    resto = p;

    molteALG = 0;
    while true
        [quoziente, resto] = deconv(resto, [1, -lambda]);
        if norm(resto, Inf) < 0.9999999
            molteALG = molteALG + 1;
            resto = quoziente; 
        else
            break;
        end
    end
end
