function molteALG = calcola_molteplicita(p, lambda)

    remainder = p;

    molteALG = 0;
    while true
       
        [quotient, remainder] = deconv(remainder, [1, -lambda]);

        if norm(remainder, Inf) < 1e-12  
            molteALG = molteALG + 1;
            remainder = quotient; 
        else
            break;
        end
    end
end
