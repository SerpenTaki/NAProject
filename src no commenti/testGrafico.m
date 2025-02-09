function [] = testGrafico(values)
    figure;
    semilogy(1:length(values), values, "-o", "LineWidth", 1.5, "MarkerSize", 6);
    xlabel("Numero di Iterazioni");
    ylabel("Valore di z");
    title("Convergenza del Metodo di Newton");
    grid on;

end