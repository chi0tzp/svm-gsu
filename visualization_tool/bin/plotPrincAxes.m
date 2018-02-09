function plotPrincAxes(MU, SIGMA, label, var_frac)
    
    line_width = 1.75;
    elpt = ellipsedata(SIGMA, MU, 100, (0:3:3));
    if (label==+1)
        for j=1:2:size(elpt,2)
            plot(elpt(:,j), elpt(:,j+1), 'g', 'LineWidth', line_width);
        end
    else
        for j=1:2:size(elpt,2)
            plot(elpt(:,j),elpt(:,j+1), 'r', 'LineWidth', line_width);
        end
    end
    
    l = 0.99;
    k = -log((1.00-l)^2);
    a = sqrt(k*SIGMA(1,1));
    b = sqrt(k*SIGMA(2,2));

    % ===
    x = MU(1)-a : 0.01 : MU(1)+a;
    y = MU(2) * ones(size(x));
    if ( label==+1 )
        plot(x, y, 'g-', 'LineWidth', line_width);
    else
        plot(x, y, 'r-', 'LineWidth', line_width);
    end
    
    % ===
    y = MU(2)-b : 0.01 : MU(2)+b;
    x = MU(1) * ones(size(y));
    
    if ( label==+1 )
        plot(x, y, 'g-', 'LineWidth', line_width);
    else
        plot(x, y, 'r-', 'LineWidth', line_width);
    end
    
    % =======
    eigenvalues = eig(SIGMA);
    tot_sum = sum(eigenvalues);
    for i=1:length(eigenvalues)
        par_sum = sum(eigenvalues(1:i));
        p = par_sum/tot_sum;
        if p > var_frac
           break; 
        end
    end

    
    
    
end