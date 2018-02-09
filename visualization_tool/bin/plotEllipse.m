function plotEllipse(MU, SIGMA, label)
    
    line_width = 1.75;
    elpt = ellipsedata(SIGMA, MU, 100, (0:3:3));
    %elpt = ellipsedata(SIGMA, MU, 100, (0:2:2));
    if ( label==+1 )
        for j=1:2:size(elpt,2)
            plot(elpt(:,j),elpt(:,j+1), 'g', 'LineWidth', line_width);
        end
    else
        for j=1:2:size(elpt,2)
            plot(elpt(:,j),elpt(:,j+1), 'r', 'LineWidth', line_width);
        end
    end
    
end