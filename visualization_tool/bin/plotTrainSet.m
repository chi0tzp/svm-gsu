function plotTrainSet(MEAN_VECTORS, COV_MATS, p)
    
    means_line_width  = 3;
    means_marker_size = 18;

    % === Plot covariance matrices ===
    for i=1:length(COV_MATS)
        y = MEAN_VECTORS(i, 3);
        centre = MEAN_VECTORS(i, 1:2);
        Sigma = COV_MATS{i};
        plotEllipse(centre, Sigma, y);
        if (p<1.00)
            plotPrincAxes(centre, Sigma, y, p);
        end
    end    
    % === Plot mean vectors ===
    for i=1:size(MEAN_VECTORS, 1)
        y = MEAN_VECTORS(i,3);
        x1 = MEAN_VECTORS(i,1);
        x2 = MEAN_VECTORS(i,2);
        if ( y==+1 )
            plot( x1, x2, '+', 'LineWidth', means_line_width, ...
                               'MarkerSize', means_marker_size, ...
                               'Color', [0 0.49804 0]);
        else
            plot( x1, x2, 'x', 'LineWidth', means_line_width, ...
                               'MarkerSize', means_marker_size, ...
                               'Color', [0.8 0 0]);
        end
    end
    
end