function h = plotSVM(svm_model, kernel_type, gamma, line_style, color)
    
    SVs = svm_model{1};
    rho = svm_model{2};
    l = size(SVs, 1);
    line_width = 2.5;
    
    % +=================================================================+ %
    % |                                                                 | %
    % |                         Linear Kernel                           | %
    % |                                                                 | %
    % +=================================================================+ %
    if strcmp( kernel_type, 'linear' )
        
        % Compute normal vector w (H: w^T*x+b=0)
        w = [0;0];
        for i=1:l
            alpha_i = SVs(i,1);
            w = w + alpha_i * [SVs(i,2); SVs(i,3)];
        end
            
        % Compute the bias term b (H: w^T*x+b=0)
        b = 0;
        for i=1:l
            alpha_i = SVs(i,1); 
            y_i     = sign( alpha_i );
            b = b + ( 1 - y_i*dot(w,[SVs(i,2);SVs(i,3)]) ) / y_i;
        end
        b = b/l;

        % Get x-, y-axis limits
        xLimits = get(gca,'XLim');
        yLimits = get(gca,'YLim');    
            
        % Plot separating hyperplane (line)
        % -- Case: Horizontal line
        if ( w(1) == 0 )
            y = -b/w(2);
            h = line( y*[1 1], xLimits);
        % -- Case: Vertical line
        elseif ( w(2) == 0 )
            x = -b/w(1);
            h = line( x*[1 1], yLimits);
        % -- Otherwise
        else
            line_cmd = sprintf( 'y=-%g*x-%g', w(1)/w(2), b/w(2) );
            h = ezplot(line_cmd, [xLimits yLimits]);
        end
        set(h , 'LineStyle' , line_style, ...
                'Color'     , color, ...
                'LineWidth' , line_width);
            
        % Remove x-, y-axis labels
        title([]);
        xlabel([]);
        ylabel([]);
        set(gca,'xtick',[]);
        set(gca,'xticklabel',[]);
        set(gca,'ytick',[]);
        set(gca,'yticklabel',[]);
            
        % Write classifier name (SVM)
        %y = 0.75 * yLimits(2);
        %x = -b/w(1) - w(2)/w(1) * y;
        %text(x, y, '\leftarrow SVM', 'FontSize', 12);
            
    % +=================================================================+ %
    % |                                                                 | %
    % |                           RBF Kernel                            | %
    % |                                                                 | %
    % +=================================================================+ %
    elseif strcmp( kernel_type, 'rbf' )
        
        % Get x-, y-axis limits
        xLimits = get(gca, 'XLim');
        yLimits = get(gca, 'YLim');

        % Create a grid of x and y points
        points_x1 = linspace(xLimits(1), xLimits(2), 50);
        points_x2 = linspace(yLimits(1), yLimits(2), 50);
        [X1, X2] = meshgrid(points_x1, points_x2);

        % Compute kernel function
        f = ones(length(points_x1),length(points_x2)) * (-rho);
        for i=1:l
            alpha_i = SVs(i,1);
            sv_i    = [SVs(i,2);SVs(i,3)];
            for j=1:length(points_x1)
                for k=1:length(points_x2)
                    x = [points_x1(k);points_x2(j)];
                    f(j,k) = f(j,k) + alpha_i * rbf_kernel(gamma, x, sv_i);
                end
            end    
        end

        % Plot contour at z = 0
        contour(X1, X2, f, [0 0], 'LineStyle' , '-', ...
                                  'Color'     , [0 0.8 1], ...
                                  'LineWidth' , 4 );
        contour(X1, X2, f, [0 0], 'LineStyle' , '--', ...
                                  'Color'     , [0 0.44706 0.74118], ...
                                  'LineWidth' , 3 );
        h = gcf;                             
    end
    
end