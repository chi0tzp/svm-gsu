function h = plotSVMGSU(svmgsu_model, kernel_type, gamma)

    line_width = 2.5;
    % +=================================================================+ %
    % |                                                                 | %
    % |                         Linear Kernel                           | %
    % |                                                                 | %
    % +=================================================================+ %
    if strcmp(kernel_type, 'linear')
        
        % Get model's w, b
        w = svmgsu_model(1:2);
        b = svmgsu_model(3);

        % Get x-, y-axis limits
        xLimits = get(gca,'XLim');
        yLimits = get(gca,'YLim');    

        % Plot separating hyperplane (line)
        % -- Case: Horizontal line
        if (w(1) == 0)
            y = -b/w(2);
            h = line( y*[1 1], xLimits);
        % -- Case: Vertical line
        elseif (w(2) == 0)
            x = -b/w(1);
            h = line( x*[1 1], yLimits);
        % -- Otherwise
        else
            line_cmd = sprintf('y=-%g*x-%g', w(1)/w(2), b/w(2));
            h = ezplot(line_cmd, [xLimits yLimits]);
        end

        set(h , 'LineStyle' , '-', ...
                'Color'     , [1     0          0], ...
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
        %text(x, y, '\leftarrow SVM-GSU', 'FontSize', 12);
        
    % +-----------------------------------------------------------------+ %
    % |                                                                 | %
    % |                           RBF Kernel                            | %
    % |                                                                 | %
    % +-----------------------------------------------------------------+ %
    elseif strcmp( kernel_type, 'rbf' )
        
        SVs = svmgsu_model{1};
        rho = svmgsu_model{2};
        l = size( SVs, 1 );

        % Get x-, y--axis limits
        xLimits = get(gca, 'XLim');
        yLimits = get(gca, 'YLim');

        % Create a grid of x and y points --
        points_x1 = linspace(xLimits(1), xLimits(2), 50);
        points_x2 = linspace(yLimits(1), yLimits(2), 50);
        [X1, X2] = meshgrid(points_x1, points_x2);

        % Compute kernel function
        %rho = 0;
        f = ones(length(points_x1),length(points_x2)) * (rho * 0);
        for i=1:l
            alpha_i = SVs(i,1);
            sv_i = [SVs(i,2);SVs(i,3)];
            for j=1:length(points_x1)
                for k=1:length(points_x2)
                    x = [points_x1(k);points_x2(j)];
                    f(j,k) = f(j,k) + alpha_i * rbf_kernel(gamma, x, sv_i);
                end
            end    
        end

        surf(X1, X2, f);
        shading interp;
        light;lighting phong;
        alpha(.75)
        contour3(X1, X2, f, 20, 'k')

        % Plot contour at z = 0
        contour(X1, X2, f, [0 0], 'LineStyle' , '-', ...
                                  'Color'     , 'r', ...
                                  'LineWidth' , 4 );
        contour(X1, X2, f, [0 0], 'LineStyle' , '--', ...
                                  'Color'     , 'm', ...
                                  'LineWidth' , 3 );
        h = gcf;
    end
    
end