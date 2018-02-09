function svmgsu_model = parseSVMGSU(kernel_type)

    fid = fopen(sprintf('../experiment/%s/svmgsu.model', kernel_type), 'r');
    % +=================================================================+ %
    % |                         Linear Kernel                           | %
    % +=================================================================+ %
    if strcmp( kernel_type, 'linear' )
        
        tline = fgetl(fid);
        w = [];
        while ischar(tline)        
            % Split line by ' '
            temp = regexp(tline, ' ', 'split');
            % Get vector w (H: w^T*x+b=0)
            if strcmp(temp{1},'w')
               for i=2:3
                  w = [w str2double(temp{i})];
               end
            end
            % Get bias term b (H: w^T*x+b=0)
            if strcmp(temp{1},'b')
               b = str2double(temp{2});
            end
            tline = fgetl(fid);
        end
        svmgsu_model = [w b];
        
    % +=================================================================+ %
    % |                           RBF Kernel                            | %
    % +=================================================================+ %
    elseif strcmp(kernel_type, 'rbf')
        
        tline = fgetl(fid);
        getSV = 0;
        line_cnt = 1;
        SVs = [];
        while ischar(tline)
            % Get support vectors
            if getSV
               temp = regexp(tline, ' ', 'split');
               temp = temp(1:end);
               alpha_i = str2double(temp(1,1));
               x_temp = regexp(temp(1,2), ':', 'split');
               x = str2double(x_temp{1}{2});
               y_temp  = regexp(temp(1,3), ':', 'split');
               y = str2double(y_temp{1}{2});
               SVs(line_cnt,:) = [alpha_i x y];
               line_cnt = line_cnt + 1;
            end
            % Get SV flag
            if strcmp(tline,'SV')
                getSV = 1;
            end
            % Get rho
            temp = regexp(tline, ' ', 'split');
            if strcmp(temp{1}, 'b')
                rho = str2double(temp{2});
            end
            tline = fgetl(fid);
        end
        svmgsu_model{1} = SVs;
        svmgsu_model{2} = rho;
    end
    fclose(fid);

end