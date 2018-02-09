function svm_model = parseSVM(kernel_type)

    fid = fopen(sprintf('../experiment/%s/svm.model', kernel_type), 'r');
    tline = fgetl(fid);
    getSV = 0;
    line_cnt = 1;
    SVs = [];
    while ischar(tline)
        % Get support vectors
        if getSV
           temp = regexp(tline, ' ', 'split');
           temp = temp(1:end-1);
           alpha_i = str2double(temp(1,1));
           x_temp = regexp(temp(1, 2), ':', 'split');
           x = str2double(x_temp{1}{2});
           y_temp = regexp(temp(1, 3), ':', 'split');
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
        if strcmp( temp{1}, 'rho' )
            rho = str2double( temp{2} );
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    svm_model{1} = SVs;
    svm_model{2} = rho;
    
end