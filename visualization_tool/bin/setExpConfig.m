% ======================================================================= % 
% Experimental configurations                                             %
% ----------------------------------------------------------------------- %
%   0: Manual settings                                                    %
%   1: Linear kernel, full covariance matrices                            %
%   2: Linear kernel, diagonal covariance matrices                        %
%   3: Linear kernel, isotropic covariance matrices                       %
%   4: RBF kernel, isotropic covariance matrices                          %
%                                                                         %
% Set the following parameters for the various scenarios:                 %
% ----------------------------------------------------------------------- %
%  -- kernel_type   : Kernel type ('linear', 'rbf')                       %
%  -- cov_mat_type  : Covariance matrices type ('full', 'diag', 'iso')    %
%  -- sgd_T         : Stochastic Gradient Descent number of iterations    %
%  -- sgd_k         : Stochastic Gradient Descent sampling size           %
%  -- svmgsu_lambda : lambda parameter for SVM-GSU (linear)               %
%  -- svmgsu_gamma  : gamma parameter for SVM-GSU (linear, RBF)           %
%  -- svm_C         : C parameter for standard SVM (linear)               %
%  -- svm_gamma     : gamma parameter for standard SVM (linear, RBF)      %
% ======================================================================= %
function setExpConfig(exp_config)
    
    KERNEL_TYPES = {'linear', 'rbf'};
    COV_MAT_TYPES = {'full', 'diag', 'iso'};
    
    % ======================== 0: Manual settings ======================= %
    if (exp_config == 0)
        kernel_type = KERNEL_TYPES{2};
        cov_mat_type = COV_MAT_TYPES{3};
        sgd_T = 6000;
        sgd_k = 4;
        p = 1.00;
        svmgsu_lambda = 2^(9);
        svmgsu_gamma = 2^(8);
        svm_C = 100;
        svm_gamma = 0.1;

    % ============ 1: Linear kernel, full covariance matrices =========== %
    elseif (exp_config == 1)
        kernel_type = KERNEL_TYPES{1};
        cov_mat_type = COV_MAT_TYPES{1};
        sgd_T = 6000;
        sgd_k = 6;
        p = 1.00;
        svmgsu_lambda = 2^(-17);
        svmgsu_gamma = 0.1;
        svm_C = 100;
        svm_gamma = 0.1;
    
    % ========== 2: Linear kernel, diagonal covariance matrices ========= %
    elseif (exp_config == 2)
        kernel_type = KERNEL_TYPES{1};
        cov_mat_type = COV_MAT_TYPES{2};
        sgd_T = 6000;
        sgd_k = 6;
        p = 1.00;
        svmgsu_lambda = 2^(-17);
        svmgsu_gamma = 0.1;
        svm_C = 100;
        svm_gamma = 0.1;
    
    % ========= 3: Linear kernel, isotropic covariance matrices ========= %
    elseif (exp_config == 3)
        kernel_type = KERNEL_TYPES{1};
        cov_mat_type = COV_MAT_TYPES{3};
        sgd_T = 6000;
        sgd_k = 4;
        p = 1.00;
        svmgsu_lambda = 2^(-17);
        svmgsu_gamma = 0.1;
        svm_C = 10;
        svm_gamma = 0.5;

    % ========== 4: RBF kernel, isotropic covariance matrices =========== %
    elseif (exp_config == 4)
        kernel_type = KERNEL_TYPES{2};
        cov_mat_type = COV_MAT_TYPES{3};
        sgd_T = 6000;
        sgd_k = 4;
        p = 1.00;
        svmgsu_lambda = 2^(9);
        svmgsu_gamma = 2^(8);
        svm_C = 100;
        svm_gamma = 2^(4);
    end
    % =================================================================== %
    
    SAVE_FIG = 0;
    SHOW_LEGEND = 0;
    legend_position = 'SouthEast';
    legend_font_size = 18;
    save('exp_config.mat');
    clear;
    
end