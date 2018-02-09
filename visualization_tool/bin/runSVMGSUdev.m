% +---------------------------------------------------------------------+ %
% | Visuialization tool for SVM-GSU                                     | %
% |                                                                     | %
% | Author : Christos Tzelepis                                          | %
% | Email  : tzelepis@iti.gr                                            | %
% +---------------------------------------------------------------------+ %
clear; close; clc;

% ======================================================================= %
% Set experimental configuration (see setExpConfig.m)                     %
% ======================================================================= %
% Experimental configurations                                             %
% ----------------------------------------------------------------------- %
%   0: Manual (edit setConfig.m)                                          %
%   1: Linear kernel, full covariance matrices                            %
%   2: Linear kernel, diagonal covariance matrices                        %
%   3: Linear kernel, isotropic covariance matrices                       %
%   4: RBF kernel, isotropic covariance matrices                          %
% ======================================================================= %
exp_config = 0;
setExpConfig(exp_config);

% Load experimental parameters:
load('exp_config.mat');

% Generate training set
[MEAN_VECTORS, COV_MATS] = genTrainSet(kernel_type, cov_mat_type);

% Save training set
saveTrainSet(MEAN_VECTORS, COV_MATS, cov_mat_type);

% Plot training set
figure(1); hold on;
plotTrainSet(MEAN_VECTORS, COV_MATS, p); axis equal;

% Create experiment dir
exp_dir = sprintf('../experiments/%s/', kernel_type);
mkdirIfNotExist(exp_dir);

% ======================================================================= %
%                              [Standard SVM]                             %
% ======================================================================= %
trainSVM(exp_dir, svm_C, svm_gamma);
svm_model = parseSVM(kernel_type);
line_style = '-';
color = [0 0.44706 0.74118];
h1 = plotSVM(svm_model, kernel_type, svm_gamma, line_style, color);

% ======================================================================= %
%                               [SVM-GSU]                                 %
% ======================================================================= %
trainSVMGSU(kernel_type, cov_mat_type, svmgsu_lambda, svmgsu_gamma, sgd_T, sgd_k, exp_dir, p);
svmgsu_model = parseSVMGSU(kernel_type);
h2 = plotSVMGSU(svmgsu_model, kernel_type, svmgsu_gamma);
axis square;

return;

if (SHOW_LEGEND == 1)
    legend([h1, h2], {'LSVM', 'LSVM-GSU'}, ...
                     'Location', legend_position, ...
                     'FontSize', legend_font_size);
end

if (SAVE_FIG == 1)
    fig_filename = sprintf('../results/SVM_vs_SVMGSU_%s', cov_mat_type);
    print(fig_filename, '-depsc');
    %close;
end
