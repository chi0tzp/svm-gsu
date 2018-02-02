clear; close; clc;

N_pos = 1000;
N_neg = 9000;
dim = 345;

EXTRACT_FULL_COV = 0;
EXTRACT_DIAG_COV = 1;
EXTRACT_ISO_COV = 1;

fprintf('Generate artifical data:\n');
fprintf(' -- Dimenionality: %d\n', dim);
fprintf(' -- Total generated data: %d\n', N_pos+N_neg);
fprintf('    -- Positive data: %d\n', N_pos);
fprintf('    -- Negative data: %d\n', N_neg);

tic;
% ======================================================================= %
% Generate negative population
% ======================================================================= %
X_mean_neg = zeros(N_neg, dim);
X_cov_full_neg = cell(N_neg, 1);
X_cov_diag_neg = cell(N_neg, 1);
X_cov_iso_neg = cell(N_neg, 1);
X_label_neg = -1 * ones(N_neg, 1);
mu_neg = -2 + rand(dim, 1);
for i=1:N_neg
    X_mean_neg(i,:) = mvnrnd(mu_neg, 0.01*eye(dim), 1);
    if (EXTRACT_FULL_COV)
        X_cov_full_neg{i} = wishrnd(0.005*eye(dim), dim);
    end
    if (EXTRACT_DIAG_COV)
        X_cov_diag_neg{i} = 0.05 * rand(dim,1);
    end
    if (EXTRACT_ISO_COV)
        X_cov_iso_neg{i} = 0.05 * rand;
    end
end

% ======================================================================= %
% Generate positive population
% ======================================================================= %
X_mean_pos = zeros(N_pos, dim);
X_cov_full_pos = cell(N_pos, 1);
X_cov_diag_pos = cell(N_pos, 1);
X_cov_iso_pos = cell(N_pos, 1);
X_label_pos = +1 * ones(N_pos, 1);
mu_pos = -1 + 2*rand(dim, 1);
for i=1:N_pos
    X_mean_pos(i,:) = mvnrnd(mu_pos, 0.01*eye(dim), 1);
    if (EXTRACT_FULL_COV)
        X_cov_full_pos{i} = wishrnd(0.005*eye(dim), dim);
    end
    if (EXTRACT_DIAG_COV)
        X_cov_diag_pos{i} = 0.05 * rand(dim,1);
    end
    if (EXTRACT_ISO_COV)
        X_cov_iso_pos{i} = 0.05 * rand;
    end
end
toc;

% +---------------------------------------------------------------------+ %
% |                                                                     | %
% |                       Write generated dataset                       | %
% |                                                                     | %
% +---------------------------------------------------------------------+ %
fprintf('Write generated data:\n');
data_dir = sprintf('./data/%dD/%d/', dim, N_neg+N_pos);
mkdirIfNotExist(data_dir);


% ======================================================================= %
% Write mean vectors
% ======================================================================= %
fprintf(' -- Mean vectors...\n');
mean_vectors = [X_mean_neg ; X_mean_pos];
mean_file = sprintf('%smean.dat', data_dir);
fid = fopen(mean_file, 'a');
for i=1:size(mean_vectors, 1)
    vector = mean_vectors(i,:);
    z = 1:length(vector);
    fprintf(fid, 'DOC%d', i);
    fprintf(fid, ' %d:%g', [z; vector(:)']);
    fprintf(fid, '\n');
end
fclose(fid);

% ======================================================================= %
% Write labels
% ======================================================================= %
fprintf(' -- Labels...\n');
labels = [X_label_neg ; X_label_pos];
labels_file = sprintf('%slabels.dat', data_dir);
fid = fopen(labels_file, 'a');
for i=1:size(labels, 1)
    fprintf(fid, 'DOC%d %d\n', i, labels(i));
end
fclose(fid);

% ======================================================================= %
% Write full covariance matrices
% ======================================================================= %
if (EXTRACT_FULL_COV)
    fprintf(' -- Full covariance matrices...\n');
    cov_full = [X_cov_full_neg ; X_cov_full_pos];
    cov_full_file = sprintf('%scov_full.dat', data_dir);
    fid = fopen(cov_full_file, 'a');
    for i=1:size(cov_full, 1)
        fprintf(fid, 'DOC%d', i);
        cov_mat = cov_full{i};
        for j=1:dim
            vector = cov_mat(j,:);
            z = 1:length(vector);
            row = j * ones(size(vector));
            fprintf(fid, ' %d,%d:%g', [row; z; vector(:)']);
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

% ======================================================================= %
% Write diagonal covariance matrices
% ======================================================================= %
if (EXTRACT_DIAG_COV)
    fprintf(' -- Diagonal covariance matrices...\n');
    cov_diag = [X_cov_diag_neg ; X_cov_diag_pos];
    cov_diag_file = sprintf('%scov_diag.dat', data_dir);
    fid = fopen(cov_diag_file, 'a');
    for i=1:size(cov_diag, 1)
        vector = cov_diag{i};
        z = 1:length(vector);
        fprintf(fid, 'DOC%d', i);
        fprintf(fid, ' %d,%d:%g', [z; z; vector(:)']);
        fprintf(fid, '\n');
    end
    fclose(fid);
end

% ======================================================================= %
% Write isotropic covariance matrices
% ======================================================================= %
if (EXTRACT_ISO_COV)
    fprintf(' -- Isotropic covariance matrices...\n');
    cov_iso = [X_cov_iso_neg ; X_cov_iso_pos];
    cov_iso_file = sprintf('%scov_iso.dat', data_dir);
    fid = fopen(cov_iso_file, 'a');
    for i=1:size(cov_iso, 1)
        fprintf(fid, 'DOC%d 1,1:%g\n', i, cov_iso{i});
    end
    fclose(fid);
end



