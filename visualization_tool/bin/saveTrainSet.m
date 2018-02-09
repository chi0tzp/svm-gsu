function saveTrainSet(MEAN_VECTORS, COV_MATS, cov_mat_type)
    
    % Training Set's filenames
    %
    train_mean_file = '../data/train_mean.dat';
    train_cov_file = '../data/train_cov.dat';
    train_gt_file = '../data/train_labels.dat';
    
    mfd = fopen(train_mean_file, 'w');
    gfd = fopen(train_gt_file, 'w');
    cfd = fopen(train_cov_file, 'w');
    
    l = size(MEAN_VECTORS,1);
    dim = 2;
    for i=1:l
        
        % === Write mean vectors file ===
        new_mean = sprintf('doc%d', i);
        for j=1:dim
            new_mean = sprintf('%s %d:%g', new_mean, j, MEAN_VECTORS(i,j));
        end
        fprintf(mfd, '%s\n', new_mean);
        
        % === Write ground truth file ===
        fprintf(gfd, 'doc%d %d\n', i, MEAN_VECTORS(i, dim+1));
        
        % === Write covariance matrices file ===
        cov_mat = COV_MATS{i};
        new_cov = sprintf( 'doc%d', i );
        if strcmp( cov_mat_type, 'iso' )
            new_cov = sprintf('%s %d,%d:%g', new_cov, 1, 1, cov_mat(1,1));
        elseif strcmp(cov_mat_type, 'diag')
            for j=1:dim
               new_cov = sprintf('%s %d,%d:%g', new_cov, j, j, cov_mat(j,j));
            end
        elseif strcmp(cov_mat_type, 'full')
            for j=1:dim
                for k=1:dim
                    new_cov = sprintf( '%s %d,%d:%g', new_cov, j, k, cov_mat(j,k) );
                end
            end
        end
        fprintf(cfd, '%s\n', new_cov);
    end
    fclose(mfd);
    fclose(gfd);
    fclose(cfd);
    
end