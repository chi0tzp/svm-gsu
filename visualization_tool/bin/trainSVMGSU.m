function trainSVMGSU(kernel_type, cov_mat_type, lambda, gamma, T, k, exp_dir, p)
    
     BUILD_SVMGSU = 1;
     
     % Build gsvm-train and copy gsvm-train to ./
     if (BUILD_SVMGSU)
         system('sh build_and_copy.sh');
     end
     
     % Training files
     train_mean_file = '../data/train_mean.dat';
     train_gt_file = '../data/train_labels.dat';
     train_cov_file = '../data/train_cov.dat';
     model_file = sprintf('%ssvmgsu.model', exp_dir);
     
     % Select covariance matrices type
     if strcmp(cov_mat_type, 'full')
         diag = 0;
     elseif strcmp(cov_mat_type, 'diag')
         diag = 1;
    elseif strcmp(cov_mat_type, 'iso')
         diag = 2;
     end
        
     % Train SVM-GSU
     if strcmp( kernel_type, 'linear' )
         train_cmd = sprintf( './gsvm-train -v 1 -t 0 -d %d -p %g -l %g -T %d -k %d %s %s %s %s', diag, p, lambda, T, k, train_mean_file, train_gt_file, train_cov_file, model_file );
     elseif strcmp( kernel_type, 'rbf' )
         train_cmd = sprintf( './gsvm-train -v 1 -t 2 -d %d -l %g -g %g  -T %d -k %d %s %s %s %s', diag, lambda, gamma, T, k, train_mean_file, train_gt_file, train_cov_file, model_file );
     end
     system( train_cmd );

end
