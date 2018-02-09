function trainSVM(exp_dir, C, gamma)
    
     % Training files
     train_data_file = '../data/train_mean.dat';
     train_labels_file = '../data/train_labels.dat';
     model_file = sprintf( '%ssvm.model', exp_dir );
     
     % Train standard SVM
     train_cmd = sprintf( './svm-train-gep -h 0 -q -b 1 -t 2 -c %g -g %g %s %s %s', C,  gamma, train_data_file, train_labels_file, model_file );
     system( train_cmd );

end
