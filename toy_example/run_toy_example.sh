#!/bin/bash
################################################################################
# svm-gsu: toy example                                                         #
#                                                                              #
# [1] C. Tzelepis, V. Mezaris, I. Patras, "Linear Maximum Margin Classifier    #
# for Learning from Uncertain Data", IEEE Transactions on Pattern Analysis and #
# Machine Intelligence, accepted for publication.                              #
# DOI:10.1109/TPAMI.2017.2772235.                                              #
################################################################################
# Define experiment's parameters
verb=1            # Set verbose mode (On: 1 / Off: 0)
cov_mat="diag"    # Choose covariance matrices type:
                  # -- "full": Full covariance matrices
                  # -- "diag": Diagonal covariance matrices
                  # -- "iso": Isotropic covariance matrices
kernel=0          # Choose kernel type:
                  # -- 0: linear kernel
                  # -- 2: RBF kernel
lambda=0.001      # Set SVM-GSU's reguralization (lambda) parameter
T=3000            # Set number of iterations for the SGD algorithm
k=100             # Set the sampling size for the SGD algorithm (1<=k<=T)
p=1.00            # Set the fraction of the total variance preserved for each
                  # input Gaussian -- 0.0<p<=1.0 (learning in linear subspaces
                  # -- See Sect. 3.2 of [1])

# Define input/output directories
data_dir="./data"
models_dir="./models"
output_dir="./output"
bin_dir_train=../build/gsvm-train
bin_dir_predict=../build/gsvm-predict

# Make output directories, if ther do not already exist
mkdir -p ${models_dir}
mkdir -p ${output_dir}

if [[ ${kernel} -eq 0 ]];
then
    kernel_type="linear"
elif [[ ${kernel} -eq 2 ]];
then
    kernel_type="rbf"
else
    echo "Invalid kernel type. Abort."
    exit;
fi

if [ "${cov_mat}" == "full" ]
then
  cov_mat_opt=0
elif [ "${cov_mat}" == "diag" ]
then
  cov_mat_opt=1
elif [ "${cov_mat}" == "iso" ]
then
  cov_mat_opt=2
fi

# Training files
train_mean_vec="${data_dir}/mean.dat"
train_labels="${data_dir}/labels.dat"
train_cov_mat="${data_dir}/cov_${cov_mat}.dat"
train_model_file="${models_dir}/svmgsu_${kernel_type}_${cov_mat}_T${T}_k${k}_l${lambda}.model"

# Testing files
test_mean_vec="${data_dir}/mean.dat"
test_labels="${data_dir}/labels.dat"
test_cov_mat="${data_dir}/cov_${cov_mat}.dat"
test_output_file="${output_dir}/svmgsu_${kernel_type}_${cov_mat}_T${T}_k${k}_l${lambda}.output"
test_metrics_file="${output_dir}/svmgsu_${kernel_type}_${cov_mat}_T${T}_k${k}_l${lambda}.metrics"

# Construct training options

opt="-v ${verb} -t ${kernel} -d ${cov_mat_opt} -l ${lambda} -T ${T} -k ${k} -p ${p}"

# Construct training command
train_cmd="${bin_dir_train}/gsvm-train ${opt} ${train_mean_vec} ${train_labels} ${train_cov_mat} ${train_model_file}"

# Train SVM-GSU
echo -n "Re-build gsvm-train..."
#make clean -C ../build/gsvm-train
make -s -C ../build/gsvm-train && echo "Done!"
echo "Run gsvm-train..."
eval ${train_cmd}

# Evaluate the above trained SVM-GSU model
echo -n "Re-build gsvm-predict..."
#make clean -C ../build/gsvm-predict
make -s -C ../build/gsvm-predict && echo "Done!"
predict_opt="-v ${verb} -t ${test_labels} -m ${test_metrics_file}"
predict_cmd="${bin_dir_predict}/gsvm-predict ${predict_opt} ${test_mean_vec} ${train_model_file} ${test_output_file}"
echo "Run gsvm-predict..."
eval ${predict_cmd}
