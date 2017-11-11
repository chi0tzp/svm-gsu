#!/bin/bash

# Define experiment's parameters
verb=1
cov_mat="diag"
kernel=0
lambda=100
T=3000
k=5000
#k=1000
p=1.00

# Define input/output directories
data_dir="./data"
models_dir="./models"
output_dir="./output"
bin_dir=../build/gsvm-train

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

# Construct training options
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
opt="-v ${verb} -t ${kernel} -d ${cov_mat_opt} -l ${lambda} -T ${T} -k ${k} -p ${p}"

# Construct training command
train_cmd="${bin_dir}/gsvm-train ${opt} ${train_mean_vec} ${train_labels} ${train_cov_mat} ${train_model_file}"

# Train SVM-GSU
echo "Training a ${kernel_type} SVM-GSU..."
eval ${train_cmd}
