# svm-gsu

A C++ framework for training/testing Support Vector Machine with Gaussian Sample Uncertainty (SVM-GSU).

This is the implementation code for the Linear SVM Gaussian Sample Uncertainty (LSVM-GSU), which was first proposed in [1], and the Kernel SVM Gaussian Sample Uncertainty (KSVM-iGSU), proposed in [2]. If you want to use one of the above classifiers, please consider citing the appropriate references (see at the end of this document).

## 0. Prerequisites and build guidelines

The code is built in C++11 using the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. In order to build the code, you need to 

```
- Eigen ??.??
- ???
```

### Linux

### Windows



## 1. Linear SVM with Gaussian Sample Uncertainty (LSVM-GSU) [1]

### Motivation

In our method we consider that our training examples are multivariate Gaussian distributions with known means and covariance matrices, each example having a different covariance matrix expressing the uncertainty around its mean. This is illustrated in the figure below

<p align="center">
  <img src="images/svmgsu_motivation.jpg" width="300" alt="SVM-GSU's motivation"/>
</p>

where the shaded regions are bounded by iso-density loci of the Gaussians, and the means of the Gaussians for examples of the positive and negative classes are located at "x" and "o" respectively. A classical linear SVM formulation (**LSVM**) would consider only the means of the Gaussians as training examples and, by optimizing the soft margin using the hinge loss and a regularization term, would arrive at the separating hyperplane depicted by the dashed line. In our formulation (**LSVM-GSU**), we optimize for the soft margin using the same regularization but the *expected* value of the hinge loss, where the expectation is taken under the given Gaussians. By doing so, we take into consideration the various uncertainties and arrive at a drastically different decision border, depicted by the solid line. 



## 2. Kernel SVM with Isotropic Gaussian Sample Uncertainty (KSVM-iGSU) [2,3]

Not available yet...



## 3. Visualization of LSVM-GSU/KSVM-iGSU

A visualization tool build in Matlab is available under XXX/





## References

[1] Tzelepis, Christos, Vasileios Mezaris, and Ioannis Patras. "Linear Maximum Margin Classifier for Learning from Uncertain Data." *IEEE Transactions on pattern analysis and machine intelligence* XX.YY (2017): pppp-pppp.

[2] Tzelepis, Christos, Vasileios Mezaris, and Ioannis Patras. "Video event detection using kernel support vector machine with isotropic gaussian sample uncertainty (KSVM-iGSU)." *International Conference on Multimedia Modeling*. Springer, Cham, 2016.

[3] Tzelepis, Christos, Eftichia Mavridaki, Vasileios Mezaris, and Ioannis Patras. "Video aesthetic quality assessment using kernel Support Vector Machine with isotropic Gaussian sample uncertainty (KSVM-iGSU)." In *Image Processing (ICIP), 2016 IEEE International Conference on*, pp. 2410-2414. IEEE, 2016.