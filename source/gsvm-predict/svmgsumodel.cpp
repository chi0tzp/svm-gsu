#include "svmgsumodel.h"

/* Default ctor */
SvmGsuModel::SvmGsuModel(){}

/* Setters */
void SvmGsuModel::setModelDim(int d){model_dim = d;}
void SvmGsuModel::setW(Eigen::VectorXd v){w = v;}
void SvmGsuModel::setAlpha(Eigen::VectorXd v){alpha = v;}
void SvmGsuModel::setSV(Eigen::MatrixXd V){SV = V;}
void SvmGsuModel::setB(double bb){b = bb;}

/* Getters */
int SvmGsuModel::getModelDim(){return model_dim;}
Eigen::VectorXd SvmGsuModel::getW(){return w;}
Eigen::VectorXd SvmGsuModel::getAlpha(){return alpha;}
Eigen::MatrixXd SvmGsuModel::getSV(){return SV;}
double SvmGsuModel::getB(){return b;}
