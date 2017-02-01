#include "svmgsumodel.h"

SvmGsuModel::SvmGsuModel()
{
    // Initialize linear model
    w.resize(1);
    w.setZero(1);
    b = 0.0;
}

/* Setters */
void SvmGsuModel::setDim(int d){dim = d;}
void SvmGsuModel::setW(Eigen::VectorXd v){w = v;}
void SvmGsuModel::setAlpha(Eigen::VectorXd v){alpha = v;}
void SvmGsuModel::setSV(Eigen::MatrixXd V){SV = V;}
void SvmGsuModel::setB(double bb){b = bb;}

/* Getters */
int SvmGsuModel::getDim(){return dim;}
Eigen::VectorXd SvmGsuModel::getW(){return w;}
Eigen::VectorXd SvmGsuModel::getAlpha(){return alpha;}
Eigen::MatrixXd SvmGsuModel::getSV(){return SV;}
double SvmGsuModel::getB(){return b;}
