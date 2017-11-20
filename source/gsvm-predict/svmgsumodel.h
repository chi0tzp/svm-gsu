#ifndef SVMGSUMODEL_H
#define SVMGSUMODEL_H

#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>


class SvmGsuModel
{
    public:
        /* Default ctor */
        SvmGsuModel();

        /* Setters */
        void setModelDim(int);
        void setW(Eigen::VectorXd);
        void setAlpha(Eigen::VectorXd);
        void setSV(Eigen::MatrixXd);
        void setB(double);

        /* Getters */
        int getModelDim();
        Eigen::VectorXd getW();
        double getB();
        Eigen::VectorXd getAlpha();
        Eigen::MatrixXd getSV();

    private:
        /* Linear Model */
        int model_dim;
        Eigen::VectorXd w;
        double b;

        /* Kernel (RBF) Model */
        Eigen::VectorXd alpha;
        Eigen::MatrixXd SV;

};

#endif // SVMGSUMODEL_H
