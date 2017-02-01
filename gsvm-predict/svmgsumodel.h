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
        void setDim(int);
        void setW(Eigen::VectorXd);
        void setAlpha(Eigen::VectorXd);
        void setSV(Eigen::MatrixXd);
        void setB(double);
        /* Getters */
        int getDim();
        Eigen::VectorXd getW();
        Eigen::VectorXd getAlpha();
        Eigen::MatrixXd getSV();
        double getB();

    private:
        int dim;
        /* Linear Model */
        Eigen::VectorXd w;
        double b;
        /* Kernel Model */
        Eigen::VectorXd alpha;
        Eigen::MatrixXd SV;

};

#endif // SVMGSUMODEL_H
