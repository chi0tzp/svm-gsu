#ifndef DIAGCOVARIANCEMATRIX_H
#define DIAGCOVARIANCEMATRIX_H
#include <string>
#include <eigen3/Eigen/Dense>

using namespace std;

class DiagCovarianceMatrix
{
    public:
        /* Ctor */
        DiagCovarianceMatrix(int, string);

        /* Setters */
        void setDim(int);
        void setDocId(string);
        void setSigma(int, double);

        /* Getters */
        int getDim();
        string getDocId();
        Eigen::VectorXd getSigma();

    private:
        int dim;
        string doc_id;
        Eigen::VectorXd Sigma;
};

#endif // DIAGCOVARIANCEMATRIX_H
