#ifndef FULLCOVARIANCEMATRIX_H
#define FULLCOVARIANCEMATRIX_H
#include <string>
#include <eigen3/Eigen/Dense>

using namespace std;

class FullCovarianceMatrix
{
    public:
        /* Ctor */
        FullCovarianceMatrix(int, string);

        /* Setters */
        void setDim(int);
        void setDocId(string);
        void setSigma(int, int, double);

        /* Getters */
        int getDim();
        string getDocId();
        Eigen::MatrixXd getSigma();

    private:
        int dim;
        string doc_id;
        Eigen::MatrixXd Sigma;
};

#endif // FULLCOVARIANCEMATRIX_H
