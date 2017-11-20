#ifndef MEANVECTOR_H
#define MEANVECTOR_H

#include <string>
#include <eigen3/Eigen/Dense>

using namespace std;

class MeanVector
{
    public:
        /* ctor */
        MeanVector(int, string);

        /* Setters */
        void setDim(int);
        void setDocId(string);
        void setXj(int,double);

        /* Getters */
        int getDim();
        string getDocId();
        Eigen::VectorXd getX();

    private:
        int dim;
        string doc_id;
        Eigen::VectorXd x;
};

#endif // MEANVECTOR_H
