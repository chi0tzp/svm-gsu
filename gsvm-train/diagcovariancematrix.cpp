#include "diagcovariancematrix.h"


/* Ctor */
DiagCovarianceMatrix::DiagCovarianceMatrix(int d, string id)
{
    dim    = d;
    doc_id = id;
    Sigma.resize(d);
    Sigma.setZero(d);
}

/* covariance_matrix: setters */
void DiagCovarianceMatrix::setDim(int d){dim = d;}
void DiagCovarianceMatrix::setDocId(string id){doc_id = id;}
void DiagCovarianceMatrix::setSigma(int idx, double val){Sigma[idx] = val;}

/* covariance_matrix: getters */
int DiagCovarianceMatrix::getDim(){return dim;}
string DiagCovarianceMatrix::getDocId(){return doc_id;}
Eigen::VectorXd DiagCovarianceMatrix::getSigma(){return Sigma;}
