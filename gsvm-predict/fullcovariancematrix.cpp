#include "fullcovariancematrix.h"

/* Ctor */
FullCovarianceMatrix::FullCovarianceMatrix(int d, string id)
{
    dim      = d;
    doc_id   = id;
    Sigma.resize(d,d);
}

/* covariance_matrix: setters */
void FullCovarianceMatrix::setDim(int d){ dim = d; }
void FullCovarianceMatrix::setDocId(string id){ doc_id = id; }
void FullCovarianceMatrix::setSigma(int row_idx, int col_idx, double val){  Sigma(row_idx, col_idx) = val; }

/* covariance_matrix: getters */
int FullCovarianceMatrix::getDim(){ return dim; }
string FullCovarianceMatrix::getDocId(){ return doc_id; }
Eigen::MatrixXd FullCovarianceMatrix::getSigma(){ return Sigma; }
