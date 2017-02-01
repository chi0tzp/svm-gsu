#include "isocovariancematrix.h"

/* Ctor */
IsoCovarianceMatrix::IsoCovarianceMatrix(string id)
{
    doc_id = id;
    Sigma  = 0.0;
}

/* Setters */
void IsoCovarianceMatrix::setDocId(string id){ doc_id = id; }
void IsoCovarianceMatrix::setSigma(double val){  Sigma = val; }

/* Getters */
string IsoCovarianceMatrix::getDocId(){ return doc_id; }
double IsoCovarianceMatrix::getSigma(){ return Sigma; }
