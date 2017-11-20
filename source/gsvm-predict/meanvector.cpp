#include "meanvector.h"
#include <string>
#include <eigen3/Eigen/Dense>

/* Ctor */
MeanVector::MeanVector(int d, string id)
{
    dim = d;
    doc_id = id;
    x.resize(d);
    x.setZero(d);
}

/* Setters */
void MeanVector::setDim(int d){dim = d;}
void MeanVector::setDocId(string id){doc_id = id;}
void MeanVector::setXj(int idx, double val){x(idx) = val;}

/* Getters */
int MeanVector::getDim(){return dim;}
string MeanVector::getDocId(){return doc_id;}
Eigen::VectorXd MeanVector::getX(){return x;}
