#ifndef ISOCOVARIANCEMATRIX_H
#define ISOCOVARIANCEMATRIX_H
#include <string>

using namespace std;

class IsoCovarianceMatrix
{
    public:
        /* Ctor */
        IsoCovarianceMatrix(string);

        /* Setters */
        void setDocId(string);
        void setSigma(double);

        /* Getters */
        string getDocId();
        double getSigma();

    private:
        string doc_id;
        double Sigma;
};

#endif // ISOCOVARIANCEMATRIX_H
