#ifndef SVMGSUPROBLEM_H
#define SVMGSUPROBLEM_H

#include "svmgsuparams.h"
#include "outputscore.h"
#include <algorithm>
#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>

using namespace std;


class SvmGsuProblem
{
    public:

        /* Ctor */
        SvmGsuProblem(const SvmGsuParams&);

        /* Getters */
        int getL();
        int getDim();
        vector<double> getLabels();

        /* Read input data (mean vectors and truth labels) */
        void readInputData(char*, char*);

        /* Predict */
        vector<OutputScore> predict();

        /* Write output file */
        void writeOutputFile(vector<OutputScore>, char*);

    private:

        /* SVMGSU's params */
        SvmGsuParams params;
        int l;
        int dim;

        /* Input data */
        vector<string> doc_id;      // Set of document ids
        vector<double> y;           // Set of truth labels
        vector<Eigen::VectorXd> x;  // Set of mean vectors

        /* Auxiliary functions */
        void tokenize(const string&, vector<string>&, const string&);
        vector<string> findDocIdsIntersection( vector<string>, vector<string>);
        int sgn(double);

};

#endif // SVMGSUPROBLEM_H
