#ifndef SVMGSUPROBLEM_H
#define SVMGSUPROBLEM_H

#include "svmgsuparams.h"
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
        vector<int> getPrincAxes(int);


        /* Read input data functions */
        void readInputDataFull(char*, char*, char*);
        void readInputDataDiag(char*, char*, char*);
        void readInputDataIso(char*, char*, char*);

        /* Eigendecomposition */
        void eigenDecomp();

        /* Solvers */
        void solveLSVMGSUxSpace();
        void solveLSVMGSUzSpace();

        /* Write linear model file */
        void writeLinearModelFile(char *);

        /* Write kernel model file */
        void writeKernelModelFile(char *);

    private:
        /* SVMGSU's params */
        SvmGsuParams params;  // SVM-GSU's parameters (see class SvmGsuParams)

        /* Problem's params */
        int l;                // Cardinality of the training set
        int dim;              // Dimensionality of original feature space
        vector<double> y;     // Set of truth labels


        /*******************************************************************************
         *                   Learning in the original feature space                    *
         *******************************************************************************/
        // -- Data --
        Eigen::MatrixXd X;                 // Matrix of mean vectors
        vector<Eigen::VectorXd> x;         // Set of mean vectors
        vector<Eigen::MatrixXd> Sigma_xf;  // Set of full covariance matrices
        vector<Eigen::VectorXd> Sigma_xd;  // Set of (diagonal) covariance matrices
        vector<double>          Sigma_xi;  // Set of (isotropic) covariance matrices

        // -- Methods --
        Eigen::VectorXd computeSumOfLossGradFullXspace(Eigen::VectorXd, double, vector<int>);
        Eigen::VectorXd computeSumOfLossGradDiagXspace(Eigen::VectorXd, double, vector<int>);
        Eigen::VectorXd computeSumOfLossGradIsoXspace(Eigen::VectorXd, double, vector<int>);

        // -- Models --


        // Linear model: w^T*x + b = 0
        Eigen::VectorXd w;     // Normal vector to H: w^T*x + b = 0
        double b;              // Bias term of H: w^T*x + b = 0
        Eigen::MatrixXd K;     // Kernel matrix
        Eigen::VectorXd alpha; //


        /*******************************************************************************
         *                        Learning in linear subspaces                         *
         *******************************************************************************/
        // -- Data --
        vector<Eigen::VectorXd> z;         // Set of mean vectors (in subspaces)
        vector<Eigen::VectorXd> Sigma_zd;  // Set of covariance matrices (in subspaces)
        vector<vector<int>> princ_axes;    // Set of principle axes

        // -- Methods --
        Eigen::VectorXd computeSumOfLossGradDiagZspace(Eigen::VectorXd, double, vector<int>);

        // -- Models --

        /* Decision values */
        Eigen::VectorXd deci;

        /* Platt scaling params */
        double sigmA;
        double sigmB;

        /* Auxiliary functions */
        double getLabel(double);
        void tokenize(const string&, vector<string>&, const string&);
        vector<string> findDocIdsIntersection( vector<string>, vector<string>);
        void initW(int, double);
        void computeDecValues();
        void plattScaling();

};

#endif // SVMGSUPROBLEM_H
