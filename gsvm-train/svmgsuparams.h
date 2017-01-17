#ifndef SVMGSUPARAMS_H
#define SVMGSUPARAMS_H
#include <string>

using namespace std;

class SvmGsuParams
{
    public:
        /* Default ctor */
        SvmGsuParams();

        /* Getters */
        int getVerbose();
        int getT();
        int getK();
        int getCovMatType();
        std::string getCovMatTypeVerb();
        int getKernelType();
        double getLambda();
        double getGamma();
        double getP();

        /* Member functions */
        void parseCommandLine( int, char**, char*, char*, char*, char*, char* );


    private:
        int verbose;           // Verbose mode flag
        int T;                 // Number of iterations for SGD
        int k;                 // Sampling size for SGD
        int cov_mat_type;      // Covariance matrices type (0: full, 1: diagonal, 2: isotropic)
        int kernel_type;       // Kernel type (0: linear, 2: rbf)
        double lambda;         // SVM-GSU's lambda parameter
        double gamma;          // KSVM-GSU's gamma parameters
        double p;              // Variance fraction to be preserved in subspace learning (in (0,1])
        void exitWithHelp();   //

};

#endif // SVMGSUPARAMS_H
