#include "svmgsumodel.h"
#include <string>
#include <vector>

#ifndef SVMGSUPARAMS_H
#define SVMGSUPARAMS_H

using namespace std;

class SvmGsuParams
{
    public:
        /* Ctor */
        SvmGsuParams(const SvmGsuModel&);

        /* SVM-GSU model */
        SvmGsuModel model;

        /* Getters */
        int getKernelType();
        double getGamma();
        double getSigmA();
        double getSigmB();
        int getVerbose();

        /* Parse command line arguments */
        void parseCommandLine(int, char**, char*, char*, char*, char*, char* );

        /* Read linear model file */
        void readModelFile(char *);

    private:
        int kernel_type;
        double gamma;
        double sigmA;
        double sigmB;
        int verbose;
        void tokenize(const string&, vector<string>&, const string&);
        void exitWithHelp();

};

#endif // SVMGSUPARAMS_H
