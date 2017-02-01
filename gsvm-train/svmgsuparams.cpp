#include "svmgsuparams.h"
#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <string>

using namespace std;

/* Default ctor */
SvmGsuParams::SvmGsuParams()
    : verbose(0),
      T(1000),
      k(1),
      cov_mat_type(0),
      kernel_type(0),
      lambda(1.0),
      gamma(0.1),
      p(1.0)
{}

/* Getters */
int SvmGsuParams::getVerbose(){return verbose;}
int SvmGsuParams::getT(){return T;}
int SvmGsuParams::getK(){return k;}
int SvmGsuParams::getKernelType(){return kernel_type;}
double SvmGsuParams::getLambda(){return lambda;}
double SvmGsuParams::getGamma(){return gamma;}
double SvmGsuParams::getP(){return p;}
int SvmGsuParams::getCovMatType(){return cov_mat_type;}
string SvmGsuParams::getCovMatTypeVerb()
{
    string cov_mat_type_str;
    switch (cov_mat_type)
    {
        case 0:
            cov_mat_type_str="full";
            break;
        case 1:
            cov_mat_type_str="diagonal";
            break;
        case 2:
            cov_mat_type_str="isotropic";
            break;
        default:
            break;
    }
    return cov_mat_type_str;
}


/* Member functions */
void SvmGsuParams::parseCommandLine(int argc,
                                    char** argv,
                                    char* mean_vectors_filename,
                                    char* labels_filename,
                                    char* covariance_matrices_filename,
                                    char* model_filename)
{
    int i;
    // -- Parse options --
    for (i=1; i<argc; i++)
    {
        if (argv[i][0]!='+')
            break;
        if (++i>=argc)
            exitWithHelp();

        switch (argv[i-1][1])
        {
            /* Select kernel type: -t {0, 2} */
            case 't':
                kernel_type = atoi(argv[i]);
                break;
            /* Diagonal covariance matrices: -d {0, 1} */
            case 'd':
                cov_mat_type = atoi(argv[i]);
                break;
            /* Select C parameter: -c <C_PARAM> */
            case 'l':
                lambda = atof(argv[i]);
                break;
            /* Select gamma parameter: -g <GAMMA_PARAM> */
            case 'g':
                gamma = atof(argv[i]);
                break;
            /* Select fraction of variance (for subspace learning) */
            case 'p':
                p =  atof(argv[i]);
                break;
            /* Select number of SGD iterations T */
            case 'T':
                T = atoi(argv[i]);
                break;
            /* Select number of sampling size k for SGD */
            case 'k':
                k = atoi(argv[i]);
                break;
            /* Verbose mode: -v {0,1} (Default: 1) */
            case 'v':
                verbose = atoi(argv[i]);
                break;
            /* Unknown option */
            default:
                cout << "### Error: unknown option: -" << argv[i-1][1] << endl;
                exitWithHelp();
                break;
        }
    }

    if (i>=argc)
        exitWithHelp();

    // -- Determine filenames --
    strcpy(mean_vectors_filename        , argv[i]);
    strcpy(labels_filename              , argv[i+1]);
    strcpy(covariance_matrices_filename , argv[i+2]);
    strcpy(model_filename               , argv[i+3]);

}




void SvmGsuParams::exitWithHelp()
{
    cout << ("Usage: gsvm-train [options] <mean_vectors> <ground_truth> <covariance_matrices> <model_file>\n"
             "Options:\n"
             "\t+v <verbose_mode>: Verbose mode (default: 0)\n"
             "\t+t <kernel_type>: Set type of kernel function (default 0)\n"
             "\t\t0 -- Linear kernel\n"
             "\t\t2 -- Radial Basis Function (RBF) kernel\n"
             "\t+d <cov_mat>: Select covariance matrices type (default: 0)\n"
             "\t\t0 -- Full covariance matrices\n"
             "\t\t1 -- Diagonal covariance matrices\n"
             "\t\t3 -- Isotropic covariance matrices\n"
             "\t+l <lambda>: Set the parameter lambda of SVM-GSU (default 1.0)\n"
             "\t+g <gamma>: Set the parameter gamma (default 1/dim)\n"
             "\t+T <iter>: Set number of SGD iterations\n"
             "\t+k <k>: Set SGD sampling size\n");
    throw std::runtime_error("exit with help");
}




