#include <iostream>
#include <chrono>
#include "svmgsuparams.h"
#include "svmgsuproblem.h"

#define MAX_FILENAME_LEN 4096

using namespace std;
using namespace std::chrono;

void getElapsedTime(int);

int main(int argc, char* argv[])
{
    /* Define filenames */
    char mean_vectors_filename[MAX_FILENAME_LEN];
    char labels_filename[MAX_FILENAME_LEN];
    char covariance_matrices_filename[MAX_FILENAME_LEN];
    char model_filename[MAX_FILENAME_LEN];

    high_resolution_clock::time_point read_data_start;
    high_resolution_clock::time_point read_data_end;

    high_resolution_clock::time_point solve_prob_start;
    high_resolution_clock::time_point solve_prob_end;

    high_resolution_clock::time_point eigdecomp_start;
    high_resolution_clock::time_point eigdecomp_end;

    high_resolution_clock::time_point data_projection_start;
    high_resolution_clock::time_point data_projection_end;

    high_resolution_clock::time_point sgd_start;
    high_resolution_clock::time_point sgd_end;


    /* Define SVM-GSU parameters and parse command line arguments */
    SvmGsuParams params;
    try
    {
        params.parseCommandLine(argc, argv,
                                mean_vectors_filename,
                                labels_filename,
                                covariance_matrices_filename,
                                model_filename);
    }
    catch(...)
    {
        return 0;
    }

    if (params.getVerbose() == 1){

        cout << ("************************************************************************\n"
                 "* svm-gsu: A framework for training/testing the Support Vector Machine *\n"
                 "*          with Gaussian Sample Uncertainty (SVM-GSU).                 *\n"
                 "*  -- gsvm-train: Train an SVM-GSU model.                              *\n"
                 "*----------------------------------------------------------------------*\n"
                 "* Version : 0.1                                                        *\n"
                 "* Author  : Christos Tzelepis                                          *\n"
                 "* Contact : tzelepis@iti.gr                                            *\n"
                 "* GitHub  : @chi0tzp                                                   *\n"
                 "************************************************************************\n");
        cout <<  ".Options selected:" << endl;
        cout <<  " \\__Kernel type: " << params.getKernelType() << endl;
        cout <<  " \\__lambda parameter: " << params.getLambda()     << endl;
        cout <<  " \\__Covariance matrices type: " << params.getCovMatTypeVerb() << endl;
        cout <<  " \\__Number of SGD iterations: T = " << params.getT() << endl;
        cout <<  " \\__SGD sampling size: k = " << params.getK() << endl;
        if (params.getP()<1.0)
            cout << " \\__Learning in linear subspaces: p=" << params.getP() << endl;
        else
            cout << " \\__Learning in the original feature space" << endl;
    }

    /* Define the SVM-GSU problem */
    SvmGsuProblem prob(params);


    /*******************************************************************************
     *                                                                             *
     *                             [ Read input data ]                             *
     *                                                                             *
     *******************************************************************************/
    if (params.getVerbose() ==  1){
        cout << ".Read input data...";
        read_data_start = high_resolution_clock::now();
    }

    switch (params.getCovMatType()){
        case 0:
            prob.readInputDataFull(mean_vectors_filename, covariance_matrices_filename, labels_filename);
            break;
        case 1:
            prob.readInputDataDiag(mean_vectors_filename, covariance_matrices_filename, labels_filename);
            break;
        case 2:
            prob.readInputDataIso(mean_vectors_filename, covariance_matrices_filename, labels_filename);
            break;
        default:
            break;
    }

    if (params.getVerbose() == 1){
        read_data_end = high_resolution_clock::now();
        cout << "Done! [Elapsed time: ";
        getElapsedTime( duration_cast<seconds>(read_data_end-read_data_start).count() );
        cout << "]" << endl;
        cout << " \\__Training set cardinality: " << prob.getL() << endl;
        cout << " \\__Feature space dimensionality: " << prob.getDim() << endl;
    }


    /*******************************************************************************
     *                                                                             *
     *                          [ Solve SVM-GSU problem ]                          *
     *                                                                             *
     *******************************************************************************/

    if (params.getVerbose() == 1){
        cout << ".Solve problem..." << endl;
        solve_prob_start = high_resolution_clock::now();
    }

    /*******************************************************************************
     *                          Learn in linear subspaces                          *
     *******************************************************************************/
    if (params.getP()<1.0){
        /* Conduct eigenanalysis on the input covariance matrices */
        if (params.getVerbose() == 1){
            cout << " \\__.Eigendecomposition...";
            eigdecomp_start = high_resolution_clock::now();
        }

        prob.eigenDecomp();

        if (params.getVerbose() == 1){
            eigdecomp_end = high_resolution_clock::now();
            cout << "Done! [Elapsed time: ";
            getElapsedTime(duration_cast<seconds>(eigdecomp_end - eigdecomp_start).count());
            cout << "]" << endl;

            vector<int> ss_dims;
            for (int i=0; i<prob.getL(); i++)
                ss_dims.push_back(prob.getPrincAxes(i).size());
            int min_dim = *std::min_element(ss_dims.begin(), ss_dims.end());
            int max_dim = *std::max_element(ss_dims.begin(), ss_dims.end());
            cout << " |   \\__min dim: " << min_dim << endl;
            cout << " |   \\__max dim: " << max_dim << endl;
        }

        /* Solve the problem using SGD */
        if (params.getVerbose() == 1){
            sgd_start = high_resolution_clock::now();
            cout << " \\__SGD...";
        }

        if (params.getKernelType() == 0)
            prob.solveLSVMGSUzSpace();
        else
            // Not Available/Implemented Yet (?)


        if (params.getVerbose() == 1){
            sgd_end = high_resolution_clock::now();
            cout << "Done! [Elapsed time: ";
            getElapsedTime(duration_cast<seconds>(sgd_end - sgd_start).count());
            cout << "]" << endl;
        }
    }
    /*******************************************************************************
     *                         Learn in the original space                         *
     *******************************************************************************/
    else{
        if (params.getVerbose() == 1){
            cout << " \\__SGD...";
            sgd_start = high_resolution_clock::now();
        }

        /* Solve the problem using SGD */
        if (params.getKernelType() == 0){
            prob.solveLSVMGSUxSpace();
        }
        else if (params.getKernelType() == 2){
            // TODO:
            prob.computeKernelMatrix();
            prob.solveKSVMiGSU();
        }


        if (params.getVerbose() == 1){
            sgd_end = high_resolution_clock::now();
            cout << "Done! [Elapsed time: ";
            getElapsedTime(duration_cast<seconds>(sgd_end - sgd_start).count());
            cout << "]" << endl;
        }
    }

    if (params.getVerbose()){
        solve_prob_end = high_resolution_clock::now();
        cout << ".Problem Solved! [Elapsed time: ";
        getElapsedTime(duration_cast<seconds>(solve_prob_end - solve_prob_start).count());
        cout << "]" << endl;
    }

    /* Write model file */

    if (params.getVerbose())
        cout << ".Write model file." << endl;

    if (params.getKernelType() == 0)
        prob.writeLinearModelFile(model_filename);
    else if (params.getKernelType() == 2)
        prob.writeKernelModelFile(model_filename);

    if (params.getVerbose())
        cout << "*** Goodbye! ***" << endl;


    return 0;
}


void getElapsedTime(int seconds)
{
    int hours, minutes;
    minutes = seconds / 60;
    hours   = minutes / 60;
    cout << int(hours) << "H:" << int(minutes%60) << "M:" << int(seconds%60) << "S";
}
