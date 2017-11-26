#include <iostream>
#include "svmgsumodel.h"
#include "svmgsuparams.h"
#include "svmgsuproblem.h"
#include "evaluationmetrics.h"
#include <chrono>


#define MAX_FILENAME_LEN 4096

using namespace std;
using namespace std::chrono;

void getElapsedTime(int);


int main(int argc, char *argv[])
{
    char mean_vectors_filename[MAX_FILENAME_LEN];
    char labels_filename[MAX_FILENAME_LEN];
    char model_filename[MAX_FILENAME_LEN];
    char output_filename[MAX_FILENAME_LEN];
    char evaluation_metrics_filename[MAX_FILENAME_LEN];

    high_resolution_clock::time_point read_data_start;
    high_resolution_clock::time_point read_data_end;
    high_resolution_clock::time_point predict_start;
    high_resolution_clock::time_point predict_end;
    high_resolution_clock::time_point write_output_start;
    high_resolution_clock::time_point write_output_end;
    high_resolution_clock::time_point comp_metrics_start;
    high_resolution_clock::time_point comp_metrics_end;


    /* Define SVM-GSU model */
    SvmGsuModel model;

    /* Define SVM-GSU parameters and parse command line arguments */
    SvmGsuParams params(model);
    try
    {
        params.parseCommandLine(argc, argv,
                                mean_vectors_filename,
                                labels_filename,
                                model_filename,
                                output_filename,
                                evaluation_metrics_filename);
    }
    catch(...)
    {
        return 0;
    }

    if (params.getVerbose() == 1)
    {
        cout << ("************************************************************************\n"
                 "* svm-gsu: A framework for training/testing the Support Vector Machine *\n"
                 "*          with Gaussian Sample Uncertainty (SVM-GSU).                 *\n"
                 "*  -- gsvm-predict: Evaluate a trained SVM-GSU model.                  *\n"
                 "*----------------------------------------------------------------------*\n"
                 "* Version : 0.1                                                        *\n"
                 "* Author  : Christos Tzelepis                                          *\n"
                 "* Contact : tzelepis@iti.gr                                            *\n"
                 "* GitHub  : @chi0tzp                                                   *\n"
                 "************************************************************************\n");
    }

    /* Read model file */
    params.readModelFile(model_filename);


    /* Define SVM-GSU problem */
    SvmGsuProblem prob(params);

    /* Read input data */
    if (params.getVerbose()){
        read_data_start = high_resolution_clock::now();
        cout << ".Read input data...";
    }

    prob.readInputData(mean_vectors_filename, labels_filename);

    if (params.getVerbose()){
        read_data_end = high_resolution_clock::now();
        cout << "Done! [Elapsed time: ";
        getElapsedTime( duration_cast<seconds>(read_data_end-read_data_start).count() );
        cout << "]" << endl;
        cout << " \\__l = " << prob.getL() << endl;
        cout << " \\__dim = " << prob.getProbDim() << endl;
    }

    /* Predict */
    if (params.getVerbose()){
        predict_start = high_resolution_clock::now();
        cout << ".Predict...";
    }

    vector<OutputScore> output = prob.predict();

    if (params.getVerbose()){
        predict_end = high_resolution_clock::now();
        cout << "Done! [Elapsed time: ";
        getElapsedTime(duration_cast<seconds>(predict_end-predict_start).count());
        cout << "]" << endl;
    }

    /* Write output file */
    if (params.getVerbose()){
        write_output_start = high_resolution_clock::now();
        cout << ".Write output file...";
    }

    prob.writeOutputFile(output, output_filename);

    if (params.getVerbose()){
        write_output_end = high_resolution_clock::now();
        cout << "Done! [Elapsed time: ";
        getElapsedTime(duration_cast<seconds>(write_output_end-write_output_start).count());
        cout << "]" << endl;
    }

    /* Compute and write evaluation metrics */
    if (params.getVerbose()){
        comp_metrics_start = high_resolution_clock::now();
        cout << ".Compute and write evaluation metrics file...";
    }

    EvaluationMetrics eval(output);
    eval.writeEvaluationMetrics(evaluation_metrics_filename);

    if (params.getVerbose()){
        comp_metrics_end = high_resolution_clock::now();
        cout << "Done! [Elapsed time: ";
        getElapsedTime(duration_cast<seconds>(comp_metrics_end-comp_metrics_start).count());
        cout << "]" << endl;
        cout << " \\__Accuracy          : " << eval.getAcc() << endl;
        cout << " \\__Average Precision : " << eval.getAP() << endl;
        cout << " \\__Precision         : " << eval.getPrec() << endl;
        cout << " \\__Recall            : " << eval.getRecall() << endl;
    }

    return 0;
}


void getElapsedTime(int seconds)
{
    int hours, minutes;
    minutes = seconds / 60;
    hours   = minutes / 60;
    cout << int(hours) << "H:" << int(minutes%60) << "M:" << int(seconds%60) << "S";
}


