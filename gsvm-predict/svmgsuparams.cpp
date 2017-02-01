#include "svmgsuparams.h"
#include "svmgsumodel.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>

using namespace std;


/* Ctor */
SvmGsuParams::SvmGsuParams(const SvmGsuModel& mod)
    : model(mod),
      kernel_type(0),
      lambda(1.0),
      gamma(0.1),
      sigmA(-1.0),
      sigmB(0.0),
      verbose(0)
{}


/* Getters */
int SvmGsuParams::getKernelType(){return kernel_type;}
double SvmGsuParams::getLambda(){return lambda;}
double SvmGsuParams::getGamma(){return gamma;}
double SvmGsuParams::getSigmA(){return sigmA;}
double SvmGsuParams::getSigmB(){return sigmB;}
int SvmGsuParams::getVerbose(){return verbose;}

/*  */
void SvmGsuParams::parseCommandLine(int argc, char** argv,
                                    char* mean_vectors_filename,
                                    char* labels_filename,
                                    char* model_filename,
                                    char* output_filename,
                                    char* evaluation_metrics_filename)
{
    int i;
    // -- Parse options --
    for (i=1; i<argc; i++)
    {
        if (argv[i][0] != '+')
            break;
        if (++i >= argc)
            exitWithHelp();

        switch (argv[i-1][1])
        {
            /* Select ground truth filename */
            case 't':
                strcpy(labels_filename, argv[i]);
                break;
            /* Select evaluation metrics filename */
            case 'm':
                strcpy(evaluation_metrics_filename, argv[i]);
                break;
            /* Verbose mode: -v {0,1} (Default: 1)  */
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
    strcpy(mean_vectors_filename, argv[i]);
    strcpy(model_filename, argv[i+1]);
    strcpy(output_filename, argv[i+2]);
}


/*  */
void SvmGsuParams::readModelFile(char* model_filename)
{
    /* Aux. variables */
    string line;
    vector<string> tokens;
    vector<string> elem;
    string first_word;
    bool read_svs = false;
    int dim = 0;
    int sv_cnt = 0;
    Eigen::VectorXd alpha;
    Eigen::MatrixXd SV;
    Eigen::VectorXd cur_sv;

    ifstream model_file( model_filename, ifstream::in );
    while ( getline(model_file,line) )
    {
        tokenize(line, tokens, " ");
        first_word = tokens[0];

        if (first_word.compare("kernel_type") == 0)
        {
            if (tokens[1].compare("linear")==0)
                kernel_type = 0;
            else if (tokens[1].compare("rbf")==0)
                kernel_type = 2;
            else if (tokens[1].compare("riemannian")==0)
                kernel_type = 3;
        }
        else if (first_word.compare("lambda") == 0)
        {
            lambda = stod(tokens[1]);
        }
        else if (first_word.compare("gamma") == 0)
        {
            gamma = stod(tokens[1]);
        }
        else if (first_word.compare("sigmA") == 0)
        {
            sigmA = stod(tokens[1]);
        }
        else if (first_word.compare("sigmB") == 0)
        {
            sigmB = stod(tokens[1]);
        }
        /**************************************************
         *               Read Linear model                *
         **************************************************/
        else if (first_word.compare("w") == 0)
        {
            dim = (int)tokens.size()-1;
            model.setDim(dim);
            Eigen::VectorXd temp_w(dim);
            temp_w.setZero(dim);
            for (int j=0; j<dim; j++)
                temp_w[j] = stod(tokens[j+1]);
            model.setW(temp_w);
        }
        else if (first_word.compare("b") == 0)
        {
            model.setB( stod( tokens[1] ) );
        }
        /**************************************************
         *               Read Kernel model                *
         **************************************************/
        else if (first_word.compare("SV") == 0)
        {
            read_svs = true;
            continue;
        }

        // Start reading SVs ...
        if (read_svs)
        {
            sv_cnt++;
            dim = (int)tokens.size()-1;
            alpha.conservativeResize(sv_cnt);
            alpha(sv_cnt-1) = stod(tokens[0]);

            cur_sv.resize(dim);
            cur_sv.setZero(dim);
            for (int j=1; j<dim+1; j++)
            {
                tokenize(tokens[j], elem, ":");
                int idx = stoi(elem[0]);
                double val = stod(elem[1]);
                cur_sv(idx-1) = val;
            }
            SV.conservativeResize(sv_cnt,dim);
            SV.row(sv_cnt-1) = cur_sv;
        }
    }
    if (read_svs)
    {
        model.setAlpha(alpha);
        model.setSV(SV);
    }
    model_file.close();
}


/*  */
void SvmGsuParams::tokenize(const string& str, vector<string>& tokens, const string& delimiters)
{
    tokens.clear();
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}


/*  */
void SvmGsuParams::exitWithHelp()
{

    cout << ("********************************************************************************\n"
             "* LIBSVMGSU: A framework for training/testing the SVM-GSU                      *\n"
             "*            (Support Vector Machine with Gaussian Sample Uncertainty).        *\n"
             "*                                                                              *\n"
             "* ** gsvm-predict **                                                           *\n"
             "*                                                                              *\n"
             "*                                                                              *\n"
             "* Version : 1.0 (February 2017)                                                *\n"
             "* Author  : Christos Tzelepis                                                  *\n"
             "* Contact : tzelepis@iti.gr                                                    *\n"
             "*           c.tzelepis@qmul.ac.uk                                              *\n"
             "********************************************************************************\n"
             "Usage: gsvm-predict [options] <mean_vectors> <model_file> <output_file>\n"
             "Options:\n"
             "\t-v <verbose_mode>       : Verbose mode (default: 0)\n"
             "\t-t <ground_truth>       : Select ground truth file\n"
             "\t-m <evaluation_metrics> : Evaluation metrics output file\n");
    throw std::runtime_error("exit with help");
}
