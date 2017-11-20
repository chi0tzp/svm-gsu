#include "svmgsuproblem.h"
#include "svmgsuparams.h"
#include "outputscore.h"
#include "meanvector.h"
#include "fullcovariancematrix.h"
#include "diagcovariancematrix.h"
#include "isocovariancematrix.h"
#include "label.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>

using namespace std;

static double str2double(const string& text)
{
    istringstream ss(text);
    double result;
    return ss >> result ? result : 0;
}

/* Ctor */
SvmGsuProblem::SvmGsuProblem(const SvmGsuParams& par)
    : params(par)
{}

/* Getters */
int SvmGsuProblem::getL(){return l;}
int SvmGsuProblem::getProbDim(){return prob_dim;}
vector<double> SvmGsuProblem::getLabels(){return y;}


/*******************************************************************************
 *                                                                             *
 *                            [ READ INPUT DATA ]                              *
 *                                                                             *
 *******************************************************************************/

/* Read input data (mean vectors and truth labels) */
void SvmGsuProblem::readInputData(char* mean_vectors_filename,
                                  char* labels_filename)
{
    /* Define vectors of labels and mean vectors */
    vector<Label> labels;
    vector<MeanVector> mean_vectors;

    /* Aux. variables */
    string line;
    vector<string> tokens;
    string doc_id_tmp;
    vector<string> elem;
    vector<string> labels_doc_ids;
    vector<string> mean_vectors_doc_ids;

    /*******************************************************************************
    *                              Read labels file                               *
    *******************************************************************************/
    ifstream labels_file( labels_filename, ifstream::in );
    Label cur_label;
    double label_tmp;
    /* Fill labels */
    while (getline(labels_file, line))
    {
        tokenize(line, tokens, " ");
        doc_id_tmp = tokens[0];
        label_tmp = str2double(tokens[1]);
        cur_label.setDocId(doc_id_tmp);
        cur_label.setY(label_tmp);
        labels.push_back(cur_label);
        labels_doc_ids.push_back(doc_id_tmp);
    }
    labels_file.close();

    /*******************************************************************************
     *                           Read mean vectors file                            *
     *******************************************************************************/
    ifstream mean_vectors_file(mean_vectors_filename, ifstream::in);
    vector<int> dims;
    int dim = 0;
    /* Find dimensionality of input space */
    while (getline(mean_vectors_file, line))
    {
        tokenize(line, tokens, " ");
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            dims.push_back(stoi(elem[0]));
        }
    }
    auto max_dim = max_element(dims.begin(), dims.end());
    dim = *max_dim;

    /* Fill mean_vectors */
    mean_vectors_file.clear();
    mean_vectors_file.seekg(0, ios::beg);
    while (getline(mean_vectors_file, line))
    {
        tokenize(line, tokens, " ");
        doc_id_tmp = tokens[0];
        MeanVector cur_mean_vector(dim, doc_id_tmp);
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            cur_mean_vector.setXj(stoi(elem[0])-1, str2double(elem[1]));
        }
        mean_vectors.push_back(cur_mean_vector);
        mean_vectors_doc_ids.push_back(doc_id_tmp);
    }
    mean_vectors_file.close();

    /* Sort labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids by doc_id */
    sort(labels_doc_ids.begin(), labels_doc_ids.end());
    sort(mean_vectors_doc_ids.begin(), mean_vectors_doc_ids.end());

    /* Find intersection between labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids based on doc_id */
    vector<string> labels_mean_intersect = findDocIdsIntersection(labels_doc_ids, mean_vectors_doc_ids);

    /* Initialize SVM-GSU problem */
    l = labels_mean_intersect.size();

    /* Populate data */
    string id;
    for (int i=0; i<l; i++)
    {
        id = labels_mean_intersect[i];
        // Get label
        for (unsigned j=0; j<labels.size(); j++)
        {
            if (labels[j].getDocId().compare(id) == 0)
            {
                doc_id.push_back(labels[j].getDocId());
                y.push_back(labels[j].getY());
                break;
            }
        }
        // Get mean vector
        for (unsigned j=0; j<mean_vectors.size(); j++)
        {
            if (mean_vectors[j].getDocId().compare(id) == 0)
            {
                x.push_back(mean_vectors[j].getX());
                break;
            }
        }
    }

    /* Set problem's dimension */
    prob_dim = dim;

    /* HERE */
    if (params.getKernelType() == 0)
    {
        int model_dim = params.model.getModelDim();

        if (model_dim != prob_dim)
        {
            Eigen::VectorXd w = params.model.getW();
            w.conservativeResizeLike(Eigen::VectorXd::Zero(prob_dim));
            params.model.setW(w);
        }
    }

}


/*  */
vector<OutputScore> SvmGsuProblem::predict()
{
    Eigen::VectorXd w = params.model.getW();
    Eigen::VectorXd alpha = params.model.getAlpha();
    Eigen::MatrixXd SV = params.model.getSV();
    double b = params.model.getB();
    double sigmA = params.getSigmA();
    double sigmB = params.getSigmB();
    double fun_val = 0.0;
    double score = 0.0;
    int label = 0;

    vector<OutputScore> output;

    if (params.getKernelType() == 0)
    {
        for (int i=0; i<l; i++)
        {
            fun_val = w.transpose()*x[i] + b;
            score = 1.0/(1.0+exp(sigmA*fun_val+sigmB));
            label = sgn(fun_val);
            OutputScore out(doc_id[i], score, label, y[i]);
            output.push_back(out);
        }
    }
    else if (params.getKernelType() == 2)
    {
        for (int i=0; i<l; i++)
        {
            fun_val = b;
            Eigen::VectorXd diff;
            for (int j=0; j<SV.rows(); j++)
            {
                diff = x[i].transpose() - SV.row(j);
                fun_val += alpha(j) * exp(-params.getGamma()*diff.squaredNorm());
            }
            score = 1.0/(1.0+exp(sigmA*fun_val+sigmB));
            label = sgn(fun_val);
            OutputScore out(doc_id[i], score, label, y[i]);
            output.push_back(out);
        }
    }

    sort(output.begin(), output.end(), [](const OutputScore & r1, const OutputScore & r2) -> bool {
                                               return r1.score > r2.score; });
    return output;
}


/*******************************************************************************
 *                                                                             *
 *                           [ Write Output File ]                             *
 *                                                                             *
 *******************************************************************************/
/* Write output file */
void SvmGsuProblem::writeOutputFile(vector<OutputScore> output, char *output_filename)
{
    ofstream f;
    f.open (output_filename);
    f << "doc_id +1 -1 label" << endl;
    for (unsigned i=0; i<output.size(); i++)
        f << output[i].doc_id << " " << output[i].score << " " << (1.0-output[i].score) << " " << output[i].predicted_label << endl;
    f.close();
}


/*******************************************************************************
 *                                                                             *
 *                          [ Auxiliary Functions ]                            *
 *                                                                             *
 *******************************************************************************/
/*  */
int SvmGsuProblem::sgn(double a){return (a<0.0)?-1:+1;}


/*  */
void SvmGsuProblem::tokenize(const string& str, vector<string>& tokens, const string& delimiters)
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
vector<string> SvmGsuProblem::findDocIdsIntersection(vector<string> A, vector<string> B)
{
    vector<string> intersection;
    int n1 = A.size();
    int n2 = B.size();
    int i = 0;
    int j = 0;

    while(i<n1 && j<n2)
    {
        if (A[i].compare(B[j]) > 0)
            j++;
        else if (B[j].compare(A[i]) > 0)
            i++;
        else
        {
            intersection.push_back(A[i]);
            i++;
            j++;
        }
    }
    return intersection;
}
