#include "svmgsuproblem.h"
#include "svmgsuparams.h"
#include "meanvector.h"
#include "label.h"
#include "fullcovariancematrix.h"
#include "diagcovariancematrix.h"
#include "isocovariancematrix.h"
#include <algorithm>
#include <random>
#include <string>
#include <eigen3/Eigen/Dense>
#include <numeric>
#include <iostream>
#include <exception>
#include <fstream>
#include <chrono>
#include <math.h>
#include <sstream>


double str2double(const string& text)
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
int SvmGsuProblem::getDim(){return dim;}
vector<int> SvmGsuProblem::getPrincAxes(int i){return princ_axes[i];}



/*******************************************************************************
 *                                                                             *
 *                            [ READ INPUT DATA ]                              *
 *                                                                             *
 *******************************************************************************/

/* Read input data (mean vectors, labels, full covariance matrices) */
void SvmGsuProblem::readInputDataFull(char* mean_vectors_filename,
                                      char* covariance_matrices_filename,
                                      char* labels_filename)
{
    /* Define vectors of labels, mean vectors, and covariance matrices */
    vector<Label> labels;
    vector<MeanVector> mean_vectors;
    vector<FullCovarianceMatrix> covariance_matrices;

    /* Aux. variables */
    string line;
    vector<string> tokens;
    string doc_id_tmp;
    vector<string> elem;
    vector<string> indices;
    vector<string> labels_doc_ids;
    vector<string> mean_vectors_doc_ids;
    vector<string> covariance_matrices_doc_ids;


    /*******************************************************************************
     *                              Read labels file                               *
     *******************************************************************************/
    ifstream labels_file(labels_filename, ifstream::in);
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
     *                       Find Input space dimensionality                       *
     *******************************************************************************/
    vector<int> dims;

    /* Mean vectors */
    ifstream mean_vectors_file(mean_vectors_filename, ifstream::in);
    while (getline(mean_vectors_file, line))
    {
        tokenize(line, tokens, " ");
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            dims.push_back(stoi(elem[0]));
        }
    }

    /* Covariance matrices */
    ifstream covariance_matrices_file(covariance_matrices_filename, ifstream::in);
    while (getline(covariance_matrices_file, line))
    {
        tokenize(line, tokens, " ");
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            tokenize(elem[0], indices, ",");
            int row_idx = stoi(indices[0]);
            int col_idx = stoi(indices[1]);
            dims.push_back(std::max(row_idx, col_idx));
        }
    }
    auto max_dim = max_element(dims.begin(), dims.end());
    dim = *max_dim;


    /*******************************************************************************
     *                           Read mean vectors file                            *
     *******************************************************************************/
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
        mean_vectors.push_back( cur_mean_vector );
        mean_vectors_doc_ids.push_back( doc_id_tmp );
    }

    /*******************************************************************************
     *                     Read full covariance matrices file                      *
     *******************************************************************************/
    covariance_matrices_file.clear();
    covariance_matrices_file.seekg(0, ios::beg);
    /* Fill covariance_matrices */
    while (getline(covariance_matrices_file, line))
    {
        tokenize(line, tokens, " ");
        doc_id_tmp = tokens[0];
        FullCovarianceMatrix cur_covariance_matrix(dim, doc_id_tmp);
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            tokenize(elem[0], indices, ",");
            cur_covariance_matrix.setSigma(stoi(indices[0])-1, stoi(indices[1])-1, str2double(elem[1]));
        }
        covariance_matrices.push_back(cur_covariance_matrix);
        covariance_matrices_doc_ids.push_back(doc_id_tmp);
    }
    mean_vectors_file.close();
    covariance_matrices_file.close();

    /* Sort labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids by doc_id */
    sort(labels_doc_ids.begin(), labels_doc_ids.end());
    sort(mean_vectors_doc_ids.begin(), mean_vectors_doc_ids.end());
    sort(covariance_matrices_doc_ids.begin(), covariance_matrices_doc_ids.end());

    /* Find intersection between labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids based on doc_id */
    vector<string> labels_mean_intersect = findDocIdsIntersection(labels_doc_ids, mean_vectors_doc_ids);
    vector<string> labels_mean_cov_intersect = findDocIdsIntersection(labels_mean_intersect, covariance_matrices_doc_ids);

    /* Initialize SVM-GSU problem */
    l = labels_mean_cov_intersect.size();

    /* Populate data */
    string doc_id;
    for (int i=0; i<l; i++)
    {
        doc_id = labels_mean_cov_intersect[i];
        // Get label
        for (unsigned j=0; j<labels.size(); j++)
        {
            if (labels[j].getDocId().compare(doc_id) == 0)
            {
                y.push_back(getLabel(labels[j].getY()));
                break;
            }
        }
        // Get mean vector
        for (unsigned j=0; j<mean_vectors.size(); j++)
        {
            if (mean_vectors[j].getDocId().compare(doc_id) == 0)
            {
                x.push_back(mean_vectors[j].getX());
                break;
            }
        }
        // Get covariance matrix
        for (unsigned j=0; j<covariance_matrices.size(); j++)
        {
            if (covariance_matrices[j].getDocId().compare(doc_id) == 0)
            {
                Sigma_xf.push_back(covariance_matrices[j].getSigma());
                break;
            }
        }
    }
}


/* Read input data (mean vectors, labels, diagonal covariance matrices) */
void SvmGsuProblem::readInputDataDiag(char* mean_vectors_filename,
                                      char* covariance_matrices_filename,
                                      char* labels_filename)
{
    /* Define vectors of labels, mean vectors, and (diagonal) covariance matrices */
    vector<Label> labels;
    vector<MeanVector> mean_vectors;
    vector<DiagCovarianceMatrix> covariance_matrices;

    /* Aux. variables */
    string line;
    vector<string> tokens;
    string doc_id_tmp;
    vector<string> elem;
    vector<string> indices;
    vector<string> labels_doc_ids;
    vector<string> mean_vectors_doc_ids;
    vector<string> covariance_matrices_doc_ids;


    /*******************************************************************************
     *                              Read labels file                               *
     *******************************************************************************/
    ifstream labels_file( labels_filename, ifstream::in );
    if (!labels_file.is_open())
        throw std::runtime_error("file not found");

    Label cur_label;
    double label_tmp;
    /* Fill labels */
    while (getline(labels_file,line))
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
     *                       Find Input space dimensionality                       *
     *******************************************************************************/
    vector<int> dims;

    /* Mean vectors */
    ifstream mean_vectors_file(mean_vectors_filename, ifstream::in);
    while (getline(mean_vectors_file, line))
    {
        tokenize(line, tokens, " ");
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            dims.push_back(stoi(elem[0]));
        }
    }

    /* Covariance matrices */
    ifstream covariance_matrices_file(covariance_matrices_filename, ifstream::in);
    while (getline(covariance_matrices_file, line))
    {
        tokenize(line, tokens, " ");
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            tokenize(elem[0], indices, ",");
            if (stoi(indices[0]) == stoi(indices[1]))
                dims.push_back(stoi(indices[0]));
            else
            {
                // TODO:
            }
        }
    }
    auto max_dim = max_element(dims.begin(), dims.end());
    dim = *max_dim;


    /*******************************************************************************
     *                           Read mean vectors file                            *
     *******************************************************************************/
    mean_vectors_file.clear();
    mean_vectors_file.seekg(0, ios::beg);
    /* Fill mean_vectors */
    while (getline(mean_vectors_file,line))
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


    /*******************************************************************************
     *                      Read covariance matrices file                          *
     *******************************************************************************/
    covariance_matrices_file.clear();
    covariance_matrices_file.seekg(0, ios::beg);
    /* Fill covariance_matrices */
    while (getline(covariance_matrices_file, line))
    {
        tokenize(line, tokens, " ");
        doc_id_tmp = tokens[0];
        DiagCovarianceMatrix cur_covariance_matrix(dim, doc_id_tmp);
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            tokenize(elem[0], indices, ",");
            cur_covariance_matrix.setSigma(stoi(indices[0])-1, str2double(elem[1]));
        }
        covariance_matrices.push_back(cur_covariance_matrix);
        covariance_matrices_doc_ids.push_back(doc_id_tmp);
    }
    mean_vectors_file.close();
    covariance_matrices_file.close();

    /* Sort labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids by doc_id */
    sort(labels_doc_ids.begin(), labels_doc_ids.end());
    sort(mean_vectors_doc_ids.begin(), mean_vectors_doc_ids.end());
    sort(covariance_matrices_doc_ids.begin(), covariance_matrices_doc_ids.end());

    /* Find intersection between labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids based on doc_id */
    vector<string> labels_mean_intersect = findDocIdsIntersection(labels_doc_ids, mean_vectors_doc_ids);
    vector<string> labels_mean_cov_intersect = findDocIdsIntersection(labels_mean_intersect, covariance_matrices_doc_ids);

    /* Initialize SVM-GSU problem */
    l = labels_mean_cov_intersect.size();

    /* Populate data */
    string doc_id;
    for (int i=0; i<l; i++)
    {
        doc_id = labels_mean_cov_intersect[i];
        // Get label
        for (unsigned j=0; j<labels.size(); j++)
        {
            if (labels[j].getDocId().compare(doc_id) == 0)
            {
                y.push_back(getLabel(labels[j].getY()));
                break;
            }
        }
        // Get mean vector
        for (unsigned j=0; j<mean_vectors.size(); j++)
        {
            if (mean_vectors[j].getDocId().compare(doc_id) == 0)
            {
                x.push_back(mean_vectors[j].getX());
                break;
            }
        }
        // Get covariance matrix
        for (unsigned j=0; j<covariance_matrices.size(); j++)
        {
            if (covariance_matrices[j].getDocId().compare(doc_id) == 0)
            {
                Sigma_xd.push_back(covariance_matrices[j].getSigma());
                break;
            }
        }
    }
}


/* Read input data (mean vectors, labels, isotropical covariance matrices) */
void SvmGsuProblem::readInputDataIso(char* mean_vectors_filename,
                                     char* covariance_matrices_filename,
                                     char* labels_filename)
{
    /* Define vectors of labels, mean vectors, and (isotropic) covariance matrices */
    vector<Label> labels;
    vector<MeanVector> mean_vectors;
    vector<IsoCovarianceMatrix>  covariance_matrices;

    /* Aux. variables */
    string line;
    vector<string> tokens;
    string doc_id_tmp;
    vector<string> elem;
    vector<string> indices;
    vector<string> labels_doc_ids;
    vector<string> mean_vectors_doc_ids;
    vector<string> covariance_matrices_doc_ids;


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


    /*******************************************************************************
     *                   Read isotropic covariance matrices file                   *
     *******************************************************************************/
    ifstream covariance_matrices_file(covariance_matrices_filename, ifstream::in);

    /* Fill covariance_matrices */
    while (getline(covariance_matrices_file,line))
    {
        tokenize(line, tokens, " ");
        doc_id_tmp = tokens[0];
        IsoCovarianceMatrix cur_covariance_matrix(doc_id_tmp);
        for (unsigned t=1; t<tokens.size(); t++)
        {
            tokenize(tokens[t], elem, ":");
            tokenize(elem[0], indices, ",");
            if (stoi(indices[0]) != stoi(indices[1]))
            {
                // TODO:
            }
            cur_covariance_matrix.setSigma(str2double(elem[1]));
        }
        covariance_matrices.push_back(cur_covariance_matrix);
        covariance_matrices_doc_ids.push_back(doc_id_tmp);
    }
    covariance_matrices_file.close();

    /* Sort labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids by doc_id */
    sort(labels_doc_ids.begin(), labels_doc_ids.end());
    sort(mean_vectors_doc_ids.begin(), mean_vectors_doc_ids.end());
    sort(covariance_matrices_doc_ids.begin(), covariance_matrices_doc_ids.end());

    /* Find intersection between labels_doc_ids, mean_vectors_doc_ids, covariance_matrices_doc_ids based on doc_id */
    vector<string> labels_mean_intersect = findDocIdsIntersection(labels_doc_ids, mean_vectors_doc_ids);
    vector<string> labels_mean_cov_intersect = findDocIdsIntersection(labels_mean_intersect, covariance_matrices_doc_ids);

    /* Initialize SVM-GSU problem */
    l = labels_mean_cov_intersect.size();

    /* Populate data */
    string doc_id;
    for (int i=0; i<l; i++)
    {
        doc_id = labels_mean_cov_intersect[i];
        // Get label
        for (unsigned j=0; j<labels.size(); j++)
        {
            if (labels[j].getDocId().compare(doc_id) == 0)
            {
                y.push_back(getLabel(labels[j].getY()));
                break;
            }
        }
        // Get mean vector
        for (unsigned j=0; j<mean_vectors.size(); j++)
        {
            if (mean_vectors[j].getDocId().compare(doc_id) == 0)
            {
                x.push_back(mean_vectors[j].getX());
                break;
            }
        }
        // Get covariance matrix
        for (unsigned j=0; j<covariance_matrices.size(); j++)
        {
            if (covariance_matrices[j].getDocId().compare(doc_id) == 0)
            {
                Sigma_xi.push_back(covariance_matrices[j].getSigma());
                break;
            }
        }
    }
}



/*******************************************************************************
 *                                                                             *
 *                           [ EIGENDECOMPOSITION ]                            *
 *                                                                             *
 *******************************************************************************/
void SvmGsuProblem::eigenDecomp()
{
    int i;
    switch (params.getCovMatType())
    {
        /*******************************************************************************
         *                          Full Covariance Matrices                           *
         *******************************************************************************/
        case 0:
            // =================== //
            // Not implemented yet //
            // =================== //
            break;

        /*******************************************************************************
         *                        Diagonal Covariance Matrices                         *
         *******************************************************************************/
        case 1:

            for (i=0; i<l; i++ )
            {
                // Get current covariance matrix (vector of variances)
                Eigen::VectorXd var_x = Sigma_xd[i];

                // Create current sorting indexes: sort_ids = (0, ..., dim-1)
                vector<int> sort_ids(dim);
                iota(sort_ids.begin(), sort_ids.end(), 0);

                // Sort indexes based on current var_xiances
                sort(sort_ids.begin(), sort_ids.end(),
                     [&var_x](size_t i1, size_t i2) {return var_x[i1] < var_x[i2];});

                // Sort current variances
                std::sort(var_x.data(), var_x.data() + var_x.size(),
                          std::greater<double>());

                int t = 0;
                double par_sum = 0.0;
                double tot_sum = var_x.sum();
                vector<int> cur_princ_axes;
                Eigen::VectorXd x_i = x[i];
                Eigen::VectorXd z_i;
                for (t=0; t<dim; t++)
                {
                    par_sum += var_x(t);
                    cur_princ_axes.push_back(sort_ids[t]);
                    z_i.conservativeResize(t+1);
                    z_i(t) = x_i(sort_ids[t]);
                    if (par_sum/tot_sum > params.getP())
                        break;
                }
                z.push_back(z_i);
                princ_axes.push_back(cur_princ_axes);

                Eigen::VectorXd var_z(t);
                var_z = var_x.block(0,0,t+1,1);
                Sigma_zd.push_back(var_z);

            }
            break;

        /*******************************************************************************
         *                       Isotropic Covariance Matrices                         *
         *******************************************************************************/
        case 2:
            // =================== //
            // Not implemented yet //
            // =================== //
            break;
        default:
            break;
    }
}



/*******************************************************************************
 *                                                                             *
 *                                [ SOLVERS ]                                  *
 *                                                                             *
 *******************************************************************************/

/* Solve the linear SVM-GSU in the original feature space */
void SvmGsuProblem::solveLSVMGSUxSpace()
{
    int T, k;

    /* Get number of iterations T, sampling size k */
    T = params.getT();
    k = params.getK();

    /* Initialize linear model w = (w_1, w_2, ..., w_dim)^T, b */
    initW(dim, params.getLambda());
    b = 0.0;

    vector<int> ids_t(l);
    iota(ids_t.begin(), ids_t.end(), 0);

    double eta_t = 1.0;
    Eigen::VectorXd sumdL(dim+1);
    sumdL.setZero(dim+1);

    // Define aux. optimization variables
    Eigen::VectorXd sumdLdw(dim);
    sumdLdw.setZero(dim);
    double sumdLdb = 0.0;

    switch (params.getCovMatType())
    {
        /* Full covariance matrices */
        case 0:
            for (int t=1; t<=T; t++)
            {
                // Get random subset of training examples' indices (ids_t)
                auto engine = std::default_random_engine{};
                std::shuffle(std::begin(ids_t), std::end(ids_t), engine);
                vector<int> ids_k(ids_t.begin(), ids_t.begin()+k);

                // Learning rate
                eta_t = 1.0/( params.getLambda()*(double)t );

                // Compute sum of gradients of loss function
                sumdL = computeSumOfLossGradFullXspace(w, b, ids_k);
                sumdLdw = sumdL.block(0,0,dim,1);
                sumdLdb = sumdL[dim];

                // Update optimization variables
                w = (1.0-1.0/(double(t)))*w - (eta_t/((double)k)) * sumdLdw;
                w = min( 1.0, 1.0/sqrt(params.getLambda()*w.squaredNorm()) ) * w;
                b = b - (eta_t/((double)k)) * sumdLdb;
            }
            break;

        /* Diagonal covariance matrices */
        case 1:
            for (int t=1; t<=T; t++)
            {
                // Get random subset of training examples' indices (ids_t)
                auto engine = std::default_random_engine{};
                std::shuffle(std::begin(ids_t), std::end(ids_t), engine);
                vector<int> ids_k(ids_t.begin(), ids_t.begin()+k);

                // Learning rate
                eta_t = 1.0/( params.getLambda()*(double)t );

                // Compute sum of gradients of loss function
                sumdL = computeSumOfLossGradDiagXspace(w, b, ids_k);
                sumdLdw = sumdL.block(0,0,dim,1);
                sumdLdb = sumdL[dim];

                // Update optimization variables
                w = (1.0-1.0/double(t))*w - (eta_t/((double)k)) * sumdLdw;
                w = min( 1.0, 1.0/sqrt(params.getLambda()*w.squaredNorm()) ) * w;
                b = b - (eta_t/((double)k)) * sumdLdb;
            }
            break;

        /* Isotropic covariance matrices */
        case 2:
            for (int t=1; t<=T; t++)
            {
                // Get random subset of training examples' indices (ids_t)
                auto engine = std::default_random_engine{};
                std::shuffle(std::begin(ids_t), std::end(ids_t), engine);
                vector<int> ids_k(ids_t.begin(), ids_t.begin()+k);

                // Learning rate
                eta_t = 1.0/( params.getLambda()*(double)t );

                // Compute sum of gradients of loss function
                sumdL = computeSumOfLossGradIsoXspace(w, b, ids_k);
                sumdLdw = sumdL.block(0,0,dim,1);
                sumdLdb = sumdL[dim];

                // Update optimization variables
                w = (1.0-1.0/double(t))*w - (eta_t/((double)k)) * sumdLdw;
                w = min( 1.0, 1.0/sqrt(params.getLambda()*w.squaredNorm()) ) * w;
                b = b - (eta_t/((double)k)) * sumdLdb;
            }
            break;

        default:
            break;
    }
    /* Compute decision values and apply Platt scaling*/
    computeDecValues();
    plattScaling();

}


/*  */
Eigen::VectorXd SvmGsuProblem::computeSumOfLossGradFullXspace(Eigen::VectorXd w,
                                                              double b,
                                                              vector<int> ids)
{
    Eigen::VectorXd sumdL(dim+1);
    Eigen::VectorXd sumdLdw(dim);
    double sumdLdb;

    sumdL.setZero(dim+1);
    sumdLdw.setZero(dim);
    sumdLdb = 0.0;

    double d_mu = 0.0;
    double d_sigma = 0.0;
    double r = 0.0;
    double erf_ = 0.0;
    double exp_ = 0.0;

    int i;
    for (unsigned t=0; t<ids.size(); t++ )
    {
        i = ids[t];
        d_mu = y[i] - w.transpose()*x[i] - b;
        d_sigma = w.transpose() * Sigma_xf[i] * w;
        r = d_mu / d_sigma;
        erf_ = erf(0.5*sqrt(2.0)*r);
        exp_ = exp(-0.5 * r * r);
        // Derivative of loss wrt w
        sumdLdw += (exp_ / (sqrt(2.0*M_PI) * d_sigma)) * Sigma_xf[i]*w - 0.5 * (erf_ + y[i]) * x[i];
        // Derivative of loss wrt b
        sumdLdb -= 0.5 * (erf_ + y[i]);
    }
    sumdL.block(0,0,dim,1) = sumdLdw;
    sumdL[dim] = sumdLdb;

    return sumdL;
}


/*  */
Eigen::VectorXd SvmGsuProblem::computeSumOfLossGradDiagXspace(Eigen::VectorXd w,
                                                              double b,
                                                              vector<int> ids)
{
    Eigen::VectorXd sumdL(dim+1);
    Eigen::VectorXd sumdLdw(dim);
    Eigen::VectorXd w_sq(dim);
    double sumdLdb;

    sumdL.setZero(dim+1);
    sumdLdw.setZero(dim);
    w_sq.setZero(dim);
    sumdLdb = 0.0;

    double d_mu = 0.0;
    double d_sigma = 0.0;
    double r = 0.0;
    double erf_ = 0.0;
    double exp_ = 0.0;

    int i;
    for (unsigned t=0; t<ids.size(); t++ )
    {
        i = ids[t];
        d_mu = y[i] - w.transpose()*x[i] - b;
        w_sq = w.array().square();
        d_sigma = sqrt(Sigma_xd[i].transpose() * w_sq);
        r = d_mu / d_sigma;
        erf_ = erf(0.5*sqrt(2.0)*r);
        exp_ = exp(-0.5 * r * r);
        // Derivative of loss wrt w
        sumdLdw += (exp_ / (sqrt(2.0*M_PI) * d_sigma))*(Sigma_xd[i].cwiseProduct(w)) - 0.5 * (erf_ + y[i]) * x[i];
        // Derivative of loss wrt b
        sumdLdb -= 0.5 * (erf_ + y[i]);
    }
    sumdL.block(0,0,dim,1) = sumdLdw;
    sumdL[dim] = sumdLdb;

    return sumdL;
}


/*  */
Eigen::VectorXd SvmGsuProblem::computeSumOfLossGradIsoXspace(Eigen::VectorXd w,
                                                             double b,
                                                             vector<int> ids)
{
    Eigen::VectorXd sumdL(dim+1);
    Eigen::VectorXd sumdLdw(dim);
    Eigen::VectorXd w_sq(dim);
    double sumdLdb;

    sumdL.setZero(dim+1);
    sumdLdw.setZero(dim);
    w_sq.setZero(dim);
    sumdLdb = 0.0;

    double d_mu = 0.0;
    double d_sigma = 0.0;
    double r = 0.0;
    double erf_ = 0.0;
    double exp_ = 0.0;

    int i;
    for (unsigned t=0; t<ids.size(); t++ )
    {
        i = ids[t];
        d_mu = y[i] - w.transpose()*x[i] - b;
        d_sigma = sqrt(Sigma_xi[i]) * w.norm();
        r = d_mu / d_sigma;
        erf_ = erf(0.5*sqrt(2.0)*r);
        exp_ = exp(-0.5 * r * r);
        // Derivative of loss wrt w
        sumdLdw += (exp_ / (sqrt(2.0*M_PI) * d_sigma)) * (Sigma_xi[i] * w) - 0.5 * (erf_ + y[i]) * x[i];
        // Derivative of loss wrt b
        sumdLdb -= 0.5 * (erf_ + y[i]);
    }
    sumdL.block(0,0,dim,1) = sumdLdw;
    sumdL[dim] = sumdLdb;

    return sumdL;
}


/* Solve the linear SVM-GSU in linear subspaces */
void SvmGsuProblem::solveLSVMGSUzSpace()
{
    int T, k;

    /* Get number of iterations T, sampling size k */
    T = params.getT();
    k = params.getK();

    /* Initialize linear model w = (w_1, w_2, ..., w_dim, b)^T, b */
    initW(dim, params.getLambda());
    b = 0.0;

    vector<int> ids_t(l);
    iota(ids_t.begin(), ids_t.end(), 0);

    double eta_t = 1.0;
    Eigen::VectorXd sumdL(dim+1);
    sumdL.setZero(dim+1);

    // Define aux. optimization variables
    Eigen::VectorXd sumdLdw(dim);
    sumdLdw.setZero(dim);
    double sumdLdb = 0.0;

    for (int t=1; t<=T; t++)
    {
        // Get random subset of training examples' indices (ids_t)
        auto engine = std::default_random_engine{};
        std::shuffle(std::begin(ids_t), std::end(ids_t), engine);
        vector<int> ids_k(ids_t.begin(), ids_t.begin()+k);

        // Learning rate
        eta_t = 1.0/( params.getLambda()*(double)t );

        // Compute sum of gradients of loss function
        sumdL = computeSumOfLossGradDiagZspace(w, b, ids_k);
        sumdLdw = sumdL.block(0,0,dim,1);
        sumdLdb = sumdL[dim];

        // Update optimization variables
        w = (1.0-1.0/double(t))*w - (eta_t/((double)k)) * sumdLdw;
        w = min( 1.0, 1.0/sqrt(params.getLambda()*w.squaredNorm()) ) * w;
        b = b - (eta_t/((double)k)) * sumdLdb;
    }

    /* Compute decision values and apply Platt scaling*/
    computeDecValues();
    plattScaling();

}


/*  */
Eigen::VectorXd SvmGsuProblem::computeSumOfLossGradDiagZspace(Eigen::VectorXd w,
                                                              double b,
                                                              vector<int> ids)
{
    Eigen::VectorXd sumdL(dim+1);
    Eigen::VectorXd sumdLdw(dim);
    Eigen::VectorXd w_sq(dim);
    double sumdLdb;

    sumdL.setZero(dim+1);
    sumdLdw.setZero(dim);
    w_sq.setZero(dim);
    sumdLdb = 0.0;

    double d_mu = 0.0;
    double d_sigma = 0.0;
    double r = 0.0;
    double erf_ = 0.0;
    double exp_ = 0.0;

    int i;
    for (unsigned t=0; t<ids.size(); t++ )
    {
        i = ids[t];

        vector<int> axes = princ_axes[i];
        Eigen::VectorXd w_z(axes.size());
        Eigen::VectorXd PSw(dim);
        Eigen::VectorXd Pzi(dim);
        Eigen::VectorXd v = Sigma_zd[i].cwiseProduct(w_z);
        PSw.setZero(dim);
        Pzi.setZero(dim);
        for (unsigned j=0; j<axes.size(); j++)
        {
            w_z(j) = w(axes[j]);
            PSw(axes[j]) = v(j);
            Pzi(axes[j]) = z[i][j];
        }

        d_mu = y[i] - z[i].transpose()*w_z - b;
        w_sq = w_z.array().square();
        d_sigma = sqrt(Sigma_zd[i].transpose() * w_sq);
        r = d_mu / d_sigma;
        erf_ = erf(0.5*sqrt(2.0)*r);
        exp_ = exp(-0.5 * r * r);

        // Derivative of loss wrt w
        sumdLdw += (exp_ / (sqrt(2.0*M_PI) * d_sigma))*PSw
                - 0.5 * (erf_ + y[i]) * Pzi;

        // Derivative of loss wrt b
        sumdLdb -= 0.5 * (erf_ + y[i]);
    }
    sumdL.block(0,0,dim,1) = sumdLdw;
    sumdL[dim] = sumdLdb;

    return sumdL;
}

/* Compute kernel matrix */
void SvmGsuProblem::computeKernelMatrix()
{
    // TODO: Implement!
    Eigen::MatrixXd D;
    D.setRandom(l, l);
    K.setZero(l, l);
    Eigen::MatrixXd X;
    X.setZero(l, dim);
    for (int i=0; i<l; i++)
        X.row(i) = x[i];

    K = exp(D.array());

    for (int i=0; i<l; i++){
        Eigen::VectorXd x_i = X.row(i);
        for (int j=0; j<l; j++){
            Eigen::VectorXd x_j = X.row(j);
            K(i, j) = exp(-params.getGamma() * (x_i-x_j).squaredNorm());
        }
    }

}


/* Solve the KSVM-iGSU */
void SvmGsuProblem::solveKSVMiGSU()
{
    int T, k;

    /* Get number of iterations T, sampling size k */
    T = params.getT();
    k = params.getK();

    /* Initialize KSVM-iGSU's alpha = (alpha_1, alpha_2, ..., alpha_l)^T, b */
    initAlpha(l, params.getLambda());
    b = 0.0;

    vector<int> ids_t(l);
    iota(ids_t.begin(), ids_t.end(), 0);

    double eta_t = 1.0;
    Eigen::VectorXd sumdLH(l + 1);
    sumdLH.setZero(l + 1);

    // Define aux. optimization variables
    Eigen::VectorXd sumdLHdalpha(l);
    sumdLHdalpha.setZero(l);
    double sumdLHdb = 0.0;

    /* SGD iteration */
    for (int t=1; t<=T; t++){

        // Get random subset of training examples' indices (ids_t)
        auto engine = std::default_random_engine{};
        std::shuffle(std::begin(ids_t), std::end(ids_t), engine);
        vector<int> ids_k(ids_t.begin(), ids_t.begin()+k);

        // Learning rate
        eta_t = 1.0/(params.getLambda() * (double)t);

        // Compute sum of gradients of loss function
        sumdLH = computeSumOfLossGradKSVMiGSU(alpha, b, ids_k);
        sumdLHdalpha = sumdLH.block(0, 0, l, 1);
        sumdLHdb = sumdLH[l];

        // Update optimization variables
        alpha = alpha - (1.0 / ((double)t)) * (K * alpha) - (eta_t / ((double)k)) * sumdLHdalpha;
        alpha = min( 1.0, 1.0/sqrt(params.getLambda() * alpha.squaredNorm()) ) * alpha;
        b = b + 0.5 * (eta_t / ((double)k)) * sumdLHdb;
    }

    /* Compute decision values and apply Platt scaling*/
    computeDecValues();
    plattScaling();
}

/*  */
Eigen::VectorXd SvmGsuProblem::computeSumOfLossGradKSVMiGSU(Eigen::VectorXd alpha,
                                                            double b,
                                                            vector<int> ids)
{
    Eigen::VectorXd K_i;
    Eigen::VectorXd sumdLH(l + 1);
    Eigen::VectorXd sumdLHdalpha(l);
    double sumdLHdb;
    double aKa;
    double d_mu;
    double d_sigma;
    double r;
    double erf_;
    double exp_;

    sumdLH.setZero(l + 1);
    sumdLHdalpha.setZero(l);
    sumdLHdb = 0.0;
    K_i.setZero(l);
    aKa = alpha.transpose() * K * alpha;
    d_mu = 0.0;
    d_sigma = 0.0;
    r = 0.0;
    erf_ = 0.0;
    exp_ = 0.0;

    int i;
    for (unsigned t=0; t<ids.size(); t++ )
    {
        i = ids[t];
        K_i = K.col(i);
        d_mu = 1.0 - y[i] * (K.col(i).transpose() * alpha + b);
        d_sigma = sqrt(2.0 * Sigma_xi[i] * aKa);
        r = d_mu / d_sigma;
        erf_ = erf(r);
        exp_ = exp(-r * r);
        // Derivative of loss wrt alpha
        sumdLHdalpha += (exp_ / (sqrt(M_PI) * d_sigma)) * (K * alpha) - 0.5 * (erf_ + 1.0) * K_i;
        // Derivative of loss wrt b
        sumdLHdb -= 0.5 * (erf_ + 1.0);
    }
    sumdLH.block(0, 0, l, 1) = sumdLHdalpha;
    sumdLH[l] = sumdLHdb;

    return sumdLH;
}



/*******************************************************************************
 *                                                                             *
 *                           [ Write Output Files ]                            *
 *                                                                             *
 *******************************************************************************/

/* Write model file for linear SVM-GSU */
void SvmGsuProblem::writeLinearModelFile(char* model_filename)
{
    ofstream m;
    m.open (model_filename);
    m << "SVM-GSU Model_File" << endl;
    m << "kernel_type linear" << endl;
    m << "lambda" << " " << params.getLambda() << endl;
    m << "sigmA " << sigmA << endl;
    m << "sigmB " << sigmB << endl;
    m << "cov_mat_type" << " " << params.getCovMatTypeVerb() << endl;
    m << "w";
    for (int j=0; j<dim; j++)
        m << " " << w[j];
    m << endl;
    m << "b " << b << endl;
    m.close();
}


/* Write model file for kernel SVM-iGSU (RBF) */
void SvmGsuProblem::writeKernelModelFile(char* model_filename)
{
    ofstream m;
    m.open (model_filename);
    m << "SVM-GSU Model_File" << endl;
    m << "kernel_type rbf" << endl;
    m << "lambda" << " "   << params.getLambda() << endl;
    m << "gamma" << " "   << params.getGamma() << endl;
    m << "b" << " " << b << endl;
    m << "sigmA " << sigmA << endl;
    m << "sigmB " << sigmB << endl;
    m << "cov_mat_type" << " " << params.getCovMatTypeVerb() << endl;
    m << "SV" << endl;

    for (int i=0; i<l; i++)
    {
        m << alpha[i];
        for (int j=0; j<dim; j++)
            m << " " << (j+1) << ":" << x[i][j];
        m << endl;
    }
    m.close();
}



/*******************************************************************************
 *                                                                             *
 *                          [ Auxiliary Functions ]                            *
 *                                                                             *
 *******************************************************************************/
/*  */
double SvmGsuProblem::getLabel(double t){return (t>0.0 ? +1.0 : -1.0);}


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

    while( i<n1 && j<n2 )
    {
        if ( A[i].compare(B[j]) > 0 )
            j++;
        else if ( B[j].compare(A[i]) > 0 )
            i++;
        else
        {
            intersection.push_back( A[i] );
            i++;
            j++;
        }
    }
    return intersection;
}


/* Initialize LSVM-GSU's w vector */
void SvmGsuProblem::initW(int d, double lambda)
{
    w.resize(d);
    w.setZero(d);

    double lower_bound = -1.0/sqrt(lambda * d);
    double upper_bound = +1.0/sqrt(lambda * d);
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    for (int j=0; j<d; j++)
        w[j] = unif(re);
}


/* Initialize KSVM-iGSU's alpha vector */
void SvmGsuProblem::initAlpha(int d, double lambda)
{
    alpha.resize(d);
    alpha.setZero(d);

    double lower_bound = -1.0/sqrt(lambda * d);
    double upper_bound = +1.0/sqrt(lambda * d);
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    for (int j=0; j<d; j++)
        alpha[j] = unif(re);
}


/*  */
void SvmGsuProblem::computeDecValues()
{
    deci.resize(l);
    deci.setZero(l);
    switch (params.getKernelType())
    {
        /* Linear kernel: Optimal variables: w, b */
        case 0:
            //deci.array() = (X*w).array() + b;
            for (int i=0; i<l; i++)
                deci(i) = w.transpose() * x[i] + b;
            break;
        /* RBF kernel: Optimal variables: a, b */
        case 2:
            // deci.array() = (K*alpha).array() + b;
            break;
        default:
            break;
    }

}


/*  */
void SvmGsuProblem::plattScaling()
{
    sigmA = -1.0;
    sigmB =  0.0;

    int i, it;

    /* Parameters Setting */
    double maxiter = 300;    // Maximum number of iterations
    double minstep = 1E-10;  // Minimum step taken in line search
    double sigma   = 1E-12;  // H' = H + sigma*I (Set to any value > 0)


    /* Construct initial values : target support in array t
                                : initial function value in fval */
    double prior0 = 0.0;        // Number of negative training examples
    double prior1 = 0.0;        // Number of positive training examples
    double len    = l;          // Total number of training examples (prior0 + prior1)

    /* Compute prior0, prior1 */
    for (i=0;i<len;i++)
        if (y[i]<0) prior0++;
        else        prior1++;

    double hiTarget = ( prior1 + 1.0 )/( prior1 + 2.0 );
    double loTarget = 1.0/( prior0 + 2.0 );


    vector<double> t(l, 0.0);

    for (i=0;i<len;i++)
        if ( y[i]>0 ) t[i] = hiTarget;
        else          t[i] = loTarget;

    double A    = 0.0;
    double B    = log((prior0+1.0)/(prior1+1.0));
    double fval = 0.0;
    double fApB = 0.0;

    double h11  = sigma;
    double h22  = sigma;
    double h21  = 0.0;
    double g1   = 0.0;
    double g2   = 0.0;
    double p    = 0.0;
    double q    = 0.0;
    double d1   = 0.0;
    double d2   = 0.0;
    double det  = 0.0;
    double dA   = 0.0;
    double dB   = 0.0;
    double gd   = 0.0;
    double newA = 0.0;
    double newB = 0.0;
    double newf = 0.0;

    double stepsize = 1.0;

    for (i=0;i<len;i++)
    {
        fApB = deci(i)*A + B;
        if (fApB>=0) fval +=     t[i]*fApB + log1p( exp(-fApB) );
        else         fval += (1-t[i])*fApB + log1p( exp(fApB)  );
    }



    for (it=0;it<maxiter;it++)
    {
        // Update Gradient and Hessian matrix (use H' = H + sigma*I)
        h11 = sigma;
        h22 = sigma;
        h21 = 0.0;
        g1  = 0.0;
        g2  = 0.0;
        for (i=0;i<len;i++)
        {
            fApB = deci(i)*A + B;
            if (fApB>=0)
            {
                p = exp(-fApB)/(1.0+exp(-fApB));
                q = 1.0/(1.0+exp(-fApB));
            }
            else
            {
                p = 1.0/(1.0+exp(fApB));
                q = exp(fApB)/(1.0+exp(fApB));
            }

            d2  = p*q;
            h11 = h11 + deci(i)*deci(i)*d2;
            h22 = h22 + d2;
            h21 = deci(i)*d2;
            d1  = t[i] - p;
            g1  = g1 + deci(i)*d1;
            g2  = g2 + d1;
        }

        /* Stoping criteria */
        if ( ( abs(g1)<1E-5 ) && ( abs(g2)<1E-5 ) )
            break;

        /* Compute modified Newton directions */
        det = h11*h22 - h21*h21;
        dA  = -( +h22*g1 - h21*g2 )/det;
        dB  = -( -h21*g1 + h11*g2 )/det;
        gd  = g1*dA + g2*dB;

        stepsize = 1;
        while (stepsize >= minstep) // Line Search
        {
            newA = A + stepsize*dA;
            newB = B + stepsize*dB;
            newf = 0.0;

            for (i=0;i<len;i++)
            {
                fApB = deci(i)*newA + newB;
                if (fApB>=0)  newf += t[i]*fApB + log1p( exp(-fApB) ) ;
                else          newf += (t[i]-1)*fApB + log1p( exp(fApB) ) ;
            }

            if ( newf<fval + 0.0001*stepsize*gd )
            {
                A    = newA;
                B    = newB;
                fval = newf;
                break; // Sufficient decrease satisfied
            }
            else
                stepsize = stepsize/2.0;

            if ( stepsize < minstep )
            {
                //cout << "[Platt Scaling Warning!] Line search failed!" << endl;
                break;
            }
        }
    }
    //if ( it>=maxiter ) cout << "[Platt Scaling Warning!] Reaching maximum iterations." << endl;
    sigmA = A;
    sigmB = B;

}
