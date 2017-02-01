#include "outputscore.h"
#include <vector>

#ifndef EVALUATIONMETRICS_H
#define EVALUATIONMETRICS_H



class EvaluationMetrics
{
    public:
        /* Ctor */
        EvaluationMetrics(vector<OutputScore>);
        /* Member functions for computing evaluation metrics */
        void countNumRel(vector<OutputScore>);
        void computeTP(vector<OutputScore>);
        void computeTN(vector<OutputScore>);
        void computeFP(vector<OutputScore>);
        void computeFN(vector<OutputScore>);
        double computeAccuracy(vector<OutputScore>);
        double computePrecision(vector<OutputScore>);
        double computeAveragePrecision(vector<OutputScore>);
        double computeRecall(vector<OutputScore>);
        void writeEvaluationMetrics(char*);

    private:

        int    num_ret;     // Total number of documents retrieved over all queries
        int    num_rel;     // Total number of relevant documents over all queries
        int    TP;          // Total number of true positives  over all queries
        int    TN;          // Total number of true negatives  over all queries
        int    FP;          // Total number of false positives  over all queries
        int    FN;          // Total number of false negatives  over all queries
        double acc;         // Accuracy over retrieved set
        double ap;          // Average Precision (AP) over retrieved set
        double ap5;         // Average Precision (AP) after    5 docs retrieved
        double ap10;        // Average Precision (AP) after   10 docs retrieved
        double ap15;        // Average Precision (AP) after   15 docs retrieved
        double ap20;        // Average Precision (AP) after   20 docs retrieved
        double ap30;        // Average Precision (AP) after   30 docs retrieved
        double ap100;       // Average Precision (AP) after  100 docs retrieved
        double ap200;       // Average Precision (AP) after  200 docs retrieved
        double ap500;       // Average Precision (AP) after  500 docs retrieved
        double ap1000;      // Average Precision (AP) after 1000 docs retrieved
        double prec;        // Precision over retrieved set
        double prec5;       // Precision after    5 docs retrieved
        double prec10;      // Precision after   10 docs retrieved
        double prec15;      // Precision after   15 docs retrieved
        double prec20;      // Precision after   20 docs retrieved
        double prec30;      // Precision after   30 docs retrieved
        double prec100;     // Precision after  100 docs retrieved
        double prec200;     // Precision after  200 docs retrieved
        double prec500;     // Precision after  500 docs retrieved
        double prec1000;    // Precision after 1000 docs retrieved
        double recall;      // Recall over retrieved set
        double recall5;     // Recall after    5 docs retrieved
        double recall10;    // Recall after   10 docs retrieved
        double recall15;    // Recall after   15 docs retrieved
        double recall20;    // Recall after   20 docs retrieved
        double recall30;    // Recall after   30 docs retrieved
        double recall100;   // Recall after  100 docs retrieved
        double recall200;   // Recall after  200 docs retrieved
        double recall500;   // Recall after  500 docs retrieved
        double recall1000;  // Recall after 1000 docs retrieved
        double fscore;      // F_1-score: the harmonic mean of
                            // precision and recall, i.e.,
                            //      F_1 = 2.0 * (prec*recall) / (prec+recall)
};


#endif // EVALUATIONMETRICS_H
