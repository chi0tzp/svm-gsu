#include "evaluationmetrics.h"
#include <fstream>
#include <iostream>

EvaluationMetrics::EvaluationMetrics(vector<OutputScore> output)
    : num_ret(0),
      num_rel(0),
      TP(0),
      TN(0),
      FP(0),
      FN(0),
      acc(0.0),
      ap(0.0),
      ap5(0.0),
      ap10(0.0),
      ap15(0.0),
      ap20(0.0),
      ap30(0.0),
      ap100(0.0),
      ap200(0.0),
      ap500(0.0),
      ap1000(0.0),
      prec(0.0),
      prec5(0.0),
      prec10(0.0),
      prec15(0.0),
      prec20(0.0),
      prec30(0.0),
      prec100(0.0),
      prec200(0.0),
      prec500(0.0),
      prec1000(0.0),
      recall(0.0),
      recall5(0.0),
      recall10(0.0),
      recall15(0.0),
      recall20(0.0),
      recall30(0.0),
      recall100(0.0),
      recall200(0.0),
      recall500(0.0),
      recall1000(0.0),
      fscore(0.0)
{
    num_ret = (int)output.size();
    countNumRel(output);

    computeTP(output);
    computeTN(output);
    computeFP(output);
    computeFN(output);

    acc    = computeAccuracy(output);
    prec   = computePrecision(output);
    ap     = computeAveragePrecision(output);
    recall = computeRecall(output);
    fscore = 2.0*(prec*recall)/(prec+recall);

    if ( num_ret>=5 )
    {
        vector<OutputScore> output5(output.begin(), output.begin()+5);
        prec5   = computePrecision(output5);
        ap5     = computeAveragePrecision(output5);
        recall5 = computeRecall(output5);
    }
    if ( num_ret>=10 )
    {
        vector<OutputScore> output10(output.begin(), output.begin()+10);
        prec10   = computePrecision(output10);
        ap10     = computeAveragePrecision(output10);
        recall10 = computeRecall(output10);
    }
    if ( num_ret>=15 )
    {
        vector<OutputScore> output15(output.begin(), output.begin()+15);
        prec15   = computePrecision(output15);
        ap15     = computeAveragePrecision(output15);
        recall15 = computeRecall(output15);
    }
    if ( num_ret>=20 )
    {
        vector<OutputScore> output20(output.begin(), output.begin()+20);
        prec20   = computePrecision(output20);
        ap20     = computeAveragePrecision(output20);
        recall20 = computeRecall(output20);
    }
    if ( num_ret>=30 )
    {
        vector<OutputScore> output30(output.begin(), output.begin()+30);
        prec30   = computePrecision(output30);
        ap30     = computeAveragePrecision(output30);
        recall30 = computeRecall(output30);
    }
    if ( num_ret>=100 )
    {
        vector<OutputScore> output100(output.begin(), output.begin()+100);
        prec100   = computePrecision(output100);
        ap100     = computeAveragePrecision(output100);
        recall100 = computeRecall(output100);
    }
    if ( num_ret>=200 )
    {
        vector<OutputScore> output200(output.begin(), output.begin()+200);
        prec200   = computePrecision(output200);
        ap200     = computeAveragePrecision(output200);
        recall200 = computeRecall(output200);
    }
    if ( num_ret>=500 )
    {
        vector<OutputScore> output500(output.begin(), output.begin()+500);
        prec500   = computePrecision(output500);
        ap500     = computeAveragePrecision(output500);
        recall500 = computeRecall(output500);
    }
    if ( num_ret>=1000 )
    {
        vector<OutputScore> output1000(output.begin(), output.begin()+1000);
        prec1000   = computePrecision(output1000);
        ap1000     = computeAveragePrecision(output1000);
        recall1000 = computeRecall(output1000);
    }

}


void EvaluationMetrics::countNumRel(vector<OutputScore> output)
{
    for (unsigned i=0; i<output.size(); i++)
        if (output[i].getPredictedLabel()==+1)
            num_rel++;
}


void EvaluationMetrics::computeTP(vector<OutputScore> output)
{
    for (unsigned i=0; i<output.size(); i++)
        if ( (output[i].getPredictedLabel()==+1) && (output[i].getActualLabel()==+1) )
            TP++;
}


void EvaluationMetrics::computeTN(vector<OutputScore> output)
{
    for (unsigned i=0; i<output.size(); i++)
        if ( (output[i].getPredictedLabel()==-1) && (output[i].getActualLabel()==-1) )
            TN++;
}


void EvaluationMetrics::computeFP(vector<OutputScore> output)
{
    for (unsigned i=0; i<output.size(); i++)
        if ( output[i].getPredictedLabel() == +1 )
            if (output[i].getActualLabel() == -1)
                FP++;
}


void EvaluationMetrics::computeFN(vector<OutputScore> output)
{
    for (unsigned i=0; i<output.size(); i++)
        if ( output[i].getPredictedLabel() == -1 )
            if (output[i].getActualLabel() == +1)
                FN++;
}




double EvaluationMetrics::computeAccuracy(vector<OutputScore> output)
{
    double accuracy = 0.0;
    double tp  = 0.0;
    double tn  = 0.0;

    for (unsigned i=0; i<output.size(); i++)
    {
        if ( (output[i].getPredictedLabel()==+1) && (output[i].getActualLabel()==+1) )
            tp++;
        if ( (output[i].getPredictedLabel()==-1) && (output[i].getActualLabel()==-1) )
            tn++;
    }
    accuracy = (tp+tn)/(double)output.size();
    return accuracy;
}


double EvaluationMetrics::computePrecision(vector<OutputScore> output)
{
    double precision = 0.0;
    double        tp = 0.0;
    double        fp = 0.0;

    for (unsigned i=0; i<output.size(); i++)
    {
        if ( output[i].getPredictedLabel() == +1 )
        {
            if (output[i].getActualLabel() == +1)
                tp++;
            else
                fp++;
        }
    }
    precision = tp/(tp+fp);
    return precision;
}


double EvaluationMetrics::computeRecall(vector<OutputScore> output)
{
    double recall = 0.0;
    double TP = 0.0;
    double  P = 0.0;

    for (unsigned i=0; i<output.size(); i++)
    {
        if (output[i].getPredictedLabel()==+1)
        {
            P++;
            if (output[i].getActualLabel()==+1)
                TP++;
        }
    }
    recall = TP/P;
    return recall;
}



double EvaluationMetrics::computeAveragePrecision(vector<OutputScore> output)
{
    double ap = 0.0;
    double tp = 0.0;

    for (unsigned i=0; i<output.size(); i++)
    {
        if (output[i].getActualLabel()==+1)
        {
            tp++;
            ap += tp/double(i+1);
        }
    }
    ap /= tp;
    return ap;
}





void EvaluationMetrics::writeEvaluationMetrics(char* evaluation_metrics_filename)
{
    ofstream f;
    f.open (evaluation_metrics_filename);
    f << "num_ret\t\t"     <<   num_ret    << endl;
    f << "num_rel\t\t"     <<   num_rel    << endl;
    f << "TP\t\t"          <<        TP    << endl;
    f << "TN\t\t"          <<        TN    << endl;
    f << "FP\t\t"          <<        FP    << endl;
    f << "FN\t\t"          <<        FN    << endl;
    f << "acc\t\t"         <<       acc    << endl;
    f << "prec\t\t"        <<      prec    << endl;
    f << "prec@5\t\t"      <<     prec5    << endl;
    f << "prec@10\t\t"     <<    prec10    << endl;
    f << "prec@15\t\t"     <<    prec15    << endl;
    f << "prec@20\t\t"     <<    prec20    << endl;
    f << "prec@30\t\t"     <<    prec30    << endl;
    f << "prec@100\t"      <<   prec100    << endl;
    f << "prec@200\t"      <<   prec200    << endl;
    f << "prec@500\t"      <<   prec500    << endl;
    f << "prec@1000\t"     <<  prec1000    << endl;
    f << "ap\t\t"          <<      ap      << endl;
    f << "ap@5\t\t"        <<     ap5      << endl;
    f << "ap@10\t\t"       <<    ap10      << endl;
    f << "ap@15\t\t"       <<    ap15      << endl;
    f << "ap@20\t\t"       <<    ap20      << endl;
    f << "ap@30\t\t"       <<    ap30      << endl;
    f << "ap@100\t\t"      <<   ap100      << endl;
    f << "ap@200\t\t"      <<   ap200      << endl;
    f << "ap@500\t\t"      <<   ap500      << endl;
    f << "ap@1000\t\t"     <<  ap1000      << endl;
    f << "recall\t\t"      <<      recall  << endl;
    f << "recall@5\t"      <<     recall5  << endl;
    f << "recall@10\t"     <<    recall10  << endl;
    f << "recall@15\t"     <<    recall15  << endl;
    f << "recall@20\t"     <<    recall20  << endl;
    f << "recall@30\t"     <<    recall30  << endl;
    f << "recall@100\t"    <<   recall100  << endl;
    f << "recall@200\t"    <<   recall200  << endl;
    f << "recall@500\t"    <<   recall500  << endl;
    f << "recall@1000\t"   <<  recall1000  << endl;
    f << "f-score\t\t"     <<      fscore  << endl;
    f.close();
}
