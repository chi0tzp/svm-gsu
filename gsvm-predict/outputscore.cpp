#include "outputscore.h"
#include <fstream>

/* Ctor */
OutputScore::OutputScore(string id, double s, int y, int actual_y)
    : doc_id(id),
      score(s),
      predicted_label(y),
      actual_label(actual_y)
{}

/* Setters */
void OutputScore::setDocId(string id){ doc_id = id; }
void OutputScore::setScore(double s){ score = s; }
void OutputScore::setPredictedLabel(int y){ predicted_label = y; }
void OutputScore::setActualLabel(int y){ actual_label = y; }

/* Getters */
string OutputScore::getDocId(){ return doc_id; }
double OutputScore::getScore(){ return score; }
int OutputScore::getPredictedLabel(){ return predicted_label; }
int OutputScore::getActualLabel(){ return actual_label; }

/*  */
bool OutputScore::operator>(const OutputScore& rhs){ return score > rhs.score; }
