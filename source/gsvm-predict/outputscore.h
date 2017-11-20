#ifndef OUTPUTSCORE_H
#define OUTPUTSCORE_H

#include <string>
#include <vector>

using namespace std;

class OutputScore
{
    public:
        string doc_id;
        double score;
        int predicted_label;
        int actual_label;

        /* Ctor */
        OutputScore(string, double, int, int);

        /* Setters */
        void setDocId(string);
        void setScore(double);
        void setPredictedLabel(int);
        void setActualLabel(int);

        /* Getters */
        string getDocId();
        double getScore();
        int getPredictedLabel();
        int getActualLabel();

        /* Comparison operators */
        bool operator>(const OutputScore&);

};

#endif // OUTPUTSCORE_H
