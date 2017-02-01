#ifndef LABEL_H
#define LABEL_H

#include <string>

using namespace std;

class Label
{
    public:
        /* Ctor */
        Label();
        /* Setters */
        void setDocId(string);
        void setY(double);
        /* Getters */
        string getDocId();
        double getY();
    private:
        string doc_id;
        double y;
};

#endif // LABEL_H
