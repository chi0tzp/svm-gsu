#include "label.h"

/* Default ctor */
Label::Label()
    : doc_id("###"),
      y(0.0)
{}

/* Setters */
void Label::setDocId(string id){doc_id = id;}
void Label::setY(double lbl){y = lbl;}

/* Getters */
string Label::getDocId(){return doc_id;}
double Label::getY(){return y;}
