TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    gsvm-predict.cpp \
    svmgsumodel.cpp \
    svmgsuparams.cpp \
    svmgsuproblem.cpp \
    fullcovariancematrix.cpp \
    diagcovariancematrix.cpp \
    isocovariancematrix.cpp \
    label.cpp \
    outputscore.cpp \
    meanvector.cpp \
    evaluationmetrics.cpp

HEADERS += \
    svmgsumodel.h \
    svmgsuparams.h \
    svmgsuproblem.h \
    fullcovariancematrix.h \
    diagcovariancematrix.h \
    isocovariancematrix.h \
    label.h \
    outputscore.h \
    meanvector.h \
    evaluationmetrics.h
