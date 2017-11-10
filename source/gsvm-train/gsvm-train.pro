TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    gsvm-train.cpp \
    svmgsuparams.cpp \
    svmgsuproblem.cpp \
    meanvector.cpp \
    label.cpp \
    isocovariancematrix.cpp \
    diagcovariancematrix.cpp \
    fullcovariancematrix.cpp

HEADERS += \
    svmgsuparams.h \
    svmgsuproblem.h \
    meanvector.h \
    label.h \
    isocovariancematrix.h

DISTFILES +=
