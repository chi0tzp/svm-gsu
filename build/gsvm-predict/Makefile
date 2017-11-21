#############################################################################
# Makefile for building: gsvm-predict
#############################################################################

DEBUG		  = 0
CXX           = g++
ifeq ($(DEBUG), 1)
CXXFLAGS      = -pipe -g -std=gnu++11 -Wall -W -fPIC
else
CXXFLAGS      = -pipe -O2 -march=x86-64 -mtune=generic -O2 -pipe -fstack-protector-strong -fno-plt -std=gnu++11 -Wall -W -fPIC $(DEFINES)
endif
INCPATH       = -I../../source/gsvm-predict -I.
LINK          = g++
LFLAGS        = -Wl,-O1 -Wl,-O1,--sort-common,--as-needed,-z,relro,-z,now

####### Files

SOURCES       = ../../source/gsvm-predict/gsvm-predict.cpp \
		../../source/gsvm-predict/svmgsuproblem.cpp \
		../../source/gsvm-predict/svmgsuparams.cpp \
		../../source/gsvm-predict/svmgsumodel.cpp \
		../../source/gsvm-predict/outputscore.cpp \
		../../source/gsvm-predict/meanvector.cpp \
		../../source/gsvm-predict/label.cpp \
		../../source/gsvm-predict/isocovariancematrix.cpp \
		../../source/gsvm-predict/fullcovariancematrix.cpp \
		../../source/gsvm-predict/evaluationmetrics.cpp \
		../../source/gsvm-predict/diagcovariancematrix.cpp 
OBJECTS       = gsvm-predict.o \
		svmgsuproblem.o \
		svmgsuparams.o \
		svmgsumodel.o \
		outputscore.o \
		meanvector.o \
		label.o \
		isocovariancematrix.o \
		fullcovariancematrix.o \
		evaluationmetrics.o \
		diagcovariancematrix.o

TARGET        = gsvm-predict

####### Build rules

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS)

all: Makefile $(TARGET)

clean:
	rm -rf *.o gsvm-predict

####### Compile

gsvm-predict.o: ../../source/gsvm-predict/gsvm-predict.cpp ../../source/gsvm-predict/svmgsumodel.h \
		../../source/gsvm-predict/svmgsuparams.h \
		../../source/gsvm-predict/svmgsuproblem.h \
		../../source/gsvm-predict/outputscore.h \
		../../source/gsvm-predict/evaluationmetrics.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o gsvm-predict.o ../../source/gsvm-predict/gsvm-predict.cpp

svmgsuproblem.o: ../../source/gsvm-predict/svmgsuproblem.cpp ../../source/gsvm-predict/svmgsuproblem.h \
		../../source/gsvm-predict/svmgsuparams.h \
		../../source/gsvm-predict/svmgsumodel.h \
		../../source/gsvm-predict/outputscore.h \
		../../source/gsvm-predict/meanvector.h \
		../../source/gsvm-predict/fullcovariancematrix.h \
		../../source/gsvm-predict/diagcovariancematrix.h \
		../../source/gsvm-predict/isocovariancematrix.h \
		../../source/gsvm-predict/label.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o svmgsuproblem.o ../../source/gsvm-predict/svmgsuproblem.cpp

svmgsuparams.o: ../../source/gsvm-predict/svmgsuparams.cpp ../../source/gsvm-predict/svmgsuparams.h \
		../../source/gsvm-predict/svmgsumodel.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o svmgsuparams.o ../../source/gsvm-predict/svmgsuparams.cpp

svmgsumodel.o: ../../source/gsvm-predict/svmgsumodel.cpp ../../source/gsvm-predict/svmgsumodel.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o svmgsumodel.o ../../source/gsvm-predict/svmgsumodel.cpp

outputscore.o: ../../source/gsvm-predict/outputscore.cpp ../../source/gsvm-predict/outputscore.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o outputscore.o ../../source/gsvm-predict/outputscore.cpp

meanvector.o: ../../source/gsvm-predict/meanvector.cpp ../../source/gsvm-predict/meanvector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o meanvector.o ../../source/gsvm-predict/meanvector.cpp

label.o: ../../source/gsvm-predict/label.cpp ../../source/gsvm-predict/label.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o label.o ../../source/gsvm-predict/label.cpp

isocovariancematrix.o: ../../source/gsvm-predict/isocovariancematrix.cpp ../../source/gsvm-predict/isocovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o isocovariancematrix.o ../../source/gsvm-predict/isocovariancematrix.cpp

fullcovariancematrix.o: ../../source/gsvm-predict/fullcovariancematrix.cpp ../../source/gsvm-predict/fullcovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fullcovariancematrix.o ../../source/gsvm-predict/fullcovariancematrix.cpp

evaluationmetrics.o: ../../source/gsvm-predict/evaluationmetrics.cpp ../../source/gsvm-predict/evaluationmetrics.h \
		../../source/gsvm-predict/outputscore.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o evaluationmetrics.o ../../source/gsvm-predict/evaluationmetrics.cpp

diagcovariancematrix.o: ../../source/gsvm-predict/diagcovariancematrix.cpp ../../source/gsvm-predict/diagcovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o diagcovariancematrix.o ../../source/gsvm-predict/diagcovariancematrix.cpp
