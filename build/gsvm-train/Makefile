#############################################################################
# Makefile for building: gsvm-train
#############################################################################

DEBUG		  = 0
CXX           = g++
ifeq ($(DEBUG), 1)
CXXFLAGS      = -pipe -g -std=gnu++11 -Wall -W -fPIC
else
CXXFLAGS      = -pipe -O2 -march=x86-64 -mtune=generic -O2 -pipe -fstack-protector-strong -fno-plt -std=gnu++11 -Wall -W -fPIC
endif
INCPATH       = -I../../source/gsvm-train -I.
LINK          = g++
LFLAGS        = 
LIBS          = $(SUBLIBS)  

####### Files

SOURCES       = gsvm-train.cpp \
		svmgsuparams.cpp \
		svmgsuproblem.cpp \
		meanvector.cpp \
		label.cpp \
		isocovariancematrix.cpp \
		diagcovariancematrix.cpp \
		fullcovariancematrix.cpp 
OBJECTS       = gsvm-train.o \
		svmgsuparams.o \
		svmgsuproblem.o \
		meanvector.o \
		label.o \
		isocovariancematrix.o \
		diagcovariancematrix.o \
		fullcovariancematrix.o
TARGET        = gsvm-train

####### Build rules

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

all: $(TARGET)

clean:
	rm -rf *.o gsvm-train

####### Compile

gsvm-train.o: ../../source/gsvm-train/gsvm-train.cpp ../../source/gsvm-train/svmgsuparams.h \
		../../source/gsvm-train/svmgsuproblem.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o gsvm-train.o ../../source/gsvm-train/gsvm-train.cpp

svmgsuparams.o: ../../source/gsvm-train/svmgsuparams.cpp ../../source/gsvm-train/svmgsuparams.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o svmgsuparams.o ../../source/gsvm-train/svmgsuparams.cpp

svmgsuproblem.o: ../../source/gsvm-train/svmgsuproblem.cpp ../../source/gsvm-train/svmgsuproblem.h \
		../../source/gsvm-train/svmgsuparams.h \
		../../source/gsvm-train/meanvector.h \
		../../source/gsvm-train/label.h \
		../../source/gsvm-train/fullcovariancematrix.h \
		../../source/gsvm-train/diagcovariancematrix.h \
		../../source/gsvm-train/isocovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o svmgsuproblem.o ../../source/gsvm-train/svmgsuproblem.cpp

meanvector.o: ../../source/gsvm-train/meanvector.cpp ../../source/gsvm-train/meanvector.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o meanvector.o ../../source/gsvm-train/meanvector.cpp

label.o: ../../source/gsvm-train/label.cpp ../../source/gsvm-train/label.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o label.o ../../source/gsvm-train/label.cpp

isocovariancematrix.o: ../../source/gsvm-train/isocovariancematrix.cpp ../../source/gsvm-train/isocovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o isocovariancematrix.o ../../source/gsvm-train/isocovariancematrix.cpp

diagcovariancematrix.o: ../../source/gsvm-train/diagcovariancematrix.cpp ../../source/gsvm-train/diagcovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o diagcovariancematrix.o ../../source/gsvm-train/diagcovariancematrix.cpp

fullcovariancematrix.o: ../../source/gsvm-train/fullcovariancematrix.cpp ../../source/gsvm-train/fullcovariancematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fullcovariancematrix.o ../../source/gsvm-train/fullcovariancematrix.cpp
