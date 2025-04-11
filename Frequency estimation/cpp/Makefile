CXXFLAGS	= -march=native -O3
CXXFLAGS	+= -std=c++11
CXXFLAGS	+= -isystem ./include
CXXFLAGS	+= -MMD
CC	= g++

TARGET	=  EstimateDistribution_Lap EstimateDistribution_RR EstimateDistribution_RAP EstimateDistribution_OLH

all: $(TARGET)

EstimateDistribution_Lap: EstimateDistribution_Lap.o
EstimateDistribution_RR: EstimateDistribution_RR.o
EstimateDistribution_RAP: EstimateDistribution_RAP.o
EstimateDistribution_OLH: EstimateDistribution_OLH.o

clean:
	rm -f *~ $(TARGET) *.d *.o

-include *.d
