CXX = g++
LDLIBS =  -larmadillo -lgsl

OBJS = Polynomial.o CubicSpline.o filters.o read.o baselineAdjustment.o peaks.o output.o
HEADERS = Polynomial.h CubicSpline.h prototypes.h structs.h legendreConstants.h


nmrAnalyzer :	main.o $(OBJS)
	$(CXX) main.o $(OBJS) $(LDLIBS) -o nmrAnalyzer

main.o : main.cpp $(HEADERS)
	$(CXX) main.cpp -c

Polynomial.o : Polynomial.cpp $(HEADERS)
	$(CXX) Polynomial.cpp -c

CubicSpline.o : CubicSpline.cpp $(HEADERS)
	$(CXX) CubicSpline.cpp -c

filters.o : filters.cpp $(HEADERS)
	$(CXX) filters.cpp -c

read.o : read.cpp $(HEADERS)
	$(CXX) read.cpp -c

baselineAdjustment.o : baselineAdjustment.cpp $(HEADERS)
	$(CXX) baselineAdjustment.cpp -c

peaks.o : peaks.cpp $(HEADERS)
	$(CXX) peaks.cpp -c

output.o : output.cpp $(HEADERS)
	$(CXX) output.cpp -c

clean:
	rm *.o

pristine:
	rm *.o
	touch *
