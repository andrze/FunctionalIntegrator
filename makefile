ODIR=obj
CXX=g++
CXXFLAGS= -std=c++14 -Wall -pedantic -O3
CFLAGS=-I. -std=c++14 -Wall -pedantic -O3 -static 

DEPS = system.h \
	   stepfunction.h \
	   regulator.h \
       realvector.h \
       plotset.h \
       plot.h \
	   integrator.h \
	   gaussquadrature.h \
	   rungekutta.h

_OBJ = main.o \
	   system.o \
	   stepfunction.o \
	   regulator.o \
	   realvector.o \
	   plotset.o \
	   plot.o \
	   integrator.o \
	   gaussquadrature.o \
	   rungekutta.o
	   
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#hello:
#	@echo $(OBJ)
#	@echo $(DEPS)

#$(ODIR)/%.o: %.c $(DEPS)
#	@echo $(DEPS)
#	$(CXX) -c -o $@ $< $(CFLAGS)


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

make: $(OBJ)
	$(CXX) -o "FunctionalIntegrator" $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core 
  