ODIR=obj
CXX=g++
CXXFLAGS= -std=c++17 -Wall -pedantic 
CFLAGS=-I. -std=c++17 -Wall -pedantic

all: executable

debug: CXXFLAGS += -g3 -fno-omit-frame-pointer -pthread
debug: CFLAGS += -g3 -fno-omit-frame-pointer -pthread
debug: executable

release: CXXFLAGS += -O3 -pthread
release: CFLAGS += -O3  -pthread
release: executable

profile: CXXFLAGS += -O3 -g3 -pthread
profile: CFLAGS += -O3 -g3 -pthread
profile: executable

static: CXXFLAGS += -static  -O3
static: CFLAGS += -static  -O3
static: executable

DEPS = system.h \
	   stepfunction.h \
	   regulator.h \
       realvector.h \
	   integrator.h \
	   gaussquadrature.h \
	   rungekutta.h \
	   terminalplot.h \
	   numericalmethods.h

_OBJ = main.o \
	   system.o \
	   stepfunction.o \
	   regulator.o \
	   realvector.o \
	   integrator.o \
	   gaussquadrature.o \
	   rungekutta.o \
	   terminalplot.o \
	   numericalmethods.o
	   
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#hello:
#	@echo $(OBJ)
#	@echo $(DEPS)

#$(ODIR)/%.o: %.c $(DEPS)
#	@echo $(DEPS)
#	$(CXX) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

executable: $(OBJ)
	$(CXX) -o "FunctionalIntegrator" $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core 
  
