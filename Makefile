CPPSOURCES=ols.cpp random.cpp mlp.cpp stats.cpp owo.cpp \
	timer.cpp setup.cpp main.cpp derivs.cpp olf.cpp\
        newton.cpp cgrad.cpp train.cpp cgutils.cpp validate.cpp matrix.cpp matops.cpp
vpath %.cpp src
vpath %.c src
INCLUDE := include
CPPSOURCES2=ols.cpp random.cpp mlp.cpp stats.cpp owo.cpp \
	timer.cpp setup.cpp main2.cpp derivs.cpp olf.cpp\
        newton.cpp cgrad.cpp train.cpp cgutils.cpp validate.cpp
PROC_SOURCES=ols.cpp random.cpp mlp.cpp stats.cpp owo.cpp \
        timer.cpp setup.cpp matops.cpp matrix.cpp derivs.cpp olf.cpp\
        newton.cpp cgrad.cpp train.cpp process.cpp
LM_SOURCES=mlp_lm.cpp allocmem.cpp conjugat.cpp get.cpp random.cpp \
           lm_conj.cpp setup.cpp lmfuncs.cpp matrix.cpp mlp.cpp matops.cpp \
           stats.cpp
LM_CSOURCES= mt19937ar.c
CSOURCES=allocmem.c mt19937ar.c readutils.c #utils.c
FSOURCES=
F77 = gfortran
CC=gcc
CXX=g++
CFLAGS=-g -c -Wall 
LD=ld
LDFLAGS= -L/home/melrobin/research/MLP -L/usr/local/atlas/lib -lf77blas -lcblas -ltatlas -lm -lgfortran -lrt -lQtGui -lQtCore 
MKLFLAGS= -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lrt
CPPOBJECTS=$(CPPSOURCES:.cpp=.o)
CPPOBJECTS2=$(CPPSOURCES2:.cpp=.o)
PROC_OBJS=$(PROC_SOURCES:.cpp=.o)
LM_OBJS=$(LM_SOURCES:.cpp=.o)
LM_COBJS=$(LM_CSOURCES:.c=.o)
UTESTOBJS = unit_tester.cpp matops.cpp matrix.cpp 
COBJS=$(CSOURCES:.c=.o)
FOBJS=$(FSOURCES:.f=.o)
EXECUTABLE=mlp

mlp: $(CPPOBJECTS) $(COBJS) $(FOBJS) 
	$(CXX) -o $@ $(CPPOBJECTS) $(COBJS) $(FOBJS) $(LDFLAGS)

mlp2: $(CPPOBJECTS2) $(COBJS) $(FOBJS) 
	$(CXX) -o $@ $(CPPOBJECTS2) $(COBJS) $(FOBJS) $(LDFLAGS)
all: mlp process

.cpp.o:
	$(CXX) -I$(INCLUDE) $(CFLAGS) $< -o $@

.c.o:
	$(CC) -I$(INCLUDE) $(CFLAGS) $< -o $@

.f.o:
	$(F77) -c $< -o $@

process: $(PROC_OBJS) $(COBJS)
	$(CXX) -o $@ $(PROC_OBJS) $(COBJS) $(LDFLAGS)

lm: $(LM_OBJS) $(LM_COBJS) 
	$(CXX) -g -o $@ $(LM_OBJS) $(LM_COBJS) $(LDFLAGS) 

unit_test:  unit_tester.cpp
	$(CXX) -g -o $@ $(UTESTOBJS) -lptlapack -lf77blas -lcblas -ltatlas -lgfortran
clean:
	rm *.o
