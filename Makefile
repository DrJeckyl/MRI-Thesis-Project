F90_SRCS +=\
Funcs.f90 \
EquationsOfMotion.f90 \
RungeKutta4.f90 \
nr.f90 \
TypesAndDefs.f90 \

OBJS += \
Funcs.o \
EquationsOfMotion.o \
RungeKutta4.o \
nr.o \
TypesAndDefs.o \

FLAGS = -L/usr/lib/X11
FLAGS += -L/usr/local/pgplot
FLAGS += -L$(MKLROOT)/lib/intel64
FLAGS += -I$(MKLROOT)/include/fftw
FLAGS += -openmp 
FLAGS += -O3
FLAGS += -fpe0
FLAGS += -march=corei7

FLAGS += -heap-arrays 4096
FLAGS += -no-prec-div
FLAGS += -xHost
FLAGS += -zero

#FLAGS += -parallel
#FLAGS += -par-num-threads=8

#Must compile with -r8 so that the fft algorthm works properly
FLAGS += -r8
#####################################################

#FLAGS += -axSSE4.1 
#FLAGS += -axSSSE3
#FLAGS += -fp-stack-check
#FLAGS += -fno-inline-functions
FLAGS += -traceback
FLAGS += -checkbounds

LIBS = -lpgplot
LIBS += -lX11
LIBS += -lpng
LIBS += -lmkl_rt

COMP = ifort
#COMP = gfortran

#all: Main


nr.o: nr.f90 $(OBJS)
Funcs.o: Funcs.f90 $(OBJS)
TypesAndDefs.o: TypesAndDefs.f90 $(OBJS)
EquationsOfMotion.o: EquationsOfMotion.f90 $(OBJS)
RungeKutta4.o: RungeKutta4.f90 $(OBJS)

%: %.f90
	$(COMP) $(FLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMP) $(FLAGS) -c -o $@ $< $(LIBS)

Main: Main.o $(OBJS)
	$(COMP) $(FLAGS) -o Main Main.o $(OBJS) $(LIBS)


clean: 
	-rm -f *.o *.mod *.MOD


