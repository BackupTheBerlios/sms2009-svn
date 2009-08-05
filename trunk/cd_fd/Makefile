
FC=gfortran
EXE=cd.exe

all: $(EXE)

$(EXE) : cd.F90 derivatives.o
	$(FC) $< derivatives.o -o $@

derivatives.o: 
	$(FC) -c derivatives.F90 -o derivatives.o

.PHONY: all clean spotless

clean:
	\rm -f derivatives.o cd.o

spotless:
	\rm -f derivatives.o cd.o $(EXE)
