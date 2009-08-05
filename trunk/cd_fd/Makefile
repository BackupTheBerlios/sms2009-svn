
FC=gfortran
EXE1=cd.exe
EXE2=cd_implicit.exe
EXES=$(EXE1) $(EXE2)
OBJS=derivatives.o solvers.o

all: $(EXES)

$(EXE1) : cd.F90 derivatives.o
	$(FC) $< derivatives.o -o $@

$(EXE2) : cd_implicit.F90 $(OBJS)
	$(FC) $< $(OBJS) -o $@

$(OBJS) : %.o: %.F90
	$(FC) -c $< -o $@ 

.PHONY: all clean spotless

clean:
	\rm -f $(OBJS)

spotless:
	\rm -f $(OBJS) $(EXES)
