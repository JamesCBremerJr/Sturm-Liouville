#-fdefault-real-8
#-freal-8-real-10

# Choose the compiler and set compiler options

CPPCOMP  = g++
CPPOPTS  = -O3 

# FCOMP    = ifort
# FOPTS    = -O3 -i8 -openmp -fno-inline
# LDOPT =  

FCOMP    = gfortran-8
FOPTS    = -march=native -w -fopenmp  -Ofast -fdefault-integer-8  -fexternal-blas
LDOPT    = -lblas -llapack 

# Set the list of programs to compile

PROGRAMS = test_chebyshev test_legendre test_tensor \
           test_odesolve test_kummer test_kummer_sturm

# Compile all of the test programs and the library

all	        : clean $(PROGRAMS) 

# List the dependencies for each module's test program


KUMMER_STURM_FILES  =     utils.o                       \
                         chebyshev.o                   \
                         tensor.o                      \
                         odesolve.o                    \
                         kummer.o                      \
                         kummer_sturm.o


KUMMER_FILES           = utils.o                       \
                         chebyshev.o                   \
                         tensor.o                      \
                         odesolve.o                    \
                         kummer.o

ODESOLVE_FILES         = utils.o                       \
                         chebyshev.o                   \
                         odesolve.o

TENSOR3D_FILES         = utils.o                       \
                         chebyshev.o                   \
                         tensor3d.o

TENSOR_FILES           = utils.o                       \
                         chebyshev.o                   \
                         tensor.o

PROLATES_FILES         = utils.o                       \
                         chebyshev.o                   \
                         prolates.o

LEGENDRE_FILES         = utils.o                       \
                         legendre.o

CHEBYSHEV_FILES        = utils.o                       \
                         chebyshev.o 





###################################################################################

test_kummer_prolates.o : $(KUMMER_PROLATES_FILES) test_kummer_prolates.f90
test_kummer_prolates   : $(KUMMER_PROLATES_FILES) test_kummer_prolates.o

generate_alegendre_roots.o: $(GENERATE_ALEGENDRE_ROOTS_FILES) generate_alegendre_roots.f90
generate_alegendre_roots  : $(GENERATE_ALEGENDRE_ROOTS_FILES) generate_alegendre_roots.o

test_kummer_alegendre.o: $(KUMMER_ALEGENDRE_FILES) test_kummer_alegendre.f90
test_kummer_alegendre  : $(KUMMER_ALEGENDRE_FILES) test_kummer_alegendre.o

test_kummer_legendre.o : $(KUMMER_LEGENDRE_FILES) test_kummer_legendre.f90
test_kummer_legendre   : $(KUMMER_LEGENDRE_FILES) test_kummer_legendre.o

test_kummer_jacobi.o   : $(KUMMER_JACOBI_FILES) test_kummer_jacobi.f90
test_kummer_jacobi     : $(KUMMER_JACOBI_FILES) test_kummer_jacobi.o


test_kummer_jacobi2.o   : $(KUMMER_JACOBI2_FILES) test_kummer_jacobi2.f90
test_kummer_jacobi2     : $(KUMMER_JACOBI2_FILES) test_kummer_jacobi2.o

test_kummer_sturm.o    : $(KUMMER_STURM_FILES) test_kummer_sturm.f90
test_kummer_sturm      : $(KUMMER_STURM_FILES) test_kummer_sturm.o

test_kummer_bessel.o   : $(KUMMER_BESSEL_FILES) test_kummer_bessel.f90
test_kummer_bessel     : $(KUMMER_BESSEL_FILES) test_kummer_bessel.o

jacobi_quad.o          : $(JACOBI_QUAD_FILES) jacobi_quad.f90
jacobi_quad            : $(JACOBI_QUAD_FILES) jacobi_quad.o

laguerre_quad.o        : $(LAGUERRE_QUAD_FILES) laguerre_quad.f90
laguerre_quad          : $(LAGUERRE_QUAD_FILES) laguerre_quad.o

test_kummer.o          : $(KUMMER_FILES) test_kummer.f90
test_kummer            : $(KUMMER_FILES) test_kummer.o

test_odesolve.o        : $(ODESOLVE_FILES) test_odesolve.f90
test_odesolve          : $(ODESOLVE_FILES) test_odesolve.o

test_tensor3d.o        : $(TENSOR3D_FILES) test_tensor3d.f90
test_tensor3d          : $(TENSOR3D_FILES) test_tensor3d.o

test_tensor.o          : $(TENSOR_FILES) test_tensor.f90
test_tensor            : $(TENSOR_FILES) test_tensor.o

test_prolates.o        : $(PROLATES_FILES) test_prolates.f90
test_prolates          : $(PROLATES_FILES) test_prolates.o

test_legendre.o        : $(LEGENDRE_FILES) test_legendre.f90
test_legendre          : $(LEGENDRE_FILES) test_legendre.o

test_chebyshev.o       : $(CHEBYSHEV_FILES) test_chebyshev.f90
test_chebyshev         : $(CHEBYSHEV_FILES) test_chebyshev.o


# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 

%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out
	rm -f $(PROGRAMS)
	rm -f gn??? gn???.dat
