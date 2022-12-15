#--------NLD----------

#CXX = g++
#CXXFLAGS = -DBLAS -O2 -Wall -I/usr/nld/gsl/include -I/usr/nld/atlas/include -I/usr/nld/netcdf-3.6.3/include
#LIBS = -L/usr/nld/gsl/lib -L/usr/nld/atlas/lib -L/usr/nld/netcdf-3.6.3/lib -lgsl -lcblas -latlas -lnetcdf_c++ -lnetcdf

#might use the additional CXXFLAG -static-libstdc++ to compile from within matlab

#------------------------------------------w/o BLAS-------------------------------------------------------

#CXX = g++
#CXXFLAGS = -O2 -Wall -I/opt/local/include
#LIBS = -L/opt/local/lib -lgsl -lgslcblas -lnetcdf_c++ -lnetcdf

#-------------------------------------------w/ BLAS by ATLAS------------------------------------------------------------

# use BLAS routines for faster orthogonalizations
#-DBLAS for BLAS2 routines in orthogonalization (matrix - vector multiplications)
# use the platform-specific optimized CBLAS routines provided by ATLAS (might need to add the include path)

CXX = g++
CXXFLAGS = -DBLAS -O2 -Wall -I/opt/local/include
LIBS = -L/opt/local/lib -lgsl -lcblas -latlas -lnetcdf_c++ -lnetcdf

#-------------------------------------------w/ BLAS by GSL------------------------------------------------------------

# or, if ATLAS is not installed use at least the CBLAS routines provided by GSL

#CXX = g++
#CXXFLAGS = -DBLAS -O2 -Wall -I/opt/local/include
#LIBS = -L/opt/local/lib -lgsl -lgslcblas -lnetcdf_c++ -lnetcdf

#-------------------------------------PARALLEL w/ BLAS----------------------------------------------------------

#use -DPAR for parallel calculation with mpi

#CXX = mpic++
#CXXFLAGS = -DPAR -DBLAS -O2 -Wall -I/opt/local/include
#LIBS = -L/opt/local/lib -lgsl -lcblas -latlas -lnetcdf_c++ -lnetcdf

#---------------------------------------------------------------------------------------------------------------



OBJS = LEanalysis.o LEons.o LEnetwork.o LEneuron.o LEinputoutput.o LEphaseneuron.o LErapidtheta.o LElif.o LEtype1type2.o LEclif.o LEtwoDlinear.o

vpath %.cpp	code

TARGET = LEquipe

all:	$(TARGET)

LEquipe: $(OBJS) LEquipe.o
		 $(CXX) $(CXXFLAGS) -o $@ $@.o $(OBJS) $(LIBS)			

clean: 
		 rm -f *.o
		 
clean_data:
		 rm -f data/*.nc

#on mac run: valgrind --dsymutil=yes ./LEquipe  for line numbers
# valgrind --dsymutil=yes --tool=cachegrind  for profiling
# performance test on nanna04 with N=500 r=10 theta neurons, K=50, J=-1, f=1, TC=5s
#			Sequential	openMPI: 1CPU	2CPU	4CPU		
# no BLAS		1507			 1524	775		402
# GSL BLAS		692				 690	372		202
# ATLAS BLAS	362				 359	197		119

