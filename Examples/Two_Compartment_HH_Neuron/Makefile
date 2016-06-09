# Makefile for two_compartment IPOPT problem

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = two_compartment_cpp

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = two_compartmentminAone_main.o \
       two_compartmentminAone_nlp.o #\

# CHANGEME: Additional libraries
#ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

# CHANGEME: Directory to the sources for the (example) problem definition
# files
#SRCDIR = 
#VPATH = 

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/usr/local/Ipopt/lib

# Directory with header files
IPOPTINCDIR = ${prefix}/include/coin

# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}
prefix = /usr/local/Ipopt

# Libraries necessary to link with IPOPT
#LIBS = -L$(IPOPTLIBDIR) -lipopt @IPADDLIBS@
LIBS = `PKG_CONFIG_PATH=/usr/local/Ipopt/lib64/pkgconfig:/usr/local/Ipopt/lib/pkgconfig:/usr/local/Ipopt/share/pkgconfig: pkg-config --libs ipopt`
# Necessary Include dirs (we use the CYGPATH_W variables to allow
# compilation with Windows compilers)
INCL =  `PKG_CONFIG_PATH=/usr/local/Ipopt/lib64/pkgconfig:/usr/local/Ipopt/lib/pkgconfig:/usr/local/Ipopt/share/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
