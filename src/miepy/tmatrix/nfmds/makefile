# Make option
compilemode := release# release or debug
objdir := .
srcdir := .
exedir := .

# Platform-specific settings
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS: statically link GCC runtime libraries (except OpenMP)
    LDFLAGS ?= -static-libgfortran -static-libgcc -static-libquadmath
else ifneq ($(findstring MINGW,$(UNAME_S)),)
    # Windows (MSYS2/MinGW): statically link GCC runtime libraries
    LDFLAGS ?= -static-libgfortran -static-libgcc -static-libquadmath
else ifneq ($(findstring MSYS,$(UNAME_S)),)
    # Windows (MSYS2): statically link GCC runtime libraries
    LDFLAGS ?= -static-libgfortran -static-libgcc -static-libquadmath
endif

sources =  Parameters.f90\
           MachParam.f90\
           Allocation.f90\
           AdditonTh.f90\
           BesLeg.f90\
           Check.f90\
           SVWF.f90\
    	   GeomLib.f90\
           GeomTrans.f90\
           IncCoeff.f90\
           InputOutput.f90\
           Integr.f90\
           Interp.f90\
           MatrixSolv.f90\
           MatrixOp.f90\
           MatrixQ.f90\
           MatrixTrans.f90\
           PostProces1.f90\
           PostProces2.f90\
           PostProces3.f90\
           Proces1.f90\
           Proces2.f90\
           Proces3.f90\
           Random.f90\
           TAXSYM.f90\
           TCOMP.f90\
           TINHOM2SPH.f90\
           TINHOM.f90\
           TINHOMSPH.f90\
           TINHOMSPHREC.f90\
           TLAY.f90\
           TMULT2SPH.f90\
           TMULT.f90\
           TMULTSPH.f90\
           TMULTSPHREC.f90\
           TNONAXSYM.f90\
           TPARTSUB.f90\
           TSPHERE.f90\
           SCTAVRGSPH.f90\
           SCT.f90\
           EFMED.f90\
           TMATRIX.f90	

objs := $(sources:%.f90=$(objdir)/%.o)

#Compiler option
flags_1 :=
outfile := tmatrix

ifeq ($(compilemode),release)
  flags   := $(flags_1) -O3
else
  flags   := $(flags_1)
endif

# Allow overriding the Fortran compiler via FC variable
FC ?= gfortran
fc = $(FC) $(flags)

# Make commands
$(exedir)/$(outfile): $(objs)
	$(fc) -o $(exedir)/$(outfile) $(objs) $(LDFLAGS)

$(objdir)/%.o: $(srcdir)/%.f90
	$(fc) -c $< -o $@

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(exedir)/$(outfile)
