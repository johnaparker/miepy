# Make option
precision   := double  # double or quad
compilemode := release # release or debug
objdir := .
srcdir := .
exedir := .

# Modules
sources = Parameters.f90\
          utils_str.f90\
          interpretator.f90\
          parser.f90\
          surface.f90\
          parser_surface.f90\
          AdditonTh.f90\
          Allocation.f90\
          BesLeg.f90\
          Check.f90\
          EFMED.f90\
          GeomLib.f90\
          GeomTrans.f90\
          IncCoeff.f90\
          InputOutput.f90\
          Integr.f90\
          Interp.f90\
          MachParam.f90\
          MatrixOp.f90\
          MatrixQ.f90\
          MatrixSolv.f90\
          MatrixTrans.f90\
          OpenMP.f90\
          Proces1.f90\
          Proces2.f90\
          Proces3.f90\
          PostProces1.f90\
          PostProces2.f90\
          PostProces3.f90\
          Random.f90\
          SCT.f90\
          SCTAVRGSPH.f90\
          SVWF.f90\
          TAXSYM.f90\
          TCOMP.f90\
          TINHOM.f90\
          TINHOM2SPH.f90\
          TINHOMSPH.f90\
          TINHOMSPHREC.f90\
          TLAY.f90\
          TMATRIX.f90\
          TMULT.f90\
          TMULT2SPH.f90\
          TMULTSPH.f90\
          TMULTSPHREC.f90\
          TNONAXSYM.f90\
          TPARTSUB.f90\
          TPARTSUB3D.f90\
          TPARTSUB3DFILM.f90\
          TPARTSUB3DMULT.f90\
          TPARTSUBFILM.f90\
          TSPHERE.f90

objs := $(sources:%.f90=$(objdir)/%.o)

#Compiler option
flags_1 := -cpp

#Compiler option
ifeq ($(precision), quad)
  flags_1 := $(flags_1) -DPRECISION_QUAD
  outfile := tmatrix_extended
else
  flags_1 := $(flags_1)
  outfile := tmatrix
endif

ifeq ($(compilemode),release)
  flags   := $(flags_1) -fopenmp -O3
else
  flags   := $(flags_1)
endif

fc = gfortran $(flags)

# Make commands
$(exedir)/$(outfile): $(objs)
	$(fc) -o $(exedir)/$(outfile) $(objs)

$(objdir)/%.o: $(srcdir)/%.f90
	$(fc) -c $< -o $@

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(exedir)/$(outfile)
