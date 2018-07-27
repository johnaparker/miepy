# Make option
precision   := double# double or quad
compilemode := release# release or debug

#Path
srcdir := ./

# Modules
objs = Parameters.o\
	utils_str.o\
	interpretator.o\
	parser.o\
	surface.o\
	parser_surface.o\
	AdditonTh.o\
	Allocation.o\
	BesLeg.o\
	Check.o\
	EFMED.o\
	GeomLib.o\
	GeomTrans.o\
	IncCoeff.o\
	InputOutput.o\
	Integr.o\
	Interp.o\
	MachParam.o\
	MatrixOp.o\
	MatrixQ.o\
	MatrixSolv.o\
	MatrixTrans.o\
	OpenMP.o\
	Proces1.o\
	Proces2.o\
	Proces3.o\
	PostProces1.o\
	PostProces2.o\
	PostProces3.o\
	Random.o\
	SCT.o\
	SCTAVRGSPH.o\
	SVWF.o\
	TAXSYM.o\
	TCOMP.o\
	TINHOM.o\
	TINHOM2SPH.o\
	TINHOMSPH.o\
	TINHOMSPHREC.o\
	TLAY.o\
	TMATRIX.o\
	TMULT.o\
	TMULT2SPH.o\
	TMULTSPH.o\
	TMULTSPHREC.o\
	TNONAXSYM.o\
	TPARTSUB.o\
	TPARTSUB3D.o\
	TPARTSUB3DFILM.o\
	TPARTSUB3DMULT.o\
	TPARTSUBFILM.o\
	TSPHERE.o
	
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
$(outfile): $(objs)
	$(fc) -o $(outfile) $(objs)

$(objs): %.o: $(srcdir)/%.f90
	$(fc) -c $< -o $@

.PHONE: clean

clean:
	rm -f *.o *.mod
