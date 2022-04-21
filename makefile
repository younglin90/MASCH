CCOMPLR = mpiicpc

CFLAGS = -std=c++17 -c -O3
# CFLAGS = -std=c++17 -c -O3 -fPIC -xCORE-AVX512 -ipo -no-prec-div
# CFLAGS = -std=c++17 -c -O3 -fPIC -ffast-math -hfp4 -ipo -no-prec-div


SHELL = /bin/bash
sp    = /-\|/
idx		=	0

LIBINCLUDE = \
             -Ilib/zlib\
             -Ilib/Metis\
             -Ilib/Scotch\
             -Ilib/amgcl\
             -Ilib\
             # -Ilib/boost\
             #-Ilib/HYPRE/include\
             #-Ilib/PETSc/include\
             # -I/home/yyl/petsc/arch-linux-c-debug/include\

LIBRARIES = \
            lib/zlib/libz.a\
            lib/Metis/libparmetis.a\
            lib/Metis/libmetis.a\
            lib/Scotch/libscotch.a\
            lib/Scotch/libscotcherr.a\
            # lib/PETSc/lib/libpetsc.so.3.14\
            #lib/HYPRE/lib/libHYPRE.a\
            #-Wl,-rpath,lib/PETSc/lib\
             # ./lib/Scotch/libscotchmetis.a\
             # -Wl,-rpath,/home/yyl/petsc/arch-linux-c-debug/lib\
             # -L/home/yyl/petsc/arch-linux-c-debug/lib\
             # -lm\
             # -lpetsc\
             # ./lib/PETSc/lib/libpetsc.so.3.14\
             # -Wl,-rpath,./lib/PETSc/lib\

SOURCES = src/mesh/mesh.cpp\
	      src/mesh/geometric.cpp\
	      src/mesh/partition.cpp\
	      src/mesh/distribute.cpp\
	      src/mesh/polyAMR.cpp\
	      src/mesh/polyRefine.cpp\
	      src/mesh/polyUnrefine.cpp\
	      src/mesh/reorder.cpp\
	      src/math/math.cpp\
	      src/math/gradient.cpp\
	      src/math/RCM.cpp\
	      src/controls/controls.cpp\
	      src/mpi/mpi.cpp\
	      src/solvers/solvers.cpp\
	      src/solvers/timestep.cpp\
	      src/solvers/norm.cpp\
	      src/solvers/incompressible/pressureBased.cpp\
	      src/solvers/incompressible/eqCoupled.cpp\
	      src/solvers/incompressible/eqMomentum.cpp\
	      src/solvers/incompressible/eqPressure.cpp\
	      src/solvers/incompressible/eqVolfrac.cpp\
	      src/solvers/compressible/densityBased.cpp\
	      src/solvers/compressible/compCoupled.cpp\
	      src/solvers/compressible/RHS.cpp\
	      src/solvers/compressible/massfrac.cpp\
	      src/solvers/compressible/pressure.cpp\
	      src/solvers/compressible/momentum.cpp\
	      src/solvers/compressible/energy.cpp\
	      src/solvers/compressible/flows.cpp\
	      src/solvers/compressible/coupled.cpp\
	      src/reconstruction/reconIncom.cpp\
	      src/reconstruction/reconComp.cpp\
	      src/reconstruction/NVD.cpp\
	      src/reconstruction/cellVarMinMax.cpp\
	      src/linearSolver/linearSolver.cpp\
	      src/linearSolver/solveAMGCL.cpp\
	      src/linearSolver/solvePETSc.cpp\
	      src/linearSolver/solveHYPRE.cpp\
	      src/convective/convective.cpp\
	      src/convective/surfaceNormalVelocity.cpp\
	      src/diffusion/diffusion.cpp\
	      src/source/gravity.cpp\
	      src/source/surfaceTension.cpp\
	      src/source/curvature.cpp\
	      src/eos/eos.cpp\
	      src/transport/transport.cpp\
	      src/turbulenceModels/LES.cpp\

SOURCES_SAVE = src/others/save.cpp\

SOURCES_LOAD = \
  src/others/load.cpp\
  src/others/load_settingFiles.cpp\
  src/others/load_vtuFiles.cpp\
  src/others/load_prim_vtuFiles.cpp\
  src/others/load_openFoam.cpp\

OBJECTS = $(SOURCES:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# Density based single time
EXE_CompDensitySingle = CompDensitySingle

SOURCES_CompDensitySingle = \
                    src/main/compressible/mainDensitySingle.cpp\
                    src/controls/controls.cpp\
                    src/mesh/mesh.cpp\
                    src/mesh/geometric.cpp\
                    src/math/math.cpp\
                    src/math/gradient.cpp\
                    src/mpi/mpi.cpp\

OBJECTS_CompDensitySingle = $(SOURCES_CompDensitySingle:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# # Density based dual time
# EXE_CompDensityDual = CompDensityDual

# OBJECTS_CompDensityDual = src/main/compressible/mainDensityDual.o $(OBJECTS)

# # Pressure based 
# EXE_IncomPressure = IncomPressure

# OBJECTS_IncomPressure = src/main/incompressible/mainPressure.o $(OBJECTS)

# # # hybrid based
# # EXE_CompHybrid = CompHybrid

# # OBJECTS_CompHybrid = src/main/compressible/hybridBased.o $(OBJECTS)

# # comp coupled based
# EXE_CompCoupled = CompCoupled

# OBJECTS_CompCoupled = src/main/compressible/mainFullCoupled.o $(OBJECTS)

# partitioning
EXE_PARTITION = Partition

SOURCES_PARTITION = src/main/partition/partitionMain.cpp\
                    src/main/partition/partition.cpp\
                    src/main/partition/parMetis.cpp\
                    src/others/controls.cpp\
                    src/others/mpi.cpp\
                    src/others/mesh.cpp\
                    src/others/log.cpp\
                    src/others/variables.cpp\
                    src/others/math.cpp\
                    src/others/setVarUDF.cpp\

OBJECTS_PARTITION = $(SOURCES_PARTITION:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# repartition
EXE_REPARTITION = Repartition

SOURCES_REPARTITION = src/main/utility/repartition.cpp\
                    src/others/controls.cpp\
                    src/others/controls_geometric.cpp\
                    src/others/mpi.cpp\
                    src/others/mesh.cpp\
                    src/others/log.cpp\
                    src/others/variables.cpp\
                    src/others/math.cpp\
                    src/others/setVarUDF.cpp\
                    src/others/repartitioning.cpp\
                    src/others/repartParMETIS.cpp\

OBJECTS_REPARTITION = $(SOURCES_REPARTITION:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# initialization
EXE_INITIAL = Initial

SOURCES_INITIAL = src/utility/initial.cpp\
                    src/controls/controls.cpp\
                    src/mesh/mesh.cpp\
                    src/mesh/geometric.cpp\
                    src/math/math.cpp\
                    src/math/gradient.cpp\
                    src/mpi/mpi.cpp\

OBJECTS_INITIAL = $(SOURCES_INITIAL:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# combine_mesh
EXE_CombineMesh = CombineMesh

SOURCES_CombineMesh = src/test/combineMesh.cpp\
                    src/controls/controls.cpp\
                    src/mesh/mesh.cpp\
                    src/mesh/geometric.cpp\
                    src/mesh/partition/partition.cpp\
                    src/mesh/partition/parMetis.cpp\
                    src/math/math.cpp\
                    src/math/gradient.cpp\
                    src/mpi/mpi.cpp\

OBJECTS_CombineMesh = $(SOURCES_CombineMesh:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)


SOURCES_Others = \
  src/others/variables.cpp\
  src/others/solvers.cpp\
  src/others/AMGCL.cpp\
  src/others/fvm.cpp\
  src/others/dpm.cpp\
  src/others/eulerianToLagrangian.cpp\
  src/others/math.cpp\
  src/others/gradient.cpp\
  src/others/curvature.cpp\
  src/others/controls.cpp\
  src/others/controls_geometric.cpp\
  src/others/mpi.cpp\
  src/others/mesh.cpp\
  src/others/log.cpp\
  src/others/polyAMR.cpp\
  src/others/polyRefine.cpp\
  src/others/polyUnrefine.cpp\
  src/others/repartitioning.cpp\
  src/others/repartParMETIS.cpp\

# test1
EXE_MultiComponent = MultiComponent

SOURCES_MultiComponent = \
  src/main/flow/compressible/multicomponent/multicomponent.cpp\
  src/main/flow/compressible/multicomponent/setAddiUDF.cpp\
  src/main/flow/compressible/multicomponent/setConvUDF.cpp\
  src/main/flow/compressible/multicomponent/setDiffUDF.cpp\
  src/main/flow/compressible/multicomponent/setFaceValUDF.cpp\
  src/main/flow/compressible/multicomponent/setGradUDF.cpp\
  src/main/flow/compressible/multicomponent/setHOUDF.cpp\
  src/main/flow/compressible/multicomponent/setOldVarUDF.cpp\
  src/main/flow/compressible/multicomponent/setSourUDF.cpp\
  src/main/flow/compressible/multicomponent/setTempUDF.cpp\
  src/main/flow/compressible/multicomponent/setVarUDF.cpp\
  src/main/flow/compressible/multicomponent/setTempStepUDF.cpp\
  src/main/flow/compressible/multicomponent/setUpdatePrimUDF.cpp\
  src/main/flow/compressible/multicomponent/setSegEq.cpp\
  src/main/flow/compressible/multicomponent/setBoundUDF.cpp\
  # src/main/flow/compressible/multicomponent/setMinMaxCellValuesUDF.cpp\

OBJECTS_MultiComponent = $(SOURCES_MultiComponent:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o) $(SOURCES_Others:.cpp=.o)



EXE_Potential = Potential

SOURCES_Potential = \
  src/main/flow/potential/potential.cpp\
  src/main/flow/potential/setAddiUDF.cpp\
  src/main/flow/potential/setConvUDF.cpp\
  src/main/flow/potential/setDiffUDF.cpp\
  src/main/flow/potential/setFaceValUDF.cpp\
  src/main/flow/potential/setGradUDF.cpp\
  src/main/flow/potential/setHOUDF.cpp\
  src/main/flow/potential/setOldVarUDF.cpp\
  src/main/flow/potential/setSourUDF.cpp\
  src/main/flow/potential/setTempUDF.cpp\
  src/main/flow/potential/setVarUDF.cpp\
  src/main/flow/potential/setTempStepUDF.cpp\
  src/main/flow/potential/setUpdatePrimUDF.cpp\
  src/main/flow/potential/setSegEq.cpp\
  src/main/flow/potential/setBoundUDF.cpp\

EXE_Laplace = Laplace

SOURCES_Laplace = \
  src/main/pde/laplace/laplace.cpp\
  src/main/pde/laplace/setAddiUDF.cpp\
  src/main/pde/laplace/setConvUDF.cpp\
  src/main/pde/laplace/setDiffUDF.cpp\
  src/main/pde/laplace/setFaceValUDF.cpp\
  src/main/pde/laplace/setGradUDF.cpp\
  src/main/pde/laplace/setHOUDF.cpp\
  src/main/pde/laplace/setOldVarUDF.cpp\
  src/main/pde/laplace/setSourUDF.cpp\
  src/main/pde/laplace/setTempUDF.cpp\
  src/main/pde/laplace/setVarUDF.cpp\
  src/main/pde/laplace/setTempStepUDF.cpp\
  src/main/pde/laplace/setUpdatePrimUDF.cpp\
  src/main/pde/laplace/setSegEq.cpp\
  src/main/pde/laplace/setBoundUDF.cpp\


OBJECTS_Laplace = $(SOURCES_Laplace:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o) $(SOURCES_Others:.cpp=.o)

EXE_INCOMP = InComp
FOLDER_INCOMP = src/main/flow/incompressible
SOURCES_INCOMP = \
  $(FOLDER_INCOMP)/main.cpp\
  $(FOLDER_INCOMP)/setVarUDF.cpp\
  $(FOLDER_INCOMP)/setSegEq.cpp\
  $(FOLDER_INCOMP)/setBoundUDF.cpp\
  $(FOLDER_INCOMP)/setGradUDF.cpp\
  $(FOLDER_INCOMP)/setAddiUDF.cpp\
  $(FOLDER_INCOMP)/setHOReconUDF.cpp\
  $(FOLDER_INCOMP)/setOldVarUDF.cpp\
  $(FOLDER_INCOMP)/setTimeStepUDF.cpp\
  $(FOLDER_INCOMP)/setTermsCellLoopUDF.cpp\
  $(FOLDER_INCOMP)/setTermsFaceLoopUDF.cpp\
  $(FOLDER_INCOMP)/setUpdatePrimUDF.cpp\

OBJECTS_INCOMP = $(SOURCES_INCOMP:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o) $(SOURCES_Others:.cpp=.o)

FOLDER_CompCoupled = src/main/flow/compCoupled

EXE_CompCoupledImplicit = CompCoupledImplicit
FOLDER_CompCoupledImplicit = $(FOLDER_CompCoupled)/implicit
SOURCES_CompCoupledImplicit = \
  $(FOLDER_CompCoupled)/main.cpp\
  $(FOLDER_CompCoupled)/setVarUDF.cpp\
  $(FOLDER_CompCoupled)/setBoundUDF.cpp\
  $(FOLDER_CompCoupled)/setGradUDF.cpp\
  $(FOLDER_CompCoupled)/setCurvatureUDF.cpp\
  $(FOLDER_CompCoupled)/setAddiUDF.cpp\
  $(FOLDER_CompCoupled)/setOldVarUDF.cpp\
  $(FOLDER_CompCoupled)/setTermsCellLoopUDF.cpp\
  $(FOLDER_CompCoupled)/setUpdatePrimUDF.cpp\
  $(FOLDER_CompCoupled)/setDPMUDF.cpp\
  $(FOLDER_CompCoupled)/setMinMaxCellValuesUDF.cpp\
  $(FOLDER_CompCoupled)/setMeanValues.cpp\
  $(FOLDER_CompCoupledImplicit)/setSegEq.cpp\
  $(FOLDER_CompCoupledImplicit)/setHOReconUDF.cpp\
  $(FOLDER_CompCoupledImplicit)/setTermsFaceLoopUDF.cpp\
  $(FOLDER_CompCoupledImplicit)/setTimeStepUDF.cpp\

OBJECTS_CompCoupledImplicit = $(SOURCES_CompCoupledImplicit:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o) $(SOURCES_Others:.cpp=.o)

EXE_CompCoupledExplicit = CompCoupledExplicit
FOLDER_CompCoupledExplicit = $(FOLDER_CompCoupled)/explicit
SOURCES_CompCoupledExplicit = \
  $(FOLDER_CompCoupled)/main.cpp\
  $(FOLDER_CompCoupled)/setVarUDF.cpp\
  $(FOLDER_CompCoupled)/setBoundUDF.cpp\
  $(FOLDER_CompCoupled)/setGradUDF.cpp\
  $(FOLDER_CompCoupled)/setCurvatureUDF.cpp\
  $(FOLDER_CompCoupled)/setAddiUDF.cpp\
  $(FOLDER_CompCoupled)/setOldVarUDF.cpp\
  $(FOLDER_CompCoupled)/setTermsCellLoopUDF.cpp\
  $(FOLDER_CompCoupled)/setUpdatePrimUDF.cpp\
  $(FOLDER_CompCoupled)/setDPMUDF.cpp\
  $(FOLDER_CompCoupled)/setMinMaxCellValuesUDF.cpp\
  $(FOLDER_CompCoupled)/setMeanValues.cpp\
  $(FOLDER_CompCoupledExplicit)/setSegEq.cpp\
  $(FOLDER_CompCoupledExplicit)/setHOReconUDF.cpp\
  $(FOLDER_CompCoupledExplicit)/setTermsFaceLoopUDF.cpp\
  $(FOLDER_CompCoupledExplicit)/setTimeStepUDF.cpp\

OBJECTS_CompCoupledExplicit = $(SOURCES_CompCoupledExplicit:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o) $(SOURCES_Others:.cpp=.o)

EXE_CompCoupledExplicitDPM = CompCoupledExplicitDPM
FOLDER_CompCoupledExplicitDPM = src/main/flow/compCoupledExplicitDPM
SOURCES_CompCoupledExplicitDPM = \
  $(FOLDER_CompCoupledExplicitDPM)/main.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setVarUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setSegEq.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setBoundUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setGradUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setCurvatureUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setAddiUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setHOReconUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setOldVarUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setTimeStepUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setTermsCellLoopUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setTermsFaceLoopUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setUpdatePrimUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setDPMUDF.cpp\
  $(FOLDER_CompCoupledExplicitDPM)/setMinMaxCellValuesUDF.cpp\

OBJECTS_CompCoupledExplicitDPM = $(SOURCES_CompCoupledExplicitDPM:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o) $(SOURCES_Others:.cpp=.o)


# # potential flow
# EXE_POTENTIAL = Potential

# SOURCES_POTENTIAL = src/utility/potential.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\
                    # src/solvers/solvers.cpp\
                    # src/eos/eos.cpp\
                    # src/linearSolver/solveAMGCL.cpp\
                    # src/linearSolver/solvePETSc.cpp\
                    # src/reconstruction/reconIncom.cpp\
                    # src/reconstruction/reconComp.cpp\
                    # src/reconstruction/NVD.cpp\

# OBJECTS_POTENTIAL = $(SOURCES_POTENTIAL:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)


# # laplace equations
# EXE_LAPLACE = Laplace

# SOURCES_LAPLACE = src/utility/laplace.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\
                    # src/solvers/solvers.cpp\
                    # src/eos/eos.cpp\
                    # src/linearSolver/solveAMGCL.cpp\
                    # src/reconstruction/reconIncom.cpp\
                    # src/reconstruction/reconComp.cpp\
                    # src/reconstruction/NVD.cpp\

# OBJECTS_LAPLACE = $(SOURCES_LAPLACE:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)


# # advection equations
# EXE_ADVECTION = Advection

# SOURCES_ADVECTION = src/utility/advection.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\
                    # src/solvers/solvers.cpp\
                    # src/eos/eos.cpp\
                    # src/linearSolver/solveAMGCL.cpp\
                    # src/reconstruction/reconIncom.cpp\
                    # src/reconstruction/reconComp.cpp\
                    # src/reconstruction/NVD.cpp\
                    # src/mesh/polyAMR.cpp\
                    # src/mesh/polyRefine.cpp\
                    # src/mesh/polyUnrefine.cpp\

# OBJECTS_ADVECTION = $(SOURCES_ADVECTION:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)


# # compression volume fraction equations
# EXE_CompressVF = CompressVF

# SOURCES_CompressVF = src/utility/compressVF.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\
                    # src/solvers/solvers.cpp\
                    # src/eos/eos.cpp\
                    # src/linearSolver/solveAMGCL.cpp\
                    # src/reconstruction/reconIncom.cpp\
                    # src/reconstruction/reconComp.cpp\
                    # src/reconstruction/NVD.cpp\
                    # src/mesh/polyAMR.cpp\
                    # src/mesh/polyRefine.cpp\
                    # src/mesh/polyUnrefine.cpp\

# OBJECTS_CompressVF = $(SOURCES_CompressVF:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)


# # calc. SMD from volume fraction
# EXE_SMD = CalcSMD

# SOURCES_SMD = src/utility/calcSMD.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\
                    # src/solvers/solvers.cpp\
                    # src/eos/eos.cpp\
                    # src/linearSolver/solveAMGCL.cpp\
                    # src/reconstruction/reconIncom.cpp\
                    # src/reconstruction/reconComp.cpp\
                    # src/reconstruction/NVD.cpp\
                    # src/mesh/polyAMR.cpp\
                    # src/mesh/polyRefine.cpp\
                    # src/mesh/polyUnrefine.cpp\

# OBJECTS_SMD = $(SOURCES_SMD:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# # gradient
# EXE_GRADIENT = Gradient

# SOURCES_GRADIENT = src/utility/gradient.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\
                    # src/solvers/solvers.cpp\
                    # src/eos/eos.cpp\
                    # src/linearSolver/solveAMGCL.cpp\
                    # src/reconstruction/reconIncom.cpp\
                    # src/reconstruction/reconComp.cpp\
                    # src/reconstruction/NVD.cpp\

# OBJECTS_GRADIENT = $(SOURCES_GRADIENT:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)



# # MapField
# EXE_MapField = MapField

# SOURCES_MapField = src/utility/mapField.cpp\
                    # src/mesh/mesh.cpp\
                    # src/mesh/geometric.cpp\
                    # src/math/math.cpp\
                    # src/math/gradient.cpp\
                    # src/mpi/mpi.cpp\
                    # src/controls/controls.cpp\

# OBJECTS_MapField = $(SOURCES_MapField:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

# # # MeshAMR
# # EXE_MeshAMR = MeshAMR

# # # SOURCES_MeshAMR = src/utility/meshAMR.cpp\
                    # # # src/mesh/build.cpp\
                    # # # src/mesh/load.cpp\
                    # # # src/mesh/save.cpp\
                    # # # src/mesh/geometric.cpp\
                    # # # src/math/math.cpp\
                    # # # src/controls/build.cpp\

# # OBJECTS_MeshAMR = src/utility/meshAMR.o $(SOURCES:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)
# # # OBJECTS_MeshAMR = $(SOURCES_MeshAMR:.cpp=.o)

# # ExtractData
# EXE_ExtractData = ExtractData

# OBJECTS_ExtractData = src/utility/extractData.o $(SOURCES:.cpp=.o) $(SOURCES_SAVE:.cpp=.o) $(SOURCES_LOAD:.cpp=.o)

COTEXT  = "\033[1;31m Compiling\033[0m\033[1m $< \033[0m"

# all : $(EXE)

# EXE_ALL = $(EXE_REPARTITION)
# OBJECTS_ALL = $(OBJECTS_REPARTITION)

EXE_ALL = $(EXE_PARTITION) $(EXE_REPARTITION) $(EXE_CompCoupledImplicit) $(EXE_CompCoupledExplicit)
OBJECTS_ALL = $(OBJECTS_PARTITION) $(OBJECTS_REPARTITION) $(OBJECTS_CompCoupledImplicit) $(OBJECTS_CompCoupledExplicit)

# EXE_ALL = $(EXE_PARTITION) $(EXE_MultiComponent) $(EXE_Potential) 

# OBJECTS_ALL = $(OBJECTS_PARTITION) $(OBJECTS_MultiComponent) $(OBJECTS_Potential)

#EXE_ALL = $(EXE_CompDensitySingle) $(EXE_PARTITION) $(EXE_test1) $(EXE_CombineMesh)

#OBJECTS_ALL = $(OBJECTS_CompDensitySingle) $(OBJECTS_PARTITION) $(OBJECTS_test1) $(OBJECTS_CombineMesh)

# EXE_ALL = $(EXE_CompDensitySingle) $(EXE_CompDensityDual) $(EXE_IncomPressure) $(EXE_CompCoupled) $(EXE_PARTITION) $(EXE_INITIAL) $(EXE_POTENTIAL) $(EXE_LAPLACE) $(EXE_ADVECTION) $(EXE_CompressVF) $(EXE_SMD) $(EXE_GRADIENT) $(EXE_MapField) $(EXE_MeshAMR) $(EXE_ExtractData)

# OBJECTS_ALL = $(OBJECTS_CompDensitySingle) $(OBJECTS_CompDensityDual) $(OBJECTS_IncomPressure) $(OBJECTS_CompCoupled) $(OBJECTS_PARTITION) $(OBJECTS_INITIAL) $(OBJECTS_POTENTIAL) $(OBJECTS_LAPLACE) $(OBJECTS_ADVECTION) $(OBJECTS_CompressVF) $(OBJECTS_SMD) $(OBJECTS_GRADIENT) $(OBJECTS_MapField) $(OBJECTS_MeshAMR)

all : $(EXE_ALL)
# all : $(EXE_LOAD)

$(EXE_CompDensitySingle) : $(OBJECTS_CompDensitySingle)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompDensitySingle) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_density_single CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_PARTITION) : $(OBJECTS_PARTITION)
	@$(CCOMPLR) -o $@ $(OBJECTS_PARTITION) $(LIBRARIES)
	@echo -e "\033[1;31m PARTITION CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_REPARTITION) : $(OBJECTS_REPARTITION)
	@$(CCOMPLR) -o $@ $(OBJECTS_REPARTITION) $(LIBRARIES)
	@echo -e "\033[1;31m REPARTITION CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_INITIAL) : $(OBJECTS_INITIAL)
	@$(CCOMPLR) -o $@ $(OBJECTS_INITIAL) $(LIBRARIES)
	@echo -e "\033[1;31m INITIAL CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CombineMesh) : $(OBJECTS_CombineMesh)
	@$(CCOMPLR) -o $@ $(OBJECTS_CombineMesh) $(LIBRARIES)
	@echo -e "\033[1;31m CombineMesh CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_MultiComponent) : $(OBJECTS_MultiComponent)
	@$(CCOMPLR) -o $@ $(OBJECTS_MultiComponent) $(LIBRARIES)
	@echo -e "\033[1;31m MultiComponent CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_Potential) : $(OBJECTS_Potential)
	@$(CCOMPLR) -o $@ $(OBJECTS_Potential) $(LIBRARIES)
	@echo -e "\033[1;31m Potential CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_Laplace) : $(OBJECTS_Laplace)
	@$(CCOMPLR) -o $@ $(OBJECTS_Laplace) $(LIBRARIES)
	@echo -e "\033[1;31m Laplace CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_INCOMP) : $(OBJECTS_INCOMP)
	@$(CCOMPLR) -o $@ $(OBJECTS_INCOMP) $(LIBRARIES)
	@echo -e "\033[1;31m InComp CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompCoupledImplicit) : $(OBJECTS_CompCoupledImplicit)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompCoupledImplicit) $(LIBRARIES)
	@echo -e "\033[1;31m " $(EXE_CompCoupledImplicit) " CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompCoupledExplicit) : $(OBJECTS_CompCoupledExplicit)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompCoupledExplicit) $(LIBRARIES)
	@echo -e "\033[1;31m " $(EXE_CompCoupledExplicit) " CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompCoupledExplicitDPM) : $(OBJECTS_CompCoupledExplicitDPM)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompCoupledExplicitDPM) $(LIBRARIES)
	@echo -e "\033[1;31m " $(EXE_CompCoupledExplicitDPM) " CODE compile/link complete \033[0m" | tee -a make.log

ifeq ("x","y")
$(EXE_CompDensityDual) : $(OBJECTS_CompDensityDual)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompDensityDual) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_density_dual CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_IncomPressure) : $(OBJECTS_IncomPressure)
	@$(CCOMPLR) -o $@ $(OBJECTS_IncomPressure) $(LIBRARIES)
	@echo -e "\033[1;31m Incom_pressure CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompCoupled) : $(OBJECTS_CompCoupled)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompCoupled) $(LIBRARIES)
	@echo -e "\033[1;31m Comp_Coupled CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_POTENTIAL) : $(OBJECTS_POTENTIAL)
	@$(CCOMPLR) -o $@ $(OBJECTS_POTENTIAL) $(LIBRARIES)
	@echo -e "\033[1;31m POTENTIAL CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_LAPLACE) : $(OBJECTS_LAPLACE)
	@$(CCOMPLR) -o $@ $(OBJECTS_LAPLACE) $(LIBRARIES)
	@echo -e "\033[1;31m LAPLACE CODE compile/link complete \033[0m" | tee -a make.log


$(EXE_ADVECTION) : $(OBJECTS_ADVECTION)
	@$(CCOMPLR) -o $@ $(OBJECTS_ADVECTION) $(LIBRARIES)
	@echo -e "\033[1;31m ADVECTION CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_CompressVF) : $(OBJECTS_CompressVF)
	@$(CCOMPLR) -o $@ $(OBJECTS_CompressVF) $(LIBRARIES)
	@echo -e "\033[1;31m ADVECTION CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_SMD) : $(OBJECTS_SMD)
	@$(CCOMPLR) -o $@ $(OBJECTS_SMD) $(LIBRARIES)
	@echo -e "\033[1;31m calcSMD CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_GRADIENT) : $(OBJECTS_GRADIENT)
	@$(CCOMPLR) -o $@ $(OBJECTS_GRADIENT) $(LIBRARIES)
	@echo -e "\033[1;31m GRADIENT CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_MapField) : $(OBJECTS_MapField)
	@$(CCOMPLR) -o $@ $(OBJECTS_MapField) $(LIBRARIES)
	@echo -e "\033[1;31m MapField CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_MeshAMR) : $(OBJECTS_MeshAMR)
	@$(CCOMPLR) -o $@ $(OBJECTS_MeshAMR) $(LIBRARIES)
	@echo -e "\033[1;31m MeshAMR CODE compile/link complete \033[0m" | tee -a make.log

$(EXE_ExtractData) : $(OBJECTS_ExtractData)
	@$(CCOMPLR) -o $@ $(OBJECTS_ExtractData) $(LIBRARIES)
	@echo -e "\033[1;31m ExtractData CODE compile/link complete \033[0m" | tee -a make.log
endif

%.o : %.cpp
	@echo -e $(COTEXT) | tee -a make.log
	@$(CCOMPLR) $(CFLAGS) $(LIBINCLUDE) $< -o $@

clean:
	@echo -e "\033[1;31m deleting objects \033[0m" | tee make.log
	@rm -fr $(OBJECTS_ALL) $(EXE_ALL) make.log *.o
