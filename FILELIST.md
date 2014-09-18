List of files organized based on their functionality {#filelist}
====================================================

Prefix keys:
------------
- [AR]: Abtin Rahimian
- [DM]: Dhairya Malhotra
- [R] : Reviewed

Device:
-------
- CPUKernels.cc
- CPUKernels.h
- StokesTest.cc
- CudaApiGlobals.cc
- CudaApiGlobals.h
- CudaKernels.cu
- CudaKernels.h
- - -
- transpose_kernel.cu
- transpose_kernel.h
- transpose_wrapper.cu
- - -
- [R] Device.h
- [R] DeviceCPU.cc
- [AR] DeviceGPU.cc
- [R] DeviceTest.cc
- [R] DeviceTest.h

Basic:
------
- [R] Enums.cc
- [R] Enums.h
- [R] EnumsTest.cc
- - -
- [R] Error.cc
- [R] Error.h
- [R] ErrorTest.cc
- - -
- [R] Logger.cc
- [R] Logger.h
- [R] LoggerTest.cc
- - -
- [R] Array.cc
- [R] Array.h
- [R] ArrayTest.cc
- - -
- [R] Scalars.cc
- [R] Scalars.h
- [R] ScalarsTest.cc
- [R] ScalarContainerTest.h
- - -
- [R] Vectors.cc
- [R] Vectors.h
- [R] VectorsTest.h
- [R] VectorsTest.cc
- - -
- [AR] SHTMats.cc
- [AR] SHTMats.h
- [AR] SHTrans.cc
- [AR] SHTrans.h
- [AR] SHTransTest.cc
- - -
- [AR] OperatorsMats.cc
- [AR] OperatorsMats.h

Core:
-----
- EvalVelocity.cc
- EvalVelocity.h
- EvalVelocityTest.cc
- - -
- EvolveSurface.cc
- EvolveSurface.h
- EvolveSurfaceTest.cc
- EvolveSurfaceTestMultThread.cc
- - -
- InterfacialForce.cc
- InterfacialForce.h
- InterfacialVelocity.cc
- InterfacialVelocity.h
- - -
- MovePole.cc
- MovePole.h
- MovePoleTest.cc
- - -
- Parameters.cc
- Parameters.h
- - -
- [DM] Repartition.h       //repartition based on FMM box, random handling of cut vesicles, could be the same class as interaction(). make abstract and subclass
- [DM] Repartition.cc
- [DM] RepartitionTest.cc
- - -
- Surface.cc
- Surface.h
- - -
- SurfaceParams.cc
- SurfaceParams.h
- SurfaceTest.cc
- - -
- TimeStepperMult.h
- - -
- [DM] VesInteraction.h    //add single and double layer interface;imporve interface. Rename it to InteractionInterface() (abstract?) and subclass this for pvfmm
- [DM] VesInteraction.cc

Util:
-----
- BenchmarkInterfacialVel.cc
- - -
- BgFlow.cc
- BgFlow.h
- BgFlowBase.h
- - -
- BiCGStab.cc
- BiCGStab.h
- BiCGStabTest.cc
- - -
- BlasToyTest.cc
- HasAtlas.h
- HasBlas.h
- HasMkl.h
- VesBlas.h
- vesblas_cygwin.h
- - -
- CurvatureFlow.h
- - -
- DataIO.cc
- DataIO.h
- DataIO_templates.cc
- - -
- GLIntegrator.cc
- GLIntegrator.h
- - -
- HelperFuns.cc
- HelperFuns.h
- - -
- Monitor.cc
- Monitor.h
- - -
- ParserLite.cc
- ParserLite.h
- ParsingTest.cc
- ParsingTestInput.in
- - -
- anyoption.cc
- anyoption.h

Comments
--------
1. Should we have one MPI per device? I think Yes.
1. For both Repartition and VesInteraction add dummy and direct evaluation classes
