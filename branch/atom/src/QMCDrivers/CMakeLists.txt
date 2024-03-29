####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(QMCDRIVERS 
  SimpleFixedNodeBranch.cpp
  BranchIO.cpp
  QMCDriver.cpp
  QMCOptimize.cpp
  QMCLinearOptimize.cpp
  QMCCorrelatedSamplingLinearOptimize.cpp
  QMCFixedSampleLinearOptimize.cpp
  QMCCSLinearOptimizeWFmanagerOMP.cpp
  QMCCostFunctionBase.cpp
  WaveFunctionTester.cpp
  WalkerControlBase.cpp
  CloneManager.cpp
  QMCUpdateBase.cpp
  VMC/VMCUpdatePbyP.cpp
  VMC/VMCUpdateAll.cpp
  VMC/VMCUpdateAllWithIons.cpp
  VMC/VMCFactory.cpp
  DMC/DMCOMP.cpp
  DMC/DMCUpdateAll.cpp
  DMC/DMCUpdateAllWithIons.cpp
  DMC/DMCUpdatePbyP.cpp
  DMC/DMCUpdatePbyPFast.cpp
  DMC/DMCFactory.cpp
  DMC/WalkerControlFactory.cpp
  DMC/WalkerReconfiguration.cpp
  RMC/RMCSingleOMP.cpp
  RMC/RMCUpdatePbyP.cpp
  RMC/RMCUpdateAll.cpp
  RMC/RMCFactory.cpp
  ../Estimators/LocalEnergyEstimator.cpp
  ../Estimators/RMCLocalEnergyEstimator.cpp
  ../Estimators/LocalEnergyEstimatorHDF.cpp
  ../Estimators/EstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  ../Estimators/CollectablesEstimator.cpp
)

IF(QMC_CUDA)
   SET(QMCDRIVERS ${QMCDRIVERS}
   VMC/VMC_CUDA.cpp
   DMC/DMC_CUDA.cpp
   QMCCostFunctionCUDA.cpp)
ENDIF(QMC_CUDA)

# IF(QMC_BUILD_LEVEL GREATER 2)
#    SET(QMCDRIVERS ${QMCDRIVERS}
#    EE/EEFactory.cpp
#    EE/VMCMultiRenyiOMP.cpp
#    EE/VMCRenyiOMP.cpp
#    EE/QMCRenyiUpdateBase.cpp
#    EE/RenyiUpdatePbyP.cpp
#    )
# ENDIF(QMC_BUILD_LEVEL GREATER 2)

  SET(QMCDRIVERS ${QMCDRIVERS} 
      VMC/VMCSingleOMP.cpp
      VMC/VMCLinearOptOMP.cpp
      QMCCostFunctionOMP.cpp
     )

IF(HAVE_MPI) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 


IF(QMC_BUILD_LEVEL GREATER 1)
  SET(QMCDRIVERS ${QMCDRIVERS}
#     PolymerEstimator.cpp
#     MultiChain.cpp
#     RQMCMultiple.cpp
    ../Estimators/CSEnergyEstimator.cpp
#     ../Estimators/CSPolymerEstimator.cpp
  )
# REMOVE broken stuff
#    ../Estimators/RQMCMultipleEstimator.cpp
#    SpaceWarp.cpp
#    RQMCMultiWarp.cpp
#    VMC/VMCMultipleWarp.cpp
#    VMC/VMCPbyPMultiWarp.cpp
#    CorrelatedSampling/CSUpdateBase.cpp
#    CorrelatedSampling/CSVMC.cpp
#    CorrelatedSampling/CSVMCUpdateAll.cpp
#    CorrelatedSampling/CSVMCUpdatePbyP.cpp
#  IF(NOT QMC_COMPLEX)
#    SET(QMCDRIVERS ${QMCDRIVERS}
#    RQMCMultiplePbyP.cpp
#    )
#  ENDIF(NOT QMC_COMPLEX)

ENDIF(QMC_BUILD_LEVEL GREATER 1)


####################################
# create libqmc 
####################################
ADD_LIBRARY(qmcdriver ${QMCDRIVERS})
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver STATIC ${QMCDRIVERS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver SHARED ${QMCDRIVERS})
#ENDIF(QMC_BUILD_STATIC)
