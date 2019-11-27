program TestEPDESolver

  ! gfortran -O3 -fdefault-real-8 -o execTest TestEPDESolver.f90 TestProblem1.f90 MultigridDataStructures.f90 preProcessingEPDE.f90 MatrixVectorProductMultigrid.f90 SmootherPML.f90 MultigridGridTransferOperators.f90 MultigridMethodsEllipticOperator.f90 extraSubroutines.f90
  use MultigridDataStructures
  use preProcessingEPDE
  use MatrixVectorProductMultigrid
  use SmootherPML
  use MultigridGridTransferOperators
  use MultigridMethodsEllipticOperator
  use TestProblem1

  character(len=*), parameter :: smoother = 'FC-GS'
  character(len=*), parameter  :: cycleType = 'W'
  integer, parameter :: numGrids = 7
  integer, parameter :: numCycles = 1
  integer, parameter :: preSweeps = 1
  integer, parameter :: postSweeps = 1
  integer, parameter :: totalIterations = 10
  integer, parameter :: orderRestriction = 2
  integer, parameter :: orderInterpolation = 2
  real(kind=8), parameter :: leftEndpoint = real(0,8)
  real(kind=8), parameter :: rightEndpoint = real(1,8)

  call MultigridSolverEllipticEquationA(evaluateParametersEPDETest1, &
    evaluatewaveNumberSquaredTest1, evaluateRHSTest1, &
    evaluateTrueSolutionTest1, smoother, cycleType, numGrids, &
    orderRestriction, orderInterpolation, totalIterations, preSweeps, &
    postSweeps, leftEndpoint, rightEndpoint)


end program TestEPDESolver
