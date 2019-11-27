program TestEPDESolver

  ! gfortran -O3 -fdefault-real-8 -o execTest TestEPDESolver.f90 TestProblem1.f90 TestProblem2.f90 TestProblem3.f90 MultigridDataStructures.f90 preProcessingEPDE.f90 MatrixVectorProductMultigrid.f90 SmootherPML.f90 MultigridGridTransferOperators.f90 MultigridMethodsEllipticOperator.f90 extraSubroutines.f90
  use MultigridDataStructures
  use preProcessingEPDE
  use MatrixVectorProductMultigrid
  use SmootherPML
  use MultigridGridTransferOperators
  use MultigridMethodsEllipticOperator
  use TestProblem1
  use TestProblem2
  use TestProblem3

  character(len=*), parameter :: smoother = 'FC-GS'
  character(len=*), parameter  :: cycleType = 'W'
  integer, parameter :: TestProblem = 2
  integer, parameter :: numGrids = 7
  integer, parameter :: numCycles = 1
  integer, parameter :: preSweeps = 1
  integer, parameter :: postSweeps = 1
  integer, parameter :: totalIterations = 10
  integer, parameter :: orderRestriction = 2
  integer, parameter :: orderInterpolation = 2
  real(kind=8), parameter :: leftEndpoint = real(0,8)
  real(kind=8), parameter :: rightEndpoint = real(1,8)

  select case ( TestProblem )

  case (1)

    call MultigridSolverEllipticEquationA(evaluateParametersEPDETest1, &
      evaluatewaveNumberSquaredTest1, evaluateRHSTest1, &
      evaluateTrueSolutionTest1, smoother, cycleType, numGrids, &
      orderRestriction, orderInterpolation, totalIterations, preSweeps, &
      postSweeps, leftEndpoint, rightEndpoint)

  case (2)

    call MultigridSolverEllipticEquationA(evaluateParametersEPDETest2, &
      evaluatewaveNumberSquaredTest2, evaluateRHSTest2, &
      evaluateTrueSolutionTest2, smoother, cycleType, numGrids, &
      orderRestriction, orderInterpolation, totalIterations, preSweeps, &
      postSweeps, leftEndpoint, rightEndpoint)

  case (3)

    call MultigridSolverEllipticEquationA(evaluateParametersEPDETest3, &
      evaluatewaveNumberSquaredTest3, evaluateRHSTest3, &
      evaluateTrueSolutionTest3, smoother, cycleType, numGrids, &
      orderRestriction, orderInterpolation, totalIterations, preSweeps, &
      postSweeps, leftEndpoint, rightEndpoint)

  case default

    call MultigridSolverEllipticEquationA(evaluateParametersEPDETest1, &
      evaluatewaveNumberSquaredTest1, evaluateRHSTest1, &
      evaluateTrueSolutionTest1, smoother, cycleType, numGrids, &
      orderRestriction, orderInterpolation, totalIterations, preSweeps, &
      postSweeps, leftEndpoint, rightEndpoint)

  end select


end program TestEPDESolver
