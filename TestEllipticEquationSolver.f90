program TestEllipticEquationSolver

  ! gfortran -fdefault-real-8 -o execTest TestEllipticEquationSolver.f90 FunctionsEllipticEquation.f90 MultigridDataStructures.f90 preProcessingEllipticOperator.f90 MatrixVectorProductMultigrid.f90 SmootherPML.f90 MultigridGridTransferOperators.f90 MultigridMethodsEllipticOperator.f90
  use FunctionsEllipticEquation
  use MultigridDataStructures
  use preProcessingEllipticOperator
  use MatrixVectorProductMultigrid
  use SmootherPML
  use MultigridGridTransferOperators
  use MultigridMethodsEllipticOperator

  character(len=*), parameter :: smoother = 'GS'
  character(len=*), parameter  :: cycleType = 'V'
  integer, parameter :: numGrids = 5
  integer, parameter :: numCycles = 1
  integer, parameter :: preSweeps = 1
  integer, parameter :: postSweeps = 1
  integer, parameter :: totalIterattions = 15
  integer, parameter :: orderRestriction = 2
  integer, parameter :: orderInterpolation = 2
  real(kind=8), parameter :: alpha = real(1.05,8)
  real(kind=8), parameter :: leftEndpoint = - real(1,8)
  real(kind=8), parameter :: rightEndpoint = real(1,8)
  real(kind=8), parameter :: PMLThickness = real(0.2,8)
  real(kind=8), parameter :: maximumSpeed = real(1,8)
  real(kind=8), parameter :: reflection = exp( - (PMLThickness**3)*alpha )
  real(kind=8), allocatable :: RHS(:,:,:), trueSolution(:,:,:), &
    augmentedRHS(:,:,:), initialGuess(:,:,:), approximateSolution(:,:,:)
  real(kind=8) :: h, preFactorDamping, norm, epsilon
  integer :: N, i, j, k, loop, count
  type(GridDataPML) :: gridDataMG(numGrids)
  type(structureRHSOperator), allocatable :: structRHS(:,:,:)

  call getMultigridStructure(getParameterEllipticEquationDx, &
    getParameterEllipticEquationDy, getParameterEllipticEquationDz, &
    getParameterEllipticEquationDxx, getParameterEllipticEquationDyy, &
    getParameterEllipticEquationDzz, getWaveNumber, gridDataMG, &
    numGrids, leftEndpoint, rightEndpoint )

  ! N = gridDataMG(1)%gridPointsPerAxis
  ! h = gridDataMG(1)%stepSize
  !
  ! do loop = 1, numGrids - 1
  !
  !   N = gridDataMG(loop)%gridPointsPerAxis
  !   h = gridDataMG(loop)%stepSize
  !
  !   allocate( RHS(0:(N+1),0:(N+1),0:(N+1)), &
  !     trueSolution(0:(N+1),0:(N+1),0:(N+1)), &
  !     augmentedRHS(0:(N+1),0:(N+1),0:(N+1)), &
  !     initialGuess(0:(N+1),0:(N+1),0:(N+1)), &
  !     approximateSolution(0:(N+1),0:(N+1),0:(N+1)))
  !
  !   call getRHSEllipticEquation(getParameterEllipticEquationDx, &
  !     getParameterEllipticEquationDy, getParameterEllipticEquationDz, &
  !     getParameterEllipticEquationDxx, getParameterEllipticEquationDyy, &
  !     getParameterEllipticEquationDzz, getWaveNumber, N, h, &
  !     leftEndpoint, RHS)
  !
  !   call getTrueSolutionEllipticEquation(N, h, leftEndpoint, trueSolution)
  !
  !   allocate( structRHS(N,N,N) )
  !
  !   call getStructureRHSOperator(getParameterEllipticEquationDx, &
  !     getParameterEllipticEquationDy, getParameterEllipticEquationDz, &
  !     getParameterEllipticEquationDxx, getParameterEllipticEquationDyy, &
  !     getParameterEllipticEquationDzz, structRHS, N, h, leftEndpoint)
  !
  !   call MATVEC_RHSOperator(structRHS, N, RHS, augmentedRHS)
  !
  !   initialGuess = real(.5,8)
  !   initialGuess(0,:,:) = real(0,8)
  !   initialGuess(N+1,:,:) = real(0,8)
  !   initialGuess(:,0,:) = real(0,8)
  !   initialGuess(:,N+1,:) = real(0,8)
  !   initialGuess(:,:,0) = real(0,8)
  !   initialGuess(:,:,N+1) = real(0,8)
  !
  !   do count = 1, 5000
  !
  !     call SMOOTH_PML(gridDataMG(loop)%LHS_OP, smoother, N, preSweeps, &
  !       augmentedRHS, initialGuess)
  !     call discreteL2Norm(h, N, initialGuess - trueSolution, norm)
  !     !call discreteL2Norm(h, N, initialGuess, norm)
  !     print *, loop, count, 'discreteL2Norm', norm
  !
  !   end do
  !
  !   deallocate( RHS, trueSolution, augmentedRHS, initialGuess, approximateSolution)
  !
  !   deallocate( structRHS )
  !
  !end do

  N = gridDataMG(1)%gridPointsPerAxis
  h = gridDataMG(1)%stepSize

  allocate( RHS(0:(N+1),0:(N+1),0:(N+1)), &
      trueSolution(0:(N+1),0:(N+1),0:(N+1)), &
      augmentedRHS(0:(N+1),0:(N+1),0:(N+1)), &
      initialGuess(0:(N+1),0:(N+1),0:(N+1)), &
      approximateSolution(0:(N+1),0:(N+1),0:(N+1)))

  call getRHSEllipticEquation(getParameterEllipticEquationDx, &
    getParameterEllipticEquationDy, getParameterEllipticEquationDz, &
    getParameterEllipticEquationDxx, getParameterEllipticEquationDyy, &
    getParameterEllipticEquationDzz, getWaveNumber, N, h, leftEndpoint, RHS)

  call getTrueSolutionEllipticEquation(N, h, leftEndpoint, trueSolution)

  allocate( structRHS(N,N,N) )

  call getStructureRHSOperator(getParameterEllipticEquationDx, &
    getParameterEllipticEquationDy, getParameterEllipticEquationDz, &
    getParameterEllipticEquationDxx, getParameterEllipticEquationDyy, &
    getParameterEllipticEquationDzz, structRHS, N, h, leftEndpoint)

  call MATVEC_RHSOperator(structRHS, N, RHS, augmentedRHS)

  initialGuess = real(0,8)
  initialGuess(0,:,:) = real(0,8)
  initialGuess(N+1,:,:) = real(0,8)
  initialGuess(:,0,:) = real(0,8)
  initialGuess(:,N+1,:) = real(0,8)
  initialGuess(:,:,0) = real(0,8)
  initialGuess(:,:,N+1) = real(0,8)

  !initialGuess = trueSolution
  !augmentedRHS = real(0,8)

  do k = 1, N
    do j = 1, N
      do i = 1, N
        print *, '=============================================================='
        print *, i, j, k, N
        print *, 'LHS: center', gridDataMG(1)%LHS_OP(i,j,k)%center
        print *, 'LHS: Side Side X(1)', gridDataMG(1)%LHS_OP(i,j,k)%SS_X(1)
        print *, 'LHS: Side Side X(2)', gridDataMG(1)%LHS_OP(i,j,k)%SS_X(2)
        print *, 'LHS: Side Side Y(1)', gridDataMG(1)%LHS_OP(i,j,k)%SS_Y(1)
        print *, 'LHS: Side Side Y(2)', gridDataMG(1)%LHS_OP(i,j,k)%SS_Y(2)
        print *, 'LHS: Side Side Z(1)', gridDataMG(1)%LHS_OP(i,j,k)%SS_Z(1)
        print *, 'LHS: Side Side Z(2)', gridDataMG(1)%LHS_OP(i,j,k)%SS_Z(2)
        print *, 'LHS: Side Corner XY(1)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XY(1)
        print *, 'LHS: Side Corner XY(2)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XY(2)
        print *, 'LHS: Side Corner XY(3)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XY(3)
        print *, 'LHS: Side Corner XY(4)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XY(4)
        print *, 'LHS: Side Corner XZ(1)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XZ(1)
        print *, 'LHS: Side Corner XZ(2)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XZ(2)
        print *, 'LHS: Side Corner XZ(3)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XZ(3)
        print *, 'LHS: Side Corner XZ(4)', gridDataMG(1)%LHS_OP(i,j,k)%SC_XZ(4)
        print *, 'LHS: Side Corner YZ(1)', gridDataMG(1)%LHS_OP(i,j,k)%SC_YZ(1)
        print *, 'LHS: Side Corner YZ(2)', gridDataMG(1)%LHS_OP(i,j,k)%SC_YZ(2)
        print *, 'LHS: Side Corner YZ(3)', gridDataMG(1)%LHS_OP(i,j,k)%SC_YZ(3)
        print *, 'LHS: Side Corner YZ(4)', gridDataMG(1)%LHS_OP(i,j,k)%SC_YZ(4)
        print *, 'RHS: Center', structRHS(i,j,k)%center/(h*h)
        print *, 'RHS: Side Side X(1)', structRHS(i,j,k)%SS_X(1)/(h*h)
        print *, 'RHS: Side Side X(2)', structRHS(i,j,k)%SS_X(2)/(h*h)
        print *, 'RHS: Side Side Y(1)', structRHS(i,j,k)%SS_Y(1)/(h*h)
        print *, 'RHS: Side Side Y(2)', structRHS(i,j,k)%SS_Y(2)/(h*h)
        print *, 'RHS: Side Side Z(1)', structRHS(i,j,k)%SS_Z(1)/(h*h)
        print *, 'RHS: Side Side Z(2)', structRHS(i,j,k)%SS_Z(2)/(h*h)
        print *, 'RHS: Side Corner', structRHS(i,j,k)%SC/(h*h)
      end do
    end do
  end do

  ! do loop = 1, 1000
  !
  !   call SMOOTH_PML(gridDataMG(1)%LHS_OP, smoother, N, preSweeps, &
  !     augmentedRHS, initialGuess)
  !   call discreteL2Norm(h, N, initialGuess - trueSolution, norm)
  !   !call discreteL2Norm(h, N, initialGuess, norm)
  !   print *, 'discreteL2Norm', norm
  !
  ! end do

  ! do loop = 1, totalIterattions
  !
  !   !initialGuess = trueSolution
  !
  !   print *, '========================', loop ,'================================'
  !   call MGCYC(gridDataMG, smoother, cycleType, N, orderRestriction, &
  !     orderInterpolation, numCycles, preSweeps, postSweeps, &
  !     numGrids, initialGuess, augmentedRHS, approximateSolution)
  !
  !   call discreteL2Norm(h, N, approximateSolution - trueSolution, norm)
  !
  !   ! do k = 0, N + 1
  !   !  do j = 0, N + 1
  !   !    do i = 0, N + 1
  !   !     ! print *, loop, i, j, k, abs( initialGuess(i,j,k) - trueSolution(i,j,k) )
  !   !      print *, loop, i, j, k, abs( approximateSolution(i,j,k) )
  !   !    end do
  !   !  end do
  !   ! end do
  !   !call discreteL2Norm(h, N, approximateSolution - trueSolution, norm)
  !   print *, norm
  !
  !   initialGuess = approximateSolution
  !
  ! end do

end program TestEllipticEquationSolver
