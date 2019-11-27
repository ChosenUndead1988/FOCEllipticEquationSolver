module preProcessingEPDE

  use MultigridDataStructures

  contains

  subroutine getMultigridStructure(evaluateParametersEPDE, &
    evaluateWaveNumberSquared, gridDataMG, numGrids, leftEndpoint, &
    rightEndpoint)
    interface
      subroutine evaluateParametersEPDE(x, y, z, a, b, c, d, e, f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE

      subroutine evaluateWaveNumberSquared(x,y,z, kappaSq)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: kappaSq
      end subroutine evaluateWaveNumberSquared
    end interface

      integer, intent(in) :: numGrids
      integer :: N, i, j, k, loop, count
      type(GridDataPML), intent(inout) :: gridDataMG(numGrids)
      real(kind=8), intent(in) :: leftEndpoint, rightEndpoint
      real(kind=8) :: h

      do loop = numGrids, 1, - 1

        count = numGrids + 1 - loop
        N = 2**loop - 1
        h = (rightEndpoint - leftEndpoint)/real(N + 1,8)
        gridDataMG(count)%gridPointsPerAxis = N
        gridDataMG(count)%stepSize = h

        allocate( gridDataMG(count)%RHS(0:(N+1),0:(N+1),0:(N+1)), &
            gridDataMG(count)%gridFunction(0:(N+1),0:(N+1),0:(N+1)), &
            gridDataMG(count)%LHS_OP(N,N,N))

        call getLHSOperator(evaluateParametersEPDE, &
          evaluateWaveNumberSquared, gridDataMG(count)%LHS_OP, &
          N, h, leftEndpoint)

      end do

  end subroutine getMultigridStructure

  ! Compute the coefficients of the LHS operator and store the entries
  ! in a data structure
  subroutine getLHSOperator(evaluateParametersEPDE, &
    evaluateWaveNumberSquared, structLHS, N, stepSize, leftEndpoint)
    interface
      subroutine evaluateParametersEPDE(x, y, z, a, b, c, d, e, f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE

      subroutine evaluateWaveNumberSquared(x,y,z, kappaSq)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: kappaSq
      end subroutine evaluateWaveNumberSquared
    end interface

    integer, intent(in) :: N
    integer :: i, j, k
    type(structureLHSOperator), intent(out) :: structLHS(N,N,N)
    real(kind=8), intent(in) :: stepSize, leftEndpoint
    real(kind=8) :: x, y, z, weightsLHS(0:18)

    do k = 1, N
      do j = 1, N
        do i = 1, N

          x = leftEndpoint + real(i,8)*stepSize
          y = leftEndpoint + real(j,8)*stepSize
          z = leftEndpoint + real(k,8)*stepSize

          call evaluateWeightsLHS_OperatorEPDE(evaluateParametersEPDE, &
            evaluateWaveNumberSquared, x, y, z, stepSize, weightsLHS)

          structLHS(i,j,k)%center = weightsLHS(0)

          structLHS(i,j,k)%SS_X(1) = weightsLHS(1)
          structLHS(i,j,k)%SS_X(2) = weightsLHS(3)
          structLHS(i,j,k)%SS_Y(1) = weightsLHS(2)
          structLHS(i,j,k)%SS_Y(2) = weightsLHS(4)
          structLHS(i,j,k)%SS_Z(1) = weightsLHS(5)
          structLHS(i,j,k)%SS_Z(2) = weightsLHS(6)

          structLHS(i,j,k)%SC_XY(1) = weightsLHS(7)
          structLHS(i,j,k)%SC_XY(2) = weightsLHS(8)
          structLHS(i,j,k)%SC_XY(3) = weightsLHS(10)
          structLHS(i,j,k)%SC_XY(4) = weightsLHS(9)

          structLHS(i,j,k)%SC_XZ(1) = weightsLHS(11)
          structLHS(i,j,k)%SC_XZ(2) = weightsLHS(13)
          structLHS(i,j,k)%SC_XZ(3) = weightsLHS(15)
          structLHS(i,j,k)%SC_XZ(4) = weightsLHS(17)

          structLHS(i,j,k)%SC_YZ(1) = weightsLHS(12)
          structLHS(i,j,k)%SC_YZ(2) = weightsLHS(14)
          structLHS(i,j,k)%SC_YZ(3) = weightsLHS(16)
          structLHS(i,j,k)%SC_YZ(4) = weightsLHS(18)

        end do
      end do
    end do

  end subroutine getLHSOperator

  ! Compute the coefficients of the RHS operator and store the entries
  ! in a data structure
  subroutine getRHSOperator(evaluateParametersEPDE, structRHS, N, stepSize, &
    leftEndpoint)
    interface
      subroutine evaluateParametersEPDE(x, y, z, a, b, c, d, e, f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE
    end interface

    integer, intent(in) :: N
    integer :: i, j, k
    type(structureRHSOperator), intent(out) :: structRHS(N,N,N)
    real(kind=8), intent(in) :: stepSize, leftEndpoint
    real(kind=8) :: x, y, z, weightsRHS(0:18)

    do k = 1, N
      do j = 1, N
        do i = 1, N

          x = leftEndpoint + real(i,8)*stepSize
          y = leftEndpoint + real(j,8)*stepSize
          z = leftEndpoint + real(k,8)*stepSize

          call evaluateWeightsRHS_OperatorEPDE(evaluateParametersEPDE, &
            x, y, z, stepSize, weightsRHS)

          structRHS(i,j,k)%center = weightsRHS(0)

          structRHS(i,j,k)%SS_X(1) = weightsRHS(1)
          structRHS(i,j,k)%SS_X(2) = weightsRHS(3)
          structRHS(i,j,k)%SS_Y(1) = weightsRHS(2)
          structRHS(i,j,k)%SS_Y(2) = weightsRHS(4)
          structRHS(i,j,k)%SS_Z(1) = weightsRHS(5)
          structRHS(i,j,k)%SS_Z(2) = weightsRHS(6)

          structRHS(i,j,k)%SC_XY(1) = weightsRHS(7)
          structRHS(i,j,k)%SC_XY(2) = weightsRHS(8)
          structRHS(i,j,k)%SC_XY(3) = weightsRHS(10)
          structRHS(i,j,k)%SC_XY(4) = weightsRHS(9)

          structRHS(i,j,k)%SC_XZ(1) = weightsRHS(11)
          structRHS(i,j,k)%SC_XZ(2) = weightsRHS(13)
          structRHS(i,j,k)%SC_XZ(3) = weightsRHS(15)
          structRHS(i,j,k)%SC_XZ(4) = weightsRHS(17)

          structRHS(i,j,k)%SC_YZ(1) = weightsRHS(12)
          structRHS(i,j,k)%SC_YZ(2) = weightsRHS(14)
          structRHS(i,j,k)%SC_YZ(3) = weightsRHS(16)
          structRHS(i,j,k)%SC_YZ(4) = weightsRHS(18)

        end do
      end do
    end do

  end subroutine getRHSOperator

  ! manually computes the coefficients of the LHS operator
  subroutine evaluateWeightsLHS_OperatorEPDE(evaluateParametersEPDE, &
    evaluateWaveNumberSquared, x, y, z, stepSize, weightsLHS)
    interface
      subroutine evaluateParametersEPDE(x,y,z,a,b,c,d,e,f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE

      subroutine evaluateWaveNumberSquared(x,y,z, kappaSq)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: kappaSq
      end subroutine evaluateWaveNumberSquared
    end interface

    real(kind=8), dimension(6), parameter :: ONES = real(1,8)
    real(kind=8), parameter :: ONE = real(1,8), TWO = real(2,8), FOUR = real(4,8)
    real(kind=8), intent(in) :: x, y, z, stepSize
    real(kind=8), intent(out) :: weightsLHS(0:18)
    real(kind=8) :: a(0:18), b(0:18), c(0:18), d(0:18), e(0:18), f(0:18), &
      s(0:18), weights(0:18), h, a_ss, b_ss, c_ss, d_ss, e_ss, f_ss, s_ss

    h = stepSize

    call evaluateLHSParametersLocally(evaluateParametersEPDE, &
      evaluateWaveNumberSquared, x, y, z, h, a, b, c, d, e, f, s)

    a_ss = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)
    b_ss = b(1) + b(2) + b(3) + b(4) + b(5) + b(6)
    c_ss = c(1) + c(2) + c(3) + c(4) + c(5) + c(6)
    d_ss = d(1) + d(2) + d(3) + d(4) + d(5) + d(6)
    e_ss = e(1) + e(2) + e(3) + e(4) + e(5) + e(6)
    f_ss = f(1) + f(2) + f(3) + f(4) + f(5) + f(6)
    s_ss = s(1) + s(2) + s(3) + s(4) + s(5) + s(6)

    weights(0) = - TWO*( FOUR*( a(0) + b(0) + c(0) ) &
      + TWO*( a_ss + b_ss + c_ss ) &
      - (a(1) - a(3))*( a(1) + b(1) + c(1) - a(3) - b(3) - c(3))/a(0) &
      - (b(2) - b(4))*( a(2) + b(2) + c(2) - a(4) - b(4) - c(4))/b(0) &
      - (c(5) - c(6))*( a(5) + b(5) + c(5) - a(6) - b(6) - c(6))/c(0) )  &
    - TWO*h*( TWO*( d(1) - d(3) + e(2) - e(4) + f(5) - f(6)) &
      - d(0)*( a(1) - a(3) - b(1) + b(3) - c(1) + c(3))/a(0) &
      + e(0)*( a(2) - a(4) - b(2) + b(4) + c(2) - c(4))/b(0) &
      + f(0)*( a(5) - a(6) + b(5) - b(6) - c(5) + c(6))/c(0) )  &
      + h*h*( TWO*s_ss - (a(1) - a(3))*(s(1) - s(3))/a(0) &
      - (b(2) - b(4))*(s(2) - s(4))/b(0) - (c(5) - c(6))*(s(5) - s(6))/c(0) &
      - FOUR*(d(0)*d(0)/a(0) + e(0)*e(0)/b(0) + f(0)*f(0)/c(0)) ) &
      + TWO*h*h*h*( d(0)*( s(1) - s(3))/a(0) + e(0)*(s(2) - s(4))/b(0) &
       + f(0)*( s(5) - s(6))/c(0) )

    weights(1) = FOUR*( a(0) - b(0) - c(0)) &
      + TWO*( a_ss - c(1) + c(3) - b(1) + b(3)) &
      - (a(1) - a(3))*( a(1) - a(3))/a(0) &
      - (a(2) - a(4))*(b(2) - b(4))/b(0)  &
      - (a(5) - a(6))*(c(5) - c(6))/c(0) &
      + TWO*( b(0) + c(0))*(a(1) - a(3))/a(0) &
    + h*( d_ss  + real(2,8)*( d(0) + d(1) - d(3)) &
        - (a(1) - a(3))*(d(1) - d(3))/(TWO*a(0)) &
        - (b(2) - b(4))*(d(2) - d(4))/(TWO*b(0)) &
        - (c(5) - c(6))*(d(5) - d(6))/(TWO*c(0)) &
        - TWO*d(0)*( b(0) + c(0))/a(0) &
        - d(0)*(a(1) - a(3))/a(0) &
        + e(0)*(a(2) - a(4))/b(0) &
        + f(0)*( a(5) - a(6))/c(0) ) &
    + h*h*( TWO*s(0) + s(1) - s(3) &
      + d(0)*(d(1) - d(3))/(TWO*a(0)) &
      + e(0)*( d(2) - d(4))/(TWO*b(0)) &
      + f(0)*( d(5) - d(6))/(TWO*c(0)) &
      - s(0)*( a(1) - a(3))/a(0) + TWO*d(0)*d(0)/a(0)) &
    + h*h*h*d(0)*s(0)/a(0)

    weights(2) = - FOUR*( a(0) - b(0) + c(0)) &
        + TWO*( b_ss - a(2) + a(4) - c(2) + c(4)) &
        - (a(1) - a(3))*(b(1) - b(3))/a(0) &
        - (b(2) - b(4))*(b(2) - b(4))/b(0) &
        - (c(5) - c(6))*(b(5) - b(6))/c(0) &
        + TWO*(a(0) + c(0))*( b(2) - b(4))/b(0) &
      + h*( e_ss + TWO*( e(0) + e(2) - e(4)) &
        - (a(1) - a(3))*(e(1) - e(3))/(TWO*a(0)) &
        - (b(2) - b(4))*(e(2) - e(4))/(TWO*b(0)) &
        - (c(5) - c(6))*(e(5) - e(6))/(TWO*c(0)) &
        - TWO*e(0)*(a(0) + c(0))/b(0) &
        + d(0)*(b(1) - b(3))/a(0) &
        - e(0)*(b(2) - b(4))/b(0) &
        + f(0)*(b(5) - b(6))/c(0)) &
      + h*h*( TWO*s(0) + s(2) - s(4) &
        + d(0)*(e(1) - e(3))/(TWO*a(0)) &
        + e(0)*(e(2) - e(4))/(TWO*b(0)) &
        + f(0)*(e(5) - e(6))/(TWO*c(0)) &
        - s(0)*(b(2) -b(4))/b(0) + TWO*e(0)*e(0)/b(0)) &
      + h*h*h*e(0)*s(0)/b(0)

    weights(3) = FOUR*( a(0) - b(0) - c(0)) &
      + TWO*( a_ss + c(1) - c(3) + b(1) - b(3)) &
      - (a(1) - a(3))*( a(1) - a(3))/a(0) &
      - (a(2) - a(4))*( b(2) - b(4))/b(0) &
      - (a(5) - a(6))*(c(5) - c(6))/c(0) &
      - TWO*( b(0) + c(0))*(a(1) - a(3))/a(0) &
      - h*( d_ss + TWO*( d(0) - d(1) + d(3)) &
        - (a(1) - a(3))*(d(1) - d(3))/(TWO*a(0)) &
        - (b(2) - b(4))*(d(2) - d(4))/(TWO*b(0)) &
        - (c(5) - c(6))*(d(5) - d(6))/(TWO*c(0)) &
        - TWO*d(0)*( b(0) + c(0))/a(0) &
        + d(0)*(a(1) - a(3))/a(0) &
        - e(0)*( a(2) - a(4)) &
        - f(0)*( a(5) - a(6)) ) &
      + h*h*( TWO*s(0) - s(1) + s(3) &
        - d(0)*(d(1) - d(3))/(TWO*a(0)) &
        - e(0)*( d(2) - d(4))/(TWO*b(0)) &
        - f(0)*( d(5) - d(6))/(TWO*c(0)) &
        + s(0)*( a(1) - a(3))/a(0) + TWO*d(0)*d(0)/a(0)) &
        - h*h*h*d(0)*s(0)/a(0)

    weights(4) = - FOUR*( a(0) - b(0) + c(0)) &
          + TWO*( b_ss + a(2) - a(4) + c(2) - c(4)) &
          - (a(1) - a(3))*(b(1) - b(3))/a(0) &
          - (b(2) - b(4))*(b(2) - b(4))/b(0) &
          - (c(5) - c(6))*(b(5) - b(6))/c(0) &
          - TWO*(a(0) + c(0))*( b(2) - b(4))/b(0) &
        - h*( e_ss + TWO*( e(0) - e(2) + e(4)) &
          - (a(1) - a(3))*(e(1) - e(3))/(TWO*a(0)) &
          - (b(2) - b(4))*(e(2) - e(4))/(TWO*b(0)) &
          - (c(5) - c(6))*(e(5) - e(6))/(TWO*c(0)) &
          - TWO*e(0)*(a(0) + c(0))/b(0) &
          - d(0)*(b(1) - b(3))/a(0) &
          + e(0)*(b(2) - b(4))/b(0) &
          - f(0)*(b(5) - b(6))/c(0)) &
        + h*h*( TWO*s(0) - s(2) + s(4) &
          - d(0)*(e(1) - e(3))/(TWO*a(0)) &
          - e(0)*(e(2) - e(4))/(TWO*b(0)) &
          - f(0)*(e(5) - e(6))/(TWO*c(0)) &
          + s(0)*(b(2) - b(4))/b(0) + TWO*e(0)*e(0)/b(0)) &
        - h*h*h*e(0)*s(0)/b(0)

    weights(5) = - FOUR*( a(0) + b(0) - c(0)) &
        + TWO*( c_ss - a(5) + a(6) - b(5) + b(6)) &
        - (a(1) - a(3))*(c(1) - c(3))/a(0) &
        - (b(2) - b(4))*(c(2) - c(4))/b(0) &
        - (c(5) - c(6))*(c(5) - c(6))/c(0) &
        + TWO*(a(0) + b(0))*( c(5) - c(6))/c(0) &
      + h*( f_ss + TWO*( f(0) + f(5) - f(6)) &
        - (a(1) - a(3))*(f(1) - f(3))/(TWO*a(0)) &
        - (b(2) - b(4))*(f(2) - f(4))/(TWO*b(0)) &
        - (c(5) - c(6))*(f(5) - f(6))/(TWO*c(0)) &
        - TWO*f(0)*(a(0) + b(0))/c(0) &
        + d(0)*( c(1) - c(3))/a(0) &
        + e(0)*( c(2) - c(4))/b(0) &
        - f(0)*( c(5) - c(6))/c(0) ) &
      + h*h*( TWO*s(0) + s(5) - s(6) &
        + d(0)*( f(1) - f(3))/(TWO*a(0)) &
        + e(0)*( f(2) - f(4))/(TWO*b(0)) &
        + f(0)*( f(5) - f(6))/(TWO*c(0)) &
        - s(0)*( c(5) - c(6))/c(0) + TWO*f(0)*f(0)/c(0)) &
      + h*h*h*f(0)*s(0)/c(0)

    weights(6) = - FOUR*( a(0) + b(0) - c(0)) &
          + TWO*( c_ss + a(5) - a(6) + b(5) - b(6)) &
          - (a(1) - a(3))*(c(1) - c(3))/a(0) &
          - (b(2) - b(4))*(c(2) - c(4))/b(0) &
          - (c(5) - c(6))*(c(5) - c(6))/c(0) &
          - TWO*(a(0) + b(0))*( c(5) - c(6))/c(0) &
        - h*( f_ss + TWO*( f(0) - f(5) + f(6)) &
          - (a(1) - a(3))*(f(1) - f(3))/(TWO*a(0)) &
          - (b(2) - b(4))*(f(2) - f(4))/(TWO*b(0)) &
          - (c(5) - c(6))*(f(5) - f(6))/(TWO*c(0)) &
          - TWO*f(0)*(a(0) + b(0))/c(0) &
          - d(0)*( c(1) - c(3))/a(0) &
          - e(0)*( c(2) - c(4))/b(0) &
          + f(0)*( c(5) - c(6))/c(0) ) &
        + h*h*( TWO*s(0) - s(5) + s(6) &
          - d(0)*( f(1) - f(3))/(TWO*a(0)) &
          - e(0)*( f(2) - f(4))/(TWO*b(0)) &
          - f(0)*( f(5) - f(6))/(TWO*c(0)) &
          + s(0)*( c(5) - c(6))/c(0) + TWO*f(0)*f(0)/c(0)) &
        - h*h*h*f(0)*s(0)/c(0)

    weights(7) = TWO*(a(0) + b(0)) + b(1) - b(3) + a(2) - a(4) &
        - b(0)*(a(1) - a(3))/a(0) - a(0)*(b(2) - b(4))/b(0) &
      + (h/TWO)*( TWO*(d(0) + e(0)) + e(1) - e(3) + d(2) - d(4) &
        - e(0)*(a(1) - a(3))/a(0) - d(0)*(b(2) - b(4))/b(0) &
        + TWO*( b(0)*d(0)/a(0) + a(0)*e(0)/b(0)) ) &
      + (h*h/TWO)*d(0)*e(0)*(ONE/a(0) + ONE/b(0))

    weights(8) = TWO*(a(0) + b(0)) - b(1) + b(3) + a(2) - a(4) &
         + b(0)*(a(1) - a(3))/a(0) - a(0)*(b(2) - b(4))/b(0) &
       + (h/TWO)*( TWO*(e(0) - d(0)) - e(1) + e(3) - d(2) + d(4) &
         + e(0)*(a(1) - a(3))/a(0) + d(0)*(b(2) - b(4))/b(0) &
         - TWO*( b(0)*d(0)/a(0) - a(0)*e(0)/b(0)) ) &
       - (h*h/TWO)*d(0)*e(0)*(ONE/a(0) + ONE/b(0))

    weights(9) = TWO*(a(0) + b(0)) - b(1) + b(3) - a(2) + a(4) &
        + b(0)*(a(1) - a(3))/a(0) + a(0)*(b(2) - b(4))/b(0) &
      - (h/TWO)*( TWO*(e(0) + d(0)) - e(1) + e(3) - d(2) + d(4) &
        + e(0)*(a(1) - a(3))/a(0) + d(0)*(b(2) - b(4))/b(0) &
          + TWO*( b(0)*d(0)/a(0) + a(0)*e(0)/b(0)) ) &
      + (h*h/TWO)*d(0)*e(0)*(ONE/a(0) + ONE/b(0))

    weights(10) = TWO*(a(0) + b(0)) + b(1) - b(3) - a(2) + a(4) &
        - b(0)*(a(1) - a(3))/a(0) + a(0)*(b(2) - b(4))/b(0) &
     + (h/TWO)*( TWO*(d(0) - e(0)) - e(1) + e(3) - d(2) + d(4) &
        + e(0)*(a(1) - a(3))/a(0) + d(0)*(b(2) - b(4))/b(0) &
        + TWO*( b(0)*d(0)/a(0) - a(0)*e(0)/b(0)) ) &
     - (h*h/TWO)*d(0)*e(0)*(ONE/a(0) + ONE/b(0) )

    weights(11) = TWO*(a(0) + c(0)) + c(1) - c(3) + a(5) - a(6) &
         - c(0)*(a(1) - a(3))/a(0) - a(0)*(c(5) - c(6))/c(0) &
      + (h/TWO)*( TWO*(f(0) + d(0)) + f(1) - f(3) + d(5) - d(6) &
          - f(0)*( a(1) - a(3))/a(0) - d(0)*(c(5) - c(6))/c(0) &
          + TWO*( c(0)*d(0)/a(0) + a(0)*f(0)/c(0)) )  &
      + (h*h/TWO)*( d(0)*f(0)*(ONE/a(0) + ONE/c(0)))

    weights(13) = TWO*(a(0) + c(0)) - c(1) + c(3) + a(5) - a(6) &
           + c(0)*(a(1) - a(3))/a(0) - a(0)*(c(5) - c(6))/c(0) &
        + (h/TWO)*( TWO*(f(0) - d(0)) - f(1) + f(3) - d(5) + d(6) &
            + f(0)*( a(1) - a(3))/a(0) + d(0)*(c(5) - c(6))/c(0) &
            - TWO*( c(0)*d(0)/a(0) - a(0)*f(0)/c(0)) )  &
        - (h*h/TWO)*( d(0)*f(0)*(ONE/a(0) + ONE/c(0)))

    weights(15) = TWO*(a(0) + c(0)) + c(1) - c(3) - a(5) + a(6) &
          - c(0)*(a(1) - a(3))/a(0) + a(0)*(c(5) - c(6))/c(0) &
      + (h/TWO)*( TWO*(d(0) - f(0)) - f(1) + f(3) - d(5) + d(6) &
           + f(0)*( a(1) - a(3))/a(0) + d(0)*(c(5) - c(6))/c(0) &
           + TWO*( c(0)*d(0)/a(0) - a(0)*f(0)/c(0)) )  &
      - (h*h/TWO)*( d(0)*f(0)*(ONE/a(0) + ONE/c(0)))

    weights(17) = TWO*(a(0) + c(0)) - c(1) + c(3) - a(5) + a(6) &
          + c(0)*(a(1) - a(3))/a(0) + a(0)*(c(5) - c(6))/c(0) &
      - (h/TWO)*( TWO*(d(0) + f(0)) - f(1) + f(3) - d(5) + d(6) &
             + f(0)*( a(1) - a(3))/a(0) + d(0)*(c(5) - c(6))/c(0) &
             + TWO*( c(0)*d(0)/a(0) + a(0)*f(0)/c(0)) )  &
      + (h*h/TWO)*( d(0)*f(0)*(ONE/a(0) + ONE/c(0)))

    weights(12) = TWO*(b(0) + c(0)) + c(2) - c(4) + b(5) - b(6) &
          - c(0)*(b(2) - b(4))/b(0) - b(0)*(c(5) - c(6))/c(0) &
      + (h/TWO)*( TWO*(e(0) + f(0)) + f(2) - f(4) + e(5) - e(6) &
            - f(0)*(b(2) - b(4))/b(0) - e(0)*(c(5) - c(6))/c(0) &
               + TWO*( c(0)*e(0)/b(0) + b(0)*f(0)/c(0)) )  &
      + (h*h/TWO)*e(0)*f(0)*(ONE/b(0) + ONE/c(0))

    weights(14) = TWO*(b(0) + c(0)) - c(2) + c(4) + b(5) - b(6) &
            + c(0)*(b(2) - b(4))/b(0) - b(0)*(c(5) - c(6))/c(0) &
        - (h/TWO)*( TWO*(e(0) - f(0)) + f(2) - f(4) + e(5) - e(6) &
              - f(0)*(b(2) - b(4))/b(0) - e(0)*(c(5) - c(6))/c(0) &
              + TWO*( c(0)*e(0)/b(0) - b(0)*f(0)/c(0)) )  &
        - (h*h/TWO)*e(0)*f(0)*(ONE/b(0) + ONE/c(0))

    weights(16) = TWO*(b(0) + c(0)) + c(2) - c(4) - b(5) + b(6) &
                - c(0)*(b(2) - b(4))/b(0) + b(0)*(c(5) - c(6))/c(0) &
        + (h/TWO)*( TWO*(e(0) - f(0)) - f(2) + f(4) - e(5) + e(6) &
              + f(0)*(b(2) - b(4))/b(0) + e(0)*(c(5) - c(6))/c(0) &
              + TWO*( c(0)*e(0)/b(0) - b(0)*f(0)/c(0)) )  &
        - (h*h/TWO)*e(0)*f(0)*(ONE/b(0) + ONE/c(0))

    weights(18) = TWO*(b(0) + c(0)) - c(2) + c(4) - b(5) + b(6) &
            + c(0)*(b(2) - b(4))/b(0) + b(0)*(c(5) - c(6))/c(0) &
       - (h/TWO)*( TWO*(e(0) + f(0)) - f(2) + f(4) - e(5) + e(6) &
              + f(0)*(b(2) - b(4))/b(0) + e(0)*(c(5) - c(6))/c(0) &
              + TWO*( c(0)*e(0)/b(0) + b(0)*f(0)/c(0)) )  &
       + (h*h/TWO)*e(0)*f(0)*(ONE/b(0) + ONE/c(0))

    weightsLHS = weights/real(24,8)

  end subroutine evaluateWeightsLHS_OperatorEPDE

  ! manually computes the coefficients of the RHS operator
  subroutine evaluateWeightsRHS_OperatorEPDE(evaluateParametersEPDE, &
    x, y, z, stepSize, weightsRHS)
    interface
      subroutine evaluateParametersEPDE(x,y,z,a,b,c,d,e,f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE

    end interface

    real(kind=8), intent(in) :: x, y, z, stepSize
    real(kind=8), intent(out) :: weightsRHS(0:18)
    real(kind=8) :: a(0:18), b(0:18), c(0:18), d(0:18), e(0:18), f(0:18), &
      s(0:18), weights(0:18), h, extraTerms(0:18)

    h = stepSize

    call evaluateRHSParametersLocally(evaluateParametersEPDE, &
      x, y, z, h, a, b, c, d, e, f)

    weights(0) = h*h*real(12,8)

    weights(1) = h*h*real(2,8) - h*h*(a(1) - a(3))/a(0) + h*h*h*d(0)/a(0)
    weights(2) = h*h*real(2,8) - h*h*(b(2) - b(4))/b(0) + h*h*h*e(0)/b(0)
    weights(3) = h*h*real(2,8)  + h*h*(a(1) - a(3))/a(0) - h*h*h*d(0)/a(0)
    weights(4) = h*h*real(2,8)  + h*h*(b(2) - b(4))/b(0) - h*h*h*e(0)/b(0)
    weights(5) = h*h*real(2,8)  - h*h*(c(5) - c(6))/c(0) + h*h*h*f(0)/c(0)
    weights(6) = h*h*real(2,8)  + h*h*(c(5) - c(6))/c(0) - h*h*h*f(0)/c(0)

    weights(7:18) = real(0,8)

    call getParametersExtraFOCTerm(a, b, c, ExtraTerms)

    weightsRHS = weights/real(24,8) + h*h*ExtraTerms

  end subroutine evaluateWeightsRHS_OperatorEPDE

  ! computes the coefficients of the ellipic partial differential equation
  ! a p_xx + b p_yy + c p_zz + d p_x + e p_y + f p_z + kappaSq p = g
  subroutine evaluateLHSParametersLocally(evaluateParametersEPDE, &
    evaluateWaveNumberSquared, x, y, z, stepSize, A, B, C, D, &
    E, F, kappaSqVec)
    interface
      subroutine evaluateParametersEPDE(x,y,z,a,b,c,d,e,f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE

      subroutine evaluateWaveNumberSquared(x,y,z, kappaSq)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: kappaSq
      end subroutine evaluateWaveNumberSquared
    end interface

    real(kind=8), intent(in) :: x, y, z, stepSize
    real(kind=8), intent(out) :: A(0:18), B(0:18), C(0:18), &
      D(0:18), E(0:18), F(0:18), kappaSqVec(0:18)
    real(kind=8) :: h, u, v, w, coordinate(0:18,1:3)
    integer :: i

    h = stepSize

    coordinate(0,1:3) = (/ x, y, z /)
    coordinate(1,1:3) = (/ x + h, y, z /)
    coordinate(2,1:3) = (/ x, y + h, z /)
    coordinate(3,1:3) = (/ x - h, y, z /)
    coordinate(4,1:3) = (/ x, y - h, z /)
    coordinate(5,1:3) = (/ x, y, z + h /)
    coordinate(6,1:3) = (/ x, y, z - h /)
    coordinate(7,1:3) = (/ x + h, y + h, z /)
    coordinate(8,1:3) = (/ x - h, y + h, z /)
    coordinate(9,1:3) = (/ x - h, y - h, z /)
    coordinate(10,1:3) = (/ x + h, y - h, z /)
    coordinate(11,1:3) = (/ x + h, y, z + h /)
    coordinate(12,1:3) = (/ x, y + h, z + h /)
    coordinate(13,1:3) = (/ x - h, y, z + h /)
    coordinate(14,1:3) = (/ x, y - h, z + h /)
    coordinate(15,1:3) = (/ x + h, y, z - h /)
    coordinate(16,1:3) = (/ x, y + h, z - h /)
    coordinate(17,1:3) = (/ x - h, y, z - h /)
    coordinate(18,1:3) = (/ x, y - h, z - h /)

    do i = 0, 18
      u = coordinate(i,1); v = coordinate(i,2); w = coordinate(i,3);
      call evaluateParametersEPDE(u, v, w, A(i),B(i),C(i),D(i),E(i),F(i))
      call evaluateWaveNumberSquared(u, v, w, kappaSqVec(i))
    end do

  end subroutine evaluateLHSParametersLocally

  subroutine evaluateRHSParametersLocally(evaluateParametersEPDE, &
    x, y, z, stepSize, A, B, C, D, E, F)
    interface
      subroutine evaluateParametersEPDE(x, y, z, a, b, c, d, e, f)
        real(kind=8), intent(in) :: x, y, z
        real(kind=8), intent(out) :: a, b, c, d, e, f
      end subroutine evaluateParametersEPDE

    end interface

    real(kind=8), intent(in) :: x, y, z, stepSize
    real(kind=8), intent(out) :: A(0:18), B(0:18), C(0:18), D(0:18), &
      E(0:18), F(0:18)
    real(kind=8) :: h, u, v, w, coordinate(0:18,1:3)
    integer :: i

    h = stepSize

    coordinate(0,1:3) = (/ x, y, z /)
    coordinate(1,1:3) = (/ x + h, y, z /)
    coordinate(2,1:3) = (/ x, y + h, z /)
    coordinate(3,1:3) = (/ x - h, y, z /)
    coordinate(4,1:3) = (/ x, y - h, z /)
    coordinate(5,1:3) = (/ x, y, z + h /)
    coordinate(6,1:3) = (/ x, y, z - h /)
    coordinate(7,1:3) = (/ x + h, y + h, z /)
    coordinate(8,1:3) = (/ x - h, y + h, z /)
    coordinate(9,1:3) = (/ x - h, y - h, z /)
    coordinate(10,1:3) = (/ x + h, y - h, z /)
    coordinate(11,1:3) = (/ x + h, y, z + h /)
    coordinate(12,1:3) = (/ x, y + h, z + h /)
    coordinate(13,1:3) = (/ x - h, y, z + h /)
    coordinate(14,1:3) = (/ x, y - h, z + h /)
    coordinate(15,1:3) = (/ x + h, y, z - h /)
    coordinate(16,1:3) = (/ x, y + h, z - h /)
    coordinate(17,1:3) = (/ x - h, y, z - h /)
    coordinate(18,1:3) = (/ x, y - h, z - h /)

    do i = 0, 18
      u = coordinate(i,1); v = coordinate(i,2); w = coordinate(i,3);
      call evaluateParametersEPDE(u, v, w, A(i),B(i),C(i),D(i),E(i),F(i))

    end do

  end subroutine evaluateRHSParametersLocally

  ! These terms increase the diagonal dominance of the RHS operator
  subroutine getParametersExtraFOCTerm(a, b, c, ExtraTerms)
    real(kind=8), intent(in) :: a(0:18), b(0:18), c(0:18)
    real(kind=8), intent(out) :: ExtraTerms(0:18)
    real(kind=8) :: phi(0:18), varphi(0:18), psi(0:18), Exy(0:18), &
      Exz(0:18), Eyz(0:18)

      phi = ( real(1,8)/a + real(1,8)/b )
      varphi = ( real(1,8)/a + real(1,8)/c )
      psi = ( real(1,8)/b + real(1,8)/c )

      Exy = real(0,8)

      Exy(0) = (a(0) + b(0))*(  real(16,8)*phi(0)  &
       - real(4,8)*( phi(1) + phi(2) + phi(3) + phi(4) )  &
       + phi(7) + phi(8) + phi(9) + phi(10) )

      Exy(1) = (a(0) + b(0))*(- real(4,8)*phi(0) + phi(2) + phi(4) &
        - real(2,8)*( phi(1) - phi(3) ) &
        + (phi(7) + phi(10) - phi(8) - phi(9))/real(2,8) )

      Exy(3) = (a(0) + b(0))*( - real(4,8)*phi(0) + phi(2) + phi(4) &
          + real(2,8)*( phi(1) - phi(3) ) &
          - (phi(7) + phi(10) - phi(8) - phi(9))/real(2,8) )

      Exy(2) = (a(0) + b(0))*( - real(4,8)*phi(0) + phi(1) + phi(3) &
          - real(2,8)*( phi(2) - phi(4) ) &
          + (phi(7) + phi(8) - phi(9) - phi(10))/real(2,8) )

      Exy(4) = (a(0) + b(0))*( - real(4,8)*phi(0) + phi(1) + phi(3) &
          + real(2,8)*( phi(2) - phi(4) ) &
          - (phi(7) + phi(8) - phi(9) - phi(10))/real(2,8) )

      Exy(7) = (a(0) + b(0))*( phi(0) &
        + (phi(1) - phi(3) + phi(2) - phi(4))/real(2,8) &
        + (phi(7) + phi(9) - phi(8) - phi(10))/real(4,8) )

      Exy(8) = (a(0) + b(0))*( phi(0) &
          + (- phi(1) + phi(3) + phi(2) - phi(4))/real(2,8) &
          - (phi(7) + phi(9) - phi(8) - phi(10))/real(4,8) )

      Exy(9) = (a(0) + b(0))*( phi(0) &
          + (- phi(1) + phi(3) - phi(2) + phi(4))/real(2,8) &
          + (phi(7) + phi(9) - phi(8) - phi(10))/real(4,8) )

      Exy(10) = (a(0) + b(0))*( phi(0) &
          + (phi(1) - phi(3) - phi(2) + phi(4))/real(2,8) &
          - (phi(7) + phi(9) - phi(8) - phi(10))/real(4,8) )

      Exz = real(0,8)

      Exz(0) = (a(0) + c(0))*(  real(16,8)*varphi(0)  &
          - real(4,8)*( varphi(1) + varphi(3) + varphi(5) + varphi(6) )  &
          + varphi(11) + varphi(13) + varphi(15) + varphi(17) )

      Exz(1) = (a(0) + c(0))*( - real(4,8)*varphi(0)  &
          + varphi(5) + varphi(6) &
          - real(2,8)*( varphi(1) - varphi(3) )  &
          + (varphi(11) - varphi(13) + varphi(15) - varphi(17))/real(2,8) )

      Exz(3) = (a(0) + c(0))*( - real(4,8)*varphi(0)  &
          + varphi(5) + varphi(6) &
          + real(2,8)*( varphi(1) - varphi(3) )  &
          - (varphi(11) - varphi(13) + varphi(15) - varphi(17) )/real(2,8) )

      Exz(5) = (a(0) + c(0))*( - real(4,8)*varphi(0)  &
          + varphi(1) + varphi(3) &
          - real(2,8)*( varphi(5) - varphi(6) )  &
          + (varphi(11) + varphi(13) - varphi(15) - varphi(17) )/real(2,8) )

      Exz(6) = (a(0) + c(0))*( - real(4,8)*varphi(0)  &
          + varphi(1) + varphi(3) &
          + real(2,8)*( varphi(5) - varphi(6) )  &
          - (varphi(11) + varphi(13) - varphi(15) - varphi(17) )/real(2,8) )

      Exz(11) = (a(0) + c(0))*( varphi(0) &
        + (varphi(1) - varphi(3) + varphi(5) - varphi(6))/real(2,8) &
        + (varphi(11) - varphi(13) - varphi(15) + varphi(17))/real(4,8))

      Exz(13) = (a(0) + c(0))*( varphi(0) &
          + (- varphi(1) + varphi(3) + varphi(5) - varphi(6))/real(2,8) &
          - (varphi(11) - varphi(13) - varphi(15) + varphi(17))/real(4,8))

      Exz(15) = (a(0) + c(0))*( varphi(0) &
          + (varphi(1) - varphi(3) - varphi(5) + varphi(6))/real(2,8) &
          - (varphi(11) - varphi(13) - varphi(15) + varphi(17))/real(4,8))

      Exz(17) = (a(0) + c(0))*( varphi(0) &
          + (- varphi(1) + varphi(3) - varphi(5) + varphi(6))/real(2,8) &
          + (varphi(11) - varphi(13) - varphi(15) + varphi(17))/real(4,8))

      Eyz = real(0,8)

      Eyz(0) = (b(0) + c(0))*(  real(16,8)*psi(0)  &
          - real(4,8)*( psi(2) + psi(4) + psi(5) + psi(6) )  &
          + psi(12) + psi(14) + psi(16) + psi(18) )

      Eyz(2) = (b(0) + c(0))*( - real(4,8)*psi(0) &
          + psi(5) + psi(6) - real(2,8)*( psi(2) - psi(4) ) &
          + (psi(12) - psi(14) + psi(16) - psi(18))/real(2,8) )

      Eyz(4) = (b(0) + c(0))*( - real(4,8)*psi(0) &
          + psi(5) + psi(6) + real(2,8)*( psi(2) - psi(4) ) &
          - (psi(12) - psi(14) + psi(16) - psi(18))/real(2,8) )

      Eyz(5) = (b(0) + c(0))*( - real(4,8)*psi(0) &
          + psi(2) + psi(4) - real(2,8)*( psi(5) - psi(6) ) &
          + (psi(12) + psi(14) - psi(16) - psi(18))/real(2,8) )

      Eyz(6) = (b(0) + c(0))*( - real(4,8)*psi(0) &
          + psi(2) + psi(4) + real(2,8)*( psi(5) - psi(6) ) &
          - (psi(12) + psi(14) - psi(16) - psi(18))/real(2,8) )

      Eyz(12) = (b(0) + c(0))*( psi(0) &
          + (psi(2) - psi(4) + psi(5) - psi(6))/real(2,8) &
          + (psi(12) - psi(14) - psi(16) + psi(18))/real(4,8) )

      Eyz(14) = (b(0) + c(0))*( psi(0) &
          + (- psi(2) + psi(4) + psi(5) - psi(6))/real(2,8) &
          - (psi(12) - psi(14) - psi(16) + psi(18))/real(4,8) )

      Eyz(16) = (b(0) + c(0))*( psi(0) &
          + (psi(2) - psi(4) - psi(5) + psi(6))/real(2,8) &
          - (psi(12) - psi(14) - psi(16) + psi(18))/real(4,8) )

      Eyz(18) = (b(0) + c(0))*( psi(0) &
          + (- psi(2) + psi(4) - psi(5) + psi(6))/real(2,8) &
          + (psi(12) - psi(14) - psi(16) + psi(18))/real(4,8) )

      ExtraTerms = (Exy + Exz + Eyz)/real(144,8)
  end subroutine getParametersExtraFOCTerm

end module preProcessingEPDE
