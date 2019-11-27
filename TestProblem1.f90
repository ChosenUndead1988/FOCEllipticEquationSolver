module TestProblem1

  real(kind=8), parameter :: pi = real(4,8)*atan(real(1,8))
  real(kind=8), parameter :: ReynoldsNumber = real(100,8)

  contains

  !=========================================================================
  !   The general elliptic partial differential equation
  !   a(x,y,z) p_{xx} + b(x,y,z)p_{yy} + c(x,y,z) p_{zz}
  !  + d(x,y,z) p_{x} + e(x,y,z) p_{y} + f(x,y,z) p_{z}
  !  + \kappa^{2}(x,y,z) p = f (x,y,z)
  !=========================================================================

  ! write explicit formulas for the elliptic PDE (w/ the exception of the waenumber)
  subroutine evaluateParametersEPDETest1(x, y, z, a, b, c, d, e, f)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8), intent(out) :: a, b, c, d, e, f

    a = real(2,8) + cos( real(2,8)*pi*x )
    b = real(2,8) + cos( real(2,8)*pi*y )
    c = real(2,8) + cos( real(2,8)*pi*z )

    d = - sin(pi*x)/real(2,8)
    e = - sin(pi*y)/real(2,8)
    f = - sin(pi*z)/real(2,8)

    !d = real(0,8)
    !e = real(0,8)
    !f = real(0,8)

  end subroutine evaluateParametersEPDETest1

  ! write explicit formula for the wavenumber squared
  subroutine evaluatewaveNumberSquaredTest1(x, y, z, kappaSq)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8), intent(out) :: kappaSq

    kappaSq = - real(1,8)
  end subroutine evaluatewaveNumberSquaredTest1

  ! the desired analytic solution. Make sure this solution has zero boundary
  ! conditions
  subroutine evaluateTrueSolutionTest1(x, y, z, analyticSolution)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: u, v, w
    real(kind=8), intent(out) :: analyticSolution

    u = pi*x; v = pi*y; w = pi*z;
    analyticSolution = sin(u)*sin(v)*sin(w)

  end subroutine evaluateTrueSolutionTest1

  ! Generate the RHS based on the above analytic solution
  subroutine evaluateRHSTest1(x, y, z, RHS)

    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: a, b, c, d, e, f, kappaSq, p, p_x, p_xx, p_y, p_yy, p_z, p_zz
    real(kind=8), intent(out) :: RHS

    u = pi*x; v = pi*y; w = pi*z

    p = sin(u)*sin(v)*sin(w)

    p_x =  pi*cos(u)*sin(v)*sin(w)
    p_y =  pi*sin(u)*cos(v)*sin(w)
    p_z =  pi*sin(u)*sin(v)*cos(w)

    p_xx = - pi*pi*p; p_yy = - pi*pi*p; p_zz = - pi*pi*p;

    call evaluateParametersEPDETest1(x, y, z, a, b, c, d, e, f)
    call evaluatewaveNumberSquaredTest1(x, y, z, kappaSq)

    RHS = a*p_xx + b*p_yy + c*p_zz + d*p_x + e*p_y + f*p_z + kappaSq*p

  end subroutine evaluateRHSTest1

end module TestProblem1
