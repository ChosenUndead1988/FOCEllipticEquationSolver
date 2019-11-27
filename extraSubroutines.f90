module extraSubroutines

  contains

    subroutine evaluateFunctionOnGrid(formulaFunction, N, &
      leftEndpoint, rightEndpoint, trueSolution)
      interface
        subroutine formulaFunction(x, y, z, trueSolution)
          real(kind=8), intent(in) :: x, y, z
          real(kind=8), intent(out) :: trueSolution
        end subroutine formulaFunction
      end interface
      integer, intent(in) :: N
      integer :: i, j, k
      real(kind=8), intent(in) :: leftEndpoint, rightEndpoint
      real(kind=8), intent(out) :: trueSolution(0:(N+1),0:(N+1),0:(N+1))
      real(kind=8) ::  x, y, z, h

      h = (rightEndpoint - leftEndpoint)/real(N + 1,8)

      do k = 0, N + 1
        do j = 0, N + 1
          do i = 0, N + 1

            x = leftEndpoint + real(i,8)*h
            y = leftEndpoint + real(j,8)*h
            z = leftEndpoint + real(k,8)*h

            call formulaFunction(x, y, z, trueSolution(i,j,k))

          end do
        end do
      end do

  end subroutine evaluateFunctionOnGrid

  subroutine discreteL2Norm(stepSize, N, gridFunction, norm)
    real(kind=8), intent(in) :: stepSize, gridFunction(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8), intent(out) :: norm
    integer, intent(in) :: N

    norm = sqrt( (stepSize**3)*SUM( gridFunction*gridFunction )  )

  end subroutine discreteL2Norm

  subroutine InfNorm(stepSize, N, gridFunction, norm)
    real(kind=8), intent(in) :: stepSize, gridFunction(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8), intent(out) :: norm
    real(kind=8) :: max
    integer :: i, j, k

    max = real(0,8)

    do k = 0, N + 1
      do j = 0, N + 1
        do i = 0, N + 1

          if ( gridFunction(i,j,k) .GE. max ) then
            max = gridFunction(i,j,k)
          end if

        end do
      end do
    end do

    norm = max

  end subroutine InfNorm

end module extraSubroutines
