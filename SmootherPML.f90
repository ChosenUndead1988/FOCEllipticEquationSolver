module SmootherPML

  use MultigridDataStructures

  contains

  ! Applies either Gauss-Seidel of Four-Color Gauss Seidel. Generally speaking
  ! Gauss-Seidel has better smoothing properties than Jacobis method. We include
  ! Four-Color Gauss Seidel for parallelization capablities

  subroutine SMOOTH_PML(structLHS, smoother, N, numSweeps, RHS, gridFunction)
    type(structureLHSOperator), intent(out) :: structLHS(N,N,N)
    character(len=*), intent(in) :: smoother
    integer, intent(in) :: N, numSweeps
    real(kind=8), intent(in) :: RHS(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8), intent(inout) :: gridFunction(0:(N+1),0:(N+1),0:(N+1))

    select case (smoother)
      case ('GS', 'GaussSeidel')
        call GaussSeidelPML(structLHS, numSweeps, N, RHS, gridFunction)
      case ('FC-GS', 'FCGS','fcgs','fc-gs', '4-GS', '4GS','4gs','4-gs', 'MultiColorGaussSeidel')
        call FourColorGaussSeidelPML(structLHS, numSweeps, N, RHS, gridFunction)
      case default
        call GaussSeidelPML(structLHS, numSweeps, N, RHS, gridFunction)
    end select
  end subroutine SMOOTH_PML

!=============================================================================
!                   Gauss-Seidel
!=============================================================================

  ! apply pointwise Gauss Seidel
  subroutine GaussSeidelPML(structLHS, numSweeps, N, RHS, gridFunction)
    type(structureLHSOperator), intent(inout) :: structLHS(N,N,N)
    integer, intent(in) :: N, numSweeps
    integer :: i, j, k, loop
    real(kind=8), intent(in) :: RHS(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8), intent(inout) :: gridFunction(0:(N+1),0:(N+1),0:(N+1))

    do loop = 1, numSweeps

      do k = 1, N
        do j = 1, N
          do i = 1, N

            call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)

          end do
        end do
      end do

    end do

  end subroutine GaussSeidelPML

  ! apply Four Color Gauss-Seidel. In this version we update one color before
  ! moving onto the next. 
  subroutine FourColorGaussSeidelPML(structLHS, numSweeps, N, RHS, gridFunction)
    type(structureLHSOperator), intent(inout) :: structLHS(N,N,N)
    integer, intent(in) :: N, numSweeps
    integer :: i, j, k, loop
    real(kind=8), intent(in) :: RHS(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8), intent(inout) :: gridFunction(0:(N+1),0:(N+1),0:(N+1))

    do loop = 1, numSweeps
      ! RED
      do k = 1, N
        if ( mod(k,2) .EQ. 1) then

          do j = 1, N, 2
            do i = 1, N, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        else

          do j = 2, N - 1, 2
            do i = 2, N - 1, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        end if

      end do

      ! BLACK
      do k = 1, N

        if ( mod(k,2) .EQ. 1) then

          do j = 2, N - 1, 2
            do i = 1, N, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        else

          do j = 1, N, 2
            do i = 2, N - 1, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        end if

      end do

      !GREEN
      do k = 1, N

        if ( mod(k,2) .EQ. 1) then

          do j = 1, N, 2
            do i = 2, N - 1, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        else

          do j = 2, N - 1, 2
            do i = 1, N, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        end if

      end do

      !ORANGE
      do k = 1, N

        if ( mod(k,2) .EQ. 1) then

          do j = 2, N - 1, 2
            do i = 2, N - 1, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        else

          do j = 1, N, 2
            do i = 1, N, 2
              call localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
            end do
          end do

        end if

      end do

    end do

  end subroutine FourColorGaussSeidelPML

  ! compute pointwise gauss seidel
  subroutine localGaussSeidelPML(structLHS, N, i, j, k, RHS, gridFunction)
    type(structureLHSOperator), intent(out) :: structLHS(N,N,N)
    integer, intent(in) :: N, i, j, k
    real(kind=8), intent(in) :: RHS(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8), intent(inout) :: gridFunction(0:(N+1),0:(N+1),0:(N+1))
    real(kind=8) :: UpperTriangular, LowerTriangular

    ! Strictly upper triangular portion of the EPDE
    UpperTriangular = structLHS(i,j,k)%SS_X(1)*gridFunction(i+1,j,k) &
      + structLHS(i,j,k)%SS_Y(1)*gridFunction(i,j+1,k) &
      + structLHS(i,j,k)%SS_Z(1)*gridFunction(i,j,k+1) &
      + structLHS(i,j,k)%SC_XY(1)*gridFunction(i+1,j+1,k) &
      + structLHS(i,j,k)%SC_XY(2)*gridFunction(i-1,j+1,k) &
      + structLHS(i,j,k)%SC_XZ(1)*gridFunction(i+1,j,k+1) &
      + structLHS(i,j,k)%SC_XZ(2)*gridFunction(i-1,j,k+1) &
      + structLHS(i,j,k)%SC_YZ(1)*gridFunction(i,j+1,k+1) &
      + structLHS(i,j,k)%SC_YZ(2)*gridFunction(i,j-1,k+1)

    ! Strictly lower triangular portion of the EPDE
    LowerTriangular = structLHS(i,j,k)%SS_X(2)*gridFunction(i-1,j,k) &
      + structLHS(i,j,k)%SS_Y(2)*gridFunction(i,j-1,k) &
      + structLHS(i,j,k)%SS_Z(2)*gridFunction(i,j,k-1) &
      + structLHS(i,j,k)%SC_XY(3)*gridFunction(i+1,j-1,k) &
      + structLHS(i,j,k)%SC_XY(4)*gridFunction(i-1,j-1,k) &
      + structLHS(i,j,k)%SC_XZ(3)*gridFunction(i+1,j,k-1) &
      + structLHS(i,j,k)%SC_XZ(4)*gridFunction(i-1,j,k-1) &
      + structLHS(i,j,k)%SC_YZ(3)*gridFunction(i,j+1,k-1) &
      + structLHS(i,j,k)%SC_YZ(4)*gridFunction(i,j-1,k-1)

    ! Applies Gauss-Seidel to the EPDE
    gridFunction(i,j,k) = (RHS(i,j,k) - UpperTriangular &
          - LowerTriangular)/structLHS(i,j,k)%center

  end subroutine localGaussSeidelPML

end module SmootherPML
