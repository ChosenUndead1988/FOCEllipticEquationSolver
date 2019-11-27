module MatrixVectorProductMultigrid

	use MultigridDataStructures

	contains

!======================================================================
!														Residual
!======================================================================

	! This subroutine computes the residual r = b - Ax
	subroutine RESIDUAL(structLHS, N, RHS, gridFunctionInput, res)
			integer, intent(in) :: N
			real(kind=8), intent(in) :: RHS(0:(N+1),0:(N+1),0:(N+1)), &
				gridFunctionInput(0:(N+1),0:(N+1),0:(N+1))
			real(kind=8), intent(out) :: res(0:(N+1),0:(N+1),0:(N+1))
			real(kind=8) :: MATVEC(0:(N+1),0:(N+1),0:(N+1))
			type(structureLHSOperator), dimension(N,N,N), intent(in) :: structLHS

			call MATVEC_LHSOperator(structLHS, N, gridFunctionInput, MATVEC)
			res = RHS - MATVEC

			! Assming zero homogeneous boundary conditions implies that the residual
			! is exact along the boundary
			res(0,:,:) = real(0,8)
			res(N+1,:,:) = real(0,8)
			res(:,0,:) = real(0,8)
			res(:,N+1,:) = real(0,8)
			res(:,:,0) = real(0,8)
			res(:,:,N+1) = real(0,8)

	end subroutine RESIDUAL

!======================================================================
!         						Discrete LHS Operator
!======================================================================

	! For the EPDE, preform the matrix vector product of the discrete LHS operator.
	subroutine MATVEC_LHSOperator(structLHS, N, gridFunctionInput, &
		gridFunctionOutput)
		integer, intent(in) :: N
		integer :: i, j, k
		real(kind=8), intent(in) :: gridFunctionInput(0:(N+1),0:(N+1),0:(N+1))
		real(kind=8), intent(out) :: gridFunctionOutput(0:(N+1),0:(N+1),0:(N+1))
		type(structureLHSOperator), dimension(N,N,N), intent(in) :: structLHS

		do k = 1, N
			do j = 1, N
				do i = 1, N

					gridFunctionOutput(i,j,k) = &
						structLHS(i,j,k)%center*gridFunctionInput(i,j,k) &
						+ structLHS(i,j,k)%SS_X(1)*gridFunctionInput(i+1,j,k) &
						+ structLHS(i,j,k)%SS_X(2)*gridFunctionInput(i-1,j,k) &
						+ structLHS(i,j,k)%SS_Y(1)*gridFunctionInput(i,j+1,k) &
						+ structLHS(i,j,k)%SS_Y(2)*gridFunctionInput(i,j-1,k) &
						+ structLHS(i,j,k)%SS_Z(1)*gridFunctionInput(i,j,k+1) &
						+ structLHS(i,j,k)%SS_Z(2)*gridFunctionInput(i,j,k-1) &
						+ structLHS(i,j,k)%SC_XY(1)*gridFunctionInput(i+1,j+1,k) &
						+ structLHS(i,j,k)%SC_XY(2)*gridFunctionInput(i-1,j+1,k) &
						+ structLHS(i,j,k)%SC_XY(3)*gridFunctionInput(i+1,j-1,k) &
						+ structLHS(i,j,k)%SC_XY(4)*gridFunctionInput(i-1,j-1,k)	&
						+ structLHS(i,j,k)%SC_XZ(1)*gridFunctionInput(i+1,j,k+1) &
						+ structLHS(i,j,k)%SC_XZ(2)*gridFunctionInput(i-1,j,k+1) &
						+ structLHS(i,j,k)%SC_XZ(3)*gridFunctionInput(i+1,j,k-1) &
						+ structLHS(i,j,k)%SC_XZ(4)*gridFunctionInput(i-1,j,k-1) &
						+ structLHS(i,j,k)%SC_YZ(1)*gridFunctionInput(i,j+1,k+1) &
						+ structLHS(i,j,k)%SC_YZ(2)*gridFunctionInput(i,j-1,k+1) &
						+ structLHS(i,j,k)%SC_YZ(3)*gridFunctionInput(i,j+1,k-1) &
						+ structLHS(i,j,k)%SC_YZ(4)*gridFunctionInput(i,j-1,k-1)

				end do
			end do
		end do

		! Assuming homogeneous boundary conditions
		gridFunctionOutput(0,:,:) = real(0,8)
		gridFunctionOutput(N+1,:,:) = real(0,8)
		gridFunctionOutput(:,0,:) = real(0,8)
		gridFunctionOutput(:,N+1,:) = real(0,8)
		gridFunctionOutput(:,:,0) = real(0,8)
		gridFunctionOutput(:,:,N+1) = real(0,8)

	end subroutine MATVEC_LHSOperator

!======================================================================
!         						Discrete RHS Operator
!======================================================================

	subroutine MATVEC_RHSOperator(structRHS, N, RHS, augmentedRHS)
		integer, intent(in) :: N
		integer :: i, j, k
		real(kind=8), intent(in) :: RHS(0:(N+1),0:(N+1),0:(N+1))
		real(kind=8), intent(out) :: augmentedRHS(0:(N+1),0:(N+1),0:(N+1))
		type(structureRHSOperator), dimension(N,N,N), intent(in) :: structRHS

		do k = 1, N
			do j = 1, N
				do i = 1, N

					augmentedRHS(i,j,k) = &
						structRHS(i,j,k)%center*RHS(i,j,k) &
						+ structRHS(i,j,k)%SS_X(1)*RHS(i+1,j,k) &
						+ structRHS(i,j,k)%SS_X(2)*RHS(i-1,j,k) &
						+ structRHS(i,j,k)%SS_Y(1)*RHS(i,j+1,k) &
						+ structRHS(i,j,k)%SS_Y(2)*RHS(i,j-1,k) &
						+ structRHS(i,j,k)%SS_Z(1)*RHS(i,j,k+1) &
						+ structRHS(i,j,k)%SS_Z(2)*RHS(i,j,k-1) &
						+ structRHS(i,j,k)%SC_XY(1)*RHS(i+1,j+1,k) &
						+ structRHS(i,j,k)%SC_XY(2)*RHS(i-1,j+1,k) &
						+ structRHS(i,j,k)%SC_XY(3)*RHS(i+1,j-1,k) &
						+ structRHS(i,j,k)%SC_XY(4)*RHS(i-1,j-1,k) &
						+ structRHS(i,j,k)%SC_XZ(1)*RHS(i+1,j,k+1) &
						+ structRHS(i,j,k)%SC_XZ(2)*RHS(i-1,j,k+1) &
						+ structRHS(i,j,k)%SC_XZ(3)*RHS(i+1,j,k-1) &
						+ structRHS(i,j,k)%SC_XZ(4)*RHS(i-1,j,k-1) &
						+ structRHS(i,j,k)%SC_YZ(1)*RHS(i,j+1,k+1) &
						+ structRHS(i,j,k)%SC_YZ(2)*RHS(i,j-1,k+1) &
						+ structRHS(i,j,k)%SC_YZ(3)*RHS(i,j+1,k-1) &
						+ structRHS(i,j,k)%SC_YZ(4)*RHS(i,j-1,k-1)

				end do
			end do
		end do

		augmentedRHS(0,:,:) = real(0,8)
		augmentedRHS(N+1,:,:) = real(0,8)
		augmentedRHS(:,0,:) = real(0,8)
		augmentedRHS(:,N+1,:) = real(0,8)
		augmentedRHS(:,:,0) = real(0,8)
		augmentedRHS(:,:,N+1) = real(0,8)

	end subroutine MATVEC_RHSOperator

end module MatrixVectorProductMultigrid
