module MultigridGridTransferOperators

  Implicit REAL*8(A-H,O-Z)

  real(kind=8), parameter :: scaling_restriction = real(4,8)
  contains

!==============================================================================
!                   Restriction Operator
!==============================================================================

  subroutine RESTRICTION_OPERATOR(order, N_coarse, N_fine, &
    gridFunctionCoarse, gridFunctionFine)
    integer, intent(in) :: order, N_coarse, N_fine
    real(kind=8), intent(in) :: gridFunctionFine(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1))
    real(kind=8):: gridFunctionCopy(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1))
    real(kind=8), dimension(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1)),  &
      intent(out) :: gridFunctionCoarse

      if (order .EQ. 0) then
        ! Injection - i.e. zeroeth order restriction operator
        call Injection(N_coarse, N_fine, gridFunctionCopy, gridFunctionFine)
      else if (order .EQ. 2) then
        ! Full Weighting - i.e. second order restriction operator
        call FullWeighting(N_coarse, N_fine, gridFunctionCopy, &
          gridFunctionFine)
      else
        ! use full weighting as the default restriction operator
        call FullWeighting(N_coarse, N_fine, gridFunctionCopy, &
          gridFunctionFine)
      end if

      ! Since we "clear the denominantor" to form the discrete linear system
      ! we must rescale the restriction operator
      gridFunctionCoarse = scaling_restriction*gridFunctionCopy

  end subroutine RESTRICTION_OPERATOR

  subroutine FullWeighting(N_coarse, N_fine, gridFunctionCoarse, &
    gridFunctionFine)
    integer, intent(in) :: N_coarse, N_fine
    integer :: i, j, k, r, s, t, UBc, UBf
    real(kind=8), intent(in) :: &
      gridFunctionFine(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1))
    real(kind=8), dimension(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1)),  &
      intent(out) :: gridFunctionCoarse
    real(kind=8) :: center, sideSide, sideCorner, cornerCorner
    real(kind=8), parameter :: c0 = real(1.0/8.0,8), cSS = real(1.0/16.0,8), &
      cSC = real(1.0/32.0,8), cCC = real(1.0/64.0,8)

    ! gridFunctionCoarse(0,:,:) = real(0,8)
    ! gridFunctionCoarse(N_coarse + 1,:,:) = real(0,8)
    ! gridFunctionCoarse(:,0,:) = real(0,8)
    ! gridFunctionCoarse(:,N_coarse + 1,:) = real(0,8)
    ! gridFunctionCoarse(:,:,0) = real(0,8)
    ! gridFunctionCoarse(:,:,N_coarse + 1) = real(0,8)

    UBc = N_coarse + 1; UBf = N_fine + 1;
    gridFunctionCoarse(0,:,:) = gridFunctionFine(0,0:UBf:2,0:UBf:2)
    gridFunctionCoarse(UBc,:,:) = gridFunctionFine(UBf,0:UBf:2,0:UBf:2)
    gridFunctionCoarse(:,0,:) = gridFunctionFine(0:UBf:2,0,0:UBf:2)
    gridFunctionCoarse(:,UBc,:) = gridFunctionFine(0:UBf:2,UBf,0:UBf:2)
    gridFunctionCoarse(:,:,0) = gridFunctionFine(0:UBf:2,0:UBf:2,0)
    gridFunctionCoarse(:,:,UBc) = gridFunctionFine(0:UBf:2,0:UBf:2,UBf)

    ! perform full-weighting operator on the interior points of the fine grid

    do k = 1, N_coarse
      do j = 1, N_coarse
        do i = 1, N_coarse
          r = 2*i; s = 2*j; t = 2*k;
          center = gridFunctionFine(r,s,t)

          sideSide = gridFunctionFine(r-1,s,t) &
            + gridFunctionFine(r+1,s,t) + gridFunctionFine(r,s-1,t) &
            + gridFunctionFine(r,s+1,t) + gridFunctionFine(r,s,t-1) &
            + gridFunctionFine(r,s,t+1)

          sideCorner = gridFunctionFine(r+1,s+1,t) &
            + gridFunctionFine(r-1,s+1,t) + gridFunctionFine(r+1,s-1,t) &
            + gridFunctionFine(r-1,s-1,t) + gridFunctionFine(r,s+1,t-1) &
            + gridFunctionFine(r,s+1,t+1) + gridFunctionFine(r,s-1,t-1) &
            + gridFunctionFine(r,s-1,t+1) + gridFunctionFine(r+1,s,t+1) &
            + gridFunctionFine(r+1,s,t-1) + gridFunctionFine(r-1,s,t-1) &
            + gridFunctionFine(r-1,s,t+1)

          cornerCorner = gridFunctionFine(r+1,s+1,t-1) &
            + gridFunctionFine(r-1,s+1,t-1) + gridFunctionFine(r+1,s-1,t-1) &
            + gridFunctionFine(r-1,s-1,t-1) + gridFunctionFine(r+1,s+1,t+1) &
            + gridFunctionFine(r-1,s+1,t+1) + gridFunctionFine(r+1,s-1,t+1) &
            + gridFunctionFine(r-1,s-1,t+1)

          gridFunctionCoarse(i,j,k) = c0*center &
            + cSS*sideSide + cSC*sideCorner + cCC*cornerCorner
        end do
      end do
    end do

  end subroutine FullWeighting

  subroutine Injection(N_coarse, N_fine, gridFunctionCoarse, &
    gridFunctionFine)
    integer, intent(in) :: N_coarse, N_fine
    integer :: Ubc, UBf
    real(kind=8), dimension(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1)),  &
      intent(in) :: gridFunctionFine
    real(kind=8), dimension(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1)),  &
      intent(out) :: gridFunctionCoarse

      ! gridFunctionCoarse(0,:,:) = real(0,8)
      ! gridFunctionCoarse(N_coarse + 1,:,:) = real(0,8)
      ! gridFunctionCoarse(:,0,:) = real(0,8)
      ! gridFunctionCoarse(:,N_coarse + 1,:) = real(0,8)
      ! gridFunctionCoarse(:,:,0) = real(0,8)
      ! gridFunctionCoarse(:,:,N_coarse + 1) = real(0,8)

      !gridFunctionCoarse(1:N_coarse, 1:N_coarse,1:N_coarse) = &
      !gridFunctionFine(2:(N_fine-1):2, 2:(N_fine-1):2, 2:(N_fine-1):2)

      Ubc = N_coarse + 1; Ubf = N_fine + 1
      gridFunctionCoarse(:,:,:) = gridFunctionFine(0:UBf:2, 0:UBf:2, 0:UBf:2)

   end subroutine Injection

!===========================================================================
!                   Prolongation Operators
!===========================================================================

   subroutine PROLONGATION_OPERATOR(order, N_coarse, N_fine, &
     gridFunctionCoarse, interpolationFine)
     integer, intent(in) :: order, N_coarse, N_fine
     real(kind=8), dimension(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1)),  &
       intent(inout) :: interpolationFine
     real(kind=8), dimension(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1)),  &
       intent(in) :: gridFunctionCoarse

       if (order .EQ. 2) then
         call TriLinearInterpolation(N_coarse, N_fine, gridFunctionCoarse, &
           interpolationFine)
       else if (order .EQ. 3 ) then
        call TriQuadraticInterpolation(N_coarse, N_fine, gridFunctionCoarse, &
          interpolationFine)
       else if (order .EQ. 4) then
         call TriCubicInterpolation(N_coarse, N_fine, gridFunctionCoarse, &
           interpolationFine)
       else
         call TriLinearInterpolation(N_coarse, N_fine, gridFunctionCoarse, &
           interpolationFine)
       end if

   end subroutine PROLONGATION_OPERATOR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             Tri Cubic Interpolation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Uses third order interpolating polynomial along each direction
  subroutine TriCubicInterpolation(N_coarse, N_fine, &
    gridFunctionCoarse, interpolationFine)
   integer, intent(in) :: N_coarse, N_fine
   integer :: i, j, k, u, v, w, flag_x, flag_y, flag_z, &
    loop_x, loop_y, loop_z, X(4), Y(4), Z(4)
   real(kind=8), intent(out) :: &
    interpolationFine(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1))
   real(kind=8), intent(in) :: &
    gridFunctionCoarse(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1))
   real(kind=8) :: lagrange(4,3)

   interpolationFine(0:(N_fine+1):2,0:(N_fine+1):2, 0:(N_fine+1):2) = &
    gridFunctionCoarse(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1))

    ! coefficients of the lagrange interpolating polynomial over the canonical
    ! interval [-1, 1]
    lagrange(1,1) = real(5,8)/real(16,8)
    lagrange(2,1) = real(15,8)/real(16,8)
    lagrange(3,1) = - real(5,8)/real(16,8)
    lagrange(4,1) =  real(1,8)/real(16,8)

    lagrange(1,2) = - real(1,8)/real(16,8)
    lagrange(2,2) = real(9,8)/real(16,8)
    lagrange(3,2) = real(9,8)/real(16,8)
    lagrange(4,2) = - real(1,8)/real(16,8)

    lagrange(1,3) = real(1,8)/real(16,8)
    lagrange(2,3) = - real(5,8)/real(16,8)
    lagrange(3,3) = real(15,8)/real(16,8)
    lagrange(4,3) = real(5,8)/real(16,8)

    ! Interpolation along the x direction
    do u = 1, N_fine, 2
      do w = 0, N_fine + 1, 2
        do v = 0, N_fine + 1, 2

          i = (u - 1)/2; j = v/2; k = w/2

          if (i .LE. N_coarse - 2 ) then
            flag_x = 1
            X = (/i, i + 1, i + 2, i + 3/)
          else if (i .EQ. N_coarse - 1) then
            flag_x = 2
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_x = 3
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_x = 1, 4
           sum = sum + gridFunctionCoarse(X(loop_x),j,k)*lagrange(loop_x, flag_x)
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

    ! Interpolation along the y direction
    do v = 1, N_fine, 2
      do w = 0, N_fine + 1, 2
        do u = 0, N_fine + 1, 2

          i = u/2; j = (v - 1)/2; k = w/2

          if (j .LE. N_coarse - 2 ) then
            flag_y = 1
            Y = (/j, j + 1, j + 2, j + 3 /)
          else if (j .EQ. N_coarse - 1) then
            flag_y = 2
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_y = 3
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_y = 1, 4
           sum = sum + gridFunctionCoarse(i,Y(loop_y),k)*lagrange(loop_y, flag_y)
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

    ! Interpolation along the z direction
    do w = 1, N_fine, 2
      do v = 0, N_fine + 1, 2
        do u = 0, N_fine + 1, 2

          i = u/2; j = v/2; k = (w - 1)/2

          if (k .LE. N_coarse - 2 ) then
            flag_z = 1
            Z = (/k, k + 1, k + 2, k + 3/)
          else if (k .EQ. N_coarse - 1) then
            flag_z = 2
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_z = 3
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_z = 1, 4
            sum = sum + gridFunctionCoarse(i,j,Z(loop_z))*lagrange(loop_z, flag_z)
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

    ! Interpolation along the xy direction
    do u = 1, N_fine, 2
      do v = 1, N_fine, 2
        do w = 0, N_fine + 1, 2

          i = (u - 1)/2; j = (v - 1)/2; k = w/2

          if (i .LE. N_coarse - 2 ) then
            flag_x = 1
            X = (/i, i + 1, i + 2, i + 3/)
          else if (i .EQ. N_coarse - 1) then
            flag_x = 2
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_x = 3
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          if (j .LE. N_coarse - 2 ) then
            flag_y = 1
            Y = (/j, j + 1, j + 2, j + 3 /)
          else if (j .EQ. N_coarse - 1) then
            flag_y = 2
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_y = 3
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_y = 1, 4
            do loop_x = 1, 4
              sum = sum + gridFunctionCoarse(X(loop_x),Y(loop_y),k)*( &
                lagrange(loop_x,flag_x)*lagrange(loop_y,flag_y) )
            end do
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

    ! Interpolation along the xz direction
    do u = 1, N_fine, 2
      do w = 1, N_fine, 2
        do v = 0, N_fine + 1, 2

          i = (u - 1)/2; j = v/2; k = (w - 1)/2

          if (i .LE. N_coarse - 2 ) then
            flag_x = 1
            X = (/i, i + 1, i + 2, i + 3/)
          else if (i .EQ. N_coarse - 1) then
            flag_x = 2
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_x = 3
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          if (k .LE. N_coarse - 2 ) then
            flag_z = 1
            Z = (/k, k + 1, k + 2, k + 3/)
          else if (k .EQ. N_coarse - 1) then
            flag_z = 2
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_z = 3
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_z = 1, 4
            do loop_x = 1, 4
              sum = sum + gridFunctionCoarse(X(loop_x),j,Z(loop_z))*( &
                lagrange(loop_x,flag_x)*lagrange(loop_z,flag_z) )
            end do
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

    ! Interpolation along the yz direction
    do w = 1, N_fine, 2
      do v = 1, N_fine, 2
        do u = 0, N_fine + 1, 2

          i = u/2; j = (v-1)/2; k = (w - 1)/2

          if (j .LE. N_coarse - 2 ) then
            flag_y = 1
            Y = (/j, j + 1, j + 2, j + 3 /)
          else if (j .EQ. N_coarse - 1) then
            flag_y = 2
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_y = 3
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          if (k .LE. N_coarse - 2 ) then
            flag_z = 1
            Z = (/k, k + 1, k + 2, k + 3/)
          else if (k .EQ. N_coarse - 1) then
            flag_z = 2
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_z = 3
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_z = 1, 4
            do loop_y = 1, 4
              sum = sum + gridFunctionCoarse(i,Y(loop_y),Z(loop_z))*( &
                lagrange(loop_y,flag_y)*lagrange(loop_z,flag_z) )
            end do
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

    ! Interpolation along the xyz direction
    do w = 1, N_fine, 2
      do v = 1, N_fine, 2
        do u = 1, N_fine , 2

          i = (u - 1)/2; j = (v - 1)/2; k = (w - 1)/2

          if (i .LE. N_coarse - 2 ) then
            flag_x = 1
            X = (/i, i + 1, i + 2, i + 3/)
          else if (i .EQ. N_coarse - 1) then
            flag_x = 2
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_x = 3
            X = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          if (j .LE. N_coarse - 2 ) then
            flag_y = 1
            Y = (/j, j + 1, j + 2, j + 3 /)
          else if (j .EQ. N_coarse - 1) then
            flag_y = 2
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_y = 3
            Y = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          if (k .LE. N_coarse - 2 ) then
            flag_z = 1
            Z = (/k, k + 1, k + 2, k + 3/)
          else if (k .EQ. N_coarse - 1) then
            flag_z = 2
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          else
            flag_z = 3
            Z = (/ N_coarse - 2, N_coarse - 1, N_coarse, N_coarse + 1/)
          end if

          sum = real(0,8)

          do loop_z = 1, 4
            do loop_y = 1, 4
              do loop_x = 1, 4

                sum = sum + gridFunctionCoarse(X(loop_x),Y(loop_y),Z(loop_z))*( &
                  lagrange(loop_x,flag_x))*(lagrange(loop_y,flag_y))*( &
                  lagrange(loop_z,flag_z))
              end do
            end do
          end do

          interpolationFine(u,v,w) = sum

        end do
      end do
    end do

  end subroutine TriCubicInterpolation

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             Tri Quadratic Interpolation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Uses second order interpolating polynomial along each direction
    subroutine TriQuadraticInterpolation(N_coarse, N_fine, &
      gridFunctionCoarse, interpolationFine)
       integer, intent(in) :: N_coarse, N_fine
       integer :: i, j, k, u, v, w, flag_x, flag_y, flag_z, &
        loop_x, loop_y, loop_z, X(3), Y(3), Z(3)
       real(kind=8), intent(out) :: &
        interpolationFine(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1))
       real(kind=8), intent(in) :: &
        gridFunctionCoarse(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1))
        real(kind=8) :: lagrange(3,2)

      interpolationFine(0:(N_fine+1):2,0:(N_fine+1):2, 0:(N_fine+1):2) = &
        gridFunctionCoarse(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1))

      lagrange(1,1) = real(3,8)/real(8,8)
      lagrange(2,1) = real(3,8)/real(4,8)
      lagrange(3,1) = - real(1,8)/real(8,8)

      lagrange(1,2) = - real(1,8)/real(8,8)
      lagrange(2,2) = real(3,8)/real(4,8)
      lagrange(3,2) = real(3,8)/real(8,8)

      ! Interpolation along the x direction
      do u = 1, N_fine, 2
        do w = 0, N_fine + 1, 2
          do v = 0, N_fine + 1, 2

            i = (u - 1)/2; j = v/2; k = w/2

            if (i .LE. N_coarse - 1 ) then
              flag_x = 1
              X = (/i, i + 1, i + 2/)
            else
              flag_x = 2
              X = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_x = 1, 3
              sum = sum + gridFunctionCoarse(X(loop_x),j,k)*lagrange(loop_x, flag_x)
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

      ! Interpolation along the y direction
      do v = 1, N_fine, 2
        do w = 0, N_fine + 1, 2
          do u = 0, N_fine + 1, 2

            i = u/2; j = (v - 1)/2; k = w/2

            if (j .LE. N_coarse - 1 ) then
              flag_y = 1
              Y = (/j, j + 1, j + 2/)
            else
              flag_y = 2
              Y = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_y = 1, 3
              sum = sum + gridFunctionCoarse(i,Y(loop_y),k)*lagrange(loop_y, flag_y)
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

      ! Interpolation along the z direction
      do w = 1, N_fine, 2
        do v = 0, N_fine + 1, 2
          do u = 0, N_fine + 1, 2

            i = u/2; j = v/2; k = (w - 1)/2

            if (k .LE. N_coarse - 1 ) then
              flag_z = 1
              Z = (/k, k + 1, k + 2/)
            else
              flag_z = 2
              Z = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_z = 1, 3
              sum = sum + gridFunctionCoarse(i,j,Z(loop_z))*lagrange(loop_z, flag_z)
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

      ! Interpolation along the xy direction
      do u = 1, N_fine, 2
        do v = 1, N_fine, 2
          do w = 0, N_fine + 1, 2

            i = (u - 1)/2; j = (v - 1)/2; k = w/2

            if (i .LE. N_coarse - 1 .AND. j .LE. N_coarse - 1) then
              flag_x = 1; flag_y = 1;
              X = (/i, i + 1, i + 2/); Y = (/j, j + 1, j + 2/);
            else if (i .LE. N_coarse - 1 .AND. j .GE. N_coarse) then
              flag_x = 1; flag_y = 2;
              X = (/i, i + 1, i + 2/)
              Y = (/N_coarse - 1, N_coarse, N_coarse + 1/)
            else if (i .GE. N_coarse .AND. j .LE. N_coarse - 1 ) then
              flag_x = 2; flag_y = 1;
              X = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
              Y = (/ j, j + 1, j + 2 /)
            else
              flag_x = 2; flag_y = 2;
              X = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
              Y = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_y = 1, 3
              do loop_x = 1, 3
                sum = sum + gridFunctionCoarse(X(loop_x),Y(loop_y),k)*( &
                  lagrange(loop_x,flag_x)*lagrange(loop_y,flag_y) )
              end do
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

      ! Interpolation along the xz direction
      do u = 1, N_fine, 2
        do w = 1, N_fine, 2
          do v = 0, N_fine + 1, 2

            i = (u - 1)/2; j = v/2; k = (w - 1)/2

            if (i .LE. N_coarse - 1 .AND. k .LE. N_coarse - 1) then
              flag_x = 1; flag_z = 1;
              X = (/i, i + 1, i + 2/); Z = (/k, k + 1, k + 2/);
            else if (i .LE. N_coarse - 1 .AND. k .GE. N_coarse) then
              flag_x = 1; flag_z = 2;
              X = (/i, i + 1, i + 2/)
              Z = (/N_coarse - 1, N_coarse, N_coarse + 1/)
            else if (i .GE. N_coarse .AND. k .LE. N_coarse - 1 ) then
              flag_x = 2; flag_z = 1;
              X = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
              Z = (/ k, k + 1, k + 2 /)
            else
              flag_x = 2; flag_z = 2;
              X = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
              Z = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_z = 1, 3
              do loop_x = 1, 3
                sum = sum + gridFunctionCoarse(X(loop_x),j,Z(loop_z))*( &
                  lagrange(loop_x,flag_x)*lagrange(loop_z,flag_z) )
              end do
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

      ! Interpolation along the yz direction
      do w = 1, N_fine, 2
        do v = 1, N_fine, 2
          do u = 0, N_fine + 1, 2

            i = u/2; j = (v-1)/2; k = (w - 1)/2

            if (j .LE. N_coarse - 1 .AND. k .LE. N_coarse - 1) then
              flag_y = 1; flag_z = 1;
              Y = (/j, j + 1, j + 2/); Z = (/k, k + 1, k + 2/);
            else if (j .LE. N_coarse - 1 .AND. k .GE. N_coarse) then
              flag_y = 1; flag_z = 2;
              Y = (/j, j + 1, j + 2/)
              Z = (/N_coarse - 1, N_coarse, N_coarse + 1/)
            else if (j .GE. N_coarse .AND. k .LE. N_coarse - 1 ) then
              flag_y = 2; flag_z = 1;
              Y = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
              Z = (/ k, k + 1, k + 2 /)
            else
              flag_y = 2; flag_z = 2;
              Y = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
              Z = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_z = 1, 3
              do loop_y = 1, 3
                sum = sum + gridFunctionCoarse(i,Y(loop_y),Z(loop_z))*( &
                  lagrange(loop_y,flag_y)*lagrange(loop_z,flag_z) )
              end do
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

      ! Interpolation along the xyz direction
      do w = 1, N_fine, 2
        do v = 1, N_fine, 2
          do u = 1, N_fine , 2

            i = (u - 1)/2; j = (v - 1)/2; k = (w - 1)/2

            if (i .LE. N_coarse - 1 ) then
              flag_x = 1
              X = (/i, i + 1, i + 2/)
            else
              flag_x = 2
              X = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            if (j .LE. N_coarse - 1 ) then
              flag_y = 1
              Y = (/j, j + 1, j + 2/)
            else
              flag_y = 2
              Y = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            if (k .LE. N_coarse - 1 ) then
              flag_z = 1
              Z = (/k, k + 1, k + 2/)
            else
              flag_z = 2
              Z = (/ N_coarse - 1, N_coarse, N_coarse + 1/)
            end if

            sum = real(0,8)

            do loop_z = 1, 3
              do loop_y = 1, 3
                do loop_x = 1, 3

                  sum = sum + gridFunctionCoarse(X(loop_x),Y(loop_y),Z(loop_z))*( &
                    lagrange(loop_x,flag_x))*(lagrange(loop_y,flag_y))*( &
                    lagrange(loop_z,flag_z))
                end do
              end do
            end do

            interpolationFine(u,v,w) = sum

          end do
        end do
      end do

    end subroutine TriQuadraticInterpolation

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             Tri Linear Interpolation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !   Linear Interpolation
    subroutine TriLinearInterpolation(N_coarse, N_fine, gridFunctionCoarse, &
      interpolationFine)
      integer, intent(in) :: N_coarse, N_fine
      integer :: i, j, k, Re, So, Se, To, Te
      real(kind=8), dimension(3) :: coeff
      real(kind=8), dimension(0:(N_fine+1),0:(N_fine+1),0:(N_fine+1)),  &
        intent(out) :: interpolationFine
      real(kind=8), dimension(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1)),  &
        intent(in) :: gridFunctionCoarse

      coeff(1) = real(1.0/2.0,8)
      coeff(2) = real(1.0/4.0,8)
      coeff(3) = real(1.0/8.0,8)

      interpolationFine(0:(N_fine+1):2,0:(N_fine+1):2,0:(N_fine+1):2) = &
       gridFunctionCoarse(0:(N_coarse+1),0:(N_coarse+1),0:(N_coarse+1))

      ! linear interpolation along the x-direction
      do k = 0, N_coarse + 1
        do j = 0, N_coarse + 1
          do i = 0, N_coarse
            interpolationFine(2*i+1,2*j,2*k) = coeff(1)*( &
              gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i+1,j,k))
          end do
        end do
      end do

      ! linear interpolation along the y-direction
      do k = 0, N_coarse + 1
        do i = 0, N_coarse + 1
          do j = 0, N_coarse
            interpolationFine(2*i,2*j+1,2*k) = coeff(1)*( &
              gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i,j+1,k))
          end do
        end do
      end do

      ! linear interpolation along the z-direction
      do j = 0, N_coarse + 1
        do i = 0, N_coarse + 1
          do k = 0, N_coarse
            interpolationFine(2*i,2*j,2*k+1) = coeff(1)*( &
              gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i,j,k+1))
          end do
        end do
      end do

      ! linear interpolation along the xy-direction
      do k = 0, N_coarse + 1
        do j = 0, N_coarse
          do i = 0, N_coarse
            interpolationFine(2*i+1 ,2*j+1 ,2*k ) = coeff(2)*( &
            gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i+1,j,k) &
            + gridFunctionCoarse(i,j+1,k) + gridFunctionCoarse(i+1,j+1,k))
          end do
       end do
      end do

      ! linear interpolation along the xz-direction
      do j = 0, N_coarse + 1
        do k = 0, N_coarse
          do i = 0, N_coarse
            interpolationFine(2*i+1 ,2*j ,2*k + 1) = coeff(2)*( &
            gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i+1,j,k) &
            + gridFunctionCoarse(i,j,k+1) + gridFunctionCoarse(i+1,j,k+1))
          end do
       end do
      end do

      ! linear interpolation along the yz-direction
      do i = 0, N_coarse + 1
        do k = 0, N_coarse
          do j = 0, N_coarse
            interpolationFine(2*i,2*j + 1 ,2*k + 1) = coeff(2)*( &
             gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i,j+1,k) &
             + gridFunctionCoarse(i,j,k+1) + gridFunctionCoarse(i,j+1,k+1))
          end do
       end do
      end do

      ! linear interpolation along the xyz-direction
      do k = 0, N_coarse
        do j = 0, N_coarse
          do i = 0, N_coarse
            interpolationFine(2*i + 1,2*j + 1,2*k + 1) = coeff(3)*( &
              gridFunctionCoarse(i,j,k) + gridFunctionCoarse(i+1,j,k) &
              + gridFunctionCoarse(i,j+1,k) + gridFunctionCoarse(i,j,k+1) &
              + gridFunctionCoarse(i+1,j+1,k) + gridFunctionCoarse(i+1,j,k+1) &
              + gridFunctionCoarse(i,j+1,k+1) + gridFunctionCoarse(i+1,j+1,k+1))

          end do
       end do
     end do

    end subroutine TriLinearInterpolation

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             Test Polynomials (for testing the degree of exactness )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine TestLinearPolynomial(N, h, a, b, gridfunction)
      integer, intent(in) :: N
      integer :: i, j, k
      real(kind=8) :: x, y, z, prod_x, prod_y, prod_z
      real(kind=8), intent(in) :: h, a, b
      real(kind=8), dimension(0:(N+1),0:(N+1),0:(N+1)), intent(out) :: gridFunction


      do k = 0, N + 1
        do j = 0, N + 1
          do i = 0, N + 1
            x = a + real(i,8)*h
            y = a + real(j,8)*h
            z = a + real(k,8)*h

            prod_x = (x - a) + real(1,8)
            prod_y = (y - a) + real(1,8)
            prod_z = (z - a) + real(1,8)

            gridFunction(i,j,k) = prod_x*prod_y*prod_z

          end do
        end do
      end do

    end subroutine TestLinearPolynomial

    subroutine TestQuadraticPolynomial(N, h, a, b, gridfunction)
      integer, intent(in) :: N
      integer :: i, j, k
      real(kind=8) :: x, y, z, prod_x, prod_y, prod_z
      real(kind=8), intent(in) :: h, a, b
      real(kind=8), dimension(0:(N+1),0:(N+1),0:(N+1)), intent(out) :: gridFunction

      do k = 0, N + 1
        do j = 0, N + 1
          do i = 0, N + 1
            x = a + real(i,8)*h
            y = a + real(j,8)*h
            z = a + real(k,8)*h

            prod_x = (x - a)*(x - b) + real(1,8)
            prod_y = (y - a)*(y - b) + real(1,8)
            prod_z = (z - a)*(z - b) + real(1,8)

            gridFunction(i,j,k) = prod_x*prod_y*prod_z

          end do
        end do
      end do

    end subroutine TestQuadraticPolynomial

    subroutine TestCubicPolynomial(N, h, a, b, gridFunction)
      integer, intent(in) :: N
      integer :: i, j, k
      real(kind=8) :: x, y, z, prod_x, prod_y, prod_z
      real(kind=8), intent(in) :: h, a, b
      real(kind=8), dimension(0:(N+1),0:(N+1),0:(N+1)), intent(out) :: gridFunction


      do k = 0, N + 1
        do j = 0, N + 1
          do i = 0, N + 1
            x = a + real(i,8)*h
            y = a + real(j,8)*h
            z = a + real(k,8)*h;

            prod_x = (x-a)*(x-a)*(x-b) + real(2,8)
            prod_y = (y-a)*(y-b)*(y-b) + real(3,8)
            prod_z = (z-a)*(z-b)*(z-b) + real(1,8)

            gridFunction(i,j,k) = prod_x*prod_y*prod_z
          end do
        end do
      end do

    end subroutine TestCubicPolynomial

    subroutine TestQuarticPolynomial(N, h, a, b, gridFunction)
      integer, intent(in) :: N
      integer :: i, j, k
      real(kind=8) :: x, y, z, prod_x, prod_y, prod_z
      real(kind=8), intent(in) :: h, a, b
      real(kind=8), dimension(0:(N+1),0:(N+1),0:(N+1)), intent(out) :: gridFunction

      gridFunction(0,:,:) = real(0,8)
      gridFunction(N+1,:,:) = real(0,8)
      gridFunction(:,0,:) = real(0,8)
      gridFunction(:,N+1,:) = real(0,8)
      gridFunction(:,:,0) = real(0,8)
      gridFunction(:,:,N+1) = real(0,8)

      do k = 1, N
        do j = 1, N
          do i = 1, N
            x = a + real(i,8)*h
            y = a + real(j,8)*h
            z = a + real(k,8)*h;

            prod_x = (x-a)*(x-a)*(x-b)*(x-b)
            prod_y = (y-a)*(y-b)*(y-a)*(y-b)
            prod_z = (z-a)*(z-b)*(z-a)*(z-b)

            gridFunction(i,j,k) = prod_x*prod_y*prod_z
          end do
        end do
      end do

    end subroutine TestQuarticPolynomial

end module MultigridGridTransferOperators
