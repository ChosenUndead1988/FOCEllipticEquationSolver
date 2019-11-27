module MultigridDataStructures

    type structureLHSOperator

      real(kind=8) :: center
      real(kind=8), dimension(2) :: SS_X, SS_Y, SS_Z
      real(kind=8), dimension(4) :: SC_XY, SC_XZ, SC_YZ

    end type structureLHSOperator

    type structureRHSOperator

      real(kind=8) :: center
      real(kind=8), dimension(2) :: SS_X, SS_Y, SS_Z
      real(kind=8), dimension(4) :: SC_XY, SC_XZ, SC_YZ

    end type structureRHSOperator

    type GridDataPML

      integer :: gridPointsPerAxis
      real(kind=8) :: stepSize
      real(kind=8), dimension(:,:,:), allocatable :: RHS, gridFunction
      type(structureLHSOperator), dimension(:,:,:), allocatable :: LHS_OP

    end type GridDataPML

end module MultigridDataStructures
