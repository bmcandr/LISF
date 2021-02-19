!------------------------------------------------------------------------------
!BOP
module LIS_create_gridMod

   use ESMF
   use LIS_logMod,     only : LIS_verify, LIS_logunit 
   use LIS_coreMod !,    only : LIS_masterproc, LIS_localPet
   !use LIS_coreMod,  only : vm=>LIS_vm, LIS_rc, LIS_Domain
   
   implicit none
   private
   public  :: createRectilinearGrid
   public  :: determine_lat_corners, determine_lon_corners
   public  :: getInteriorGrid
   public  :: create_rectilinear_grid
   public  :: create_regular_grid
   public  :: create_regular_grid_2D
!
! !DESCRIPTION:
! This module contains functions for creating ESMF grids.
   
!EOP
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
   function createRectilinearGrid(lon_center_points, lat_center_points, &
                                 lon_corner_points, lat_corner_points, &
                                 grid_name, Nx, Ny, &
                                 coordSys, xgrids, ygrids, using_model, periodic) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: lat_center_points(:)
      real,             intent(in)    :: lon_center_points(:)
      real,             intent(in)    :: lat_corner_points(:)
      real,             intent(in)    :: lon_corner_points(:)
      type(ESMF_CoordSys_Flag)        :: coordSys ! Coordinate system: cartesian, spherical
      integer,          intent(in)    :: xgrids(Nx*Ny) ! number of x grid points in each processor
      integer,          intent(in)    :: ygrids(Nx*Ny) ! number of y grid points in each processor
      logical,          intent(in)    :: using_model   ! are we trying to create the model grid?
      logical,               optional :: periodic
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create a generic ESMF 2D rectilinear grid given arbitrary latitude and
! longitude grid points.
!
! !LOCAL VARIABLES:
      integer                         :: i, j, num_lats, num_lons
      integer                         :: rc, STATUS
      logical                         :: periodic_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = "createRectilinearGrid"
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF Rectilinear Grid: '// TRIM(grid_name)

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      IF (using_model) THEN
         countsPerDEDimX(:) = xgrids(1:Nx)
         countsPerDEDimY(:) = ygrids(1:(Ny-1)*Nx:Nx)
      ELSE
         num_lons = SIZE(lon_center_points)
         call decomposeDim(num_lons, countsPerDEDimX, Nx )

         num_lats = SIZE(lat_center_points)
         call decomposeDim(num_lats, countsPerDEDimY, Ny )
      ENDIF

      ! Create the ESMF Grid
      !---------------------
      periodic_ = .FALSE.
      if (PRESENT(periodic)) periodic_ = periodic

      if (periodic_) then
         new_grid = ESMF_GridCreate1PeriDim( &
                         name            = TRIM(grid_name), &
                         countsPerDEDim1 = countsPerDEDimX, &
                         countsPerDEDim2 = countsPerDEDimY, &
                         indexflag       = ESMF_INDEX_DELOCAL, & !ESMF_INDEX_GLOBAL, &
                         gridEdgeLWidth  = [0,0], &
                         gridEdgeUWidth  = [0,1], &
                         coordDep1       = [1,2], &
                         coordDep2       = [1,2], &
                         coordSys        = coordSys, rc=rc)

         call LIS_verify(rc, 'ESMF_GridCreate1PeriDim failed in '//TRIM(IAm))
      else
        new_grid = ESMF_GridCreateNoPeriDim( &
                         name            = TRIM(grid_name), &
                         countsPerDEDim1 = countsPerDEDimX, &
                         countsPerDEDim2 = countsPerDEDimY, &
                         indexflag       = ESMF_INDEX_DELOCAL, &! ESMF_INDEX_GLOBAL, &
                         gridEdgeLWidth  = [0,0], &
                         gridEdgeUWidth  = [1,1], &
                         coordDep1       = [1,2], &
                         coordDep2       = [1,2], &
                         coordSys        = coordSys, rc=rc)

         call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed in '//TRIM(IAm))
      endif

      ! Allocate coords at default stagger location
      call ESMF_GridAddCoord(new_grid, rc=rc)
      call LIS_verify(rc, 'Gneric ESMF_GridAddCoord failed in '//TRIM(IAm))

      call ESMF_GridAddCoord(new_grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed in '//TRIM(IAm))

      call ESMF_AttributeSet(new_grid, 'GridType', 'LatLon', rc=rc)
      call LIS_verify(rc, 'ESMF_AttributeSet failed in '//TRIM(IAm))

      if (.NOT. periodic_) then
         call ESMF_AttributeSet(new_grid, 'Global', .FALSE., rc=rc)
      endif

      call add_horz_coordinates(new_grid, periodic_, &
                                lat_center_points, lon_center_points, &
                                lat_corner_points, lon_corner_points)

      write(LIS_logunit,*) '[INFO] Done creating the ESMF Rectilinear Grid: '// TRIM(grid_name)

   end function createRectilinearGrid
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine add_horz_coordinates(grid, periodic, &
                                   lat_center_points, lon_center_points, &
                                   lat_corner_points, lon_corner_points, rc)
!
! !INPUT PARAMETERS:
      logical,          intent(in)    :: periodic
      real,             intent(in)    :: lat_center_points(:)
      real,             intent(in)    :: lon_center_points(:)
      real,             intent(in)    :: lat_corner_points(:)
      real,             intent(in)    :: lon_corner_points(:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Grid), intent(inout) :: grid
      integer, optional, intent(out) :: rc
!
! ! LOCAL VARIABLES:
      integer :: i_1, i_n, j_1, j_n ! regional array bounds
      integer :: ic_1,ic_n,jc_1,jc_n ! regional corner bounds
      real(kind=ESMF_KIND_R8), pointer :: centers(:,:)
      real(kind=ESMF_KIND_R8), pointer :: corners(:,:)
      integer :: status, num_lats, num_lons
      integer :: i, j, ij(4)
!EOP
!------------------------------------------------------------------------------
!BOC
      call getInteriorGrid(grid, i_1, i_n, j_1, j_n)
      ij(1)=i_1
      ij(2)=i_n
      ij(3)=j_1
      ij(4)=j_n

      if (.not. any(ij == -1)) then
         if (periodic) then
            ic_1=i_1
            ic_n=i_n
         else
            ic_1=i_1
            if (i_n == num_lons) then
               ic_n=i_n+1
            else
               ic_n=i_n
            end if
         end if

         jc_1=j_1
         if (j_n == num_lats) then
            jc_n=j_n+1
         else
            jc_n=j_n
         end if

         ! First we handle longitudes:
         call ESMF_GridGetCoord(grid, coordDim=1, localDE=0, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                       farrayPtr=centers, rc=status)
         call LIS_verify(status, 'ESMF_GridGetCoord for lon center failed')

         call ESMF_GridGetCoord(grid, coordDim=1, localDE=0, &
                       staggerloc=ESMF_STAGGERLOC_CORNER, &
                       farrayPtr=corners, rc=status)
         call LIS_verify(status, 'ESMF_GridGetCoord for lon corner failed')

         do j = 1, size(centers,2)
            centers(:,j) = lon_center_points(i_1:i_n)
         end do

         do j = 1, size(corners,2)
            corners(:,j) = lon_corner_points(ic_1:ic_n)
         end do

         ! Now latitudes
         call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                       farrayPtr=centers, rc=status)
         call LIS_verify(status, 'ESMF_GridGetCoord for lat center failed')

         call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
                       staggerloc=ESMF_STAGGERLOC_CORNER, &
                       farrayPtr=corners, rc=status)
         call LIS_verify(status, 'ESMF_GridGetCoord for lat corner failed')

         do i = 1, size(centers,1)
            centers(i,:) = lat_center_points(j_1:j_n)
         end do
         do i = 1, size(corners,1)
            corners(i,:) = lat_corner_points(jc_1:jc_n)
         end do
      end if

      end subroutine add_horz_coordinates
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine getInteriorGrid(GRID, I1, IN, J1, JN)
!
! INPUT PARAMETERS:
      type (ESMF_Grid), intent(IN) :: grid
!
! OUTPUT PARAMETERS:
      integer, intent(OUT)         :: I1, IN, J1, JN

! !LOCAL VARIABLES:
      integer                               :: status
      type (ESMF_DistGrid)                  :: distGrid
      type(ESMF_DELayout)                   :: LAYOUT
      integer,               allocatable    :: AL(:,:)
      integer,               allocatable    :: AU(:,:)
      integer                               :: nDEs,localDECount
      integer                               :: deId
      integer                               :: gridRank
      integer,               allocatable    :: localDeToDeMap(:)
      !character(len=ESMF_MAXSTR)            :: IAm='getInteriorGrid'
!EOP
!------------------------------------------------------------------------------
!BOC
      i1=-1
      j1=-1
      in=-1
      jn=-1
      call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=STATUS)
      call LIS_verify(status, 'ESMF_GridGet failed')

      call ESMF_DistGridGet(distGRID, delayout=layout, rc=STATUS)
      call LIS_verify(status, 'ESMF_DistGridGet failed')

      call ESMF_DELayoutGet(layout, deCount = nDEs, localDeCount=localDeCount,rc=status)
      call LIS_verify(status, 'ESMF_DELayoutGet failed')

      if (localDeCount > 0) then
         allocate(localDeToDeMap(localDeCount),stat=status)
         call ESMF_DELayoutGet(layout, localDEtoDeMap=localDeToDeMap,rc=status)
         call LIS_verify(status, 'ESMF_DELayoutGet failed')

         deId=localDeToDeMap(1)

         allocate (AL(gridRank,0:nDEs-1),  stat=status)
         allocate (AU(gridRank,0:nDEs-1),  stat=status)

         call ESMF_DistGridGet(distGrid,  &
                               minIndexPDe = AL, &
                               maxIndexPDe = AU, rc=status)
         call LIS_verify(status, 'ESMF_DistGridGet failed')

         I1 = AL(1, deId)
         IN = AU(1, deId)
         !    _ASSERT(gridRank > 1, 'tilegrid is 1d (without RC this only for info')
         J1 = AL(2, deId)
         JN = AU(2, deId)
         deallocate(AU, AL, localDeToDeMap)
      end if

      end subroutine getInteriorGrid
!EOC
!------------------------------------------------------------------------------
!BOP
      function determine_lon_corners(center_lons, periodic) result(corner_lons)
!
! !INPUT PARAMETERS:
      real,    intent(IN) :: center_lons(:)
      logical, intent(in) :: periodic
!
! !RETURNED VALUE:
      real, allocatable :: corner_lons(:)
!
! !DESCRIPTION:
! Given longitude grid points at center locations,
! this function deteremine the corresponding grid points at corners.
! The center_points array can be:
! \begin{itemize}
! \item Increasing values without periodicity
! \item Increasing values with periodicity
! \item Longitude values (from 0 to 180-dx and from -180+dx to -dx) 
!       from a Gaussian grid
! \end{itemize}
!
! !LOCAL VARIABLES:
      integer :: num_lons, i, m, k
      real    :: dx, dy
      real    :: min_point, max_point
      real, allocatable :: temp_corners(:), temp_centers(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      num_lons = size(center_lons)
      allocate(corner_lons(num_lons+1))
      corner_lons(:) = -999.0

      if ( (center_lons(1) == minval(center_lons)) .AND. &
           (center_lons(num_lons) == maxval(center_lons)) )then
         ! Values are in increasing order
         ! This applies to latitudes and longitudes in lat/lon grid

         ! Longitudes from -180 to 180
         dx = center_lons(2)-center_lons(1)
         corner_lons(  1) = max(-180.0, center_lons(1) - 0.5*dx)
         dx = center_lons(num_lons)-center_lons(num_lons-1)
         corner_lons(num_lons+1) = min( 180.0, center_lons(num_lons) + 0.5*dx)
         do i=2,num_lons
            corner_lons(i) = 0.5*(center_lons(i) + center_lons(i-1))
          enddo
      else
         ! This is for the longitudes in a Gaussian grids
         ! Values are from 0 to 180-dx and from -180+dx to -dx

         ! Determine the first time a negative value is encountered
         do i = 1, num_lons
            if (center_lons(i) < 0) then
               m = i
               EXIT
            endif
         end do

         ALLOCATE(temp_centers(num_lons))
         ALLOCATE(temp_corners(num_lons))

         temp_centers(1:num_lons-m+1) = center_lons(m:num_lons)
         temp_centers(num_lons-m+2:) = center_lons(1:m-1)
         dx = temp_centers(2) - temp_centers(1)
         temp_corners(1) = max(-180.0, temp_centers(1) - 0.5*dx)
         do i = 2, num_lons
            temp_corners(i) = 0.5*(temp_centers(i) + temp_centers(i-1))
         enddo

         do i = 1, num_lons
            if (temp_corners(i) > 0) then
               k = i
               EXIT
            endif
         end do
         corner_lons(1:num_lons-k+1) = temp_corners(k:num_lons)
         corner_lons(num_lons-k+2:)  = temp_corners(1:k-1)

         DEALLOCATE(temp_corners, temp_centers)

!         ! ---> from 0 to 180-dx
!         corner_lons(1) = 0.5*(center_lons(1) + center_lons(num_lons))
!         do i = 1, m-1
!            corner_lons(i) = 0.5*(center_lons(i) + center_lons(i-1))
!         enddo
!
!         ! ---> from -180+dx to -dx
!         corner_lons(m) = 0.5*(-180.0 + center_lons(m))
!         do i = m+1, num_lons
!            corner_lons(i) = 0.5*(center_lons(i-1) + center_lons(i))
!         enddo
!         corner_lons(num_lons+1) = 0.5*(center_lons(num_lons) + center_lons(1))
         !PRINT*, "LONS: ", num_lons,m, center_lons(m-1), center_lons(m),corner_lons(m)
      endif

      end function determine_lon_corners
!EOC
!------------------------------------------------------------------------------
!BOP
      function determine_lat_corners(center_lats, periodic) result(corner_lats)
!
! !INPUT PARAMETERS:
      real,    intent(IN) :: center_lats(:)
      logical, intent(in) :: periodic
!
! !RETURNED VALUE:
      real, allocatable :: corner_lats(:)
!
! !DESCRIPTION:
! Given (latitude or longitude) grid points at center locations,
! this function deteremine the corresponding grid points at corners.
! The center_lats array can be:
! \begin{itemize}
! \item Increasing values without periodicity
! \item Increasing ovalues with periodicity
! \item Latitude values (from 90-dy to -90+dy) from a Gaussian grid
! \end{itemize}
!
! !LOCAL VARIABLES:
      integer :: num_lats, i, m
      real    :: dx, dy
      real    :: min_point, max_point
!EOP
!------------------------------------------------------------------------------
!BOC
      num_lats = size(center_lats)
      allocate(corner_lats(num_lats+1))
      corner_lats(:) = -999.0

      if ((maxval(center_lats) == center_lats(1)) .AND. &
          (minval(center_lats) == center_lats(num_lats))) then
         ! Values are in decreasing order
         ! This applies to latitude values in a Gaussian grid
         ! Values go from 90-dy to -90+dy
         if (periodic) then
            corner_lats(  1) = 0.5*( 90.0 + center_lats(1))
            corner_lats(num_lats+1) = 0.5*(-90.0 + center_lats(num_lats))
         else
            dy = center_lats(1) - center_lats(2)
            corner_lats(1) = min(90.0, center_lats(1) + 0.5*dy)

            dy = center_lats(num_lats-1) - center_lats(num_lats)
            corner_lats(num_lats+1) = max(-90.0, center_lats(num_lats) - 0.5*dy)
         endif
         do i=2,num_lats
            corner_lats(i) = 0.5*(center_lats(i) + center_lats(i-1))
         enddo
      elseif ((minval(center_lats) == center_lats(1)) .AND. &
              (maxval(center_lats) == center_lats(num_lats))) then
         ! Values are in increasing order
         ! This applies to latitudes in lat/lon grid
         ! Latitudes between -90 and 90
         dy = center_lats(2)-center_lats(1)
         corner_lats( 1)         = max(-90.0, center_lats(1) - 0.5*dy)
         corner_lats(num_lats+1) = min( 90.0, center_lats(num_lats) + 0.5*dy)
         do i=2,num_lats
            corner_lats(i) = 0.5*(center_lats(i) + center_lats(i-1))
         enddo
      endif

      end function determine_lat_corners
!EOC
!------------------------------------------------------------------------------
!BOP
   function create_rectilinear_grid(lon_points, lat_points, grid_name, Nx, Ny, &
                                 coordSys, periodic, staggerloc) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: lat_points(:)
      real,             intent(in)    :: lon_points(:)
      type(ESMF_CoordSys_Flag)        :: coordSys ! Coordinate system: cartesian, spherical
      logical,               optional :: periodic
      type(ESMF_STAGGERLOC), optional :: staggerloc
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create a generic ESMF 2D rectilinear grid given arbitrary latitude and
! longitude grid points.
!
! !LOCAL VARIABLES:
      integer                         :: i, j, num_lats, num_lons
      integer                         :: tlb(2), tub(2)
      integer                         :: rc, STATUS
      real(ESMF_KIND_R8), pointer     :: coordX(:,:)
      real(ESMF_KIND_R8), pointer     :: coordY(:,:)
      type(ESMF_STAGGERLOC)           :: staggerloc_
      logical                         :: periodic_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = "create_rectilinear_grid"
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF Rectilinear Grid: '// TRIM(grid_name)

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      num_lons = SIZE(lon_points)
      call decomposeDim(num_lons, countsPerDEDimX, Nx )

      num_lats = SIZE(lat_points)
      call decomposeDim(num_lats, countsPerDEDimY, Ny )

      ! Create the ESMF Grid
      !---------------------
      periodic_ = .FALSE.
      if (PRESENT(periodic)) periodic_ = periodic
      
      if (periodic_) then
         new_grid = ESMF_GridCreate1PeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1,2/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/1,2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = coordSys, &
                          name            = TRIM(grid_name), rc=rc)

         call LIS_verify(rc, 'ESMF_GridCreate1PeriDim failed in '//TRIM(IAm))
      else
         new_grid = ESMF_GridCreateNoPeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1,2/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/1,2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = coordSys, & 
                          name            = TRIM(grid_name), rc=rc)

         call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed in '//TRIM(IAm))
      endif

      !-------------------------------------------------------------------
      ! Allocate coordinate storage and associate it with the center
      ! stagger location.  Since no coordinate values are specified in
      ! this call no coordinate values are set yet.
      !-------------------------------------------------------------------
      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Get the pointer to the first coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for x-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Get the pointer to the second coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 2, &
                                       farrayPtr           = coordY, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for y-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates 
      !-------------------------------------------------------------------
      do j = tlb(2), tub(2)
         do i = tlb(1), tub(1)
            coordX(i,j) = lon_points(i)
            coordY(i,j) = lat_points(j)
         enddo
      enddo

      write(LIS_logunit,*) '[INFO] Done creating the ESMF Rectilinear Grid: '// TRIM(grid_name)

   end function create_rectilinear_grid
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
   function create_regular_grid_2D(grid_desc, grid_name, Nx, Ny, staggerloc) result(new_grid)
      integer, intent(in) :: Nx, Ny
      character(len=*) :: grid_name
      real,    intent(in) :: grid_desc(:)
      type(ESMF_STAGGERLOC), optional :: staggerloc

      type(ESMF_Grid)                 :: new_grid
      integer                         :: i, j, i1, i2, j1, j2
      integer                         :: numGridPoints(2), tlb(2), tub(2)
      integer                         ::  rc, STATUS
      real(ESMF_KIND_R8)              ::  dx, dy
      real(ESMF_KIND_R8), pointer     :: coordX(:,:)
      real(ESMF_KIND_R8), pointer     :: coordY(:,:)
      type(ESMF_STAGGERLOC)           :: staggerloc_

      real, parameter                 :: D2R = 4.0*ATAN(1.0) / 180.

      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF lat/lon Grid: '// TRIM(grid_name)

      !if (LIS_masterproc) print *, "Create the ESMF Grid: ", TRIM(grid_name)

      numGridPoints(1) = int(grid_desc( 2))   ! Along longitudes
      numGridPoints(2) = int(grid_desc( 3))   ! Along latitudes

      dx = grid_desc( 9) ! (grid_desc(8) - grid_desc(5))/(numGridPoints(1)-1)
      dy = grid_desc(10) ! (grid_desc(7) - grid_desc(4))/(numGridPoints(2)-1)

      !if (LIS_masterproc) print *, "KNJR->"//TRIM(grid_name)//"dx/dy:", dx, dy, grid_desc( 9),grid_desc(10)
      !  -------------------------------------------------------
      new_grid = ESMF_GridCreateNoPeriDim(minIndex       = (/1,1/), &
                                          maxIndex       = numGridPoints, &
                                          gridEdgeLWidth = (/0,0/), &
                                          gridEdgeUWidth = (/0,0/), &
                                          coordSys=ESMF_COORDSYS_CART, & ! Cartesian coordinates (x,y not lat,lon)
                                          indexflag      = ESMF_INDEX_GLOBAL, &
                                          regDecomp      = (/NX, NY /), &
                                          name           = TRIM(grid_name), rc=rc)

	call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed')

      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed')

      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for x-axis')

      call ESMF_GridGetCoord(new_grid, localDE    = 0, &
                                       staggerLoc = staggerloc_, &
                                       coordDim   = 2, &
                                       farrayPtr  = coordY, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for y-axis')

      do j = tlb(2), tub(2)
         do i = tlb(1), tub(1)
            coordX(i,j) = (i-1)*dx  + grid_desc(5)    ! longitude grid points
            coordY(i,j) = (j-1)*dy  + grid_desc(4)    ! latitude  grid points
         enddo
      enddo
      !if (LIS_masterproc) print *, "Done creating the ESMF Grid: ", TRIM(grid_name)
      write(LIS_logunit,*) '[INFO] Done creating the ESMF lat/lon Grid: '// TRIM(grid_name)

   end function create_regular_grid_2D
!EOC
!------------------------------------------------------------------------------
!BOP
   function create_regular_grid(grid_desc, grid_name, Nx, Ny, staggerloc) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: grid_desc(:)
      type(ESMF_STAGGERLOC), optional :: staggerloc
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create an ESMF 2D lat/lon grid
!
! !LOCAL VARIABLES:
      integer                         :: i, j, i1, i2, j1, j2
      integer                         :: numGridPoints(2)
      integer                         :: tlb(1), tub(1)
      integer                         ::  rc, STATUS
      real(ESMF_KIND_R8)              ::  dx, dy
      real(ESMF_KIND_R8), pointer     :: coordX(:)
      real(ESMF_KIND_R8), pointer     :: coordY(:)
      type(ESMF_STAGGERLOC)           :: staggerloc_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = 'create_regular_grid'
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      write(LIS_logunit,*) '[INFO] Create the ESMF lat/lon Grid: '// TRIM(grid_name)

      numGridPoints(1) = int(grid_desc(2))   ! Along longitudes
      numGridPoints(2) = int(grid_desc(3))   ! Along latitudes

      dx = grid_desc( 9) 
      dy = grid_desc(10) 

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      !ALLOCATE(countsPerDEDimX(Nx))
      call decomposeDim(numGridPoints(1), countsPerDEDimX, Nx )

      !ALLOCATE(countsPerDEDimY(Ny))
      call decomposeDim(numGridPoints(2), countsPerDEDimY, Ny )

      ! Create the ESMF Grid
      !---------------------
      new_grid = ESMF_GridCreateNoPeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = ESMF_COORDSYS_CART, & ! use cartesian coordinates
                          name            = TRIM(grid_name), rc=rc)
      call LIS_verify(rc, 'ESMF_GridCreateNoPeriDim failed in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Allocate coordinate storage and associate it with the center
      ! stagger location.  Since no coordinate values are specified in
      ! this call no coordinate values are set yet.
      !-------------------------------------------------------------------
      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call LIS_verify(rc, 'ESMF_GridAddCoord failed in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Get the pointer to the first coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for x-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the first dimension.
      !-------------------------------------------------------------------
      do i = tlb(1), tub(1)
         coordX(i) = (i-1)*dx  + grid_desc(5)    ! longitude grid points
      enddo

      !-------------------------------------------------------------------
      ! Get the pointer to the second coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 2, &
                                       farrayPtr           = coordY, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call LIS_verify(rc, 'ESMF_GridGetCoord failed for y-axis in '//TRIM(IAm))

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the second dimension.
      !-------------------------------------------------------------------
      do j = tlb(1), tub(1)
         coordY(j) = (j-1)*dy  + grid_desc(4)    ! latitude  grid points
      enddo

      write(LIS_logunit,*) '[INFO] Done creating the ESMF lat/lon Grid: '// TRIM(grid_name)

   end function create_regular_grid
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine decomposeDim(dim_world, dim_array, NDEs )
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: dim_world ! total number of grid points
      integer, intent(in)  :: NDEs      ! number of DEs
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: dim_array(0:NDEs-1) 
!
! !DESCRIPTION:
! For a given number of grid points and a number of available processors,
! this subroutines determines the number of grid points assigned to each
! processor.
!
! !LOCAL VARIABLES:
      integer ::   n, im, rm
!EOP
!------------------------------------------------------------------------------
!BOC
      im = dim_world/NDEs
      rm = dim_world-NDEs*im
      do n=0,NDEs-1
                      dim_array(n) = im
      if( n.le.rm-1 ) dim_array(n) = im+1
      enddo
      end subroutine decomposeDim
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: gaussian_find_row
! \label{gaussian_find_row}
!
! !INTERFACE:
function gaussian_find_row(jmax, gaussian_lat_array, lat)
! !USES:
   use LIS_logMod , only : LIS_logunit, LIS_endrun

   implicit none
! !ARGUMENTS:
   integer, intent(in)               :: jmax
   real, dimension(jmax), intent(in) :: gaussian_lat_array
   real, intent(in)                  :: lat
!
! !DESCRIPTION: 
!  This function computes the row number within the global quasi-regular 
!  Gaussian grid corresponding to the latitude given by lat.
!
!   \begin{description}
!   \item [jmax]
!      the total number of latitude circles
!   \item [gaussian\_lat\_array]
!      array of computed latitudes, in degrees
!   \item [lat]
!      latitude to search for
!   \end{description}
!EOP

   integer :: gaussian_find_row

   real    :: eps
   integer :: r

   eps = 180./(2*jmax)

   gaussian_find_row = -9999
   do r = 1, jmax
      if ( abs(gaussian_lat_array(r) - lat) < eps ) then
      !if ( gaussian_latitudes(1,r) == lat ) then
         gaussian_find_row = r
      endif
   enddo

   if ( gaussian_find_row == -9999 ) then
      write(LIS_logunit,fmt=*) '[ERR] gaussian_find_row -- '// &
                               'Could not find lat',lat
      call LIS_endrun
   endif

end function gaussian_find_row
!------------------------------------------------------------------------------

end module LIS_create_gridMod
