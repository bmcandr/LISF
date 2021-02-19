!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdas_forcingMod
!BOP
! !MODULE: gdas_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (GDAS) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. GDAS forcing variables are produced
!  on a quadratic gaussian grid. LIS uses the 00, 03, 06 and as
!  needed, the 09 forecasts. The forecasts are produced at 6 hr intervals.
!  The resolution of GDAS forcing varies as follows:
!
!   upto 2000/1/24          :   T126 (384x190)  grid \newline
!   2001/01/24 - 2002/10/29 :   T170 (512x256)  grid \newline
!   2002/10/29 - 2005/05/31 :   T254 (768x384)  grid \newline
!   2005/05/31 - 2010/07/27 :   T382 (1152x576) grid \newline
!   2010/07/28 - 2015/01/14 :   T574 (1760x880) grid \newline
!   2015/01/14 - onwards    :  T1534 (3072x1536) grid
!
!  On 2019/06/12 12Z, GDAS removed precipitation fields from the f00 data
!  files. The data fields in these files are now all instantaneous values.
!  When the reader is using data files after this time, a new subroutine will
!  be used that excludes precipitation as well as reads in instantaneous radiation
!  data. For data files prior to the switch, the reader will use the old subroutine.
!
!  The implementation in LIS has the derived data type {\tt gdas\_struc} that
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GDAS data
!  \item[gdastime1]
!    The nearest, previous 3 hour instance of the incoming
!    data (as a real time).
!  \item[gdastime2]
!    The nearest, next 3 hour instance of the incoming
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch GDAS resolution to T126
!  \item[griduptime2]
!    The time to switch GDAS resolution to T170
!  \item[griduptime3]
!    The time to switch GDAS resolution to T254
!  \item[griduptime4]
!    The time to switch GDAS resolution to T382
!  \item[griduptime5]
!    The time to switch GDAS resolution to T574
!  \item[griduptime6]
!    The time to switch GDAS resolution to T1534
!  \item[datastructime1]
!    The time to switch to new data structure for f00 files
!    that removed precipitation fields.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \item[gdasdir]
!    Directory containing the input data
!  \item[elevfile]
!    File with the elevation definition for the input data.
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES:
      USE ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDAS      !defines the native resolution of
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdas_struc
  public :: create_gdasForcing_ESMFobjects
  public :: create_gdas_ESMFroutehandle
!EOP

  type, public :: gdas_type_dec
     real          :: ts
     integer       :: ncold, nrold   !AWIPS 212 dimensions
     integer       :: nmif
     character*100 :: gdasdir   !GDAS Forcing Directory
     character*50  :: met_interp

     real*8        :: gdastime1, gdastime2
     real*8        :: griduptime1, griduptime2, griduptime3
     real*8        :: griduptime4, griduptime5, griduptime6
     real*8        :: datastructime1
     logical       :: gridchange1, gridchange2, gridchange3
     logical       :: gridchange4, gridchange5, gridchange6
     logical       :: dstrucchange1
     integer       :: mi
     integer       :: findtime1, findtime2

     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

     ! For ESMF regridding
     type(ESMF_RouteHandle)       :: routehandle_bilinear
     type(ESMF_DynamicMask)       :: dynamicMask_bilinear
     type(ESMF_RegridMethod_Flag) :: regridMethod_bilinear = ESMF_REGRIDMETHOD_BILINEAR
     type(ESMF_RouteHandle)       :: routehandle_conserve
     type(ESMF_DynamicMask)       :: dynamicMask_conserve
     type(ESMF_RegridMethod_Flag) :: regridMethod_conserve = ESMF_REGRIDMETHOD_CONSERVE
     type(ESMF_RouteHandle)       :: routehandle_neighbor
     type(ESMF_DynamicMask)       :: dynamicMask_neighbor
     type(ESMF_RegridMethod_Flag) :: regridMethod_neighbor = ESMF_REGRIDMETHOD_NEAREST_STOD
     type(ESMF_TypeKind_Flag)     :: type_kind = ESMF_TYPEKIND_R4
     type(ESMF_Grid)              :: forcing_grid
     type(ESMF_Grid)              :: model_grid
     type(ESMF_Field)             :: forcing_field
     type(ESMF_Field)             :: model_field
     real                         :: undefined_value ! for missing value
     integer                      :: i_min, i_max
     integer                      :: j_min, j_max

  end type gdas_type_dec

  type(gdas_type_dec), allocatable :: gdas_struc(:)

      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)

contains

!BOP
!
! !ROUTINE: init_GDAS
!  \label{init_GDAS}
!
! !REVISION HISTORY:
! 11Dec2003: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine init_GDAS(findex)
! !USES:
    use LIS_coreMod !,    only : LIS_rc, LIS_domain
    use LIS_logMod !,     only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep

      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod

    implicit none
! !ARGUMENTS:
    integer,  intent(in)    :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for GDAS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}). Based on the GDAS data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the GDAS resolution switches are also
!  defined in this routine.
!
!  The arguments are:
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_gdas](\ref{readcrd_gdas}) \newline
!     reads the runtime options specified for GDAS data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!  \end{description}
!EOP
    real :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    integer          :: ic, rc
    real             :: forcing_gridDesc(10)

    character(len=100)         :: name_tracer
    integer                    :: fieldcount

    write(LIS_logunit,*) '[INFO] Defines the native resolution '
    write(LIS_logunit,*) '       of the input forcing for GDAS data'

    allocate(gdas_struc(LIS_rc%nnest))

    do n=1, LIS_rc%nnest
       gdas_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, gdas_struc(n)%ts)
    enddo

    call readcrd_gdas()

    gdas_struc(:)%nmif    = 9
    LIS_rc%met_nf(findex) = 9 !number of met variables in GDAS forcing

    gdas_struc(:)%ncold = 192
    gdas_struc(:)%nrold = 94

    do n=1,LIS_rc%nnest

       allocate(gdas_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gdas_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gdas_struc(n)%metdata1 = 0
       gdas_struc(n)%metdata2 = 0

       gdas_struc(n)%findtime1 = 0
       gdas_struc(n)%findtime2 = 0

       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
       gridDesci = 0
       gridDesci(1) = 4
       gridDesci(2) = 192
       gridDesci(3) = 94
       gridDesci(4) = 88.542
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) =  -88.542
       gridDesci(8) = -1.875
       gridDesci(9) = 1.875
       gridDesci(10) = 47
       gridDesci(20) = 0
       gdas_struc(n)%mi = gdas_struc(n)%ncold*gdas_struc(n)%nrold

       ! This grid is good for some time in the 1990's.
       ! Look up the exact dates.
       yr1 = 1991
       mo1 = 01
       da1 = 01
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time( gdas_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2000
       mo1 = 01
       da1 = 24
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time( gdas_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2002     !grid update time ~ 0.469
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2005     !grid update time ~ 0.313
       mo1 = 05
       da1 = 31
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2010     !grid update time ~ 0.205
       mo1 = 07
       da1 = 28
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2015     !grid update time ~ 0.117
       mo1 = 01
       da1 = 14
       hr1 = 6
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       ! Set time for f00 data structure change
       yr1 = 2019
       mo1 = 06
       da1 = 12
       hr1 = 9 !09Z is when the reader reads in the 12Zf00 file
       mn1 = 0; ss1 = 0
       call LIS_date2time(gdas_struc(n)%datastructime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)
 
       gdas_struc(n)%gridchange1 = .true.
       gdas_struc(n)%gridchange2 = .true.
       gdas_struc(n)%gridchange3 = .true.
       gdas_struc(n)%gridchange4 = .true.
       gdas_struc(n)%gridchange5 = .true.
       gdas_struc(n)%gridchange6 = .true.

       gdas_struc(n)%dstrucchange1 = .true.

       ! Initialize ESMF objects
       !------------------------
       IF (LIS_rc%do_esmfRegridding) THEN
          !gdas_struc(n)%type_kind       = ESMF_TYPEKIND_R4
          gdas_struc(n)%undefined_value = LIS_rc%udef

          call create_gdasModel_ESMFobjects(n, findex)
       ELSE
          ! Setting up weights for Interpolation
          call gdas_reset_interp_input(n, findex, gridDesci)
       ENDIF

    enddo

  end subroutine init_GDAS

!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasForcing_ESMFobjects(n, findex, forcing_gridDesc)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod !,      only :create_gaussian_grid
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
      real,    intent(in) :: forcing_gridDesc(10)
!
! !DESCRIPTION:
! Create the GDAS forcing ESMF grid and field.
! This subroutine might be called several times depending on the integration date.
! Howver, it will be called once for any period when the GDAS forcing resolution
! does not change.
!
! !LOCAL VARIABLES:
      character(len=2) :: num_st
      integer          :: rc, ic
      type(ESMF_Grid)  :: gdas_grid
      real             :: dummy_array(1,1) = 0.0
      real, pointer    :: lat_centers(:), slat(:), lat_weights(:)
      real, pointer    :: lon_centers(:)
      real, pointer    :: lat_corners(:)
      real, pointer    :: lon_corners(:)
      real             :: dx, min_lon, min_lat, max_lat
      integer          :: num_lons, num_lats
      real, parameter  :: pi=3.14159265358979
      real, parameter  :: dpr=180.0/pi
      real(kind=4), pointer        :: PTR4(:,:)
      type(ESMF_ArraySpec)         :: arrayspec
      logical               :: periodic
      logical               :: using_model = .FALSE.
      integer  :: xgrids(LIS_rc%npesx*LIS_rc%npesy)
      integer  :: ygrids(LIS_rc%npesx*LIS_rc%npesy)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize GDAS ESMF object.'

      periodic = .TRUE.               ! this grid is periodic
      num_lons = forcing_gridDesc(2)
      min_lat  = forcing_gridDesc(4)
      min_lon  = forcing_gridDesc(5)
      max_lat  = forcing_gridDesc(7)
      dx       = forcing_gridDesc(9)

      ALLOCATE(lon_centers(num_lons))
      do ic = 1, num_lons
         lon_centers(ic) = (ic-1)*dx + min_lon
         ! values must be between -180 to 180
         IF (lon_centers(ic) > 180.0) lon_centers(ic) = lon_centers(ic) - 360.0
      enddo

      ! Determine the global Gaussian latitude grid points
      num_lats = forcing_gridDesc(3)

      ALLOCATE(lat_centers(num_lats))
      ALLOCATE(slat(num_lats))
      ALLOCATE(lat_weights(num_lats))

      call gausslat(num_lats, slat, lat_weights)
      do ic = 1, num_lats
         lat_centers(ic) = dpr*asin(slat(ic))
      enddo

      ALLOCATE(lon_corners(num_lons+1))
      ALLOCATE(lat_corners(num_lats+1))
      lon_corners = determine_lon_corners(lon_centers, periodic)
      lat_corners = determine_lat_corners(lat_centers, periodic)

      if ( (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") .OR. &
           (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") )THEN
         gdas_struc(n)%forcing_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "GDAS Grid", LIS_rc%npesx, LIS_rc%npesy,  &
                                  ESMF_COORDSYS_SPH_DEG, &
                                  xgrids, ygrids, using_model, &
                                  periodic   = periodic)
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         gdas_struc(n)%forcing_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "GDAS Grid", LIS_rc%npesx, LIS_rc%npesy,  &
                                  ESMF_COORDSYS_CART, &
                                  xgrids, ygrids, using_model, &
                                  periodic   = periodic)
      endif

      DEALLOCATE(lon_centers, lat_centers)
      DEALLOCATE(lon_corners, lat_corners)
      DEALLOCATE(slat, lat_weights)

      ! Interior grid indices
      call getInteriorGrid(gdas_struc(n)%forcing_grid, &
                           gdas_struc(n)%i_min, gdas_struc(n)%i_max, &
                           gdas_struc(n)%j_min, gdas_struc(n)%j_max)
                          

      ! Create the forcing field
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=gdas_struc(n)%type_kind)

      gdas_struc(n)%forcing_field = ESMF_FieldCreate(gdas_struc(n)%forcing_grid, arrayspec, &
                                 indexflag=ESMF_INDEX_DELOCAL,  &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "Merra2 Forcing Field", rc=rc)
      call LIS_verify(rc, 'ESMF_FieldCreate failed ')

      call ESMF_FieldGet(gdas_struc(n)%forcing_field, farrayPtr=PTR4, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldGet failed ')
      PTR4 = 0.0

      end subroutine create_gdasForcing_ESMFobjects
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasModel_ESMFobjects(n, findex)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod 
!
! !INPUT PRAMETERS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
!
! !DESCRIPTION:
! Create the model ESMF grid and field.
! This subroutine is called once as the model resolution does not change.
!
! !LOCAL VARIABLES:
      integer          :: rc, ic, num_lats, num_lons
      real, pointer    :: lat_centers(:)
      real, pointer    :: lon_centers(:)
      real, pointer    :: lat_corners(:)
      real, pointer    :: lon_corners(:)
      real(kind=4), pointer :: PTR4(:,:)
      type(ESMF_ArraySpec)  :: arrayspec
      real                  :: dx, dy
      logical               :: periodic
      logical               :: using_model = .TRUE.
      integer  :: xgrids(LIS_rc%npesx*LIS_rc%npesy)
      integer  :: ygrids(LIS_rc%npesx*LIS_rc%npesy)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize GDAS and model ESMF objects.'

      xgrids(:) = LIS_ewe_halo_ind(n,:) - LIS_ews_halo_ind(n,:) + 1 ! LIS_rc%lnc(n)
      ygrids(:) = LIS_nse_halo_ind(n,:) - LIS_nss_halo_ind(n,:) + 1 ! LIS_rc%lnr(n)

      periodic = .FALSE.    ! the grid is not periodic

      ! ---> ESMF grid for for the model
      num_lons = LIS_rc%gnc(n)
      num_lats = LIS_rc%gnr(n)

      ALLOCATE(lon_centers(num_lons))
      lon_centers(:) = LIS_domain(n)%glon(1:num_lons)

      ALLOCATE(lat_centers(num_lats))
      lat_centers(:) = LIS_domain(n)%glat(1:(num_lats-1)*num_lons:num_lons)

      ALLOCATE(lon_corners(num_lons+1))
      ALLOCATE(lat_corners(num_lats+1))
      lon_corners = determine_lon_corners(lon_centers, periodic)
      lat_corners = determine_lat_corners(lat_centers, periodic)

      if ( (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") .OR. &
           (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") )THEN
         gdas_struc(n)%model_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                   lon_corners, lat_corners, &
                                   "Model Grid", LIS_rc%npesx, LIS_rc%npesy, &
                                   ESMF_COORDSYS_SPH_DEG, &
                                   xgrids, ygrids, using_model, &
                                   periodic = periodic)
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         gdas_struc(n)%model_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                   lon_corners, lat_corners, &
                                   "Model Grid", LIS_rc%npesx, LIS_rc%npesy, &
                                   ESMF_COORDSYS_CART, &
                                   xgrids, ygrids, using_model, &
                                   periodic = periodic)
      endif


      DEALLOCATE(lon_centers, lat_centers)
      DEALLOCATE(lon_corners, lat_corners)

      ! Create the model field
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=gdas_struc(n)%type_kind)

      gdas_struc(n)%model_field = ESMF_FieldCreate(gdas_struc(n)%model_grid, arrayspec, &
                                 indexflag=ESMF_INDEX_DELOCAL,  &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "Model Field", rc=rc)
      call LIS_verify(rc, 'ESMF_FieldCreate failed ')

      call ESMF_FieldGet(gdas_struc(n)%model_field, farrayPtr=PTR4, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldGet failed ')
      PTR4 = 0.0

      end subroutine create_gdasModel_ESMFobjects
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdas_ESMFroutehandle(n, findex)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun

      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_ESMF_Regrid_Utils, only : createESMF_RouteHandle
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
! 
! !DESCRIPTION:
! Determine the ESMF routehandle needed for thr regridding between 
! the GDAS forcing and the model.
!
! !LOCAL VARIABLES:
      integer          :: ftc(2), ftlb(2), ftub(2)
      real, pointer    :: farray2dd(:,:)
      integer          :: rc
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Determine the ESMF routehandle.'

      if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") THEN
         call createESMF_RouteHandle(gdas_struc(n)%forcing_field, &
                         gdas_struc(n)%model_field,  &
                         gdas_struc(n)%regridMethod_bilinear, &
                         gdas_struc(n)%undefined_value, &
                         gdas_struc(n)%routehandle_bilinear, &
                         gdas_struc(n)%dynamicMask_bilinear, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         !write(LIS_logunit,*) '[INFO] Done with the bilinear routehandle.'
      else if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         call createESMF_RouteHandle(gdas_struc(n)%forcing_field, &
                         gdas_struc(n)%model_field,  &
                         gdas_struc(n)%regridMethod_bilinear, &
                         gdas_struc(n)%undefined_value, &
                         gdas_struc(n)%routehandle_bilinear, &
                         gdas_struc(n)%dynamicMask_bilinear, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         !write(LIS_logunit,*) '[INFO] Done with the bilinear routehandle.'

         call createESMF_RouteHandle(gdas_struc(n)%forcing_field, &
                         gdas_struc(n)%model_field,  &
                         gdas_struc(n)%regridMethod_conserve, &
                         gdas_struc(n)%undefined_value, &
                         gdas_struc(n)%routehandle_conserve, &
                         gdas_struc(n)%dynamicMask_conserve, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         !write(LIS_logunit,*) '[INFO] Done with the conservative routehandle.'
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         call createESMF_RouteHandle(gdas_struc(n)%forcing_field, &
                         gdas_struc(n)%model_field,  &
                         gdas_struc(n)%regridMethod_neighbor, &
                         gdas_struc(n)%undefined_value, &
                         gdas_struc(n)%routehandle_neighbor, &
                         gdas_struc(n)%dynamicMask_neighbor)
         !write(LIS_logunit,*) '[INFO] Done with the nearest neighbor routehandle.'
      endif

      end subroutine create_gdas_ESMFroutehandle
!EOC
!------------------------------------------------------------------------------

end module gdas_forcingMod

