!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdasT1534_forcingMod
!BOP
! !MODULE: gdasT1534_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (GDAST1534) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. GDAST1534 forcing variables are produced
!  on a quadratic gaussian grid. 
!  The implementation in LIS has the derived data type {\tt gdasT1534\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GDAST1534 data
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \item[gdasT1534dir]
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
      use ESMF

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDAST1534      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdasT1534_struc
!EOP

  type, public :: gdasT1534_type_dec
     real          :: ts
     integer       :: ncold, nrold   !AWIPS 212 dimensions
     integer       :: nmif
     character*100 :: gdasT1534dir   !GDAST1534 Forcing Directory
     real*8        :: gdasT1534time1,gdasT1534time2
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
     
     integer, allocatable   :: n113(:)
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

     ! For ESMF regridding
     type(ESMF_RouteHandle)       :: routehandle_bilinear
     type(ESMF_DynamicMask)       :: dynamicMask_bilinear
     type(ESMF_RegridMethod_Flag) :: regridMethod_bilinear = ESMF_REGRIDMETHOD_BILINEAR
     type(ESMF_RouteHandle)       :: routehandle_conserve
     type(ESMF_DynamicMask)       :: dynamicMask_conserve
     type(ESMF_RegridMethod_Flag) :: regridMethod_conserve = ESMF_REGRIDMETHOD_CONSERVE_2ND
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

  end type gdasT1534_type_dec
  
  type(gdasT1534_type_dec), allocatable :: gdasT1534_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_GDAST1534
!  \label{init_GDAST1534}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GDAST1534(findex)
! !USES: 
    use LIS_coreMod !,    only : LIS_rc, LIS_domain
    use LIS_logMod !,     only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep

      use LIS_FORC_AttributesMod
    
    implicit none
! !ARGUMENTS: 
    integer,  intent(in)    :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GDAST1534
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}. Based on the GDAST1534 data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the GDAST1534 resolution switches are also 
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
!   \item[readcrd\_gdasT1534](\ref{readcrd_gdasT1534}) \newline
!     reads the runtime options specified for GDAST1534 data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    integer          :: ic, rc

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GDAS-T1534 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(gdasT1534_struc(LIS_rc%nnest))
    
    do n=1, LIS_rc%nnest
       gdasT1534_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, gdasT1534_struc(n)%ts)
    enddo

    call readcrd_gdasT1534()
    
    gdasT1534_struc(:)%nmif    = 16 
    LIS_rc%met_nf(findex) = 16 !number of met variables in GDAST1534 forcing

    do n=1,LIS_rc%nnest

       allocate(gdasT1534_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gdasT1534_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gdasT1534_struc(n)%metdata1 = 0
       gdasT1534_struc(n)%metdata2 = 0

       gdasT1534_struc(n)%findtime1 = 0 
       gdasT1534_struc(n)%findtime2 = 0 

       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
       gridDesci = 0 
       gridDesci(1) = 4
       gridDesci(2) = 3072
       gridDesci(3) = 1536
       gridDesci(4) = 89.90935 ! 89.91153
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) = -89.90935 ! -89.91153
       gridDesci(8) = -0.117187
       gridDesci(9) = 0.117187
       gridDesci(10) = 768.0
       gridDesci(20) = 0
       gdasT1534_struc(n)%ncold = 3072
       gdasT1534_struc(n)%nrold = 1536
       gdasT1534_struc(n)%mi = gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold
       
       IF (LIS_rc%do_esmfRegridding) THEN
          gdasT1534_struc(n)%undefined_value = LIS_rc%udef

          CALL create_gdasT1534Forcing_ESMFobjects(n, findex, gridDesci(1:10))
          CALL create_gdasT1534Model_ESMFobjects(n, findex)
          CALL create_gdasT1534_ESMFroutehandle(n, findex)
       ELSE
          !Setting up weights for Interpolation
          if(LIS_rc%met_interp(findex).eq."bilinear") then 

             allocate(gdasT1534_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
                  gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,&
                  gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
                  gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221)
          elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 
             allocate(gdasT1534_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(gdasT1534_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
                  gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,&
                  gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
                  gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221)
             allocate(gdasT1534_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(gdasT1534_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             call conserv_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n112,gdasT1534_struc(n)%n122,&
                  gdasT1534_struc(n)%n212,gdasT1534_struc(n)%n222,&
                  gdasT1534_struc(n)%w112,gdasT1534_struc(n)%w122,&
                  gdasT1534_struc(n)%w212,gdasT1534_struc(n)%w222)
          elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
             allocate(gdasT1534_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call neighbor_interp_input(n,gridDesci,&
                  gdasT1534_struc(n)%n113)
          else
             write(LIS_logunit,*) 'The specified spatial interpolation option not'
             write(LIS_logunit,*) 'supported for GDAST1534..'
             write(LIS_logunit,*) 'program stopping..'
             call LIS_endrun()
          endif

       ENDIF

    enddo
  end subroutine init_GDAST1534

!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasT1534Forcing_ESMFobjects(n, findex, forcing_gridDesc)
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
      integer          :: rc, ic
      real, pointer    :: lat_centers(:), slat(:), lat_weights(:)
      real, pointer    :: lon_centers(:)
      real, pointer    :: lat_corners(:)
      real, pointer    :: lon_corners(:)
      real             :: dx, min_lon
      integer          :: num_lons, num_lats, m
      real, parameter  :: pi=3.14159265358979
      real, parameter  :: dpr=180.0/pi
      real(ESMF_KIND_R4), pointer :: PTR4(:,:)
      type(ESMF_ArraySpec)        :: arrayspec
      logical                     :: periodic
      logical               :: using_model = .FALSE.
      integer  :: xgrids(LIS_rc%npesx*LIS_rc%npesy)
      integer  :: ygrids(LIS_rc%npesx*LIS_rc%npesy)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize GDAS T1534 ESMF object.'

      periodic = .TRUE.             ! this grid is periodic
      num_lons = forcing_gridDesc(2)
      min_lon  = forcing_gridDesc(5)
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
         !lat_centers(ic) = dpr*asin(slat(num_lats-ic+1))
      enddo

      ALLOCATE(lon_corners(num_lons+1))
      ALLOCATE(lat_corners(num_lats+1))
      lon_corners = determine_lon_corners(lon_centers, periodic)
      lat_corners = determine_lat_corners(lat_centers, periodic)

      if ( (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") .OR. &
           (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") )THEN
         gdasT1534_struc(n)%forcing_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "GDAS T1534 Grid", LIS_rc%npesx, LIS_rc%npesy,  &
                                   ESMF_COORDSYS_SPH_DEG, xgrids, ygrids, using_model, periodic = periodic)
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         gdasT1534_struc(n)%forcing_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "GDAS T1534 Grid", LIS_rc%npesx, LIS_rc%npesy,  &
                                   ESMF_COORDSYS_CART, xgrids, ygrids, using_model, periodic = periodic)
      endif

      DEALLOCATE(lon_centers, lat_centers)
      DEALLOCATE(lon_corners, lat_corners)
      DEALLOCATE(slat, lat_weights)

      ! Interior grid indices
      call getInteriorGrid(gdasT1534_struc(n)%forcing_grid, &
                           gdasT1534_struc(n)%i_min, gdasT1534_struc(n)%i_max, &
                           gdasT1534_struc(n)%j_min, gdasT1534_struc(n)%j_max)

      ! Create the forcing field
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R4)

      gdasT1534_struc(n)%forcing_field = ESMF_FieldCreate(gdasT1534_struc(n)%forcing_grid, arrayspec, &
                              indexflag=ESMF_INDEX_DELOCAL,  &
                              staggerloc=ESMF_STAGGERLOC_CENTER, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "Merra2 Forcing Field", rc=rc)
      call LIS_verify(rc, 'ESMF_FieldCreate failed ')

      call ESMF_FieldGet(gdasT1534_struc(n)%forcing_field, farrayPtr=PTR4, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldGet failed ')
      PTR4 = 0.0

      end subroutine create_gdasT1534Forcing_ESMFobjects
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasT1534Model_ESMFobjects(n, findex)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_field_bundleMod
      use LIS_create_gridMod !,       only : create_gaussian_grid
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
      integer          :: rc, ic, num_lats, num_lons, m
      real             :: model_gridDesc(10)
      real             :: dx
      real, pointer    :: lat_centers(:)
      real, pointer    :: lon_centers(:)
      real, pointer    :: lat_corners(:)
      real, pointer    :: lon_corners(:)
      real(ESMF_KIND_R4), pointer :: PTR4(:,:)
      type(ESMF_ArraySpec)        :: arrayspec
      logical                     :: periodic
      logical               :: using_model = .TRUE.
      integer  :: xgrids(LIS_rc%npesx*LIS_rc%npesy)
      integer  :: ygrids(LIS_rc%npesx*LIS_rc%npesy)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize model ESMF objects.'

      xgrids(:) = LIS_ewe_halo_ind(n,:) - LIS_ews_halo_ind(n,:) + 1 ! LIS_rc%lnc(n)
      ygrids(:) = LIS_nse_halo_ind(n,:) - LIS_nss_halo_ind(n,:) + 1 ! LIS_rc%lnr(n)

      periodic = .TRUE.         ! this grid is periodic
      num_lons = LIS_rc%gnc(n)
      num_lats = LIS_rc%gnr(n)

      ALLOCATE(lon_centers(num_lons))
      ALLOCATE(lat_centers(num_lats))
      lon_centers(:) = LIS_domain(n)%glon(1:num_lons)
      lat_centers(:) = LIS_domain(n)%glat(1:(num_lats-1)*num_lons:num_lons)

      ALLOCATE(lon_corners(num_lons+1))
      ALLOCATE(lat_corners(num_lats+1))
      lon_corners = determine_lon_corners(lon_centers, periodic)
      lat_corners = determine_lat_corners(lat_centers, periodic)

      if ( (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") .OR. &
           (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") )THEN
         gdasT1534_struc(n)%model_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "Model Grid", LIS_rc%npesx, LIS_rc%npesy, &
                                   ESMF_COORDSYS_SPH_DEG, xgrids, ygrids, using_model, periodic = periodic)
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         gdasT1534_struc(n)%model_grid = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "Model Grid", LIS_rc%npesx, LIS_rc%npesy, &
                                   ESMF_COORDSYS_CART, xgrids, ygrids, using_model, periodic = periodic)
      endif

      DEALLOCATE(lon_centers, lat_centers)
      DEALLOCATE(lon_corners, lat_corners)

      ! Create the model field
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R4)

      gdasT1534_struc(n)%model_field = ESMF_FieldCreate(gdasT1534_struc(n)%model_grid, arrayspec, &
                              indexflag=ESMF_INDEX_DELOCAL, &
                              staggerloc=ESMF_STAGGERLOC_CENTER, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "Model Field", rc=rc)
      call LIS_verify(rc, 'ESMF_FieldCreate failed ')

      call ESMF_FieldGet(gdasT1534_struc(n)%model_field, farrayPtr=PTR4, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldGet failed ')
      PTR4 = 0.0

      end subroutine create_gdasT1534Model_ESMFobjects
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_gdasT1534_ESMFroutehandle(n, findex)
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
      type(ESMF_FIELD) :: model_field, gdasT1534_field
      integer          :: ftc(2), ftlb(2), ftub(2)
      real, pointer    :: farray2dd(:,:)
      integer          :: rc
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Determine the ESMF routehandle.'

      if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") THEN
         call createESMF_RouteHandle(gdasT1534_struc(n)%forcing_field, &
                         gdasT1534_struc(n)%model_field,  &
                         gdasT1534_struc(n)%regridMethod_bilinear, &
                         gdasT1534_struc(n)%undefined_value, &
                         gdasT1534_struc(n)%routehandle_bilinear, &
                         gdasT1534_struc(n)%dynamicMask_bilinear, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         write(LIS_logunit,*) '[INFO] Done with the bilinear routehandle.'
      else if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         call createESMF_RouteHandle(gdasT1534_struc(n)%forcing_field, &
                         gdasT1534_struc(n)%model_field,  &
                         gdasT1534_struc(n)%regridMethod_bilinear, &
                         gdasT1534_struc(n)%undefined_value, &
                         gdasT1534_struc(n)%routehandle_bilinear, &
                         gdasT1534_struc(n)%dynamicMask_bilinear, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         write(LIS_logunit,*) '[INFO] Done with the bilinear routehandle.'

         call createESMF_RouteHandle(gdasT1534_struc(n)%forcing_field, &
                         gdasT1534_struc(n)%model_field,  &
                         gdasT1534_struc(n)%regridMethod_conserve, &
                         gdasT1534_struc(n)%undefined_value, &
                         gdasT1534_struc(n)%routehandle_conserve, &
                         gdasT1534_struc(n)%dynamicMask_conserve, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         write(LIS_logunit,*) '[INFO] Done with the conservative routehandle.'
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         call createESMF_RouteHandle(gdasT1534_struc(n)%forcing_field, &
                         gdasT1534_struc(n)%model_field,  &
                         gdasT1534_struc(n)%regridMethod_neighbor, &
                         gdasT1534_struc(n)%undefined_value, &
                         gdasT1534_struc(n)%routehandle_neighbor, &
                         gdasT1534_struc(n)%dynamicMask_neighbor, &
                         lineType=ESMF_LINETYPE_CART)
         write(LIS_logunit,*) '[INFO] Done with the nearest neighbor routehandle.'
      endif

      end subroutine create_gdasT1534_ESMFroutehandle
!EOC
!------------------------------------------------------------------------------

end module gdasT1534_forcingMod

