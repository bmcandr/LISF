!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module merra2_forcingMod
!BOP
! !MODULE: merra2_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the MERRA2 forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt merra2\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[merra2time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2dir]
!    Directory containing the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
! !USES:
      USE ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_merra2      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: merra2_struc

!EOP
  type, public ::  merra2_type_dec
     real         :: ts
     integer      :: ncold, nrold
     character*40 :: merra2dir   !MERRA2 Forcing Directory
     real*8       :: merra2time1,merra2time2
     logical      :: reset_flag

     integer                :: mi
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
     integer                :: findtime1, findtime2
     logical                :: startFlag, dayFlag
     real, allocatable      :: merraforc1(:,:,:,:), merraforc2(:,:,:,:)

     integer            :: nvars
     integer            :: uselml
     integer            :: usecorr

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

     integer                 :: usescalef
     integer                 :: usepcpsampling
     integer                 :: pcpscal_cmo
     integer                 :: use2mwind
     character*100           :: scaleffile
     integer                 :: nbins
     real, allocatable       :: refxrange(:,:,:,:)
     real, allocatable       :: refcdf(:,:,:,:)
     real, allocatable       :: refmean(:,:,:)
     real, allocatable       :: refmean_ip(:)
     real, allocatable       :: refstdev(:,:,:)
     real, allocatable       :: refstdev_ip(:)
     real, allocatable       :: merraxrange(:,:,:,:)
     real, allocatable       :: merracdf(:,:,:,:)
     integer, allocatable    :: rseed(:,:)

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
     type(ESMF_Grid)              :: forcing_grid, forcing_gridCS
     type(ESMF_Grid)              :: model_grid, model_gridCS
     type(ESMF_Field)             :: forcing_field, forcing_fieldCS
     type(ESMF_Field)             :: model_field, model_fieldCS
     real                         :: undefined_value ! for missing value
     integer                      :: i_min, i_max
     integer                      :: j_min, j_max

  end type merra2_type_dec

  type(merra2_type_dec), allocatable :: merra2_struc(:)

      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)

contains

!BOP
!
! !ROUTINE: init_merra2
! \label{init_merra2}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
!
! !INTERFACE:
  subroutine init_merra2(findex)

! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_forecastMod

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MERRA2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_merra2](\ref{readcrd_merra2}) \newline
!     reads the runtime options specified for MERRA2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(LIS_rc%nnest,50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    integer          :: ic, rc
    real             :: forcing_gridDesc(10)

    allocate(merra2_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       merra2_struc(n)%ncold = 576
       merra2_struc(n)%nrold = 361
    enddo

    call readcrd_merra2()
    LIS_rc%met_nf(findex) = 14

    merra2_struc%reset_flag = .false.

    do n=1, LIS_rc%nnest
       merra2_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, merra2_struc(n)%ts)
    enddo

    gridDesci = 0

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0
       gridDesci(n,2) = merra2_struc(n)%ncold
       gridDesci(n,3) = merra2_struc(n)%nrold
       gridDesci(n,4) = -90.000
       gridDesci(n,5) = -180.000
       gridDesci(n,6) = 128
       gridDesci(n,7) = 90.000
       gridDesci(n,8) = 179.375
       gridDesci(n,9) = 0.625
       gridDesci(n,10) = 0.5
       gridDesci(n,20) = 0

       merra2_struc(n)%mi = merra2_struc(n)%ncold*merra2_struc(n)%nrold

       IF (LIS_rc%do_esmfRegridding) THEN
          !merra2_struc(n)%type_kind       = ESMF_TYPEKIND_R4
          merra2_struc(n)%undefined_value = LIS_rc%udef

          forcing_gridDesc(:) = 0
          forcing_gridDesc(2) = gridDesci(n, 2) ! num points along x
          forcing_gridDesc(3) = gridDesci(n, 3) ! num points along y
          forcing_gridDesc(4) = gridDesci(n, 4) ! lower lat
          forcing_gridDesc(5) = gridDesci(n, 5) ! lower lon
          forcing_gridDesc(7) = gridDesci(n, 7) ! upper lat
          forcing_gridDesc(8) = gridDesci(n, 8) ! upper lon
          forcing_gridDesc(9) = gridDesci(n, 9) ! x-grid size
          forcing_gridDesc(10)= gridDesci(n, 10) ! y-grid size

          CALL create_merra2_Forcing_ESMFobjects(n, findex, forcing_gridDesc(1:10))
          CALL create_merra2_Model_ESMFobjects(n, findex)
          CALL create_merra2_ESMFroutehandle(n, findex)
       ELSE
          ! Setting up weights for Interpolation
          if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
             allocate(merra2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n111,merra2_struc(n)%n121,&
                  merra2_struc(n)%n211,merra2_struc(n)%n221,&
                  merra2_struc(n)%w111,merra2_struc(n)%w121,&
                  merra2_struc(n)%w211,merra2_struc(n)%w221)

          elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
             allocate(merra2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(merra2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call bilinear_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n111,merra2_struc(n)%n121,&
                  merra2_struc(n)%n211,merra2_struc(n)%n221,&
                  merra2_struc(n)%w111,merra2_struc(n)%w121,&
                  merra2_struc(n)%w211,merra2_struc(n)%w221)

             allocate(merra2_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(merra2_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             call conserv_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n112,merra2_struc(n)%n122,&
                  merra2_struc(n)%n212,merra2_struc(n)%n222,&
                  merra2_struc(n)%w112,merra2_struc(n)%w122,&
                  merra2_struc(n)%w212,merra2_struc(n)%w222)

          elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
             allocate(merra2_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             call neighbor_interp_input(n, gridDesci(n,:),&
                  merra2_struc(n)%n113)
   
          else
             write(LIS_logunit,*) '[ERR] Interpolation option '// &
                  trim(LIS_rc%met_interp(findex))//&
                  ' for MERRA2 forcing is not supported'
             call LIS_endrun()
          endif
       ENDIF

       call LIS_registerAlarm("MERRA2 forcing alarm", 86400.0, 86400.0)
       merra2_struc(n)%startFlag = .true.
       merra2_struc(n)%dayFlag = .true.

       merra2_struc(n)%nvars = 14

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then 
          
          if(mod(LIS_rc%nensem(n),&
            LIS_forecast_struc(1)%niterations).ne.0) then 
             write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
             write(LIS_logunit,*) '[ERR] of the number of iterations '
             write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
             write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif

          allocate(merra2_struc(n)%merraforc1(&
               LIS_forecast_struc(1)%niterations,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2_struc(n)%merraforc2(&
               LIS_forecast_struc(1)%niterations,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          merra2_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          merra2_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          merra2_struc(n)%nIter = LIS_forecast_struc(1)%niterations
       
          allocate(merra2_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(merra2_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       ! Regular retrospective or non-forecast mode:
       else
          allocate(merra2_struc(n)%merraforc1(1,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2_struc(n)%merraforc2(1,&
               merra2_struc(n)%nvars, 24, &
               LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          merra2_struc(n)%st_iterid = 1
          merra2_struc(n)%en_iterId = 1
          merra2_struc(n)%nIter = 1
       
          allocate(merra2_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(merra2_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
       
       endif

       merra2_struc(n)%metdata1 = 0
       merra2_struc(n)%metdata2 = 0

       merra2_struc(n)%merraforc1 = LIS_rc%udef
       merra2_struc(n)%merraforc2 = LIS_rc%udef

       if ( LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
            LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" ) then
          call read_merra2_elev(n,findex)
       endif

       ! Set up precipitation climate downscaling:
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               merra2_struc(n)%ncold,&
               merra2_struc(n)%nrold)
       endif

    enddo   ! End nest loop

  end subroutine init_merra2

!------------------------------------------------------------------------------
!BOP
      subroutine create_merra2_Forcing_ESMFobjects(n, findex, forcing_gridDesc)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_create_gridMod !,      only : create_regular_grid
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
      real,    intent(in) :: forcing_gridDesc(10)
!
! !DESCRIPTION:
! Create the MERRA2 forcing ESMF grid and field.
! This subroutine might be called several times depending on the integration date.
! Howver, it will be called once for any period when the MERRA2 forcing resolution
! does not change.
!
! !LOCAL VARIABLES:
      integer          :: rc, ic
      real, pointer    :: lat_centers(:)
      real, pointer    :: lon_centers(:)
      real, pointer    :: lat_corners(:)
      real, pointer    :: lon_corners(:)
      real             :: min_lon, min_lat, dx, dy
      integer          :: num_lons, num_lats
      real(kind=sp), pointer :: PTR4(:,:)
      type(ESMF_ArraySpec)  :: arrayspec
      logical               :: periodic
      logical               :: using_model = .FALSE.
      integer  :: xgrids(LIS_rc%npesx*LIS_rc%npesy)
      integer  :: ygrids(LIS_rc%npesx*LIS_rc%npesy)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize MERRA2 ESMF objects.'

      periodic = .TRUE.         ! this grid is periodic
      num_lons = forcing_gridDesc(2)
      min_lon  = forcing_gridDesc(5)
      dx       = forcing_gridDesc(9)

      num_lats = forcing_gridDesc(3)
      min_lat  = forcing_gridDesc(4)
      dy       = forcing_gridDesc(10)

      ALLOCATE(lon_centers(num_lons))
      do ic = 1, num_lons
         lon_centers(ic) = (ic-1)*dx + min_lon
      enddo

      ALLOCATE(lat_centers(num_lats))
      do ic = 1, num_lats
         lat_centers(ic) = (ic-1)*dy + min_lat
      enddo

      ALLOCATE(lon_corners(num_lons+1))
      ALLOCATE(lat_corners(num_lats+1))
      lon_corners = determine_lon_corners(lon_centers, periodic)
      lat_corners = determine_lat_corners(lat_centers, periodic)

      merra2_struc(n)%forcing_grid = createRectilinearGrid(lon_centers, lat_centers, &
                               lon_corners, lat_corners, &
                               "MERRA2 Grid", LIS_rc%npesx, LIS_rc%npesy,  &
                                ESMF_COORDSYS_CART, xgrids, ygrids, using_model, periodic   = periodic)

      if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         merra2_struc(n)%forcing_gridCS = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "MERRA2 Grid Conservative", LIS_rc%npesx, LIS_rc%npesy,  &
                                   ESMF_COORDSYS_SPH_DEG, xgrids, ygrids, using_model, periodic   = periodic)
      endif

      DEALLOCATE(lon_centers, lat_centers)
      DEALLOCATE(lon_corners, lat_corners)

      ! Interior grid indices
      call getInteriorGrid(merra2_struc(n)%forcing_grid, &
                           merra2_struc(n)%i_min, merra2_struc(n)%i_max, &
                           merra2_struc(n)%j_min, merra2_struc(n)%j_max)

      ! Create the forcing field
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=merra2_struc(n)%type_kind)

      merra2_struc(n)%forcing_field = ESMF_FieldCreate(merra2_struc(n)%forcing_grid, arrayspec, &
                                 indexflag=ESMF_INDEX_DELOCAL,  &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "Merra2 Forcing Field", rc=rc)
      call LIS_verify(rc, 'ESMF_FieldCreate failed ')

      if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         merra2_struc(n)%forcing_fieldCS = ESMF_FieldCreate(merra2_struc(n)%forcing_gridCS, &
                                    arrayspec, &
                                    indexflag=ESMF_INDEX_DELOCAL,  &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                                 name = "Merra2 Forcing Field Conservative", rc=rc)
         call LIS_verify(rc, 'ESMF_FieldCreate failed ')
      endif

      call ESMF_FieldGet(merra2_struc(n)%forcing_field, farrayPtr=PTR4, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldGet failed ')
      PTR4 = 0.0
      if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         call ESMF_FieldGet(merra2_struc(n)%forcing_fieldCS, farrayPtr=PTR4, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldGet failed ')
         PTR4 = 0.0
      endif

      write(LIS_logunit,*) '[INFO] Done Initialize MERRA2 ESMF objects.'
      end subroutine create_merra2_Forcing_ESMFobjects
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_merra2_Model_ESMFobjects(n, findex)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun
      use LIS_FORC_AttributesMod
      use LIS_create_gridMod !,       only : create_regular_grid
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
      real             :: model_gridDesc(10)
      real, pointer    :: lat_centers(:)
      real, pointer    :: lon_centers(:)
      real, pointer    :: lat_corners(:)
      real, pointer    :: lon_corners(:)
      real             :: dx, dy
      real(kind=sp), pointer  :: PTR4(:,:)
      type(ESMF_ArraySpec)   :: arrayspec
      logical                :: periodic
      logical               :: using_model = .TRUE.
      integer  :: xgrids(LIS_rc%npesx*LIS_rc%npesy)
      integer  :: ygrids(LIS_rc%npesx*LIS_rc%npesy)
!EOP
!------------------------------------------------------------------------------
!BOC
      write(LIS_logunit,*) '[INFO] Initialize model ESMF objects.'

      xgrids(:) = LIS_ewe_halo_ind(n,:) - LIS_ews_halo_ind(n,:) + 1 ! LIS_rc%lnc(n)
      ygrids(:) = LIS_nse_halo_ind(n,:) - LIS_nss_halo_ind(n,:) + 1 ! LIS_rc%lnr(n)

      periodic = .FALSE.         ! this grid is not periodic
      num_lons = LIS_rc%gnc(n)
      num_lats = LIS_rc%gnr(n)

      ! Get the latitudes and longitudes at centers from the model
      ALLOCATE(lon_centers(num_lons))
      ALLOCATE(lat_centers(num_lats))
      lon_centers(:) = LIS_domain(n)%glon(1:num_lons)
      lat_centers(:) = LIS_domain(n)%glat(1:(num_lats-1)*num_lons:num_lons)

      ! Determine the latitudes and longitudes at corners
      ALLOCATE(lon_corners(num_lons+1))
      ALLOCATE(lat_corners(num_lats+1))
      lon_corners = determine_lon_corners(lon_centers, periodic)
      lat_corners = determine_lat_corners(lat_centers, periodic)

      merra2_struc(n)%model_grid = createRectilinearGrid(lon_centers, lat_centers, &
                               lon_corners, lat_corners, &
                               "Model Grid", LIS_rc%npesx, LIS_rc%npesy, &
                                ESMF_COORDSYS_CART, xgrids, ygrids, using_model, periodic = periodic)

      if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         merra2_struc(n)%model_gridCS = createRectilinearGrid(lon_centers, lat_centers, &
                                  lon_corners, lat_corners, &
                                  "Model Grid Conservative", LIS_rc%npesx, LIS_rc%npesy, &
                                   ESMF_COORDSYS_SPH_DEG, xgrids, ygrids, using_model, periodic = periodic)
      endif

      DEALLOCATE(lon_centers, lat_centers)
      DEALLOCATE(lon_corners, lat_corners)

      ! Create the model field
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=merra2_struc(n)%type_kind)

      merra2_struc(n)%model_field = ESMF_FieldCreate(merra2_struc(n)%model_grid, arrayspec, &
                                 indexflag=ESMF_INDEX_DELOCAL,  &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "Model Field", rc=rc)
      call LIS_verify(rc, 'ESMF_FieldCreate failed ')

      if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         merra2_struc(n)%model_fieldCS = ESMF_FieldCreate(merra2_struc(n)%model_gridCS, &
                                 arrayspec, &
                                 indexflag=ESMF_INDEX_DELOCAL,  &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                                 name = "Model Field Conservative", rc=rc)
         call LIS_verify(rc, 'ESMF_FieldCreate failed ')
      endif

      call ESMF_FieldGet(merra2_struc(n)%model_field, farrayPtr=PTR4, rc=rc)
      call LIS_verify(rc, 'ESMF_FieldGet failed ')
      PTR4 = 0.0
      if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         call ESMF_FieldGet(merra2_struc(n)%model_fieldCS, farrayPtr=PTR4, rc=rc)
         call LIS_verify(rc, 'ESMF_FieldGet failed ')
         PTR4 = 0.0
      endif

      end subroutine create_merra2_Model_ESMFobjects
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine create_merra2_ESMFroutehandle(n, findex)
!
! !USES:
      use ESMF
      use LIS_coreMod !,    only : LIS_rc, LIS_domain
      use LIS_logMod !,     only : LIS_logunit, LIS_endrun

      use LIS_FORC_AttributesMod
      use LIS_ESMF_Regrid_Utils, only : createESMF_RouteHandle
!
! !INPUT PARAMETERS:
      integer, intent(in) :: n
      integer, intent(in) :: findex
! 
! !DESCRIPTION:
! Determine the ESMF routehandle needed for thr regridding between 
! the MERRA2 forcing and the model.
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
         call createESMF_RouteHandle(merra2_struc(n)%forcing_field, &
                         merra2_struc(n)%model_field,  &
                         merra2_struc(n)%regridMethod_bilinear, &
                         merra2_struc(n)%undefined_value, &
                         merra2_struc(n)%routehandle_bilinear, &
                         merra2_struc(n)%dynamicMask_bilinear, &
                         lineType=ESMF_LINETYPE_CART)
         write(LIS_logunit,*) '[INFO] Done with the bilinear routehandle.'
      else if (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") THEN
         call createESMF_RouteHandle(merra2_struc(n)%forcing_field, &
                         merra2_struc(n)%model_field,  &
                         merra2_struc(n)%regridMethod_bilinear, &
                         merra2_struc(n)%undefined_value, &
                         merra2_struc(n)%routehandle_bilinear, &
                         merra2_struc(n)%dynamicMask_bilinear, &
                         lineType=ESMF_LINETYPE_CART)
         write(LIS_logunit,*) '[INFO] Done with the bilinear routehandle.'

         call createESMF_RouteHandle(merra2_struc(n)%forcing_fieldCS, &
                         merra2_struc(n)%model_fieldCS,  &
                         merra2_struc(n)%regridMethod_conserve, &
                         merra2_struc(n)%undefined_value, &
                         merra2_struc(n)%routehandle_conserve, &
                         merra2_struc(n)%dynamicMask_conserve, &
                         lineType=ESMF_LINETYPE_GREAT_CIRCLE)
         write(LIS_logunit,*) '[INFO] Done with the conservative routehandle.'
      else if (trim(LIS_rc%met_interp(findex)) .eq. "neighbor") THEN
         call createESMF_RouteHandle(merra2_struc(n)%forcing_field, &
                         merra2_struc(n)%model_field,  &
                         merra2_struc(n)%regridMethod_neighbor, &
                         merra2_struc(n)%undefined_value, &
                         merra2_struc(n)%routehandle_neighbor, &
                         merra2_struc(n)%dynamicMask_neighbor, &
                         lineType=ESMF_LINETYPE_CART)
         write(LIS_logunit,*) '[INFO] Done with the nearest neighbor routehandle.'
      endif

      end subroutine create_merra2_ESMFroutehandle
!EOC
!------------------------------------------------------------------------------
end module merra2_forcingMod
