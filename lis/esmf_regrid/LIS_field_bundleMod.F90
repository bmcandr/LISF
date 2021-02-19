!------------------------------------------------------------------------------
! NASA GSFC 
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: LIS_field_bundleMod
!
! !INTERFACE:
!
      module LIS_field_bundleMod
!
! !USES:
      use ESMF
      use LIS_logMod,     only : LIS_verify

      implicit none

! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: addTracerToBundle
      public  :: updateTracerToBundle
      public  :: getPointerFromBundle
      public  :: ESMF_recordBundle
      public  :: copy_esmfBundle
      public  :: getPointerFromField

      interface addTracerToBundle
         module procedure addTracerToBundle3D
         module procedure addTracerToBundle2D
      end interface

      interface updateTracerToBundle
         module procedure updateTracerToBundle_ByIndex3D_r4
         module procedure updateTracerToBundle_ByIndex3D_r8
         module procedure updateTracerToBundle_ByName3D_r4
         module procedure updateTracerToBundle_ByName3D_r8
         module procedure updateTracerToBundle_ByIndex2D_r4
         module procedure updateTracerToBundle_ByIndex2D_r8
         module procedure updateTracerToBundle_ByName2D_r4
         module procedure updateTracerToBundle_ByName2D_r8
      end interface

      interface getPointerFromBundle
         module procedure getPointerFromBundle_ByIndex3D_r4
         module procedure getPointerFromBundle_ByIndex3D_r8
         module procedure getPointerFromBundle_ByName3D_r4
         module procedure getPointerFromBundle_ByName3D_r8
         module procedure getPointerFromBundle_ByIndex2D_r4
         module procedure getPointerFromBundle_ByIndex2D_r8
         module procedure getPointerFromBundle_ByName2D_r4
         module procedure getPointerFromBundle_ByName2D_r8
      end interface

      interface getPointerFromField
         module procedure getPointerFromField_r4_3D
         module procedure getPointerFromField_r8_3D
         module procedure getPointerFromField_r4_2D
         module procedure getPointerFromField_r8_2D
      end interface
      

      type ESMF_recordBundle
          type(ESMF_FieldBundle) :: prevBundle ! for the previous record
          type(ESMF_FieldBundle) :: nextBundle ! for the next     record
      end type ESMF_recordBundle

      integer, PARAMETER :: MAXSTR = ESMF_MAXSTR
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)

! !DESCRIPTION:
! Basic routines for manipulating ESMF bundles.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
      subroutine copy_esmfBundle(dstBundle, srcBundle)
!
! !INPUT PARAMETERS:
      type(ESMF_FieldBundle),  intent(in)    :: srcBundle
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_FieldBundle),  intent(inOut)    :: dstBundle
!
! !DESCRIPTION:
! ! Copy the source ESMF bundle into the destination ESMF bundle.
!
! !LOCAL VARIABLES:
      integer :: fieldCount, fieldRank, j, rc
      !real(kind=sp), pointer :: var2d_dstR4(:,:), var2d_srcR4(:,:)
      !real(kind=sp), pointer :: var2d_dstR4(:,:,:), var2d_srcR4(:,:,:)
      !real(kind=dp), pointer :: var2d_dstR8(:,:), var2d_srcR8(:,:)
      !real(kind=dp), pointer :: var2d_dstR8(:,:,:), var2d_srcR8(:,:,:)
      !type(ESMF_Array)       :: array
      !type(ESMF_TypeKind_Flag)      :: kind_type
      type (ESMF_Field)     :: field_dst, field_src
      character(len=MAXSTR), allocatable :: fields_names(:)   
      character(len=MAXSTR) :: IAm='copy_esmfBundle'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Get the number of fields in the destination bundle
      call ESMF_FieldBundleGet(dstBundle, fieldCount = fieldCount, rc=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed ')

      ! Get the list of field names in the destination bundle
      allocate(fields_names(fieldCount))
      call ESMF_FieldBundleGet(dstBundle, fieldNameList = fields_names, rc=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed ')

      ! Loop over the fields to copy them
      do j = 1,fieldCount
         call ESMF_FieldBundleGet(srcBundle, TRIM(fields_names(j)), field=field_src, rc=rc)
         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fields_names(j)))

         call ESMF_FieldBundleGet(dstBundle, TRIM(fields_names(j)), field=field_dst, rc=rc)
         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fields_names(j)))

         call ESMF_FieldCopy(field_dst, field_src, rc=rc)
         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldCopy failed for '//TRIM(fields_names(j)))

         !call ESMF_FieldGet  (field_dst, array=array, name=fields_names(j), RC=STATUS)
         !call ESMF_ArrayGet(array,typekind=kind_type, rc=status )

         !call ESMF_FieldGet(field_dst, dimCount=fieldRank, rc=rc)
         !call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed to get field rank ')

         !if (fieldRank == 2) then
         !   if (kind_type == ESMF_TYPEKIND_R4) THEN
         !      call ESMF_FieldGet(field_src, localDE=0, farrayPtr=var2d_srcR4, rc=rc)
         !      call ESMF_FieldGet(field_dst, localDE=0, farrayPtr=var2d_dstR4, rc=rc)
         !      var2d_dstR4 = var2d_srcR4
         !   else
         !      call ESMF_FieldGet(field_src, localDE=0, farrayPtr=var2d_srcR8, rc=rc)
         !      call ESMF_FieldGet(field_dst, localDE=0, farrayPtr=var2d_dstR8, rc=rc)
         !      var2d_dstR8 = var2d_srcR8
         !   endif
         !else if (fieldRank == 3) then
         !   if (kind_type == ESMF_TYPEKIND_R4) THEN
         !      call ESMF_FieldGet(field_src, localDE=0, farrayPtr=var3d_srcR4, rc=rc)
         !      call ESMF_FieldGet(field_dst, localDE=0, farrayPtr=var3d_dstR4, rc=rc)
         !      var3d_dstR4 = var3d_srcR4
         !   else
         !      call ESMF_FieldGet(field_src, localDE=0, farrayPtr=var3d_srcR8, rc=rc)
         !      call ESMF_FieldGet(field_dst, localDE=0, farrayPtr=var3d_dstR8, rc=rc)
         !      var3d_dstR8 = var3d_srcR8
         !   endif
         !endif
      enddo

      deallocate(fields_names)

      end subroutine copy_esmfBundle
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine addTracerToBundle3D(bundle, grid, fieldName, type_kind, dummy_array)
!
! !INPUT PARAMETERS:
      type (ESMF_Grid),           intent(in) :: grid
      type (ESMF_TypeKind_Flag),  intent(in) :: type_kind
      character(len=*),           intent(in) :: fieldName
      real,                       intent(in) :: dummy_array(1,1,1)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Adds tracers to a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=dp), pointer :: PTR8(:,:,:)
      real(kind=sp), pointer :: PTR4(:,:,:)
      type(ESMF_ArraySpec)  :: arrayspec
      type (ESMF_Field)     :: field
      character(len=MAXSTR) :: IAm='addTracerToBundle3D'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_ArraySpecSet(arrayspec, rank=3, typekind=type_kind)

      ! Create the field and have it create the array internally
      field = ESMF_FieldCreate(grid, arrayspec, &
                                    totalLWidth=(/0,0,0/), totalUWidth=(/0,0,0/), &
                                    name = TRIM(fieldName), rc=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldCreate failed for '//TRIM(fieldName))

!      if (type_kind == ESMF_TYPEKIND_R4) then
!         call ESMF_FieldGet(field, farrayPtr=PTR4, rc=rc)
!         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))
!         PTR4 = 0.0
!      elseif (type_kind == ESMF_TYPEKIND_R8) then
!         call ESMF_FieldGet(field, farrayPtr=PTR8, rc=rc)
!         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))
!         PTR8 = 0.0d0
!      endif

      call ESMF_FieldBundleAdd ( bundle, (/field/), rc=rc )
      call LIS_verify(rc, &
                    TRIM(IAm)//': ESMF_FieldBundleAdd failed for '//TRIM(fieldName))

      return

      end subroutine addTracerToBundle3D
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName3D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=dp), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByName3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByName3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName3D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=sp), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByName3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByName3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex3D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=dp), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByIndex3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex3D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=sp), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByIndex3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByName3D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByName3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      return

      end subroutine getPointerFromBundle_ByName3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByName3D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByName3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      return

      end subroutine getPointerFromBundle_ByName3D_r4
!!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByIndex3D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByIndex3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromBundle_ByIndex3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByIndex3D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByIndex3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromBundle_ByIndex3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine addTracerToBundle2D(bundle, grid, fieldName, type_kind, dummy_array)
!
! !INPUT PARAMETERS:
      type (ESMF_Grid),           intent(in) :: grid
      type (ESMF_TypeKind_Flag),  intent(in) :: type_kind
      character(len=*),           intent(in) :: fieldName
      real,                       intent(in) :: dummy_array(1,1)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Adds tracers to a bundle.
!
! !LOCAL VARIABLES:
      real(kind=dp), pointer        :: PTR8(:,:)
      real(kind=sp), pointer        :: PTR4(:,:)
      type(ESMF_ArraySpec)         :: arrayspec
      integer                      :: rc, status
      type (ESMF_Field)            :: field
      character(len=MAXSTR)        :: IAm='addTracerToBundle2D'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=type_kind)

      ! Create the field and have it create the array internally
      field = ESMF_FieldCreate(grid, arrayspec, &
                                    totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                                    name = TRIM(fieldName), rc=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldCreate failed for '//TRIM(fieldName))

!      if (type_kind == ESMF_TYPEKIND_R4) then
!         call ESMF_FieldGet(field, farrayPtr=PTR4, rc=rc)
!         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))
!         PTR4 = 0.0
!      elseif (type_kind == ESMF_TYPEKIND_R8) then
!         call ESMF_FieldGet(field, farrayPtr=PTR8, rc=rc)
!         call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))
!         PTR8 = 0.0d0
!      endif

      call ESMF_FieldBundleAdd ( bundle, (/field/), rc=rc )
      call LIS_verify(rc, &
                   TRIM(IAm)//': ESMF_FieldBundleAdd failed for '//TRIM(fieldName))

      return

      end subroutine addTracerToBundle2D
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName2D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=dp), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByName2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      !call ESMF_FieldGet  (field, array=array,                 RC=rc)
      !call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      !call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=rc)
      !call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, farrayptr=ptr2D,   RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByName2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName2D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=sp), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByName2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      !call ESMF_FieldGet  (field, array=array,                 RC=rc)
      !call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      !call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=rc)
      !call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, localDE=0, farrayPtr=ptr2D,   RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByName2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex2D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=dp), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByIndex2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex2D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(kind=sp), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='updateTracerToBundle_ByIndex2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByName2D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByName2D_r8'
!
!EOP
!---------------------------------------------------------------------------
!BOC
     call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      return

      end subroutine getPointerFromBundle_ByName2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByName2D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByName2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, TRIM(fieldName), field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed for '//TRIM(fieldName))

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed for '//TRIM(fieldName))

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed for '//TRIM(fieldName))

      return

      end subroutine getPointerFromBundle_ByName2D_r4
!!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByIndex2D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByIndex2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromBundle_ByIndex2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromBundle_ByIndex2D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromBundle_ByIndex2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldBundleGet failed')

      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromBundle_ByIndex2D_r4
!EOC
!------------------------------------------------------------------------------
      function createEmptyField(field_name, grid, staggerloc) result(field)
        character(len=*),       intent(IN   ) :: field_name
        type (ESMF_STAGGERLOC), intent(in   ) :: staggerloc
        type (ESMF_Grid),       intent(INout) :: grid
        type (ESMF_Field)                     :: FIELD
  
        character(len=ESMF_MAXSTR),parameter  :: IAm=" createEmptyField"
        integer                               :: rc
  
        field = ESMF_FieldEmptyCreate(name=field_name, rc=rc)
        call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldEmptyCreate failed')
  
        call ESMF_FieldEmptySet(field, &
                  grid       = grid, &
                  staggerloc = staggerloc, rc = rc)
        call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldEmptySet failed')
  
      end function createEmptyField
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromField_r8_2D(field, PTR)
!
! !OUTPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Field), intent(inOut) :: field
!
! !DESCRIPTION:
! Obtains data pointer from a field.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromField_r8_2D'
!
!EOP
!---------------------------------------------------------------------------
!BOC
      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromField_r8_2D
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromField_r4_2D(field, PTR)
!
! !OUTPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Field), intent(inOut) :: field
!
! !DESCRIPTION:
! Obtains data pointer from a field.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromField_r4_2D'
!
!EOP
!---------------------------------------------------------------------------
!BOC
      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromField_r4_2D
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromField_r8_3D(field, PTR)
!
! !OUTPUT PARAMETERS:
      real(kind=dp), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Field), intent(inOut) :: field
!
! !DESCRIPTION:
! Obtains data pointer from a field.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromField_r8_3D'
!
!EOP
!---------------------------------------------------------------------------
!BOC
      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromField_r8_3D
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine getPointerFromField_r4_3D(field, PTR)
!
! !OUTPUT PARAMETERS:
      real(kind=sp), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Field), intent(inOut) :: field
!
! !DESCRIPTION:
! Obtains data pointer from a field.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_ARRAY)           :: array
      character(len=MAXSTR) :: IAm='getPointerFromField_r4_3D'
!
!EOP
!---------------------------------------------------------------------------
!BOC
      call ESMF_FieldGet  (field, array=array,                 RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_FieldGet failed')

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=rc)
      call LIS_verify(rc, TRIM(IAm)//': ESMF_ArrayGet failed')

      return

      end subroutine getPointerFromField_r4_3D
!EOC
!---------------------------------------------------------------------------

      end module LIS_field_bundleMod
