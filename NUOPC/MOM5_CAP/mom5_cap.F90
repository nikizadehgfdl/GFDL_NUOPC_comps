!--------------- MOM5 Ocean solo model -----------------
! This is the MOM5 ocean model cap
!
! Author:  Fei.Liu@gmail.com
!
! 5/10/13
!

module mom_cap_mod
  use constants_mod,            only: constants_init
  use data_override_mod,        only: data_override_init, data_override
  use diag_manager_mod,         only: diag_manager_init, diag_manager_end
  use field_manager_mod,        only: field_manager_init, field_manager_end
  use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
  use fms_mod,                  only: close_file, file_exist, uppercase
  use fms_io_mod,               only: fms_io_exit
  use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain, mpp_get_compute_domains
  use mpp_domains_mod,          only: mpp_get_ntile_count, mpp_get_pelist, mpp_get_global_domain
  use mpp_domains_mod,          only: mpp_get_domain_npes, mpp_global_field
  use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
  use mpp_mod,                  only: input_nml_file, mpp_error, FATAL, NOTE, mpp_pe, mpp_npes, mpp_set_current_pelist
  use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, MAXPES
  use time_interp_external_mod, only: time_interp_external_init
  use time_manager_mod,         only: set_calendar_type, time_type, increment_date
  use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
  use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )
  use time_manager_mod,         only: operator( * ), operator( /= ), operator( > )
  use time_manager_mod,         only: date_to_string
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
  use ocean_domains_mod,        only: get_local_indices, get_global_indices, get_domain_offsets
  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end, get_ocean_grid
  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_model_mod,          only: ocean_model_data_get
  use ocean_types_mod,          only: ice_ocean_boundary_type, ocean_grid_type

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  use time_utils_mod

  implicit none
  private
  public SetServices

  type ocean_internalstate_type
    type(ocean_public_type),       pointer :: ocean_public_type_ptr
    type(ocean_state_type),        pointer :: ocean_state_type_ptr
    type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
  end type

  type ocean_internalstate_wrapper
    type(ocean_internalstate_type), pointer :: ptr
  end type

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToOcn_num = 0
  type (fld_list_type) :: fldsToOcn(fldsMax)
  integer :: fldsFrOcn_num = 0
  type (fld_list_type) :: fldsFrOcn(fldsMax)

  integer   :: import_slice = 1
  integer   :: export_slice = 1
  character(len=256) :: tmpstr
  integer   :: dbrc

  type(ESMF_Grid), save   :: mom_grid_i
  logical                 :: write_diagnostics = .true.
  logical                 :: profile_memory = .true.
  integer(ESMF_KIND_I8)   :: restart_interval

  contains
  !-----------------------------------------------------------------------
  !------------------- Solo Ocean code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(mom_cap:SetServices)'

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=ocean_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    ! No need to change clock settings
    !call ESMF_MethodAdd(gcomp, label=model_label_SetClock, &
    !  userRoutine=SetClock, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    
  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    character(len=10)                         :: value

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write_diagnostics=(trim(value)=="true")

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    profile_memory=(trim(value)/="false")

    ! Retrieve restart_interval in (seconds)
    ! A restart_interval value of 0 means no restart will be written.
    call ESMF_AttributeGet(gcomp, name="restart_interval", value=value, defaultValue="0", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    restart_interval = ESMF_UtilString2Int(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(restart_interval < 0) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="MOM5_CAP: OCN attribute: restart_interval cannot be negative.", &
        line=__LINE__, &
        file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_LogWrite('MOM5_CAP:restart_interval = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    type(ESMF_VM)                          :: vm
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT
    
    type (ocean_public_type),      pointer :: Ocean_sfc   => NULL()
    type (ocean_state_type),       pointer :: Ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate

    type(time_type)                        :: Run_len      ! length of experiment 
    type(time_type)                        :: Time        
    type(time_type)                        :: Time_restart
    type(time_type)                        :: DT
    integer                                :: DT_OCEAN
    integer                                :: isc,iec,jsc,jec
    integer                                :: dt_cpld  = 86400
    integer                                :: year=0, month=0, day=0, hour=0, minute=0, second=0
    integer                                :: mpi_comm_mom

    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut

    integer                                :: npet, npet_x, npet_y
    character(len=*),parameter  :: subname='(mom_cap:InitializeAdvertise)'

    rc = ESMF_SUCCESS

    allocate(Ice_ocean_boundary)
    !allocate(Ocean_state) ! ocean_model_init allocate this pointer
    allocate(Ocean_sfc)
    allocate(ocean_internalstate%ptr)
    ocean_internalstate%ptr%ice_ocean_boundary_type_ptr => Ice_ocean_boundary
    ocean_internalstate%ptr%ocean_public_type_ptr       => Ocean_sfc
    ocean_internalstate%ptr%ocean_state_type_ptr        => Ocean_state

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(VM, mpiCommunicator=mpi_comm_mom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HOUR,    M =MINUTE,    S =SECOND,  &
                                        RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fms_init(mpi_comm_mom)
    call constants_init
    call field_manager_init
    call diag_manager_init
    call set_calendar_type (JULIAN                )
    ! this ocean connector will be driven at set interval
    dt_cpld = DT_OCEAN
    DT = set_time (DT_OCEAN, 0)         
    Time = set_date (YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)

    call ocean_model_init(Ocean_sfc, Ocean_state, Time, Time)
    call data_override_init(Ocean_domain_in = Ocean_sfc%domain)
    call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)

    allocate ( Ice_ocean_boundary% u_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% v_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% t_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% q_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% salt_flux (isc:iec,jsc:jec),       &
               Ice_ocean_boundary% lw_flux (isc:iec,jsc:jec),         &
               Ice_ocean_boundary% sw_flux_vis_dir (isc:iec,jsc:jec), &
               Ice_ocean_boundary% sw_flux_vis_dif (isc:iec,jsc:jec), &
               Ice_ocean_boundary% sw_flux_nir_dir (isc:iec,jsc:jec), &
               Ice_ocean_boundary% sw_flux_nir_dif (isc:iec,jsc:jec), &
               Ice_ocean_boundary% lprec (isc:iec,jsc:jec),           &
               Ice_ocean_boundary% fprec (isc:iec,jsc:jec),           &
               Ice_ocean_boundary% runoff (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% calving (isc:iec,jsc:jec),         &
               Ice_ocean_boundary% runoff_hflx (isc:iec,jsc:jec),     &
               Ice_ocean_boundary% calving_hflx (isc:iec,jsc:jec),    &
               Ice_ocean_boundary% mi (isc:iec,jsc:jec),              &
               Ice_ocean_boundary% p (isc:iec,jsc:jec))

    Ice_ocean_boundary%u_flux          = 0.0
    Ice_ocean_boundary%v_flux          = 0.0
    Ice_ocean_boundary%t_flux          = 0.0
    Ice_ocean_boundary%q_flux          = 0.0
    Ice_ocean_boundary%salt_flux       = 0.0
    Ice_ocean_boundary%lw_flux         = 0.0
    Ice_ocean_boundary%sw_flux_vis_dir = 0.0
    Ice_ocean_boundary%sw_flux_vis_dif = 0.0
    Ice_ocean_boundary%sw_flux_nir_dir = 0.0
    Ice_ocean_boundary%sw_flux_nir_dif = 0.0
    Ice_ocean_boundary%lprec           = 0.0
    Ice_ocean_boundary%fprec           = 0.0
    Ice_ocean_boundary%runoff          = 0.0
    Ice_ocean_boundary%calving         = 0.0
    Ice_ocean_boundary%runoff_hflx     = 0.0
    Ice_ocean_boundary%calving_hflx    = 0.0
    Ice_ocean_boundary%mi              = 0.0
    Ice_ocean_boundary%p               = 0.0

    call external_coupler_sbc_init(Ocean_sfc%domain, dt_cpld, Run_len)

    ocean_internalstate%ptr%ocean_state_type_ptr => Ocean_state
    call ESMF_GridCompSetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call MOM_FieldsSetup(ice_ocean_boundary, ocean_sfc)

    call MOM_AdvertiseFields(importState, fldsToOcn_num, fldsToOcn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call MOM_AdvertiseFields(exportState, fldsFrOcn_num, fldsFrOcn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(*,*) '----- MOM5 initialization phase Advertise completed'

  end subroutine InitializeAdvertise
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut
    type(ESMF_DeLayout)                    :: delayout
    type(ESMF_Distgrid)                    :: Distgrid
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    type (ocean_public_type),      pointer :: Ocean_sfc   => NULL()
    type (ocean_state_type),       pointer :: Ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    integer                                :: npet, ntiles
    integer                                :: nxg, nyg, cnt
    integer                                :: isc,iec,jsc,jec
    integer, allocatable                   :: xb(:),xe(:),yb(:),ye(:),pe(:)
    integer, allocatable                   :: deBlockList(:,:,:), &
                                              petMap(:),deLabelList(:), &
                                              indexList(:)
    integer                                :: ioff, joff
    integer                                :: i, j, n, i1, j1, n1
    integer                                :: lbnd1,ubnd1,lbnd2,ubnd2
    integer                                :: lbnd3,ubnd3,lbnd4,ubnd4
    integer                                :: nblocks_tot
    logical                                :: found
    real(ESMF_KIND_R8), allocatable        :: ofld(:,:), gfld(:,:)
    real(ESMF_KIND_R8), pointer            :: t_surf(:,:)
    integer(ESMF_KIND_I4), pointer         :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_area(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_xcen(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_ycen(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_xcor(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_ycor(:,:)
    type(ESMF_Field)                       :: field_t_surf
    character(len=*),parameter  :: subname='(mom_cap:InitializeRealize)'
    
    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, petCount=npet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !---------------------------------
    ! global mom grid size
    !---------------------------------

    call mpp_get_global_domain(Ocean_sfc%domain, xsize=nxg, ysize=nyg)
    write(tmpstr,'(a,2i6)') subname//' nxg,nyg = ',nxg,nyg
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  

    !---------------------------------
    ! number of tiles per PET, assumed to be 1, and number of pes (tiles) total
    !---------------------------------

    ntiles=mpp_get_ntile_count(Ocean_sfc%domain) ! this is tiles on this pe
    if (ntiles /= 1) then
      rc = ESMF_FAILURE
      call ESMF_LogWrite(subname//' ntiles must be 1', ESMF_LOGMSG_ERROR, rc=dbrc)  
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    ntiles=mpp_get_domain_npes(Ocean_sfc%domain)
    write(tmpstr,'(a,1i6)') subname//' ntiles = ',ntiles
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  

    !---------------------------------
    ! get start and end indices of each tile and their PET
    !---------------------------------

    allocate(xb(ntiles),xe(ntiles),yb(ntiles),ye(ntiles),pe(ntiles))
    call mpp_get_compute_domains(Ocean_sfc%domain, xbegin=xb, xend=xe, ybegin=yb, yend=ye)
    call mpp_get_pelist(Ocean_sfc%domain, pe)
    do n = 1,ntiles
      write(tmpstr,'(a,6i6)') subname//' tiles ',n,pe(n),xb(n),xe(n),yb(n),ye(n)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  
    enddo

    !---------------------------------
    ! create delayout and distgrid
    !---------------------------------

    allocate(deBlockList(2,2,ntiles))
    allocate(petMap(ntiles))
    allocate(deLabelList(ntiles))

    do n = 1, ntiles
       deLabelList(n) = n
       deBlockList(1,1,n) = xb(n)
       deBlockList(1,2,n) = xe(n)
       deBlockList(2,1,n) = yb(n)
       deBlockList(2,2,n) = ye(n)
       petMap(n) = pe(n)
       ! write(tmpstr,'(a,3i8)') subname//' iglo = ',n,deBlockList(1,1,n),deBlockList(1,2,n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       ! write(tmpstr,'(a,3i8)') subname//' jglo = ',n,deBlockList(2,1,n),deBlockList(2,2,n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       ! write(tmpstr,'(a,2i8)') subname//' pe  = ',n,petMap(n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       !--- assume a tile with starting index of 1 has an equivalent wraparound tile on the other side
    enddo

    delayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    
    allocate(connectionList(2))
    ! bipolar boundary condition at top row: nyg
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
      tileIndexB=1, positionVector=(/nxg+1, 2*nyg+1/), &
      orientationVector=(/-1, -2/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! periodic boundary condition along first dimension
    call ESMF_DistGridConnectionSet(connectionList(2), tileIndexA=1, &
      tileIndexB=1, positionVector=(/nxg, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nxg,nyg/), &
!        indexflag = ESMF_INDEX_DELOCAL, &
        deBlockList=deBlockList, &
!        deLabelList=deLabelList, &
        delayout=delayout, &
        connectionList=connectionList, &
        rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    deallocate(xb,xe,yb,ye,pe)
    deallocate(connectionList)
    deallocate(deLabelList)
    deallocate(deBlockList)
    deallocate(petMap)

    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
    allocate(indexList(cnt))
    write(tmpstr,'(a,i8)') subname//' distgrid cnt= ',cnt
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
    write(tmpstr,'(a,4i8)') subname//' distgrid list= ',indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    deallocate(IndexList)

    !---------------------------------
    ! create grid
    !---------------------------------

    gridIn = ESMF_GridCreate(distgrid=distgrid, &
       gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
       coordSys = ESMF_COORDSYS_SPH_DEG, &
       rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    mom_grid_i = gridIn

    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGetCoord(gridIn, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_xcen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_ycen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGetCoord(gridIn, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=dataPtr_xcor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=dataPtr_ycor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_mask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_AREA, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_area, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !---------------------------------
    ! load up area, mask, center and corner values
    ! area, mask, and centers should be same size in mom and esmf grid
    ! corner points may not be, need to offset corner points by 1 in i and j
    !   for esmf and also need to "make up" j=1 values.  use wraparound in i
    !---------------------------------

    call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)

    lbnd1 = lbound(dataPtr_area,1)
    ubnd1 = ubound(dataPtr_area,1)
    lbnd2 = lbound(dataPtr_area,2)
    ubnd2 = ubound(dataPtr_area,2)

    lbnd3 = lbound(dataPtr_xcor,1)
    ubnd3 = ubound(dataPtr_xcor,1)
    lbnd4 = lbound(dataPtr_xcor,2)
    ubnd4 = ubound(dataPtr_xcor,2)

    write(tmpstr,*) subname//' iscjsc = ',isc,iec,jsc,jec
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' lbub12 = ',lbnd1,ubnd1,lbnd2,ubnd2
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' lbub34 = ',lbnd3,ubnd3,lbnd4,ubnd4
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    if (iec-isc /= ubnd1-lbnd1 .or. jec-jsc /= ubnd2-lbnd2) then
       rc=ESMF_FAILURE
       call ESMF_LogWrite(subname//' fld and grid not same size', ESMF_LOGMSG_ERROR, rc=dbrc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    allocate(ofld(isc:iec,jsc:jec))
    allocate(gfld(nxg,nyg))

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld mask = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld mask = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'area', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld area = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld area = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_area(i,j) = ofld(i1,j1)
    enddo
    enddo

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'tlon', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld xt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld xt = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_xcen(i,j) = ofld(i1,j1)
       dataPtr_xcen(i,j) = mod(dataPtr_xcen(i,j)+720.0_ESMF_KIND_R8,360.0_ESMF_KIND_R8)
    enddo
    enddo

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'tlat', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld yt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld yt = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_ycen(i,j) = ofld(i1,j1)
    enddo
    enddo

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'ulon', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld xu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld xu = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd4, ubnd4
    do i = lbnd3, ubnd3
       j1 = j - lbnd4 + jsc - 1
       i1 = mod(i - lbnd3 + isc - 2 + nxg, nxg) + 1
       if (j1 == 0) then
          dataPtr_xcor(i,j) = 2*gfld(i1,1) - gfld(i1,2)
!          if (dataPtr_xcor(i,j)-dataPtr_xcen(i,j) > 180.) dataPtr_xcor(i,j) = dataPtr_xcor(i,j) - 360.
!          if (dataPtr_xcor(i,j)-dataPtr_xcen(i,j) < 180.) dataPtr_xcor(i,j) = dataPtr_xcor(i,j) + 360.
       elseif (j1 >= 1 .and. j1 <= nyg) then
          dataPtr_xcor(i,j) = gfld(i1,j1)
       else
          rc=ESMF_FAILURE
          call ESMF_LogWrite(subname//' error in xu j1', ESMF_LOGMSG_ERROR, rc=dbrc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       endif
       dataPtr_xcor(i,j) = mod(dataPtr_xcor(i,j)+720.0_ESMF_KIND_R8,360.0_ESMF_KIND_R8)
       ! write(tmpstr,*) subname//' ijfld xu = ',i,i1,j,j1,dataPtr_xcor(i,j)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    enddo

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'ulat', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld yu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld yu = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd4, ubnd4
    do i = lbnd3, ubnd3
       j1 = j - lbnd4 + jsc - 1
       i1 = mod(i - lbnd3 + isc - 2 + nxg, nxg) + 1
       if (j1 == 0) then
          dataPtr_ycor(i,j) = 2*gfld(i1,1) - gfld(i1,2)
       elseif (j1 >= 1 .and. j1 <= nyg) then
          dataPtr_ycor(i,j) = gfld(i1,j1)
       else
          rc=ESMF_FAILURE
          call ESMF_LogWrite(subname//' error in xu j1', ESMF_LOGMSG_ERROR, rc=dbrc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       endif
       ! write(tmpstr,*) subname//' ijfld yu = ',i,i1,j,j1,dataPtr_ycor(i,j)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    enddo

    write(tmpstr,*) subname//' mask = ',minval(dataPtr_mask),maxval(dataPtr_mask)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' area = ',minval(dataPtr_area),maxval(dataPtr_area)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' xcen = ',minval(dataPtr_xcen),maxval(dataPtr_xcen)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' ycen = ',minval(dataPtr_ycen),maxval(dataPtr_ycen)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' xcor = ',minval(dataPtr_xcor),maxval(dataPtr_xcor)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' ycor = ',minval(dataPtr_ycor),maxval(dataPtr_ycor)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    deallocate(gfld)

    gridOut = gridIn ! for now out same as in

    !---------------------------------
    ! realize fields on grid
    !---------------------------------

    call MOM_RealizeFields(importState, gridIn , fldsToOcn_num, fldsToOcn, "Ocn import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call MOM_RealizeFields(exportState, gridOut, fldsFrOcn_num, fldsFrOcn, "Ocn export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(exportState, itemName='sea_surface_temperature', field=field_t_surf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldGet(field_t_surf, localDe=0, farrayPtr=t_surf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)

    lbnd1 = lbound(t_surf,1)
    ubnd1 = ubound(t_surf,1)
    lbnd2 = lbound(t_surf,2)
    ubnd2 = ubound(t_surf,2)

    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       if (ofld(i1,j1) == 0.) t_surf(i,j) = 0.0
    enddo
    enddo

    deallocate(ofld)

    call NUOPC_Write(exportState, fileNamePrefix='init_field_ocn_export_', &
      timeslice=1, relaxedFlag=.true., rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(*,*) '----- MOM5 initialization phase Realize completed'

  end subroutine InitializeRealize
  
  !-----------------------------------------------------------------------------

  ! Ocean model uses same clock as parent gridComp
  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep, timestep
    character(len=*),parameter  :: subname='(mom_cap:SetClock)'

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    call ESMF_TimeIntervalSet(stabilityTimeStep, m=120, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine SetClock

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
    character(len=64)                      :: timestamp

    type (ocean_public_type),      pointer :: Ocean_sfc          => NULL()
    type (ocean_state_type),       pointer :: Ocean_state        => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate

    ! define some time types 
    type(time_type)                        :: Time        
    type(time_type)                        :: Time_step_coupled
    type(time_type)                        :: Time_restart_current

    integer :: dth, dtm, dts, dt_cpld  = 86400
    integer :: isc,iec,jsc,jec,lbnd1,ubnd1,lbnd2,ubnd2
    integer :: i,j,i1,j1
    real(ESMF_KIND_R8), allocatable        :: ofld(:,:), ocz(:,:), ocm(:,:)
    real(ESMF_KIND_R8), allocatable        :: mmmf(:,:), mzmf(:,:)
    integer :: nc
    real(ESMF_KIND_R8), pointer :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_mmmf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_mzmf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocz(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocm(:,:) 
!    real(ESMF_KIND_R8), pointer :: dataPtr_oci(:,:)
!    real(ESMF_KIND_R8), pointer :: dataPtr_ocj(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_frazil(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sensi(:,:)
    type(ocean_grid_type), pointer :: Ocean_grid
    character(240)              :: msgString
    character(len=*),parameter  :: subname='(mom_cap:ModelAdvance)'

    rc = ESMF_SUCCESS
    if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM5 Model_ADVANCE: ")
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
      timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", &
      unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeIntervalGet(timeStep, h=dth, m=dtm, s=dts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Time = esmf2fms_time(currTime)
    Time_step_coupled = esmf2fms_time(timeStep)
    dt_cpld = dth*3600+dtm*60+dts

    call ice_ocn_bnd_from_data(Ice_ocean_boundary, Time, Time_step_coupled)

    call external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

    if(write_diagnostics) then
      call NUOPC_Write(importState, fileNamePrefix='field_ocn_import_', &
        timeslice=import_slice, relaxedFlag=.true., rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      import_slice = import_slice + 1
    endif

    ! rotate the lat/lon wind vector (CW) onto local tripolar coordinate system

    call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)

    call State_getFldPtr(exportState,'ocean_mask',dataPtr_mask,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    lbnd1 = lbound(dataPtr_mask,1)
    ubnd1 = ubound(dataPtr_mask,1)
    lbnd2 = lbound(dataPtr_mask,2)
    ubnd2 = ubound(dataPtr_mask,2)

    call get_ocean_grid(Ocean_grid)

    call State_getFldPtr(importState,'mean_zonal_moment_flx',dataPtr_mzmf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_merid_moment_flx',dataPtr_mmmf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_evap_rate',dataPtr_evap,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mean_sensi_heat_flx',dataPtr_sensi,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    dataPtr_evap = - dataPtr_evap
    dataPtr_sensi = - dataPtr_sensi

    allocate(mzmf(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(mmmf(lbnd1:ubnd1,lbnd2:ubnd2))
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc  ! work around local vs global indexing
        i1 = i - lbnd1 + isc
        mzmf(i,j) = Ocean_grid%cos_rot(i1,j1)*dataPtr_mzmf(i,j) &
                  + Ocean_grid%sin_rot(i1,j1)*dataPtr_mmmf(i,j)
        mmmf(i,j) = Ocean_grid%cos_rot(i1,j1)*dataPtr_mmmf(i,j) &
                  - Ocean_grid%sin_rot(i1,j1)*dataPtr_mzmf(i,j)
      enddo
    enddo
    dataPtr_mzmf = mzmf
    dataPtr_mmmf = mmmf
    deallocate(mzmf, mmmf)

    !Optionally write restart files when currTime-startTime is integer multiples of restart_interval
    if(restart_interval > 0 ) then
      time_elapsed = currTime - startTime
      call ESMF_TimeIntervalGet(time_elapsed, s_i8=time_elapsed_sec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
      n_interval = time_elapsed_sec / restart_interval
      if((n_interval .gt. 0) .and. (n_interval*restart_interval == time_elapsed_sec)) then
          time_restart_current = esmf2fms_time(currTime)
          timestamp = date_to_string(time_restart_current)
          call ESMF_LogWrite("MOM5: Writing restart at "//trim(timestamp), ESMF_LOGMSG_INFO, rc=dbrc)
          write(*,*) 'calling ocean_model_restart'
          call ocean_model_restart(Ocean_state, timestamp)
      endif
    endif

    if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM5 update_ocean_model: ")
    call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, Time, Time_step_coupled)
    if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM5 update_ocean_model: ")

    allocate(ofld(isc:iec,jsc:jec))

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo
    deallocate(ofld)

    ! Now rotate ocn current from tripolar grid back to lat/lon grid (CCW)
    allocate(ocz(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(ocm(lbnd1:ubnd1,lbnd2:ubnd2))

    call State_getFldPtr(exportState,'ocn_current_zonal',dataPtr_ocz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'ocn_current_merid',dataPtr_ocm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
!    call State_getFldPtr(exportState,'ocn_current_idir',dataPtr_oci,rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
!    call State_getFldPtr(exportState,'ocn_current_jdir',dataPtr_ocj,rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'freezing_melting_potential',dataPtr_frazil,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    dataPtr_frazil = dataPtr_frazil/dt_cpld !convert from J/m^2 to W/m^2 for CICE coupling

    ocz = dataPtr_ocz
    ocm = dataPtr_ocm
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc  ! work around local vs global indexing
        i1 = i - lbnd1 + isc
        dataPtr_ocz(i,j) = Ocean_grid%cos_rot(i1,j1)*ocz(i,j) &
                         - Ocean_grid%sin_rot(i1,j1)*ocm(i,j)
        dataPtr_ocm(i,j) = Ocean_grid%cos_rot(i1,j1)*ocm(i,j) &
                         + Ocean_grid%sin_rot(i1,j1)*ocz(i,j)
!        dataPtr_oci(i,j) = ocz(i,j)
!        dataPtr_ocj(i,j) = ocm(i,j)
      enddo
    enddo
    deallocate(ocz, ocm)

    if(write_diagnostics) then
      call NUOPC_Write(exportState, fileNamePrefix='field_ocn_export_', &
        timeslice=export_slice, relaxedFlag=.true., rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      export_slice = export_slice + 1
    endif

    call external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

    !write(*,*) 'MOM: --- run phase called ---'
    call dumpMomInternal(mom_grid_i, import_slice, "mean_zonal_moment_flx", "will provide", Ice_ocean_boundary%u_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_merid_moment_flx", "will provide", Ice_ocean_boundary%v_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_sensi_heat_flx"  , "will provide", Ice_ocean_boundary%t_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_evap_rate"       , "will provide", Ice_ocean_boundary%q_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_salt_rate"       , "will provide", Ice_ocean_boundary%salt_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_lw_flx"      , "will provide", Ice_ocean_boundary%lw_flux  )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_vis_dir_flx", "will provide", Ice_ocean_boundary%sw_flux_vis_dir)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_vis_dif_flx", "will provide", Ice_ocean_boundary%sw_flux_vis_dif)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_ir_dir_flx" , "will provide", Ice_ocean_boundary%sw_flux_nir_dir)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_ir_dif_flx" , "will provide", Ice_ocean_boundary%sw_flux_nir_dif)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_prec_rate"       , "will provide", Ice_ocean_boundary%lprec  )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_fprec_rate"      , "will provide", Ice_ocean_boundary%fprec  )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_runoff_rate"     , "will provide", Ice_ocean_boundary%runoff )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_calving_rate"    , "will provide", Ice_ocean_boundary%calving)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_runoff_heat_flx" , "will provide", Ice_ocean_boundary%runoff_hflx )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_calving_heat_flx", "will provide", Ice_ocean_boundary%calving_hflx)
    call dumpMomInternal(mom_grid_i, import_slice, "inst_pres_height_surface" , "will provide", Ice_ocean_boundary%p )
    call dumpMomInternal(mom_grid_i, import_slice, "mass_of_overlying_sea_ice", "will provide", Ice_ocean_boundary%mi)

!--------- export fields -------------

    call dumpMomInternal(mom_grid_i, export_slice, "ocean_mask", "will provide", dataPtr_mask)
    call dumpMomInternal(mom_grid_i, export_slice, "sea_surface_temperature", "will provide", Ocean_sfc%t_surf)
    call dumpMomInternal(mom_grid_i, export_slice, "s_surf"    , "will provide", Ocean_sfc%s_surf )
    call dumpMomInternal(mom_grid_i, export_slice, "ocn_current_zonal", "will provide", Ocean_sfc%u_surf )
    call dumpMomInternal(mom_grid_i, export_slice, "ocn_current_merid", "will provide", Ocean_sfc%v_surf )
!    call dumpMomInternal(mom_grid_i, export_slice, "ocn_current_idir", "will provide", dataPtr_oci )
!    call dumpMomInternal(mom_grid_i, export_slice, "ocn_current_jdir", "will provide", dataPtr_ocj )
    call dumpMomInternal(mom_grid_i, export_slice, "sea_lev"   , "will provide", Ocean_sfc%sea_lev)

    if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM5 Model_ADVANCE: ")
  end subroutine ModelAdvance

  subroutine ocean_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type (ocean_public_type),      pointer :: Ocean_sfc          
    type (ocean_state_type),       pointer :: Ocean_state
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    type(TIME_TYPE)                        :: Time        
    type(ESMF_Clock)                       :: clock
    type(ESMF_Time)                        :: currTime
    character(len=64)                      :: timestamp
    character(len=*),parameter  :: subname='(mom_cap:ocean_model_finalize)'

    write(*,*) 'MOM: --- finalize called ---'
    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time = esmf2fms_time(currTime)

    call ocean_model_end (Ocean_sfc, Ocean_State, Time)
    call diag_manager_end(Time )
    call field_manager_end

    call fms_io_exit
    call fms_end

    write(*,*) 'MOM: --- completed ---'

  end subroutine ocean_model_finalize

!====================================================================
! get forcing data from data_overide 
  subroutine ice_ocn_bnd_from_data(x, Time, Time_step_coupled)

      type (ice_ocean_boundary_type) :: x
      type(Time_type), intent(in)    :: Time, Time_step_coupled

      type(Time_type)                :: Time_next
      character(len=*),parameter  :: subname='(mom_cap:ice_ocn_bnd_from_data)'

      Time_next = Time + Time_step_coupled

      !call data_override('OCN', 't_flux',          x%t_flux         , Time_next)
      !call data_override('OCN', 'u_flux',          x%u_flux         , Time_next)
      !call data_override('OCN', 'v_flux',          x%v_flux         , Time_next)
      !call data_override('OCN', 'q_flux',          x%q_flux         , Time_next)
      !call data_override('OCN', 'salt_flux',       x%salt_flux      , Time_next)
      !call data_override('OCN', 'lw_flux',         x%lw_flux        , Time_next)
      !call data_override('OCN', 'sw_flux_vis_dir', x%sw_flux_vis_dir, Time_next)
      !call data_override('OCN', 'sw_flux_vis_dif', x%sw_flux_vis_dif, Time_next)
      !call data_override('OCN', 'sw_flux_nir_dir', x%sw_flux_nir_dir, Time_next)
      !call data_override('OCN', 'sw_flux_nir_dif', x%sw_flux_nir_dif, Time_next)
      !call data_override('OCN', 'lprec',           x%lprec          , Time_next)
      !call data_override('OCN', 'fprec',           x%fprec          , Time_next)
      !call data_override('OCN', 'runoff',          x%runoff         , Time_next)
      !call data_override('OCN', 'calving',         x%calving        , Time_next)
      !call data_override('OCN', 'p',               x%p              , Time_next)
            
  end subroutine ice_ocn_bnd_from_data


!-----------------------------------------------------------------------------------------
! 
! Subroutines  for enabling coupling to external programs through a third party coupler
! such as OASIS/PRISM.
! If no external coupler then these will mostly be dummy routines.
! These routines can also serve as spots to call other user defined routines
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

! Dummy subroutines.

  subroutine external_coupler_mpi_init(mom_local_communicator, external_initialization)
  implicit none
  integer, intent(out) :: mom_local_communicator
  logical, intent(out) :: external_initialization
  external_initialization = .false.
  mom_local_communicator = -100         ! Is there mpp_undefined parameter corresponding to MPI_UNDEFINED?
                                        ! probably wouldn't need logical flag.
  return
  end subroutine external_coupler_mpi_init

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_sbc_init(Dom, dt_cpld, Run_len)
  implicit none
  type(domain2d) :: Dom
  integer :: dt_cpld
  type(time_type) :: Run_len
  return
  end  subroutine external_coupler_sbc_init

  subroutine external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
  implicit none
  type (ice_ocean_boundary_type), intent(INOUT) :: Ice_ocean_boundary
  type (ocean_public_type) , intent(INOUT)        :: Ocean_sfc
  integer , intent(IN)                       :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_before


  subroutine external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
  type (ice_ocean_boundary_type) :: Ice_ocean_boundary
  type (ocean_public_type)         :: Ocean_sfc
  integer                        :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_after

  subroutine external_coupler_restart( dt_cpld, num_cpld_calls )
  implicit none
  integer, intent(in)               :: dt_cpld, num_cpld_calls
  return
  end subroutine external_coupler_restart

  subroutine external_coupler_exit
  return
  end subroutine external_coupler_exit

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_mpi_exit(mom_local_communicator, external_initialization)
  implicit none
  integer, intent(in) :: mom_local_communicator
  logical, intent(in) :: external_initialization
  return
  end subroutine external_coupler_mpi_exit
!-----------------------------------------------------------------------------------------
    subroutine writeSliceFields(state, filename_prefix, slice, rc)
      type(ESMF_State)                :: state
      character(len=*)                :: filename_prefix
      integer                         :: slice
      integer, intent(out), optional  :: rc

      integer                         :: n, nfields
      type(ESMF_Field)                :: field
      type(ESMF_StateItem_Flag)       :: itemType
      character(len=40)               :: fileName
      character(len=64),allocatable   :: fieldNameList(:)
      character(len=*),parameter :: subname='(mom_cap:writeSliceFields)'

      if (present(rc)) rc = ESMF_SUCCESS
      
      if (ESMF_IO_PIO_PRESENT .and. &
        (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT)) then

        call ESMF_StateGet(state, itemCount=nfields, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allocate(fieldNameList(nfields))
        call ESMF_StateGet(state, itemNameList=fieldNameList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        do n=1, size(fieldNameList)
          call ESMF_StateGet(state, itemName=fieldNameList(n), &
            itemType=itemType, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          if (itemType /= ESMF_STATEITEM_NOTFOUND) then
            ! field is available in the state
            call ESMF_StateGet(state, itemName=fieldNameList(n), field=field, &
              rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
            ! -> output to file
            write (fileName,"(A)") &
              filename_prefix//trim(fieldNameList(n))//".nc"
            call ESMF_FieldWrite(field, fileName=trim(fileName), &
              timeslice=slice, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
        enddo

        deallocate(fieldNameList)

      endif


    end subroutine writeSliceFields

  !-----------------------------------------------------------------------------

  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(mom_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr

  !-----------------------------------------------------------------------------
  subroutine MOM_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(mom_cap:MOM_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine MOM_AdvertiseFields

  !-----------------------------------------------------------------------------

  subroutine MOM_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=*),parameter  :: subname='(mom_cap:MOM_RealizeFields)'
 
    rc = ESMF_SUCCESS

    do i = 1, nfields

      if (field_defs(i)%assoc) then
        write(tmpstr, *) subname, tag, ' Field ', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2)
        call ESMF_LogWrite(tmpstr, ESMF_LOGMSG_INFO, rc=dbrc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
!          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_GLOBAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
!        call ESMF_FieldPrint(field=field, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      else
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo

  end subroutine MOM_RealizeFields

  !-----------------------------------------------------------------------------

  subroutine MOM_FieldsSetup(ice_ocean_boundary,ocean_sfc)
    type(ice_ocean_boundary_type), intent(in)   :: Ice_ocean_boundary
    type(ocean_public_type), intent(in)         :: Ocean_sfc
    character(len=*),parameter  :: subname='(mom_cap:MOM_FieldsSetup)'

  !!! fld_list_add(num, fldlist, stdname, transferOffer, data(optional), shortname(optional))

!--------- import fields -------------

! tcraig, don't point directly into mom data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_zonal_moment_flx", "will provide", data=Ice_ocean_boundary%u_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_merid_moment_flx", "will provide", data=Ice_ocean_boundary%v_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_sensi_heat_flx"  , "will provide", data=Ice_ocean_boundary%t_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_evap_rate"       , "will provide", data=Ice_ocean_boundary%q_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_salt_rate"       , "will provide", data=Ice_ocean_boundary%salt_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_lw_flx"      , "will provide", data=Ice_ocean_boundary%lw_flux  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dir_flx", "will provide", data=Ice_ocean_boundary%sw_flux_vis_dir)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dif_flx", "will provide", data=Ice_ocean_boundary%sw_flux_vis_dif)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dir_flx" , "will provide", data=Ice_ocean_boundary%sw_flux_nir_dir)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dif_flx" , "will provide", data=Ice_ocean_boundary%sw_flux_nir_dif)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_prec_rate"       , "will provide", data=Ice_ocean_boundary%lprec  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_fprec_rate"      , "will provide", data=Ice_ocean_boundary%fprec  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_rate"     , "will provide", data=Ice_ocean_boundary%runoff )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_rate"    , "will provide", data=Ice_ocean_boundary%calving)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_heat_flx" , "will provide", data=Ice_ocean_boundary%runoff_hflx )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_heat_flx", "will provide", data=Ice_ocean_boundary%calving_hflx)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "inst_pres_height_surface" , "will provide", data=Ice_ocean_boundary%p )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mass_of_overlying_sea_ice", "will provide", data=Ice_ocean_boundary%mi)

!--------- export fields -------------

    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocean_mask", "will provide")
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_temperature", "will provide", data=Ocean_sfc%t_surf)
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "s_surf"    , "will provide", data=Ocean_sfc%s_surf )
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_zonal", "will provide", data=Ocean_sfc%u_surf )
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_merid", "will provide", data=Ocean_sfc%v_surf )
!    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_idir", "will provide")
!    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_jdir", "will provide")
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_lev"   , "will provide", data=Ocean_sfc%sea_lev)
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "freezing_melting_potential"   , "will provide", data=Ocean_sfc%frazil)

  end subroutine MOM_FieldsSetup

  !-----------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, transferOffer, data, shortname)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: transferOffer
    real(ESMF_KIND_R8), dimension(:,:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(mom_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif
    fldlist(num)%transferOffer  = trim(transferOffer)
    if (present(data)) then
      fldlist(num)%assoc        = .true.
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add

  subroutine dumpMomInternal(grid, slice, stdname, nop, farray)

    type(ESMF_Grid)          :: grid
    integer, intent(in)      :: slice
    character(len=*)         :: stdname
    character(len=*)         :: nop
    real(ESMF_KIND_R8), dimension(:,:), target :: farray

    type(ESMF_Field)         :: field
    real(ESMF_KIND_R8), dimension(:,:), pointer  :: f2d
    integer                  :: rc

    if(.not. write_diagnostics) return ! nop in production mode

    field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
      indexflag=ESMF_INDEX_DELOCAL, &
      name=stdname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldGet(field, farrayPtr=f2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    f2d(:,:) = farray(:,:)

    call ESMF_FieldWrite(field, fileName='field_ocn_internal_'//trim(stdname)//'.nc', &
      timeslice=slice, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldDestroy(field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

end module mom_cap_mod
