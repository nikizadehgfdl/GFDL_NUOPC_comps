!--------------- SIS2 sea-ice model -----------------
! This is the SIS2 sea-ice model cap
!
! Author:  Niki.Zadeh
!
! 3/10/17
!

module sis2_cap_mod
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

  use ice_model_mod,            only: ice_model_init, ice_model_end
  use ice_model_mod,            only: update_ice_model_slow_up
  use ice_model_mod,            only: update_ice_model_fast
  use ice_model_mod,            only: update_ice_model_slow_dn
  use ice_model_mod,            only: ice_data_type, land_ice_boundary_type
  use ice_model_mod,            only: ocean_ice_boundary_type, atmos_ice_boundary_type
  use ice_model_mod,            only: ice_model_restart

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  implicit none
  private
  public SetServices

  type ice_internalstate_type
    type(ice_data_type),           pointer :: ice_data_type_ptr
    type(atmos_ice_boundary_type), pointer :: atmos_ice_boundary_type_ptr
    type(ocean_ice_boundary_type), pointer :: ocean_ice_boundary_type_ptr
  end type

  type ice_internalstate_wrapper
    type(ice_internalstate_type), pointer :: ptr
  end type

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:), pointer :: farrayPtr
  end type fld_list_type

!  interface fld_list_add
!        module procedure fld_list_add_2d, fld_list_add_3d
!  end interface fld_list_add

  integer,parameter :: fldsMax = 100
  integer :: fldsToIce_num = 0
  type (fld_list_type) :: fldsToIce(fldsMax)
  integer :: fldsFrIce_num = 0
  type (fld_list_type) :: fldsFrIce(fldsMax)

  integer   :: import_slice = 1
  integer   :: export_slice = 1
  character(len=256) :: tmpstr
  integer   :: dbrc

  type(ESMF_Grid), save   :: mom_grid_i
  logical                 :: write_diagnostics = .true.
  logical                 :: profile_memory = .true.
  logical                 :: ocean_solo = .true.
  integer(ESMF_KIND_I8)   :: restart_interval

  contains
  !-----------------------------------------------------------------------
  !------------------- Solo Ice code starts here -----------------------
  !-----------------------------------------------------------------------

  !> NUOPC SetService method is the only public entry point.
  !! SetServices registers all of the user-provided subroutines
  !! in the module with the NUOPC layer.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param rc return code  

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(sis2_cap:SetServices)'

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
!    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
!      specRoutine=ModelAdvance, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
!      specRoutine=ocean_model_finalize, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

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

  !> First initialize subroutine called by NUOPC.  The purpose
  !! is to set which version of the Initialize Phase Definition (IPD)
  !! to use.
  !!
  !! For this MOM cap, we are using IPDv01.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code
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

    call ESMF_AttributeGet(gcomp, name="SIS2Solo", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ocean_solo=(trim(value)=="true")

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
        msg="SIS2_CAP: SIS2 attribute: restart_interval cannot be negative.", &
        line=__LINE__, &
        file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_LogWrite('SIS2_CAP:restart_interval = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  !> Called by NUOPC to advertise import and export fields.  "Advertise"
  !! simply means that the standard names of all import and export
  !! fields are supplied.  The NUOPC layer uses these to match fields
  !! between components in the coupled system.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    type(ESMF_VM)                          :: vm
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT
    
    type(ice_data_type), pointer           :: ice_data          
    type(atmos_ice_boundary_type), pointer :: atmos_ice_boundary => NULL()
    type(ocean_ice_boundary_type), pointer :: ocean_ice_boundary => NULL()
    type(ice_internalstate_wrapper)        :: ice_internalstate

    type(time_type)                        :: Run_len      ! length of experiment 
    type(time_type)                        :: Time        
    type(time_type)                        :: Time_restart
    type(time_type)                        :: Time_step_slow, Time_step_fast
    integer                                :: DT_ICE
    integer                                :: isc,iec,jsc,jec,kd
    integer                                :: dt_cpld  = 86400
    integer                                :: year=0, month=0, day=0, hour=0, minute=0, second=0
    integer                                :: mpi_comm_mom

    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut

    integer                                :: npet, npet_x, npet_y
    character(len=*),parameter  :: subname='(sis2_cap:InitializeAdvertise)'

    rc = ESMF_SUCCESS

    allocate(atmos_ice_boundary)
    allocate(ocean_ice_boundary) 
    !allocate(ice_data) ! ice_model_init allocates this
    allocate(ice_internalstate%ptr)
    ice_internalstate%ptr%atmos_ice_boundary_type_ptr => atmos_ice_boundary
    ice_internalstate%ptr%ocean_ice_boundary_type_ptr => ocean_ice_boundary
    ice_internalstate%ptr%ice_data_type_ptr         => ice_data

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

    CALL ESMF_TimeIntervalGet(TINT, S=DT_ICE, RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !
    !GFDL/FMS subroutines,  setup calls
    !
    call fms_init(mpi_comm_mom)
    call constants_init
    call field_manager_init
    call diag_manager_init
    call set_calendar_type (JULIAN                )
    ! this ocean connector will be driven at set interval
    dt_cpld = DT_ICE
    Time_step_slow = set_time (dt_cpld ,0)
    Time_step_fast = set_time (dt_cpld,0)

    Time = set_date (YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
    allocate(ice_data)
    ice_data%pe = .true.
    ice_data%fast_ice_pe = .true.
    ice_data%slow_ice_pe = .true.

    call data_override_init
    call ice_model_init( ice_data, Time, Time, Time_step_fast, Time_step_slow )

    call data_override_init(Ice_domain_in = ice_data%domain)
    call mpp_get_compute_domain(ice_data%domain, isc, iec, jsc, jec)
    !allocate and initialize ocean_ice_boundary
    allocate ( ocean_ice_boundary%u(isc:iec,jsc:jec),          &
               ocean_ice_boundary%v(isc:iec,jsc:jec),          &
               ocean_ice_boundary%t(isc:iec,jsc:jec),          &
               ocean_ice_boundary%s(isc:iec,jsc:jec),          &
               ocean_ice_boundary%frazil(isc:iec,jsc:jec),     &
               ocean_ice_boundary%sea_level (isc:iec,jsc:jec)  )

    ocean_ice_boundary%u=0.0
    ocean_ice_boundary%v=0.0
    ocean_ice_boundary%t=273.0
    ocean_ice_boundary%s=0.0
    ocean_ice_boundary%frazil=0.0
    ocean_ice_boundary%sea_level=0.0

    !allocate and initialize atmos_ice_boundary
    kd = size(ice_data%part_size,3)
    allocate( atmos_ice_boundary%u_flux(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%v_flux(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%u_star(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%t_flux(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%q_flux(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%lw_flux(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%sw_flux_vis_dir(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%sw_flux_vis_dif(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%sw_flux_nir_dir(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%sw_flux_nir_dif(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%lprec(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%fprec(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%dhdt(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%dedt(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%drdt(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%coszen(isc:iec,jsc:jec,kd) )
    allocate( atmos_ice_boundary%p(isc:iec,jsc:jec,kd) )

    atmos_ice_boundary%u_flux=0.0
    atmos_ice_boundary%v_flux=0.0
    atmos_ice_boundary%u_star=0.0
    atmos_ice_boundary%t_flux=0.0
    atmos_ice_boundary%q_flux=0.0
    atmos_ice_boundary%lw_flux=0.0
    atmos_ice_boundary%sw_flux_vis_dir=0.0
    atmos_ice_boundary%sw_flux_vis_dif=0.0
    atmos_ice_boundary%sw_flux_nir_dir=0.0
    atmos_ice_boundary%sw_flux_nir_dif=0.0
    atmos_ice_boundary%lprec=0.0
    atmos_ice_boundary%fprec=0.0
    atmos_ice_boundary%dhdt=0.0
    atmos_ice_boundary%dedt=0.0
    atmos_ice_boundary%drdt=0.0
    atmos_ice_boundary%coszen=0.0
    atmos_ice_boundary%p=0.0

ice_data%t_surf = 274.0

    ice_internalstate%ptr%ice_data_type_ptr => ice_data
    call ESMF_GridCompSetInternalState(gcomp, ice_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call SIS_FieldsSetup(atmos_ice_boundary, ocean_ice_boundary, ice_data)

    call SIS_AdvertiseFields(importState, fldsToIce_num, fldsToIce, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call SIS_AdvertiseFields(exportState, fldsFrIce_num, fldsFrIce, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(*,*) '----- SIS2 initialization phase Advertise completed'

  end subroutine InitializeAdvertise
  
  !-----------------------------------------------------------------------------
  !> Called by NUOPC to realize import and export fields.  "Realizing" a field
  !! means that its grid has been defined and an ESMF_Field object has been
  !! created and put into the import or export State.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code
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
    
    type(ice_data_type), pointer           :: ice_data          => NULL()
    type(atmos_ice_boundary_type), pointer :: atmos_ice_boundary => NULL()
    type(ocean_ice_boundary_type), pointer :: ocean_ice_boundary => NULL()
    type(ice_internalstate_wrapper)        :: ice_internalstate

    integer                                :: npet, ntiles
    integer                                :: nxg, nyg, cnt
    integer                                :: isc,iec,jsc,jec
    integer, allocatable                   :: xb(:),xe(:),yb(:),ye(:),pe(:)
    integer, allocatable                   :: deBlockList(:,:,:), &
                                              petMap(:),deLabelList(:), &
                                              indexList(:)
    integer                                :: ioff, joff
    integer                                :: i, j, n, i1, j1, n1, icount
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
    character(len=*),parameter  :: subname='(sis2_cap:InitializeRealize)'
    
    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, ice_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ocean_ice_boundary => ice_internalstate%ptr%ocean_ice_boundary_type_ptr
    atmos_ice_boundary => ice_internalstate%ptr%atmos_ice_boundary_type_ptr
    ice_data           => ice_internalstate%ptr%ice_data_type_ptr

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

    call mpp_get_global_domain(ice_data%domain, xsize=nxg, ysize=nyg)
    write(tmpstr,'(a,2i6)') subname//' nxg,nyg = ',nxg,nyg
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  

    !---------------------------------
    ! number of tiles per PET, assumed to be 1, and number of pes (tiles) total
    !---------------------------------

    ntiles=mpp_get_ntile_count(ice_data%domain) ! this is tiles on this pe
    if (ntiles /= 1) then
      rc = ESMF_FAILURE
      call ESMF_LogWrite(subname//' ntiles must be 1', ESMF_LOGMSG_ERROR, rc=dbrc)  
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    ntiles=mpp_get_domain_npes(ice_data%domain)
    write(tmpstr,'(a,1i6)') subname//' ntiles = ',ntiles
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  

    !---------------------------------
    ! get start and end indices of each tile and their PET
    !---------------------------------

    allocate(xb(ntiles),xe(ntiles),yb(ntiles),ye(ntiles),pe(ntiles))
    call mpp_get_compute_domains(ice_data%domain, xbegin=xb, xend=xe, ybegin=yb, yend=ye)
    call mpp_get_pelist(ice_data%domain, pe)
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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGetCoord(gridIn, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_xcen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_ycen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGetCoord(gridIn, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=dataPtr_xcor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=dataPtr_ycor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_mask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_AREA, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_area, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !---------------------------------
    ! load up area, mask, center and corner values
    ! area, mask, and centers should be same size in mom and esmf grid
    ! corner points may not be, need to offset corner points by 1 in i and j
    !   for esmf and also need to "make up" j=1 values.  use wraparound in i
    !---------------------------------

    call mpp_get_compute_domain(ice_data%domain, isc, iec, jsc, jec)

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
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
         msg=SUBNAME//": fld and grid do not have the same size.", &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    allocate(ofld(isc:iec,jsc:jec))
    allocate(gfld(nxg,nyg))

!   call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)
    ofld = ice_data%sCS%G%mask2dT
    write(tmpstr,*) subname//' ofld mask = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ice_data%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld mask = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo

!    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'area', ofld, isc, jsc)
    ofld = ice_data%area
    write(tmpstr,*) subname//' ofld area = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ice_data%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld area = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_area(i,j) = ofld(i1,j1)
    enddo
    enddo

!    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'tlon', ofld, isc, jsc)
    ofld = ice_data%sCS%G%geoLonT
    write(tmpstr,*) subname//' ofld xt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ice_data%domain, ofld, gfld)
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

!    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'tlat', ofld, isc, jsc)
    ofld = ice_data%sCS%G%geoLatT
    write(tmpstr,*) subname//' ofld yt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ice_data%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld yt = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_ycen(i,j) = ofld(i1,j1)
    enddo
    enddo

!    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'ulon', ofld, isc, jsc)
    ofld = ice_data%sCS%G%geoLonCu
    write(tmpstr,*) subname//' ofld xu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ice_data%domain, ofld, gfld)
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

!    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'ulat', ofld, isc, jsc)
    ofld = ice_data%sCS%G%geoLatCu
    write(tmpstr,*) subname//' ofld yu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ice_data%domain, ofld, gfld)
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

    call SIS_RealizeFields(importState, gridIn , fldsToIce_num, fldsToIce, "Ice import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call SIS_RealizeFields(exportState, gridOut, fldsFrIce_num, fldsFrIce, "Ice export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(exportState, itemSearch="mean_net_lw_flx", itemCount=icount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! Do sst initialization if it's part of export state
!    if(icount /= 0) then
!      call ESMF_StateGet(exportState, itemName='sea_surface_temperature', field=field_t_surf, rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out
!      call ESMF_FieldGet(field_t_surf, localDe=0, farrayPtr=t_surf, rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out
!      call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)
!      lbnd1 = lbound(t_surf,1)
!      ubnd1 = ubound(t_surf,1)
!      lbnd2 = lbound(t_surf,2)
!      ubnd2 = ubound(t_surf,2)
!      do j = lbnd2, ubnd2
!      do i = lbnd1, ubnd1
!         j1 = j - lbnd2 + jsc
!         i1 = i - lbnd1 + isc
!         if (ofld(i1,j1) == 0.) t_surf(i,j) = 0.0
!      enddo
!      enddo
!    endif
    deallocate(ofld)

    call NUOPC_Write(exportState, fileNamePrefix='init_field_ice_export_', &
      timeslice=1, relaxedFlag=.true., rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    write(*,*) '----- SIS2 initialization phase Realize completed'

  end subroutine InitializeRealize

  subroutine SIS_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(sis2_cap:SIS_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine SIS_AdvertiseFields
  !-----------------------------------------------------------------------------

  subroutine SIS_RealizeFields(state, grid, nfields, field_defs, tag, rc)

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
    character(len=*),parameter  :: subname='(sis2_cap:SIS_RealizeFields)'

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

  end subroutine SIS_RealizeFields

  subroutine SIS_FieldsSetup(atmos_ice_boundary, ocean_ice_boundary, ice_data)
    type(ice_data_type),           intent(in) :: ice_data
    type(atmos_ice_boundary_type), intent(in) :: atmos_ice_boundary
    type(ocean_ice_boundary_type), intent(in) :: ocean_ice_boundary

    character(len=*),parameter  :: subname='(sis2_cap:SIS_FieldsSetup)'

!--------- import fields to Sea Ice -------------
    !--from ATM boundary
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_net_lw_flx", data=atmos_ice_boundary%lw_flux(:,:,1))
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_net_sw_vis_dir_flx",data=atmos_ice_boundary%sw_flux_vis_dir(:,:,1))
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_net_sw_vis_dif_flx",data=atmos_ice_boundary%sw_flux_vis_dif(:,:,1))
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_net_sw_ir_dir_flx",data=atmos_ice_boundary%sw_flux_nir_dir(:,:,1))
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_net_sw_ir_dif_flx",data=atmos_ice_boundary%sw_flux_nir_dif(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%u_flux(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%v_flux(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%u_star(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%t_flux(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%q_flux(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%lprec(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%fprec(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%dhdt(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%dedt(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%drdt(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%coszen(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=atmos_ice_boundary%p(:,:,1))
    !--from LND boundary
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=land_ice_boundary%runoff(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=land_ice_boundary%calving(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=land_ice_boundary%runoff_hflx(:,:,1))
!    call fld_list_add(fldsToIce_num, fldsToIce, "", data=land_ice_boundary%calving_hflx(:,:,1))
    !--From OCN boundary. These should correspond to fields in OCN export fields.
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_zonal",          data=ocean_ice_boundary%u(:,:))
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_merid",          data=ocean_ice_boundary%v(:,:))
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_temperature",    data=ocean_ice_boundary%t(:,:))
    call fld_list_add(fldsToIce_num, fldsToIce, "s_surf",                     data=ocean_ice_boundary%s(:,:))
    call fld_list_add(fldsToIce_num, fldsToIce, "freezing_melting_potential", data=ocean_ice_boundary%frazil(:,:))

!--------- export fields from Sea Ice -------------
    call fld_list_add(fldsFrIce_num, fldsFrIce, "sea_ice_temperature"    ,data=ice_data%t_surf(:,:,1)         )
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dir_albedo",data=ice_data%albedo_vis_dir(:,:,1))
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dir_albedo", data=ice_data%albedo_nir_dir(:,:,1))
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dif_albedo",data=ice_data%albedo_vis_dif(:,:,1))
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dif_albedo", data=ice_data%albedo_nir_dif(:,:,1))

  end subroutine SIS_FieldsSetup

  subroutine fld_list_add(num, fldlist, stdname, data, shortname)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    real(ESMF_KIND_R8), dimension(:,:),optional, target :: data
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(sis2_cap:fld_list_add)'
    character(len=*), parameter :: transferOffer='can provide'

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

end module sis2_cap_mod
