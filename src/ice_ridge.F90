!> A partially implmented ice ridging parameterizations
!! that does not yet work with an arbitrary number of vertical layers in the ice
module ice_ridging_mod

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! replaced T. Martin code with wrapper for Icepack ridging function - mw 1/18  !
!                                                                              !
! Prioritized to do list as of 6/4/19 (mw):                                    !
!                                                                              !
! 1) implement new snow_to_ocean diagnostic to record this flux.               !
! 2) implement ridging_rate diagnostics: ridging_shear, ridging_conv           !
! 3) implement "do_j" style optimization as in "compress_ice" or               !
!    "adjust_ice_categories" (SIS_transport.F90) if deemed necessary           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_domains,       only : pass_var, pass_vector, BGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, read_param, log_version, param_file_type
use MOM_unit_scaling,  only : unit_scale_type
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_types,         only : ice_state_type
use fms_io_mod,        only : register_restart_field, restart_file_type
use SIS_tracer_registry, only : SIS_tracer_registry_type, SIS_tracer_type
use SIS2_ice_thm,    only : get_SIS2_thermo_coefs
use ice_grid,          only : ice_grid_type
!Icepack modules
use icepack_kinds
use icepack_itd, only: icepack_init_itd, cleanup_itd
use icepack_mechred, only: ridge_ice
use icepack_warnings, only: icepack_warnings_flush, icepack_warnings_aborted, &
                            icepack_warnings_setabort
use icepack_tracers, only: icepack_init_tracer_indices

implicit none ; private

#include <SIS2_memory.h>

public :: ice_ridging, ridge_rate

real, parameter :: hlim_unlim = 1.e8   !< Arbitrary huge number used in ice_ridging
real    :: s2o_frac       = 0.5        !< Fraction of snow dumped into ocean during ridging [nondim]
logical :: rdg_lipscomb = .true.       !< If true, use the Lipscomb ridging scheme
                                       !! TODO: These parameters belong in a control structure

! future namelist parameters?
integer (kind=int_kind), parameter :: &
     krdg_partic = 0  , & ! 1 = new participation, 0 = Thorndike et al 75
     krdg_redist = 0      ! 1 = new redistribution, 0 = Hibler 80

! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
real(kind=dbl_kind), parameter ::  mu_rdg = 3.0

contains



!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ridge_rate finds the deformation rate or total energy dissipation rate due
!!   to ridging (Flato and Hibler, 1995, JGR) or the net area loss in riding
!!   (CICE documentation) depending on the state of the ice drift
function ridge_rate(del2, div) result (rnet)
  real, intent(in)  :: del2 !< The magnitude squared of the shear rates [T-2 ~> s-2].
  real, intent(in)  :: div  !< The ice flow divergence [T-1 ~> s-1]

  ! Local variables
  real :: rnet ! The net rate of area loss or energy dissipation rate due to ridging [T-1 ~> s-1]
  real :: del, rconv, rshear ! Local variables with rates [T-1 ~> s-1]
  !TOM> cs is now set in namelist:
  !Niki: this was commented out
  real, parameter   :: cs=0.25 !(CICE documentation)

  del=sqrt(del2)

  rconv  = -min(div,0.0)           ! energy dissipated by convergence ...
  rshear = 0.5*(del-abs(div))      ! ... and by shear
  rnet   = rconv + cs*rshear       ! net energy contains only part of the
                                   !  shear energy as only a fraction is
           !  dissipated in the ridging process
end function ridge_rate

!
! ice_ridging is a wrapper for the icepack ridging routine ridge_ice
!
subroutine ice_ridging(IST, G, IG, mca_ice, mca_snow, mca_pond, TrReg, US, dt)
  type(ice_state_type),              intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type),                       intent(inout) :: G !<  G The ocean's grid structure.
  type(ice_grid_type),                           intent(inout) :: IG !<   The sea-ice-specific grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),   intent(inout) :: mca_ice, mca_snow, mca_pond
  type(SIS_tracer_registry_type),                pointer       :: TrReg  ! TrReg - The registry of registered SIS ice and snow tracers.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors.
  real (kind=dbl_kind),                          intent(in)    :: dt !<   The amount of time over which the ice dynamics are to be.
                                                                     !    advanced in seconds.
!  logical,                                       intent(in)    :: dyn_Cgrid !<  True if using C-gridd velocities, B-grid if False.



  ! these strain metrics are calculated here from the velocities used for advection
  real :: sh_Dt ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                ! all metric terms, in s-1.
  real :: sh_Dd ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                ! metric terms, in s-1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                ! including all metric terms, in s-1.

  integer :: i, j, k ! loop vars
  integer isc, iec, jsc, jec ! loop bounds
  integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.

  integer (kind=int_kind) :: &
       ncat  , & ! number of thickness categories
       nilyr , & ! number of ice layers
       nslyr  ! number of snow layers


  real (kind=dbl_kind), dimension(0:IG%CatIce) :: hin_max   ! category limits (m)

  logical (kind=log_kind) :: &
       closing_flag, &! flag if closing is valid
       tr_brine       ! if .true., brine height differs from ice thickness

  ! optional history fields
  real (kind=dbl_kind) :: &
       dardg1dt   , & ! rate of fractional area loss by ridging ice (1/s)
       dardg2dt   , & ! rate of fractional area gain by new ridges (1/s)
       dvirdgdt   , & ! rate of ice volume ridged (m/s)
       opening    , & ! rate of opening due to divergence/shear (1/s)
       closing    , & ! rate of closing due to divergence/shear (1/s)
       fpond      , & ! fresh water flux to ponds (kg/m^2/s)
       fresh      , & ! fresh water flux to ocean (kg/m^2/s)
       fhocn          ! net heat flux to ocean (W/m^2)

  real (kind=dbl_kind), dimension(IG%CatIce) :: &
       dardg1ndt  , & ! rate of fractional area loss by ridging ice (1/s)
       dardg2ndt  , & ! rate of fractional area gain by new ridges (1/s)
       dvirdgndt  , & ! rate of ice volume ridged (m/s)
       aparticn   , & ! participation function
       krdgn      , & ! mean ridge thickness/thickness of ridging ice
       araftn     , & ! rafting ice area
       vraftn     , & ! rafting ice volume
       aredistn   , & ! redistribution function: fraction of new ridge area
       vredistn       ! redistribution function: fraction of new ridge volume

  real (kind=dbl_kind), dimension(0) :: &
       faero_ocn      ! aerosol flux to ocean (kg/m^2/s)

  real (kind=dbl_kind), dimension(0) :: &
       fiso_ocn       ! isotope flux to ocean (kg/m^2/s)

  integer (kind=int_kind) :: &
       ndtd = 1  , & ! number of dynamics subcycles
       n_aero = 0, & ! number of aerosol tracers
       ntrcr = 0     ! number of tracer level


  real(kind=dbl_kind) :: &
       del_sh        , & ! shear strain measure
       rdg_conv = 0.0, & ! normalized energy dissipation from convergence (1/s)
       rdg_shear= 0.0    ! normalized energy dissipation from shear (1/s)

  real(kind=dbl_kind), dimension(IG%CatIce) :: &
       aicen, & ! concentration of ice
       vicen, & ! volume per unit area of ice          (m)
       vsnon    ! volume per unit area of snow         (m)

  ! ice tracers; ntr*(NkIce+NkSnow) guaranteed to be enough for all (intensive)
  real(kind=dbl_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),IG%CatIce) :: trcrn

  real(kind=dbl_kind) :: aice0          ! concentration of open water

  integer (kind=int_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow)) :: &
       trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon (weighting to use)
       n_trcr_strata  ! number of underlying tracer layers

  real(kind=dbl_kind), dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),3) :: &
       trcr_base      ! = 0 or 1 depending on tracer dependency
                    ! argument 2:  (1) aice, (2) vice, (3) vsno

  integer, dimension(TrReg%ntr*(IG%NkIce+IG%NkSnow),IG%CatIce) :: &
       nt_strata      ! indices of underlying tracer layers

  type(SIS_tracer_type), dimension(:), pointer :: Tr=>NULL() ! SIS2 tracers


  real :: rho_ice, rho_snow
  integer :: m, n ! loop vars for tracer; n is tracer #; m is tracer layer

  nSlyr = IG%NkSnow
  nIlyr = IG%NkIce
  nCat  = IG%CatIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  call get_SIS2_thermo_coefs(IST%ITV, rho_snow=rho_snow)

  ! copy strain calculation code from SIS_C_dynamics; might be a more elegant way ...
  !
  halo_sh_Ds = min(isc-G%isd, jsc-G%jsd, 2)
!  if (dyn_Cgrid) then
     do J=jsc-halo_sh_Ds,jec+halo_sh_Ds-1 ; do I=isc-halo_sh_Ds,iec+halo_sh_Ds-1
       ! This uses a no-slip boundary condition.
       sh_Ds(I,J) = (2.0-G%mask2dBu(I,J)) * &
            (G%dxBu(I,J)*G%IdyBu(I,J)*(IST%u_ice_C(I,j+1)*G%IdxCu(I,j+1) - IST%u_ice_C(I,j)*G%IdxCu(I,j)) + &
            G%dyBu(I,J)*G%IdxBu(I,J)*(IST%v_ice_C(i+1,J)*G%IdyCv(i+1,J) - IST%v_ice_C(i,J)*G%IdyCv(i,J)))
     enddo; enddo
 ! else
 !    call SIS_error(FATAL,"Icepack ridging only implemented for C-grid versions of SIS2")
 ! endif

  ! set category limits; Icepack has a max on the largest, unlimited, category (why?)
  do k=1,nCat
     hin_max(k) = IG%mH_cat_bound(k)/(Rho_ice*IG%kg_m2_to_H)
  end do
  hin_max(nCat+1) = 1e5; ! not sure why this is needed, set big

  trcr_base = 0.0; n_trcr_strata = 0; nt_strata = 0; ! init some tracer vars
  ! When would we use icepack tracer "strata"?

  ! set icepack tracer index "nt_lvl" to (last) pond tracer so it gets dumped when
  ! ridging in ridge_ice (this is what happens to "level" ponds); first add up ntrcr;
  ! then set nt_lvl to ntrcr+1; could move this to an initializer - mw
  ntrcr = 0
  if (TrReg%ntr>0) then ! sum tracers
    Tr => TrReg%Tr_snow
    do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
      ntrcr = ntrcr + 1
    enddo ; enddo
    Tr => TrReg%Tr_ice
    do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
          ntrcr = ntrcr + 1
    enddo ; enddo
  endif
  call icepack_init_tracer_indices(nt_vlvl_in=ntrcr+1); ! pond will be last tracer

  do j=jsc,jec; do i=isc,iec;
  if ((G%mask2dT(i,j) .gt. 0.0) .and. (sum(IST%part_size(i,j,1:nCat)) .gt. 0.0)) then
  ! feed locations to Icepack's ridge_ice

    ! start like we're putting ALL the snow in the ocean
    IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) + sum(mca_snow(i,j,:))
    IST%enth_snow_to_ocn(i,j) = IST%enth_snow_to_ocn(i,j) + sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1));
    IST%water_to_ocn(i,j) = IST%water_to_ocn(i,j) + sum(mca_pond(i,j,:));
    aicen = IST%part_size(i,j,1:nCat);

    if (sum(aicen) .le. 0.0) then ! no ice -> no ridging
      IST%part_size(i,j,0) = 1.0
    else

      ! set up ice and snow volumes
      vicen = mca_ice(i,j,:) /Rho_ice
      vsnon = mca_snow(i,j,:)/Rho_snow

      sh_Dt = (G%dyT(i,j)*G%IdxT(i,j)*(G%IdyCu(I,j) * IST%u_ice_C(I,j) - &
                                       G%IdyCu(I-1,j)*IST%u_ice_C(I-1,j)) - &
               G%dxT(i,j)*G%IdyT(i,j)*(G%IdxCv(i,J) * IST%v_ice_C(i,J) - &
                                       G%IdxCv(i,J-1)*IST%v_ice_C(i,J-1)))
      sh_Dd = (G%IareaT(i,j)*(G%dyCu(I,j) * IST%u_ice_C(I,j) - &
                              G%dyCu(I-1,j)*IST%u_ice_C(I-1,j)) + &
               G%IareaT(i,j)*(G%dxCv(i,J) * IST%v_ice_C(i,J) - &
                              G%dxCv(i,J-1)*IST%v_ice_C(i,J-1)))

      del_sh = sqrt(sh_Dd**2 + 0.25 * (sh_Dt**2 + &
                   (0.25 * ((sh_Ds(I-1,J-1) + sh_Ds(I,J)) + &
                            (sh_Ds(I-1,J) + sh_Ds(I,J-1))))**2 ) ) ! H&D eqn 9
      rdg_conv  = -min(sh_Dd,0.0)              ! energy dissipated by convergence ...
      rdg_shear = 0.5*(del_sh-abs(sh_Dd))      ! ... and by shear

!rdg_shear = 86400.0/100.0 ! ... for column model testing

      aice0 = 1.0+dt*sh_Dd-sum(aicen)

      ntrcr = 0
      if (TrReg%ntr>0) then ! load tracer array

        Tr => TrReg%Tr_snow
        do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
          ntrcr = ntrcr + 1
          trcrn(ntrcr,:) = Tr(n)%t(i,j,:,m)
          trcr_depend(ntrcr) = 2; ! 2 means snow-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,3) = 1.0; ! 3rd index for snow
        enddo ; enddo

        Tr => TrReg%Tr_ice
        do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
          ntrcr = ntrcr + 1
          trcrn(ntrcr,:) = Tr(n)%t(i,j,:,m)
          trcr_depend(ntrcr) = 1; ! 1 means ice-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,2) = 1.0; ! 2nd index for ice
        enddo ; enddo

      endif ! have tracers to load

      ! load pond on top of stack
      ntrcr = ntrcr + 1
      trcrn(ntrcr,:) = IST%mH_pond(i,j,:)
      trcr_depend(ntrcr) = 0; ! 1 means ice area-based tracer
      trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,1) = 1.0; ! 1st index for ice area

      tr_brine = .false.
      dardg1dt=0.0
      dardg2dt=0.0
      opening=0.0
      fpond=0.0
      fresh=0.0
      fhocn=0.0
      faero_ocn=0.0
      fiso_ocn=0.0
      aparticn=0.0
      krdgn=0.0
      aredistn=0.0
      vredistn=0.0
      dardg1ndt=0.0
      dardg2ndt=0.0
      dvirdgndt=0.0
      araftn=0.0
      vraftn=0.0
      closing_flag=.false.

      ! call Icepack routine; how are ponds treated?
      call ridge_ice (dt,           ndtd,           &
                      ncat,         n_aero,         &
                      nilyr,        nslyr,          &
                      ntrcr,        hin_max,        &
                      rdg_conv,     rdg_shear,      &
                      aicen,                        &
                      trcrn,                        &
                      vicen,        vsnon,          &
                      aice0,                        &
                      trcr_depend,                  &
                      trcr_base,                    &
                      n_trcr_strata,                &
                      nt_strata,                    &
                      krdg_partic,  krdg_redist,    &
                      mu_rdg,       tr_brine,       &
                      dardg1dt,     dardg2dt,       &
                      dvirdgdt,     opening,        &
                      fpond,                        &
                      fresh,        fhocn,          &
                      krdgn=krdgn,             &
                      aredistn=aredistn,       &
                      vredistn=vredistn,       &
                      dardg1ndt=dardg1ndt,     &
                      dardg2ndt=dardg2ndt,      &
                      dvirdgndt=dvirdgndt,                    &
                      araftn=araftn,       &
                      vraftn=vraftn,         &
                      closing_flag=closing_flag )


      if ( icepack_warnings_aborted() ) then
        call icepack_warnings_flush(0);
        call icepack_warnings_setabort(.false.)
        call SIS_error(WARNING,'icepack ridge_ice error');
      endif

      ! pop pond off top of stack
      IST%mH_pond(i,j,:) = trcrn(ntrcr,:)
      mca_pond(i,j,:) = IST%mH_pond(i,j,:)*aicen

      if (TrReg%ntr>0) then
        ! unload tracer array reversing order of load -- stack-like fashion

        Tr => TrReg%Tr_ice
        do n=TrReg%ntr,1,-1 ; do m=Tr(n)%nL,1,-1
          ntrcr = ntrcr - 1
          Tr(n)%t(i,j,:,m) = trcrn(ntrcr,:)
        enddo ; enddo

        Tr => TrReg%Tr_snow
        do n=TrReg%ntr,1,-1 ; do m=Tr(n)%nL,1,-1
          ntrcr = ntrcr - 1
          Tr(n)%t(i,j,:,m) = trcrn(ntrcr,:)
        enddo ; enddo

      endif ! have tracers to unload

      ! output: snow/ice masses/thicknesses
      do k=1,nCat
        if (aicen(k) > 0.0) then
          IST%part_size(i,j,k)  = aicen(k)
          mca_ice(i,j,k)  = vicen(k)*Rho_ice
          IST%mH_ice(i,j,k)   = vicen(k)*Rho_ice/aicen(k)
          mca_snow(i,j,k) = vsnon(k)*Rho_snow
          IST%mH_snow(i,j,k)  = vsnon(k)*Rho_snow/aicen(k)
        else
          IST%part_size(i,j,k) = 0.0
          mca_ice(i,j,k)  = 0.0
          IST%mH_ice(i,j,k) = 0.0
          mca_snow(i,j,k) = 0.0
          IST%mH_snow(i,j,k) = 0.0
        endif
        ! How to treat ponds?
      enddo

      ! negative part_sz(i,j,0) triggers compress_ice clean_up later
      IST%part_size(i,j,0) = 1.0 - sum(IST%part_size(i,j,1:nCat))

    endif
    ! subtract new snow/pond mass and energy on ice to sum net fluxes to ocean
    IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) - sum(mca_snow(i,j,:));
    IST%enth_snow_to_ocn(i,j) = IST%enth_snow_to_ocn(i,j) - sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1));
    IST%water_to_ocn(i,j) = IST%water_to_ocn(i,j) - sum(mca_pond(i,j,:));

  endif; enddo; enddo ! part_sz, j, i

end subroutine ice_ridging


! !TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! !> ice_ridging parameterizes mechanical redistribution of thin (undeformed) ice
! !! into thicker (deformed/ridged) ice categories
! subroutine ice_ridging(km, cn, hi, hs, rho_ice, t1, t2, age, snow_to_ocn, enth_snow_to_ocn, rdg_rate, hi_rdg, &
!                        dt, hlim_in, rdg_open, vlev, US)
!   !  Subroutine written by T. Martin, 2008
!   integer,             intent(in)    :: km  !< The number of ice thickness categories
!   real, dimension(0:km), intent(inout) :: cn  !< Fractional concentration of each thickness category,
!                                             !! including open water fraction [nondim]
!   real, dimension(km), intent(inout) :: hi  !< ice mass in each category [R Z ~> kg m-2]
!   real, dimension(km), intent(inout) :: hs  !< snow mass in each category [R Z ~> kg m-2]
!   real,                intent(in)    :: rho_ice !< Nominal ice density [R ~> kg m-3]
!    ! CAUTION: these quantities are extensive here,
!   real, dimension(km), intent(inout) :: t1  !< Volume integrated upper layer temperature [degC m3]?
!   real, dimension(km), intent(inout) :: t2  !< Volume integrated upper layer temperature [degC m3]?
!   real, dimension(km), intent(inout) :: age !< Volume integrated ice age [m3 years]?
!   real,                intent(out)   :: enth_snow_to_ocn !< average of enthalpy of the snow dumped into
!                                             !! ocean due to this ridging event [Q ~> J kg-1]
!   real,                intent(out)   :: snow_to_ocn !< total snow mass dumped into ocean due to this
!                                             !! ridging event [R Z ~> kg m-2]
!   real,                intent(in)    :: rdg_rate  !< Ridging rate from subroutine ridge_rate [T-1 ~> s-1]
!   real, dimension(km), intent(inout) :: hi_rdg    !< A diagnostic of the ridged ice volume in each
!                                                   !! category [R Z ~> kg m-2].
!   real,                intent(in)    :: dt        !< time step [T ~> s]
!   real, dimension(km), intent(in)    :: hlim_in   !< ice thickness category limits
!   real,                intent(out)   :: rdg_open  !< Rate of change in open water area due to
!                                                   !! newly formed ridges [T-1 ~> s-1]
!   real,                intent(out)   :: vlev      !< mass of level ice participating in ridging [R Z T-1 ~> kg m-2 s-1]
!   type(unit_scale_type), intent(in)  :: US  !< A structure with unit conversion factors

!   ! Local variables
!   integer :: k, kd, kr, n_iterate
!   integer, parameter :: n_itermax = 10 ! maximum number of iterations for redistribution
! !    real, parameter :: frac_hs_rdg = 0.5 ! fraction of snow that remains on ridged ice;
!   real :: frac_hs_rdg  ! fraction of snow that remains on ridged ice [nondim];
!                        !  (1.0-frac_hs_rdg)*hs falls into ocean
!   real, dimension(0:km) :: part_undef  ! fraction of undeformed ice or open water participating in ridging
!   real                  :: area_undef  ! fractional area of parent and ...
!   real, dimension(km) :: area_def    ! ... newly ridged ice, respectively [nondim]
!   real, dimension(km) :: vol_def     ! fractional volume of newly ridged ice [nondim]
!   real, dimension(km) :: cn_old      ! concentrations at beginning of iteration loop [nondim]
!   real, dimension(km) :: hi_old      ! thicknesses at beginning of iteration loop [R Z ~> kg m-2]
!   real, dimension(km) :: rdg_frac    ! ratio of ridged and total ice volume [nondim]
!   real                :: alev        ! area of level ice participating in ridging [nondim]
!   real                :: ardg, vrdg  ! area and volume of newly formed rdiged (vlev=vrdg!!!)
!   real, dimension(km) :: hmin, hmax, efold ! [R Z ~> kg m-2]
!   real, dimension(km) :: rdg_ratio ! [nondim]
!   real, dimension(km) :: hlim      ! [R Z ~> kg m-2]
!   real :: hl, hr
!   real :: snow_dump, enth_dump
!   real :: cn_tot, part_undef_sum
!   real :: div_adv, Rnet, Rdiv, Rtot ! [T-1 ~> s-1]
!   real :: rdg_area, rdgtmp, hlimtmp
!   real                  :: area_frac
!   real, dimension(km) :: area_rdg ! [nondim]
!   real, dimension(km) :: &
!     frac_hi, frac_hs, &    ! Portion of the ice and snow that are being ridged [R Z ~> kg m-2]
!     frac_t1, frac_t2, &
!     frac_age
!   logical               :: rdg_iterate
!   !-------------------------------------------------------------------
!   ! some preparations
!   !-------------------------------------------------------------------
!   hlimtmp = hlim_in(km)
!   hlim(km) = hlim_unlim   ! ensures all ridged ice is smaller than thickest ice allowed
!   frac_hs_rdg = 1.0 - s2o_frac
!   snow_to_ocn = 0.0 ; enth_snow_to_ocn = 0.0
!   alev=0.0 ; ardg=0.0 ; vlev=0.0 ; vrdg=0.0
!   !
!   call ice_ridging_init(km, cn, hi, rho_ice, part_undef, part_undef_sum, &
!                         hmin, hmax, efold, rdg_ratio, US)

!   !-------------------------------------------------------------------
!   ! opening and closing rates of the ice cover
!   !-------------------------------------------------------------------

!   ! update total area fraction as this may exceed 1 after transportation/advection
!   ! (cn_tot <= 1 after ridging!)
!   cn_tot = sum(cn(0:km))

!   ! dissipated energy in ridging from state of ice drift
!   !  after Flato and Hibler (1995, JGR)
!   !  (see subroutine ridge_rate in ice_dyn_mod),
!   !  equals net closing rate times ice strength
!   ! by passing to new, local variable rdg_rate is saved for diagnostic output
!   Rnet = rdg_rate
!   ! the divergence rate given by the advection scheme ...
!   div_adv = (1.-cn_tot) / dt
!   ! ... may exceed the rate derived from the drift state (Rnet)
!   if (div_adv < 0.) Rnet = max(Rnet, -div_adv)
!   ! non-negative opening rate that ensures cn_tot <=1 after ridging
!   Rdiv = Rnet + div_adv
!   ! total closing rate
!   Rtot = Rnet / part_undef_sum

!   !-------------------------------------------------------------------
!   ! iteration of ridging redistribution
!   do n_iterate=1, n_itermax
!   !-------------------------------------------------------------------

!     ! save initial state of ice concentration, total and ridged ice volume
!     !  at beginning of each iteration loop
!     do k=1,km
!       cn_old(k) = cn(k)
!       hi_old(k) = hi(k)

!       rdg_frac(k) = 0.0 ; if (hi(k)>0.0) rdg_frac(k) = hi_rdg(k)/hi(k)
!     enddo

!     ! reduce rates in case more than 100% of any category would be removed
!     do k=1,km ; if (cn(k)>1.0e-10 .and. part_undef(k)>0.0) then
!       rdg_area = part_undef(k) * Rtot * dt   ! area ridged in category k
!       if (rdg_area > cn(k)) then
!         rdgtmp = cn(k)/rdg_area
!         Rtot = Rtot * rdgtmp
!         Rdiv = Rdiv * rdgtmp
!       endif
!     endif ; enddo

!     !-------------------------------------------------------------------
!     ! redistribution of ice
!     !-------------------------------------------------------------------

!     ! changes in open water area
!     cn(0) = max(cn(0) + (Rdiv - part_undef(1)*Rtot) * dt, 0.0)

!     if (Rtot>0.0) then

!       ! area, volume and energy changes in each category
!       do kd=1,km   ! donating category
!         area_undef = min(part_undef(kd)*Rtot*dt, cn_old(kd))   ! area that experiences ridging in category k,
!                                                                ! make sure that not more than 100% are used
!         if (cn_old(kd) > 1.0e-10) then
!           area_frac    = area_undef / cn_old(kd)              ! fraction of level ice area involved in ridging
!           area_rdg(kd) = area_undef / rdg_ratio(kd)           ! area of new ridges in category k
!         else
!           area_frac    = 0.0
!           area_rdg(kd) = 0.0
!         endif
!         !if (rdg_ratio(kd) > 0.0) then     ! distinguish between level and ridged ice in
!         !else                              !  each category: let only level ice ridge;
!         !endif                             !  sea also change of hi_rdg below

!         ! reduce area, volume and energy of snow and ice in source category
!         frac_hi(kd)  = hi(kd)   * area_frac
!         frac_hs(kd)  = hs(kd)   * area_frac
!         frac_t1(kd)  = t1(kd)   * area_frac
!         frac_t2(kd)  = t2(kd)   * area_frac
!         frac_age(kd) = age(kd)  * area_frac

!         cn(kd)  = cn(kd)  - area_undef
!         hi(kd)  = hi(kd)  - frac_hi(kd)
!         hs(kd)  = hs(kd)  - frac_hs(kd)
!         t1(kd)  = t1(kd)  - frac_t1(kd)
!         t2(kd)  = t2(kd)  - frac_t2(kd)
!         age(kd) = age(kd) - frac_age(kd)

!         alev = alev + area_undef   ! diagnosing area of level ice participating in ridging
!         vlev = vlev + frac_hi(kd)  ! diagnosing total ice volume moved due to ridging
!                                    !  (here donating categories)

!         !    Here it is assumed that level and ridged ice
!         !  of a category participate in ridging in equal
!         !  measure; this also means that ridged ice may be ridged again
!         hi_rdg(kd) = hi_rdg(kd) - rdg_frac(kd)*frac_hi(kd)
!         hi_rdg(kd) = max(hi_rdg(kd), 0.0)      ! ensure hi_rdg >= 0

!         ! dump part of the snow in ocean (here, sum volume, transformed to flux in update_ice_model_slow)
!         snow_dump = frac_hs(kd)*(1.0-frac_hs_rdg)
!         if (snow_to_ocn > 0.0) then
!           enth_dump = t1(kd)  !### THIS IS WRONG, BUT IS A PLACEHOLDER FOR NOW.
!           enth_snow_to_ocn = (enth_snow_to_ocn*snow_to_ocn + enth_dump*snow_dump) / (snow_to_ocn + snow_dump)
!           snow_to_ocn = snow_to_ocn + snow_dump
!         endif

!       enddo

!       ! split loop in order to derive frac_... variables with initial status (before ridging redistribution)
!       do kd=1,km

!         !----------------------------------------------------------------------------------------
!         ! add area, volume and energy in receiving category :
!         ! A) after Lipscomb, 2007 (negative exponential distribution)
!         ! B) after Hibler, 1980, Mon. Weather Rev. (uniform distribution)
!         !----------------------------------------------------------------------------------------
!         if (rdg_lipscomb) then
!           ! ************
!           ! *   A      *
!           ! ************
!           if (efold(kd)>0.0) then
!             do kr=1,km-1   ! receiving categories
!               if (hmin(kd) >= hlim(kr)) then
!                 area_def(kr) = 0.0
!                 vol_def(kr)  = 0.0
!               else
!                 hl = max(hmin(kd), hlim(kr-1))
!                 hr = hlim(kr)
!                 area_def(kr) = exp((hmin(kd)-hl)/efold(kd))   &
!                  -             exp((hmin(kd)-hr)/efold(kd))
!                 vol_def(kr)  = ( (hl+efold(kd))*exp((hmin(kd)-hl)/efold(kd))   &
!                  -   (hr+efold(kd))*exp((hmin(kd)-hr)/efold(kd)) ) &
!                  / (hmin(kd)+efold(kd))
!               endif
!             enddo   ! k receiving
!             ! thickest categery is a special case:
!             hl = max(hmin(kd), hlim(km-1))
!             area_def(km) =                  exp((hmin(kd)-hl)/efold(kd))
!             vol_def(km)  = ( (hl+efold(kd))*exp((hmin(kd)-hl)/efold(kd)) ) &
!              / (hmin(kd)+efold(kd))
!           else
!             do kr=1,km
!               area_def(kr) = 0.0
!               vol_def(kr)  = 0.0
!             enddo
!           endif
!         !----------------------------------------------------------------------------------------
!         else ! not rdg_lipscomb
!           ! ************
!           ! *   B      *
!           ! ************
!           if (hmax(kd)==hmin(kd)) then
!             do kr=1,km ; area_def(kr) = 0.0 ; vol_def(kr)  = 0.0 ; enddo
!           else
!             do kr=1,km   ! receiving categories
!               if (hmin(kd) >= hlim(kr) .or. hmax(kd) <= hlim(kr-1)) then
!                 hl = 0.0
!                 hr = 0.0
!               else
!                 hl = max(hmin(kd), hlim(kr-1))
!                 hr = min(hmax(kd), hlim(kr)  )
!               endif
!               area_def(kr) = (hr   -hl   ) / (hmax(kd)   -hmin(kd)   )
!               !vol_def(kr) = (hr**2-hl**2) / (hmax(kd)**2-hmin(kd)**2)
!               vol_def(kr)  = area_def(kr) * (hr+hl) / (hmax(kd)+hmin(kd))
!             enddo   ! k receiving
!           endif

!         endif
!         !----------------------------------------------------------------------------------------

!         ! update ice/snow area, volume, energy for receiving categories
!         do kr=1,km   ! receiving categories
!           cn(kr)  = cn(kr)  + area_def(kr) * area_rdg(kd)
!           rdgtmp  = vol_def(kr)  * frac_hi(kd)
!           hi(kr)  = hi(kr)  + rdgtmp
!           hs(kr)  = hs(kr)  + vol_def(kr)  * frac_hs(kd) * frac_hs_rdg
!           t1(kr)  = t1(kr)  + vol_def(kr)  * frac_t1(kd)
!           t2(kr)  = t2(kr)  + vol_def(kr)  * frac_t2(kd)
!           age(kr) = age(kr) + vol_def(kr)  * frac_age(kd)

!           ardg = ardg + area_def(kr) * area_rdg(kd) ! diagnosing area of newly ridged ice
!           vrdg = vrdg + rdgtmp                      ! diagnosing total ice volume moved due to ridging
!                                                     !  (here receiving categories, cross check with vlev)

!           ! add newly ridged ice volume to total ridged ice in each category
!           hi_rdg(kr) = hi_rdg(kr) + rdgtmp
!         enddo

!       enddo   ! kd loop over donating categories

!     endif ! Rtot>0.0

!     ! update total area fraction and check if this is now <= 1
!     ! and rerun the ice redistribution when necessary
!     cn_tot = sum(cn(0:km))
!     rdg_iterate = .false.
!     if (abs(cn_tot-1.) > 1.0e-10) then
!        rdg_iterate = .true.
!        div_adv = (1.-cn_tot) / dt
!        Rnet    = max(0.0, -div_adv)
!        Rdiv    = max(0.0,  div_adv)
!        call ice_ridging_init(km, cn, hi, rho_ice, part_undef, part_undef_sum, &
!                              hmin, hmax, efold, rdg_ratio, US)
!        Rtot    = Rnet / part_undef_sum
!     endif

!     !-------------------------------------------------------------------
!     if (.not. rdg_iterate) exit
!   enddo   ! n_iterate
!   !-------------------------------------------------------------------

!   ! check ridged ice volume for natural limits
!   do k=1,km
!     hi_rdg(k) = max(hi_rdg(k),0.0)   ! ridged ice volume positive
!     hi_rdg(k) = min(hi_rdg(k),hi(k)) ! ridged ice volume not to exceed total ice volume
!   enddo

!   ! calculate opening rate of ridging
!   rdg_open = (alev - ardg) / dt

!   ! cross check ice volume transferred from level to ridged ice
!   if (abs(vlev - vrdg) > 1.0e-10*US%kg_m3_to_R*US%m_to_Z) then
!     print *,'WARNING: subroutine ice_ridging: parent ice volume does not equal ridge volume', vlev, vrdg
!   endif
!   ! turn vlev into a rate for diagnostics
!   vlev = vlev / dt

!   ! return to true upper most ice thickness category limit
!   !hlim(km) = hlimtmp

! end subroutine ice_ridging

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_ridging_end deallocates the memory associated with this module.
subroutine ice_ridging_end()

end subroutine ice_ridging_end

end module ice_ridging_mod
