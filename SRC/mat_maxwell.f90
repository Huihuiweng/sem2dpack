module mat_maxwell

  use prop_mat

  implicit none
  private

 !-- maxwell
  type matwrk_maxwell_type
    private
    double precision, pointer :: eta(:,:) => null()
    double precision, pointer :: sigmaIntegral(:,:,:) => null()
    double precision          :: dt
  end type matwrk_maxwell_type
  
  integer, save :: isMaxwell = 0
  logical, save :: NormalizeByDT=.true.

  ! for memory report
  integer, save :: MAT_Maxwell_mempro = 0
  integer, save :: MAT_Maxwell_memwrk = 0

  public :: matwrk_maxwell_type &
          , MAT_isMaxwell, MAT_Maxwell_read, MAT_Maxwell_init_elem_prop &
          , MAT_Maxwell_init_elem_work, MAT_Maxwell_add_etav &
          , MAT_Maxwell_mempro, MAT_Maxwell_memwrk

contains

!=======================================================================
  logical function MAT_isMaxwell(m)

  type(matpro_elem_type), intent(in) :: m

  MAT_isMaxwell = MAT_isKind(m,isMaxwell)

  end function MAT_isMaxwell

!=======================================================================
! BEGIN INPUT BLOCK
!
! NAME   : MAT_Maxwell
! GROUP  : MATERIALS
! PURPOSE:  XXX
! SYNTAX : &MAT_Maxwell eta, ETAxDT /
!          &MAT_Maxwell etaH, ETAxDT / followed by a DIST_XXX input block
!
! ARG: eta      [dble][0d0] Viscosity coefficient
! ARG: ETAxDT   [log][T] If eta is given in units of dt (timestep)
!
! NOTE   :  XXX
!
! END INPUT BLOCK

 
  subroutine MAT_Maxwell_read(input,iin)

  use echo, only : echo_input, iout
  use stdio, only : IO_abort

  type (matpro_input_type), intent(inout) :: input
  integer, intent(in) :: iin

  double precision :: eta
  character(20) :: etaH
  logical :: ETAxDT

  NAMELIST / MAT_Maxwell / eta, etaH, ETAxDT

  eta = 0d0
  etaH = ' '
  ETAxDT = .true.

  read(iin, MAT_Maxwell, END=100)

  call MAT_setKind(input,isMaxwell)

  call MAT_setProp(input,'eta',eta,etaH,iin,etaH)
  NormalizeByDT = ETAxDT

  if (echo_input) write(iout,200) etaH,ETAxDT

  return
  
  100 call IO_abort('MAT_Maxwell_read: MAT_Maxwell input block not found')

  200   format(5x, &
    'Viscosity coefficient . . . . .(eta/etaH) = ',a,/5x, &
    'Normalized by dt. . . . . . . . .(ETAxDT) = ',L1)

  end subroutine MAT_Maxwell_read

!=======================================================================
  subroutine MAT_Maxwell_init_elem_prop(elem,ecoord)

  type(matpro_elem_type), intent(inout) :: elem
  double precision, intent(in) :: ecoord(:,:,:)

  call MAT_setProp(elem,'eta',ecoord,MAT_Maxwell_mempro)

  end subroutine MAT_Maxwell_init_elem_prop

!=======================================================================
  subroutine MAT_Maxwell_init_elem_work(matwrk,matpro,ngll,dt,ndof)

  type(matwrk_maxwell_type), intent(inout) :: matwrk
  type(matpro_elem_type), intent(in) :: matpro
  double precision, intent(in) :: dt
  integer, intent(in) :: ngll,ndof

  if (.not. MAT_isMaxwell(matpro)) return
  allocate( matwrk%eta(ngll,ngll) )
  call MAT_getProp(matwrk%eta, matpro,'eta')
  if (NormalizeByDT) matwrk%eta = dt * matwrk%eta

  allocate( matwrk%sigmaIntegral(ngll,ngll,ndof) )
  matwrk%sigmaIntegral(:,:,:) = 0.0
  matwrk%dt =dt

  MAT_Maxwell_memwrk = MAT_Maxwell_memwrk &
                + size( transfer(matwrk, (/ 0d0 /) )) &
                + size(matwrk%eta)

  end subroutine MAT_Maxwell_init_elem_work

!=======================================================================

  subroutine MAT_Maxwell_add_etav(d,v,m,ngll,ndof)

  integer, intent(in) :: ngll,ndof
  double precision, intent(inout) :: d(ngll,ngll,ndof)
  double precision, intent(in) :: v(ngll,ngll,ndof)
  type(matwrk_maxwell_type), intent(in) :: m

  integer :: i

  do i=1,ndof
    d(:,:,i) = d(:,:,i) - m%sigmaIntegral(:,:,i) / m%eta

    ! Update the integral of the stress function: int_0^t sigma dt. Here d() is for the stress sigma
    m%sigmaIntegral(:,:,i) = m%sigmaIntegral(:,:,i) + d(:,:,i) * m%dt
  enddo

  end subroutine MAT_Maxwell_add_etav

end module mat_maxwell
