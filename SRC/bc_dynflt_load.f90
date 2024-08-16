module bc_dynflt_load

  implicit none
  private

  type load_type
  private
  integer :: kind_ampl, kind_load
  double precision :: x0, R, T, Delta_t
  end type load_type

  public :: load_type, ld_read, load_form 

contains
!=====================================================================
! BEGIN INPUT BLOCK
!
! NAME   : BC_DYNFLT_LOAD
! GROUP  : DYNAMIC_FAULT
! PURPOSE: Load stress perturbation depends on time and space
! SYNTAX : &BC_DYNFLT_LOAD kind, x0, R, T, Sti /
!
! ARG: kind_ampl     [int] [1] amplitude of the perturbations, kind =
!                        2 = smooth (as in TPV105 description eq.17)
!                        2 = ellipse
!                        ...                   
! ARG: kind_load      [dble] [1] Load perturbations form:
!                        1 = smooth (as in TPV105 description eq.18)
!                        2 = step
!                        ...
! ARG: x0       [dble] [0d0] central point of x
! ARG: R        [dble] [None] Friction coefficient at the hypocenter at time=0
! ARG: T        [dble] [None] Load time
! ARG: Delta_t      [dble] [None] Initial delta_stress
!
! NOTE   : Time-weakening is usually applied as an artificial nucleation procedure.
!          The maximum size of the nucleation region is 2*V*T if kind=1, V*T/2 if kind=2
        
  subroutine ld_read(ld,iin)

  use echo, only : echo_input,iout
  use stdio, only: IO_abort

  type(load_type), intent(out) :: ld
  integer, intent(in) :: iin
  double precision :: x0, R, T, Delta_t
  integer :: kind_load, kind_ampl

  NAMELIST / BC_DYNFLT_LOAD / kind_ampl, kind_load, x0, R, T, Delta_t

  kind_load = 1
  kind_ampl = 1
  x0 = 0d0
  read(iin,BC_DYNFLT_LOAD,END=300)
  300 continue


  ld%kind_ampl = kind_ampl
  ld%kind_load = kind_load
  ld%x0 = x0
  ld%R = R
  ld%T = T
  ld%Delta_t = Delta_t

  if (echo_input) write(iout,200) kind_ampl,kind_load,R,T,Delta_t

  return
  200 format(5x,'Perturbations with time and space:                     ', &
            /5x,'  Type of asperity amplitude  . . . . (kind_ampl) = ',I0, &
            /5x,'  Type of perturbation form   . . . . (kind_load) = ',I0, &
            /5x,'  Perturbation Radius . . . . .  . . . . .(R) = ',EN13.3, &
            /5x,'  Perturbation Time . . . . . . .  . . . .(T) = ',EN13.3, &
            /5x,'  Initial perturbation value. . . . (Delta_t) = ',EN13.3  )

  end subroutine ld_read
        
  function load_form(ld, field, coord, time) result(field_new)
       ! Load perturbation with time
         type(load_type), intent(in) :: ld
         double precision, intent(in) :: coord(:,:),time
         double precision, intent(in), dimension(size(coord(1,:))) :: field

         integer :: j
         double precision :: t
         double precision, dimension(size(field)) :: ampl, field_new
         
         t = time
         write(*,*) 'time =', t
         ampl = load_ampl(ld,coord)
         select case(ld%kind_load)
         case(1)
                  do j=1,size(ampl)
                        if (t < ld%T .and. t > 0d0) then
                                field_new(j) = field(j) + ld%Delta_t* ampl(j) * exp((t - ld%T)**2 / (t * (t - 2*ld%T)))
                        elseif( t >= ld%T) then
                                field_new(j) = field(j) + ld%Delta_t * ampl(j) * 1.0
                        end if 
                  end do
          end select
                         
  end function load_form
         





  function load_ampl(ld,coord) result(ampl)
     ! amplitude of the perturbations form
        type(load_type), intent(in) :: ld
        double precision, intent(in) :: coord(:,:)
        double precision, dimension(size(coord(1,:))) :: ampl, rad_x

        integer :: i

        select case(ld%kind_ampl)
        case(1)
           do i=1,size(rad_x)
                rad_x(i) =  abs(coord(1,i)-ld%x0)
                if (rad_x(i) < ld%R) then
                        ampl(i) = exp(rad_x(i)**2 / (rad_x(i)**2 - ld%R**2))
                else
                        ampl(i) = 0
                end if
            end do
         end select
        !case(2) ..new ampl to be added       

  end function load_ampl


end module bc_dynflt_load
