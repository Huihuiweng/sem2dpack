module distribution_boxer

  implicit none
  private
 
  type boxer_dist_type
    private
    integer :: NumZon ! number of zones
    double precision :: CenPnt   &! centraL point
                       ,InBox_W  &! interior box half-width
                       ,Tran_w     ! increment at Transition zone
    double precision :: Inter, Incre  ! Value of inside and increment at transition zone

  end type boxer_dist_type

  public :: boxer_dist_type,read_boxer_dist,generate_boxer_dist ,&
            destroy_boxer_dist

  contains

!
!***********************************************************************
!
! BEGIN INPUT BLOCK
!
! NAME   : DIST_BOXER
! GROUP  : DISTRIBUTIONS_2D
! PURPOSE: Piecewise constant radial (2D) distribution.
!          This distribution defines a  box zones with three zones: 
!          Interior zone        Transition zone          exterior zone
!          Transition zone distribution follow the boxer function, Interior/exterior zone keep constant value.
! SYNTAX : &DIST_PWCONR central_at, W, w /
!
! ARG: x0            [dble] [0d0] Central of the box zones
! ARG: W             [dble] [0d0] Half length of the interior box zone
! ARG: w             [dble] [0d0] Transition zone width
!                    first zone x <= W, 
!                    second W < x <= W + w 
!                    last x > W + w 
! ARG: Inter             [dble(num)] [none] Value inside              
! ARG: Incre         [dble][0d0] Increment of transition zone
! END INPUT BLOCK

  subroutine read_boxer_dist (data, file)

  use stdio, only: IO_abort
  type(boxer_dist_type), intent(out) :: data
  integer, intent(in) :: file

  double precision :: x0, W, T, Inter, Incre
  NAMELIST / DIST_BOXER / x0, W, T, Inter, Incre

  x0 = 0
  read (file,DIST_BOXER)
  data%CenPnt = x0
  data%InBox_W = W
  data%Tran_w = T
  data%Inter = Inter
  data%Incre = Incre
 
 
  end subroutine read_boxer_dist
!
!***********************************************************************
!

  subroutine generate_boxer_dist(field, coord, par)

  double precision, intent(in), dimension(:,:) :: coord
  type(boxer_dist_type), intent(in) :: par

  double precision, intent(out), dimension(:) :: field 
  double precision :: rad_x, tmp, B
  integer :: i,izone

  do i =1,size(field)
    rad_x =  abs(coord(1,i)-par%CenPnt) 
    if (rad_x <= par%Inbox_W ) then
                field(i) = par%Inter
    else if(par%InBox_W < rad_x .and. rad_x < par%InBox_W + par%Tran_w) then
                tmp = par%Tran_w / (rad_x - par%InBox_W - par%Tran_w) + par%Tran_w / (rad_x - par%InBox_W)
                B = 0.5 * (1 + tanh(tmp))
                field(i) = par%Inter + par%Incre * (1 - B) 
    else
                field(i) = par%Inter + par%Incre
    endif
  end do
 
  end subroutine generate_boxer_dist

!
!***********************************************************************
!
!! boxer_dist_type destructor
subroutine destroy_boxer_dist(d)
  type(boxer_dist_type), pointer :: d 
  deallocate(d)
end subroutine destroy_boxer_dist

end module distribution_boxer
