! verlet_list_module.f90
! Verlet list handling routines for MD simulation
MODULE verlet_list_module

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor

  IMPLICIT NONE
  PRIVATE

  ! Public routines
  PUBLIC :: initialize_list, finalize_list, make_list

  ! The initialize_list routine sets the value of r_list_factor. If you wish to read it
  ! from standard input using a NAMELIST nml_list, just uncomment the indicated statements
  ! It is assumed that all positions and displacements are divided by box
  ! r_list_box is set to r_cut_box*r_list_factor

  ! Public data
  INTEGER*8, DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: point ! index to neighbour list (n)
  INTEGER*8, DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: list  ! Verlet neighbour list (nl)


  ! Private data
  INTEGER*8                           :: nl         ! Size of list
  REAL*8                              :: r_list     ! List range parameter
  REAL*8                              :: r_skin     ! List skin parameter
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: r_save     ! Saved positions for list (3,n)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: dr         ! Displacements (3,n)

CONTAINS

  SUBROUTINE initialize_list (n,r_cut)
    IMPLICIT NONE
    INTEGER*8, INTENT(in) :: n         ! number of particles
    REAL*8, INTENT(in) :: r_cut         ! cut off distance
    
    REAL*8    :: r_list_factor

    REAL*8, PARAMETER :: pi = 4.0d0*DATAN(1.0d0)

    ! Uncomment the following statements if you wish to read in the value of r_list_factor
!!$    INTEGER :: ioerr
!!$    NAMELIST /nml_list/ r_list_factor

    ! Sensible default for r_list_factor
    r_list_factor = 1.2d0

    ! Uncomment the following statements if you wish to read in the value of r_list_factor
!!$    READ ( unit=input_unit, nml=nml_list, iostat=ioerr ) ! namelist input
!!$    IF ( ioerr /= 0 ) THEN
!!$       WRITE ( unit=error_unit, fmt='(a,i15)' ) 'Error reading namelist nml_list from standard input', ioerr
!!$       IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
!!$       IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
!!$       STOP 'Error in initialize_list'
!!$    END IF

!    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list factor = ', r_list_factor
    IF ( r_list_factor <= 1.0d0 ) THEN
       WRITE ( unit=error_unit, fmt='(a,f15.6)') 'r_list_factor must be > 1', r_list_factor
       STOP 'Error in initialize_list'
    END IF
    r_list = r_cut * r_list_factor
    r_skin = r_list - r_cut
!    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list range  = ', r_list
!    WRITE ( unit=output_unit, fmt='(a,t40,f15.6)') 'Verlet list skin = ', r_skin
    
    ! Estimate list size based on density + 10 per cent
!    nl = CEILING ( 1.1d0*(4.0d0*pi/3.0d0)*(r_list**3.0d0)*DBLE(n*(n-1)) / 2.0d0 ,kind=8)
    nl = n*30

    WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'Verlet list size = ', nl

    ALLOCATE ( r_save(3,n), dr(3,n), point(n), list(nl) )

  END SUBROUTINE initialize_list

  SUBROUTINE finalize_list
    IMPLICIT NONE

    DEALLOCATE ( r_save, dr, point, list )

  END SUBROUTINE finalize_list

  SUBROUTINE resize_list ! reallocates list array, somewhat larger
    IMPLICIT NONE
    INTEGER*8, DIMENSION(:), ALLOCATABLE :: tmp
    INTEGER*8                            :: nl_new

    nl_new = CEILING ( 1.25d0 * DBLE(nl),kind=8 )
    WRITE( unit=output_unit, fmt='(a)', advance='no' ) 'Warning: new Verlet list array size = '

    ALLOCATE ( tmp(nl_new) ) ! new size for list
    tmp(1:nl) = list(:)      ! copy elements across

    CALL MOVE_ALLOC ( tmp, list )
    nl = SIZE(list)
    WRITE( unit=error_unit, fmt='(t60,i15)') nl

  END SUBROUTINE resize_list

  SUBROUTINE make_list ( n, r, list_flag )
    IMPLICIT NONE
    INTEGER*8,                 INTENT(in) :: n
    REAL*8,    DIMENSION(3,n), INTENT(in) :: r
    LOGICAL, OPTIONAL                     :: list_flag

    INTEGER*8            :: i, j, k
    REAL*8               :: r_list_sq, rij_sq, dr_sq_max
    REAL*8, DIMENSION(3) :: rij

    LOGICAL, SAVE :: first_call = .TRUE.
    
    IF(PRESENT(list_flag))THEN
        first_call = .TRUE.
    END IF
    
    IF ( .NOT. first_call ) THEN
        IF(PRESENT(list_flag)) THEN
            OPEN(unit=10,file='.\Output\check.dat',status='old',position='append')
            WRITE(unit=10,fmt='(a)')'Shouldn''t be here!'
            WRITE(10,*)point
            WRITE(unit=10,fmt='(a)')'End Shouldn''t be here'
            CLOSE(10)
        END IF
       dr = r - r_save                             ! Displacement since last list update
       dr_sq_max = MAXVAL ( SUM(dr**2.0d0,dim=1) )     ! Squared maximum displacement
       IF ( 4.0d0*dr_sq_max < r_skin**2.0d0 ) RETURN ! No need to make list
    END IF
    
    first_call = .FALSE.

    k = 0
    r_list_sq = r_list ** 2.0d0

    
    DO i = 1, n - 1 ! Begin outer loop over atoms

       point(i) = k + 1

       DO j = i + 1, n ! Begin inner loop over partner atoms

          rij(:) = r(:,i) - r(:,j)
          rij_sq = SUM ( rij**2.0d0 )
          
          IF ( rij_sq < r_list_sq ) THEN

             k = k + 1
             IF ( k > nl ) CALL resize_list
             list(k) = j

          END IF

       END DO ! End inner loop over partner atoms

    END DO ! End outer loop over atoms
    
    point(n) = k + 1

    r_save = r
    
  END SUBROUTINE make_list

END MODULE verlet_list_module
