!+
! Subroutine tao_hook_lattice_calc (calc_ok)
!
! This hook is used to do custom lattice calculations before any of the standard calculations
! in the tao_lattice_calc routine.
!
! See tao/code/tao_lattice_calc.f90 for how the standard lattice 
! calculation is performed and how this routine is called.
!
! Also consider:
!   tao_hook_calc_calc_post_process   ! Called at the end of tao_lattice_calc's universe loop.
!   tao_hook_branch_calc
!
!
! Input:
!   calc_ok    -- logical: Current state of the lattice calculation.
!
! Output:
!   s%u(i)%universe_recalc 
!              -- logical: Set this to False to suppress tao_lattice_recalc.
!   calc_ok    -- logical: Set False if there was an error in the calculation like a particle 
!                   was lost or a lattice is unstable. Note to programmers: Do not set True.
!-

subroutine tao_hook_lattice_calc (calc_ok)

use tao_interface

implicit none

logical calc_ok

! 

end subroutine tao_hook_lattice_calc
