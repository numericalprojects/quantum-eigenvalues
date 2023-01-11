!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Constructs the diagonal and off diagonal elements of the hamiltonian
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! construct_T, construct_square_well_V, construct_QHO_V, construct_woods_V,  
!! construct_hamiltonian
!!----------------------------------------------------------------------
module hamiltonian
use types
implicit none
real(dp), parameter :: hbar = 197.3_dp, mass = 939.0_dp, v_0 = 50.0_dp, a = 0.2e-15_dp
! This is a good place to define subroutines and functions to calculate 
! the different contributions to the hamiltonian; the diagonal and 
! off-diagonal elements of the kinetic energy, and the diagonal elements
! for the different potentials

! Since more than one subroutine in this module will use the value for
! hbar and mass, it would be a good idea to define them here as 
! parameters


private
public :: construct_T, construct_square_well_V, construct_QHO_V, construct_woods_V, & 
construct_hamiltonian  
contains


!-----------------------------------------------------------------------
!! Subroutine: construct_T
!-----------------------------------------------------------------------
!! by: Nathan Crawford 
!!
!! Constructs the diagonal and off diagonal kinetic energy terms
!!----------------------------------------------------------------------
!! Input: n_points      number of intervals between the edges of the 
!!                      potential
!!        dx            x increment value between each interval
!!
!!
!-----------------------------------------------------------------------
!! Output: kinetic_diag     diagonal elements of kinetic energy 
!!         kinetic_offdiag  off diagonal elements of kinetic energy
!!
!! 
!-----------------------------------------------------------------------
subroutine construct_T(dx, n_points, kinetic_diag, kinetic_offdiag) 
implicit none 
real(dp), allocatable, intent(out) :: kinetic_diag(:), kinetic_offdiag(:) 
real(dp), intent(in) :: dx 
integer, intent(in) :: n_points 
integer :: n_offdiag, i 

!Let's declare a variable to hold the number of offdiagonal 
!elements in a matrix which should be the number of diagonal elements 
!minus 1. 

n_offdiag = n_points - 1

allocate(kinetic_diag(1:n_points)) 
allocate(kinetic_offdiag(1:n_offdiag)) 

!Contribution to diagonal elements of kinetic energy operator is 
!hbar^2/mass*dx^2 and the contribution to the offdiagonal is 
!-1/2 * hbar^2/mass*dx^2 

do i = 1, n_points 
    kinetic_diag(i) = hbar**2/(mass * (dx**2)) 
end do 

do i = 1, n_offdiag 
    kinetic_offdiag(i) = -0.5_dp * (hbar**2)/(mass * (dx**2)) 
end do 

end subroutine construct_T 

!-----------------------------------------------------------------------
!! Subroutine: construct_square_well_V
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Constructs the potential for the infinite square well
!!----------------------------------------------------------------------
!! Input: n_points      number of interval points from -L to L
!!
!!
!-----------------------------------------------------------------------
!! Output: potential    array containing values for the potential from 
!!                      -L to L. In this case all 0. 
!!
!-----------------------------------------------------------------------
subroutine construct_square_well_V(n_points, potential) 
implicit none 
real(dp), allocatable, intent(out) :: potential(:) 
integer :: i 
integer, intent(in) :: n_points

!Potential of infinite square well is 0
allocate(potential(1:n_points))
do i = 1, n_points 
    potential(i) = 0.0_dp 
end do 
end subroutine construct_square_well_V

!-----------------------------------------------------------------------
!! Subroutine: Construct_QHO_V
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Constructs potential for the QHO
!!----------------------------------------------------------------------
!! Input: n_points      number of interval points from each edge of the box
!!        x             array containing the sampled points from -L to L
!!
!-----------------------------------------------------------------------
!! Output: potential    array containing values of the potential from -L to L
!!
!! 
!-----------------------------------------------------------------------
subroutine construct_QHO_V(x, n_points, potential)  
implicit none  
real(dp), allocatable, intent(out) :: potential(:) 
real(dp), intent(in) :: x(:)
integer :: i 
integer, intent(in) :: n_points  

!Potential of QHO is hbar^2 * x^2/2*mass
allocate(potential(1:n_points))
do i = 1, n_points 
    potential(i) = (hbar**2) * (x(i)**2)/(2 * mass) 
end do 
end subroutine construct_QHO_V 

!-----------------------------------------------------------------------
!! Subroutine: construct_woods_V
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Constructs the potential for the woods saxon box
!!----------------------------------------------------------------------
!! Input: radius        radius of the woods saxon potential
!!        x             sampled points between each edge of the potential
!!        n_points      number of interval points between -L and L
!-----------------------------------------------------------------------
!! Output: potential    array containing the different values
!!
!! 
!-----------------------------------------------------------------------
subroutine construct_woods_V(x, n_points, radius, potential) 
implicit none 
real(dp), allocatable, intent(out) :: potential(:) 
real(dp), intent(in) :: x(:), radius
integer, intent(in) :: n_points 
integer :: i 

!Potential of the woods saxon potential is 
!-V_0/1+e^|x-R|/a
allocate(potential(1:n_points)) 

do i = 1, n_points 
    potential(i) = -v_0/(1 + exp((abs(x(i)) - radius)/a)) 
end do 
end subroutine construct_woods_V

!-----------------------------------------------------------------------
!! Subroutine: construct_hamiltonian
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Constructs hamiltonian operator
!!----------------------------------------------------------------------
!! Input: kinetic_diag      array containing diagonal elements of kinetic 
!!                          energy
!!        potential         array containing diagonal elements of potential 
!!                          energy
!!
!-----------------------------------------------------------------------
!! Output: hamiltonian      array containing diagonal elements of hamiltonian 
!!                          operator
!!
!! 
!-----------------------------------------------------------------------
subroutine construct_hamiltonian(kinetic_diag, potential, n_points, hamiltonian) 
implicit none 
real(dp), allocatable, intent(out) :: hamiltonian(:) 
real(dp) :: kinetic_diag(:), potential(:) 
integer :: n_points, i  
    
    !To optimize this slighty we know the potential energy 
    !operator is 0 everywhere except the diagonal elements 
    !So when adding the 2 matrices we only add the diagonals 
    
    allocate(hamiltonian(1:n_points))
    do i = 1, n_points 
        hamiltonian(i) = kinetic_diag(i) + potential(i) 
    end do 

end subroutine construct_hamiltonian


    
end module hamiltonian