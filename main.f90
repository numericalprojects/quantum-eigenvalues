!-----------------------------------------------------------------------------
!! Program: schrodinger_solution
!! By: Nathan Crawford 
!!
!! This program solves the time independent Schrodinger Equation for the 3 
!! cases of the infinite square well, quantum harmonic oscillator(will be 
!! abbreviated QHO) and the Woods Saxon problem. 

!! It does this with the lapack subroutine dstev that takes the 
!! diagonals and off diagonals of the hamiltonian and returns 
!! the ordered eigenvalues and eigenvectors. 
!!
!! This program will then normalize the eigenvectors/wavefunctions 
!! and print to screen the analytical and numerical energies in the 
!! particle in a box and QHO case for comparison and write the probability 
!! densities in all 3 cases as a function of x. 
!! X is some point along the potential from -L to L where L is the length of 
!! the potential from the edge to 0. 
!!
!! This program will then write the ground state, first excited state and second 
!! excited state energies as a function of the Woods Saxon potential radius R. 
!-----------------------------------------------------------------------------
program schrodinger_solution 

use types
use read_write, only : read_input, write_probability_density, print_energies, write_woods_saxon_energies
use qm_solver, only: sample_box, solve_infinite_well, analytic_infinite_well, solve_harmonic_oscillator,&
    analytic_harmonic_oscillator, solve_woods_saxon
implicit none

integer :: n_points
real(dp) :: length, radius, dx

integer, parameter :: n_energies = 3
real(dp) :: energies(1:n_energies), analytic_energies(1:n_energies)
real(dp), allocatable :: wave_functions(:,:), lowest_wave(:, :)
real(dp), allocatable :: x_vector(:)
real(dp), parameter :: r_min = 2._dp, r_max = 10._dp

call read_input(length, radius, n_points)

call sample_box(n_points, x_vector, length, dx)

! Solving particle in a box
call solve_infinite_well(dx, wave_functions, n_points, energies, lowest_wave) 

call analytic_infinite_well(analytic_energies, length)
call print_energies('Infinite Well', energies, analytic_energies) 
print * 
print *, "Writing probability densities as a function of x in infinite_well_wf.dat"
call write_probability_density('infinite_well_wf.dat', x_vector, lowest_wave)

! Solving harmonic oscillator 

call solve_harmonic_oscillator(x_vector, dx, wave_functions, n_points, energies, lowest_wave) 

call analytic_harmonic_oscillator(analytic_energies)
call print_energies('Harmonic oscillator', energies, analytic_energies) 
print * 
print *, "Writing probability densities as a function of x in harmonic_oscillator_wf.dat" 
call write_probability_density('harmonic_oscillator_wf.dat', x_vector, lowest_wave)

! Solving Woods Saxon 

call solve_woods_saxon(x_vector, dx, wave_functions, n_points, energies, radius, lowest_wave)  
print * 
print *, "Writing probability densities as a function of x in woods_saxon_wf.dat"
call write_probability_density('woods_saxon_wf.dat', x_vector, lowest_wave)

! Woods Saxon Energies as a function of radius 
print * 
print *, "Writing 3 lowest energies of woods saxon as a function of radius in woods_saxon_ener.dat"
call write_woods_saxon_energies('woods_saxon_ener.dat', n_points, dx, wave_functions, x_vector, energies, r_min, r_max)

end program schrodinger_solution