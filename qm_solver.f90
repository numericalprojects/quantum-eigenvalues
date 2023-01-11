!-----------------------------------------------------------------------
!Module: qm_solver
!-----------------------------------------------------------------------
!! By: Nathan Crawford 
!!
!! This module finds the 3 lowest energies and the corresponding wave functions 
!! of the particle in a box, quantum harmonic oscillator and woods saxon problems. 
!! It then normalizes the wave functions and returns the lowest energies and wave 
!! functions. 
!! 
!! The general solution process is setting up the hamiltonian operator and sending 
!! it to another subroutine to solve the eigenvalue problem. Documentation is provided 
!! in those modules
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_infinite_well, sample_box, analytic_infinite_well, solve_harmonic_oscillator, 
!! analytic_harmonic_oscillator, solve_woods_saxon
!!----------------------------------------------------------------------
module qm_solver
use types
use hamiltonian, only : construct_T, construct_square_well_V, construct_QHO_V, construct_woods_V, & 
construct_hamiltonian
use eigen_solver, only : solve_evp 
implicit none

private
public solve_infinite_well, sample_box, analytic_infinite_well, solve_harmonic_oscillator, &
    analytic_harmonic_oscillator, solve_woods_saxon

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Samples x coordinates from -L to L where L is the length of the box. 
!! Every potential is discretized in this range. 
!!----------------------------------------------------------------------
!! Input: n_points      Number of points from -L to L
!!        
!!        length        length of the box if the length is measured 
!!                      from 0 to one of the boundaries 
!!
!!        x             array that contains sampled points from -L to L  
!!
!!        dx            increment value between each interval point from -L 
!!                      to L
!-----------------------------------------------------------------------
!! Output: dx           how much x is incremented at each n_points interval 
!!                      from 0 to L 
!!         x            array that contains sampled points from -L to L 
!!
!-----------------------------------------------------------------------
subroutine sample_box(n_points, x, length, dx)
    implicit none 
    real(dp), allocatable, intent(out):: x(:) 
    real(dp), intent(out) :: dx 
    real(dp) :: length
    integer, intent(in) :: n_points
    integer :: i 
    !Since our potentials in this program are 
    !defined from -L to L we use that and the number 
    !of points to define our increment value. 
    dx = (2 * length)/(n_points - 1)
    !We want the potential to start at -L and end at L so let 
    !x at any point i be x(i) = -L + dx * (i-1)
    allocate(x(1:n_points)) 
    do i = 1, n_points 
        x(i) = -length + dx * (i - 1) 
    end do 
end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: solve_infinite_well
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Solves the infinite square well problem and returns the lowest energies 
!! and the corresponding normalized wave functions
!!----------------------------------------------------------------------
!! Input: n_points      Number of points from -L to L
!!        
!!        dx            increment value between each interval point from -L 
!!                      to L
!!        
!!        wave_functions 2 dimensional array to hold the eigenvectors 
!!                       from the solved eigenvalue problem
!!        
!!        lowest_energies array containing 3 lowest energies found 
!!                        from the time independent Schrodinger Equation
!-----------------------------------------------------------------------
!! Output: lowest_wave      3 normalized wave functions corresponding to 
!!                          the 3 lowest energies
!!         lowest_energies  3 lowest energies of the infinite square well
!!
!-----------------------------------------------------------------------
subroutine solve_infinite_well(dx, wave_functions, n_points, lowest_energies, lowest_wave)
    implicit none 
    real(dp), intent(in) :: dx 
    real(dp), intent(out) :: lowest_energies(:) 
    real(dp), allocatable, intent(out) :: wave_functions(:, :), lowest_wave(:, :)
    real(dp), allocatable :: kinetic_diag(:), kinetic_offdiag(:), potential(:), hamiltonian(:), energies(:) 
    integer :: n_points, i
    
    !Our process for any problem should be this 
    !Make the Kinetic energy and potential energy 
    !matrices and add them to construct the hamiltonian. 
    
    !Then solve the eigenvalue problem to get our wavefunctions 
    !and energies then normalize the wave functions and return the 
    !3 lowest energies and corresponding eigenvalues. 
    
    allocate(lowest_wave(1:n_points, 1:3)) 
    
    
    call construct_T(dx, n_points, kinetic_diag, kinetic_offdiag) 
  
    call construct_square_well_v(n_points, potential) 
    
    call construct_hamiltonian(kinetic_diag, potential, n_points, hamiltonian) 
   
    call solve_evp(hamiltonian, kinetic_offdiag, wave_functions, energies) 
    
    call normalize_wave_functions(wave_functions, dx, n_points) 
    
    do i = 1, size(lowest_energies) 
        lowest_energies(i) = energies(i) 
        lowest_wave(:, i) = wave_functions(:, i) 
    end do 
    
end subroutine solve_infinite_well 


!-----------------------------------------------------------------------
!! Subroutine: normalize_wave_functions
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Normalizes the wave functions 
!!----------------------------------------------------------------------
!! Input: n_points      Number of points from -L to L
!!        
!!        dx            increment value between each interval point from -L 
!!                      to L       
!-----------------------------------------------------------------------
!! Output: wave_functions  normalized wave functions  
!!
!-----------------------------------------------------------------------
subroutine normalize_wave_functions(wave_functions, dx, n_points) 
implicit none 
real(dp), intent(in) :: dx 
real(dp), intent(out) :: wave_functions(:, :) 
real(dp) :: summation 
real(dp), allocatable :: normalization_constant(:)  
integer, intent(in) :: n_points 
integer :: i, j 

!Each wave function will have a different normalization 
!constant(it doesn't for a particle in a box but we should be 
!more general) so lets allocate an array to hold all of them. 
allocate(normalization_constant(1:n_points)) 
do j = 1, n_points 
    summation = 0.0_dp
    do i = 1, n_points
        !Let's make a simple way to evaluate the normalization integral 
        !Let's assume the wave function is not complex and just square each 
        !evaluation and multiply it by a dx and sum it up.
        summation = summation + (wave_functions(i, j)**2 * dx)  
    end do 
    normalization_constant(j) = sqrt(1/summation) 
end do 

 
!Apply the normalization constant to each wave function to make them 
!normalized
do j = 1, n_points 
    do i = 1, n_points 
        wave_functions(i, j) = wave_functions(i, j) * normalization_constant(j)  
    end do 
end do 

end subroutine normalize_wave_functions

!-----------------------------------------------------------------------
!! Subroutine: analytic_infinite_well
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! calculates the analytic energies according to the analytic formula
!!----------------------------------------------------------------------
!! Input: length        length of the box
!!
!-----------------------------------------------------------------------
!! Output: analytic_energies        energy of eigenstate in MeV
!!
!-----------------------------------------------------------------------
subroutine analytic_infinite_well(analytic_energies, length)
    implicit none 
    real(dp) :: hbar, mass, length 
    real(dp), intent(out) :: analytic_energies(:) 
    integer :: i 
    
        hbar = 197.3_dp 
        mass = 939.0_dp 
        !this is more straightforward. Analytic solutions for the energy already exist 
        !So let's simply apply the formula. 
        do i = 1, size(analytic_energies) 
            analytic_energies(i) = (i**2.0_dp) * (hbar**2 * pi**2)/(8 * mass * length**2) 
        end do 
end subroutine analytic_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: solve_harmonic_oscillator
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Solves the quantum harmonic oscillator problem and returns the 
!! 3 lowest energies and corresponding wave functions
!!----------------------------------------------------------------------
!! Input: n_points      Number of points from -L to L
!!        
!!        dx            increment value between each interval point from -L 
!!                      to L
!!        
!!        wave_functions 2 dimensional array to hold the eigenvectors 
!!                       from the solved eigenvalue problem
!!        
!!        lowest_energies array containing 3 lowest energies found 
!!                        from the time independent Schrodinger Equation
!!
!-----------------------------------------------------------------------
!! Output: lowest_wave      3 normalized wave functions corresponding to 
!!                          the 3 lowest energies
!!         lowest_energies  3 lowest energies of the infinite square well
!!
!-----------------------------------------------------------------------
subroutine solve_harmonic_oscillator(x, dx, wave_functions, n_points, lowest_energies, lowest_wave)
    implicit none 
    real(dp), intent(in) :: x(:), dx 
    real(dp), allocatable, intent(out) :: lowest_wave(:, :)
    real(dp), intent(out) :: lowest_energies(:) 
    real(dp), allocatable :: kinetic_diag(:), kinetic_offdiag(:), potential(:), hamiltonian(:), wave_functions(:, :), energies(:)
    integer :: n_points, i 
    
    !Our process for any problem should be this 
    !Make the Kinetic energy and potential energy 
    !matrices and add them to construct the hamiltonian. 
    
    !Then solve the eigenvalue problem to get our wavefunctions 
    !and energies then normalize the wave functions and return the 
    !3 lowest energies and corresponding eigenvalues.  
    allocate(lowest_wave(1:n_points, 1:size(lowest_energies)))  
    
    call construct_T(dx, n_points, kinetic_diag, kinetic_offdiag) 
    call construct_qho_v(x, n_points, potential) 
    call construct_hamiltonian(kinetic_diag, potential, n_points, hamiltonian)
    call solve_evp(hamiltonian, kinetic_offdiag, wave_functions, energies) 
    call normalize_wave_functions(wave_functions, dx, n_points) 
    
    do i = 1, size(lowest_energies) 
        lowest_energies(i) = energies(i) 
        lowest_wave(:, i) = wave_functions(:, i) 
    end do 
end subroutine solve_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: analytic_harmonic_oscillator
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Calculates energies of quantum harmonic oscillator analytically
!!----------------------------------------------------------------------
!! Input: hbar      quantum mechanical constant equal to 197.3 MeV
!!        
!!        mass      mass of the particle equal to 939 MeV
!-----------------------------------------------------------------------
!! Output: analytic_energies    energies of eigenstates of QHO
!!
!-----------------------------------------------------------------------
subroutine analytic_harmonic_oscillator(analytic_energies)
    implicit none 
    real(dp), intent(out) :: analytic_energies(:) 
    real(dp) :: hbar, mass 
    integer :: i 
    hbar = 197.3_dp 
    mass = 939.0_dp 
    !same as before. We simply use a formula to evaluate energies 
    !analytically.
    do i = 1, size(analytic_energies) 
        analytic_energies(i) = (i - 0.5_dp) * (hbar**2)/mass 
    end do 
end subroutine analytic_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: solve_woods_saxon
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Solves the wood saxon problem and returns the normalized wave functions 
!! that correspond to the 3 lowest energies
!!----------------------------------------------------------------------
!! Input: n_points      Number of points from -L to L
!!        
!!        dx            increment value between each interval point from -L 
!!                      to L
!!        
!!        wave_functions 2 dimensional array to hold the eigenvectors 
!!                       from the solved eigenvalue problem
!!        
!!        lowest_energies array containing 3 lowest energies found 
!!                        from the time independent Schrodinger Equation
!!
!-----------------------------------------------------------------------
!! Output: lowest_wave      3 normalized wave functions corresponding to 
!!                          the 3 lowest energies
!!         lowest_energies  3 lowest energies of the infinite square well
!!
!!
!-----------------------------------------------------------------------
subroutine solve_woods_saxon(x, dx, wave_functions, n_points, lowest_energies, radius, lowest_wave)
    implicit none 
    real(dp), intent(in) :: x(:), dx 
    real(dp), allocatable, intent(out) :: lowest_wave(:, :) 
    real(dp), intent(out) :: lowest_energies(:)
    real(dp) :: radius
    real(dp), allocatable :: kinetic_diag(:), kinetic_offdiag(:), potential(:), hamiltonian(:),  wave_functions(:,:), energies(:) 
    integer :: n_points, i
    !Our process for any problem should be this 
    !Make the Kinetic energy and potential energy 
    !matrices and add them to construct the hamiltonian. 
    
    !Then solve the eigenvalue problem to get our wavefunctions 
    !and energies then normalize the wave functions and return the 
    !3 lowest energies and corresponding eigenvalues. 
    allocate(lowest_wave(1:n_points, 1:size(lowest_energies))) 
    
    
    
    call construct_T(dx, n_points, kinetic_diag, kinetic_offdiag) 
    call construct_woods_v(x, n_points, radius, potential) 
    call construct_hamiltonian(kinetic_diag, potential, n_points, hamiltonian)
    call solve_evp(hamiltonian, kinetic_offdiag, wave_functions, energies) 
    call normalize_wave_functions(wave_functions, dx, n_points) 
    
    do i = 1, size(lowest_energies) 
        lowest_energies(i) = energies(i) 
        lowest_wave(:, i) = wave_functions(:, i)
    end do 
end subroutine solve_woods_saxon

end module qm_solver