!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! Nathan Crawford
!!
!! Contains subroutines and functions related to reading input from the
!! user and  writing output into a text file
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types 
use qm_solver, only : solve_woods_saxon
implicit none

private
public :: read_input, write_probability_density, print_energies, write_woods_saxon_energies

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! Nathan Crawford
!!
!! Displays a message describing what the program does and the expected
!! input. After that it uses the `read_real` and `read_integer`
!! functions to assign values to the different parameters.
!!----------------------------------------------------------------------
!! Output:
!!
!! n_points     integer     number of grid points the discretized wave function
!! length       real        length of the box
!! radius       real        radius of the Woods-Saxon potential
!-----------------------------------------------------------------------
subroutine read_input(length, radius, n_points)
    implicit none 
    integer :: n_points 
    real(dp) :: length, radius 
    print *, "Greetings user!" 
    print * 
    print *, "This program will solve the Schrodinger Equation for the " 
    print *, "Infinite square well, quantum harmonic oscillator and " 
    print *, "the Woods Saxon potential." 
    print *
    print *, "In order to do this we need your input as to the values of" 
    print *, "The number of points discretized between the 2 edges of the" 
    print *, "potentials, the length of the potential measured from 0 to one edge" 
    print *, "and the radius of the Woods Saxon Potential." 
    print *
    n_points = read_integer('Enter grid points of wave function n_points') 
    length = read_real('Enter length of the box of infinite square well') 
    radius = read_real('Enter radius of the Woods-Saxon potential') 
    

end subroutine read_input 

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Tests if user input was a real number
!!----------------------------------------------------------------------
!! Input: name      String that was printed
!!
!!
!-----------------------------------------------------------------------
!! Output: x        user input
!!
!! 
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none 
    
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print * 
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    ! Here we are using a do loop to check if the user entered a real number greater than 0. 
    ! If the user entered something non empty and if it is a real number greater than 0 then the loop 
    ! exits because the input is valid. 
    do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .AND. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a number or is a negative number, please provide a non negative number'
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo
    
end function read_real 

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Tests if user input is an integer
!!----------------------------------------------------------------------
!! Input: name      user input
!!
!!
!-----------------------------------------------------------------------
!! Output: x        user input value that is tested
!!
!! 
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len = 120) :: string
    integer :: ierror

    print *
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    ! Here we are using a do loop to check if the user entered an integer greater than 0. 
    ! If the user entered something non empty and if it is an integer greater than 0 then the loop 
    ! exits because the input is valid.

    ! Remember to declare above what type of variables string and ierror are

    do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .and. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not an integer or is a negative number, please provide a non negative number'
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!! Subroutine: write_probability_density
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! writes probability densities as a function of x in a file
!!----------------------------------------------------------------------
!! Input: x_vector      array containing sampled points from -L to L
!!        wave_functions 2-D array containing 3 lowest normalized wave 
!!                       functions
!!
!-----------------------------------------------------------------------
!! Output: No output. written in a file.
!!
!! 
!-----------------------------------------------------------------------
subroutine write_probability_density(file_name, x_vector, wave_functions)
    implicit none
    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: x_vector(:), wave_functions(:,:)
    integer :: unit, n, i  
    n = size(x_vector) 
    
    
    open(newunit=unit,file=trim(file_name))
    write(unit,'(4a28)') 'x', 'ground state', '1st excited', '2nd excited'
        
        !The probability density of a wave function value is that value squared, assuming 
        !the wave function is real.
        
        do i = 1, n 
            write(unit, *) x_vector(i), wave_functions(i, 1)**2, wave_functions(i, 2)**2, wave_functions(i, 3)**2 
        end do 
    close(unit)
    print *, "Probability densities were written in ", file_name
end subroutine write_probability_density

!-----------------------------------------------------------------------
!! Subroutine: print_energies
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Prints to screen numerical and analytical energies of ground state 
!! first excited and second excited
!!----------------------------------------------------------------------
!! Input:
!!       numerical      array containing 3 lowest energies calculated 
!!                      numerically
!!       analytic       array containing 3 lowest energies calculated 
!!                      analytically
!-----------------------------------------------------------------------
!! Output: No output. Printed to screen
!!
!! 
!-----------------------------------------------------------------------
subroutine print_energies(name, numerical, analytic)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: numerical(:), analytic(:) 
    integer :: i 

    if (size(numerical) /= size(analytic)) then
        print*, "arrays size don't match in print_energies"
        stop
    endif

    print*, 'Comparing numerical and anlytic solutions in'
    print*, trim(name)
    print*
    print'(a9,2a15)', 'number', 'numerical', 'analytic'
    
    !Let i be a particular eigenstate and the numerical is the 
    !numerical result and analytic is the analytic result.
    do i = 1, size(numerical) 
        print *, i, numerical(i), analytic(i)
    end do 
    
end subroutine print_energies

!-----------------------------------------------------------------------
!! Subroutine: write_woods_saxon_energies
!-----------------------------------------------------------------------
!! by: Nathan Crawford
!!
!! Writes in a file the 3 lowest woods saxon energies as a function of 
!! the radius from 2 to 10
!!----------------------------------------------------------------------
!! Input: r_min         starting radius of 2 
!!        r_max         ending radius of 10 
!!        energies      3 lowest energies of the woods saxon problem
!-----------------------------------------------------------------------
!! Output: No output. Writes output into a file.
!!
!! 
!-----------------------------------------------------------------------
subroutine write_woods_saxon_energies(file_name, n_points, dx, wave_functions, x_vector, energies, r_min, r_max)
    implicit none
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: n_points
    real(dp) :: r_min, r_max, dx  
    real(dp), allocatable :: wave_functions(:, :), lowest_wave(:, :) 
    real(dp) :: r_increment, r, energies(:), x_vector(:)
    integer :: unit 
    !Value for how much to increment r by from 2 to 10.
    r_increment = 0.05_dp  
    !r_min and r_max are parameter values. So let's define a new variable 
    !to hold r_min and all current or future radius values
    r = r_min
    open(newunit=unit,file=trim(file_name))
    write(unit, *) 'radius ', 'ground state energy ', 'first excited state ', 'second excited state' 
    do 
    !We have to solve woods saxon for many different values of radius so lets do that in an 
    !infinite do loop and exit once r is at the maximum
    call solve_woods_saxon(x_vector, dx, wave_functions, n_points, energies, r, lowest_wave) 
    
    write(unit, *) r, energies(1), energies(2), energies(3) 
    
        if(r >= r_max) then 
            exit 
        end if 
   
    r = r + r_increment 
    
    end do 
    close(unit)
    print *, "Woods Saxon ground state, first excited state, second excited state energies written in ", file_name
end subroutine write_woods_saxon_energies


end module read_write
