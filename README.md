## Goal of the Program 
 The goal of this program is to find numerical solutions to the eigenvalue problem, the time-independent Schrodinger Equation. 
 We will work in 1 dimension so the equation becomes: 
 $[\frac{-ℏ^2}{2m}\frac{d^2}{dx^2} + V(x)]Ψ(x) = EΨ(x)$ 
 
 We want to solve this for different potentials V(x). The potentials will be infinite square well, harmonic oscillator 
 and Woods Saxon potential. We can compare analytically for infinite square well and the harmonic oscillator but 
 the Woods Saxon has no analytical solution and requires numerical methods. 
 
We will discretize the wavefunction  on a lattice, turning it into a column vector. On the same lattice, the potential 
energy will be a diagonal matrix, while the kinetic energy will be a tridiagonal matrix. Therefore the Hamiltonian will 
be another tridiagonal matrix. We can then use the lapack library to find the eigenvalues (energies) and eigenvectors (wavefunctions). 
If you can't install lapack on your system you can use the tqli and eigsrt subroutines provided in the code.

Remember the ultimate goal is to find the eigenenergies and eigenfunctions for the Woods Saxon potential: 
$V(x)= \frac{-V_0}{1 + exp((|x|-R)/a)}$ 
where $V_0$ is the potential's depth, $R$ is the radius, and $a$ is the surface thickness. We will let $ℏ = 197.3 MeV fm, c = 1$ and 
$m = 939MeV$ (this would be $m = 939MeV/c^2$ normally but $c = 1$ ). Let $V_0 = 50 MeV$ and $a = 0.2 fm$ just for some typical values. 

The wavefunction will be discretized in a “box” from `–L` to `+L`, with `N`
points and  a step size `dx = (2L)/(N-1)`. Hence the sampled points for x will
be `x(i) = -L + dx*(i-1)`.

The kinetic energy is slightly subtler, but not much. Remember that it is  $\frac{-ℏ^2}{2m}\frac{d^2}{dx^2}Ψ(x)$ , and on a lattice of equally
spaced points we can discretize the second derivate with a symmetric
three-point formula: $(Ψ(x-dx) + Ψ(x+dx) - 2Ψ(x))/dx^2$  Hence the kinetic
energy contributes to the diagonal a constant `+hbar**2/mass/dx**2` and along
the off-diagonal also a constant `-0.5*hbar**2/mass/dx**2`.

Adding the potential and kinetic energy matrices and we have a discretized
representation of the Hamiltonian. 

With all that being said, these are the parts of the program. 

### First part. Infinite well 
 The first potential is the particle in a box scenario in which $V(x) = 0$. The program will construct the hamiltonian then calculate 
 the eigenvalues and eigenvectors(as well as making sure the eigenvalues are in ascending order). Then the program normalizes the eigenfunctions. 
 Then returns the 3 lowest energies and eigenfunctions. 
 
 The particle in a box scenario is a very well understood problem and is a common homework problem for any physics undergrad in a quantum 
 mechanics course. These have analytical closed form solutions that we can compare to. In fact if we were to plot these eigenfunctions we expect 
 sines and cosines. And we can compare the eigenvalues using this analytical formula: $E_i = \frac{i^2 ℏ^2 π^2}{8mL^2}$. 
 
 Our numerical solutions may not exactly match which is fine. We have a second order derivative with only a first order approximation after all. The  program will print the 3 lowest analytical and numerical energies on the screen then write in a file the three normalized probability densities and the sampled points. 
 
 ### Second part. Harmonic Oscillator 
 The next part does the same thing but with the quantum harmonic oscillator where $V(x) = \frac{ℏ^2 x^2}{2m}$. This is also a well understood problem with 
 an analytical solution(granted it's a little harder). The eigenfunctions here should look like Gaussians. We can compare the eigenvalues to the analytical 
 expression: 
 $E_i = (i- \frac{1}{2})ℏ^2/m$. The program will do the same things as before. 
 
 ### Third part. Woods Saxon 
 The last part of this program is the Woods Saxon potential $V(x)= \frac{-V_0}{1 + exp((|x|-R)/a)}$. The program will do the same thing as the previous 2, 
 but this time will also write the 3 lowest energies as a function of $R$ going from 2 to 10. 
 
 ### Jupyter Notebook 
 Once the program is finished we can head over to jupyter and there should be 4 different plots: 
 1. The ground state normalized probability density for the 3 different problems
 2.The 1st excited state normalized probability density for the 3 different problems
 3.The 2nd excited state normalized probability density for the 3 different problems
 4.The three lowest energies for the Woods-Saxon potential as a function of the radius 
 
 ### Compiling the program 
 It would be quite beneficial to you if you had a Linux system because it would enable you to use the makefile included. 
If you don't then you'll need to adjust the source code itself to solve the eigenvalue problem.

If this is the case then what you do is open a terminal, use the cd command to change to this directory. 

Then type make. 

You'll see some gfortran commands being executed. All of this has created an exectuable file called woods_saxon

Then the program will ask for input. Traps have been set in case you accidentally enter invalid input. 

Results will be displayed on screen and written into several files. Now you can head over to Jupyter Notebook 
to analyze the results. 
