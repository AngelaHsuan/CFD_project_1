# CFD_project_1
10/25/2017 written

Trying to use numerical method to solve the poisson eqqation in the channel flow which only suffers from the pressure gradient.
"A" for the cross section aspect ratio, "w" for relaxation coefficient, iteration method is the initail condition given by user.

file: CFD-project1.f90 is the main program
file: PointSOR.f90 is subroutine using point SOR
file: LineSOR-row.f90 is subroutine using line SOR
file: LineSOR-column.f90 is subroutine using line SOR
file: ADI.f90 is subroutine using ADI
file: printout.f90 is subroutine for information output
file: TDMA.f90 is subroutine for TDMA
