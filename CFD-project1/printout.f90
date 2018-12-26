!------------------------------------------------------------
! Result output
subroutine printout(delta_y,delta_z,Ny,Nz,u_new,printoutPLA)
    implicit NONE
    
    integer :: j,k,Ny,Nz,Pad
    real :: delta_z,delta_y
    real*8, dimension(Ny,Nz) :: u_new
    character :: printoutPLA
    
    Pad = 100
    if (printoutPLA == 'P') then
        open (unit=Pad, file="PSOR.txt", status="UNKNOWN")
    else if (printoutPLA == 'L') then
        open (unit=Pad, file="LSOR.txt", status="UNKNOWN")
    else
        open (unit=Pad, file="ADI.txt", status="UNKNOWN")
    end if
    
    write(Pad,*) 'variables="x","y","Velocity"'
    write(Pad,*) 'zone j=', Ny,'k=',Nz,'DATAPACKING=POINT'
    
    do k = 1,Nz,1
        do j = 1,Ny,1
            write(Pad,'(3F20.12)') 0.d0+(j-1)*delta_y, 0.d0+(k-1)*delta_z, u_new(j,k)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine printout