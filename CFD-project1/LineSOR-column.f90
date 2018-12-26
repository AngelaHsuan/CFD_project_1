!---------------------------------------------------------------
! Line SOR of column
subroutine LineSORcol(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
    implicit NONE
    
    integer :: j,k,m,n,Ny,Nz,N_iter,i
    real :: delta_y,delta_z,err,u_errMAX,w,conv
    real*8, dimension(Nz) :: b,C_p,C_N,C_S,C_E,C_W,u_fnew
    real*8, dimension(Ny,Nz) :: u_old,u_new,u_err
    
    u_old = 0.d0
    u_errMAX = 1.d0
    N_iter = 0.d0
    
    do while (u_errMAX > conv)
        do j = 1,Ny,1
        
            if (j == 1) then
                do k = 1,Nz,1
                    if (k == 1) then
                        C_N(k) = 0
                        C_P(k) = 1
                        C_S(k) = 0
                        C_E(k) = 0
                        b(k) = 0
                        u_new(j,k) = 0
                        
                    else if (k == Nz) then      !Corner
                        C_N(k) = 0
                        C_P(k) = -1/delta_z**2 - 1/delta_y**2
                        C_S(k) = 1/delta_z**2
                        C_E(k) = - 1/delta_y**2
                        C_W(k) = 0
                        b(k) = -1/2 + C_E(k)*u_old(j+1,k)
                        
                    else                        !Left
                        C_N(k) = 1/delta_z**2
                        C_P(k) = -2/delta_z**2 - 2/delta_y**2
                        C_S(k) = 1/delta_z**2
                        C_E(k) = - 2/delta_y**2
                        C_W(k) = 0
                        b(k) = -1 + C_E(k)*u_old(j+1,k)
                        
                    end if
                end do
                
            else if (j == Ny) then
                do k = 1,Nz,1
                    C_N(k) = 0
                    C_P(k) = 1
                    C_S(k) = 0
                    b(k) = 0
                    u_new(j,k) = 0
                end do
                
            else
                do k = 1,Nz,1
                    if (k == 1) then
                        C_N(k) = 0
                        C_P(k) = 1
                        C_S(k) = 0
                        b(k) = 0
                        u_new(j,k) = 0
                    
                    else if (k == Nz) then      !Top
                        C_N(k) = 0
                        C_P(k) = -2.0/delta_z**2 - 2.0/delta_y**2
                        C_S(k) = 2.0/delta_z**2
                        C_E(k) = - 1/delta_y**2
                        C_W(k) = - 1/delta_y**2
                        b(k) = -1 + C_E(k)*u_old(j+1,k) + C_W(k)*u_old(j-1,k)
                        
                        !write(*,*)'hi'
                        !pause
                        
                    else                        !Center
                        C_N(k) = 1/delta_z**2
                        C_P(k) = -2/delta_z**2 - 2/delta_y**2
                        C_S(k) = 1/delta_z**2
                        C_E(k) = - 1/delta_y**2
                        C_W(k) = - 1/delta_y**2
                        b(k) = -1 + C_E(k)*u_old(j+1,k) + C_W(k)*u_old(j-1,k)
                        
                    end if
                
                end do
                
            end if
            
            
            call TDMA(Nz,C_S,C_P,C_N,b,u_fnew)
            
            do i = 1,Nz,1
                u_new(j,i) = u_fnew(i)
                u_err(j,i) = abs(u_new(j,i)-u_old(j,i))
                u_old(j,i) = w*u_new(j,i) + (1-w)*u_old(j,i)
            end do
            
            
        end do
        
        
        
    u_errMax = maxval(u_err)
    
    write(*,*) 'number of iteration', N_iter,'    error    ', u_errMax
    u_old = w*u_new + (1-w)*u_old
    N_iter = N_iter + 1
        
    end do
    
    
end subroutine LineSORcol