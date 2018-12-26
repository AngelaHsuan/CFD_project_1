!---------------------------------------------------------
! ADI
subroutine ADI(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
    implicit NONE
    
    integer :: j,k,m,n,Ny,Nz,N_iter,i
    real :: delta_y,delta_z,err,u_errMAX,w,conv
    real*8, dimension(Ny) :: br,C_pr,C_Nr,C_Sr,C_Er,C_Wr,u_newr
    real*8, dimension(Nz) :: bc,C_pc,C_Nc,C_Sc,C_Ec,C_Wc,u_newc
    real*8, dimension(Ny,Nz) :: u_old,u_new,u_err
    character :: circle
        
    
    u_old = 0.d0
    u_errMAX = 1.d0
    N_iter = 0.d0
    circle = 'o'
    
    do while (u_errMAX > conv)
        if (circle == 'o') then
            
                        
            !開始row運算--------------
            do k = 1,Nz,1
        
                if (k == 1) then
                    do i = 1,Ny,1
                        u_new(i,k) = 0.d0
                    
                    end do
                
                else if (k == Nz) then
                    do j = 1,Ny,1
                        if (j == 1) then        !Corner
                        
                            C_Pr(j) = -1/delta_z**2 - 1/delta_y**2
                            C_Er(j) = 1/delta_y**2
                            C_Sr(j) = -1/delta_z**2
                            br(j) = -1/2 + C_Sr(j)*u_old(j,k-1)
                    
                        else if (j == Ny) then  !右上
                            C_Wr(j) = 0
                            C_Pr(j) = 1
                            br(j) = 0
                            u_new(j,k) = 0
                        else                    !Top
                            C_Wr(j) = 1/delta_y**2
                            C_Pr(j) = -2/delta_z**2 - 2/delta_y**2
                            C_Er(j) = 1/delta_y**2
                            C_Sr(j) = -2/delta_z**2
                            br(j) = -1+C_Sr(j)*u_old(j,k-1)
                    
                        end if
                
                    end do
                    call TDMA(Ny,C_Wr,C_Pr,C_Er,br,u_newr)
                    do i = 1,Ny,1
                        u_new(i,k) = u_newr(i)
                        u_err(i,k) = abs(u_new(i,k)-u_old(i,k))
                        u_old(i,k) = w*u_new(i,k) + (1-w)*u_old(i,k)
                    end do
                else
                    do j = 1,Ny,1
                        if (j == 1) then        !Left
                        
                            C_Pr(j) = -2/delta_z**2 - 2/delta_y**2
                            C_Er(j) = 2/delta_y**2
                            C_Nr(j) = -1/delta_z**2
                            C_Sr(j) = -1/delta_z**2
                            br(j) = -1 + C_Sr(j)*u_old(j,k-1) + C_Nr(j)*u_old(j,k+1)
                    
                        else if (j == Ny) then  !Right
                            C_Wr(j) = 0
                            C_Pr(j) = 1
                            br(j) = 0
                            u_new(j,k) = 0
                        else                    !Center
                            C_Wr(j) = 1/delta_y**2
                            C_Pr(j) = -2/delta_y**2 - 2/delta_z**2
                            C_Er(j) = 1/delta_y**2
                            C_Nr(j) = -1/delta_z**2
                            C_Sr(j) = -1/delta_z**2
                            br(j) = -1 + C_Sr(j)*u_old(j,k-1) + C_Nr(j)*u_old(j,k+1)
                        
                        end if
                
                    end do
                    call TDMA(Ny,C_Wr,C_Pr,C_Er,br,u_newr)
                    do i = 1,Ny,1
                        u_new(i,k) = u_newr(i)
                    end do
                end if
        
            end do

            circle = 'e'
    
        else
            
            !開始col運算--------------
            do j = 1,Ny,1
        
                if (j == 1) then
                    do k = 1,Nz,1
                        if (k == 1) then
                            C_Nc(k) = 0
                            C_Pc(k) = 1
                            C_Sc(k) = 0
                            C_Ec(k) = 0
                            bc(k) = 0
                            u_new(j,k) = 0
                        
                        else if (k == Nz) then      !Corner
                            C_Nc(k) = 0
                            C_Pc(k) = -1/delta_z**2 - 1/delta_y**2
                            C_Sc(k) = 1/delta_z**2
                            C_Ec(k) = - 1/delta_y**2
                            C_Wc(k) = 0
                            bc(k) = -1/2 + C_Ec(k)*u_old(j+1,k)
                        
                        else                        !Left
                            C_Nc(k) = 1/delta_z**2
                            C_Pc(k) = -2/delta_z**2 - 2/delta_y**2
                            C_Sc(k) = 1/delta_z**2
                            C_Ec(k) = - 2/delta_y**2
                            C_Wc(k) = 0
                            bc(k) = -1 + C_Ec(k)*u_old(j+1,k)
                        
                        end if
                    end do
                
                else if (j == Ny) then
                    do k = 1,Nz,1
                        C_Nc(k) = 0
                        C_Pc(k) = 1
                        C_Sc(k) = 0
                        bc(k) = 0
                        u_new(j,k) = 0
                    end do
                
                else
                    do k = 1,Nz,1
                        if (k == 1) then
                            C_Nc(k) = 0
                            C_Pc(k) = 1
                            C_Sc(k) = 0
                            bc(k) = 0
                            u_new(j,k) = 0
                    
                        else if (k == Nz) then      !Top
                            C_Nc(k) = 0
                            C_Pc(k) = -2.0/delta_z**2 - 2.0/delta_y**2
                            C_Sc(k) = 2.0/delta_z**2
                            C_Ec(k) = - 1/delta_y**2
                            C_Wc(k) = - 1/delta_y**2
                            bc(k) = -1 + C_Ec(k)*u_old(j+1,k) + C_Wc(k)*u_old(j-1,k)
                        
                                                    
                        else                        !Center
                            C_Nc(k) = 1/delta_z**2
                            C_Pc(k) = -2/delta_z**2 - 2/delta_y**2
                            C_Sc(k) = 1/delta_z**2
                            C_Ec(k) = - 1/delta_y**2
                            C_Wc(k) = - 1/delta_y**2
                            bc(k) = -1 + C_Ec(k)*u_old(j+1,k) + C_Wc(k)*u_old(j-1,k)
                        
                        end if
                
                    end do
                
                end if
            
            
                call TDMA(Nz,C_Sc,C_Pc,C_Nc,bc,u_newc)
            
                do i = 1,Nz,1
                    u_new(j,i) = u_newc(i)
                    u_err(j,i) = abs(u_new(j,i)-u_old(j,i))
                    u_old(j,i) = w*u_new(j,i) + (1-w)*u_old(j,i)
                end do
            
            
            end do
    
    
    
        
            circle = 'o'
        end if
    
    u_errMax = maxval(u_err)
    
    write(*,*) 'nember of iteration', N_iter,'    error    ', u_errMax
    
    N_iter = N_iter + 1
    
    end do
    
end subroutine ADI