!---------------------------------------------------------------
subroutine Initial(N, del_x, del_y, Pe_m, conv, w,beta,G,printout, DCM)
    implicit NONE
    
    integer :: N, Nv, N_select, Pe_select, method_select, G, i, j,printout, DCM_select, DCM
    real*8 :: Pe_m, del_x, del_y, conv, w, beta, pi
    
    conv = 1D-6
    
    write(*,*)'How many grid points do you want?'
    write(*,*)'(1) 61 by 61'
    write(*,*)'(2) 81 by 81'
    read(*,*)N_select
    select case (N_select)
        case(1)
            N = 61
        case(2)
            N = 81
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
    
    del_x = 1.d0/(N-1)
    del_y = 1.d0/(N-1)
    
    
    
    !-------------------------------------------------
    write(*,*)'Choose the number you want for Pe?'
    write(*,*)'(1) 10'
    write(*,*)'(2) 100'
    write(*,*)'(3) 1000'
    write(*,*)'(4) 1500'
    read(*,*)Pe_select
    select case (Pe_select)
        case(1)
            Pe_m = 10.d0
        case(2)
            Pe_m = 100.d0
        case(3)
            Pe_m = 1000.d0
        case(4)
            Pe_m = 1500.d0
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
    
    
    !-------------------------------------------------
    write(*,*)'Choose a method to solve the problem.'
    write(*,*)'(1) FOU'
    write(*,*)'(2) SOCD'
    write(*,*)'(3) QUICK'
    write(*,*)'(4) TOU'
    write(*,*)'(5) SOU'
    read(*,*)method_select
    select case (method_select)
        case(1)
            G = 0
            beta = 0
            printout = 1
        case(2)
            G = 1
            beta = 0
            printout = 2
        case(3)
            G = 1
            beta = 1.d0/8.d0
            printout = 3
        case(4)
            G = 1
            beta = 1.d0/6.d0
            printout = 4
        case(5)
            G = 1
            beta = 1.d0/2.d0
            printout = 5
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
    
    !----------------------------------------------------
    write(*,*)'Please enter the relaxation coefficient.'
    read(*,*)w
    !----------------------------------------------------
    
    write(*,*)'Do you want to use Deferred Correction Method?'
    write(*,*)'(1) Yes'
    write(*,*)'(2) No'
    read(*,*)DCM_select
    select case (DCM_select)
        case(1)
            DCM = 1
        case(2)
            DCM = 0
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select    
    
end subroutine Initial