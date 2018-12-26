!----------------------------------------------------------------------------------
subroutine deposition_rate(N,divide_1, divide_2,del_x, del_y, CVD,J_x, J_avg_x)
    implicit NONE
    
    integer :: i,N
    real*8 :: divide_1, divide_2, del_y,del_x,sum
    real*8, dimension(N) :: J_x, J_avg_x
    real*8, dimension(N,N) :: CVD
    
    !calculate J
    J_x = 0.d0
    J_avg_x = 0.d0
    
    do i = divide_1,divide_2,1
        J_x(i) = -(CVD(i,N) - CVD(i,N-1))/del_y
    end do
    
    !calculate J_bar
    sum = 0.d0
    J_avg_x(divide_1) = J_x(divide_1)
    do i = divide_1+1,divide_2,1
        J_avg_x(i) = (sum + (J_x(i-1)+J_x(i))*del_x/2.d0)/((i-divide_1)*del_x)
        sum = sum + (J_x(i-1)+J_x(i))*del_x/2.d0
    end do
    
end subroutine deposition_rate