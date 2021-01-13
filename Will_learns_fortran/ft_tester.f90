program ft_tester
    implicit none

    integer :: i_step, n_step, flag = 1, i
    integer, parameter :: neqn = 3
    real :: t, t_out, t_start, t_stop, abserr, relerr
    external gLV
    real, dimension(neqn):: x, xdot, mu
    real, dimension(neqn,neqn):: A
    
    ! Test Variables
    x = (/0.1, 0.1, 0.1/)
    print *, ""
    print '(a, 20f6.2)', 'x0 =', x
    xdot = (/0.1, 0.1, 0.1/) 
    mu = (/0.05, 0.1, 0.2/)
    A(1,1) = -0.02
    A(1,2) = 0
    A(1,3) = -0.05
    A(2,1) = -0.5
    A(2,2) = -0.2
    A(2,3) = -0.2
    A(3,1) = 0.3
    A(3,2) = -0.3
    A(3,3) = -0.3

    print *, "A="
    do i=1,neqn; print '(20f6.2)', A(i,:); enddo
    
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    t_start = 0.0E+00
    t_stop = 200.0E+00
    n_step = 200
    t_out = 0.0E+00
    t = t_out

    print *, ""
    print *, " Flag    Time       x1       x2       x3      xd1      xd2      xd3"
    call gLV ( neqn, t, x, xdot, mu, A )
    
    ! Iterative calls to RK45
    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start  &
            + real (          i_step - 1, kind = 4 ) * t_stop ) & 
            / real ( n_step,              kind = 4 )
    
        t_out = ( real ( n_step - i_step, kind = 4 ) * t_start  &
                + real (          i_step, kind = 4 ) * t_stop ) & 
                / real ( n_step,          kind = 4 )
        
        ! print '(a, 20f6.1)', "t = ", t
        ! print '(a, 20f6.1)', "t_out = ", t_out
        ! call r4_rkf45 ( r4_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )
        call r4_rkf45_glv ( gLV, neqn, x, xdot, t, t_out, relerr, abserr, flag, mu, A )

        ! print '(i5, f9.1, 20f9.3, 20f9.3)', flag, t, x, xdot
        print *, x, ";"

        if (flag == 7) then
            flag = 2
        end if
        
    end do

end program ft_tester

subroutine gLV ( n, t, x, xdot, mu, A )
    ! Calculation of the derivative of a gLV system with bounds x>=0
    ! willsharpless@berkeley.edu 

    ! n - size of community
    ! t - time
    ! x - current state array
    ! xdot - current state derivative array
    ! mu - array of innate growth rates
    ! A - matrix of interaction coefficients
    implicit none

    integer :: n
    real :: t, x(n), xdot(n), mu(n), A(n,n)

    ! Population is dead (apply bounds)
    where (x < 0.0) x = 0.0
    xdot = x * ( mu + matmul(A,x) )

    return
end subroutine gLV