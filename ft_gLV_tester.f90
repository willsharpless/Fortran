program ft_tester
    implicit none

    integer :: n = 0,i
    real :: t = 0.0
    real, dimension(:), allocatable :: x, xdot, mu
    real, dimension(:,:), allocatable :: A
    external gLV

    print *, "Testing the gLV subroutine"
    print *, ""
    print *, "Size of community?"
    read *, n
    allocate(x(1:n))
    allocate(xdot(1:n))
    xdot = 0 ! initialized value irrelevant
    allocate(mu(1:n))
    allocate(A(1:n,1:n))

    print *, "Current state, x?"
    read *, x

    print *, "Growth Rates, mu?"
    read *, mu

    print *, "Interaction matrix, A?"
    read *, A
    print *, "A="
    do i=1,n; print '(20f6.2)', A(i,:); enddo
    
    call gLV (n, t, x, xdot, mu, A)
    print *, 'xdot = '
    print '(20f6.2)', xdot
    print *, "End of Program!"

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