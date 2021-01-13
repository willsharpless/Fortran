program Motif_Evolution_GA
    use M_scramble, only : scramble
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    implicit none
    external DGEEV
    external DGETRF
    external DGETRI
    integer :: res, i, ii, is, ig, ic, iara, n, unstable, singular, display, maxidx
    integer, parameter :: rows_nc = 1000, cols_tot = 132, cols_sc = 40
    integer, parameter :: N_tot = 11, N_sc = 5, N_nc = 6
    integer, parameter :: S = 40, G = 300, parents = 10, seed_n = 5, choice = 1
    integer, parameter :: Total_NC = rows_nc * seed_n
    double precision, dimension(Total_NC, cols_tot) :: nc_params ! stores all NC parameters, over written and used for scoring full comms
    double precision, dimension(seed_n, cols_sc) :: sc_params ! stores 5 SC parameters
    double precision, dimension(S, cols_sc) :: Lambda ! stores the evolved communities
    double precision:: Pi(S,4), SB(G,cols_sc+4+8) ! scoring mat for Lambda, and genertaional score board
    integer :: seed(5), sc_idx(cols_sc)
    character(100) :: filename
    double precision, allocatable :: r(:,:)
    integer :: r_int(S/2)
    ! integer :: test(4), test2(4,4)
    ! real, allocatable :: f(:,:)
    double precision, allocatable :: xfp(:,:), mu(:,:), A(:,:)
    real :: t_stop = 10000.0E+00
    integer :: n_step = 100
    integer :: skip(Total_NC,1) ! For inf tending communities (not sure why they exist in csv but checked w matlab)
    integer :: supported_naturally(Total_NC, 1)
    integer :: nan_penalty, nan_flag, scis
    double precision :: ra_nc(Total_NC, 1), ra_growth(Total_NC, 1), avg_ra_growth
    double precision :: standard ! score of the the seed Support Community
    real(real32) nan
    double precision :: holder_pi(4), holder_lambda(cols_sc)
    double precision :: avg_pts(4), avg_pop(4) 
    double precision :: phi(cols_sc), child(cols_sc)
    !r(S - 2*parents, cols_sc), r_int(S - 2*parents)

    print *, ""
    print *, ""
    print *, ""
    print *, ""
    print *, ""
    print *, "Beginning of Evolution Script"
    print *, "Will Sharpless - 12/22/20"
    print *, ""
    print *, ""
    print *, ""

    ! Initialize some stuff
    nan = ieee_value(nan, ieee_quiet_nan)
    Skip = 0
    Pi = 0
    supported_naturally = 0
    ra_nc = 0
    display = 0
    standard = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Open and Read the native community and motif parameter .csvs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Potential SC ID's to seed evolution:
    seed = (/ 678, 799, 1056, 1189, 1278 /)
    ! choice = 1 corresponds to SC 678

    ! other possislbe: 11, 331, 1402, 1763, 1943

    ! For each SC, the # of evolved params (Lambda cols)should be,
    ! 5 growth rates + 
    ! 5 self interactions (excludes target self) +
    ! 30 interactions for 6 mem comm 
    ! sc_params/Lambda should thus have 40 columns
    sc_idx = (/ 7,8,9,10,11,73,74,75,76,77,83,84,85,86,87,88,94,95,96,97,98,99,105 &
    ,106,107,108,109,110,116,117,118,119,120,121,127,128,129,130,131,132 /)

    ! Note I will use the 5 sets of interactions between the NC and SC from the 5 seed motifs to assay any evolved solutions (to limit dependence on those relationships)

    !!! Extracting Native Community and Support Community data from csv's
    print '(a)', "Reading the Native Communities, and their NC-SC interactions, from ..."
    do ii = 1,5
        write (filename, fmt='(a,i0,a)') 'Ten_Top_Scoring_SC/KillerComm_', seed(ii), '_stode_final_parameters.csv'
        print *, filename
        open(10, file=filename, access='sequential',form="formatted",iostat=res)

        do i=1,1000
            read(10,*) nc_params(i + rows_nc*(ii-1),:)
            if (i == 1) then
                sc_params(ii, :) = nc_params(i + rows_nc*(ii-1), sc_idx)
            end if
        end do
    end do

    print*,""

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Generating Lambda with SC 678
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Choice Evolution Method
    ! - Take each vector and vary within 5-30% (?) range of the value
    ! - Should it be all values or just a couple of them?
    ! - This will produce less random/diverse values, but more likely stable communities

    !! Alternatively
    ! - Randomly seed within the entire range of the motif
    ! - Will produce more diversity in the pool but also some wildly different, likely unstable values

    ! print '(a)', ""
    ! print '(a)', "Max and Min Value of the NC"
    ! print *, maxval(nc_params)
    ! print *, minval(nc_params)
    ! print '(a)', ""
    ! print '(a)', "Max and Min Value of the SC_1"
    ! print *, maxval(sc_params(1,:))
    ! print *, minval(sc_params(1,:))
    ! print '(a)', "Max and Min Value of the SC_2"
    ! print *, maxval(sc_params(2,:))
    ! print *, minval(sc_params(2,:))
    ! print '(a)', "Max and Min Value of the SC_3"
    ! print *, maxval(sc_params(3,:))
    ! print *, minval(sc_params(3,:))
    ! print '(a)', "Max and Min Value of the SC_4"
    ! print *, maxval(sc_params(4,:))
    ! print *, minval(sc_params(4,:))
    ! print '(a)', "Max and Min Value of the SC_5"
    ! print *, maxval(sc_params(5,:))
    ! print *, minval(sc_params(5,:))
    ! print*,""

    ! Seeding Lambda with Motif 678
    Lambda(1,:) = sc_params(choice,:) ! First row is SC itself
    ! next 20 rows are +-10% variation
    allocate (r(20, cols_sc))
    call random_number(r)
    r = 1 + ((r - 0.5)/5)
    Lambda(2:21,:) = r*transpose(spread(sc_params(1,:), 2, 20))
    deallocate (r)
    ! last 19 rows are with 50% variation given to 50% of the parameters
    Lambda(22:40,:) = transpose(spread(sc_params(1,:), 2, 19)) !base SC
    allocate (r(19, cols_sc))
    call random_number(r)
    r = (3*r - 1.5) ! +-150% range, so it can switch signs
    do i=22,40
        !r_int(1,:) = scramble(40)
        r_int(:) = scramble(S)
        !Lambda(i,r_int(1,:)) = r(i-21,r_int(1,:))*Lambda(i,r_int(1,:))
        Lambda(i,r_int(:)) = r(i-21,r_int(:))*Lambda(i,r_int(:))
    enddo
    ! Ensure positive growth rates (150% mightve fucked it up)
    Lambda(:,1:5) = abs(Lambda(:,1:5))
    deallocate (r)
    allocate (r(S - 2*parents, cols_sc))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Score baseline survival without Support Community
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n = N_nc
    allocate (mu(n,1))
    allocate (A(n,n))
    allocate (xfp(n,1))
    do ic=1,Total_NC

        mu(:,1) = nc_params(ic,1:n)
        do i =1,n
            A(i,:) = nc_params(ic, (i*N_tot + 1) : (i*N_tot + n))
            ! print*,A(i,:)
        enddo
        xfp = 0.001 ! initial conc
        nan_flag = 0
        call gLV_solution (n, xfp, mu, A, t_stop, n_step, display, nan_flag)
        xfp = dble(xfp)
        ! print*,xfp(6,1)
        ! nan_flag will catch some but not all nans
        if (nan_flag == 1) then; xfp = nan; endif
        if (isnan(xfp(6,1))) then; skip(ic,1)=1; cycle; endif

        ! if it represents less than 1% before then not supported naturally
        if (xfp(6,1)/sum(xfp) > 0.01) then
            supported_naturally(ic,1) = 1 
            ra_nc(ic,1) = xfp(6,1)/sum(xfp)
        endif
        
        ! Store xfp(6,1) total abundance to guage improvement with SC later?
        ! is xfp(6,1) close to/ == 0? store in binary?

    enddo
    deallocate (mu)
    deallocate (A)
    deallocate (xfp)

    print*, "Number of infinity tending NC:"
    print*, (sum(skip)/5)
    print*, "Number of NC which support MOI naturally:"
    print*, (sum(supported_naturally)/5)
    print*,""
    ! print*,ra_nc(1:1000,1)

    ! print*, "Skip"
    ! print*, skip(1:25,1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lambda Generations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n = N_tot
    allocate (mu(n,1))
    allocate (A(n,n))
    allocate (xfp(n,1))


    do ig = 1,G

        ! Print generation % and store top scoring Lambda
        print*, "Lambda Generation ", ig, "/", G
        print*, ""
        call timestamp

        if (.not. ig == 1 ) then
            print*, "Best Score thus far:"
            print*, SB(ig-1,cols_sc+4)
        endif

        do is = 1,S !!! Iterate through the vectors of Lambda

            ! Skip already scored parents
            if ((ig > 1) .and. (is <= parents)) then; cycle; endif
            
            ! Print generation % and store top scoring Lambda
            ! print*, "Lambda ", is, ":"
            ! call timestamp
            ! print*, ""

            ra_growth = 0
            nan_penalty = 0
            avg_ra_growth = 0
            iara = 0
            scis = 0

            do ic = 1,Total_NC !!! Iterate through the native communities

                ! Skip the inf tending native communities
                if (skip(ic,1) == 1) then;  cycle; endif

                ! Stash Lambda in nc_params(ic, sc_idx) for simple porting
                nc_params(ic, sc_idx) = Lambda(is,:)
                ! Distribute community to temporary work variables
                mu(:,1) = nc_params(ic,1:n)
                do i =1,n
                    A(i,:) = nc_params(ic, (i*n + 1) : (i*n + n))
                enddo
                ! numerically derive xfp, after 10000 hours
                xfp = 0.001
                nan_flag = 0
                ! if (is == 2) then; display = 1; endif
                call gLV_solution (n, xfp, mu, A, t_stop, n_step, display, nan_flag)

                if (nan_flag == 1) then; xfp = nan; endif

                if (.not. isnan(xfp(6,1))) then
                    ra_growth(ic,1) = xfp(6,1)/sum(xfp) - ra_nc(ic,1)

                    ! nan_flag catches some, but not all nans
                    if (isnan(ra_growth(ic,1))) then
                        nan_penalty = nan_penalty + 1
                        ra_growth(ic,1) = 0
                    else
                        ! Summing for the Average Relative Abundance Change
                        avg_ra_growth = avg_ra_growth + ra_growth(ic,1)
                        iara = iara + 1

                        ! Summing the # of c where xfp rel abu jumps from <1% to > 10% (survival)
                        if ((supported_naturally(ic,1) == 0) .and. (xfp(6,1)/sum(xfp) > 0.1)) then
                            scis = scis + 1
                        endif
                    endif
                else
                    nan_penalty = nan_penalty + 1
                endif
                
                ! Reset dummy variables for stability analysis
                ! unstable = 0
                ! singular = 0
                ! ! gLV Stability computes the nontrivial fixed point and its stability
                ! call gLV_stability (n, xfp, mu, A, unstable, singular)
                ! ! Filter out unstable, singular and negative stabilizing communities
                ! if ((unstable == 1) .or. (singular == 1) .or. (minval(xfp) < 0)) then
                !     ! make score for Lambda(s,:) huge
                !     ! skip loop step
                ! endif
                ! If semipositive and stable xfp, then we score based on abundances

            enddo

            avg_ra_growth = avg_ra_growth/iara
            ! print*,"Average Relative Abundance Growth"
            ! print*,avg_ra_growth
            ! print*,"Support Community Induced Survival"
            ! print*,scis
            ! print*,"Nan Penalty (excludes nc nans)"
            ! print*,nan_penalty
            ! print*, ""

            Pi(is, 1) = avg_ra_growth
            Pi(is, 2) = scis
            Pi(is, 3) = nan_penalty
            Pi(is, 4) = 1000*avg_ra_growth + scis + (-2)*nan_penalty - standard ! Overall score of Lambda

            ! Standardize scores to original motif score
            if ((ig == 1) .and. (is == 1)) then
                standard = Pi(is, 4)
                Pi(is, 4) = Pi(is, 4) - standard ! standardizing original sc as well
            endif
            
        enddo

        ! call vector_breeding(Lambda, Pi, SB, sc_params, seed_n, choice, parents, ig, S, G, cols_sc )
        ! I am concerned there is a copy-in/copy-out bug in gfortran that originates from this subroutine
        ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58310
        ! https://groups.google.com/g/comp.lang.fortran/c/VuFvOsLs6hE
        ! at least others on the internet were seeing the similar error
        ! I need to get valgrind or gdb to get in depth debuggin but for now I am just going to move the subroutine mechanics here

        !!!  Sort Pi and Lambda
        do i=1,S
            maxidx = maxloc(Pi(i:S,4), dim=1)+i-1

            holder_pi(:) = Pi(i,:)
            holder_lambda(:) = Lambda(i,:)

            Pi(i,:) = Pi(maxidx,:)
            Lambda(i,:) = Lambda(maxidx,:)

            Pi(maxidx,:) = holder_pi(:)
            Lambda(maxidx,:) = holder_lambda(:)

            !! Running averages after sorting
            if (i <= parents) then
                avg_pts(:) = avg_pts(:) + Pi(i,:)
            endif
            avg_pop(:) = avg_pop(:) + Pi(i,:)
        enddo

        avg_pts = avg_pts/parents
        avg_pop = avg_pop/S

        !!! Store best lambda, best score, avg parents, avg pop in SB
        SB(ig, 1:cols_sc) = Lambda(1,:)
        SB(ig, cols_sc+1:cols_sc+4) = Pi(1,:)
        SB(ig, cols_sc+5:cols_sc+8) = avg_pts(:)
        SB(ig, cols_sc+9:cols_sc+12) = avg_pop(:)

        ! print*, "Scoreboard------------------" 
        ! print*, "Generation: "
        ! print*, ig
        ! print*, ""
        ! print*, SB(ig,:)
        ! print*,""

        !!! Iterate through parents, breeding new children vectors
        do i=1,parents,2
            ! Child 1
            call random_number(phi)
            phi = 2*phi + 0.5 ! 150% range between (0.5x to 1.5y, no sign change)
            child(:) = Lambda(i,:)*phi + Lambda(i+1,:)*(1-phi)
            Lambda(parents+i,:) = child(:)
            Pi(parents+i,:) = 0

            ! Child 2
            call random_number(phi)
            phi = 2*phi + 0.5 ! 150% range between (0.5x to 1.5y, no sign change)
            child(:) = Lambda(i,:)*phi + Lambda(i+1,:)*(1-phi)
            Lambda(parents+i+1,:) = child(:)
            Pi(parents+i+1,:) = 0

            ! Might try larger ranges, sign flipping, or more children/parents
        enddo

        ! print*, "Did Breeding"

        !!!  Generate fresh vectors at bottom
        ! New pool generated with 50% variation of seed SC given to 50% of the parameters
        Lambda(2*parents+1:S,:) = transpose(spread(sc_params(choice,:), 2, 20))
        ! If stuck, might be worth populating with best vector in lambda not seed SC
        call random_number(r)
        r = (3*r - 1.5) ! +-150% range, so it can switch signs
        do i=2*parents+1,S
            r_int(:) = scramble(S)
            Lambda(i,r_int(:)) = r(i-2*parents, r_int(:))*Lambda(i,r_int(:))
        enddo
        Lambda(:,1:5) = abs(Lambda(:,1:5)) ! Ensure positive growth rates (150% mightve fucked it up)
        Pi(2*parents+1:S, :) = 0 ! Reset score matrix

        print*,""
    enddo

    ! Plot the Best, Average-Parents, and Average-Pop scores over generations

    ! Simulate and Plot the Best Scoring Motif

    ! Vary the Scoring Metric
    ! Vary the Evolution breeding, particularly the randomization
    ! Try with the other 5 motifs

    print*, ""
    print*, "---------------------------------------------------------------------------------------------------------------"
    print*, "-----------------------------------------------   Score Board   -----------------------------------------------"
    print*, "---------------------------------------------------------------------------------------------------------------"

    do i=1,G
        print*, SB(i,:)
    enddo

    print*, ""
    print*, ""
    print*, ""
    print*, "He terminado"
    call timestamp
    print *, ""
    print *, ""
    print *, ""
    print *, ""
    print *, ""
end program Motif_Evolution_GA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    real :: x(n), xdot(n), mu(n), A(n,n), t

    ! Population is dead (apply bounds)
    where (x < 0.0) x = 0.0
    xdot = x * ( mu + matmul(A,x) )

    return
end subroutine gLV

subroutine gLV_solution (n, xinput, muinput, Ainput, t_stop, n_step, display, nan_flag)
    implicit none
    integer :: i, n_step, flag = 1, n, display, nan_flag
    real :: t, t_out, t_start, t_stop, abserr, relerr
    external gLV
    double precision :: xinput(n,1), muinput(n,1), Ainput(n,n)
    real :: x(n), xdot(n), mu(n), A(n,n)
    
    x = real(reshape(xinput,shape(x)))
    mu = real(reshape(muinput,shape(mu)))
    A = real(Ainput)
    xdot = 0
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )
    t_start = 0.0E+00
    t_out = 0.0E+00
    t = t_out

    ! print*,mu
    ! do i=1,n; print*,A(i,:); enddo

    ! Initialize with call to gLV
    if (display == 1) then
        print *, ""
        print *, " Flag    Time       x1       x2       x3      x4      x5      x6"
    endif

    call gLV ( n, t, x, xdot, mu, A )

    if (display == 1) then; print '(i5, f9.1, 20f9.3, 20f9.3)', flag, t, x; endif

    ! Iterative calls to RK45
    do i = 1, n_step
        t = ( real ( n_step - i + 1, kind = 4 ) * t_start  &
            + real (          i - 1, kind = 4 ) * t_stop ) & 
            / real ( n_step,              kind = 4 )
        t_out = ( real ( n_step - i, kind = 4 ) * t_start  &
                + real (          i, kind = 4 ) * t_stop ) & 
                / real ( n_step,          kind = 4 )
        
        ! print '(a, 20f6.1)', "t = ", t
        ! print '(a, 20f6.1)', "t_out = ", t_out
        ! print '(a,i5)', "Flag:", flag
        ! print '(a,i5)', "Nan_Flag:", nan_flag 
        call r4_rkf45_glv ( gLV, n, x, xdot, t, t_out, relerr, abserr, flag, mu, A)

        ! The GLV systems that start at 0.001 and exceed 1 tend to infinity
        if (sum(x) >= 1) then
            nan_flag = 1
            flag = 1
            ! print*, "NAN FLAG" 
            return
        endif

        if (display == 1) then
            print '(i5, f9.1, 20f9.3, 20f9.3)', flag, t, x
        endif

        if (flag == 7) then
            abserr = 10*abserr
            relerr = 10*relerr
            flag = 2
        elseif (flag == 6) then
            abserr = 10*abserr
            relerr = 10*relerr
            flag = 2
        end if

    end do

    xinput = dble(reshape(x,shape(xinput)))

end subroutine gLV_solution

subroutine gLV_stability (n, xfp, mu, A, unstable, singular)
    ! Calculation of the fixed point of a glv system and its stability
    ! willsharpless@berkeley.edu 

    ! n - size of community
    ! x - state array at the fixed point
    ! mu - innate growth rate array
    ! A - matrix of interaction coefficients
    ! Jc - Jacobian array

    ! The Jacobian of the GLV system is equivalent to:
    ! Diag(r + Ax) + (x' hadamard A)
    implicit none
    ! external DGEEV
    ! external DGETRF
    ! external DGETRI
    integer :: n,i,j,ok,info,unstable,singular
    integer :: ipiv(n)
    double precision :: xfp(n,1), mu(n,1), A(n,n), Ainv(n,n), Jc(n,n), tmp(1), eigR(n), eigI(n)
    double precision :: dummy(n,n), work(n)
    Jc = 0

    ! Compute nontrivial fixed point
    Ainv = A
    call DGETRF(n, n, Ainv, n, ipiv, info) ! LU Factorization of A
    if (info /= 0) then
        print*, 'Matrix is numerically singular!'
        singular = 1
        return
    end if
    call DGETRI(n, Ainv, n, ipiv, work, n, info) ! Inversion of LU Fact of A

    ! print*, "Ainv:"
    ! do i =1,n
    !     print*,Ainv(i,:)
    ! enddo
    ! print*,""

    xfp = matmul(Ainv,-mu)
    print*, "xfp"
    print *, xfp
    stop
    
    ! Compute Jacobian at nontrivial fixed point
    do i=1,n
        do j=1,n
            if (i==j) then 
                tmp = matmul(A(i,:),xfp)
                Jc(i,j) = mu(i,1) + tmp(1)
            endif
            Jc(i,j) = Jc(i,j) + A(i,j)*xfp(i,1)
        enddo
    enddo

    ! Get Eigenvalues of (real) Jacobian
    call DGEEV('N', 'N', n, Jc, n, eigR, eigI, dummy, n, dummy, n, work, 3*n, work, ok)

    do i=1,n
        if (eigR(i) >= 0) then
            unstable = 1
        endif
    enddo

    return
end subroutine gLV_stability

    ! do i=1,10
    !     r_int(i,:) = (/(ii, ii=1,20, 1)/)
    !     ! r_int(i,:) = (/1,2,3,4,5/)
    ! enddo

subroutine vector_breeding (Lambda, Pi, SB, sc_params, seed_n, choice, parents, ig, S, G, cols_sc )
    use M_scramble, only : scramble
    implicit none
    integer :: seed_n, choice, parents, ig, S, G, cols_sc
    integer :: i, maxidx
    double precision :: Lambda(S, cols_sc), Pi(S,4), SB(G,cols_sc+4+8), sc_params(seed_n, cols_sc)
    double precision :: holder_pi(4), holder_lambda(cols_sc)
    double precision :: avg_pts(4), avg_pop(4) 
    double precision :: phi(cols_sc), child(cols_sc), r(S - 2*parents, cols_sc), r_int(S - 2*parents)


    ! print*, "Scores (arag, scis, np):"
    ! do i=1,S; print*,Pi(i,:);enddo
    ! print*,""

    ! print*, "First Lambda pre sorting"
    ! print*, Lambda(1,:)
    ! print*, ""

    !!!  Sort Pi and Lambda
    do i=1,S
        maxidx = maxloc(Pi(i:S,4), dim=1)+i-1

        holder_pi(:) = Pi(i,:)
        holder_lambda(:) = Lambda(i,:)

        Pi(i,:) = Pi(maxidx,:)
        Lambda(i,:) = Lambda(maxidx,:)

        Pi(maxidx,:) = holder_pi(:)
        Lambda(maxidx,:) = holder_lambda(:)

        !! Running averages after sorting
        if (i <= parents) then
            avg_pts(:) = avg_pts(:) + Pi(i,:)
        endif
        avg_pop(:) = avg_pop(:) + Pi(i,:)
    enddo

    ! print*, "Scores (arag, scis, np):"
    ! do i=1,S; print*,Pi(i,:);enddo

    avg_pts = avg_pts/parents
    avg_pop = avg_pop/S

    ! print*, "First Lambda post sorting (Best):"
    ! print*, Lambda(1,:)
    ! print*, ""

    !!! Store best lambda, best score, avg parents, avg pop in SB
    SB(ig, 1:cols_sc) = Lambda(1,:)
    SB(ig, cols_sc+1:cols_sc+4) = Pi(1,:)
    SB(ig, cols_sc+5:cols_sc+8) = avg_pts(:)
    SB(ig, cols_sc+9:cols_sc+12) = avg_pop(:)

    print*, "Scoreboard------------------" 
    print*, "Generation: "
    print*, ig
    print*, ""
    print*, SB(ig,:)
    print*,""

    !!! Iterate through parents, breeding new children vectors
    do i=1,parents,2
        ! Child 1
        call random_number(phi)
        phi = 2*phi + 0.5 ! 150% range between (0.5x to 1.5y, no sign change)
        child(:) = Lambda(i,:)*phi + Lambda(i+1,:)*(1-phi)
        Lambda(parents+i,:) = child(:)
        Pi(parents+i,:) = 0

        ! Child 2
        call random_number(phi)
        phi = 2*phi + 0.5 ! 150% range between (0.5x to 1.5y, no sign change)
        child(:) = Lambda(i,:)*phi + Lambda(i+1,:)*(1-phi)
        Lambda(parents+i+1,:) = child(:)
        Pi(parents+i+1,:) = 0

        ! Might try larger ranges, sign flipping, or more children/parents
    enddo

    ! print*, "Did Breeding"

    !!!  Generate fresh vectors at bottom
    ! New pool generated with 50% variation of seed SC given to 50% of the parameters
    Lambda(2*parents+1:S,:) = transpose(spread(sc_params(choice,:), 2, 20))
    ! If stuck, might be worth populating with best vector in lambda not seed SC
    call random_number(r)
    r = (3*r - 1.5) ! +-150% range, so it can switch signs
    do i=2*parents+1,S
        r_int(:) = scramble(S)
        Lambda(i,r_int(:)) = r(i-2*parents, r_int(:))*Lambda(i,r_int(:))
    enddo
    Lambda(:,1:5) = abs(Lambda(:,1:5)) ! Ensure positive growth rates (150% mightve fucked it up)
    Pi(2*parents+1:S, :) = 0 ! Reset score matrix

    ! print*, "Regenerated Lambda and Reset Scores"
    ! print*, ""

    ! print*, " Regen Scores (arag, scis, np):"
    ! do i=1,S; print*,Pi(i,:);enddo
    ! print*, ""
    return
end subroutine vector_breeding