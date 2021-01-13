program Motif_Evolution_GA_validation
    ! A genetic algorithm for evolving support communities that promote the survival of a member of interest
    ! This is a special validation script
    ! Will Sharpless - 1/12/21
    ! willsharpless@berkeley.edu

    use M_scramble, only : scramble
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    implicit none
    external DGEEV
    external DGETRF
    external DGETRI
    integer :: res, i, ii, is, ig, ic, iara, n, display, maxidx, nol,lam_count, file_count !unstable, singular
    integer, parameter :: rows_nc = 1000, cols_tot = 132, cols_sc = 40
    integer, parameter :: N_tot = 11, N_sc = 5, N_nc = 6
    integer, parameter :: S = 40, G = 300, parents = 10, seed_n = 5, choice = 1
    integer, parameter :: Total_NC = rows_nc * seed_n
    integer, parameter :: N_tb = 17
    double precision, dimension(Total_NC, cols_tot) :: nc_params ! stores all NC parameters, over written and used for scoring full comms
    double precision, dimension(seed_n, cols_sc) :: sc_params ! stores 5 SC parameters
    double precision, dimension(N_tb, cols_sc) :: Lambda ! stores the evolved communities
    double precision:: Pi(S,4), SB(G,cols_sc+4+8) ! scoring mat for Lambda, and genertaional score board
    integer :: seed(5), sc_idx(cols_sc)
    character(100) :: filename
    character(30), dimension(N_tb) :: lambda_names
    double precision, allocatable :: r(:,:)
    integer :: r_int(S/2)
    double precision, allocatable :: xfp(:,:), mu(:,:), A(:,:)
    double precision :: NC_xfp(rows_nc, N_nc), seedSC_xfp(Total_nc, N_tot), evolSC_xfp(Total_nc, N_tot)
    real :: t_stop = 10000.0E+00
    integer :: n_step = 100
    integer :: skip(Total_NC,1) ! For inf tending communities (not sure why they exist in csv but checked w matlab)
    integer :: supported_naturally(rows_nc, 1)
    integer :: nan_penalty, nan_flag, scis, scis_aorig
    double precision :: ra_nc(Total_NC, 1), ra_growth(Total_NC, 1), avg_ra_growth
    double precision :: standard ! score of the the seed Support Community
    real(real32) nan
    double precision :: holder_pi(4), holder_lambda(cols_sc)
    double precision :: avg_pts(4), avg_pop(4) 
    double precision :: phi(cols_sc), child(cols_sc)
    double precision :: load(cols_sc)

    ! This script will read several top scoring lambda stored in "thoroughbred_011221.csv" and compute the xfp of them, the original seed support communities and native communities
    ! and write to files called "thoroughbred_performance_011221.csv"

    ! v2 edits:
! 	- Remove avg_ra_growth from score, still tracking though for interest
!   - Bounded all generated values by the NC max's and min's
! 	- Create a binary cut off >< 1E-8 to discern "survival" (instead of <1% to >10% relative abundance jump)
! 	- Lift the death cut off in the simulations to 1E-8 (instead of <=0)
!   - Lower the community total abundance to 0.1 (instead of 1)

    ! v3 edits:
!   - Instead of generating new lambda by % variations of the seed, populating with totally random values within the max and min of the corresponding NC values

    print *, ""
    print *, ""
    print *, ""
    print *, ""
    print *, ""
    print *, "Beginning of Evolution Validation Script"
    print *, "Will Sharpless"
    call timestamp
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
    NC_xfp = 0
    file_count = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Open and Read the native community and motif parameter .csvs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Potential SC ID's to seed evolution:
    seed = (/ 678, 799, 1056, 1189, 1278 /)
    ! choice = 1 corresponds to SC 678

    ! other possislbe: 11, 331, 1402, 1763, 1943

    sc_idx = (/ 7,8,9,10,11,73,74,75,76,77,83,84,85,86,87,88,94,95,96,97,98,99,105 &
    ,106,107,108,109,110,116,117,118,119,120,121,127,128,129,130,131,132 /)

    ! Note I will use the 5 sets of interactions between the NC and SC from the 5 seed motifs to assay any evolved solutions (to limit dependence on those relationships)

    !!! Extracting Native Community and Support Community data from csv's
    print '(a)', "Reading the Native Communities, and their NC-SC interactions, from ..."
    do ii = 1,5
        write (filename, fmt='(a,i0,a)') 'Ten_Top_Scoring_SC/KillerComm_', seed(ii), '_stode_final_parameters.csv'
        print *, filename
        open(file_count, file=filename, access='sequential',form="formatted",iostat=res)

        do i=1,1000
            read(file_count,*) nc_params(i + rows_nc*(ii-1),:)
            if (i == 1) then
                sc_params(ii, :) = nc_params(i + rows_nc*(ii-1), sc_idx)
            end if
        end do
    end do
    print*,""

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Reading Lambda from Thoroughbred
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Lambda = 0
    file_count = file_count + 1
    print*,""
    print '(a)', "Reading the Throughbred Lambdas, and their NC-SC interactions"
    write (filename, fmt='(a)') 'Evolution_runs/thoroughbred_011221.csv'
    open(file_count, file=filename, access='sequential',form="formatted",iostat=res)
    read(file_count,*)
    read(file_count,*)
    do ii = 1,15 !!!NEEDS TO BE CHANGED WHEN SCRIPT FINISHED N_tb
        read(file_count,*) lambda_names(ii), Lambda(ii,:)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Native Comm xfp Derivation & Naturally Supported Score
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n = N_nc
    allocate (mu(n,1))
    allocate (A(n,n))
    allocate (xfp(n,1))
    do ic=1,rows_nc

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

        NC_xfp(ic,:) = xfp(:,1) !store xfp for validation scoring
        
        if (isnan(xfp(6,1))) then
            skip(ic,1)=1 
            skip(ic+1000,1)=1
            skip(ic+2000,1)=1
            skip(ic+3000,1)=1
            skip(ic+4000,1)=1
            cycle 
        endif

        ! if MOI abundance > 1.0E-008 then its "surviving"
        if (xfp(6,1) > 1.0E-008) then ! old req: (xfp(6,1)/sum(xfp) > 0.01)
            supported_naturally(ic,1) = 1 
        endif

    enddo
    deallocate (mu)
    deallocate (A)
    deallocate (xfp)

    print*, "Number of NaN-tending NC (10k hrs):"
    print*, (sum(skip)/5), " / 1000"
    print*, "Number of NC where xfp(MOI) > 1E-8 naturally:"
    print*, (sum(supported_naturally)), " / 1000"

    ! Write Native Comm xfp to csv
    file_count = file_count + 1
    write (filename, fmt='(a)') 'Evolution_runs/xfp_csvs/xfp_NC.csv'
    open(file_count, file=filename, access='sequential',form="formatted",iostat=res)
    do ic=1,rows_nc
        write(file_count,'(*(G0.6,:,","))') NC_xfp(ic,:)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Support Comm xfp Derivation & Validation Scoring
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n = N_tot
    allocate (mu(n,1))
    allocate (A(n,n))
    allocate (xfp(n,1))
    lam_count = 0

    do ig = 1,seed_n !!! Iterate through seed support comms

        if (ig == 1) then
            nol = 5 + 1 !runs + 1 for seed
        else
            nol = 3 + 1 !runs + 1 for seed
        endif

        do is = 1, nol !!! Iterate through the vectors of Lambda

            if (is == 1) then ! seed round
                load = sc_params(ig,:)
                seedSC_xfp = 0
                print*,""
                print*, "Original Support Community:", seed(ig)
            else ! lambda round
                lam_count = lam_count + 1
                load = Lambda(lam_count,:)
                evolSC_xfp = 0
                print*,""
                print*, "Evolved Support Community:    ", trim(lambda_names(lam_count))
            endif

            ra_growth = 0
            avg_ra_growth = 0
            iara = 0
            nan_penalty = 0
            scis = 0
            scis_aorig = 0

            do ic = 1,Total_NC !!! Iterate through the native communities

                ! Skip the inf tending native communities
                if (skip(ic,1) == 1) then;  cycle; endif

                ! Stash Lambda in nc_params(ic, sc_idx) for simple porting
                ! nc_params(ic, sc_idx) = Lambda(is,:)
                nc_params(ic, sc_idx) = load(:)
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

                ! Storing the xfp for export
                if (is == 1) then ! seed round
                    seedSC_xfp(ic,:) = xfp(:,1)
                else ! lambda round
                    evolSC_xfp(ic,:) = xfp(:,1)
                endif

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

                        ! Summing the # of c where xfp rel abu jumps from <1.0E-008 to >1.0E-008 (=:survival)
                        if (xfp(6,1) > 1.0E-008) then ! old requirement (xfp(6,1)/sum(xfp) > 0.1)) >10%
                            scis = scis + 1

                            ! Summing the # for the specific alpha-nc-sc derived from the seed community
                            if ((ic >= (rows_nc*(ig-1) + 1)) .and. (ic < (rows_nc*(ig) + 1))) then
                                scis_aorig = scis_aorig + 1
                            endif
                        endif
                    endif
                else
                    nan_penalty = nan_penalty + 1
                endif
            enddo

            avg_ra_growth = avg_ra_growth/iara
            ! print*,"Average Relative Abundance Growth"
            ! print*,avg_ra_growth
            ! print*,"Support Community Induced Survival"
            ! print*,scis
            ! print*,"Nan Penalty (excludes nc nans)"
            ! print*,nan_penalty
            ! print*, ""

            print*, "- # of comm where xfp(MOI) > 1E-8, with original SC {a_nc-sc}: ", scis_aorig, "/1000"
            print*, "- # of comm where xfp(MOI) > 1E-8, with 5 different {a_nc-sc}: ", scis, "/5000"
            print*, "(NaNs: ", (nan_penalty + (sum(skip)/5))," /5000 )"

            !!! Export to csv's
            file_count = file_count + 1
            if (is == 1) then ! seed round
                write (filename, fmt='(a,i0,a)') 'Evolution_runs/xfp_csvs/xfp_SC_', seed(ig), '.csv'
                open(file_count, file=filename, access='sequential',form="formatted",iostat=res)
                do ic=1,Total_NC
                    write(file_count,'(*(G0.6,:,","))') seedSC_xfp(ic,:)
                enddo
            else ! lambda round
                write (filename, fmt='(a,a,a)') 'Evolution_runs/xfp_csvs/xfp_eSC_', trim(lambda_names(lam_count)), '.csv'
                open(file_count, file=filename, access='sequential',form="formatted",iostat=res)
                do ic=1,Total_NC
                    write(file_count,'(*(G0.6,:,","))') evolSC_xfp(ic,:)
                enddo
            endif
        enddo
    enddo


    ! Plot the Best, Average-Parents, and Average-Pop scores over generations

    ! Simulate and Plot the Best Scoring Motif

    ! Vary the Scoring Metric
    ! Vary the Evolution breeding, particularly the randomization
    ! Try with the other 5 motifs

    ! print*, ""
    ! print*, "---------------------------------------------------------------------------------------------------------------"
    ! print*, "-----------------------------------------------   Score Board   -----------------------------------------------"
    ! print*, "---------------------------------------------------------------------------------------------------------------"

    ! do i=1,G
    !     print*, SB(i,:)
    ! enddo

    ! print*, ""
    ! print*, ""
    ! print*, "Best evolved Support Community:"
    ! print*, SB(G,1:cols_sc)
    ! print*, ""
    ! print*, "Which induces MOI survival in (out of 5000),"
    ! print*, SB(G,1:cols_sc)
    ! print*, "excluding the number of naturally supporting Native Communities (out of 5000),"
    ! print*, sum(supported_naturally)
    print *, ""
    print *, ""
    print*, "He terminado"
    call timestamp
    print *, ""
    print *, ""
    print *, ""
end program Motif_Evolution_GA_validation

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
    where (x < 1.0E-008) x = 0.0
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
        ! Reduced to 0.1 in v2
        if (sum(x) >= 0.1) then
            nan_flag = 1
            flag = 1
            ! print*, "NAN FLAG" 
            return
        endif

        if (display == 1) then
            print '(i5, f9.1, 20f9.3, 20f9.3)', flag, t, x
        endif

        if ((flag == 7) .or. (flag == 6)) then
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

! subroutine vector_breeding (Lambda, Pi, SB, sc_params, seed_n, choice, parents, ig, S, G, cols_sc )
!     use M_scramble, only : scramble
!     implicit none
!     integer :: seed_n, choice, parents, ig, S, G, cols_sc
!     integer :: i, maxidx
!     double precision :: Lambda(S, cols_sc), Pi(S,4), SB(G,cols_sc+4+8), sc_params(seed_n, cols_sc)
!     double precision :: holder_pi(4), holder_lambda(cols_sc)
!     double precision :: avg_pts(4), avg_pop(4) 
!     double precision :: phi(cols_sc), child(cols_sc), r(S - 2*parents, cols_sc), r_int(S - 2*parents)


!     ! print*, "Scores (arag, scis, np):"
!     ! do i=1,S; print*,Pi(i,:);enddo
!     ! print*,""

!     ! print*, "First Lambda pre sorting"
!     ! print*, Lambda(1,:)
!     ! print*, ""

!     !!!  Sort Pi and Lambda
!     do i=1,S
!         maxidx = maxloc(Pi(i:S,4), dim=1)+i-1

!         holder_pi(:) = Pi(i,:)
!         holder_lambda(:) = Lambda(i,:)

!         Pi(i,:) = Pi(maxidx,:)
!         Lambda(i,:) = Lambda(maxidx,:)

!         Pi(maxidx,:) = holder_pi(:)
!         Lambda(maxidx,:) = holder_lambda(:)

!         !! Running averages after sorting
!         if (i <= parents) then
!             avg_pts(:) = avg_pts(:) + Pi(i,:)
!         endif
!         avg_pop(:) = avg_pop(:) + Pi(i,:)
!     enddo

!     ! print*, "Scores (arag, scis, np):"
!     ! do i=1,S; print*,Pi(i,:);enddo

!     avg_pts = avg_pts/parents
!     avg_pop = avg_pop/S

!     ! print*, "First Lambda post sorting (Best):"
!     ! print*, Lambda(1,:)
!     ! print*, ""

!     !!! Store best lambda, best score, avg parents, avg pop in SB
!     SB(ig, 1:cols_sc) = Lambda(1,:)
!     SB(ig, cols_sc+1:cols_sc+4) = Pi(1,:)
!     SB(ig, cols_sc+5:cols_sc+8) = avg_pts(:)
!     SB(ig, cols_sc+9:cols_sc+12) = avg_pop(:)

!     print*, "Scoreboard------------------" 
!     print*, "Generation: "
!     print*, ig
!     print*, ""
!     print*, SB(ig,:)
!     print*,""

!     !!! Iterate through parents, breeding new children vectors
!     do i=1,parents,2
!         ! Child 1
!         call random_number(phi)
!         phi = 2*phi + 0.5 ! 150% range between (0.5x to 1.5y, no sign change)
!         child(:) = Lambda(i,:)*phi + Lambda(i+1,:)*(1-phi)
!         Lambda(parents+i,:) = child(:)
!         Pi(parents+i,:) = 0

!         ! Child 2
!         call random_number(phi)
!         phi = 2*phi + 0.5 ! 150% range between (0.5x to 1.5y, no sign change)
!         child(:) = Lambda(i,:)*phi + Lambda(i+1,:)*(1-phi)
!         Lambda(parents+i+1,:) = child(:)
!         Pi(parents+i+1,:) = 0

!         ! Might try larger ranges, sign flipping, or more children/parents
!     enddo

!     !!!  Generate fresh vectors at bottom
!     ! New pool generated with 50% variation of seed SC given to 50% of the parameters
!     Lambda(2*parents+1:S,:) = transpose(spread(sc_params(choice,:), 2, 20))
!     ! If stuck, might be worth populating with best vector in lambda not seed SC
!     call random_number(r)
!     r = (3*r - 1.5) ! +-150% range, so it can switch signs
!     do i=2*parents+1,S
!         r_int(:) = scramble(S)
!         Lambda(i,r_int(:)) = r(i-2*parents, r_int(:))*Lambda(i,r_int(:))
!     enddo
!     Lambda(:,1:5) = abs(Lambda(:,1:5)) ! Ensure positive growth rates (150% mightve fucked it up)
    
!     Pi(2*parents+1:S, :) = 0 ! Reset score matrix

!     ! print*, "Regenerated Lambda and Reset Scores"
!     ! print*, ""

!     ! print*, " Regen Scores (arag, scis, np):"
!     ! do i=1,S; print*,Pi(i,:);enddo
!     ! print*, ""
!     return
! end subroutine vector_breeding