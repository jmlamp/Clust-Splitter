        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                          THE CLUST-SPLITTER METHOD                               | | 
        !| |                                (version 1.1)                                     | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |                             by Jenni Lampainen                                   | |
        !| |                                                                                  | |
        !| |                         (Last modified 29.4.2025)                                | |
        !| |                                                                                  | |       
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| |     The software is free for academic teaching and research purposes but I       | |
        !| |     ask you to refer the reference given below, if you use it.                   | |
        !| |                                                                                  | |
        !| |                                                                                  | |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..*     
        !|                                                                                      |
        !|                                                                                      |
        !|   Codes include:                                                                     |
        !|                                                                                      |
        !|   clustsplitter.f95       - Main program for Clust-Splitter (this file)              |
        !|   Makefile                - Makefile                                                 |
        !|                                                                                      |
        !|                                                                                      |
        !|                                                                                      |
        !|   To USE the software SET the following PARAMETER VALUES                             |
		!|      at the end of the clustsplitter.f95 file                                        |
        !|                                                                                      |
        !|   used_method                                                                        |
        !|   max_cl                                                                             |
        !|   n_outlier                                                                          |
        !|   delete_outlier                                                                     |
        !|   opt_startingpoint                                                                  |
        !|   allstart                                                                           |
        !|   ncenter1, ncenter2                                                                 |
        !|   index1, index2                                                                     |
        !|   opt1, opt2, opt3                                                                   |
        !|   noopt1, noopt2, noopt3, noopt4, noopt5, noopt6                                     |
        !|                                                                                      |
        !|   Also MODIFY the loop: DO task=... at the end of the clustsplitter.f95 file		    |
        !|      to obtain the correct data set.                                                 |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..* 
        !|                                                                                      |
        !|                                                                                      |
        !|   References:                                                                        |
        !|                                                                                      |
        !|   [1] Lampainen, J., Joki, K., Karmitsa, N., & Mäkelä, M. M. "Clust-Splitter − an    |
		!|       efficient nonsmooth optimization-based algorithm for clustering large          |
		!|       datasets", arXiv:2505.????1 [math.OC], 2025.                                   |
        !|                                                                                      |
        !|                                                                                      |
        !|   [2] N. Haarala, K. Miettinen, M.M. Mäkelä, "Globally Convergent Limited Memory     |
		!|       Bundle Method for Large-Scale Nonsmooth Optimization", Mathematical            |
		!|       Programming, Vol. 109, No. 1, pp. 181-205, 2007. DOI 10.1007/s10107-006-0728-2.|
        !|                                                                                      |
        !|                                                                                      |
        !|   [3] M. Haarala, K. Miettinen, M.M. Mäkelä, "New Limited Memory Bundle Method for   |
        !|       Large-Scale Nonsmooth Optimization", Optimization Methods and Software,        |
		!|       Vol. 19, No. 6, pp. 673-692, 2004. DOI 10.1080/10556780410001689225.	        |
        !|                                                                                      |
        !|                                                                                      |
        !| .**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**. |
        !*..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..**..* 



!======================================================
!
!                  MODULE constants
!
! Module that specifies the precision of real numbers
!======================================================

MODULE constants
    IMPLICIT NONE
    
    ! Double precision (i.e. accuracy)
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
    
END MODULE constants

!======================================================
!
!                  MODULE functions
!
! Module of all functions
!======================================================

MODULE functions

    USE constants
    
    IMPLICIT NONE
    
    INTEGER, SAVE :: nrecords_original                          ! Number of rows in original data a                                                    
    INTEGER, SAVE :: nrecords_a                                 ! Number of rows in modified data a                                                    
    INTEGER, SAVE :: nfeatures_a                                ! Number of colums in data a (original and modified)                                       
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: a       ! Whole data matrix a
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: a_clust         ! Which cluster each point in matrix a belongs to
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, SAVE :: a_error   ! Distances between cluster center and data points (in data a)
          
    INTEGER, SAVE :: nrecords_b                                 ! Number of rows in subdata b                                                 
    INTEGER, SAVE :: nfeatures_b                                ! Number of colums in subdata b                                            
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: b       ! Subdata matrix b
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: b_clust         ! Which cluster each point in submatrix b belongs to
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, SAVE :: b_error   ! Distances between cluster center and data points (in data b)
          
    INTEGER, SAVE :: nclust               ! Number of clusters
    REAL(KIND=dp), SAVE :: coefficient    ! Coefficient for calculating indices
    INTEGER, SAVE :: better_indices       ! 1 = we want to calculate better indeces in 2-clustering auxiliary problem, 0 = we don't want to
    INTEGER, SAVE :: better_indices2      ! 1 = we want to calculate better indeces in k-clustering problem, 0 = we don't want to

    CONTAINS
    
    !---------------------------------------------------
    ! Subroutine where validity indices are calculated
    !---------------------------------------------------
    
        SUBROUTINE check(kl,x_in,db,dn)

            IMPLICIT NONE
            
            ! Arguments
            INTEGER, INTENT(IN) :: kl                                           ! Number of clusters
            REAL(KIND=dp), DIMENSION(nfeatures_a*nclust), INTENT(IN) :: x_in    ! Cluster centers
            REAL(KIND=dp), INTENT(OUT) :: db                                    ! Davies-Bouldin validity index (DBI)
            REAL(KIND=dp), INTENT(OUT) :: dn                                    ! Dunn validity index (DI)
                
            ! Local variables
            REAL(KIND=dp), DIMENSION(nfeatures_a*nclust) :: x   ! Modified cluster centers
            INTEGER, DIMENSION(nclust,nrecords_a) :: nk         ! Help variable
            INTEGER, DIMENSION(nclust) :: nel                   ! Help variable
            INTEGER, DIMENSION(nrecords_a) :: lcand             ! Help variable
            REAL(KIND=dp), DIMENSION(nrecords_a) :: dminim      ! Help variable
            INTEGER :: same_clusters(nclust)                    ! 1 = there are identical clusters, 0 = there are not identical clusters
            REAL(KIND=dp) :: sum1                               ! Help variable
            INTEGER :: modify                                   ! modify > 0 means that there are two or more identical clusters
            INTEGER :: empty_place                              ! Help variable for checking if there are identical clusters

            ! Local real variables
            REAL(KIND=dp) ::      &
                f2,               &
                fdb,              &
                fk2,              &
                fm,               &
                f1(nclust),       &
                fk(nclust),       &
                rad(nclust),      &
                radmax(nclust),   &
                dc(nclust,nclust),&
                dc1,              &
                dn1,              &
                dn2,              &
                radm

            ! Local integer variables
            INTEGER ::  &
                j,      &
                i,      &
                k,      &
                k1,     &
                imin,   &
                mf,     &
                nrecord,&
                ncand,  &
                nc,     &
                list1(nrecords_a)
            
            !--------------------------
            ! Initializations
            !--------------------------
            
            mf = nfeatures_a    
            nrecord = nrecords_a
            nc = kl
            
            ! Check if there are identical clusters
            modify = 0
            DO i=1,nc ! iterate through clusters
                same_clusters(i) = 0
                DO j=i+1,nc ! iterate through other clusters than i
                    IF (same_clusters(i) == 0) THEN
                        sum1 = 0.0_dp
                        DO k=1,nfeatures_a ! iterate through features
                            sum1 = sum1 + (x_in(k+(i-1)*nfeatures_a)-x_in(k+(j-1)*nfeatures_a))**2 ! Calculate the distance between two cluster centers
                        END DO
                        IF (sum1 <= 0.00001_dp) THEN ! clusters i and j are the same
                            PRINT*, 'Clusters', i, ' and', j, ' are identical!'
                            same_clusters(i) = 1
                            modify = modify + 1
                        END IF
                    END IF
                END DO
            END DO
            
            IF (modify > 0) THEN ! Two or more clusters were identical
                empty_place = 0
                DO i=1,nc
                    IF (same_clusters(i)==0) THEN
                        empty_place = empty_place + 1
                        DO j=1,nfeatures_a
                            x(j+(empty_place-1)*nfeatures_a) = x_in(j+(i-1)*nfeatures_a)
                        END DO
                    END IF
                END DO
                nc = nc - modify
            ELSE ! There were no identical clusters
                x = x_in
            END IF

            DO j=1,nclust
                nel(j)=0
                rad(j)=0.0_dp
                radmax(j)=0.0_dp
            END DO

            outerloop: DO k=1,nrecord
                DO j=1,nc
                    f1(j)=0.0_dp
                    DO k1=1,mf
                        f1(j)=f1(j)+(a(k1,k)-x(k1+(j-1)*mf))**2
                    END DO
                END DO
                f2=f1(1)
                imin = 1
                DO  j=2,nc
                    IF (f1(j) < f2) THEN
                        f2=f1(j)
                        imin = j
                    END IF
                END DO
                
                dminim(k)=f2
                nel(imin)=nel(imin)+1
                nk(imin,nel(imin))=k
                list1(k)=imin
                rad(imin)=rad(imin)+f2
                radmax(imin)=MAX(radmax(imin),f2)
            END DO outerloop

            DO k=1,nc
                IF (nel(k) > 0) THEN
                    rad(k)=rad(k)/REAL(nel(k),dp)
                ELSE
                    rad(k)=0.0_dp
                END IF
            END DO

            ncand=0
            DO k=1,nc
                DO j=1,nel(k)
                    i=nk(k,j)
                    IF(dminim(i) > rad(k)) THEN
                        ncand=ncand+1
                        lcand(ncand)=i
                    END IF
                END DO
            END DO
            
            !=====================================================
            ! Calculation of distances between cluster centers
            !=====================================================

            DO i=1,nc
                dc(i,i) = 0.0_dp
            END DO

            DO i=1,nc
                DO j=i+1,nc
                    dc1 = 0.0_dp
                    DO k=1,mf
                        dc1=dc1+(x(k+(i-1)*mf)-x(k+(j-1)*mf))**2                        
                    END DO
                    dc(i,j)=SQRT(dc1)
                    dc(j,i)=dc(i,j)
                END DO
            END DO

            !=====================================================
            ! Calculation of Davies-Bouldin validity index (DBI)
            !=====================================================
            
            DO i=1,nc
                fk(i)=0.0_dp
            END DO

            fdb=0.0_dp

            DO i=1,nc
                fk(i)=SQRT(rad(i))
            END DO

            DO k=1,nc
                fm=0.0_dp
                DO i=1,nc
                    IF (i.ne.k) THEN
                        fk2=fk(i)+fk(k)             
                        f2=fk2/dc(i,k)
                        fm=MAX(fm,f2)
                    END IF
                END DO
                fdb=fdb+fm
                
            END DO
            db=fdb/REAL(nc,dp)
            
            
            !=====================================================
            ! Calculation of Dunn validity index (DI)
            !=====================================================
                        
            radm=0.0_dp
            DO i=1,nc
                radm=MAX(radm,radmax(i)) ! maximum radius of all clusters
            END DO
            radm = SQRT(radm)
            
            dn = 3.40282347*10.**38
            DO i=1,nc
                dn1 = 3.40282347*10.**38
                DO j=1,nc
                    IF (j.NE.i) THEN                       
                        dn2=dc(i,j)/radm
                        dn1=MIN(dn2,dn1) ! distance of cluster i to the closest cluster j
                    END IF
                END DO
                dn = MIN(dn,dn1)
            END DO
            
            IF (nc == 1) dn=0.0_dp

            RETURN

        END SUBROUTINE check
    
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to set whether indices should be considered in the objective function or not
    !------------------------------------------------------------------------------------------
    
        SUBROUTINE consider_indices(ind, ind2)
        
            IMPLICIT NONE
            
            ! Arguments
            INTEGER, INTENT(IN) :: ind          ! 1 = yes, 0 = no
            INTEGER, INTENT(IN) :: ind2         ! 1 = yes, 0 = no
            
            better_indices = ind
            better_indices2 = ind2
        
        END SUBROUTINE consider_indices
    
    
    !---------------------------------------------
    ! Subroutine where coefficient one is set
    !---------------------------------------------
    
        SUBROUTINE set_coefficient(x, new_coefficient, old_coefficient)
        
            IMPLICIT NONE
            
            ! Arguments
            REAL(KIND=dp), DIMENSION(2*nfeatures_a), INTENT(IN) :: x    ! Cluster centers
            INTEGER, INTENT(IN) :: new_coefficient                      ! 1 = The coefficient from the previous iteration was negative;
                                                                            ! A new coefficient is now calculated
                                                                        ! 0 = First iteration; 
                                                                            ! A coefficient is now calculated
            REAL(KIND=dp), INTENT(INOUT) :: old_coefficient             ! Calculated coefficient
            
            ! Local variables
            REAL(KIND=dp) :: func_value             ! Function value
            REAL(KIND=dp) :: cluster_distance       ! Distance between two cluster centers
            INTEGER :: i                            ! Help variable
            
            ! The coefficient from the previous iteration was negative. A new coefficient is now calculated.
            IF (new_coefficient == 1) THEN
                PRINT*, 'Coefficient was negative'
                coefficient = old_coefficient/10
                old_coefficient = coefficient
            ! The first iteration. A coefficient is now calculated.
            ELSE
                ! Calculate the function value
                CALL func_b(x, 2*nfeatures_a, func_value)
                
                ! Calculate distance between two cluster centers
                cluster_distance = 0.0_dp
                DO i=1,nfeatures_a
                    cluster_distance = cluster_distance + (x(i)-x(i+nfeatures_a))**2
                END DO
                
                ! Set the coefficient
                IF (func_value > cluster_distance) THEN
                    IF (cluster_distance <= 0.001_dp) THEN
                        coefficient = func_value
                    ELSE
                        coefficient = func_value / cluster_distance
                        IF (func_value < 10**7) THEN
                            coefficient = coefficient / 1000000
                        ELSE IF (func_value < 10**8) THEN
                            coefficient = coefficient / 100000
                        ELSE
                            coefficient = coefficient / 10000
                        END IF
                    END IF 
                ELSE
                    IF (cluster_distance <= 0.001_dp) THEN
                        coefficient = 1
                    ELSE
                        coefficient = func_value / cluster_distance
                        coefficient = coefficient / 10000
                    END IF
                END IF
                
                ! Store the coefficient
                old_coefficient = coefficient
                
            END IF
        
        END SUBROUTINE set_coefficient
        
        
    !---------------------------------------------
    ! Subroutine where coefficient two is set
    !---------------------------------------------
    
        SUBROUTINE set_coefficient2(x, new_coefficient, old_coefficient)
        
            IMPLICIT NONE
            
            ! Arguments
            REAL(KIND=dp), DIMENSION(nclust*nfeatures_a), INTENT(IN) :: x   ! Cluster centers
            INTEGER, INTENT(IN) :: new_coefficient                          ! 1 = The coefficient from the previous iteration was negative;
                                                                                ! A new coefficient is now calculated
                                                                            ! 0 = First iteration; 
                                                                                ! A coefficient is now calculated
            REAL(KIND=dp), INTENT(INOUT) :: old_coefficient                 ! Calculated coefficient
            
            ! Local variables
            REAL(KIND=dp) :: func_value             ! Function value
            REAL(KIND=dp) :: cluster_distance       ! Distance between two cluster centers
            INTEGER :: i, j, k                      ! Help variables
            REAL(KIND=dp) :: value1                 ! Help variable
            INTEGER :: divisor                      ! Help variable
            
            ! Set the divisor to the number of cluster pairs
            divisor = 0
            DO i=1,nclust
                DO j=i+1,nclust
                    divisor = divisor + 1
                END DO
            END DO
            
            ! Initialize cluster distance to zero
            cluster_distance = 0.0_dp
            
            !The coefficient from the previous iteration was negative. A new coefficient is now calculated.
            IF (new_coefficient == 1) THEN
                PRINT*, 'Coefficient was negative'
                coefficient = old_coefficient/10
                old_coefficient = coefficient
            ! First iteration. A coefficient is now calculated
            ELSE
                ! Calculate the function value
                CALL func_a(x, nclust*nfeatures_a, func_value)
                
                ! Calculate distances between cluster centers
                DO i=1,nclust   ! iterate through clusters
                    DO j=i+1,nclust  ! iterate through other clusters than i
                        value1 = 0.0_dp
                        DO k=1,nfeatures_a  ! iterate through features
                            value1 = value1 + (x(k+nfeatures_a*(i-1))-x(k+nfeatures_a*(j-1)))**2
                        END DO
                        cluster_distance = cluster_distance + value1
                    END DO
                END DO
                
                ! Set the coefficient
                coefficient = func_value / cluster_distance
                coefficient = coefficient / 100
                coefficient = coefficient / divisor      ! scale the coefficient to the number of cluster pairs
                
                ! Store the coefficient
                old_coefficient = coefficient
                
            END IF
        
        END SUBROUTINE set_coefficient2
        
        
    !----------------------------------------------------
    ! Subroutine to calculate radius of cluster (data a)
    !----------------------------------------------------
    
        SUBROUTINE calculate_radius(center, radius)
        
            IMPLICIT NONE
            
            ! Arguments
            REAL(KIND=dp), DIMENSION(nfeatures_a), INTENT(IN) :: center     ! Cluster center
            REAL(KIND=dp), INTENT(OUT) :: radius                            ! Radius of the cluster
            
            ! Local variables
            INTEGER :: i,j
            REAL(KIND=dp) :: help_radius
            
            radius = 0.0_dp
            DO i=1,nrecords_a
                help_radius = 0.0_dp
                DO j=1,nfeatures_a
                    help_radius = help_radius + (center(j)-a(j,i))**2
                END DO
                IF (help_radius > radius) THEN
                    radius = help_radius
                END IF
            END DO
        
        END SUBROUTINE calculate_radius


    !----------------------------------------------------
    ! Subroutine to calculate radius of cluster (data b)
    !----------------------------------------------------
    
        SUBROUTINE calculate_radius2(center, radius)
        
            IMPLICIT NONE
            
            ! Arguments
            REAL(KIND=dp), DIMENSION(nfeatures_a), INTENT(IN) :: center     ! Cluster center
            REAL(KIND=dp), INTENT(OUT) :: radius                            ! Radius of the cluster
            
            ! Local variables
            INTEGER :: i,j
            REAL(KIND=dp) :: help_radius
            
            radius = 0.0_dp
            DO i=1,nrecords_b
                help_radius = 0.0_dp
                DO j=1,nfeatures_b
                    help_radius = help_radius + (center(j)-b(j,i))**2
                END DO
                IF (help_radius > radius) THEN
                    radius = help_radius
                END IF
            END DO
        
        END SUBROUTINE calculate_radius2
    
    
    !--------------------------------------------------
    ! Subroutine where we calculate errors in data a
    !--------------------------------------------------
    
        SUBROUTINE set_errors_a(x)
        
            IMPLICIT NONE
            
            ! Argument
            REAL(KIND=dp), DIMENSION(nfeatures_a), INTENT(IN) :: x   ! Cluster center
            
            ! Local variables
            REAL(KIND=dp) :: error1         ! Help variable
            INTEGER :: i,j                  ! Help variables
            
            ! Set the errors i.e. distance between data points (in data a) and cluster center
            DO i=1,nrecords_a
                error1 = 0.0_dp
                DO j=1,nfeatures_a
                    error1 = error1 + (x(j)-a(j,i))**2
                END DO
                a_error(i) = error1
            END DO
        
        END SUBROUTINE set_errors_a
        
        
    !--------------------------------------------------
    ! Subroutine where we calculate errors in data b
    !--------------------------------------------------
    
        SUBROUTINE set_errors_b(x)
        
            IMPLICIT NONE
            
            ! Argument
            REAL(KIND=dp), DIMENSION(nfeatures_b), INTENT(IN) :: x   ! Cluster center
            
            ! Local variables
            REAL(KIND=dp) :: error1         ! Help variable
            INTEGER :: i,j                  ! Help variables
            
            ! Set the errors i.e. distance between data points (in data b) and cluster center
            DO i=1,nrecords_b
                error1 = 0.0_dp
                DO j=1,nfeatures_b
                    error1 = error1 + (x(j)-b(j,i))**2
                END DO
                b_error(i) = error1
            END DO
        
        END SUBROUTINE set_errors_b
    
    
    !-----------------------------------------------------------
    ! Subroutine where nclust (number of clusters) is set
    !-----------------------------------------------------------
        
        SUBROUTINE set_nclust(nclust0) 
        
            IMPLICIT NONE

            ! Argument
            INTEGER, INTENT(IN) :: nclust0      ! Number of clusters
            
            nclust = nclust0
                                         
        END SUBROUTINE set_nclust
        
        
    !-------------------------------------------------
    ! Subroutine where data matrix a is modified
    !-------------------------------------------------

        SUBROUTINE modify_data(nclmax, examined_ind, n_a_sep, indexes)
        
            IMPLICIT NONE
            
            ! Arguments
            INTEGER, INTENT(IN) :: nclmax                                   ! Maximum number of clusters
            INTEGER, INTENT(IN) :: examined_ind                             ! Number of cluster we want modify
            INTEGER, DIMENSION(nclmax), INTENT(IN) :: n_a_sep               ! Number of points in clusters
            INTEGER, DIMENSION(nclmax,nrecords_a), INTENT(IN) :: indexes    ! Indices that tell which data point belongs to which cluster
            
            ! Local variables
            INTEGER, DIMENSION(n_a_sep(examined_ind)) :: emptyplaces        ! Empty place indexes
            REAL(KIND=dp), DIMENSION(nfeatures_a) :: help_point             ! Help variable for storing a point
            INTEGER :: i, j, ind                                            ! Help variables
            INTEGER :: nchange, empty, npoints                              ! Help variables
            INTEGER :: old_place, new_place                                 ! Help variables for places
            
            ! Set the number of outliers
            npoints = n_a_sep(examined_ind)
            
            ! Determine the indices of the outlier points
            nchange = npoints
            DO i=1,npoints
                empty = 0
                DO j=1,npoints
                    IF (indexes(examined_ind,j) == nrecords_a+i-npoints) THEN
                        empty = empty+1
                        nchange = nchange-1
                    END IF
                END DO
                IF (empty == 0) THEN
                    emptyplaces(i) = nrecords_a+i-npoints
                ELSE
                    emptyplaces(i) = 0  ! Place is full
                END IF
            END DO
        
            ! Rearrange the data so that outlier points are moved to the end
            ind = 1
            DO i=1,nchange
                old_place = indexes(examined_ind,i)
                new_place = emptyplaces(ind)
                DO WHILE (new_place == 0)
                    ind = ind + 1
                    new_place = emptyplaces(ind)
                END DO
                help_point = a(:,old_place)
                a(:,old_place)=a(:,new_place)
                a(:,new_place)=help_point
                ind = ind + 1
            END DO
            
            ! Update records in data set a
            nrecords_a = nrecords_a - npoints
            
        END SUBROUTINE modify_data
        
        
    !----------------------------------------------------------------------------
    ! Subroutine where data is read into the matrix a and matrix b is allocated
    !----------------------------------------------------------------------------

        SUBROUTINE allocate_whole_data(infile, nrecords, nfeatures)
        
            IMPLICIT NONE
            
            ! Arguments
            CHARACTER(LEN=*), INTENT(IN) :: infile      ! Name of the dataset file
            INTEGER, INTENT(IN) :: nrecords             ! Number of rows in data
            INTEGER, INTENT(IN) :: nfeatures            ! Number of colums in data
            
            ! Local variables
            REAL(KIND=dp), DIMENSION(nrecords,nfeatures) :: a_in
            INTEGER :: i, j     
               
            nrecords_original = nrecords
            nrecords_a = nrecords
            nfeatures_a = nfeatures
            nrecords_b = nrecords
            nfeatures_b = nfeatures
            
            ! Allocate data matrix a
            ALLOCATE(a(nfeatures, nrecords))   
            
            ! Allocate a_clust
            ALLOCATE(a_clust(nrecords))
            
            ! Allocate a_error
            ALLOCATE(a_error(nrecords))
            
            ! Allocate data matrix b
            ALLOCATE(b(nfeatures, nrecords))  
            
            ! Allocate b_clust
            ALLOCATE(b_clust(nrecords))
            
            ! Allocate b_error
            ALLOCATE(b_error(nrecords))
              
            ! Open and read data from the file
            OPEN(10, file=infile, status='old', form='formatted')
            DO i = 1, nrecords
                READ(10, *) (a_in(i, j), j = 1, nfeatures)
            END DO
            CLOSE(10)

            ! Fill matrix a
            DO i = 1, nrecords
                DO j = 1, nfeatures
                    a(j,i) = a_in(i,j)
                END DO
            END DO
                                         
        END SUBROUTINE allocate_whole_data
        
        
    !-----------------------------------------------------------------------------------------------------
    ! Subroutine where submatrix b is formed (b is matrix containing the data points of a single cluster)
    !-----------------------------------------------------------------------------------------------------

        SUBROUTINE set_submatrix(ind_clust)
        
            IMPLICIT NONE
            
            ! Arguments
            INTEGER, INTENT(IN) :: ind_clust                ! Cluster number to focus on

            ! Local variables
            INTEGER :: i, j, k, row                         ! Loop variables and row counter

            ! Count how many rows have the target cluster
            k = 0
            DO i = 1, nrecords_a
                IF (a_clust(i) == ind_clust) THEN
                    k = k + 1
                END IF
            END DO
            
            ! Set the number of records in data b
            nrecords_b = k

            ! Fill submatrix b
            row = 0
            DO i = 1, nrecords_a
                IF (a_clust(i) == ind_clust) THEN
                    row = row + 1
                    DO j = 1, nfeatures_a
                        b(j, row) = a(j, i) ! Copy all columns from matrix a to submatrix b
                    END DO
                END IF
            END DO
            
        END SUBROUTINE set_submatrix
        
        
    !------------------------------------------------------------------------------
    ! Subroutine to deallocate data matrices a and b along with related resources
    !------------------------------------------------------------------------------
        
        SUBROUTINE deallocate_whole_data() 
        
            IMPLICIT NONE
            
            DEALLOCATE(a)
            DEALLOCATE(a_clust)
            DEALLOCATE(a_error)
            DEALLOCATE(b)
            DEALLOCATE(b_clust)
            DEALLOCATE(b_error)
                                         
        END SUBROUTINE deallocate_whole_data
        
    
    !-------------------------------------------------------------------------------------------------------------------
    ! Subroutine where the clustering function f value is calculated from the data a (each cluster gets its own value)
    !-------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_a_sep(x, n, ncl, nclmax, f_a_sep, n_a_sep, indexes)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                    ! Cluster centers
        INTEGER, INTENT(IN) :: n                                        ! Number of clusters * number of features
        INTEGER, INTENT(IN) :: ncl                                      ! Number of clusters
        INTEGER, INTENT(IN) :: nclmax                                   ! Maximum number of clusters
        REAL(KIND=dp), DIMENSION(nclmax), INTENT(OUT) :: f_a_sep        ! Function value
        INTEGER, DIMENSION(nclmax), INTENT(OUT) :: n_a_sep              ! Number of points in clusters
        INTEGER, DIMENSION(nclmax,nrecords_a), INTENT(OUT) :: indexes   ! Data point indexes in each cluster
        
        ! Local variables
        INTEGER :: i, j, k, num                     ! Help variables
        REAL(KIND=dp) :: smallest_value, value1     ! Help variables
    
        ! Initialize the values of all clustering functions to 0 and number of points in clusters to 0
        f_a_sep = 0.0_dp
        n_a_sep = 0
        
        ! Calculating the value of the function
        DO i=1,nrecords_a  ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            num = 1
            DO j=1,ncl  ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a  ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    a_clust(i) = j  ! Determine which cluster the data point belongs to
                    num = j
                END IF
            END DO
            n_a_sep(num) = n_a_sep(num) + 1 
            f_a_sep(num) = f_a_sep(num)+smallest_value
            indexes(num,n_a_sep(num)) = i
        END DO

    END SUBROUTINE func_a_sep
    
    !--------------------------------------------------------------------------------------------------------------------
    ! Subroutine where the clustering function f_b value is calculated from the data b (each cluster gets its own value)
    !--------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_b_sep(x, n, ncl, f_b_sep)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        INTEGER, INTENT(IN) :: ncl                                  ! Number of clusters
        REAL(KIND=dp), DIMENSION(ncl), INTENT(OUT) :: f_b_sep       ! Function value
        
        ! Local variables
        INTEGER :: i, j, k, num                                     ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_b_sep = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            num = 1
            DO j=1,ncl                                          ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    b_clust(i) = j  ! Determine which cluster the data point belongs to
                    num = j
                END IF
            END DO
            f_b_sep(num) = f_b_sep(num)+smallest_value
        END DO

    END SUBROUTINE func_b_sep
    

    !-----------------------------------------------------------------------------------------------------------------------
    ! Subroutine where the value of the clustering function f is calculated and a_clust is utilized (one value is obtained)
    !-----------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_a_clust(x, n, f_a)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_a                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k, num                                     ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_a = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_a                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            num = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    a_clust(i) = j  ! Determine which cluster the data point belongs to
                    num = j
                END IF
            END DO
            f_a = f_a+smallest_value
        END DO

    END SUBROUTINE func_a_clust
    
    
    !-----------------------------------------------------------------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the value of the clustering function f is calculated and a_clust is utilized (one value is obtained)
    !-----------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_a_clust_index(x, n, f_a)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_a                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k, num                                     ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_a = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_a                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            num = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    a_clust(i) = j  ! Determine which cluster the data point belongs to
                    num = j
                END IF
            END DO
            f_a = f_a+smallest_value
        END DO
        
        ! For better indices
        value1 = 0.0_dp
        DO i=1,nclust   ! iterate through clusters
            DO j=i+1,nclust  ! iterate through other clusters than i
                DO k=1,nfeatures_a  ! iterate through features
                    value1 = value1 + (x(k+nfeatures_a*(i-1))-x(k+nfeatures_a*(j-1)))**2
                END DO
            END DO
        END DO
        f_a = f_a - coefficient*value1

    END SUBROUTINE func_a_clust_index
    
    
    !-------------------------------------------------------------------------------------------------------------------------
    ! Subroutine where the value of the clustering function f_b is calculated and b_clust is utilized (one value is obtained)
    !-------------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_b_clust(x, n, f_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_b                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k, num                                     ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_b = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            num = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    b_clust(i) = j  ! Determine which cluster the data point belongs to
                    num = j
                END IF
            END DO
            f_b = f_b+smallest_value
        END DO

    END SUBROUTINE func_b_clust
    

    !-------------------------------------------------------------------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the value of the clustering function f_b is calculated and b_clust is utilized (one value is obtained)
    !-------------------------------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_b_clust_index(x, n, f_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_b                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k, num                                     ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                               ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_b = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            num = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    b_clust(i) = j  ! Determine which cluster the data point belongs to
                    num = j
                END IF
            END DO
            f_b = f_b+smallest_value
        END DO
        
        ! For better indices
        value1 = 0.0_dp
        DO k=1,nfeatures_b
            value1 = value1 - ((x(k)-x(k+nfeatures_b))**2)
        END DO
        f_b = f_b+value1*coefficient

    END SUBROUTINE func_b_clust_index
    
    
    !-------------------------------------------------------------------------------------------------
    ! Subroutine where the value of the clustering function f is calculated (one value is obtained)
    !-------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_a(x, n, f_a)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_a                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k                                          ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_a = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_a                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                END IF
            END DO
            f_a = f_a+smallest_value
        END DO

    END SUBROUTINE func_a
    
    
    !-------------------------------------------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the value of the clustering function f is calculated (one value is obtained)
    !-------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_a_index(x, n, f_a)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_a                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k                                          ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_a = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_a                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                END IF
            END DO
            f_a = f_a+smallest_value
        END DO
        
        ! For better indices
        value1 = 0.0_dp
        DO i=1,nclust   ! iterate through clusters
            DO j=i+1,nclust  ! iterate through other clusters than i
                DO k=1,nfeatures_a  ! iterate through features
                    value1 = value1 + (x(k+nfeatures_a*(i-1))-x(k+nfeatures_a*(j-1)))**2
                END DO
            END DO
        END DO
        f_a = f_a - coefficient*value1
        

    END SUBROUTINE func_a_index
    
    
    !---------------------------------------------------------------------------------------------
    ! Subroutine to calculate function value at the starting point auxiliary problem (for data a)
    !---------------------------------------------------------------------------------------------
    
    SUBROUTINE func_help_a(x_new, n, f_help)
    
        IMPLICIT NONE
        
        ! Arguments
        INTEGER, INTENT(IN) :: n                              ! Dimension, nfeatures_a
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x_new      ! New cluster center
        REAL(KIND=dp), INTENT(OUT) :: f_help                  ! Function value
        
        ! Local variables
        INTEGER :: i, j                                       ! Help variables
        REAL(KIND=dp) :: f, small1                            ! Help variables
    
        ! Initialize function value to 0
        f_help = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_a ! Iterate through the data points
            f = 0.0_dp
            DO j=1,nfeatures_a  ! Iterate through the features
                f = f + (x_new(j)-a(j,i))**2  ! Calculate the distance between cluster center and data point
            END DO
            ! Determine which distance is smaller
                ! the distance between the data point and the original cluster center or
                ! the distance between the data point and the new cluster center
            small1 = MIN(a_error(i),f)
            f_help = f_help + small1 ! Add the smaller distance to the function value
        END DO

    END SUBROUTINE func_help_a
    
    
    !---------------------------------------------------------------------------------------------
    ! Subroutine to calculate function value at the starting point auxiliary problem (for data b)
    !---------------------------------------------------------------------------------------------
    
    SUBROUTINE func_help_b(x_new, n, f_help)
    
        IMPLICIT NONE
        
        ! Arguments
        INTEGER, INTENT(IN) :: n                              ! Dimension, nfeatures_b
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x_new      ! New cluster center
        REAL(KIND=dp), INTENT(OUT) :: f_help                  ! Function value
        
        ! Local variables
        INTEGER :: i, j                                       ! Help variables
        REAL(KIND=dp) :: f, small1                            ! Help variables
     
        ! Initialize function value to 0
        f_help = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_b  ! Iterate through the data points
            f = 0.0_dp
            DO j=1,nfeatures_b  ! Iterate through the features
                f = f + (x_new(j)-b(j,i))**2  ! Calculate the distance between cluster center and data point
            END DO
            ! Determine which distance is smaller
                ! the distance between the data point and the original cluster center or
                ! the distance between the data point and the new cluster center
            small1 = MIN(b_error(i),f)
            f_help = f_help + small1  ! Add the smaller distance to the function value
        END DO

    END SUBROUTINE func_help_b
    
    
    !-------------------------------------------------------------------------------------------------
    ! Subroutine to calculate subgradient value at the starting point auxiliary problem (for data a)
    !-------------------------------------------------------------------------------------------------  
    
    SUBROUTINE subgrad_help_a(x_new, n, subgrad)
    
        IMPLICIT NONE
    
        ! Arguments
        INTEGER, INTENT(IN) :: n                                    ! Dimension, nfeatures_a
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x_new            ! New cluster center
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: subgrad         ! Subgradient for x_new

        ! Local variables
        INTEGER :: i, j                             ! Iteration variables
        REAL(KIND=dp) :: f                          ! Function value
        REAL(KIND=dp), DIMENSION(n) :: grad         ! Gradient

        ! Initialize subgradient value to zero
        subgrad = 0.0_dp

        ! Loop over each data point
        DO i = 1, nrecords_a
            ! Initializations
            f = 0.0_dp
            grad = 0.0_dp
            ! Calculate function value and gradient
            DO j = 1, nfeatures_a
                f = f + (x_new(j) - a(j,i))**2
                grad(j) = 2.0_dp * (x_new(j) - a(j,i))
            END DO
            ! Add gradient value based on the minimum condition
            IF (f <= a_error(i)) THEN
                subgrad = subgrad + grad
            END IF
        END DO

    END SUBROUTINE subgrad_help_a
    
    
    !-------------------------------------------------------------------------------------------------
    ! Subroutine to calculate subgradient value at the starting point auxiliary problem (for data b)
    !-------------------------------------------------------------------------------------------------  
    
    SUBROUTINE subgrad_help_b(x_new, n, subgrad)
    
        IMPLICIT NONE

        ! Arguments
        INTEGER, INTENT(IN) :: n                                    ! Dimension, nfeatures_b
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x_new            ! New cluster center
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: subgrad         ! Subgradient for x_new

        ! Local variables
        INTEGER :: i, j                             ! Iteration variables
        REAL(KIND=dp) :: f                          ! Function value
        REAL(KIND=dp), DIMENSION(n) :: grad         ! Gradient

        ! Initialize subgradient values to zero
        subgrad = 0.0_dp

        ! Loop over each data point
        DO i = 1, nrecords_b
            ! Initializations
            f = 0.0_dp
            grad = 0.0_dp
            ! Calculate function value and gradient
            DO j = 1, nfeatures_b
                f = f + (x_new(j) - a(j,i))**2
                grad(j) = 2.0_dp * (x_new(j) - b(j,i))
            END DO
            ! Add gradient value based on the minimum condition
            IF (f <= b_error(i)) THEN
                subgrad = subgrad + grad
            END IF
        END DO

    END SUBROUTINE subgrad_help_b
    
    
    !-------------------------------------------------------------------------------------------------
    ! Subroutine where the value of the clustering function f is calculated (one value is obtained)
    !-------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_a_original(x, n, f_a)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_a                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k                                          ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_a = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_original                                ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                END IF
            END DO
            f_a = f_a+smallest_value
        END DO

    END SUBROUTINE func_a_original
    
    
    !-------------------------------------------------------------------------------------------------
    ! Subroutine where the value of the clustering function f_b is calculated (one value is obtained)
    !-------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_b(x, n, f_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_b                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k                                          ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_b = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                END IF
            END DO
            f_b = f_b+smallest_value
        END DO

    END SUBROUTINE func_b
    
    
    !-------------------------------------------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the value of the clustering function f_b is calculated (one value is obtained)
    !-------------------------------------------------------------------------------------------------
    
    SUBROUTINE func_b_index(x, n, f_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), INTENT(OUT) :: f_b                           ! Function value
        
        ! Local variables
        INTEGER :: i, j, k                                          ! Help variables
        REAL(KIND=dp) :: smallest_value, value1                     ! Help variables
    
        ! Initialize the values of all clustering functions to 0
        f_b = 0.0_dp
        
        ! Calculating the value of the function
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                END IF
            END DO
            f_b = f_b+smallest_value
        END DO
        
        ! For better indices
        value1 = 0.0_dp
        DO k=1,nfeatures_b
            value1 = value1 - ((x(k)-x(k+nfeatures_b))**2)
        END DO
        f_b = f_b+value1*coefficient

    END SUBROUTINE func_b_index
    

    !--------------------------------------------------------------------------------------
    ! Subroutine where the subgradient g is calculated (from data a and utilizing a_clust)
    !--------------------------------------------------------------------------------------
    
    SUBROUTINE subgrad_a_clust(x, n, subgrad_a)
    
        IMPLICIT NONE
    
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: subgrad_a       ! Subgradient
        
        ! Local variables
        INTEGER :: i, j, t                                          ! Help variables
        REAL(KIND=dp) :: gvalue1                                    ! Help variable
    
        ! Initialize the values of all subgradients to 0
        subgrad_a = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_a                                       ! Iterate through the data points
            DO j=1,nfeatures_a                                  ! Iterate through the features
                t = (a_clust(i)-1)*nfeatures_a
                gvalue1 = 2*(x(j+t)-a(j,i))                     ! Calculate the subgradient
                subgrad_a(j+t) = subgrad_a(j+t)+gvalue1         ! Store the value of the subgradient
            END DO
        END DO

    END SUBROUTINE subgrad_a_clust
    
    
    !--------------------------------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the subgradient g is calculated (from data a and utilizing a_clust)
    !--------------------------------------------------------------------------------------
    
    SUBROUTINE subgrad_a_clust_index(x, n, subgrad_a)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: subgrad_a       ! Subgradient
        
        ! Local variables
        REAL(KIND=dp), DIMENSION(n) :: subgrad_a2                   ! Help subgradient
        INTEGER :: i, j, k, t                                       ! Help variables
        REAL(KIND=dp) :: gvalue1                                    ! Help variable
    
        ! Initialize the values of all subgradients to 0
        subgrad_a = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_a                                       ! Iterate through the data points
            DO j=1,nfeatures_a                                  ! Iterate through the features
                t = (a_clust(i)-1)*nfeatures_a
                gvalue1 = 2*(x(j+t)-a(j,i))                     ! Calculate the subgradient
                subgrad_a(j+t) = subgrad_a(j+t)+gvalue1         ! Store the value of the subgradient
            END DO
        END DO
        
        ! For better indices
        subgrad_a2 = 0.0_dp
        DO i=1,nclust   ! iterate through clusters
            DO j=i+1,nclust  ! iterate through other clusters than i
                DO k=1,nfeatures_a  ! iterate through features
                    gvalue1 = -2*(x(k+nfeatures_a*(i-1))-x(k+nfeatures_a*(j-1)))
                    subgrad_a2(k+nfeatures_a*(i-1)) = subgrad_a2(k+nfeatures_a*(i-1)) + gvalue1
                    gvalue1 = -2*(x(k+nfeatures_a*(j-1))-x(k+nfeatures_a*(i-1)))
                    subgrad_a2(k+nfeatures_a*(j-1)) = subgrad_a2(k+nfeatures_a*(j-1)) + gvalue1
                END DO
            END DO
        END DO
        
        subgrad_a2 = subgrad_a2 * coefficient
        subgrad_a = subgrad_a + subgrad_a2

    END SUBROUTINE subgrad_a_clust_index
    
    
    !--------------------------------------------------------------------------------------
    ! Subroutine where the subgradient g is calculated (from data b and utilizing b_clust)
    !--------------------------------------------------------------------------------------
    
    SUBROUTINE subgrad_b_clust(x, n, subgrad_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: subgrad_b       ! Subgradient
        
        ! Local variables
        INTEGER :: i, j, t                                          ! Help variables
        REAL(KIND=dp) :: gvalue1                                    ! Help variable
    
        ! Initialize the values of all subgradients to 0
        subgrad_b = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_b                                       ! Iterate through the data points
            DO j=1,nfeatures_b                                  ! Iterate through the features
                t = (b_clust(i)-1)*nfeatures_b
                gvalue1 = 2*(x(j+t)-b(j,i))                     ! Calculate the subgradient
                subgrad_b(j+t) = subgrad_b(j+t)+gvalue1         ! Store the value of the subgradient
            END DO
        END DO

    END SUBROUTINE subgrad_b_clust
    
    
    !--------------------------------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the subgradient g is calculated (from data b and utilizing b_clust)
    !--------------------------------------------------------------------------------------
    
    SUBROUTINE subgrad_b_clust_index(x, n, subgrad_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x                ! Cluster centers
        INTEGER, INTENT(IN) :: n                                    ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: subgrad_b       ! Subgradient
        
        ! Local variables
        INTEGER :: i, j, k, t                                       ! Help variables
        REAL(KIND=dp) :: gvalue1                                    ! Help variable
    
        ! Initialize the values of all subgradients to 0
        subgrad_b = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_b                                       ! Iterate through the data points
            DO j=1,nfeatures_b                                  ! Iterate through the features
                t = (b_clust(i)-1)*nfeatures_b
                gvalue1 = 2*(x(j+t)-b(j,i))                     ! Calculate the subgradient
                subgrad_b(j+t) = subgrad_b(j+t)+gvalue1         ! Store the value of the subgradient
            END DO
        END DO
        
        ! For better indices
        DO k=1,nfeatures_b
            gvalue1 = -2*coefficient*(x(k)-x(k+nfeatures_b))
            subgrad_b(k) = subgrad_b(k) + gvalue1
            gvalue1 = -2*coefficient*(x(k+nfeatures_b)-x(k))
            subgrad_b(k+nfeatures_b) = subgrad_b(k+nfeatures_b) + gvalue1
        END DO

    END SUBROUTINE subgrad_b_clust_index
    
    
    !----------------------------------------------------------------
    ! Subroutine where the subgradient g is calculated (from data a)
    !----------------------------------------------------------------
    
    SUBROUTINE subgrad_a(x, n, grad_a)
    
        IMPLICIT NONE

        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x            ! Cluster centers
        INTEGER, INTENT(IN) :: n                                ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: grad_a      ! Subgradient
        
        ! Local variables
        INTEGER :: i, j, k, ii, ind1, ind2                      ! Help variables
        REAL(KIND=dp) :: gvalue1, value1, smallest_value        ! Help variables
    
        ! Initialize the values of all subgradients to 0
        grad_a = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_a                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            ind1 = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    ind1 = j
                END IF
            END DO
            DO ii=1,nfeatures_a                         ! Iterate through the features
                ind2 = ii+nfeatures_a*(ind1-1)
                gvalue1 = 2*(x(ind2)-a(ii,i))           ! Calculate the subgradient
                grad_a(ind2) = grad_a(ind2)+gvalue1     ! Store the value of the subgradient
            END DO
        END DO

    END SUBROUTINE subgrad_a
    
    
    !----------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the subgradient g is calculated (from data a)
    !----------------------------------------------------------------
    
    SUBROUTINE subgrad_a_index(x, n, grad_a)
    
        IMPLICIT NONE

        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x            ! Cluster centers
        INTEGER, INTENT(IN) :: n                                ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: grad_a      ! Subgradient
        
        ! Local variables
        REAL(KIND=dp), DIMENSION(n) :: grad_a2                  ! Help subgradient
        INTEGER :: i, j, k, ii, ind1, ind2                      ! Help variables
        REAL(KIND=dp) :: gvalue1, value1, smallest_value        ! Help variables
    
        ! Initialize the values of all subgradients to 0
        grad_a = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_a                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            ind1 = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_a                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_a*(j-1))-a(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    ind1 = j
                END IF
            END DO
            DO ii=1,nfeatures_a                         ! Iterate through the features
                ind2 = ii+nfeatures_a*(ind1-1)
                gvalue1 = 2*(x(ind2)-a(ii,i))           ! Calculate the subgradient
                grad_a(ind2) = grad_a(ind2)+gvalue1     ! Store the value of the subgradient
            END DO
        END DO
        
        ! For better indices
        grad_a2 = 0.0_dp
        DO i=1,nclust   ! iterate through clusters
            DO j=i+1,nclust  ! iterate through other clusters than i
                DO k=1,nfeatures_a  ! iterate through features
                    gvalue1 = -2*(x(k+nfeatures_a*(i-1))-x(k+nfeatures_a*(j-1)))
                    grad_a2(k+nfeatures_a*(i-1)) = grad_a2(k+nfeatures_a*(i-1)) + gvalue1
                    gvalue1 = -2*(x(k+nfeatures_a*(j-1))-x(k+nfeatures_a*(i-1)))
                    grad_a2(k+nfeatures_a*(j-1)) = grad_a2(k+nfeatures_a*(j-1)) + gvalue1
                END DO
            END DO
        END DO
        
        grad_a2 = grad_a2 * coefficient
        grad_a = grad_a + grad_a2

    END SUBROUTINE subgrad_a_index


    !----------------------------------------------------------------
    ! Subroutine where the subgradient g is calculated (from data b)
    !----------------------------------------------------------------

    SUBROUTINE subgrad_b(x, n, grad_b)
    
        IMPLICIT NONE
        
        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x            ! Cluster centers
        INTEGER, INTENT(IN) :: n                                ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: grad_b      ! Subgradient
        
        ! Loca variables
        INTEGER :: i, j, k, ii, ind1, ind2                      ! Help variables
        REAL(KIND=dp) :: gvalue1, value1, smallest_value        ! Help variables
    
        ! Initialize the values of all subgradients to 0
        grad_b = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            ind1 = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    ind1 = j
                END IF
            END DO
            DO ii=1,nfeatures_b                         ! Iterate through the features
                ind2 = ii+nfeatures_b*(ind1-1)
                gvalue1 = 2*(x(ind2)-b(ii,i))           ! Calculate the subgradient
                grad_b(ind2) = grad_b(ind2)+gvalue1     ! Store the value of the subgradient
            END DO
        END DO

    END SUBROUTINE subgrad_b
    
    
    !----------------------------------------------------------------
    ! FOR BETTER INDICES
    ! Subroutine where the subgradient g is calculated (from data b)
    !----------------------------------------------------------------

    SUBROUTINE subgrad_b_index(x, n, grad_b)
    
        IMPLICIT NONE

        ! Arguments
        REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x            ! Cluster centers
        INTEGER, INTENT(IN) :: n                                ! Number of clusters * number of features
        REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: grad_b      ! Subgradient
        
        ! Local variables
        INTEGER :: i, j, k, ii, ind1, ind2                      ! Help variables
        REAL(KIND=dp) :: gvalue1, value1, smallest_value        ! Help variables
    
        ! Initialize the values of all subgradients to 0
        grad_b = 0.0_dp
        
        ! Calculating the subgradient
        DO i=1,nrecords_b                                       ! Iterate through the data points
            smallest_value = 3.40282347*10.**38
            ind1 = 1
            DO j=1,nclust                                       ! Iterate through the clusters
                value1 = 0.0_dp
                DO k=1,nfeatures_b                              ! Iterate through the features
                    value1 = value1 + (x(k+nfeatures_b*(j-1))-b(k,i))**2
                END DO
                IF (value1 < smallest_value) THEN
                    smallest_value = value1
                    ind1 = j
                END IF
            END DO
            DO ii=1,nfeatures_b                         ! Iterate through the features
                ind2 = ii+nfeatures_b*(ind1-1)
                gvalue1 = 2*(x(ind2)-b(ii,i))           ! Calculate the subgradient
                grad_b(ind2) = grad_b(ind2)+gvalue1     ! Store the value of the subgradient
            END DO
        END DO
        
        ! For better indices
        DO k=1,nfeatures_b
            gvalue1 = -2*coefficient*(x(k)-x(k+nfeatures_b))
            grad_b(k) = grad_b(k) + gvalue1
            gvalue1 = -2*coefficient*(x(k+nfeatures_b)-x(k))
            grad_b(k+nfeatures_b) = grad_b(k+nfeatures_b) + gvalue1
        END DO

    END SUBROUTINE subgrad_b_index
    
    
END MODULE functions


!=======================================================================================
!
!                          *************  LMBM  *************
!
!=======================================================================================

MODULE param  ! Parameters
  USE constants
  IMPLICIT NONE

! Parameters
  INTEGER, PARAMETER, PUBLIC :: maxeps = 20, maxnrs = 2000
  REAL(KIND=dp), PARAMETER, PUBLIC :: &
       zero    = 0.0_dp,    & ! 
       half    = 0.5_dp,    & ! 
       one     = 1.0_dp,    & ! 
       large   = 3.40282347*10.**38,  & !
       small   = 1.17549435*10.**(-38)     ! 

END MODULE param


MODULE exe_time  ! Execution time
  USE constants
  IMPLICIT NONE

  PUBLIC :: getime

CONTAINS
  SUBROUTINE getime(tim)  ! Execution time.
  IMPLICIT NONE
      
! Scalar Arguments
    REAL(KIND=dp), INTENT(OUT):: tim  ! Current time, REAL argument.

! Intrinsic Functions
    INTRINSIC CPU_TIME

    CALL CPU_TIME(tim)

  END SUBROUTINE getime

END MODULE exe_time


MODULE lmbm_sub  ! Subprograms for lmbm

  USE constants
  IMPLICIT NONE

! MODULE lmbm_sub includes the following subroutines (S) and functions (F).
  PUBLIC :: &
       vdot, &   ! F Dot product of two vectors.
       vneg, &   ! S Change the signs of vector elements.
       copy, &   ! S Copying a vector.
       copy2, &  ! S Copying two vectors.
       xdiffy, & ! S Difference of two vectors z:= x - y.
       xdiffy2, &! S Difference of two vectors z:= x - y. (Variant)
       xsumy, &  ! S Sum of two vectors z:= x + y.
       xsumy2, & ! S Sum of two vectors z:= x + y. (Variant)
       scdiff, & ! S Difference of the scaled vector and a vector z:= a*x - y.
       scdiff2, &! S Difference of the scaled vector and a vector z:= a*x - y. (Variant)
       scsum, &  ! S Sum of a vector and the scaled vector z:= y + a*x.
       scsum2, & ! S Sum of a vector and the scaled vector y:= y + a*x. (variant INOUT for y)
       vxdiag, & ! S Vector is multiplied by a diagonal matrix y:=d*x.
       symax, &  ! S Multiplication of a dense symmetric matrix A by a vector x.
       cwmaxv, & ! S Multiplication of a columnwise stored dense rectangular matrix by a vector.
       rwaxv2, & ! S Multiplication of two rowwise stored dense rectangular  
                 !   matrices A and B by vectors X and Y.
       trlieq, & ! S Solving x from linear equation u*x=y or u'*x=y, 
                 !   where u is an upper triangular matrix.
       trlieq2, &! S Solving x from linear equation u*x=y or u'*x=y, (variant)
                 !   where u is an upper triangular matrix.
       lineq, &  ! S Solver from linear equation.
       calq      ! S Solving x from linear equation A*x=y. Contains:
                 !     S mxdpgf   Gill-Murray decomposition of a dense symmetric matrix.

CONTAINS

  FUNCTION vdot(n,x,y) RESULT(xty) ! Dot product of two vectors.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    REAL(KIND=dp) xty
    INTEGER :: i

    xty = zero
    DO i = 1,n
       xty = xty + x(i)*y(i)
    END DO

  END FUNCTION vdot

  SUBROUTINE vneg(n,x,y)  ! Change the signs of vector elements.
    IMPLICIT NONE
      
! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= -x.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER :: i

    DO i = 1,n
       y(i) = -x(i)
    END DO

  END SUBROUTINE vneg

  SUBROUTINE scalex(n,a,x,y)  ! Scaling a vector y:= a*x.
    IMPLICIT NONE
      
! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= a*x.

! Scalar Arguments
    REAL(KIND=dp), INTENT(IN) :: &
         a           ! Scaling parameter.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER :: i

    DO i = 1,n
       y(i) = a*x(i)
    END DO
      
  END SUBROUTINE scalex

  SUBROUTINE xdiffy(n,x,y,z)  ! Difference of two vectors z:= x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= x - y.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = x(i) - y(i)
    END DO
 
  END SUBROUTINE xdiffy

  ! Variant with x:= x - y
!  SUBROUTINE xdiffy(n,x,y,z)  ! Difference of two vectors z:= x - y.
  SUBROUTINE xdiffy2(n,x,y)  ! Difference of two vectors x:= x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
!         x,y         ! Input vectors.
         y         ! Input vectors.
!    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
!         z           ! Output vector z:= x - y.
         x           ! Output vector x:= x - y.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
!       z(i) = x(i) - y(i)
       x(i) = x(i) - y(i)
    END DO
 
  END SUBROUTINE xdiffy2

  SUBROUTINE xsumy(n,x,y,z)  ! Sum of two vectors z:= x + y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= x + y.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = x(i) + y(i)
    END DO

  END SUBROUTINE xsumy

  ! Variant with x = z
  SUBROUTINE xsumy2(n,x,y)  ! Sum of two vectors x:= x + y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
!         x,y         ! Input vectors.
         y           ! Input vectors.
!    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
!         z           ! Output vector z:= x + y.
         x           ! Output vector x:= x + y.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
!       z(i) = x(i) + y(i)
       x(i) = x(i) + y(i)
    END DO

  END SUBROUTINE xsumy2

  SUBROUTINE scdiff(n,a,x,y,z)  ! Difference of the scaled vector and a vector z:= a*x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= a*x - y.

! Scalar Arguments
    REAL(KIND=dp), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = a*x(i) - y(i)
    END DO
 
  END SUBROUTINE scdiff

  ! Variant subroutine for cases where x = z when calling the subroutine
!  SUBROUTINE scdiff(n,a,x,y,z)  ! Difference of the scaled vector and a vector z:= a*x - y.
  SUBROUTINE scdiff2(n,a,x,y)  ! Difference of the scaled vector and a vector x:= a*x - y.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
!         x,y         ! Input vectors.
         y           ! Input vectors.
!    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         x           ! Output vector x:= a*x - y. (variant)

! Scalar Arguments
    REAL(KIND=dp), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       x(i) = a*x(i) - y(i)
    END DO
 
  END SUBROUTINE scdiff2


  SUBROUTINE scsum(n,a,x,y,z)  ! Sum of a vector and the scaled vector z:= y + a*x.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,y         ! Input vectors.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         z           ! Output vector z:= a*x + y.

! Scalar Arguments
    REAL(KIND=dp), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       z(i) = a*x(i) + y(i)
    END DO

  END SUBROUTINE scsum

! Variant for subroutine scsum with one paramameter as INOUT

  SUBROUTINE scsum2(n,a,x,y)  ! Sum of a vector and the scaled vector y:= y + a*x. 
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x         ! Input vectors.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         y           ! Output vector z:= a*x + y.

! Scalar Arguments
    REAL(KIND=dp), INTENT(IN) :: &
         a           ! Scaling factor.
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       y(i) = a*x(i) + y(i)
    END DO

  END SUBROUTINE scsum2

  SUBROUTINE copy(n,x,y)  ! Copying of vector Y:= X.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x           ! Input vector.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER :: i

    DO i = 1,n
       y(i) = x(i)
    END DO
      
  END SUBROUTINE copy

  SUBROUTINE copy2(n,x,y,z,v)  ! Copying of two vectors: y:=x, v:=z.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,z         ! Input vectors.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y,v         ! Output vectors.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.

! Local Scalars
    INTEGER :: i

    DO i = 1,n
       y(i) = x(i)
       v(i) = z(i)
    END DO

  END SUBROUTINE copy2

  SUBROUTINE vxdiag(n,d,x,y)  ! Vector is multiplied by a diagonal matrix y:=d*x.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector.
         d           ! Diagonal matrix stored as a vector with n elements.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= d*x.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n           ! Vectors dimension.
      
! Local Scalars
    INTEGER :: i

    DO  i = 1,n
       y(i) = x(i)*d(i)
    END DO
 
  END SUBROUTINE vxdiag

  SUBROUTINE symax(n,m,iold,a,x,y)  ! Multiplication of a dense symmetric matrix A by a vector x.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector stored in a circular order.
         a           ! Dense symmetric matrix stored in the packed form: a(n*(n+1)/2).
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector y:= a*x. Vector y has the same circular order than x.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &         ! Order of matrix A.
         m, &         ! Length of vector x, m >= n, note that only n
                      ! components from vector x are used.
         iold         ! Index, which controlls the circular order of
                      ! the vector x.

! Local Scalars
    INTEGER :: i,j,k,l

    DO j=1,n
       l=j+iold-1
       IF (l > m) l=l-m
       y(l) = zero
       k=l
       DO i=j,n
          y(l) = a((i-1)*i/2+j)*x(k)+y(l)
          k=k+1
          IF (k > m) k=k-m
       END DO
    END DO

    DO j=2,n
       l=j+iold-1
       IF (l > m) l=l-m
       k=iold
       DO i=1,j-1
          IF (k > m) k=k-m
          y(l) = a((j-1)*j/2+i)*x(k)+y(l)
          k=k+1
       END DO
    END DO
      
  END SUBROUTINE symax

  SUBROUTINE cwmaxv(n,m,a,x,y)  ! Multiplication of a columnwise stored dense 
                                ! rectangular matrix A by a vector x.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x, &        ! Input vector (dimension m).
         a           ! Rectangular matrix stored columnwise in the
                     ! one-dimensional array (dimension n*m).
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         y           ! Output vector equal to s*a*x. If m = 0 y is a zero vector. 

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Number of rows of the matrix A.
         m           ! Number of columns of the matrix A.

! Local Scalars
    INTEGER :: i,j,k
      
    DO i = 1,n
       y(i) = zero
    END DO

    k = 1
    DO j = 1,m
       CALL scsum2(n,x(j),a(k:),y)
       k = k + n
    END DO

  END SUBROUTINE cwmaxv

  SUBROUTINE rwaxv2(n,m,a,b,x,y,v,w)  ! Multiplication of two rowwise stored dense rectangular  
                                      ! matrices A and B by vectors X and Y.
    USE param, ONLY : zero
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x,y, &      ! Input vectors (dimension n).
         a,b         ! Rectangular matrices stored rowwise in the
                     ! one-dimensional array (dimension n*m).
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         v,w         ! Output vectors v=a*x and w=b*y. 

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Number of columns of the matrices A and B.
         m           ! Number of rows of the matrices A and B.

! Local Scalars
    REAL(KIND=dp) :: tmp1,tmp2
    INTEGER :: i,j,k
      
    k = 0
    DO i = 1,m
       tmp1 = zero
       tmp2 = zero
       DO j = 1,n
          tmp1 = tmp1 + a(k+j)*x(j)
          tmp2 = tmp2 + b(k+j)*y(j)
       END DO
       v(i) = tmp1
       w(i) = tmp2
       k = k + n
    END DO

  END SUBROUTINE rwaxv2


  SUBROUTINE trlieq(n,m,iold,u,x,y,job,ierr)  ! Solving x from linear equation u*x=y or u'*x=y, 
                                              ! where u is an upper triangular matrix.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         y, &        ! Input vector stored in a circular order.
         u           ! Triangular matrix.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector y has the same circular order than x.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &         ! Order of matrix U.
         m, &         ! Length of vectors x and y, m >= n, note that only n
                      ! components from vectors are used.
         iold, &      ! Index, which controlls the circular order of
                      ! the vectors x and y.
         job          ! Option:
                      !   0  - x:=(u')**(-1)*y, u upper triangular.
                      !   1  - x:=u**(-1)*y, u upper triangular.
    INTEGER, INTENT(OUT) :: &
         ierr         ! Error indicador: 
                      !   0   - Everything is ok.
                      !  -3   - Error; 0 at diagonal.

! Local Scalars
    INTEGER :: i,ii,ij,j,k,l,ji

! Intrinsic Functions
    INTRINSIC ABS
      
    ierr = -3
      
    DO i=1,m
       x(i)=y(i)
    END DO
      
    IF (job == 0) THEN
     
! x=u'**(-1)*y, u' = [u1         ] is lower triangular.
!                    [u2 u3      ]
!                    [u4 u5 u6   ]
!                    [.  .  .  . ]
         
       ii = 0
       DO  i = 1,n
          ii=ii+i
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          x(l) = x(l)/u(ii)
          DO j = i+1,n
             ji = (j-1)*j/2+i
             k=j+iold-1
             IF (k > m) k=k-m
             x(k) = x(k) - u(ji)*x(l)
          END DO
       END DO
             
         
    ELSE IF (job == 1) THEN
     
! x=u**(-1)*y, u = [u1 u2 u4 . ] is upper triangular.
!                  [   u3 u5 . ]
!                  [      u6 . ]
!                  [         . ]
         
       ii = n* (n+1)/2
       DO i = n,1,-1
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          ij = ii
          DO j = i + 1,n
             k=j+iold-1
             IF (k > m) k=k-m
             ij = ij + j - 1
             x(l) = x(l) - u(ij)*x(k)
          END DO
          x(l)=x(l)/u(ii)
          ii = ii - i
       END DO
         
         
    ELSE
         
       RETURN
    END IF
      
    ierr = 0

  END SUBROUTINE trlieq

  ! VARIANT SUBROUTINE where trlieq(n,m,iold,u,x,y,job,ierr) -> trlieq2(n,m,iold,u,y,job,ierr) for inputs where x = y
  SUBROUTINE trlieq2(n,m,iold,u,y,job,ierr)    ! Solving x from linear equation u*x=y or u'*x=y, 
                                               ! where u is an upper triangular matrix.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
!   REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
!        y, &        ! Input vector stored in a circular order.
         y           ! Input vector stored in a circular order.
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         u           ! Triangular matrix.
!   REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
!        x           ! Output vector x_new:= a*x. Vector x_new has the same circular order than x. (variant: x_new replaces y)
                     ! Note that x may be equal to y in calling sequence. ! DEBUGGING: Causes compiler errors!

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &         ! Order of matrix U.
         m, &         ! Length of vectors x and y, m >= n, note that only n
                      ! components from vectors are used.
         iold, &      ! Index, which controlls the circular order of
                      ! the vectors x and y.
         job          ! Option:
                      !   0  - x:=(u')**(-1)*y, u upper triangular.
                      !   1  - x:=u**(-1)*y, u upper triangular.
    INTEGER, INTENT(OUT) :: &
         ierr         ! Error indicador: 
                      !   0   - Everything is ok.
                      !  -3   - Error; 0 at diagonal.

! Local Scalars
    INTEGER :: i,ii,ij,j,k,l,ji

! Local Array
    REAL(KIND=dp), DIMENSION(m) :: &
         x            ! DEBUGGING: For variant trlieq2 to avoid compiler errors a local x is used

! Intrinsic Functions
    INTRINSIC ABS
      
    ierr = -3
      
    DO i=1,m
       x(i)=y(i)
    END DO
      
    IF (job == 0) THEN
     
! x=u'**(-1)*y, u' = [u1         ] is lower triangular.
!                    [u2 u3      ]
!                    [u4 u5 u6   ]
!                    [.  .  .  . ]
         
       ii = 0
       DO  i = 1,n
          ii=ii+i
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          x(l) = x(l)/u(ii)
          DO j = i+1,n
             ji = (j-1)*j/2+i
             k=j+iold-1
             IF (k > m) k=k-m
             x(k) = x(k) - u(ji)*x(l)
          END DO
       END DO
             
         
    ELSE IF (job == 1) THEN
     
! x=u**(-1)*y, u = [u1 u2 u4 . ] is upper triangular.
!                  [   u3 u5 . ]
!                  [      u6 . ]
!                  [         . ]
         
       ii = n* (n+1)/2
       DO i = n,1,-1
          l=i+iold-1
          IF (l > m) l=l-m
          IF (ABS(u(ii)) <= small) RETURN
          ij = ii
          DO j = i + 1,n
             k=j+iold-1
             IF (k > m) k=k-m
             ij = ij + j - 1
             x(l) = x(l) - u(ij)*x(k)
          END DO
          x(l)=x(l)/u(ii)
          ii = ii - i
       END DO
         
         
    ELSE
         
       RETURN
    END IF
      
    ierr = 0

  ! trlieq2 variant; copying x to y so a single identical parameter can be used for calling the subroutine
    DO i=1,m
       y(i)=x(i)
    END DO

  END SUBROUTINE trlieq2
      
  SUBROUTINE lineq(n,m,iold,a,x,y,ierr)  ! Solving X from linear equation A*X=Y. 
                                         ! Positive definite matrix A+E is given using 
                                         ! the factorization A+E=L*D*L' obtained by the
                                         ! subroutine mxdpgf.
    USE param, ONLY : small
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         y, &        ! Input vector stored in a circular order (dimension m).
         a           ! Factorization a+e=l*d*l' obtained by the subroutine mxdpgf.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold        ! Index, which controlls the circular order of
                     ! the vectors x and y.
    INTEGER, INTENT(OUT) :: &
         ierr        ! Error indicador: 
                     !   0   - Everything is ok.
                     !  -2   - Error; indefinite matrix.

! Local Scalars
    INTEGER :: i,ii,ij,j,k,l


    ierr = -2
      
! Phase 1: x=l**(-1)*x

    ij = 0
    DO i = 1,n
       l=i+iold-1
       IF (l > m) l=l-m
       x(l) = y(l)
         
       DO j = 1,i - 1
          ij = ij + 1
          k=j+iold-1
          IF (k > m) k=k-m
          x(l) = x(l) - a(ij)*x(k)
       END DO
       ij = ij + 1
    END DO

! Phase 2 : x:=d**(-1)*x

    ii = 0
    DO i = 1,n
       ii = ii + i
       IF (a(ii) <= small) RETURN
       l=i+iold-1
       IF (l > m) l=l-m
       x(l) = x(l)/a(ii)
    END DO

! Phase 3 : x:=trans(l)**(-1)*x

    ii = n* (n-1)/2
    DO i = n - 1,1,-1
       ij = ii
       l=i+iold-1
       IF (l > m) l=l-m
       DO j = i + 1,n
          k=j+iold-1
          IF (k > m) k=k-m
          ij = ij + j - 1
          x(l) = x(l) - a(ij)*x(k)
       END DO
       ii = ii - i
    END DO

    ierr = 0

  END SUBROUTINE lineq

  SUBROUTINE calq(n,m,iold,a,x,y,iprint)  ! Solving x from linear equation A*x=y.
    USE param, ONLY : zero,small,one
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         y           ! Input vector stored in a circular order (dimension m).
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         a           ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         x           ! Output vector y:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold, &     ! Index, which controlls the circular order of
                     ! the vectors x and y.
         iprint      ! Printout specification.
!    INTEGER, INTENT(OUT) :: &
!         ierr        ! Error indicador: 
!                     !   0   - Everything is ok.
!                     !  -2   - Error; indefinite matrix.

! Local Scalars
    REAL(KIND=dp) :: eta,bet
    INTEGER :: inf,ierr

      
    eta = small+small
      
    CALL mxdpgf(n,a,inf,eta,bet)


    IF (iprint == 2) THEN
       IF (inf < 0) THEN
 !         WRITE (6,FMT='(1X,''Warning: Insufficiently positive'' &
 !              '' definite matrix detected. '')')
 !         WRITE (6,FMT='(1X,''Correction added.'')')
 !        
       ELSE IF (inf > 0) THEN
 !         WRITE (6,FMT='(1X,''Warning: Indefinite'' &
 !           '' matrix detected. Correction added.'')')
       END IF
    END IF
      
    CALL lineq(n,m,iold,a,x,y,ierr)
    IF (ierr /= 0) THEN
       IF (iprint == 2) THEN
!          WRITE (6,FMT='(1X,''Warning: Indefinite matrix detected. '')')
       END IF
    END IF

  CONTAINS

    SUBROUTINE mxdpgf(n,a,inf,alf,tau)  ! Factorization A+E=L*D*trans(L) of a dense symmetric positive
                                        ! definite matrix A+E, where D and E are diagonal positive 
                                        ! definite matrices and L is a lower triangular matrix. 
                                        ! If A is sufficiently positive definite then E=0.
      
! Array Arguments
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
           a         ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).

! Scalar Arguments
      REAL(KIND=dp), INTENT(INOUT) :: &
           alf       ! On input a desired tolerance for positive definiteness. 
                     ! On output the most negative diagonal element used in the factorization
                     ! process (if inf>0).
      REAL(KIND=dp), INTENT(OUT) :: &
           tau       ! Maximum diagonal element of matrix E.

      INTEGER, INTENT(IN) :: &
           n         ! Order of matrix a.
      INTEGER, INTENT(OUT) :: &
           inf       ! An information obtained in the factorization process:
                     !    inf=0  - A is sufficiently positive definite and E=0. 
                     !    inf<0  - A is not sufficiently positive definite and E>0.
                     !    inf>0  - A is indefinite and inf is an index of the most negative 
                     !             diagonal element used in the factorization process.

! Local Scalars
      REAL(KIND=dp) :: bet,del,gam,rho,sig,tol
      INTEGER :: i,ij,ik,j,k,kj,kk,l

! Intrinsic Functions
      INTRINSIC ABS,MAX

      l = 0
      inf = 0
      tol = alf
      

! Estimation of the matrix norm

      alf = zero
      bet = zero
      gam = zero
      tau = zero
      kk = 0

      DO k = 1,n
         kk = kk + k
         bet = MAX(bet,ABS(a(kk)))
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
      END DO
      bet = MAX(tol,bet,gam/n)

      del = tol*MAX(bet,one)
      kk = 0
      DO k = 1,n
         kk = kk + k

!     Determination of a diagonal correction

         sig = a(kk)
         IF (alf > sig) THEN
            alf = sig
            l = k
         END IF

         gam = zero
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
         gam = gam*gam
         rho = MAX(ABS(sig),gam/bet,del)
         IF (tau < rho-sig) THEN
            tau = rho - sig
            inf = -1
         END IF
         
! Gaussian elimination

         a(kk) = rho
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = a(kj)
            a(kj) = gam/rho
            ik = kk
            ij = kj
            DO i = k + 1,j
               ik = ik + i - 1
               ij = ij + 1
               a(ij) = a(ij) - a(ik)*gam
            END DO
         END DO
      END DO
      IF (l > 0 .AND. ABS(alf) > del) inf = l
      
    END SUBROUTINE mxdpgf
  END SUBROUTINE calq

! Variant calq2 where x = y

!  SUBROUTINE calq2(n,m,iold,a,x,y,iprint)  ! Solving x from linear equation A*x=y.
  SUBROUTINE calq2(n,m,iold,a,y,iprint)  ! Solving x from linear equation A*x=y.
    USE param, ONLY : zero,small,one
    IMPLICIT NONE

! Array Arguments
!    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
!         y           ! Input vector stored in a circular order (dimension m).
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         a           ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).
!    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         y           ! Output vector x:= a*x. Vector x has the same circular order than y.
                     ! Note that x may be equal to y in calling sequence.

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         n, &        ! Order of matrix a.
         m, &        ! Length of vectors x and y, m >= n, note that only n
                     ! components from vectors are used.
         iold, &     ! Index, which controlls the circular order of
                     ! the vectors x and y.
         iprint      ! Printout specification.
!    INTEGER, INTENT(OUT) :: &
!         ierr        ! Error indicador: 
!                     !   0   - Everything is ok.
!                     !  -2   - Error; indefinite matrix.

! Local Array
    REAL(KIND=dp), DIMENSION(m) :: & ! Variant
         x

! Local Scalars
    REAL(KIND=dp) :: eta,bet
    INTEGER :: inf,ierr,i ! i added for variant


  ! Variant; copy y to x
  DO i = 1,m
     x(i) = y(i)
  END DO
      
    eta = small+small
      
    CALL mxdpgf(n,a,inf,eta,bet)


    IF (iprint == 2) THEN
       IF (inf < 0) THEN
 !         WRITE (6,FMT='(1X,''Warning: Insufficiently positive'' &
 !              '' definite matrix detected. '')')
 !         WRITE (6,FMT='(1X,''Correction added.'')')
 !        
       ELSE IF (inf > 0) THEN
 !         WRITE (6,FMT='(1X,''Warning: Indefinite'' &
 !           '' matrix detected. Correction added.'')')
       END IF
    END IF

    CALL lineq(n,m,iold,a,x,y,ierr)
    IF (ierr /= 0) THEN
       IF (iprint == 2) THEN
!          WRITE (6,FMT='(1X,''Warning: Indefinite matrix detected. '')')
       END IF
    END IF

  ! Variant: Copy x to back to y
  DO i = 1,m
     y(i) = x(i)
  END DO


  CONTAINS

    SUBROUTINE mxdpgf(n,a,inf,alf,tau)  ! Factorization A+E=L*D*trans(L) of a dense symmetric positive
                                        ! definite matrix A+E, where D and E are diagonal positive 
                                        ! definite matrices and L is a lower triangular matrix. 
                                        ! If A is sufficiently positive definite then E=0.
      
! Array Arguments
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
           a         ! On input: Dense symmetric matrix stored in the packed form. 
                     ! On output: factorization A+E=L*D*trans(L).

! Scalar Arguments
      REAL(KIND=dp), INTENT(INOUT) :: &
           alf       ! On input a desired tolerance for positive definiteness. 
                     ! On output the most negative diagonal element used in the factorization
                     ! process (if inf>0).
      REAL(KIND=dp), INTENT(OUT) :: &
           tau       ! Maximum diagonal element of matrix E.

      INTEGER, INTENT(IN) :: &
           n         ! Order of matrix a.
      INTEGER, INTENT(OUT) :: &
           inf       ! An information obtained in the factorization process:
                     !    inf=0  - A is sufficiently positive definite and E=0. 
                     !    inf<0  - A is not sufficiently positive definite and E>0.
                     !    inf>0  - A is indefinite and inf is an index of the most negative 
                     !             diagonal element used in the factorization process.

! Local Scalars
      REAL(KIND=dp) :: bet,del,gam,rho,sig,tol
      INTEGER :: i,ij,ik,j,k,kj,kk,l

! Intrinsic Functions
      INTRINSIC ABS,MAX

      l = 0
      inf = 0
      tol = alf
      

! Estimation of the matrix norm

      alf = zero
      bet = zero
      gam = zero
      tau = zero
      kk = 0

      DO k = 1,n
         kk = kk + k
         bet = MAX(bet,ABS(a(kk)))
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
      END DO
      bet = MAX(tol,bet,gam/n)

      del = tol*MAX(bet,one)
      kk = 0
      DO k = 1,n
         kk = kk + k

!     Determination of a diagonal correction

         sig = a(kk)
         IF (alf > sig) THEN
            alf = sig
            l = k
         END IF

         gam = zero
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = MAX(gam,ABS(a(kj)))
         END DO
         gam = gam*gam
         rho = MAX(ABS(sig),gam/bet,del)
         IF (tau < rho-sig) THEN
            tau = rho - sig
            inf = -1
         END IF
         
! Gaussian elimination

         a(kk) = rho
         kj = kk
         DO j = k + 1,n
            kj = kj + j - 1
            gam = a(kj)
            a(kj) = gam/rho
            ik = kk
            ij = kj
            DO i = k + 1,j
               ik = ik + i - 1
               ij = ij + 1
               a(ij) = a(ij) - a(ik)*gam
            END DO
         END DO
      END DO
      IF (l > 0 .AND. ABS(alf) > del) inf = l
      
    END SUBROUTINE mxdpgf
  END SUBROUTINE calq2


END MODULE lmbm_sub


MODULE initializat  ! Initialization of parameters and x_var for LDGBM and LMBM

  USE constants
  USE functions                  ! Contains INFORMATION from the USER

  IMPLICIT NONE

    ! Parameters 
    INTEGER, PARAMETER :: &
    
        na     =    2, &        ! Size of the bundle na >= 2.
        mcu    =    15, &        ! Upper limit for maximum number of stored corrections, mcu >= 3.
        mcinit =    7, &        ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                      ! If mcinit <= 0, the default value mcinit = 3 will be used. 
                                      ! However, the value mcinit = 7 is recommented.
        inma    = 3, &          ! Selection of line search method:
                                      !   inma = 0, Armijo line search, Not here
                                      !   inma = 1, nonmonotone Armijo line search. Not here
                                      !   inma = 2, weak Wolfe line search.
                                      !   inma = 3, nonmonotone  weak Wolfe line search.
        mnma    = 10, &         ! Maximum number of function values used in nonmonotone line search.
        maxnin  = 20            ! Maximum number of interpolations, maxnin >= 0.
                                      ! The value maxnin = 2-20 is recommented with inma=0,
                                      ! maxnin >= 20 with inma=1 and 3, and maxnin =200 with inma=2.
                                      ! For example:
                                      !   inma = 0, maxin = 20.
                                      !   inma = 1, mnma = 20, maxin = 30.
                                      !   inma = 2, maxnin = 200.
                                      !   inma = 3, mnma=10, maxnin = 20.

    REAL(KIND=dp), PARAMETER :: &
        time =  1000.0_dp        !  Maximum CPU-time in seconds. If time <= 0.0 the maximum time 
                                       !  is ignored. REAL argument.

    ! Real parameters (if parameter value <= 0.0 the default value of the parameter will be used).
    REAL(KIND=dp), SAVE :: &
        tolf = 0.0_dp, &                          ! Tolerance for change of function values (default = 1.0E-8).
        tolf2 = -10.0_dp, &                         ! Second tolerance for change of function values.
                                         !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                         !   criterion will be ignored. 
                                         !   - If tolf2 = 0 the default value 1.0E+4 will be used. 
        tolb  =  0.0_dp , &             ! Tolerance for the function value (default = -large).
        
        tolg  = 1.0E-5_dp, &      ! Tolerance for the first termination criterion (default = 1.0E-6).
        tolg2 = 1.0E-3_dp, &      ! Tolerance for the second termination criterion (default = tolg). clustering code small data
    
        eta   =    1.0E-4_dp, &   ! Distance measure parameter, eta >= 0. 
                                         !   - If eta < 0  the default value 0.5 will be used. 
          
        epsl  =    0.24E+00_dp, &  ! Line search parameter, 0 < epsl < 0.25 (default = 1.0E-4.

        xmax  =    1000.0_dp     ! Maximum stepsize, 1 < XMAX (default = 1.5).

    INTEGER, SAVE :: n        ! The dimension

    INTEGER, SAVE :: nproblem ! The solved problem 


    ! Integer parameters (if value <= 0 the default value of the parameter will be used).
    INTEGER, SAVE :: &
    
        mittt   = 500, &          ! Maximun number of iterations (default = 10000).
        mfe     = 500, &          ! Maximun number of function evaluations (default = n*mittt).
        mtesf   =     0, &        ! Maximum number of iterations with changes of
                                  ! function values smaller than tolf (default = 10).
        iiprint =      0, &       ! Printout specification:
                                  !     0  - Only the error messages.
                                  !     1  - The final values of the objective function 
                                  !          (default used if iiprint < -1).
                                  !     2  - The final values of the objective function and the 
                                  !          most serious warning messages.
                                  !     3  - The whole final solution. 
                                  !     4  - At each iteration values of the objective function.
                                  !     5  - At each iteration the whole solution
        method  =      0, &       ! Selection of the method:
                                  !     0  - Limited memory bundle method (default).
                                  !     1  - L-BFGS bundle method.
                                  !     2  - Limited memory discrete gradient bundle method.
        iscale  =      0          ! Selection of the scaling:
                                  !     0  - Scaling at every iteration with STU/UTU (default).
                                  !     1  - Scaling at every iteration with STS/STU.
                                  !     2  - Interval scaling with STU/UTU.
                                  !     3  - Interval scaling with STS/STU.
                                  !     4  - Preliminary scaling with STU/UTU.
                                  !     5  - Preliminary scaling with STS/STU.
                                  !     6  - No scaling.      


    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, SAVE :: x_var   ! Vector of variables (initialize x in subroutine init_x)          

CONTAINS


    SUBROUTINE init_lmbmpar(kl)  ! User supplied subroutine for further initialization of parameters
                               ! (when needed) for LMBM. May be left empty.

        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: kl

        ! For big data sets use larger values
        tolg =1.0E+00
        tolg2 = 1.0E+3
        IF (nproblem == 1) THEN
            tolg = 1.0E+3
            IF(kl == 5) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 10) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 15) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 20) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 25) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 30) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 35) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 40) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 45) THEN
                tolg = 1.0E-3
            ELSE IF(kl == 50) THEN
                tolg = 1.0E-3
            END IF
        END IF

    ! For small data sets
    !        tolg2 = 1.0E-1
    !        IF(ns == 1) tolg = 1.0E-5
    !        IF(ns == 2) tolg = 1.0E-4

    END SUBROUTINE init_lmbmpar

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!-------------------------------------------------------------------------------------------
    SUBROUTINE allocate_x_var(user_N) 
        
        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        INTEGER, INTENT(IN) :: user_N                 ! the of features/inputs in a data point   
        !************************************************************************************
        n=user_N
   
        ALLOCATE(x_var(n))
  
    END SUBROUTINE allocate_x_var
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
!-------------------------------------------------------------------------------------------  
    SUBROUTINE deallocate_x_var() 
         
        IMPLICIT NONE
  
        DEALLOCATE(x_var)
  
    END SUBROUTINE deallocate_x_var  
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
  
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
!-------------------------------------------------------------------------------------------  
    SUBROUTINE init_x_var(x_help)  ! User supplied subroutine for initialization of vector x(n) and the solved problem

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x_help 
        !***********************************************************************************
        INTEGER :: i
    
        DO i = 1, n
            x_var(i) = x_help(i)
        END DO   
    
    END SUBROUTINE init_x_var
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
   
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    
!-------------------------------------------------------------------------------------------  
    SUBROUTINE init_problem(problem)

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        INTEGER, INTENT(IN) :: problem       ! 1 = whole data a, 2 = subdata b     
        !***********************************************************************************
    
        nproblem = problem
    
    END SUBROUTINE init_problem
!-------------------------------------------------------------------------------------------
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  

!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  
!------------------------------------------------------------------------------------------- 
    SUBROUTINE copy_x_var(x_help)  ! User supplied subroutine to copy the vector x(n).

        IMPLICIT NONE
        !**************************** NEEDED FROM USER *************************************
        REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: x_help 
        !***********************************************************************************
        INTEGER :: i
    
        DO i = 1, n
            x_help(i) = x_var(i)
        END DO   
    
    END SUBROUTINE copy_x_var  
!-------------------------------------------------------------------------------------------    
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/  

END MODULE initializat


MODULE obj_fun  ! Computation of the value and the subgradient of the objective function

  USE constants
  IMPLICIT NONE

  PUBLIC :: &
       myf, &    ! Computation of the value of the objective. User spesified function. 
       myg, &    ! Computation of the subgradient of the objective. User spesified function. 
       myf2, &
       myg2

CONTAINS
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE myf *                                               *
!*                                                                      *
!*     Computation of the value of the objective. User spesified        *
!*     function.                                                        *
!*                                                                      *
!************************************************************************
     
  SUBROUTINE myf(n,x,f,iterm,nproblem)

    USE functions       ! Computation of the value for problem next.                                

    IMPLICIT NONE

! Scalar Arguments
    REAL(KIND=dp), INTENT(OUT) :: f    ! Value of the objective function. 
    INTEGER, INTENT(IN) :: n           ! Number of variables.
    INTEGER, INTENT(IN) :: nproblem    ! The number of the problem solved   
    INTEGER, INTENT(OUT) :: iterm      ! Cause of termination:
                                                   !   0  - Everything is ok.
                                                   !  -3  - Failure in function calculations 
                                                   !        (assigned by the user).

! Array Arguments
    REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: &
         x  ! Vector of variables.

    iterm = 0
    
! Function evaluation (give your function here).
    IF (nproblem == 1) THEN             ! data a
        IF (better_indices2 == 1) THEN
            CALL func_a_index(x, n, f)
        ELSE
            CALL func_a(x, n, f)
        END IF
    ELSE IF (nproblem == 2) THEN        ! data b
        IF (better_indices == 1) THEN
            CALL func_b_index(x, n, f)
        ELSE
            CALL func_b(x, n, f)
        END IF
    ELSE IF (nproblem == 3) THEN        ! auxiliary problem, data a
        CALL func_help_a(x, n, f)
    ELSE
        CALL func_help_b(x, n, f)       ! auxiliary problem, data b
    END IF 
      
  END SUBROUTINE myf
  
  
    SUBROUTINE myf2(n,x,f,iterm,nproblem)

    USE functions       ! Computation of the value for problem next.                             

    IMPLICIT NONE

! Scalar Arguments
    REAL(KIND=dp), INTENT(OUT) :: f    ! Value of the objective function. 
    INTEGER, INTENT(IN) :: n           ! Number of variables.
    INTEGER, INTENT(IN) :: nproblem    ! The number of the problem solved   
    INTEGER, INTENT(OUT) :: iterm      ! Cause of termination:
                                                   !   0  - Everything is ok.
                                                   !  -3  - Failure in function calculations 
                                                   !        (assigned by the user).

! Array Arguments
    REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: &
         x  ! Vector of variables.


    iterm = 0

! Function evaluation (give your function here).
    IF (nproblem == 1) THEN             ! data a
        IF (better_indices2 == 1) THEN
            CALL func_a_clust_index(x, n, f)
        ELSE
            CALL func_a_clust(x, n, f)
        END IF    
    ELSE IF (nproblem == 2) THEN        ! data b
        IF (better_indices == 1) THEN
            CALL func_b_clust_index(x, n, f)
        ELSE
            CALL func_b_clust(x, n, f)
        END IF
    ELSE IF (nproblem == 3) THEN        ! auxiliary problem, data a
        CALL func_help_a(x, n, f)
    ELSE
        CALL func_help_b(x, n, f)       ! auxiliary problem, data b
    END IF 
      
  END SUBROUTINE myf2

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE MYG *                                               *
!*                                                                      *
!*     Computation of the subgradient of the objective function. User   *
!*     spesified function. If discrete gradients are used, this         *
!*     subroutine may be left empty.                                    *
!*                                                                      *
!************************************************************************
     
  SUBROUTINE myg(n,x,g,iterm,nproblem)

    USE functions       ! Computation of the subgradient for problem next.

    IMPLICIT NONE

! Scalar Arguments
    INTEGER, INTENT(IN) :: n            ! Number of variables.
    INTEGER, INTENT(IN) :: nproblem     ! The number of the problem solved
    INTEGER, INTENT(OUT) :: iterm       ! Cause of termination:
                                                    !   0  - Everything is ok.
                                                    !  -3  - Failure in subgradient calculations 
                                                    !        (assigned by the user).
! Array Arguments
    REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
    REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

    iterm = 0

! Subgradient evaluation (give your function here).
    IF (nproblem == 1) THEN             ! data a
        IF (better_indices2 == 1) THEN
            CALL subgrad_a_index(x, n, g)
        ELSE
            CALL subgrad_a(x, n, g)
        END IF    
    ELSE IF (nproblem == 2) THEN        ! data b
        IF (better_indices == 1) THEN
            CALL subgrad_b_index(x, n, g)
        ELSE
            CALL subgrad_b(x, n, g)
        END IF
    ELSE IF (nproblem == 3) THEN        ! auxiliary problem, data a
        CALL subgrad_help_a(x, n, g)
    ELSE
        CALL subgrad_help_b(x, n, g)    ! auxiliary problem, data b
    END IF
      
  END SUBROUTINE myg
  
  
    SUBROUTINE myg2(n,x,g,iterm,nproblem)

    USE functions       ! Computation of the subgradient for problem next.

    IMPLICIT NONE

! Scalar Arguments
    INTEGER, INTENT(IN) :: n            ! Number of variables.
    INTEGER, INTENT(IN) :: nproblem     ! The number of the problem solved
    INTEGER, INTENT(OUT) :: iterm       ! Cause of termination:
                                                    !   0  - Everything is ok.
                                                    !  -3  - Failure in subgradient calculations 
                                                    !        (assigned by the user).
! Array Arguments
    REAL(KIND=dp), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
    REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

    iterm = 0

! Subgradient evaluation (give your function here).
    IF (nproblem == 1) THEN             ! data a
        IF (better_indices2 == 1) THEN
            CALL subgrad_a_clust_index(x, n, g)
        ELSE
            CALL subgrad_a_clust(x, n, g)
        END IF     
    ELSE IF (nproblem == 2) THEN        ! data b
        IF (better_indices == 1) THEN
            CALL subgrad_b_clust_index(x, n, g)
        ELSE
            CALL subgrad_b_clust(x, n, g)
        END IF
    ELSE IF (nproblem == 3) THEN        ! auxiliary problem, data a
        CALL subgrad_help_a(x, n, g)
    ELSE
        CALL subgrad_help_b(x, n, g)    ! auxiliary problem, data b
    END IF
      
  END SUBROUTINE myg2
  

END MODULE obj_fun


MODULE lmbm_mod  ! Limited memory bundle method

  USE constants
  IMPLICIT NONE

! MODULE lmbm_mod includes the following subroutines (S) and functions (F).
  PUBLIC :: &
       init_lmbm, &  ! S Initialization for limited memory bundle subroutine.
       lmbm          ! S Limited memory bundle subroutine for nonsmooth 
                     !   large-scale optimization. Contains:
                     !     S restar  Initialization and reinitialization.
                     !     S dobun   Bundle construction.
  PRIVATE :: &
       indic1, &     ! S Initialization of indices.
       tinit, &      ! S Calculation of initial step size. Contains:
                     !     S destep  Stepsize selection using polyhedral 
                     !                 approximation for descent steps.
                     !     S nulstep Stepsize selection using polyhedral 
                     !                 approximation for null steps.
       lls, &        ! S Line search. Contains:
                     !     F qint    Quadratic interpolation.
       nmlls, &      ! S Nonmonotonic Weak Wolfe line search. Contains:
                      !     F qint    Quadratic interpolation.
       dlbfgs, &     ! S Computing the search direction by limited memory 
                     !   BFGS update. Contains:
                     !     F sclpar  Selection of scaling parameter.
       dlskip, &     ! S Skipping the updates and computing the search
                     !   direction by limited memory BFGS update.    
       dlsr1, &      ! S Computing the search direction by limited
                     !   memory SR1 update.    
       agbfgs, &     ! S Simplified subgradient aggregation.
       aggsr1, &     ! S Subgradient aggregation.
       agskip        ! S Subgradient aggregation using BFGS update.
       !, &     ! S Subgradient aggregation using BFGS update.
!       wprint, &     ! S Printout the error and warning messages.
!       rprint        ! S Printout the results.

CONTAINS
!***********************************************************************
!*                                                                     *
!*     * SUBROUTINE init_lmbm *                                        *
!*                                                                     *
!*     Initialization for limited memory bundle subroutine for         *
!*     large-scale unconstrained nonsmooth optimization.               *
!*                                                                     *
!***********************************************************************
      
  SUBROUTINE init_lmbm(mc,iterm)
    USE param, ONLY : zero,half,small,large
    USE initializat, ONLY : &
         n, &         ! Number of variables.
         na, &        ! Maximum bundle dimension, na >= 2.
         mcu, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
         iiprint, &    ! Printout specification (see initializat for more details).
         !method, &    ! Selection of the method (see initializat for more details).
         iscale, &    ! Selection of the scaling (see initializat for more details).
         tolf, &      ! Tolerance for change of function values.
         tolf2, &     ! Second tolerance for change of function values.
         tolb, &      ! Tolerance for the function value.
         tolg, &      ! Tolerance for the first termination criterion. 
         tolg2, &     ! Tolerance for the second termination criterion.
         xmax, &      ! Maximum stepsize, 1 < XMAX.
         eta, &       ! Distance measure parameter, eta >= 0. 
         epsl, &      ! Line search parameter, 0 < epsl < 0.25.
         mtesf, &     ! Maximum number of iterations with changes of
                      ! function values smaller than tolf.
         mittt, &       ! Maximun number of iterations.
         mfe          ! Maximun number of function evaluations.

    IMPLICIT NONE

! Scalar Arguments
    INTEGER, INTENT(INOUT) :: mc     ! Initial maximum number of stored corrections.
    INTEGER, INTENT(OUT) :: iterm    ! Cause of termination:
                                                 !   0  - Everything is ok.
                                                 !  -5  - Invalid input parameters.
 
! Initialization and error checking

    IF (iiprint < -1) iiprint  = 1     ! Printout specification. 
    iterm = 0
      
    !IF (n <= 0) THEN
       !iterm = -5
!       CALL wprint(iterm,iiprint,1)
       !RETURN 
    !END IF

    IF (na < 2) THEN
       iterm = -5
!       CALL wprint(iterm,iiprint,3)
       RETURN
    END IF

    IF (epsl >= 0.25_dp) THEN
       iterm = -5
!       CALL wprint(iterm,iiprint,4)
       RETURN
    END IF

    IF (mcu <= 3) THEN
       iterm = -5
!       CALL wprint(iterm,iiprint,2)
       RETURN
    END IF

    IF (mc > mcu) THEN
       mc = mcu
!       CALL wprint(iterm,iiprint,-1)
    END IF

! Default values

    IF (mc    <= 0) mc       = 3               ! Initial maximum number of corrections.
    IF (mittt   <= 0) mittt      = 500  !oli ennen 10000           ! Maximum number of iterations.
    IF (mfe   <= 0) mfe      = n*mittt           ! Maximum number of function evaluations.
    IF (tolf  <= zero) tolf  = 1.0E-08_dp    ! Tolerance for change of function values.
    IF (tolf2 == zero) tolf2 = 1.0E+04_dp    ! Second tolerance for change of function values.
    IF (tolb  == zero) tolb  = -large + small  ! Tolerance for the function value.
    IF (tolg  <= zero) tolg  = 1.0E-05_dp    ! Tolerance for the first termination criterion.
    IF (tolg2 <= zero) tolg = 1.0E-03_dp            ! Tolerance for the second termination criterion.
    IF (xmax  <= zero) xmax  = 1.5_dp        ! Maximum stepsize.
    IF (eta   <  zero) eta   = half            ! Distance measure parameter
    IF (epsl  <= zero) epsl  = 1.0E-04_dp    ! Line search parameter,
    IF (mtesf <= 0) mtesf    = 10              ! Maximum number of iterations with changes 
                                               ! of function values smaller than tolf.
    !IF (method > 1 .OR. method < 0) method = 0 ! Selection of the method.
    IF (iscale > 6 .OR. iscale < 0) iscale = 0 ! Selection of the scaling.

     
  END SUBROUTINE init_lmbm

!***********************************************************************
!*                                                                     *
!*     * SUBROUTINE lmbm *                                             *
!*                                                                     *
!*     Limited memory bundle subroutine for nonsmooth optimization.    *
!*                                                                     *
!***********************************************************************
      
  SUBROUTINE lmbm(mc,f,nit,nfe,nge,iterm,strtim)
    USE param, ONLY : small,large,zero,half,one
    USE exe_time, ONLY : getime  ! Execution time.
    USE initializat, ONLY : &
         n, &         ! Number of variables.
         x_var, &         ! Vector of variables
         nproblem, &  ! The solved problem
         na, &        ! Maximum bundle dimension, na >= 2.
         mcu, &       ! Upper limit for maximum number of stored corrections, mcu >= 3.
         inma, &      ! Selection of line search method.
         mnma, &      ! Maximum number of function values used in nonmonotone line search.
         iiprint, &    ! Printout specification (see initializat for more details).
         method, &    ! Selection of the method (see initializat for more details).
         iscale, &    ! Selection of the scaling (see initializat for more details).
         tolf, &      ! Tolerance for change of function values.
         tolf2, &     ! Second tolerance for change of function values.
         tolb, &      ! Tolerance for the function value.
         tolg, &      ! Tolerance for the first termination criterion. 
         tolg2, &     ! Tolerance for the second termination criterion.
         xmax, &      ! Maximum stepsize, 1 < XMAX.
         eta, &       ! Distance measure parameter, eta >= 0. 
         epsl, &      ! Line search parameter, 0 < epsl < 0.25.
         mtesf, &     ! Maximum number of iterations with changes of
                      ! function values smaller than tolf.
         mittt, &     ! Maximun number of iterations.
         mfe, &       ! Maximun number of function evaluations.
         time!, &         ! Maximum time

    USE obj_fun, ONLY : &
         myf, &       ! Computation of the value of the objective function. 
         myg, &          ! Computation of the subgradient of the objective function. 
         myf2, &
         myg2
    USE lmbm_sub, ONLY : &
         vdot, &      ! Dot product.
         vneg, &      ! Copying of a vector with change of the sign.
         xdiffy, &    ! Difference of two vectors.
         xdiffy2, &   ! Difference of two vectors. (Variant)
         copy, &      ! Copying of a vector.
         copy2        ! Copying of two vectors.


    IMPLICIT NONE

! Scalar Arguments
    REAL(KIND=dp), INTENT(OUT) :: &
         f            ! Value of the objective function.
    INTEGER, INTENT(OUT) :: & 
         nit,nfe,nge  ! Number of iterations, and function and subgradient evaluations.
    INTEGER, INTENT(INOUT) :: &
         mc           ! Maximum number of stored corrections.
    INTEGER, INTENT(OUT) :: &
         iterm        ! Cause of termination:
                      !   1  - The problem has been solved with desired accuracy.
                      !   2  - Changes in function values < tolf in mtesf
                      !        subsequent iterations.
                      !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                      !        where small is the smallest positive number such that 
                      !        1.0 + small > 1.0.
                      !   4  - Number of function calls > mfe.
                      !   5  - Number of iterations > mittt.
                      !   6  - Time limit exceeded. 
                      !   7  - f < tolb.
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.
                      !  -6  - Unspecified error.


! Local Arrays
    REAL(KIND=dp), DIMENSION(n) :: &
         xo, &        ! Previous vector of variables.
         g, &         ! Subgradient of the objective function.
         gp, &        ! Previous subgradient of the ohjective function.
         ga, &        ! Aggregate subgradient.
         s, &         ! Difference of current and previous variables.
         u, &         ! Difference of current and previous subgradients.
         d, &         ! Direction vector.
         tmpn1        ! Auxiliary array.
    REAL(KIND=dp), DIMENSION(n*(mcu+1)) :: &
         sm,um        ! Matrises whose columns are stored differences of
                      ! variables (sm) and subgradients (um).
    REAL(KIND=dp), DIMENSION((mcu+2)*(mcu+1)/2) :: &
         rm, &        ! pper triangular matrix stored columnwise in the 
                      ! one-dimensional array.
         umtum        ! Matrix whose columns are stored subgradient differences.
    REAL(KIND=dp), DIMENSION(mcu+1) :: &
         cm, &        ! Diagonal matrix.
         smtgp, &     ! smtgp = trans(sm)*gp.
         umtgp, &     ! umtgp = trans(um)*gp.
         tmpmc1, &    ! Auxiliary array.
         tmpmc2       ! Auxiliary array.
    REAL(KIND=dp), DIMENSION(n*na) :: &
         ax,ag        ! Matrix whose columns are bundle points and subgradients.
    REAL(KIND=dp), DIMENSION(na) :: &
         af           ! Vector of bundle values.
    REAL(KIND=dp), DIMENSION(mnma) :: &
         fold         ! Old function values.
         
! Local Scalars
    REAL(KIND=dp) :: &
         alfn, &      ! Locality measure.
         alfv, &      ! Aggregate locality measure.
         epsr, &      ! Line search parameter.
         dnorm, &     ! Euclidean norm of the direction vector.
         gnorm, &     ! Euclidean norm of the aggregate subgradient.
         xnorm, &     ! Stopping criterion.
         pxnorm, &    ! Previous stopping criterion.
         p, &         ! Directional derivative.
         tmax, &      ! Maximum stepsize.
         t, &         ! Stepsize.
         theta, &     ! Correction parameter for stepsize.
         fo, &        ! Previous value of the objective.
         gamma!, &     ! Scaling parameter.
!         rmse_best    ! The best RMSE for the validation set.
    INTEGER :: i, &
         mcinit, &    ! Initial maximum number of stored corrections.
         mcc, &       ! Current number of stored corrections.
         inew, &      ! Index for the circular arrays.
         inewnma, &   ! Index for the circular arrays containing function values.
         ibfgs, &     ! Index of the type of BFGS update.
         isr1, &      ! Index of the type of SR1 update.
         iters, &     ! Null step indicator.
                      !   0  - Null step.
                      !   1  - Serious step.
         nnk, &       ! Consecutive null steps counter.
         ibun, &      ! Index for the circular arrays in bundle updating.
         mal, &       ! Current size of the bundle.
         ic, &        ! Correction indicator.
         icn, &       ! Correction indicator for null steps.
         iflag, &     ! Index for adaptive version.
         neps, &      ! Number of consecutive equal stopping criterions.
         ntesf, &     ! Number of tests on function decrease.
         ncres, &     ! Number of restarts.
         nres, &      ! Number of consecutive restarts.
         nress, &     ! Number of consecutive restarts in case of tmax < tmin.
         nout         ! Auxilary printout specification.
! INTEGER :: ww !debug
! Intrinsic Functions
    INTRINSIC ABS,MAX,SQRT

! Computational Time
    REAL(KIND=dp) strtim,elapstim
       
! Parameters
    REAL(KIND=dp), PARAMETER :: &
         fmin    = -large, &        ! Smallest acceptable value of the function.
         tmin    = 1.0E-12_dp, &  ! Minimum stepsize.
         lengthd = 1.0E+20_dp, &  ! Direction vector length.
         rho     = 1.0E-12_dp     ! Correction parameter.
    INTEGER, PARAMETER :: &
         maxeps = 20, &             ! Maximum number of consecutive equal stopping criterions.
         maxnrs = 2000              ! Maximum number of restarts.
      
      
!    IF (iiprint > 3) THEN
!       IF (method == 0) WRITE (6,FMT='(1X,''Entry to LMBM:'')')
!       IF (method == 1) WRITE (6,FMT='(1X,''Entry to LBB:'')')
!    END IF
      
     
! Initialization

    inewnma = 1
    nout   = 0
    nit    = 0
    nfe    = 0
    nge    = 0
    ntesf  = 0
    nres   = 1
    ncres  = -1
    nress  = 0
    neps   = 0
    iterm  = 0
    iters  = 1
    nnk    = 0
    isr1   = 0
    alfn   = zero
    alfv   = zero
    mcinit = mc
      
    tmax   = xmax
    xnorm  = large
    fold   = -large
    

    !epsr   = 0.25_dp+small
    epsr   = 1.0E-05_dp+small
    IF (epsl+epsl >= epsr) THEN
       epsr = epsl+epsl + small
       IF (epsr >= half) THEN
          CALL wprint(iterm,iiprint,-2)
       END IF
    END IF
        
  
! Computation of the value and the subgradient of the objective
! function and the search direction for the first iteration

    CALL myf2(n,x_var,f,iterm,nproblem)
    CALL myg2(n,x_var,g,iterm,nproblem)

    nfe = nfe + 1
    nge = nge + 1

    IF (iterm /= 0) GO TO 900
    
    CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
         alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
      
    CALL dobun(n,na,mal,x_var,g,f,ax,ag,af,iters,ibun)

     
! Start of the iteration
            
    iteration: DO
    
! Computational time

       IF (time > 0.0E+00) THEN
          CALL getime(elapstim)
          IF (elapstim-strtim > time) THEN
             iterm = 6
             EXIT iteration
          END IF
       END IF


! Computation of norms

       IF (iters > 0) THEN
          gnorm = vdot(n,g,g)
          dnorm = SQRT(vdot(n,d,d))

          p = vdot(n,g,d)
          
       ELSE
          gnorm = vdot(n,ga,ga)
          dnorm = SQRT(vdot(n,d,d))

          p = vdot(n,ga,d)
       END IF
       

   
! Test on descent direction

       IF (p+small*SQRT(gnorm)*dnorm <= zero) THEN
          nres = 0
          
       ELSE
          nres = nres + 1
          IF (nres == 1) THEN
             CALL wprint(iterm,iiprint,-3)
             
             CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
                  alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
             IF (ncres > maxnrs) THEN
                nout = maxnrs
                iterm = -2
                EXIT iteration
             END IF
             
             CALL dobun(n,na,mal,x_var,g,f,ax,ag,af,iters,ibun)
             
             CYCLE iteration
          END IF
          nout = -1
          iterm = -1
          EXIT iteration
       END IF
       
 
 
! Stopping criterion

       nit = nit + 1
       pxnorm = xnorm
       xnorm = -p + 2.0_dp*alfv

     
! Tests for termination

       IF (xnorm <= 1.0E+03_dp*tolg .AND. &
            (mcc > 0 .OR. ibfgs == 2)) THEN
          
          IF(half*gnorm + alfv <= tolg2 .AND. &
               xnorm <= tolg) THEN
          
             iterm = 1
             EXIT iteration
          END IF
       
          IF (mc < mcu .AND. iflag == 0) THEN
             mc = mc+1
             iflag = 1
          END IF
       END IF



       IF (nfe >= mfe) THEN
          nout = mfe
          iterm = 4
          EXIT iteration
       END IF

      
       IF (nit >= mittt) THEN
          nout = mittt
          iterm = 5
          EXIT iteration
       END IF

      
       IF (f <= tolb) THEN
          iterm = 7
          EXIT iteration
       END IF
    
      
       IF (iters == 0) THEN
          IF (ABS(xnorm - pxnorm) <= small) THEN
             neps = neps + 1
          
             IF (neps > maxeps) THEN
                iterm = -4
                EXIT iteration
             END IF

          ELSE
             neps = 0
          END IF

       ELSE
          neps = 0
       END IF


! Correction

       IF (-p < rho*gnorm .OR. icn == 1) THEN

          xnorm = xnorm + rho*gnorm
          dnorm = SQRT(dnorm*dnorm-2.0_dp*rho*p+rho*rho*gnorm)
         
          IF (iters > 0) THEN
             DO i=1,n
                d(i) = d(i)-rho*g(i)
             END DO

          ELSE
             DO i=1,n
                d(i) = d(i)-rho*ga(i)
             END DO
             icn = 1
          END IF

          ic = 1
            
       ELSE
          ic = 0
       END IF


       IF (pxnorm < xnorm .AND. nnk > 2) THEN
          CALL wprint(iterm,iiprint,-4)
       END IF
      

       CALL rprint(n,nit,nfe,nge,x_var,f,xnorm,half*gnorm+alfv,iterm,iiprint)
      
     
!     Preparation of line search

       fo = f
       
       IF (iters > 0) THEN
          CALL copy2(n,x_var,xo,g,gp)
       END IF

       IF (dnorm > zero) tmax = xmax/dnorm
       
       IF (tmax > tmin) THEN
          nress = 0

       ELSE
          nress = nress + 1
          IF (nress == 1) THEN
             CALL wprint(iterm,iiprint,-5)
             
             CALL restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
                  alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)
             
             IF (ncres > maxnrs) THEN
                nout = maxnrs
                iterm = -2
                EXIT iteration
             END IF
             
             CALL dobun(n,na,mal,x_var,g,f,ax,ag,af,iters,ibun)
             CYCLE iteration
          END IF
          iterm = -1
          EXIT iteration
       END IF
 
 

! Initial step size

       CALL tinit(n,na,mal,x_var,af,ag,ax,ibun,d,f,p,t,tmax,tmin, &
            eta,iters,iterm)

       IF (iterm /= 0) EXIT iteration
     

! Line search with directional derivatives which allows null steps

       theta = one
       IF (dnorm > lengthd) THEN
           theta=lengthd/dnorm
       END IF
       

       IF (inma == 2) THEN! Line search with directional derivatives which allows null steps
                 ! With this the global convergence can be guaranteed.
      
            CALL lls(n,x_var,g,d,xo,t,fo,f,p,alfn,tmin,dnorm,xnorm,theta, &
                epsl,epsr,eta,iters,nfe,nge,nnk,iterm)
                
        ELSE ! Nonmonotone line search with directional derivatives which allows null steps
                 ! With this the global convergence can be guaranteed.

            CALL nmlls(n,x_var,g,d,xo,t,fo,f,fold,p,alfn,tmin,dnorm,xnorm,theta, &
                epsl,epsr,eta,iters,inewnma,nfe,nge,nnk,iterm)

            END IF
          
       IF (iterm /= 0) EXIT iteration
 
       IF (tolf2 >= 0) THEN
          IF (ABS(fo-f) <= tolf2*small*MAX(ABS(f),ABS(fo),one) &
               .AND. iters == 1) THEN
             
             iterm = 3
             EXIT iteration
          END IF
       END IF

       IF (ABS(fo-f) <= tolf) THEN
          ntesf = ntesf + 1
          
          if (ntesf >= mtesf .AND. iters == 1) THEN
             iterm = 2
             EXIT iteration
          END IF
          
       ELSE
          ntesf = 0
       END IF
      

! Bundle updating

       CALL dobun(n,na,mal,x_var,g,f,ax,ag,af,iters,ibun)
  

! Computation of variables difference 

       CALL xdiffy(n,x_var,xo,s)
  

! Computation of aggregate values and gradients difference

       IF (iters == 0) THEN
          nnk = nnk + 1

          IF (nnk == 1) THEN
             CALL copy(n,gp,tmpn1)
             CALL xdiffy(n,g,gp,u)
             
             CALL agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um, &
                  rm,cm,umtum,alfn,alfv,gamma,ic,rho)
             
          ELSE
             IF (method == 0) THEN
                CALL copy(n,ga,tmpn1)
                CALL aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv &
                     ,umtum,rm,gamma,smtgp,umtgp,tmpmc1,tmpmc2,sm &
                     ,um,icn,rho)
                CALL xdiffy(n,g,gp,u)
             ELSE
                CALL copy(n,ga,tmpn1)
                CALL xdiffy(n,g,gp,u)
                CALL agskip(n,mc,mcc,inew,iflag,g,gp,ga,d,u,alfn,alfv &
                     ,umtum,rm,cm,gamma,smtgp,umtgp,tmpmc1,tmpmc2,sm &
                     ,um,icn,rho)
             END IF
          END IF
          
          CALL copy(n,xo,x_var)
          f = fo
          
       ELSE
          IF (nnk /= 0) THEN
             CALL copy(n,ga,tmpn1)
          ELSE
             CALL copy(n,gp,tmpn1)
          END IF
          nnk = 0
          CALL xdiffy(n,g,gp,u)
       END IF

     
! Serious step initialization

       IF (iters > 0) THEN
          icn = 0
          alfn = zero
          alfv = zero
       END IF

      
! Direction finding
    
       IF (iters > 0) THEN
         
     
! BFGS update and direction determination

          CALL dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm, &
               umtum,cm,smtgp,umtgp,gamma,tmpn1,method,iscale)

       ELSE
          IF (method == 0) THEN
            
     
! SR1 update and direction determination
             CALL dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm, &
                  umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk,iiprint)
             ibfgs=0
          ELSE
             
     
! BFGS skipping and direction determination

             CALL dlskip(n,mc,mcc,inew,ibfgs,iflag,d,ga,sm,um,rm, &
                  umtum,cm,tmpmc1,tmpmc2,gamma,iscale)
          END IF
       END IF
    
    END DO iteration
      

900 CONTINUE
        

! Printout the final results

    IF (iiprint > 3) THEN
       IF (method == 0) WRITE (6,FMT='(1X,''Exit from LMBM:'')')
!       IF (method == 1) WRITE (6,FMT='(1X,''Exit from LBB:'')')
    END IF
    CALL wprint(iterm,iiprint,nout)
    CALL rprint(n,nit,nfe,nge,x_var,f,xnorm,half*gnorm+alfv,iterm,iiprint)

  CONTAINS

    SUBROUTINE restar(n,mc,mcc,mcinit,inew,ibun,ibfgs,iters,gp,g,nnk, &
         alfv,alfn,gamma,d,ic,icn,mal,ncres,iflag)  ! Initialization
      
      USE param, ONLY : zero,one  ! given in host
      USE lmbm_sub, ONLY : vneg,copy  ! given in host
      IMPLICIT NONE

! Array Arguments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
           gp        ! Basic subgradient of the objective function.
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
           g         ! Current (auxiliary) subgradient of the objective function.
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
           d         ! Search direction.

! Scalar Arguments
      INTEGER, INTENT(IN) :: & 
           n, &      ! Number of variables
           mcinit    ! Initial maximum number of stored corrections.
      INTEGER, INTENT(OUT) :: &
           mc, &     ! Current maximum number of stored corrections.
           mcc, &    ! Current number of stored corrections.
           inew, &   ! Index for the circular arrays.
           ibun, &   ! Index for the circular arrays in bundle updating.
           ibfgs, &  ! Index of the type of BFGS update.
           nnk, &    ! Consecutive null steps counter.
           ic, &     ! Correction indicator.
           icn, &    ! Correction indicator for null steps.
           mal, &    ! Current size of the bundle.
           iflag     ! Index for adaptive version.
      INTEGER, INTENT(INOUT) :: &
           iters, &  ! Null step indicator.
                     !   0  - Null step.
                     !   1  - Serious step.
           ncres     ! Number of restarts.
      REAL(KIND=dp), INTENT(OUT) :: & 
           alfn, &   ! Locality measure.
           alfv, &   ! Aggregate locality measure.
           gamma     ! Scaling parameter.


! Restart
      mc    = mcinit
      mcc   = 0
      inew  = 1
      ibun  = 1
      ibfgs = 0
      ic    = 0
      icn   = 0
      mal   = 0
      ncres = ncres + 1
      iflag = 0
        
      IF (iters == 0) THEN
         CALL copy(n,gp,g)
         iters = 1
         nnk = 0
         alfv=zero
         alfn=zero
      END IF

      gamma = one
      CALL vneg(n,g,d)
      
    END SUBROUTINE restar
    
      
    SUBROUTINE dobun(n,ma,mal,x,g,f,ay,ag,af,iters,ibun)  ! Bundle construction

      USE lmbm_sub, ONLY : copy2 ! given in host
      IMPLICIT NONE

! Array Arguments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
           x, &      ! Vector of variables
           g         ! Subgradient of the objective function.
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
           ay, &     ! Matrix whose columns are bundle points 
                     ! (stored in one-dimensional n*ma -array).
           ag, &     ! Matrix whose columns are bundle subgradients 
                     ! (stored in one-dimensional n*ma -array).
           af        ! Vector of values of bundle functions 
                     ! (stored in one-dimensional ma -array).

! Scalar Arguments
      INTEGER, INTENT(IN) :: & 
           n, &      ! Number of variables
           iters, &  ! Null step indicator.
                     !   0  - Null step.
                     !   1  - Serious step.
           ma        ! Maximum size of the bundle.
      INTEGER, INTENT(INOUT) :: &
           ibun, &   ! Index for the circular arrays in bundle updating.
           mal       ! Current size of the bundle.
      REAL(KIND=dp), INTENT(IN) :: & 
           f         ! Value of the objective function.

! Local Scalars
      INTEGER :: i,j
      
      IF (iters == 1) THEN
     
! Serious step

         af(ibun) = f
         i = (ibun-1)*n + 1
         CALL copy2(n,g,ag(i:),x,ay(i:))
         
      ELSE

! Null step

         IF (mal < ma) THEN
            
            af(ibun) = af(mal)
            af(mal) = f
            
            i = mal*n + 1
            CALL copy2(n,ag(i-n:),ag(i:),ay(i-n:),ay(i:))
            CALL copy2(n,g,ag(i-n:),x,ay(i-n:))
            
         ELSE
            i = ibun-1
            IF (i < 1) i = mal
            af(ibun) = af(i)
            af(i) = f
            
            i = (ibun-2)*n + 1
            IF (i < 1) i = (mal-1)*n + 1
            j = (ibun-1)*n + 1
            CALL copy2(n,ag(i:),ag(j:),ay(i:),ay(j:))
            CALL copy2(n,g,ag(i:),x,ay(i:))
         END IF
         
      END IF
      
      mal = mal + 1
      IF (mal > ma) mal = ma
      
      ibun = ibun + 1
      IF (ibun > ma) ibun = 1
      
    END SUBROUTINE dobun
    
  END SUBROUTINE lmbm

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE tinit *                                             *
!*                                                                      *
!*     Initial stepsize selection for limited memory bundle method      *
!*                                                                      *
!************************************************************************
  SUBROUTINE tinit(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,tmax,tmin,eta,iters,iterm)

    USE param, ONLY : zero,half,one,large
!    USE param, ONLY : zero,one  ! these are the one needed in tinit itself 
                                 ! half and large are used in destep and nulstep
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x, &      ! Vector of variables (n array).
         d, &      ! Direction vector (n array).
         ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
         ag        ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Scalar Arguments
    REAL(KIND=dp), INTENT(OUT) :: & 
         t         ! Initial stepsize
    REAL(KIND=dp), INTENT(IN) :: & 
         p, &      ! Directional derivative.
         eta, &    ! Distance measure parameter.
         f, &      ! Value of the objective function.
         tmax, &   ! Upper limit for stepsize parameter.
         tmin      ! Lower limit for stepsize parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         na, &     ! Maximum size of the bundle.
         mal, &    ! Current size of the bundle.
         ibun, &   ! Index for the circular arrays in bundle updating.
         iters     ! Null step indicator.
                   !   0  - Null step.
                   !   1  - Serious step.
    INTEGER, INTENT(OUT) :: &
         iterm     ! Cause of termination:
                   !   0  - Everything is ok.
                   !  -6  - Error.


! Intrinsic Functions
    INTRINSIC MAX,MIN

    t = MIN(one,tmax)

    IF (p == zero) RETURN

    IF (iters == 1) THEN
       CALL destep(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,eta,iterm)
    ELSE
       CALL nulstep(n,na,mal,x,af,ag,ay,ibun,d,f,p,t,eta,iterm)
    END IF

    t = MIN(MAX(t,tmin),tmax)

  CONTAINS

    SUBROUTINE destep(n,ma,mal,x,af,ag,ay,ibun,d,f,df,t,eta,iterm)  ! Stepsize selection

      USE param, ONLY : zero,half,one,large  ! given in host
      IMPLICIT NONE

! Scalar Arguments
      REAL(KIND=dp), INTENT(INOUT) :: & 
           t         ! Initial stepsize
      REAL(KIND=dp), INTENT(IN) :: & 
           df, &     ! Directional derivative.
           eta, &    ! Distance measure parameter.
           f         ! Value of the objective function.
      INTEGER, INTENT(IN) :: & 
           n, &      ! Number of variables
           ma, &     ! Maximum size of the bundle.
           mal, &    ! Current size of the bundle.
           ibun      ! Index for the circular arrays in bundle updating.
      INTEGER, INTENT(OUT) :: &
           iterm     ! Cause of termination:
                     !   0  - Everything is ok.
                     !  -6  - Error.

! Array Arguments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
           x, &      ! Vector of variables (n array).
           d, &      ! Direction vector (n array).
           ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
           ag, &     ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
           af        ! Vector of values of bundle functions (stored in one-dimensional ma -array).

! Local Arrays
      REAL(KIND=dp), DIMENSION(2*ma) :: &
           tmparray  ! Auxiliary array.

! Local Scalars
      REAL(KIND=dp) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
      INTEGER :: i,j,jn,k,l,lq,ib

! Intrinsic Functions
      INTRINSIC ABS,REAL,MAX,MIN,SQRT

      iterm = 0
      alfl = zero
      betl = zero
      
      w = df*t* (one - t*half)
      
     
! Initial choice of possibly active lines
      
      k = 0
      l = -1
      jn = (ibun-1)*n
      betr = - large
      DO j=1,mal-1
         ib = ibun - 1 + j
         IF (ib > mal) ib = ib - mal
         IF (jn >= mal*n) jn = jn - mal*n
         r = zero
         bet = zero
         alfl = af(ib) - f
         DO i=1,n
            dx = x(i) - ay(jn+i)
            q = ag(jn+i)
            r = r + dx*dx
            alfl = alfl + dx*q
            bet = bet + d(i)*q
         END DO
         alf = MAX(ABS(alfl),eta*r)
         r = one - bet/df
         IF (r*r + (alf+alf)/df > 1.0E-6_dp) THEN
            k = k + 1
            tmparray(k) = alf
            tmparray(ma+k) = bet
            r = t*bet - alf
            IF (r > w) THEN
               w = r
               l = k
            END IF
         END IF
         
         betr = MAX(betr,bet-alf)
         jn = jn + n
      END DO
      lq = -1
      IF (betr <= df*half) RETURN
      lq = 1
      IF (l <= 0) RETURN
      betr = tmparray(ma+l)
      IF (betr <= zero) THEN
         IF (t < one .OR. betr == zero) RETURN
         lq = 2
      END IF
      
      alfr = tmparray(l)


! Iteration loop

      ds_iteration: DO 
         
         IF (lq >= 1) THEN
            q = one - betr/df
            r = q + SQRT(q*q + (alfr+alfr)/df)
            IF (betr >= zero) r = - (alfr+alfr)/ (df*r)
            r = MIN(1.95_dp,MAX(zero,r))
         ELSE
            IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_dp*df) RETURN
            IF (betr-betl  == zero) THEN
               iterm = -6
               RETURN
            END IF
            r = (alfr-alfl)/ (betr-betl)
         END IF
       
         IF (ABS(t-r) < 1.0E-4_dp) RETURN
         t = r
         tmparray(l) = - one
         w = t*betr - alfr
         l = -1
         DO j = 1,k
            alf = tmparray(j)
            IF (alf < zero) EXIT
            bet = tmparray(ma+j)
            r = t*bet - alf
            IF (r > w) THEN
               w = r
               l = j
            END IF
         END DO
         IF (l < 0) RETURN
       
         bet = tmparray(ma+l)
         IF (bet == zero) RETURN
     

!     New interval selection

         alf = tmparray(l)
         IF (bet < zero) THEN
            IF (lq == 2) THEN
               alfr = alf
               betr = bet
               
            ELSE
               alfl = alf
               betl = bet
               lq = 0
            END IF
       
         ELSE
            IF (lq == 2) THEN
               alfl = alfr
               betl = betr
               lq = 0
            END IF
       
            alfr = alf
            betr = bet
         END IF
       
      END DO ds_iteration
      
    END SUBROUTINE destep

     
    SUBROUTINE nulstep(n,ma,mal,x,af,ag,ay,ibun,d,f,df,t,eta,iterm)  ! Stepsize selection

      USE param, ONLY : zero,half,one,large  ! given in host
      IMPLICIT NONE

! Scalar Arguments
      REAL(KIND=dp), INTENT(INOUT) :: & 
           t         ! Initial stepsize
      REAL(KIND=dp), INTENT(IN) :: & 
           df, &     ! Directional derivative.
           eta, &    ! Distance measure parameter.
           f         ! Value of the objective function.
      INTEGER, INTENT(IN) :: & 
           n, &      ! Number of variables
           ma, &     ! Maximum size of the bundle.
           mal, &    ! Current size of the bundle.
           ibun      ! Index for the circular arrays in bundle updating.
      INTEGER, INTENT(OUT) :: &
           iterm     ! Cause of termination:
                     !   0  - Everything is ok.
                     !  -6  - Error.

! Array Arguments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
           x, &      ! Vector of variables (n array).
           d, &      ! Direction vector (n array).
           ay, &     ! Matrix whose columns are bundle points (stored in one-dimensional n*ma -array).
           ag, &     ! Matrix whose columns are bundle subgradients (stored in one-dimensional n*ma -array).
           af        ! Vector of values of bundle functions (stored in one-dimensional 4*ma -array).

! Local Arrays
      REAL(KIND=dp), DIMENSION(2*ma) :: &
           tmparray  ! Auxiliary array.

! Local Scalars
      REAL(KIND=dp) :: alf,alfl,alfr,bet,betl,betr,dx,q,r,w
      INTEGER :: i,j,jn,k,l,ib

! Intrinsic Functions
      INTRINSIC ABS,REAL,MAX,MIN,SQRT

      iterm = 0
      w = df*t
     
! Initial choice of possibly active parabolas

      k = 0
      l = -1
      jn = (ibun-1)*n
      betr = - large
      DO j = 1,mal - 1
         ib = ibun - 1 + j
         IF (ib > mal) ib = ib - mal
         IF (jn >= mal*n) jn = jn - mal*n
         bet = zero
         r = zero
         alfl = af(ib) - f
         DO i = 1,n
            dx = x(i) - ay(jn+i)
            r = r + dx*dx
            q = ag(jn+i)
            alfl = alfl + dx*q
            bet = bet + d(i)*q
         END DO
         alf = MAX(ABS(alfl),eta*r)
         betr = MAX(betr,bet-alf)
         IF (alf < bet-df) THEN
            k = k + 1
            r = t*bet - alf
            tmparray(k) = alf
            tmparray(ma+k) = bet
            IF (r > w) THEN
               w = r
               l = k
            END IF
         END IF
         jn = jn + n
      END DO
      IF (l <= 0) RETURN
      betr = tmparray(ma+l)
      alfr = tmparray(l)
      alf = alfr
      bet = betr
      alfl = zero
      betl = df
    
     
!     Iteration loop
    
      ns_iteration: DO

         w = bet/df
         IF (ABS(betr-betl)+ABS(alfr-alfl)< -1.0E-4_dp*df) RETURN
         IF (betr-betl  == zero) THEN
            iterm = -6
            RETURN
         END IF
         r = (alfr-alfl)/ (betr-betl)
         IF (ABS(t-w) < ABS(t-r)) r = w
         q = t
         t = r
         IF (ABS(t-q) < 1.0E-3_dp) RETURN
         tmparray(l) = - one
         w = t*bet - alf
         l = -1
         DO j=1,k
            alf = tmparray(j)
            IF (alf < zero) EXIT
            bet = tmparray(ma+j)
            r = t*bet - alf
            IF (r > w) THEN
               w = r
               l = j
            END IF
         END DO

         IF (l <= 0) RETURN
         bet = tmparray(ma+l)
         q = bet - t*df
         IF (Q == zero) RETURN

     
!     New interval selection

         alf = tmparray(l)
         IF (q < zero) THEN
            alfl = alf
            betl = bet
         ELSE
            alfr = alf
            betr = bet
         END IF

      END DO ns_iteration

    END SUBROUTINE nulstep

  END SUBROUTINE tinit
      
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE lls *                                               *
!*                                                                      *
!*     Special line search for limited memory bundle method             *
!*                                                                      *
!************************************************************************

  SUBROUTINE lls(n,x,g,d,xo,t,fo,f,p,alfn,tmin,dnorm,wk,theta,epsl,epsr,&
       eta,iters,nfe,nge,nnk,iterm)

    USE param, ONLY : zero,half,one
    USE obj_fun, ONLY : &
         myf,myg   ! Computation of the value and the subgradient of
                   ! the objective function. 
    USE lmbm_sub, ONLY : &
         scsum, &  ! Sum of a vector and the scaled vector.
         vdot      ! Dot product of vectors.
    USE initializat, ONLY : &
        nproblem, &      ! the solved problem
        maxint => maxnin ! Maximum number of interpolations.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         xo        ! Previous vector of variables.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         x         ! Vector of variables.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         g         ! Subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: & 
         epsl,&    ! Linesearch parameter.
         epsr,&    ! Linesearch parameter.
         t, &      ! Stepsize
         p         ! Directional derivative.
    REAL(KIND=dp), INTENT(IN) :: & 
         theta, &  ! Scaling parameter.
         eta, &    ! Distance measure parameter.
         fo, &     ! Previous value of the objective function.
         wk, &     ! Stopping parameter.
         dnorm, &  ! Euclidean norm of the direction vector.
         tmin      ! Lower limit for stepsize parameter.
    REAL(KIND=dp), INTENT(OUT) :: & 
         f, &      ! Value of the objective function.
         alfn      ! Locality measure.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         nnk       ! Number of consequtive null steps.
    INTEGER, INTENT(INOUT) :: & 
         nfe, &    ! Number of function evaluations.
         nge       ! Number of subgradient evaluations.
    INTEGER, INTENT(OUT) :: &
         iters, &  ! Null step indicator.
                   !   0  - Null step.
                   !   1  - Serious step.
         iterm     ! Cause of termination:
                   !   0  - Everything is ok.
                   !  -3  - Failure in function or subgradient calculations.

! Local Scalars
    INTEGER :: nin ! Number of interpolations.
    REAL(KIND=dp) :: &
         tl,tu, &  ! Lower and upper limits for t used in interpolation.
         fl,fu, &  ! Values of the objective function for t=tl and t=tu.
         epsa,epst,epslk,epsrk, & ! Line search parameters.
         thdnorm,epsawk,epstwk,epslwk,epsrwk ! Auxiliary scalars.

! Parameters
!    INTEGER, PARAMETER :: maxint = 200  ! Maximum number of interpolations.

! Intrinsic Functions
    INTRINSIC ABS,MAX
!    INTEGER :: i
      

! Initialization

    nin = 0

    epst = epsl+epsl
    epsa = half*(epsr - epst)
    thdnorm = theta*dnorm

    tl = zero
    tu = t
    fl = fo

    IF (theta < one) THEN
       epst  = theta*epst
       epsa  = theta*epsa
       epslk = epsl
       epsl  = theta*epsl
       epsrk = epsr
       epsr  = theta*epsr
    END IF
 
    epsawk   = epsa*wk
    epslwk   = epsl*wk
    epsrwk   = epsr*wk
    epstwk   = epst*wk


! Function evaluation at a new point

    lls_iteration: DO
       
       CALL scsum(n,theta*t,d,xo,x)

       CALL myf(n,x,f,iterm,nproblem)
       nfe = nfe + 1
       IF (iterm /= 0) RETURN

    
! Null/descent step test (ITERS=0/1)

       iters = 1
       IF (f <= fo - t*epstwk) THEN
          tl = t
          fl = f
       ELSE
          tu = t
          fu = f
       END IF
      

! Additional interpolation

       IF (f > fo .AND. tu-tl >= tmin*0.1_dp &
            .AND. nnk >= 1 .AND. nin < maxint) THEN

          nin=nin+1
          IF (tl == zero .AND. wk > zero) THEN
             t = qint(tu,fl,fu,wk,one-half/(one-epst))
          ELSE
             t = half*(tu+tl)
          END IF
          CYCLE lls_iteration
       END IF

       CALL myg(n,x,g,iterm,nproblem)
       
       !WRITE(*,*)
       !WRITE(*,*) 'x:'
       !DO i =1, 31        
       !WRITE(*,*) 'x(',i,')=', x(i)      
       !WRITE(*,*) 'x(',i,')=', x(i)      
       !END DO     
       !WRITE(*,*)  
       !WRITE(*,*) 'g='            
       !DO i =1, 31        
       !WRITE(*,*) 'g(',i,')=', g(i)      
       !END DO  
       !WRITE(*,*) 
       nge = nge + 1
       IF (iterm /= 0) RETURN

       p = theta*vdot(n,g,d)
       alfn = MAX(ABS(fo-f+p*t),eta*(t*thdnorm)**2)


! Serious step

       IF (f <= fo - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lls_iteration


! Null step

       IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_dp .OR. &
            nin >= maxint) THEN
          ITERS = 0
          EXIT lls_iteration
       END IF


! Interpolation

       nin=nin+1
       IF (tl == zero .AND. wk > zero) THEN
          t = qint(tu,fl,fu,wk,one-half/(one-epst))
       ELSE
          t = half*(tu+tl)
       END IF

    END DO lls_iteration

    IF (theta /= one) THEN
       epsl = epslk
       epsr = epsrk
    END IF

  CONTAINS

    FUNCTION qint(tu,fl,fu,wk,kappa) RESULT(t) ! Quadratic interpolation

!      USE param, ONLY : half,one  ! given in host
      IMPLICIT NONE
      
! Scalar Arguments
      REAL(KIND=dp), INTENT(IN) :: & 
           fl, &  ! Value of the objective function.
           fu, &  ! Value of the objective function for t=tu.
           wk, &  ! Directional derivative.
           tu, &  ! Upper value of the stepsize parameter.
           kappa  ! Interpolation parameter.
      REAL(KIND=dp) :: & 
           t      ! Stepsize.

! Local Scalars
      REAL(KIND=dp) :: tmp1,tmp2

! Intrinsic Functions
      INTRINSIC MAX


      tmp1 = (fu-fl)/ (-wk*tu)

     
! Quadratic interpolation with one directional derivative

      tmp2 = 2.0_dp * (one - tmp1)

      IF (tmp2 > one) THEN
   
! Interpolation accepted
       
         t = MAX(kappa*tu,tu/tmp2)
         RETURN
      END IF
      
     
! Bisection
    
      t = half*tu
      
    END FUNCTION qint

  END SUBROUTINE lls

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE nmlls *                                             *
!*                                                                      *
!*     Nonmonotonic weak Wolfe-type line search                         *
!*                                                                      *
!************************************************************************

SUBROUTINE nmlls(n,x,g,d,xo,t,fo,f,fold,p,alfn,tmin,dnorm,wk,theta,epsl,epsr,&
    eta,iters,inewnma,nfe,nge,nnk,iterm)

    USE param, ONLY : zero,half,one
    USE obj_fun, ONLY : &
        myf,myg    ! Computation of the value and the subgradient of
                   ! the objective function.
    USE lmbm_sub, ONLY : &
        scsum      ! Sum of a vector and the scaled vector.
    USE initializat, ONLY : &
        nproblem, &  ! The solved problem
        mnma, &      ! Maximum number of function values used in line search.
        maxint => maxnin   ! Maximum number of interpolations.

    IMPLICIT NONE

    ! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
        d, &      ! Direction vector.
        xo        ! Previous vector of variables.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
        x         ! Vector of variables.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
        g         ! Subgradient of the objective function.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
        fold           ! Old function values.

    ! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: &
        epsl,&    ! Linesearch parameter.
        epsr,&    ! Linesearch parameter.
        t, &      ! Stepsize
        p         ! Directional derivative.
    REAL(KIND=dp), INTENT(IN) :: &
        theta, &  ! Scaling parameter.
        eta, &    ! Distance measure parameter.
        fo, &     ! Previous value of the objective function.
        wk, &     ! Stopping parameter.
        dnorm, &  ! Euclidean norm of the direction vector.
        tmin      ! Lower limit for stepsize parameter.
    REAL(KIND=dp), INTENT(OUT) :: &
        f, &      ! Value of the objective function.
        alfn      ! Locality measure.
    INTEGER, INTENT(IN) :: &
        n, &      ! Number of variables
        nnk       ! Number of consequtive null steps.
    INTEGER, INTENT(INOUT) :: &
        inewnma, &! index for array.
        nfe, &    ! Number of function evaluations.
        nge       ! Number of subgradient evaluations.
    INTEGER, INTENT(OUT) :: &
        iters, &  ! Null step indicator.
                  !   0  - Null step.
                  !   1  - Serious step.
        iterm     ! Cause of termination:
                  !   0  - Everything is ok.
                  !  -3  - Failure in function or subgradient calculations.

    ! Local Scalars
    INTEGER :: nin ! Number of interpolations.
    REAL(KIND=dp) :: &
        tl,tu, &  ! Lower and upper limits for t used in interpolation.
        fl,fu, &  ! Values of the objective function for t=tl and t=tu.
        ftmp, &   ! Maximum function value from mnma last iterations.
        epsa,epst,epslk,epsrk, & ! Line search parameters.
        thdnorm,epsawk,epstwk,epslwk,epsrwk ! Auxiliary scalars.

    ! Intrinsic Functions
    INTRINSIC ABS,MAX,MAXVAL,DOT_PRODUCT



    ! Initialization

    nin = 0

    epst = epsl+epsl
    epsa = half*(epsr - epst)
    thdnorm = theta*dnorm

    tl = zero
    tu = t
    fl = fo
    fu = fo

    ! Debugging: epslk and epsrk moved out from condition theta < one
    epslk = epsl
    epsrk = epsr
    IF (theta < one) THEN
        epst  = theta*epst
        epsa  = theta*epsa
        epsl  = theta*epsl
        epsr  = theta*epsr
    END IF

    epsawk   = epsa*wk
    epslwk   = epsl*wk
    epsrwk   = epsr*wk
    epstwk   = epst*wk

    ! Updating circular array fold

    IF (iters==1) THEN
        fold(inewnma) = fo
        inewnma = inewnma + 1
        IF (inewnma > mnma) inewnma = 1
    END IF
    ftmp = MAXVAL(fold)


    ! Function evaluation at a new point

    lls_iteration: DO

        CALL scsum(n,theta*t,d,xo,x)

        CALL myf(n,x,f,iterm,nproblem)
        nfe = nfe + 1
        IF (iterm /= 0) RETURN
        

        ! Null/descent step test (ITERS=0/1)

        iters = 1
        IF (f <= ftmp - t*epstwk) THEN
            tl = t
            fl = f
        ELSE
            tu = t
            fu = f
        END IF


        ! Additional interpolation

        IF (f > ftmp .AND. tu-tl >= tmin*0.1_dp &
            .AND. nnk >= 1 .AND. nin < maxint) THEN

            nin=nin+1
            IF (tl == zero .AND. wk > zero) THEN
                t = qint(tu,fl,fu,wk,one-half/(one-epst))
            ELSE
                t = half*(tu+tl)
            END IF
            CYCLE lls_iteration
        END IF

        CALL myg(n,x,g,iterm,nproblem)
        nge = nge + 1
        IF (iterm /= 0) RETURN

        p = theta*DOT_PRODUCT(g,d)
        alfn = MAX(ABS(fo-f+p*t),eta*(t*thdnorm)**2)


        ! Serious step

        IF (f <= ftmp - t*epslwk .AND. (t >= tmin .OR. alfn > epsawk)) EXIT lls_iteration


        ! Null step

        IF (p-alfn >= -epsrwk .OR. tu-tl < tmin*0.1_dp .OR. &
            nin >= maxint) THEN
            ITERS = 0
            EXIT lls_iteration
        END IF


        ! Interpolation

        nin=nin+1
        IF (tl == zero .AND. wk > zero) THEN
            t = qint(tu,fl,fu,wk,one-half/(one-epst))
        ELSE
            t = half*(tu+tl)
        END IF

    END DO lls_iteration

    IF (theta /= one) THEN
        epsl = epslk
        epsr = epsrk
    END IF

CONTAINS

    FUNCTION qint(tu,fl,fu,wk,kappa) RESULT(t) ! Quadratic interpolation

        USE param, ONLY : half,one  ! given in host
        IMPLICIT NONE

        ! Scalar Arguments
        REAL(KIND=dp), INTENT(IN) :: &
            fl, &  ! Value of the objective function.
            fu, &  ! Value of the objective function for t=tu.
            wk, &  ! Directional derivative.
            tu, &  ! Upper value of the stepsize parameter.
            kappa  ! Interpolation parameter.
        REAL(KIND=dp) :: &
            t      ! Stepsize.

        ! Local Scalars
        REAL(KIND=dp) :: tmp1,tmp2

        ! Intrinsic Functions
        INTRINSIC MAX


        tmp1 = (fu-fl)/ (-wk*tu)


        ! Quadratic interpolation with one directional derivative

        tmp2 = 2.0_dp * (one - tmp1)

        IF (tmp2 > one) THEN

            ! Interpolation accepted

            t = MAX(kappa*tu,tu/tmp2)
            RETURN
        END IF


        ! Bisection

        t = half*tu

    END FUNCTION qint

END SUBROUTINE nmlls



!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlbfgs *                                            *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*g  *
!*     by the limited memory BFGS update.                               *
!*                                                                      *
!************************************************************************
    
  SUBROUTINE dlbfgs(n,mc,mcc,inew,ibfgs,iflag,d,g,gp,s,u,sm,um,rm, &
       umtum,cm,smtgp,umtgp,gamma,tmpn1,method,iscale)
      
    USE param, ONLY : zero,small,one,half ! half is used at sclpar
    USE lmbm_sub, ONLY : &
         vdot, &   ! Dot product.
         vneg, &   ! Copying of a vector with change of the sign.
         xdiffy, & ! Difference of two vectors.
         xdiffy2, &! Difference of two vectors. (Variant)
         xsumy, &  ! Sum of two vectors.
         xsumy2, & ! Sum of two vectors. (Variant)
         scdiff, & ! Difference of the scaled vector and a vector.
         scdiff2, &! Difference of the scaled vector and a vector. (Variant)
         scsum, &  ! Sum of a vector and the scaled vector.
         scsum2, & ! Sum of a vector and the scaled vector. (Variant)
         vxdiag, & ! Multiplication of a vector and a diagonal matrix.
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.         
         trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
         trlieq2, &! Solving x from linear equation L*x=y or trans(L)*x=y. (Variant constructed for one fewer parameters in debugging)
         copy2     ! Copying of two vectors.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         g, &      ! Current subgradient of the objective function.
         gp        ! Previous subgradient of the objective function.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         s, &      ! Difference of current and previous variables.
         u, &      ! Difference of current and previous subgradients.
         tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         d         ! Direction vector.

! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: & 
         gamma     ! Scaling parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         method, & ! Selection of the method:
                   !   0  - Limited memory bundle method.
                   !   1  - L-BFGS bundle method.
         iscale    ! Selection of the scaling:
                   !   0  - Scaling at every iteration with STU/UTU.
                   !   1  - Scaling at every iteration with STS/STU.
                   !   2  - Interval scaling with STU/UTU.
                   !   3  - Interval scaling with STS/STU.
                   !   4  - Preliminary scaling with STU/UTU.
                   !   5  - Preliminary scaling with STS/STU.
                   !   6  - No scaling.      
    INTEGER, INTENT(INOUT) :: & 
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
    INTEGER, INTENT(OUT) :: & 
         ibfgs     ! Index of the type of BFGS update:
                   !   1  - BFGS update: the corrections are stored.
                   !   2  - BFGS update: the corrections are not stored.
                   !   3  - BFGS update is skipped.
      
! Local Arrays
    REAL(KIND=dp), DIMENSION(mcc+1) :: &
         tmpmc1,tmpmc2,tmpmc3,tmpmc4

! Local Scalars
    REAL(KIND=dp) :: &
         stu, &    ! stu = trans(s)*u. 
         sts       ! sts = trans(s)*s. 
    INTEGER :: i,j,k, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         iflag2, & ! Index for adaptive version.
         ierr      ! Error indicator

! Intrinsic Functions
    INTRINSIC SQRT,MIN,MAX


    ierr = 0
    ibfgs = 0
    iflag2 = 0
    stu = vdot(n,s,u)
    sts = vdot(n,s,s)


! Positive definiteness

    IF (stu > zero) THEN
       IF (-vdot(n,d,u)-vdot(n,tmpn1,s) < -small .OR. method == 1) THEN
          
     
! Update matrices
         
          ibfgs = 1

! Initialization of indices.

          CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
            
     
! Update sm and um

          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))
           
     
! Computation of trans(sm)*g and trans(um)*g

          IF (inew >= mcnew) THEN
             CALL rwaxv2(n,mcnew,sm((inew-mcnew)*n+1:),&
                  um((inew-mcnew)*n+1:),g,g,tmpmc1(iold:),tmpmc2(iold:))
          ELSE
             CALL rwaxv2(n,inew,sm,um,g,g,tmpmc1,tmpmc2)
             CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),&
                  um((iold-1)*n+1:),g,g,tmpmc1(iold:),tmpmc2(iold:))
          END IF
            
            
! Computation of trans(sm)*u and trans(um)*u

          IF (inew >= mcnew) THEN
             DO i=iold,inew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          ELSE
             DO i=1,inew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
             DO i=iold,mcnew+1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          END IF
          tmpmc3(inew) = tmpmc1(inew) - vdot(n,s,gp)
          smtgp(inew)  = tmpmc1(inew)
          tmpmc4(inew) = tmpmc2(inew) - vdot(n,u,gp)
          umtgp(inew)  = tmpmc2(inew)
            
         
! Update rm and umtum

          IF (mcc >= mc .AND. iflag2 /= 1) THEN
             DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
             END DO
          END IF
                      
          IF (inew >= mcnew) THEN
             CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),&
                  tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
          ELSE
             CALL copy2(mcnew-inew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:)&
                  ,tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
             CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),&
                  tmpmc4,umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
          END IF
            

! Update CM

          cm(inew) = stu
            
! Computation of gamma

          gamma = sclpar(mcc,iscale,method,sts,stu,tmpmc4(inew))
            
          inew = inew + 1
          IF (inew > mc + 1) inew = 1
          IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
          
       ELSE

            
! BFGS update, corrections are not saved.
     
          ibfgs = 2

! Initialization of indices.

          CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
          
! Update sm and um
          
          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))
            
! Computation of trans(sm)*g and trans(um)*g

          CALL rwaxv2(n,mcnew,sm,um,g,g,tmpmc1,tmpmc2)

! Computation of trans(sm)*u and trans(um)*u
          
          IF (iold /= 1) THEN
             DO i=1,inew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
             DO i=iold,mcnew
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          ELSE
             DO i=1,mcnew-1
                tmpmc3(i) = tmpmc1(i) - smtgp(i)
                smtgp(i)  = tmpmc1(i)
                tmpmc4(i) = tmpmc2(i) - umtgp(i)
                umtgp(i)  = tmpmc2(i)
             END DO
          END IF
          tmpmc3(inew) = tmpmc1(inew) - vdot(n,s,gp)
          smtgp(inew)  = tmpmc1(inew)
          tmpmc4(inew) = tmpmc2(inew) - vdot(n,u,gp)
          umtgp(inew)  = tmpmc2(inew)


! Update rm and umtum

          IF (iold /= 1) THEN
             CALL copy2(mcnew-inew,tmpmc3(iold:),&
                  rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                  umtum((mcnew-1)*mcnew/2+1:))
             CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                  umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
          ELSE
             CALL copy2(mcnew,tmpmc3,rm((mcnew-1)*mcnew/2+1:)&
                  ,tmpmc4,umtum((mcnew-1)*mcnew/2+1:))
          END IF
            

! Update cm

          cm(inew) = stu
            

! Computation of gamma

          gamma = sclpar(mcc,iscale,method,sts,stu,tmpmc4(inew))
               
       END IF
    ELSE
         
     
! BFGS update is skipped

       ibfgs = 3

       IF (mcc == 0) THEN
          iflag = 0
          CALL vneg(n,g,d)
          RETURN
       END IF
         

!     Initialization of indices.

       CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)


!     Computation of gamma

       IF (iscale >= 4) gamma = one
               
         
!     Computation of trans(sm)*g and trans(um)*g and the two intermediate values

       IF (iold <= 2) THEN
          CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),g,g,&
               smtgp(iold:),umtgp(iold:))
       ELSE
          CALL rwaxv2(n,inew-1,sm,um,g,g,smtgp,umtgp)
          CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),&
               um((iold-1)*n+1:),g,g,smtgp(iold:),umtgp(iold:))
       END IF
    END IF


! Computation of two intermediate values tmpmc1 and tmpmc2

    IF (iold == 1 .OR. ibfgs == 2) THEN
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,smtgp,1,ierr)
       CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc3)
       CALL vxdiag(mcnew,cm,tmpmc1,tmpmc2)
!       CALL scsum(mcnew,gamma,tmpmc3,tmpmc2,tmpmc2)
       CALL scsum2(mcnew,gamma,tmpmc3,tmpmc2) ! Variant
       CALL scsum(mcnew,-gamma,umtgp,tmpmc2,tmpmc3)
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc2,tmpmc3,0,ierr)

    ELSE IF (iflag == 0) THEN
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc1,smtgp,1,ierr)
       CALL symax(mcnew,mc+1,iold,umtum,tmpmc1,tmpmc3)
       CALL vxdiag(mc+1,cm,tmpmc1,tmpmc2)
!       CALL scsum(mc+1,gamma,tmpmc3,tmpmc2,tmpmc2)
       CALL scsum2(mc+1,gamma,tmpmc3,tmpmc2) ! Variant
       CALL scsum(mc+1,-gamma,umtgp,tmpmc2,tmpmc3)
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc2,tmpmc3,0,ierr)

    ELSE
       CALL trlieq(mcnew,mc,iold,rm,tmpmc1,smtgp,1,ierr)
       CALL symax(mcnew,mc,iold,umtum,tmpmc1,tmpmc3)
       CALL vxdiag(mc,cm,tmpmc1,tmpmc2)
!       CALL scsum(mc,gamma,tmpmc3,tmpmc2,tmpmc2)
       CALL scsum2(mc,gamma,tmpmc3,tmpmc2) ! Variant
       CALL scsum(mc,-gamma,umtgp,tmpmc2,tmpmc3)
       CALL trlieq(mcnew,mc,iold,rm,tmpmc2,tmpmc3,0,ierr)
    END IF

      
! Computation of the search direction d

    IF (iold == 1 .OR. ibfgs == 2) THEN
       CALL cwmaxv(n,mcnew,um,tmpmc1,d)
!       CALL xdiffy(n,d,g,d)
       CALL xdiffy2(n,d,g) ! Variant
       CALL cwmaxv(n,mcnew,sm,tmpmc2,tmpn1)
!       CALL scdiff(n,gamma,d,tmpn1,d)
       CALL scdiff2(n,gamma,d,tmpn1) ! Variant
    ELSE 
       CALL cwmaxv(n,inew-1,um,tmpmc1,d)
       CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc1(iold:),tmpn1)
!       CALL xsumy(n,d,tmpn1,d)
       CALL xsumy2(n,d,tmpn1)
!       CALL xdiffy(n,d,g,d)
       CALL xdiffy2(n,d,g)
       CALL cwmaxv(n,inew-1,sm,tmpmc2,tmpn1)
!       CALL scdiff(n,gamma,d,tmpn1,d)
       CALL scdiff2(n,gamma,d,tmpn1) ! Variant
       CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc2(iold:),tmpn1)
!       CALL xdiffy(n,d,tmpn1,d)
       CALL xdiffy2(n,d,tmpn1)
    END IF

  CONTAINS

    FUNCTION sclpar(mcc,iscale,method,sts,stu,utu) RESULT(spar) ! Calculation of the scaling parameter.
 !!!! ALLA KOMMENTOITU TAKAISIN  DEBUG testi   
      USE param, ONLY : small,one,half ! given in host
      IMPLICIT NONE

! Scalar Arguments
      REAL(KIND=dp), INTENT(IN) :: & 
           sts, &    ! sts = trans(s)*s. 
           stu, &    ! stu = trans(s)*u. 
           utu       ! utu = trans(u)*u. 
      REAL(KIND=dp) :: & 
           spar      ! Scaling parameter.
      INTEGER, INTENT(IN) :: & 
           mcc, &    ! Current number of stored corrections.
           method, & ! Selection of the method:
                     !   0  - Limited memory bundle method.
                     !   1  - L-BFGS bundle method.
           iscale    ! Selection of the scaling:
                     !   0  - Scaling at every iteration with STU/UTU.
                     !   1  - Scaling at every iteration with STS/STU.
                     !   2  - Interval scaling with STU/UTU.
                     !   3  - Interval scaling with STS/STU.
                     !   4  - Preliminary scaling with STU/UTU.
                     !   5  - Preliminary scaling with STS/STU.
                     !   6  - No scaling.      

! Intrinsic Functions
      INTRINSIC SQRT


! Computation of scaling parameter.

      SELECT CASE(iscale)
     
! Scaling parameter = STU/UTU

      CASE(0,2,4)
         IF (utu < SQRT(small)) THEN
            spar = one
            RETURN
         ELSE
            spar = stu/utu
         END IF
    
! Scaling parameter = STS/STU

      CASE(1,3,5)
         IF (stu < SQRT(small)) THEN
            spar = one
            RETURN
         ELSE
            spar = sts/stu
         END IF

! No scaling

      CASE DEFAULT
         spar = one
         RETURN
      END SELECT

               
!     Scaling
            
      IF (MCC == 0) THEN
         IF (spar < 0.01_dp) spar=0.01_dp
         IF (spar > 100.0_dp) spar=100.0_dp

      ELSE 

         SELECT CASE(iscale)

! Interval scaling
         CASE(2)
            IF (method == 0) THEN
               IF (spar < 0.6_dp .OR. spar > 6.0_dp) spar = one
            ELSE
               IF (spar < 0.01_dp .OR. spar > 100.0_dp) spar = one
            END IF

         CASE(3)
            IF (spar < half .OR. spar > 5.0_dp) spar = one
               
! Preliminary scaling
         CASE(4,5)
            spar = one

! Scaling at every iteration
         CASE DEFAULT
            CONTINUE
         END SELECT

      END IF

      IF (spar < 1.0E+03_dp*small) spar = 1.0E+03_dp*small
         
    END FUNCTION sclpar

  END SUBROUTINE dlbfgs

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlskip *                                            *
!*                                                                      *
!*     Matrix skipping and computation of the search direction          *
!*     d = -dm*ga by the limited memory BFGS update.                    *
!*                                                                      *
!************************************************************************

  SUBROUTINE dlskip(n,mc,mcc,inew,ibfgs,iflag,d,ga,sm,um,rm,&
       umtum,cm,tmpmc1,tmpmc2,gamma,iscale)
   
    USE param, ONLY : zero,small,one
    USE lmbm_sub, ONLY : &
         vneg, &   ! Copying of a vector with change of the sign.
         xsumy, &  ! Sum of two vectors.
         xsumy2, & ! Sum of two vectors. (Variant)
         xdiffy, & ! Difference of two vectors.
         xdiffy2, &! Difference of two vectors. (Variant)
         scsum, &  ! Sum of a vector and the scaled vector.
         scsum2, & ! Sum of a vector and the scaled vector. (Variant)
         scdiff, & ! Difference of the scaled vector and a vector.
         scdiff2, &! Difference of the scaled vector and a vector. (Variant)
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.         
         trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
         trlieq2, &! Solving x from linear equation L*x=y or trans(L)*x=y. (Variant constructed for one fewer parameters in debugging)
         vxdiag    ! Multiplication of a vector and a diagonal matrix.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         ga, &     ! Current aggregate subgradient of the objective function.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum     ! Matrix umtum = trans(um) * um.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         tmpmc1, & ! Auxiliary array. On output: trans(sm)*ga.
         tmpmc2, & ! Auxiliary array. On output: trans(um)*ga.
         d         ! Direction vector.

! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: & 
         gamma     ! Scaling parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         mc        ! Declared number of stored corrections.
    INTEGER, INTENT(INOUT) :: & 
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         iflag, &  ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
         iscale    ! Selection of the scaling:
                   !   0   - Scaling at every iteration with STU/UTU.
                   !   1  - Scaling at every iteration with STS/STU.
                   !   2  - Interval scaling with STU/UTU.
                   !   3  - Interval scaling with STS/STU.
                   !   4  - Preliminary scaling with STU/UTU.
                   !   5  - Preliminary scaling with STS/STU.
                   !   6  - No scaling.      
    INTEGER, INTENT(OUT) :: & 
         ibfgs     ! Index of the type of L-SR1 update:
                   !   3  - BFGS update is skipped.
      
! Local Arrays
    REAL(KIND=dp), DIMENSION(n) :: tmpn1
    REAL(KIND=dp), DIMENSION(mcc+1) :: tmpmc3,tmpmc4,tmpmc5
 
! Local Scalars
    INTEGER :: &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr,&    ! Error indicator.
         iflag2    ! Index for adaptive version.


! BFGS update is skipped

    ierr = 0
    ibfgs = 3

    IF (mcc == 0) THEN
       iflag = 0
       CALL vneg(n,ga,d)
       RETURN
    END IF
      
   
! Initialization of indices.
 
    CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,ibfgs)
      
         
! Computation of trans(sm)*ga and trans(um)*ga
    
    IF (iold <= 2) THEN
       CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),ga,ga,&
            tmpmc1(iold:),tmpmc2(iold:))
    ELSE
       CALL rwaxv2(n,inew-1,sm,um,ga,ga,tmpmc1,tmpmc2)
       CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
            ga,ga,tmpmc1(iold:),tmpmc2(iold:))
    END IF


! Computation of GAMMA

    IF (iscale >= 4) gamma = one


! Computation of two intermediate values tmpmc3 and tmpmc4

    IF (iold == 1) THEN
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc1,1,ierr)
       CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc5)
       CALL vxdiag(mcnew,cm,tmpmc3,tmpmc4)
!       CALL scsum(mcnew,gamma,tmpmc5,tmpmc4,tmpmc4)
       CALL scsum2(mcnew,gamma,tmpmc5,tmpmc4)
       CALL scsum(mcnew,-gamma,tmpmc2,tmpmc4,tmpmc5)
       CALL trlieq(mcnew,mcnew,iold,rm,tmpmc4,tmpmc5,0,ierr)

    ELSE IF (iflag == 0) THEN
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc3,tmpmc1,1,ierr)
       CALL symax(mcnew,mc+1,iold,umtum,tmpmc3,tmpmc5)
       CALL vxdiag(mc+1,cm,tmpmc3,tmpmc4)
!       CALL scsum(mc+1,gamma,tmpmc5,tmpmc4,tmpmc4)
       CALL scsum2(mc+1,gamma,tmpmc5,tmpmc4) ! Variant
       CALL scsum(mc+1,-gamma,tmpmc2,tmpmc4,tmpmc5)
       CALL trlieq(mcnew,mc+1,iold,rm,tmpmc4,tmpmc5,0,ierr)
    ELSE
       CALL trlieq(mcnew,mc,iold,rm,tmpmc3,tmpmc1,1,ierr)
       CALL symax(mcnew,mc,iold,umtum,tmpmc3,tmpmc5)
       CALL vxdiag(mc,cm,tmpmc3,tmpmc4)
!       CALL scsum(mc,gamma,tmpmc5,tmpmc4,tmpmc4)
       CALL scsum2(mc,gamma,tmpmc5,tmpmc4) ! Variant
       CALL scsum(mc,-gamma,tmpmc2,tmpmc4,tmpmc5)
       CALL trlieq(mcnew,mc,iold,rm,tmpmc4,tmpmc5,0,ierr)
    END IF
      

! Computation of the search direction d
    
    IF (iold == 1) THEN
       CALL cwmaxv(n,mcnew,um,tmpmc3,d)
    ELSE 
       CALL cwmaxv(n,inew-1,um,tmpmc3,d)
       CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
!       CALL xsumy(n,d,tmpn1,d)
       CALL xsumy2(n,d,tmpn1)
    END IF

!    CALL xdiffy(n,d,ga,d)
    CALL xdiffy2(n,d,ga) ! Variant
      
    IF (iold == 1) THEN
       CALL cwmaxv(n,mcnew,sm,tmpmc4,tmpn1)
!       CALL scdiff(n,gamma,d,tmpn1,d)
       CALL scdiff2(n,gamma,d,tmpn1) ! Variant
    ELSE
       CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
!       CALL scdiff(n,gamma,d,tmpn1,d)
       CALL scdiff2(n,gamma,d,tmpn1) ! Variant
       CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
!       CALL xdiffy(n,d,tmpn1,d)
       CALL xdiffy2(n,d,tmpn1) ! Variant
    END IF

  END SUBROUTINE dlskip
      
!************************************************************************
!*                                                                      *
!*     * SUBROUTINE dlsr1 *                                             *
!*                                                                      *
!*     Matrix update and computation of the search direction d = -dm*ga *
!*     by the limited memory SR1 update.                                *
!*                                                                      *
!************************************************************************

  SUBROUTINE dlsr1(n,mc,mcc,inew,isr1,iflag,d,gp,ga,s,u,sm,um,rm,&
       umtum,cm,smtgp,umtgp,gamma,tmpmc1,tmpmc2,tmpn1,nnk,iiprint)
      
    USE param, ONLY : zero,small,one
    USE lmbm_sub, ONLY : &
         vdot, &   ! Dot product.   
         vneg, &   ! Copying of a vector with change of the sign.
         scalex, & ! Scaling a vector.
         xdiffy, & ! Difference of two vectors.
         xdiffy2, &! Difference of two vectors. (Variant)
         scdiff, & ! Difference of the scaled vector and a vector.
         scdiff2, &! Difference of the scaled vector and a vector. (Variant)
         xsumy, &  ! Sum of two vectors.
         xsumy2, & ! Sum of two vectors. (Variant)
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.
         calq, &   ! Solving x from linear equation A*x=y.
         calq2, &  ! Solving x from linear equation A*x=y. (Variant)
         copy, &   ! Copying of a vector.
         copy2     ! Copying of two vectors.
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         ga, &     ! Current aggregate subgradient of the objective function.
         gp        ! Basic subgradient of the objective function.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         s, &      ! Difference of current and previous variables.
         u, &      ! Difference of current and previous subgradients.
         tmpn1     ! Auxiliary array. On input: previous aggregate subgradient.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         tmpmc1, & ! Auxiliary array. On output: trans(sm)*ga.
         tmpmc2, & ! Auxiliary array. On output: trans(um)*ga.
         d         ! Direction vector.

! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: & 
         gamma     ! Scaling parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         nnk, &    ! Consecutive null steps counter.
         iiprint    ! Printout specification.
    INTEGER, INTENT(INOUT) :: & 
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
    INTEGER, INTENT(OUT) :: & 
         isr1      ! Index of the type of L-SR1 update:
                   !   1  - SR1 update: the corrections are stored.
                   !   3  - SR1 update is skipped.
      
! Local Arrays
    REAL(KIND=dp), DIMENSION(n) :: tmpn2
    REAL(KIND=dp), DIMENSION((mcc+1)*(mcc+2)/2) :: tmpmat
    REAL(KIND=dp), DIMENSION(mcc+1) :: tmpmc3,tmpmc4,tmpmc5,tmpmc6

! Local Scalars
    REAL(KIND=dp) :: &
         stu, &    ! stu = trans(s)*u. 
         a, &      ! a = trans(ga) dm_(k-1) ga.
         b         ! b = trans(ga) dm_k ga.
    INTEGER :: i,j,k, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         iflag2    ! Index for adaptive version.


    iflag2 = 0
    isr1 = 0 
      

! Initialization of indices
      
    CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,3)
      
         
! Computation of gamma 

    gamma = one


! Computation of trans(sm)*ga and trans(um)*ga

    IF (iold <= 2) THEN
       CALL rwaxv2(n,mcnew,sm((iold-1)*n+1:),um((iold-1)*n+1:),ga,ga,&
            tmpmc1(iold:),tmpmc2(iold:))
    ELSE
       CALL rwaxv2(n,inew-1,sm,um,ga,ga,tmpmc1,tmpmc2)
       CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
            ga,ga,tmpmc1(iold:),tmpmc2(iold:))
    END IF


! Positive definiteness

    IF (-vdot(n,d,u) - vdot(n,tmpn1,s) >= -small) THEN
      
     
! SR1 update is skipped

       isr1 = 3
       
       IF (mcc == 0) THEN
          iflag = 0
          CALL vneg(n,ga,d)
          RETURN
       END IF

    ELSE
      
       stu = vdot(n,s,u)
        
       tmpmc1(inew) = vdot(n,s,ga)
       tmpmc2(inew) = vdot(n,u,ga)
      

!     Convergence conditions

       IF ((nnk == 1 .OR. mcc < mc) .OR. &
            (iflag == 1 .AND. (inew == 1 .OR. inew == mc))) THEN

! SR1 Update

          isr1 = 1
 

! Initialization of indices

          CALL indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,1)        

          if (iflag2 == 1 .and. iold == 2) then
             tmpmc1(inew) = tmpmc1(1)
             tmpmc2(inew) = tmpmc2(1)
          end if


! Update sm and um

          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))

         
! Update trans(sm)*gp and trans(um)*gp

          smtgp(inew) = vdot(n,s,gp)
          umtgp(inew) = vdot(n,u,gp)

     
! Computation of trans(sm)*u and trans(um)*u

          IF (iold <= 2) THEN
             CALL rwaxv2(n,mcnew-1,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                  u,u,tmpmc3(iold:),tmpmc4(iold:))
          ELSE
             CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
             CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
                  u,u,tmpmc3(iold:),tmpmc4(iold:))
          END IF

          tmpmc3(inew) = stu
          tmpmc4(inew) = vdot(n,u,u)

         
! Update rm and umtum

          IF (mcc >= mc .AND. iflag2 /= 1) THEN
             DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
             END DO
          END IF


          IF (inew >= mcnew) THEN
             CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),&
                  tmpmc4(iold:),umtum((mcnew-1)*mcnew/2+1:))
          ELSE
             CALL copy2(mcnew-inew,tmpmc3(iold:),&
                  rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                  umtum((mcnew-1)*mcnew/2+1:))
             CALL COPY2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                  umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
          END IF

      
! Update CM

          cm(inew) = stu

          inew = inew + 1
          IF (inew > mc + 1) inew = 1
          IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1

       ELSE


! Calculation of matrix (umtum-rm-trans(rm)+cm) from previous iteration
         
          DO  i=1,mcnew*(mcnew+1)/2
             tmpmat(i)= gamma * umtum(i) - rm(i)
          END DO

     
! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

          IF (iold == 1) THEN
             CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc5,0)

          ELSE IF (iflag == 0) THEN
             CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc5,0)

          ELSE
             CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
             CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc5,0)
          END IF


! Computation of a = -trans(ga)*dm_(k-1)*ga
      
          IF (iold <= 2) THEN
             CALL scalex(mcnew,gamma,tmpmc4(iold:),tmpmc3(iold:))
             CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
             CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
!             CALL xsumy(n,tmpn2,tmpn1,tmpn2)
             CALL xsumy2(n,tmpn2,tmpn1)
             
          ELSE
             CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
             CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,tmpn2)
             CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),&
                  tmpn1)
!             CALL xdiffy(n,tmpn2,tmpn1,tmpn2)
             CALL xdiffy2(n,tmpn2,tmpn1)
             CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
!             CALL xsumy(n,tmpn2,tmpn1,tmpn2)
             CALL xsumy2(n,tmpn2,tmpn1)
             CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),&
                  tmpn1)
!             CALL xsumy(n,tmpn2,tmpn1,tmpn2)
             CALL xsumy2(n,tmpn2,tmpn1) ! Variant
          END IF
          
          a = vdot(n,ga,tmpn2)
          
          IF (iflag == 0) THEN
             mcnew = mc
             iold = inew + 2
             IF (iold > mc+1) iold = iold - mc - 1
          ELSE
             mcnew = mc - 1
             iold = inew + 2
             IF (iold > mc) iold = iold - mc
          END IF
      

! Calculation of the new canditate for search direction
! Updates are not necessarily saved

! Update sm and um

          CALL copy2(n,s,sm((inew-1)*n+1:),u,um((inew-1)*n+1:))

     
! Computation of trans(sm)*u and trans(um)*u

          IF (iold == 1 .OR. iold == 2) THEN
             CALL rwaxv2(n,mcnew-1,sm((iold-1)*n+1:),um((iold-1)*n+1:),u,u,&
                  tmpmc3(iold:),tmpmc4(iold:))
          ELSE
             CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
             CALL rwaxv2(n,mcnew-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),u,u,&
                  tmpmc3(iold:),tmpmc4(iold:))
          END IF
       
          tmpmc3(inew) = stu
          tmpmc4(inew) = vdot(n,u,u)


! Calculation of matrix (umtum-rm-trans(rm)+cm) without updating
! matrices rm, umtum and cm
      
          DO i=1,mcnew*(mcnew+1)/2
             tmpmat(i)= gamma * umtum(i) - rm(i)
          END DO

          DO i=1,mcnew-1
             j=(i-1)*i/2+1
             k=i*(i+1)/2+2
             CALL copy(i,tmpmat(k:),tmpmat(j:))
          END DO
         
          CALL scdiff(mcnew+1,gamma,tmpmc4,tmpmc3,tmpmc5)
         
          IF (inew >= mcnew) THEN
             CALL copy(mcnew,tmpmc5(iold:),tmpmat((mcnew-1)*mcnew/2+1:))
          ELSE
             CALL copy(mcnew-inew,tmpmc5(iold:),tmpmat((mcnew-1)*mcnew/2+1:))
             CALL copy(inew,tmpmc5,tmpmat((mcnew-1)*mcnew/2+mcnew-inew+1:))
          END IF
      
          IF (iflag == 0) THEN
             CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc5)
!             CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc5,tmpmc5,iiprint)
             CALL calq2(mcnew,mc+1,iold,tmpmat,tmpmc5,iiprint) ! Variant

          ELSE
             CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc5)
!             CALL calq(mcnew,mc,iold,tmpmat,tmpmc5,tmpmc5,iiprint)
             CALL calq2(mcnew,mc,iold,tmpmat,tmpmc5,iiprint) ! Variant
          END IF


! Calculation of the new canditate for search direction d = -dm_k*ga
! and computation of b = -trans(ga)*dm_k*ga
      
          IF (iold <= 2) THEN
             CALL scalex(mcnew,gamma,tmpmc5(iold:),tmpmc6(iold:))
             CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc5(iold:),tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,d)
             CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc6(iold:),tmpn1)
!             CALL xsumy(n,d,tmpn1,d)
             CALL xsumy2(n,d,tmpn1) ! Variant
         
          ELSE
             CALL scalex(mcnew+1,gamma,tmpmc5,tmpmc6)
             CALL cwmaxv(n,inew,sm,tmpmc5,tmpn1)
             CALL scdiff(n,-gamma,ga,tmpn1,d)
             CALL cwmaxv(n,mcnew-inew,sm((iold-1)*n+1:),tmpmc5(iold:),&
                  tmpn1)
!             CALL xdiffy(n,d,tmpn1,d)
             CALL xdiffy2(n,d,tmpn1) ! Variant
             CALL cwmaxv(n,inew,um,tmpmc6,tmpn1)
!             CALL xsumy(n,d,tmpn1,d)
             CALL xsumy2(n,d,tmpn1) ! Variant
             CALL cwmaxv(n,mcnew-inew,um((iold-1)*n+1:),tmpmc6(iold:),&
                  tmpn1)
!             CALL xsumy(n,d,tmpn1,d)
             CALL xsumy2(n,d,tmpn1) ! Variant
          END IF

          b = vdot(n,ga,d)


! Checking the convergence conditions

          IF (b - a < zero) THEN
             isr1 = 3
             CALL copy(n,tmpn2,d)
            
          ELSE

             isr1 = 1
         
     
! Update trans(sm)*gp and trans(um)*gp
     
             smtgp(inew) = vdot(n,s,gp)
             umtgp(inew) = vdot(n,u,gp)

                     
! Update rm and umtum

             DO i=1,mcnew-1
                j=(i-1)*i/2+1
                k=i*(i+1)/2+2
                CALL copy2(i,rm(k:),rm(j:),umtum(k:),umtum(j:))
             END DO

             IF (inew >= mcnew) THEN
                CALL copy2(mcnew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                     umtum((mcnew-1)*mcnew/2+1:))
             ELSE
                CALL copy2(mcnew-inew,tmpmc3(iold:),rm((mcnew-1)*mcnew/2+1:),tmpmc4(iold:),&
                     umtum((mcnew-1)*mcnew/2+1:))
                CALL copy2(inew,tmpmc3,rm((mcnew-1)*mcnew/2+mcnew-inew+1:),tmpmc4,&
                     umtum((mcnew-1)*mcnew/2+mcnew-inew+1:))
             END IF
            

! Update cm

             cm(inew) = stu
                     
             inew = inew + 1
             IF (inew > mc + 1) inew = 1
             IF (iflag == 0 .AND. mcc < mc + 1) mcc = mcc + 1
            
          END IF
         
          RETURN
         
       END IF
    END IF
      
    DO i=1,mcnew*(mcnew+1)/2
       tmpmat(i)= gamma * umtum(i) - rm(I)
    END DO
      
     
! Computation of tmpmat*tmpmc4 = gamma*trans(um)*ga-trans(sm)*ga

    IF (iold == 1) THEN
       CALL scdiff(mcnew,gamma,tmpmc2,tmpmc1,tmpmc4)
!       CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc4,tmpmc4,iiprint)
       CALL calq2(mcnew,mcnew,iold,tmpmat,tmpmc4,iiprint) ! Variant
    ELSE IF (iflag == 0) THEN
       CALL scdiff(mc+1,gamma,tmpmc2,tmpmc1,tmpmc4)
!       CALL calq(mcnew,mc+1,iold,tmpmat,tmpmc4,tmpmc4,iiprint)
       CALL calq2(mcnew,mc+1,iold,tmpmat,tmpmc4,iiprint) ! Variant
    ELSE
       CALL scdiff(mc,gamma,tmpmc2,tmpmc1,tmpmc4)
!       CALL calq(mcnew,mc,iold,tmpmat,tmpmc4,tmpmc4,iiprint)
       CALL calq2(mcnew,mc,iold,tmpmat,tmpmc4,iiprint) ! Variant
    END IF
      

! Computation of the search direction d
      
    IF (iold <= 2) THEN
       CALL scalex(mcnew,gamma,tmpmc4(iold:),tmpmc3(iold:))
       CALL cwmaxv(n,mcnew,sm((iold-1)*n+1:),tmpmc4(iold:),tmpn1)
       CALL scdiff(n,-gamma,ga,tmpn1,d)
       CALL cwmaxv(n,mcnew,um((iold-1)*n+1:),tmpmc3(iold:),tmpn1)
!       CALL xsumy(n,d,tmpn1,d)
       CALL xsumy2(n,d,tmpn1) ! Variant
    ELSE
       CALL scalex(mcc,gamma,tmpmc4,tmpmc3)
       CALL cwmaxv(n,inew-1,sm,tmpmc4,tmpn1)
       CALL scdiff(n,-gamma,ga,tmpn1,d)
       CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc4(iold:),&
            tmpn1)
!       CALL xdiffy(n,d,tmpn1,d)
       CALL xdiffy2(n,d,tmpn1) ! Variant
       CALL cwmaxv(n,inew-1,um,tmpmc3,tmpn1)
!       CALL xsumy(n,d,tmpn1,d)
       CALL xsumy2(n,d,tmpn1) ! Variant
       CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc3(iold:),&
            tmpn1)
!       CALL xsumy(n,d,tmpn1,d)
       CALL xsumy2(n,d,tmpn1) ! Variant
    END IF
      
  END SUBROUTINE dlsr1

!************************************************************************
!*                                                                      *
!*     * SUBROUTINE indic1 *                                            *
!*                                                                      *
!*     Initialization of indices.                                       *
!*                                                                      *
!************************************************************************

  SUBROUTINE indic1(mc,mcc,mcnew,inew,iold,iflag,iflag2,itype)

    IMPLICIT NONE

! Scalar Arguments

    INTEGER, INTENT(IN) :: & 
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         itype     ! Type of Initialization:
                   !   1  - corrections are stored,
                   !   2  - corrections are not stored,
                   !   3  - update is skipped.
    INTEGER, INTENT(INOUT) :: & 
         inew, &   ! Index for circular arrays.
         iflag, &  ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
         iflag2    ! Index for adaptive version.
                   !   0  - iflag has not been changed.
                   !   1  - iflag has been changed.
    INTEGER, INTENT(OUT) :: & 
         mcnew, &  ! Current size of vectors.
         iold      ! Index of the oldest corrections.
      
    IF (itype == 1) THEN
       IF (mcc < mc) THEN
          mcnew = mcc + 1
          iold = 1
          iflag = 0
       ELSE
          IF (iflag == 0) THEN
             mcnew = mc
             iold = inew + 2
             IF (iold > mc+1) iold = iold - mc - 1
          ELSE
             IF (inew == 1) THEN
                inew = mc + 1
                mcnew = mc
                iold = 2
                iflag = 0
                iflag2 = 1
             ELSE IF (inew == mc) THEN
                mcnew = mc
                iold = 1
                iflag = 0
                iflag2 = 1
             ELSE
                mcnew = mc - 1
                iold = inew + 2
                IF (iold > mc) iold = iold - mc
             END IF
          END IF
       END IF
      
    ELSE IF (itype == 2) THEN
       IF (mcc < mc) THEN
          mcnew = mcc + 1
          iold = 1
          iflag = 0
       ELSE
          IF (iflag == 0) THEN
             mcnew = mc + 1
             iold = inew + 1
             IF (iold > mc + 1) iold = 1
          ELSE
             mcnew = mc
             iold = inew + 1
             IF (iold > mc) iold = 1
          END IF
       END IF

    ELSE 
       IF (mcc < mc) THEN
          mcnew = mcc
          iold = 1
          iflag = 0
       ELSE
          IF (iflag == 0) THEN
             mcnew = mc
             iold = inew + 1
             IF (iold > mc + 1) iold = 1
          ELSE
             mcnew = mc - 1
             iold = inew + 1
             IF (iold > mc) iold = 1
          END IF
       END IF
    END IF
  END SUBROUTINE indic1
  
!************************************************************************
!*
!*     * SUBROUTINE agbfgs *
!*
!*     Computation of aggregate values by the limited memory BFGS update.
!*
!************************************************************************

  SUBROUTINE agbfgs(n,mc,mcc,inew,ibfgs,iflag,g,gp,ga,u,d,sm,um, &
       rm,cm,umtum,alfn,alfv,gamma,ic,rho)

    USE param, ONLY : zero,half,one
    USE lmbm_sub, ONLY : &
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.
         trlieq, & ! Solving x from linear equation L*x=y or trans(L)*x=y.
         trlieq2, &! Solving x from linear equation L*x=y or trans(L)*x=y. (Variant built in debugging with one fewer parameter)
         vdot      ! Dot product.

    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         g, &      ! Current (auxiliary) subgradient of the objective function.
         gp, &     ! Previous subgradient of the objective function.
         u, &      ! Difference of trial and aggregate gradients.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         cm        ! Diagonal matrix.
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: &
         ga        ! Next aggregate subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=dp), INTENT(OUT) :: & 
         alfv      ! Aggregate locality measure.
    REAL(KIND=dp), INTENT(IN) :: & 
         gamma, &  ! Scaling parameter.
         alfn, &   ! Locality measure.
         rho       ! Correction parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         ibfgs, &  ! Index of the type of BFGS update.
         ic        ! Correction indicator.
    INTEGER, INTENT(INOUT) :: & 
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.

! Local arrays
    REAL(KIND=dp), DIMENSION(mcc+1) :: tmpmc1, tmpmc2

! Local Scalars
    REAL(KIND=dp) :: &
         p, &      ! p = trans(d)*u - alfn.
         q, &      ! q = trans(u)*dm*u, where dm is the inverse approximation of 
                   ! the Hessian calculated by using the L-BFGS formula.
         lam, &    ! Multiplier used to calculate aggregate values.
         w         ! Correction.
    INTEGER :: i, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr      ! Error indicador.

! Intrinsic Functions
    INTRINSIC MAX,MIN,SIGN

    ierr = 0

    IF (mcc < mc) THEN
       IF (ibfgs == 2) THEN
          mcnew = mcc + 1
       ELSE
          mcnew = mcc
       END IF
       iold = 1

    ELSE
       IF (iflag == 0) THEN
          IF (ibfgs == 2) THEN
             mcnew = mc + 1
          ELSE
             mcnew = mc
          END IF
          iold = inew + 1
          IF (iold > mc+1) iold = 1

       ELSE
          IF (ibfgs == 2) THEN
             mcnew = mc
          ELSE
             mcnew = mc - 1
          END IF
          iold = inew + 1
          IF (iold > mc) iold = 1
       END IF
    END IF
      
      
! Computation of trans(d) * u - alfn

    p = vdot(n,d,u) - alfn
    q = vdot(n,u,u)

    IF (ic == 1) THEN
       w = rho * q
    ELSE
       w = zero
    END IF
         
     
! Computation of the product trans(u)*dm*u

    IF (mcc > 0 .OR. ibfgs == 2) THEN

       IF (iold == 1 .OR. ibfgs == 2) THEN
          CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc1,tmpmc2)
!          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc1,tmpmc1,1,ierr)
          CALL trlieq2(mcnew,mcnew,iold,rm,tmpmc1,1,ierr) ! Debugged version with no compiler warnings

          q = q - 2.0_dp*vdot(mcnew,tmpmc2,tmpmc1)
          q = gamma*q
            
          DO i=1,mcnew
             tmpmc2(i) = cm(i)*tmpmc1(i)
          END DO

          q = q + vdot(mcnew,tmpmc1,tmpmc2)

          CALL symax(mcnew,mcnew,iold,umtum,tmpmc1,tmpmc2)

          q = q + gamma*vdot(mcnew,tmpmc1,tmpmc2)

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc1,tmpmc2)
          CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:), &
               u,u,tmpmc1(iold:),tmpmc2(iold:))
!          CALL trlieq(mcnew,mcc,iold,rm,tmpmc1,tmpmc1,1,ierr)
          CALL trlieq2(mcnew,mcc,iold,rm,tmpmc1,1,ierr) ! Debugged version with no compiler warnings

          q = q - 2.0_dp*(vdot(mcc-inew,tmpmc2(iold:),tmpmc1(iold:)) + &
               vdot(inew-1,tmpmc2,tmpmc1))
          q = gamma*q

          DO i=1,mcc
             tmpmc2(i) = cm(i)*tmpmc1(i)
          END DO

          q = q + vdot(mcc-inew,tmpmc1(iold:),tmpmc2(iold:)) + &
               vdot(inew-1,tmpmc1,tmpmc2)

          CALL symax(mcnew,mcc,iold,umtum,tmpmc1,tmpmc2)

          q = q + gamma*(vdot(mcc-inew,tmpmc1(iold:),tmpmc2(iold:)) + &
               vdot(inew-1,tmpmc1,tmpmc2))
       END IF

    END IF
    
    q = q + w
    
    lam = half + SIGN(half,p)

    IF (q > zero) lam = MIN(one,MAX(zero,p/q))
      

! Computation of the aggregate values

    p = one - lam
    DO i=1,n
       ga(i)=lam*g(i) + p*gp(i)
    END DO
      
    alfv = lam*alfn
      
  END SUBROUTINE agbfgs

!************************************************************************
!*
!*     * SUBROUTINE aggsr1 *
!*
!*     Computation of aggregate values by the limited memory SR1 update.
!*
!************************************************************************
      
  SUBROUTINE aggsr1(n,mc,mcc,inew,iflag,g,gp,ga,d,alfn,alfv, &
       umtum,rm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)
      
    USE param, ONLY : zero,one,small
    USE lmbm_sub, ONLY : &
         vdot, &   ! Dot product.
         scalex, & ! Scaling a vector.
         xsumy, &  ! Sum of two vectors.
         xsumy2, & ! Sum of two vectors. (Variant)
         xdiffy, & ! Difference of two vectors.
         xdiffy2, &! Difference of two vectors. (Variant)
         scsum, &  ! Sum of a vector and the scaled vector.
         scsum2, & ! Sum of a vector and the scaled vector. (Variant)
         scdiff, & ! Difference of the scaled vector and a vector.
         scdiff2, &! Difference of the scaled vector and a vector. ! Variant with INOUT
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.   
         cwmaxv, & ! Multiplication of a vector by a dense rectangular matrix.
         lineq, &  ! Solver for linear equation.
         calq, &   ! Solving x from linear equation A*x=y.
         calq2     ! Solving x from linear equation A*x=y. (Variant)
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         g, &      ! Current (auxiliary) subgradient of the objective function.
         gp, &     ! Previous subgradient of the objective function.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         smtga, &  ! vector smtga = trans(sm)*ga.
         umtga     ! vector umtga = trans(um)*ga.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         ga        ! Aggregate subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: & 
         alfv      ! Aggregate locality measure.
    REAL(KIND=dp), INTENT(IN) :: & 
         gamma, &  ! Scaling parameter.
         alfn, &   ! Locality measure.
         rho       ! Correction parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         icn       ! Correction indicator.
    INTEGER, INTENT(INOUT) :: & 
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
      
! Local arrays
    REAL(KIND=dp), DIMENSION(n) :: tmpn2,tmpn3,tmpn4
    REAL(KIND=dp), DIMENSION((mcc+1)*(mcc)/2) :: tmpmat
    REAL(KIND=dp), DIMENSION(mcc+1) :: tmpmc3, tmpmc4

! Local Scalars
    REAL(KIND=dp) :: &
         pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm
                   ! presents the L-SR1- approximation of Hessian.
         rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
         prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
         rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
         qr, &     ! qr = trans(g-ga) dm (g-ga).
         pq, &     ! pq = trans(g-gp) dm (g-gp).
         qqp, &    ! qqp = trans(g-gp) dm g + alfn.
         lam1, &   ! Multiplier used to calculate aggregate values.
         lam2, &   ! Multiplier used to calculate aggregate values.
         w, &      ! Correction.
         tmp1, &   ! Auxiliary scalar.
         tmp2      ! Auxiliary scalar.
    INTEGER :: i, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr      ! Error indicador.

      
! Intrinsic Functions
    INTRINSIC MIN,MAX


    ierr = 0
      
    IF (mcc < mc) THEN
       iold = 1
       mcnew = mcc
    ELSE IF (iflag == 0) THEN
       mcnew = mc
       iold = inew + 1
       IF (iold > mc+1) iold = 1
    ELSE
       mcnew = mc - 1
       iold = inew + 1
       IF (iold > mc) iold = 1
    END IF
      
    CALL xdiffy(n,gp,ga,tmpn2)
     
      
! Calculation of tmpn3 = trans(gp - ga)dm

    IF (mcc > 0) THEN

       DO i=1,mcnew*(mcnew+1)/2
          tmpmat(i)= gamma * umtum(i) - rm(i)
       END DO

       IF (iold == 1) THEN
          CALL xdiffy(mcnew,umtgp,umtga,tmpmc4)
!          CALL scdiff(mcnew,gamma,tmpmc4,smtgp,tmpmc4)
          CALL scdiff2(mcnew,gamma,tmpmc4,smtgp) ! Variant with INOUT
!          CALL xsumy(mcnew,tmpmc4,smtga,tmpmc4)
          CALL xsumy2(mcnew,tmpmc4,smtga) ! Variant
          
          CALL calq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,0)
          CALL scalex(mcnew,gamma,tmpmc3,tmpmc4)
          
          CALL cwmaxv(n,mcnew,sm,tmpmc3,tmpn4)
          CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
          CALL cwmaxv(n,mcnew,um,tmpmc4,tmpn4)
!          CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
          CALL xdiffy2(n,tmpn3,tmpn4) ! Variant

       ELSE
          CALL xdiffy(mcc,umtgp,umtga,tmpmc4)
!          CALL scdiff(mcc,gamma,tmpmc4,smtgp,tmpmc4)
          CALL scdiff2(mcc,gamma,tmpmc4,smtgp) ! Variant with INOUT
!          CALL xsumy(mcc,tmpmc4,smtga,tmpmc4)
          CALL xsumy2(mcc,tmpmc4,smtga) ! Variant
            
          CALL calq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,0)
          CALL scalex(mcc,gamma,tmpmc3,tmpmc4)
          
          CALL cwmaxv(n,inew-1,sm,tmpmc3,tmpn4)
          CALL scsum(n,gamma,tmpn2,tmpn4,tmpn3)
          CALL cwmaxv(n,mcnew-inew+1,sm((iold-1)*n+1:),tmpmc3(iold:)&
               ,tmpn4)
!          CALL xsumy(n,tmpn3,tmpn4,tmpn3)
          CALL xsumy2(n,tmpn3,tmpn4) ! Variant
          CALL cwmaxv(n,inew-1,um,tmpmc4,tmpn4)
!          CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
          CALL xdiffy2(n,tmpn3,tmpn4)
          CALL cwmaxv(n,mcnew-inew+1,um((iold-1)*n+1:),tmpmc4(iold:)&
               ,tmpn4)
!          CALL xdiffy(n,tmpn3,tmpn4,tmpn3)
          CALL xdiffy2(n,tmpn3,tmpn4)
       END IF
       
       IF (icn == 1) THEN
!          CALL scsum(n,rho,tmpn2,tmpn3,tmpn3)
          CALL scsum2(n,rho,tmpn2,tmpn3)
       END IF
         
       pr = vdot(n,tmpn3,tmpn2)
       rrp = vdot(n,tmpn3,ga) 
       CALL xdiffy(n,g,ga,tmpn4)
       prqr = vdot(n,tmpn3,tmpn4)
       rrq = -vdot(n,tmpn4,d)

    ELSE

       pr = vdot(n,tmpn2,tmpn2)
       rrp = vdot(n,tmpn2,ga) 
       CALL xdiffy(n,g,ga,tmpn4)
       prqr = vdot(n,tmpn2,tmpn4)
       rrq = -vdot(n,tmpn4,d)
    END IF

! calculation of qr = trans(g - ga) dm (g - ga)

    qr = vdot(n,tmpn4,tmpn4)
    IF (icn == 1) THEN
       w = rho*qr
    ELSE
       w = zero
    END IF
      
    IF (mcc > 0) THEN
       qr = gamma*qr

       IF (iold == 1) THEN
          CALL rwaxv2(n,mcnew,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
!          CALL scsum(mcnew,-gamma,tmpmc3,tmpmc4,tmpmc4)
          CALL scsum2(mcnew,-gamma,tmpmc3,tmpmc4)
          CALL lineq(mcnew,mcnew,iold,tmpmat,tmpmc3,tmpmc4,ierr)
            
          qr = qr - vdot(mcnew,tmpmc4,tmpmc3) + w

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,tmpn4,tmpn4,tmpmc4,tmpmc3)
          CALL rwaxv2(n,mcnew-inew+1,sm((iold-1)*n+1:),&
               um((iold-1)*n+1:),tmpn4,tmpn4,tmpmc4(iold:),tmpmc3(iold:))
!          CALL scsum(mcc,-gamma,tmpmc3,tmpmc4,tmpmc4)
          CALL scsum2(mcc,-gamma,tmpmc3,tmpmc4)
          CALL lineq(mcnew,mcc,iold,tmpmat,tmpmc3,tmpmc4,ierr)
          
          qr = qr - vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) -&
               vdot(inew-1,tmpmc4,tmpmc3) + w
       END IF

    END IF
      
    pq = qr - prqr - prqr + pr
    qqp = pq + prqr + rrq - pr - rrp + alfn
    rrp = rrp - alfv
    rrq = rrq + alfn - alfv

! computation of multipliers lam1 and lam2

    IF (pr > zero .AND. qr > zero) THEN
       tmp1 = rrq/qr
       tmp2 = prqr/qr
       w = pr - prqr*tmp2
       IF (w /= zero) THEN
          lam1 = (tmp1*prqr - rrp)/w
          lam2 = -tmp1 - lam1*tmp2
          IF (lam1*(lam1 - one) < zero .AND. &
               lam2*(lam1 + lam2 - one) < zero) GO TO 200
       END IF
    END IF

! Minimum on the boundary

!100 continue
    continue
    lam1 = zero
    lam2 = zero
    IF (alfn <= alfv) lam2 = one
    IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
    w = (lam2*qr + rrq+rrq)*lam2
!    w = (lam2*qr + 2.0_dp*rrq)*lam2
    tmp1 = zero
    IF (alfv >= zero) tmp1 = one
    IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
!    tmp2 = (tmp1*pr + 2.0_dp*rrp)*tmp1
    tmp2 = (tmp1*pr + rrp+rrp)*tmp1
    IF (tmp2 < w) THEN
       w = tmp2
       lam1 = tmp1
       lam2 = zero
    END IF
    
    IF (qqp*(qqp - pq) < zero) THEN
       IF (qr + rrq + rrq - qqp*qqp/pq < W) THEN
          lam1 = qqp/pq
          lam2 = one - lam1
       END IF
    END IF
    
200 CONTINUE
    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
         .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
         lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)

! Computation of the aggregate values
      
    tmp1 = one - lam1 - lam2
    DO i=1,n
       ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
    END DO
    
    alfv = lam2*alfn + tmp1*alfv
    
  END SUBROUTINE aggsr1
      
!************************************************************************
!*
!*     * SUBROUTINE agskip *
!*
!*     Computation of aggregate values after consecutive null steps
!*     by the limited memory BFGS update.
!*
!************************************************************************
      
  SUBROUTINE agskip(n,mc,mcc,inew,iflag,g,gp,ga,d,u,alfn,alfv, &
       umtum,rm,cm,gamma,smtgp,umtgp,smtga,umtga,sm,um,icn,rho)
      
    USE param, ONLY : zero,half,one,small
    USE lmbm_sub, ONLY : &
         xdiffy, & ! Difference of two vectors.
         xdiffy2, &! Difference of two vectors. (Variant)
         symax, &  ! Multiplication of a dense symmetric matrix by a vector.
         rwaxv2, & ! Multiplication of two rowwise stored dense rectangular 
                   ! matrices A and B by vectors x and y.         
         trlieq, & ! Solving x from linear equation l*x=y or trans(l)*x=y.
         trlieq2, &! Solving x from linear equation l*x=y or trans(l)*x=y. (Debugging variant)
         vdot      ! Dot product
    IMPLICIT NONE

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         d, &      ! Direction vector.
         g, &      ! Current (auxiliary) subgradient of the objective function.
         gp, &     ! Previous subgradient of the objective function.
         u, &      ! Difference of trial and aggregate gradients.
         sm, &     ! Matrix whose columns are stored corrections.
         um, &     ! Matrix whose columns are stored subgradient differences.
         rm, &     ! Upper triangular matrix.
         cm, &     ! Diagonal matrix.
         umtum, &  ! Matrix umtum = trans(um) * um.
         smtgp, &  ! Vector smtgp = trans(sm)*gp.
         umtgp, &  ! vector umtgp = trans(um)*gp.
         smtga, &  ! vector smtga = trans(sm)*ga.
         umtga     ! vector umtga = trans(um)*ga.
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: &
         ga        ! Aggregate subgradient of the objective function.

! Scalar Arguments
    REAL(KIND=dp), INTENT(INOUT) :: & 
         alfv      ! Aggregate locality measure.
    REAL(KIND=dp), INTENT(IN) :: & 
         gamma, &  ! Scaling parameter.
         alfn, &   ! Locality measure.
         rho       ! Correction parameter.
    INTEGER, INTENT(IN) :: & 
         n, &      ! Number of variables
         mc, &     ! Declared number of stored corrections.
         mcc, &    ! Current number of stored corrections.
         inew, &   ! Index for circular arrays.
         icn       ! Correction indicator.
    INTEGER, INTENT(INOUT) :: & 
         iflag     ! Index for adaptive version:
                   !   0  - Maximum number of stored corrections
                   !        has not been changed at this iteration.
                   !   1  - Maximum number of stored corrections
                   !        has been changed at this iteration.
      
! Local arrays
    REAL(KIND=dp), DIMENSION(n) :: tmpn2
    REAL(KIND=dp), DIMENSION(mcc+1) :: tmpmc3, tmpmc4

! Local Scalars
    REAL(KIND=dp) :: &
         pr, &     ! pr = trans(gp-ga) dm (gp-ga), where dm
                   ! presents the L-SR1- approximation of Hessian.
         rrp, &    ! rrp = trans(gp-ga) dm ga - alfv.
         prqr, &   ! prqr = trans(gp-ga) dm (g-ga).
         rrq, &    ! rrq = trans(g-ga) dm ga - alfv + alfn.
         qr, &     ! qr = trans(g-ga) dm (g-ga).
         pq, &     ! pq = trans(g-gp) dm (g-gp).
         qqp, &    ! qqp = trans(g-gp) dm g + alfn.
         lam1, &   ! Multiplier used to calculate aggregate values.
         lam2, &   ! Multiplier used to calculate aggregate values.
         w, &      ! Correction.
         tmp1, &   ! Auxiliary scalar.
         tmp2      ! Auxiliary scalar.
    INTEGER :: i, &
         mcnew, &  ! Current size of vectors.
         iold, &   ! Index of the oldest corrections.
         ierr      ! Error indicador.
    
! Intrinsic Functions
    INTRINSIC MIN,MAX

    ierr = 0
           
    IF (mcc < mc) THEN
       iold = 1
       mcnew = mcc
    ELSE
       IF (iflag == 0) THEN
          mcnew = mc
          iold = inew + 1
          IF (iold > mc+1) iold = 1
       ELSE
          mcnew = mc - 1
          iold = inew + 1
          IF (iold > mc) iold = 1
       END IF
    END IF

      
! Calculation of pq = trans(g-gp) dm (g-gp) = trans(u) dm u.

    pq = vdot(n,u,u)

    IF (icn == 1) THEN
       w = rho * pq
    ELSE
       w = zero
    END IF
    
    IF (mcc > 0) THEN

       IF (iold == 1) THEN
          CALL rwaxv2(n,mcnew,sm,um,u,u,tmpmc3,tmpmc4)
!          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
          CALL trlieq2(mcnew,mcnew,iold,rm,tmpmc3,1,ierr) ! Debugged variant with no compiler warnings
          
          pq = pq - 2.0_dp*vdot(mcnew,tmpmc4,tmpmc3)
          pq = gamma*pq
          
          DO i=1,mcnew
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pq = pq + vdot(mcnew,tmpmc3,tmpmc4)
          
          CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)
          
          pq = pq + gamma*vdot(mcnew,tmpmc3,tmpmc4)

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,u,u,tmpmc3,tmpmc4)
          CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
               u,u,tmpmc3(iold:),tmpmc4(iold:))
!          CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)
          CALL trlieq2(mcnew,mcc,iold,rm,tmpmc3,1,ierr) ! Debugged variant with no compiler warnings

          pq = pq - 2.0_dp*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
               vdot(inew-1,tmpmc4,tmpmc3))
          pq = gamma*pq
          
          DO i=1,mcc
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pq = pq + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
               vdot(inew-1,tmpmc3,tmpmc4)

          CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)
          
          pq = pq + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
               + vdot(inew-1,tmpmc3,tmpmc4))
       END IF

    END IF

    pq = pq + w
      

! Calculation of pr = trans(gp-ga) dm (gp-ga).
      
    CALL xdiffy(n,gp,ga,tmpn2)
    pr = vdot(n,tmpn2,tmpn2)

    IF (icn == 1) THEN
       w = rho * pr
    ELSE
       w = zero
    END IF

    IF (mcc > 0) THEN
       
       IF (iold == 1) THEN 
          DO i=1, mcnew
             tmpmc3(i)=smtgp(i)-smtga(i)
             tmpmc4(i)=umtgp(i)-umtga(i)
          END DO
!          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
          CALL trlieq2(mcnew,mcnew,iold,rm,tmpmc3,1,ierr) ! Debugged variant with no compiler warnings
             
          pr = pr - 2.0_dp*vdot(mcnew,tmpmc4,tmpmc3)
          pr = gamma*pr
            
          DO i=1,mcnew
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pr = pr + vdot(mcnew,tmpmc3,tmpmc4)

          CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)

          pr = pr + gamma*vdot(mcnew,tmpmc3,tmpmc4)

       ELSE
          DO i=1, mcc
             tmpmc3(i)=smtgp(i)-smtga(i)
             tmpmc4(i)=umtgp(i)-umtga(i)
          END DO
!          CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)
          CALL trlieq2(mcnew,mcc,iold,rm,tmpmc3,1,ierr) ! Debugged variant with no compiler warnings

          pr = pr - 2.0_dp*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
               vdot(inew-1,tmpmc4,tmpmc3))
          pr = gamma*pr

          DO  i=1,mcc
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          pr = pr + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
               vdot(inew-1,tmpmc3,tmpmc4)

          CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)

          pr = pr + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
               + vdot(inew-1,tmpmc3,tmpmc4))
       END IF

    END IF

    pr = pr + w

      
! Calculation of rrp = trans(gp-ga) dm ga - alfv.
      
    rrp = - vdot(n,tmpn2,d) - alfv
      

! Calculation of qr = trans(g-ga) dm (g-ga).

    CALL xdiffy(n,g,ga,tmpn2)
    qr = vdot(n,tmpn2,tmpn2)

    IF (icn == 1) THEN
       w = rho * qr
    ELSE
       w = zero
    END IF

    IF (mcc > 0) THEN

       IF (iold == 1) THEN
          CALL rwaxv2(n,mcnew,sm,um,tmpn2,tmpn2,tmpmc3,tmpmc4)
!          CALL trlieq(mcnew,mcnew,iold,rm,tmpmc3,tmpmc3,1,ierr)
          CALL trlieq2(mcnew,mcnew,iold,rm,tmpmc3,1,ierr) ! Debugged variant with no compiler warnings

          qr = qr - 2.0_dp*vdot(mcnew,tmpmc4,tmpmc3)
          qr = gamma*qr
            
          DO i=1,mcnew
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO

          qr = qr + vdot(mcnew,tmpmc3,tmpmc4)
          
          CALL symax(mcnew,mcnew,iold,umtum,tmpmc3,tmpmc4)
          
          qr = qr + gamma*vdot(mcnew,tmpmc3,tmpmc4)

       ELSE
          CALL rwaxv2(n,inew-1,sm,um,tmpn2,tmpn2,tmpmc3,tmpmc4)
          CALL rwaxv2(n,mcc-inew,sm((iold-1)*n+1:),um((iold-1)*n+1:),&
               tmpn2,tmpn2,tmpmc3(iold:),tmpmc4(iold:))
!          CALL trlieq(mcnew,mcc,iold,rm,tmpmc3,tmpmc3,1,ierr)
          CALL trlieq2(mcnew,mcc,iold,rm,tmpmc3,1,ierr) ! Debugged variant with no compiler warnings

          qr = qr - 2.0_dp*(vdot(mcc-inew,tmpmc4(iold:),tmpmc3(iold:)) + &
               vdot(inew-1,tmpmc4,tmpmc3))
          qr = gamma*qr
          
          DO i=1,mcc
             tmpmc4(i) = cm(i)*tmpmc3(i)
          END DO
          
          qr = qr + vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) + &
               vdot(inew-1,tmpmc3,tmpmc4)
          
          CALL symax(mcnew,mcc,iold,umtum,tmpmc3,tmpmc4)
          
          qr = qr + gamma*(vdot(mcc-inew,tmpmc3(iold:),tmpmc4(iold:)) &
               +vdot(inew-1,tmpmc3,tmpmc4))
       END IF
       
    END IF
    
    qr = qr + w
      

! Calculation of rrq = trans(g-ga) dm ga - alfv + alfn.

    rrq = - vdot(n,tmpn2,d) - alfv + alfn

     
! Calculation of prqr = trans(gp-ga) dm (g-ga).
      
    prqr = half*(qr - pq + pr)

     
! Calculation of qqp = trans(g-gp) dm g + alfn.

    qqp = pq + prqr + rrq - pr - rrp

     
! Computation of multipliers lam1 and lam2

    IF (pr > zero .AND. qr > zero) THEN
       tmp1 = rrq/qr
       tmp2 = prqr/qr
       w = pr - prqr*tmp2
       IF (w /= zero) THEN

          lam1 = (tmp1*prqr - rrp)/w
          lam2 = -tmp1 - lam1*tmp2

          IF (lam1*(lam1 - one) < zero .AND. &
               lam2*(lam1 + lam2 - one) < zero) GO TO 200
       END IF
    END IF


! Minimum on the boundary

    lam1 = zero
    lam2 = zero
    IF (alfn <= alfv) lam2 = one
    IF (qr > zero) lam2 = MIN(one,MAX(zero,-rrq/qr))
    w = (lam2*qr + rrq + rrq)*lam2
    tmp1 = zero
    IF (alfv >= zero) tmp1 = one
    IF (pr > zero) tmp1 = MIN(one,MAX(zero,-rrp/pr))
    tmp2 = (tmp1*pr + rrp + rrp)*tmp1
    IF (tmp2 < w) THEN
       w = tmp2
       lam1 = tmp1
       lam2 = zero
    END IF
      
    IF (qqp*(qqp - pq) < zero) THEN
       IF (qr + rrq + rrq - qqp*qqp/pq < w) THEN
          lam1 = qqp/pq
          lam2 = one - lam1
       END IF
    END IF

200 CONTINUE
    IF (lam1 == zero .AND. lam2*(lam2 - one) < zero &
         .AND. -rrp - lam2*prqr > zero .AND. pr > zero) &
         lam1 = MIN(one - lam2, (-rrp-lam2*prqr)/pr)
      

! Computation of the aggregate values
      
    tmp1 = one - lam1 - lam2
    DO i=1,n
       ga(i)=lam1*gp(i)+lam2*g(i)+tmp1*ga(i)
    END DO
    
    alfv = lam2*alfn + tmp1*alfv
      
  END SUBROUTINE agskip

!************************************************************************
!*
!*     * SUBROUTINE wprint *
!*
!*     Printout the warning and error messages.
!*
!************************************************************************
      
  SUBROUTINE wprint(iterm,iiprint,nout) 
    IMPLICIT NONE

! Scalar Arguments
    INTEGER, INTENT(IN) :: &
         iiprint, &    ! Printout specification:
                      !  -1  - No printout.
                      !   0  - Only the error messages.
                      !   1  - The final values of the objective
                      !        function.
                      !   2  - The final values of the objective
                      !        function and the most serious
                      !        warning messages.
                      !   3  - The whole final solution. 
                      !   4  - At each iteration values of the
                      !        objective function.
                      !   5  - At each iteration the whole
                      !        solution
         nout, &      ! Auxilary printout specification.
         iterm        ! Cause of termination:
                      !   1  - The problem has been solved with desired accuracy.
                      !   2  - Changes in function values < tolf in mtesf
                      !        subsequent iterations.
                      !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                      !        where small is the smallest positive number such that 
                      !        1.0 + small > 1.0.
                      !   4  - Number of function calls > mfe.
                      !   5  - Number of iterations > mittt.
                      !   6  - Time limit exceeded. 
                      !   7  - f < tolb.
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.


    IF (iiprint >= 0) THEN

! Initial error messages

       IF (iterm == -5) THEN
          IF (nout == 1) WRITE (6,FMT='(1X,''Error: '' &
               & ''Number of variables (n) is too small, iterm='',I3)') iterm
          IF (nout == 2) WRITE (6,FMT='(1X,''Error: '' &
               &''The maximum number of stored corrections (mcu) '' &
               &''is too small, iterm='',I3)') iterm
          IF (nout == 3) WRITE (6,FMT='(1X,''Error: '' &
               &''The size of the bundle (na) is too small, iterm='' &
               &,I3)') iterm
          IF (nout == 4) WRITE (6,FMT='(1X,''Error: '' &
               &''Line search parameter epsl >= 0.25, iterm='',I3)') iterm
          RETURN
       END IF

        
! Warning messages

       IF (iiprint >= 2) THEN
          IF (iterm == 0) THEN
             IF (nout == -1) WRITE (6,FMT='(1X,''Warning: '' &
                  &''mc > mcu. Assigned mc = mcu.'')')
             IF (nout == -2) WRITE (6,FMT='(1X,''Warning: '' &
                  &''A line search parameter epsr >= 0.5.'')')
             IF (nout == -3) WRITE (6,FMT='(1X,''Warning: '' &
                  &''A nondescent search direction occured. Restart.'')')
             IF (nout == -4) WRITE (6,FMT='(1X,''Warning: '' &
                  &''Does not converge.'')')
             IF (nout == -5) WRITE (6,FMT='(1X,''Warning: '' &
                  &''tmax < tmin. Restart.'')')
             RETURN
          END IF
         

! Printout the final results
            
          IF (iterm == 6) WRITE (6,FMT='(1X,''Abnormal exit: Time is up.'')')
          IF (iterm == 7) WRITE (6,FMT='(1X,''Abnormal exit: f < tolb.'')')
          IF (iterm == 2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''Too many steps without significant progress.'')')
          IF (iterm == 3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''The value of the function does not change.'')')
          IF (iterm == 5) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''Number of iterations > '',I5)') nout
          IF (iterm == 4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''Number of function evaluations > '',I5)') nout
          IF (iterm == -1) THEN
             IF (nout == -1) THEN
                WRITE (6,FMT='(1X,''Abnormal exit: Two consecutive restarts.'')')
             ELSE
                WRITE (6,FMT='(1X,''Abnormal exit: tmax < tmin in two subsequent iterations.'')')
             END IF
          END IF
          IF (iterm == -2) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''Number of restarts > '',I5''.'')') nout
          IF (iterm == -3) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''Failure in function or subgradient calculations.'')')
          IF (iterm == -4) WRITE (6,FMT='(1X,''Abnormal exit: '' &
               &''Failure in attaining the demanded accuracy.'')')
       END IF
    END IF
      
  END SUBROUTINE wprint

            
!************************************************************************
!*
!*     * SUBROUTINE rprint *
!*      
!*     Printout the (final) results.
!*
!************************************************************************
      
  SUBROUTINE rprint(n,nit,nfe,nge,x,f,wk,qk,iterm,iiprint) 
    IMPLICIT NONE

! Scalar Arguments
    INTEGER, INTENT(IN) :: & 
         n, &         ! Number of variables 
         nit, &       ! Number of used iterations.
         nfe, &       ! Number of used function evaluations.
         nge, &       ! Number of used subgradient evaluations.
         iiprint, &    ! Printout specification:
                      !  -1  - No printout.
                      !   0  - Only the error messages.
                      !   1  - The final values of the objective
                      !        function.
                      !   2  - The final values of the objective
                      !        function and the most serious
                      !        warning messages.
                      !   3  - The whole final solution. 
                      !   4  - At each iteration values of the
                      !        objective function.
                      !   5  - At each iteration the whole
                      !        solution
         iterm        ! Cause of termination:
                      !   1  - The problem has been solved with desired accuracy.
                      !   2  - Changes in function values < tolf in mtesf
                      !        subsequent iterations.
                      !   3  - Changes in function value < tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                      !        where small is the smallest positive number such that 
                      !        1.0 + small > 1.0.
                      !   4  - Number of function calls > mfe.
                      !   5  - Number of iterations > mittt.
                      !   6  - Time limit exceeded. 
                      !   7  - f < tolb.
                      !  -1  - Two consecutive restarts.
                      !  -2  - Number of restarts > maximum number of restarts.
                      !  -3  - Failure in function or subgradient calculations 
                      !        (assigned by the user).
                      !  -4  - Failure in attaining the demanded accuracy.
                      !  -5  - Invalid input parameters.


    REAL(KIND=dp), INTENT(IN) :: &
         f, &         ! Value of the objective function.
         wk, &        ! Value of the first stopping criterion.
         qk           ! Value of the second stopping criterion.

! Array Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: &
         x            ! Vector of variables
         
! Local Scalars
    INTEGER :: i

! Intermediate results
    
    IF (iterm == 0) THEN
       IF (iiprint > 3) WRITE (6,FMT='(1X,''nit='',I5,2X, &
            &''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
            &''qk='',D11.4,2X)') nit,nfe,nge,f,wk,qk
       IF (iiprint == 5) WRITE (6,FMT='(1X,''x='', &
            &5D15.7:/(4X,5D15.7))')(x(i),i=1,n)
       RETURN
    END IF
         

! Final results

    IF (iiprint > 0) WRITE (6,FMT='(1X,''nit='',I5,2X, &
         &''nfe='',I5,2X,''nge='',I5,2X,''f='',D15.8,2X,''wk='',D11.4,2X, &
         &''qk='',D11.4,2X,''iterm='',I3)') nit,nfe,nge,f,wk,qk,iterm
    IF (iiprint .EQ. 3 .OR. iiprint .EQ. 5) &
         WRITE (6,FMT='(1X,''x='',5D15.7:/(4X,5D15.7))')(x(i),i=1,n)
      
  END SUBROUTINE rprint
      
END MODULE lmbm_mod   !LMBM ENDS

!=======================================================================================
!
!                     *************  CLUST-SPLITTER  *************
!
!=======================================================================================

MODULE our_method

    USE constants
    USE functions
    USE param
    USE exe_time
    USE lmbm_sub
    USE initializat
    USE obj_fun
    USE lmbm_mod
    
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE clustsplitter(infile, usedmethod, max_cl, delete_outlier, n_outlier, opt1, opt2, opt3,&
            noopt1, noopt2, noopt3, noopt4, noopt5, noopt6, optimize_startingpoint, ncenter1, ncenter2, all_start,&
            index1, index2, features, records, x, f, f_method4, db, dn, nf, nsub, cpu, datareadingtime)
            
            IMPLICIT NONE
        
            !---------------------
            ! Arguments
            !---------------------
            CHARACTER(LEN=*), INTENT(IN) :: infile                  ! Infile data, for example 'iris.txt'
            INTEGER, INTENT(IN) :: features                         ! Number of features in the infile data                                 
            INTEGER, INTENT(IN) :: records                          ! Number of records in the infile data  
            INTEGER, INTENT(IN) :: usedmethod                       ! Method to use
                                                                        ! 4 = Clust-Splitter
                                                                        ! 2 = Clust-Splitter without k-clustering problem
            INTEGER, INTENT(INOUT) :: max_cl                        ! Maximum number of clusters   
            INTEGER, INTENT(IN) :: delete_outlier                   ! Do you want to delete outlier points from the data
                                                                        ! 0=yes, 1=no
            INTEGER, INTENT(INOUT) :: n_outlier                     ! Limit for outlier points  
            ! Starting points for starting point auxiliary problem when we optimize starting point (optimize_startingpoint = 1)
            INTEGER, INTENT(IN) :: opt1                             ! Average center of random data points (when the distance of the average point from the cluster center does not matter)
            INTEGER, INTENT(IN) :: opt2                             ! Average center of random data points (when the distance of the average point from the cluster center matters)
            INTEGER, INTENT(IN) :: opt3                             ! The cluster center
            ! Starting points for 2-clustering problem when we do not optimize starting point (optimize_startingpoint = 0)
            INTEGER, INTENT(IN) :: noopt1                           ! The cluster center + random data point
            INTEGER, INTENT(IN) :: noopt2                           ! The cluster center + the average of random data points
            INTEGER, INTENT(IN) :: noopt3                           ! The cluster center + The cluster center
            INTEGER, INTENT(IN) :: noopt4                           ! Random data point + random data point
            INTEGER, INTENT(IN) :: noopt5                           ! Random data point + the average of random data points
            INTEGER, INTENT(IN) :: noopt6                           ! The average of random data points + the average of random data points
            INTEGER, INTENT(IN) :: all_start                        ! Calculate 2-clustering auxiliry problem from all starting points or just from the best
                                                                        ! 1=best, 0=all
            INTEGER, INTENT(IN) :: index1                           ! Do you want to make DBI indices better using 2-clustering auxiliary problem (at the cost of the function value)
                                                                        ! 1=yes 0=no
            INTEGER, INTENT(IN) :: index2                           ! Do you want to make DBI indices better using k-clustering problem (at the cost of the function value)
                                                                        ! 1=yes 0=no
            INTEGER, INTENT(IN) :: optimize_startingpoint           ! Do you want to optimize starting points using the starting point auxiliary problem
                                                                        ! 1=yes, 0=no
            INTEGER, INTENT(IN) :: ncenter1, ncenter2               ! Number of random points used to calculate the average center in auxiliary problems 
            REAL(KIND=dp), INTENT(INOUT) :: datareadingtime         ! Time for reading the data in
            REAL(KIND=dp), DIMENSION(features*max_cl,max_cl), INTENT(OUT) :: x  ! Cluster centers
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(OUT) :: f                  ! Function values for clusters of different sizes  
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(OUT) :: f_method4          ! Function values for clusters of different sizes for method 4 when we want to delete outliers
                                                                                    ! f_method4 differs from f only when delete_outlier=0
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(INOUT) :: nf               ! Number of function evaluations
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(INOUT) :: nsub             ! Number of subgradient evaluations 
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(INOUT) :: cpu              ! CPU times (without data reading time)
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(INOUT) :: db               ! Davies-Bouldin indexes (DBI)
            REAL(KIND=dp), DIMENSION(max_cl), INTENT(INOUT) :: dn               ! Dunn indexes (DI)

            !------------------------
            ! Local variables
            !------------------------
            
            ! Real variables
            REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: xpoint              ! Cluster centers
            REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: new_xpoint          ! Help variable for cluster centers
            REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: newpoint            ! Help variable for cluster centers
            REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: bestpoint           ! Best result so far
            REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f_table             ! Table for function values
            REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: x_table           ! Table for points
            REAL(KIND=dp), DIMENSION(features * max_cl) :: y                ! Variable where intermediate points are stored
            REAL(KIND=dp), DIMENSION(features * max_cl) :: help_y           ! Help variable where intermediate points are stored
            REAL(KIND=dp), DIMENSION(features) :: centroid                  ! Cluster center
            REAL(KIND=dp), DIMENSION(features,ncenter1) :: center_table1    ! Table for storing data points for which we want to calculate the centroid
            REAL(KIND=dp), DIMENSION(features,ncenter2) :: center_table2    ! Table for storing data points for which we want to calculate the centroid
            REAL(KIND=dp), DIMENSION(features) :: helppoint                 ! Help point
            REAL(KIND=dp), DIMENSION(max_cl) :: f_a_sep                     ! Function value of each cluster
            REAL(KIND=dp), DIMENSION(features) :: old_center                ! Old cluster center
            REAL(KIND=dp) :: LMBMstart                                      ! The starting time in LMBM
            REAL(KIND=dp) :: LMBMend                                        ! The ending time in LMBM
            REAL(KIND=dp) :: s_time, f_time                                 ! Times for calculating CPU times
            REAL(KIND=dp) :: first_time, second_time                        ! Times for calculating data reading time
            REAL(KIND=dp) :: davbou                                         ! Davies-Bouldin index (DBI)
            REAL(KIND=dp) :: dunn                                           ! Dunn index (DI)
            REAL(KIND=dp) :: f_subp                                         ! The objective function value at the solution of the reduced problem
            REAL(KIND=dp) :: f_a                                            ! Function value
            REAL(KIND=dp) :: f_a_method4                                    ! Function value in method 4
            REAL(KIND=dp) :: f_smallest                                     ! Help variable, smallest function value so far
            REAL(KIND=dp) :: f_smallest_new                                 ! Help variable, smallest function value so far
            REAL(KIND=dp) :: beta                                           ! Random number between 0 and 1
            REAL(KIND=dp) :: largestvalue                                   ! Help variable, large value
            REAL(KIND=dp) :: element_sum                                    ! Variable for summing up elements
            REAL(KIND=dp) :: oldcoef                                        ! Old coefficient
            REAL(KIND=dp) :: value1                                         ! Help variable for storing value
            REAL(KIND=dp) :: smallest_value                                 ! Help variable, small value
            REAL(KIND=dp) :: difference                                     ! Difference between two points
            REAL(KIND=dp) :: radius                                         ! Cluster radius
            
            ! Integer variables
            INTEGER, ALLOCATABLE, DIMENSION(:) :: predicted_labels      ! Vector of predicted labels
            INTEGER, DIMENSION(max_cl) :: n_a_sep                       ! Number of points in each cluster
            INTEGER, DIMENSION(max_cl) :: n_a_sep_help                  ! Help variable, number of points in each cluster
            INTEGER, DIMENSION(max_cl,records) :: indexes               ! Help variable, index of cluster that datapoint belongs to
            INTEGER, DIMENSION(4) :: iout                               ! Output integer parameters.
                                                                         !   iout(1)   Number of used iterations.
                                                                         !   iout(2)   Number of used function evaluations.
                                                                         !   iout(3)   Number of used subgradient evaluations (LMBM)
                                                                         !               or number of outer iterations (LDGBM).
                                                                         !   iout(4)   Cause of termination:
                                                                         !               1  - The problem has been solved
                                                                         !                    with desired accuracy.
                                                                         !               2  - Changes in function values < tolf in mtesf
                                                                         !                    subsequent iterations.
                                                                         !               3  - Changes in function value < 
                                                                         !                    tolf*small*MAX(|f_k|,|f_(k-1)|,1),
                                                                         !                    where small is the smallest positive number
                                                                         !                    such that 1.0 + small > 1.0.
                                                                         !               4  - Number of function calls > mfe.
                                                                         !               5  - Number of iterations > mit.
                                                                         !               6  - Time limit exceeded. 
                                                                         !               7  - f < tolb.
                                                                         !              -1  - Two consecutive restarts.
                                                                         !              -2  - Number of restarts > maximum number
                                                                         !                    of restarts.
                                                                         !              -3  - Failure in function or subgradient
                                                                         !                    calculations (assigned by the user).
                                                                         !              -4  - Failure in attaining the demanded
                                                                         !                    accuracy.
                                                                         !              -5  - Invalid input parameters.
                                                                         !              -6  - Unspecified error.
                                                                         
            INTEGER :: nproblem             ! 1 = whole data a, 2 = subdata b
            INTEGER :: kl                   ! Number of clusters, iteration counter
            INTEGER :: i, j, k              ! Help variables
            INTEGER :: l_limit, u_limit     ! Help variables, lower and upper limits
            INTEGER :: largest_ind          ! Largest index so far
            INTEGER :: y_dim                ! Dimension of y
            INTEGER :: old_largest_ind      ! Largest index of previous iteration round
            INTEGER :: mc                   ! Number of corrections (for lmbm)
            INTEGER :: continuing           ! Help variable to check if we should continue or not
            INTEGER :: loop_counter         ! Loop counter
            INTEGER :: action_taken         ! Variable for checking if we did something
            INTEGER :: are_outliers         ! Variable for checking if there are outliers or not
            INTEGER :: change_ind           ! Index we are modifying
            INTEGER :: n_outlier_2          ! Help variable for storing the number of outliers
            INTEGER :: beta2                ! Random integer number
            INTEGER :: nrandom              ! Loop counter
            INTEGER :: smallest_index       ! Smallest index so far
            INTEGER :: n_randompoints       ! Total number of random points
            INTEGER :: neg_fvalue           ! Is the function value negative; 0=no, 1=yes
            INTEGER :: divisor              ! Divisor
            
            ! Logical variables
            LOGICAL :: should_continue     ! Help variable to check if we should continue of not
            
            !===================================
            ! Beginning of the actual program
            !===================================
            
            CALL random_seed
            
            mc = 7    ! Initialize parameter mc
            
            nf = 0.0_dp         ! Initialize the number of function evaluations to 0
            nsub = 0.0_dp       ! Initialize the number of subgradient evaluations to 0
            cpu = 0.0_dp        ! Initialize all CPU times to 0
            
            ! Time for reading the data
            CALL cpu_time(first_time)
            CALL allocate_whole_data(infile, records, features)  ! Read the data in and initialize all parameters
            CALL cpu_time(second_time)
            datareadingtime = second_time-first_time
            PRINT*, 'DATA READING: ', datareadingtime
            
            CALL cpu_time(s_time)  ! Time for starting the computations
            
            ! Initialize variables x and y to zero
            y = 0.0_dp
            x = 0.0_dp
            
            ! Loop for cluster size 1 and 2
            DO kl = 1,2
            
                ! Check if we want to make DBI indices better (at the cost of function value)
                IF (kl==1) THEN
                    CALL consider_indices(0, 0) ! When kl=1, we do not consider indices
                ELSE
                    CALL consider_indices(index1, index2)
                END IF
            
                PRINT*, 'NUMBER OF CLUSTERS ', kl
                
                ! Initialize the number of starting points
                IF (kl == 1) THEN
                    n_randompoints = 1  ! When kl=1, the problem is convex so only one starting point is enough
                ELSE IF (optimize_startingpoint == 1) THEN
                    n_randompoints = opt1+opt2+opt3  ! We want to optimize starting points using starting point auxiliary problem
                ELSE
                    n_randompoints = noopt1+noopt2+noopt3+noopt4+noopt5+noopt6  ! We do not want to optimize starting points
                END IF
            
                ! Allocate the starting point for LMBM and other necessary variables
                ALLOCATE(xpoint(kl*features))
                ALLOCATE(bestpoint(kl*features))
                ALLOCATE(f_table(n_randompoints))
                ALLOCATE(x_table(features*kl,n_randompoints))
                
                ! Initializations
                xpoint = 0.0_dp
                f_smallest = 3.40282347*10.**38
                f_smallest_new = 3.40282347*10.**38
                smallest_index = 1
                
                DO nrandom = 1,n_randompoints
                    
                    IF (kl==1) THEN  ! This loop when kl=1
                            
                        ! Take amount of 'ncenter1' random data points
                        DO j=1,ncenter1
                            CALL random_number(beta)
                            ! Check that random number was not exactly zero
                            DO WHILE (beta==0)
                                CALL random_number(beta)
                            END DO
                            beta2 = ceiling(records*beta)
                            ! Add random data point into table
                            DO k=1,features
                                center_table1(k,j)=a(k,beta2)
                            END DO
                        END DO
                        ! Calculate the average of the random data points
                        DO j=1,features
                            element_sum = 0.0_dp
                            DO k=1,ncenter1
                                element_sum = element_sum + center_table1(j,k)
                            END DO
                            ! Calculate the average (i.e. center)
                            centroid(j) = element_sum / ncenter1
                        END DO
                        DO i=1,features
                            ! Starting point is the calculated center
                            xpoint(i)=centroid(i)
                        END DO

                    ELSE  ! This loop when kl=2
                    
                        !-------------------------------------
                        ! Starting point auxiliary problem
                        !-------------------------------------
                        IF (optimize_startingpoint == 1) THEN ! This loop when kl=2 and we want to optimize starting point
                            
                            ! Starting point is the center of random data points
                            IF (nrandom<=opt1) THEN
                                DO j=1,ncenter1
                                    CALL random_number(beta)
                                    DO WHILE (beta==0)
                                        CALL random_number(beta)
                                    END DO
                                    beta2 = ceiling(nrecords_a*beta)
                                    DO k=1,features
                                        center_table1(k,j)=a(k,beta2)
                                    END DO
                                END DO
                                DO j=1,features
                                    element_sum = 0.0_dp
                                    DO k=1,ncenter1
                                        element_sum = element_sum + center_table1(j,k)
                                    END DO
                                    centroid(j) = element_sum / ncenter1
                                END DO
                                DO i=1,features
                                    xpoint(i)=centroid(i) ! Save starting point
                                END DO
                                
                            ! Starting point is the center of random data points
                            ! We also ensure that the calculated center is sufficiently far from the original cluster center
                            ELSE IF (nrandom>opt1 .AND. nrandom /= n_randompoints) THEN
                                difference = 0.0_dp
                                DO i=1,features
                                    old_center(i)=y(i) ! Save the old cluster center
                                END DO
                                CALL calculate_radius(old_center, radius)
                                divisor = 2
                                DO WHILE (difference < radius/divisor) ! Check the difference between new center and old center
                                    DO j=1,ncenter2
                                        CALL random_number(beta)
                                        DO WHILE (beta==0)
                                            CALL random_number(beta)
                                        END DO
                                        beta2 = ceiling(nrecords_a*beta)
                                        DO k=1,features
                                            center_table2(k,j)=a(k,beta2)
                                        END DO
                                    END DO
                                    DO j=1,features
                                        element_sum = 0.0_dp
                                        DO k=1,ncenter2
                                            element_sum = element_sum + center_table2(j,k)
                                        END DO
                                        centroid(j) = element_sum / ncenter2
                                    END DO
                                    difference = 0.0_dp
                                    DO i=1,features
                                        difference = difference + (old_center(i)-centroid(i))**2
                                    END DO
                                    divisor = divisor * 2
                                END DO ! Do While ends
                                DO i=1,features
                                    xpoint(i)=centroid(i) ! Save starting point
                                END DO
                                
                            ! Starting point is the cluster center
                            ELSE
                                DO i=1,features
                                    xpoint(i)=y(i) ! Save starting point
                                END DO
                            END IF
                            
                            ! Calculate errors
                            DO i=1,features
                                helppoint(i)=y(i)
                            END DO
                            CALL set_errors_a(helppoint)
                            
                            ! Solve starting point auxiliary problem with LMBM
                            nproblem = 3
                            CALL set_nclust(1)
                            ! Calls for LMBM
                            CALL allocate_x_var(1*features)
                            CALL init_x_var(xpoint)
                            CALL init_problem(nproblem)
                            CALL init_lmbmpar(1)                    
                            CALL cpu_time(LMBMstart)   ! Start CPU timining
                            CALL lmbm(mc,f_subp,iout(1),iout(2),iout(3),iout(4),LMBMstart)
                            CALL copy_x_var(helppoint)
                            CALL deallocate_x_var()
                            
                            ! Update the counts of function and subgradient evaluations
                            DO i=kl,max_cl
                                nf(i) = nf(i)+iout(2)
                                nsub(i) = nsub(i)+iout(3)
                            END DO
                            
                            ! Save the old cluster center
                            DO i=1,features
                                xpoint(i)=y(i)
                            END DO
                            ! Copy the result of starting point auxiliary problem
                            DO i=features+1,2*features
                                xpoint(i)=helppoint(i-features)
                            END DO
                            
                            ! Check if the result point gave smaller function value than previous points
                            IF (f_subp < f_smallest) THEN
                                f_smallest = f_subp
                                DO i=1,2*features
                                    bestpoint(i)=xpoint(i)  ! Save the best starting points
                                END DO
                            END IF
                            
                        !-----------------------------------------
                        ! Starting point auxiliary problem ends
                        !-----------------------------------------
                            
                        ELSE ! This loop when kl=2 and we do not want to optimize starting point
                        
                            ! Starting points: cluster center + random data point
                            IF (nrandom <= noopt1) THEN
                                DO i=1,features
                                    xpoint(i)=y(i) ! Save the first starting point
                                END DO
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_a*beta)
                                DO i=features+1,2*features
                                    xpoint(i)=a(i-features,beta2) ! Save the second starting point
                                END DO
                                
                            ! Starting points: cluster center + center of random data points
                            ELSE IF (nrandom <= noopt1+noopt2) THEN
                                DO i=1,features
                                    xpoint(i)=y(i) ! Save the first starting point
                                END DO
                                DO j=1,ncenter1
                                    CALL random_number(beta)
                                    DO WHILE (beta==0)
                                        CALL random_number(beta)
                                    END DO
                                    beta2 = ceiling(nrecords_a*beta)
                                    DO k=1,features
                                        center_table1(k,j)=a(k,beta2)
                                    END DO
                                END DO
                                DO j=1,features
                                    element_sum = 0.0_dp
                                    DO k=1,ncenter1
                                        element_sum = element_sum + center_table1(j,k)
                                    END DO
                                    centroid(j) = element_sum / ncenter1
                                END DO
                                DO i=features+1,2*features
                                    xpoint(i)=centroid(i-features) ! Save the second starting point
                                END DO
                               
                            ! Starting points: cluster center + cluster center
                            ELSE IF (nrandom <= noopt1+noopt2+noopt3) THEN
                                DO i=1,features
                                    xpoint(i)=y(i) ! Save the first starting point
                                END DO
                                DO i=features+1,2*features
                                    xpoint(i)=y(i-features) ! Save the second starting point
                                END DO
                            
                            ! Starting points: random data point + random data point
                            ELSE IF (nrandom <= noopt1+noopt2+noopt3+noopt4) THEN
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_a*beta)
                                DO i=1,features
                                    xpoint(i)=a(i,beta2) ! Save the first starting point
                                END DO
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_a*beta)
                                DO i=features+1,2*features
                                    xpoint(i)=a(i-features,beta2) ! Save the second starting point
                                END DO
                            
                            ! Starting points: random data point + center of random data points
                            ELSE IF (nrandom <= noopt1+noopt2+noopt3+noopt4+noopt5) THEN
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_a*beta)
                                DO i=1,features
                                    xpoint(i)=a(i,beta2) ! Save the first starting point
                                END DO
                                DO j=1,ncenter1
                                    CALL random_number(beta)
                                    DO WHILE (beta==0)
                                        CALL random_number(beta)
                                    END DO
                                    beta2 = ceiling(nrecords_a*beta)
                                    DO k=1,features
                                        center_table1(k,j)=a(k,beta2)
                                    END DO
                                END DO
                                DO j=1,features
                                    element_sum = 0.0_dp
                                    DO k=1,ncenter1
                                        element_sum = element_sum + center_table1(j,k)
                                    END DO
                                    centroid(j) = element_sum / ncenter1
                                END DO
                                DO i=features+1,2*features
                                    xpoint(i)=centroid(i-features) ! Save the second starting point
                                END DO
                            
                            ! Starting points: center of random data points + center of random data points
                            ELSE
                                DO j=1,ncenter1
                                    CALL random_number(beta)
                                    DO WHILE (beta==0)
                                        CALL random_number(beta)
                                    END DO
                                    beta2 = ceiling(nrecords_a*beta)
                                    DO k=1,features
                                        center_table1(k,j)=a(k,beta2)
                                    END DO
                                END DO
                                DO j=1,features
                                    element_sum = 0.0_dp
                                    DO k=1,ncenter1
                                        element_sum = element_sum + center_table1(j,k)
                                    END DO
                                    centroid(j) = element_sum / ncenter1
                                END DO
                                DO i=1,features
                                    xpoint(i)=centroid(i) ! Save the first starting point
                                END DO
                                DO j=1,ncenter1
                                    CALL random_number(beta)
                                    DO WHILE (beta==0)
                                        CALL random_number(beta)
                                    END DO
                                    beta2 = ceiling(nrecords_a*beta)
                                    DO k=1,features
                                        center_table1(k,j)=a(k,beta2)
                                    END DO
                                END DO
                                DO j=1,features
                                    element_sum = 0.0_dp
                                    DO k=1,ncenter1
                                        element_sum = element_sum + center_table1(j,k)
                                    END DO
                                    centroid(j) = element_sum / ncenter1
                                END DO
                                DO i=features+1,2*features
                                    xpoint(i)=centroid(i-features) ! Save the second starting point
                                END DO
                            END IF
                                
                        END IF
                    END IF
                    
                    ! This loop when:
                        ! we want to solve 2-clustering auxiliary problem from all starting points OR
                        ! we want to solve 2-clustering auxiliary problem only from the best point and we already know that point
                    IF (all_start == 0 .OR. nrandom == n_randompoints) THEN
                        
                        ! This loop when kl=2 and we want to solve 2-clustering auxiliary problem only from the best starting point
                        IF (kl==2 .AND. all_start == 1) THEN
                            DO i=1,kl*features
                                xpoint(i) = bestpoint(i) ! Set the best starting point into variable xpoint
                            END DO
                        END IF
                        
                        nproblem = 1
                        CALL set_nclust(kl) ! Set the number of clusters
                        
                        ! Allocate new_xpoint
                        ALLOCATE(new_xpoint(kl*features))
                        
                        ! Initializations
                        neg_fvalue = 0
                        oldcoef = 1
                        
                        444 CONTINUE
                        
                        ! Check if we want to calculate better DBI indices (at the cost of function value)
                        IF (kl==2 .AND. index2 == 1) THEN
                            CALL set_coefficient2(xpoint, neg_fvalue, oldcoef)
                        END IF
                        
                        !----------------------------------
                        ! 2-clustering auxiliary problem
                        !----------------------------------
                        
                        ! Calls for LMBM
                        CALL allocate_x_var(kl*features)
                        CALL init_x_var(xpoint)
                        CALL init_problem(nproblem)
                        CALL init_lmbmpar(kl)                    
                        CALL cpu_time(LMBMstart) ! Start CPU timining
                        CALL lmbm(mc,f_subp,iout(1),iout(2),iout(3),iout(4),LMBMstart)
                        CALL copy_x_var(new_xpoint)
                        CALL deallocate_x_var()
                        
                        ! Update the counts of function and subgradient evaluations
                        DO i=kl,max_cl
                            nf(i) = nf(i)+iout(2)
                            nsub(i) = nsub(i)+iout(3)
                        END DO
                        
                        neg_fvalue = 1
                            
                        IF (f_subp < 0.0_dp) GO TO 444  ! We do not want negative f_subp value
                        
                        ! Calculate f value (without indexes)
                        IF (kl==2 .AND. index2 == 1) THEN
                            CALL func_a(new_xpoint, kl*features, f_subp)
                        END IF
                        
                        xpoint = new_xpoint ! Save the result point into variable xpoint

                        ! Save result values (the function value and the cluster centers)
                        f_table(nrandom) = f_subp
                        DO i=1,features*kl
                            x_table(i,nrandom) = xpoint(i)
                        END DO
                        
                        ! Save smallest function value and smallest index
                        IF (f_subp < f_smallest_new) THEN
                            f_smallest_new = f_subp
                            smallest_index = nrandom
                        END IF
                        
                    END IF
                    
                END DO
                
                ! Update the best solution (cluster centers)
                DO i=1,features*kl
                    y(i)=x_table(i,smallest_index)
                    x(i,kl)=x_table(i,smallest_index)
                END DO
                
                ! Store the obtained function value
                f(kl)=f_smallest_new  ! The "normal" function value
                f_method4(kl)=f_smallest_new  ! The function value for data with outlier points removed
                
                ! Calculate validity indices
                CALL check(kl,y,davbou,dunn)
                db(kl) = davbou ! Davies-Bouldin index (DBI)
                dn(kl) = dunn ! Dunn index (DI)
                
                ! Calculate predicted labels if max_cl=2
                IF (kl == max_cl) THEN
                    ALLOCATE(predicted_labels(nrecords_a))
                    predicted_labels = 100 ! Initialize predicted labels
                    DO i=1,nrecords_a ! Iterate through the data points
                        smallest_value = 3.40282347*10.**38
                        DO j=1,kl ! Iterate through the clusters
                            value1 = 0.0_dp
                            DO k=1,nfeatures_a ! Iterate through the features
                                value1 = value1 + (y(k+nfeatures_a*(j-1))-a(k,i))**2
                            END DO
                            IF (value1 < smallest_value) THEN
                                smallest_value = value1
                                predicted_labels(i) = j-1 ! Determine which cluster the data point belongs to
                            END IF
                        END DO
                    END DO
                    ! Print predicted labels
                    !PRINT*, 'Predicted labels:'
                    !PRINT*, predicted_labels
                    DEALLOCATE(predicted_labels)
                END IF
                
                ! Deallocations
                DEALLOCATE(xpoint)
                DEALLOCATE(new_xpoint)
                DEALLOCATE(bestpoint)
                DEALLOCATE(f_table)
                DEALLOCATE(x_table)
                
                ! Finish time of the iteration
                CALL cpu_time(f_time)
                
                !PRINT*, 'TIME', f_time-s_time
                !PRINT*, ''
                
                ! The CPU of the iteration
                cpu(kl) = f_time-s_time
                
            END DO
            
            ! Loop for cluster size > 2
            DO kl = 3,max_cl
            
                PRINT*, 'NUMBER OF CLUSTERS ', kl
                
                ! Initializations
                loop_counter = 0
                old_largest_ind = 1
                y_dim = features * max_cl
                
                ! Find the cluster to be split next
                CALL func_a_sep(y, y_dim, kl-1, max_cl, f_a_sep, n_a_sep, indexes)
                
                ! When k=3, determine the index of the cluster with the highest function value
                IF (kl == 3) THEN
                    largestvalue = -3.40282347*10.**38
                    DO i=1,2
                        IF (f_a_sep(i) > largestvalue) THEN
                            old_largest_ind = i
                            largestvalue = f_a_sep(i)
                        END IF
                    END DO
                END IF
                
                ! Check if there are outlier clusters (clusters containing less than n_outlier points)
                are_outliers = 0
                DO i=1,kl-1
                    IF (n_a_sep(i) < n_outlier .AND. n_a_sep(i) > 0) THEN
                        are_outliers = 1 ! Cluster i is outlier cluster
                    END IF
                END DO
                
                ! Initializations
                action_taken = 0
                largestvalue = -3.40282347*10.**38
                
                IF (are_outliers==0) THEN ! There are no outlier clusters
                    DO i = 1,kl-1
                        IF (f_a_sep(i)>largestvalue) THEN
                            largestvalue = f_a_sep(i)
                            largest_ind = i ! Check the index of the cluster with the highest function value
                            action_taken = 1 ! There are non-outlier clusters
                        END IF
                    END DO
                ELSE  ! There are outlier cluster(s)
                    DO i = 1,kl-1
                        ! Check the largest function value only within non-outlier clusters
                        IF (n_a_sep(i) > n_outlier-1 .AND. i /= old_largest_ind) THEN
                            IF (f_a_sep(i)>largestvalue) THEN
                                largestvalue = f_a_sep(i)
                                largest_ind = i ! Store the index of the cluster with the highest function value
                                action_taken = 1 ! There is non-outlier cluster(s)
                            END IF
                        END IF
                    END DO
                END IF
                
                ! All the clusters were outlier clusters
                IF (action_taken == 0) THEN
                    largest_ind = old_largest_ind ! Examine the same cluster as in the previous iteration
                END IF
                
                ! Store the index of the cluster with highest function value for the next iteration
                old_largest_ind = largest_ind
                
                ! Form the subdata b containing the points of the cluster under consideration
                CALL set_submatrix(largest_ind)
                
                ! Initializations
                f_smallest = 3.40282347*10.**38
                f_smallest_new = 3.40282347*10.**38
                
                ! Allocations
                ALLOCATE(xpoint(2*features))
                ALLOCATE(newpoint(2*features))
                ALLOCATE(bestpoint(2*features))
                ALLOCATE(f_table(n_randompoints))
                ALLOCATE(x_table(features*kl,n_randompoints))
                
                ! Store the cluster centers
                DO i=1,n_randompoints
                    DO j=1,features*(kl-1)
                        x_table(j,i)=y(j)
                    END DO
                END DO
                
                ! Loop over the number of starting points
                DO nrandom = 1,n_randompoints
                
                    !-----------------------------------
                    ! Starting point auxiliary problem
                    !-----------------------------------
                    IF (optimize_startingpoint == 1) THEN ! This loop when we want to optimize starting point
                        
                        ! Starting point is the center of random data points
                        IF (nrandom<=opt1) THEN
                            DO j=1,ncenter1
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_b*beta)
                                DO k=1,features
                                    center_table1(k,j)=b(k,beta2)
                                END DO
                            END DO
                            DO j=1,features
                                element_sum = 0.0_dp
                                DO k=1,ncenter1
                                    element_sum = element_sum + center_table1(j,k)
                                END DO
                                centroid(j) = element_sum / ncenter1
                            END DO
                            DO i=1,features
                                xpoint(i)=centroid(i) ! Save the starting point
                            END DO
                            
                        ! Starting point is the center of random data points
                        ! We also ensure that the calculated center is sufficiently far from the original cluster center
                        ELSE IF (nrandom>opt1 .AND. nrandom /= n_randompoints) THEN
                            difference = 0.0_dp
                            DO i=1,features
                                old_center(i)=y(i+features*(largest_ind-1))
                            END DO
                            CALL calculate_radius2(old_center, radius)
                            divisor = 2
                            DO WHILE (difference < radius/divisor) ! Check the difference between new center and old center
                                DO j=1,ncenter2
                                    CALL random_number(beta)
                                    DO WHILE (beta==0)
                                        CALL random_number(beta)
                                    END DO
                                    beta2 = ceiling(nrecords_b*beta)
                                    DO k=1,features
                                        center_table2(k,j)=b(k,beta2)
                                    END DO
                                END DO
                                DO j=1,features
                                    element_sum = 0.0_dp
                                    DO k=1,ncenter2
                                        element_sum = element_sum + center_table2(j,k)
                                    END DO
                                    centroid(j) = element_sum / ncenter2
                                END DO
                                difference = 0.0_dp
                                DO i=1,features
                                    difference = difference + (old_center(i)-centroid(i))**2
                                END DO
                                divisor = divisor * 2
                            END DO ! Do while ends  
                            DO i=1,features
                                xpoint(i)=centroid(i) ! Save the starting point
                            END DO
                            
                        ! Starting point is the cluster center
                        ELSE
                            DO i=1,features
                                xpoint(i)=y(i+features*(largest_ind-1)) ! Save the starting point
                            END DO
                        END IF
                        
                        ! Calculate errors
                        DO i=1,features
                            helppoint(i)=y(i+features*(largest_ind-1))
                        END DO
                        CALL set_errors_b(helppoint)
                        
                        ! Solve starting point auxiliary problem with LMBM
                        nproblem = 4
                        CALL set_nclust(1)
                        ! Calls for LMBM
                        CALL allocate_x_var(1*features)
                        CALL init_x_var(xpoint)
                        CALL init_problem(nproblem)
                        CALL init_lmbmpar(1)                    
                        CALL cpu_time(LMBMstart)   ! Start CPU timining
                        CALL lmbm(mc,f_subp,iout(1),iout(2),iout(3),iout(4),LMBMstart)
                        CALL copy_x_var(helppoint)
                        CALL deallocate_x_var()
                        CALL cpu_time(LMBMend)
                        !PRINT*, 'TIME LMBM', LMBMend-LMBMstart
                        
                        ! Update the counts of function and subgradient evaluations
                        DO i=kl,max_cl
                            nf(i) = nf(i)+iout(2)
                            nsub(i) = nsub(i)+iout(3)
                        END DO
                        
                        ! Save the old cluster center
                        DO i=1,features
                            xpoint(i)=y(i+features*(largest_ind-1))     ! Tallennetaan vanha centroid
                        END DO
                        
                        ! Copy the result of starting point auxiliary problem
                        DO i=features+1,2*features
                            xpoint(i)=helppoint(i-features)
                        END DO
                        
                        ! Check if the result point gave smaller function value than previous points
                        IF (f_subp < f_smallest_new) THEN
                            f_smallest_new = f_subp
                            DO i=1,2*features
                                bestpoint(i)=xpoint(i) ! Save the best starting point
                            END DO
                        END IF
                        
                        !----------------------------------------
                        ! Starting point auxiliary problem ends
                        !----------------------------------------
                        
                    ELSE ! This loop when we do not want to optimize starting point
                    
                        ! Starting points: cluster center + random data point
                        IF (nrandom <= noopt1) THEN ! VAIHTOEHTO 1
                            DO i=1,features
                                xpoint(i)=y(i+features*(largest_ind-1)) ! Save the first starting point
                            END DO
                            CALL random_number(beta)
                            DO WHILE (beta==0)
                                CALL random_number(beta)
                            END DO
                            beta2 = ceiling(nrecords_b*beta)
                            DO i=features+1,2*features
                                xpoint(i)=b(i-features,beta2) ! Save the second starting point
                            END DO
                            
                        ! Starting points: cluster center + center of random data points
                        ELSE IF (nrandom <= noopt1+noopt2) THEN
                            DO i=1,features
                                xpoint(i)=y(i+features*(largest_ind-1)) ! Save the first starting point
                            END DO
                            DO j=1,ncenter1
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_b*beta)
                                DO k=1,features
                                    center_table1(k,j)=b(k,beta2)
                                END DO
                            END DO
                            DO j=1,features
                                element_sum = 0.0_dp
                                DO k=1,ncenter1
                                    element_sum = element_sum + center_table1(j,k)
                                END DO
                                centroid(j) = element_sum / ncenter1
                            END DO
                            DO i=features+1,2*features
                                xpoint(i)=centroid(i-features) ! Save the second starting point
                            END DO
                        ! Starting points: cluster center + cluster center
                        ELSE IF (nrandom <= noopt1+noopt2+noopt3) THEN
                            DO i=1,features
                                xpoint(i)=y(i+features*(largest_ind-1)) ! Save the first starting point
                            END DO
                            DO i=features+1,2*features
                                xpoint(i)=y(i+features*(largest_ind-2)) ! Save the second starting point
                            END DO
                        
                        ! Starting points: random data point + random data point
                        ELSE IF (nrandom <= noopt1+noopt2+noopt3+noopt4) THEN
                            CALL random_number(beta)
                            DO WHILE (beta==0)
                                CALL random_number(beta)
                            END DO
                            beta2 = ceiling(nrecords_b*beta)
                            DO i=1,features
                                xpoint(i)=b(i,beta2) ! Save the first starting point
                            END DO
                            CALL random_number(beta)
                            DO WHILE (beta==0)
                                CALL random_number(beta)
                            END DO
                            beta2 = ceiling(nrecords_b*beta)
                            DO i=features+1,2*features
                                xpoint(i)=b(i-features,beta2) ! Save the second starting point
                            END DO
                        
                        ! Starting points: random data point + center of random data points
                        ELSE IF (nrandom <= noopt1+noopt2+noopt3+noopt4+noopt5) THEN
                            CALL random_number(beta)
                            DO WHILE (beta==0)
                                CALL random_number(beta)
                            END DO
                            beta2 = ceiling(nrecords_b*beta)
                            DO i=1,features
                                xpoint(i)=b(i,beta2) ! Save the first starting point
                            END DO
                            DO j=1,ncenter1
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_b*beta)
                                DO k=1,features
                                    center_table1(k,j)=b(k,beta2)
                                END DO
                            END DO
                            DO j=1,features
                                element_sum = 0.0_dp
                                DO k=1,ncenter1
                                    element_sum = element_sum + center_table1(j,k)
                                END DO
                                centroid(j) = element_sum / ncenter1
                            END DO
                            DO i=features+1,2*features
                                xpoint(i)=centroid(i-features) ! Save the second starting point
                            END DO
                        
                        ! Starting points: center of random data points + center of random data points
                        ELSE
                            DO j=1,ncenter1
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_b*beta)
                                DO k=1,features
                                    center_table1(k,j)=b(k,beta2)
                                END DO
                            END DO
                            DO j=1,features
                                element_sum = 0.0_dp
                                DO k=1,ncenter1
                                    element_sum = element_sum + center_table1(j,k)
                                END DO
                                centroid(j) = element_sum / ncenter1
                            END DO
                            DO i=1,features
                                xpoint(i)=centroid(i) ! Save the first starting point
                            END DO
                            DO j=1,ncenter1
                                CALL random_number(beta)
                                DO WHILE (beta==0)
                                    CALL random_number(beta)
                                END DO
                                beta2 = ceiling(nrecords_b*beta)
                                DO k=1,features
                                    center_table1(k,j)=b(k,beta2)
                                END DO
                            END DO
                            DO j=1,features
                                element_sum = 0.0_dp
                                DO k=1,ncenter1
                                    element_sum = element_sum + center_table1(j,k)
                                END DO
                                centroid(j) = element_sum / ncenter1
                            END DO
                            DO i=features+1,2*features
                                xpoint(i)=centroid(i-features) ! Save the second starting point
                            END DO
                        END IF
                    END IF

                    ! This loop when:
                        ! we want to solve 2-clustering auxiliary problem from all starting points OR
                        ! we want to solve 2-clustering auxiliary problem only from the best point and we already know that point
                    IF (all_start == 0 .OR. nrandom == n_randompoints) THEN 
                        
                        ! This loop when we want to solve 2-clustering auxiliary problem only from the best starting point
                        IF (all_start == 1) THEN
                            DO i=1,2*features
                                xpoint(i) = bestpoint(i) ! Set the best starting point into variable xpoint
                            END DO
                        END IF
                        
                        ! Initializations
                        neg_fvalue = 0
                        oldcoef = 1
                        
                        222 CONTINUE
                        
                        ! Check if we want to calculate better DBI indices (at the cost of function value)
                        IF (index1 == 1) THEN
                            CALL set_coefficient(bestpoint, neg_fvalue, oldcoef)
                        END IF
                        
                        ! Initializations
                        nproblem = 2
                        CALL set_nclust(2) ! Set the number of clusters into to 2 (next problem is 2-clustering auxiliary problem)
                        CALL init_problem(nproblem)
                        should_continue = .TRUE.
                        
                        DO WHILE (should_continue)
                        
                            !----------------------------------
                            ! 2-clustering auxiliary problem
                            !----------------------------------

                            ! Calls for LMBM
                            CALL allocate_x_var(2*features)
                            CALL init_x_var(xpoint)
                            CALL init_lmbmpar(kl)    
                            CALL cpu_time(LMBMstart)   ! Start CPU timining
                            CALL lmbm(mc,f_subp,iout(1),iout(2),iout(3),iout(4),LMBMstart)  
                            CALL copy_x_var(newpoint)
                            CALL deallocate_x_var()
                            CALL cpu_time(LMBMend)
                            !PRINT*, 'TIME LMBM', LMBMend-LMBMstart
                            
                            ! Update the counts of function and subgradient evaluations
                            DO i=kl,max_cl
                                nf(i) = nf(i)+iout(2)
                                nsub(i) = nsub(i)+iout(3)
                            END DO

                            neg_fvalue = 1
                            
                            IF (f_subp < 0.0_dp) GO TO 222  ! We do not want negative f_subp value
                            
                            ! Store the function value generated by LMBM in the table
                            f_table(nrandom) = f_subp
                            
                            ! Compute the lower and upper index limits in the cluster center vector for the cluster under consideration
                            l_limit = features*(largest_ind-1)+1
                            u_limit = features*largest_ind
                            
                            ! Save the first result point from 2-clustering auxiliary problem into the correct position in the cluster center table
                            DO i=l_limit,u_limit
                                x_table(i,nrandom)=newpoint(i-features*(largest_ind-1))
                            END DO
                            
                            ! Compute the lower and upper index limits in the cluster center vector for the new cluster center
                            l_limit = features*(kl-1)+1
                            u_limit = features*kl
                            
                            ! Save the second result point from 2-clustering auxiliary problem into the correct position in the cluster center table
                            DO i=l_limit,u_limit
                                x_table(i,nrandom)=newpoint(i-features*(kl-2))
                            END DO
                            
                            ! Copy result table into variable help_y
                            DO i=1,features*kl
                                help_y(i) = x_table(i,nrandom)
                            END DO
                            
                            ! Save smallest function value and smallest index
                            IF (f_subp < f_smallest) THEN
                                f_smallest = f_subp
                                smallest_index = nrandom
                            END IF
                            
                            ! Check if we want to delete outliers (if there are any)
                            IF (delete_outlier == 0) THEN
                                
                                ! Calculate cluster function values for every cluster
                                CALL func_a_sep(help_y, y_dim, kl, max_cl, f_a_sep, n_a_sep_help, indexes)
                                
                                ! Initializations
                                continuing = 0
                                n_outlier_2 = n_outlier
                                
                                ! Find the cluster with smallest number of outlier points
                                DO i=1,kl
                                    IF (n_a_sep_help(i) > 0) THEN
                                        IF (n_a_sep_help(i) < n_outlier_2) THEN
                                            continuing = 1 ! We found an outlier cluster
                                            n_outlier_2 = n_a_sep_help(i)
                                            change_ind = i ! The index of the outlier cluster
                                        END IF
                                    END IF
                                END DO
                                
                                IF (continuing == 0) THEN ! This loop when there are no outlier clusters
                                    should_continue = .FALSE.
                                ELSE ! Whis loop when there are outlier clusters
                                    loop_counter = loop_counter+1
                                    CALL modify_data(max_cl, change_ind, n_a_sep_help, indexes) ! Modify the data (delete outliers)
                                    CALL func_a_sep(y, y_dim, kl-1, max_cl, f_a_sep, n_a_sep, indexes) ! Compute new function values
                                    largestvalue = -3.40282347*10.**38
                                    ! Find the cluster with the highest cluster function value
                                    DO i = 1,kl-1
                                        IF (f_a_sep(i)>largestvalue) THEN
                                            largestvalue = f_a_sep(i)
                                            largest_ind = i
                                        END IF
                                    END DO
                                    CALL set_submatrix(largest_ind) ! Create a subdata set b containing the points of the cluster under consideration

                                END IF
                            
                            ELSE ! We do not want to delete outliers
                                should_continue = .FALSE.
                            END IF
                            
                        END DO  ! DO WHILE ends
                        
                    END IF
                    
                END DO
                
                ! Store the best result (cluster centers) into variable y
                DO i=1,features*kl
                    y(i)=x_table(i,smallest_index)
                END DO
                
                IF (usedmethod == 4) THEN ! We use 'method4' i.e. we calculate k-clustering problem
                
                    ! Deallocations
                    DEALLOCATE(xpoint)
                    DEALLOCATE(bestpoint)

                    nproblem = 1 ! The whole data a is in use
                    ALLOCATE(xpoint(kl*features)) ! Allocate xpoint
                    
                    ! Store the cluster centers into variable xpoint
                    DO i=1,kl*features
                        xpoint(i)=y(i)
                    END DO
                    
                    ! Initializations
                    CALL set_nclust(kl) ! The number of clusters is kl
                    neg_fvalue = 0
                    oldcoef = 1
                        
                    333 CONTINUE
                    
                    ! Check if we want to make DBI indices better (at the cost of function value)
                    IF (index2 == 1) THEN
                        CALL set_coefficient2(xpoint, neg_fvalue, oldcoef)
                    END IF
                    
                    !------------------------
                    ! k-clustering problem
                    !------------------------

                    ! Calls for LMBM
                    CALL allocate_x_var(kl*features)
                    CALL init_x_var(xpoint)
                    CALL init_problem(nproblem)  
                    CALL init_lmbmpar(kl)    
                    CALL cpu_time(LMBMstart)   ! Start CPU timining
                    CALL lmbm(mc,f_subp,iout(1),iout(2),iout(3),iout(4),LMBMstart)
                    CALL copy_x_var(xpoint)
                    CALL deallocate_x_var()
                    CALL cpu_time(LMBMend)
                    !PRINT*, 'TIME LMBM', LMBMend-LMBMstart

                    ! Update the counts of function and subgradient evaluations
                    DO i=kl,max_cl
                        nf(i) = nf(i)+iout(2)
                        nsub(i) = nsub(i)+iout(3)
                    END DO
                    
                    neg_fvalue = 1
                            
                    IF (f_subp < 0.0_dp) GO TO 333  ! We do not want negative f_subp value
                    
                    ! Save the result cluster centers into variable y
                    DO i=1,features*kl
                        y(i)=xpoint(i)
                    END DO

                END IF ! Extra step of 'method4' ends
                
                ! Store the cluster centers into table x
                DO i=1,features*kl
                    x(i,kl)=y(i)
                END DO 
                
                ! Calculate and store cluster function values
                CALL func_a_original(y, y_dim, f_a) ! Function value for original data a
                CALL func_a(y, y_dim, f_a_method4) ! Function value for data set a after removing outliers
                f(kl)=f_a
                f_method4(kl)=f_a_method4
                
                ! Calculate validity indices
                CALL check(kl,y,davbou,dunn)
                db(kl) = davbou ! Davies-Bouldin index (DBI)
                dn(kl) = dunn ! Dunn index (DI)
                
                ! Calculate predicted labels if kl=max_cl
                IF (kl == max_cl) THEN
                    ALLOCATE(predicted_labels(nrecords_a))
                    predicted_labels = 100 ! Initialize predicted labels
                    DO i=1,nrecords_a ! Iterate through the data points
                        smallest_value = 3.40282347*10.**38
                        DO j=1,kl ! Iterate through the clusters
                            value1 = 0.0_dp
                            DO k=1,nfeatures_a ! Iterate through the features
                                value1 = value1 + (y(k+nfeatures_a*(j-1))-a(k,i))**2
                            END DO
                            IF (value1 < smallest_value) THEN
                                smallest_value = value1
                                predicted_labels(i) = j-1 ! Determine which cluster the data point belongs to
                            END IF
                        END DO
                    END DO
                    ! Print predicted labels
                    !PRINT*, 'Predicted labels:'
                    !PRINT*, predicted_labels
                    DEALLOCATE(predicted_labels)
                END IF
                
                ! Deallocations
                DEALLOCATE(xpoint)
                DEALLOCATE(newpoint)
                DEALLOCATE(f_table)
                DEALLOCATE(x_table)
                
                ! Finish time of the iteration
                CALL cpu_time(f_time)
            
                !PRINT*, 'TIME', f_time-s_time
                !PRINT*, ''
                
                ! The CPU of the iteration              
                cpu(kl) = f_time-s_time

            END DO
            
            ! Final deallocations
            CALL deallocate_whole_data()

    END SUBROUTINE clustsplitter

END MODULE our_method

    
    
PROGRAM clust_splitter

    USE our_method
    
    IMPLICIT NONE

    CHARACTER(LEN=50) :: infile                             ! Infile data, for example 'iris.txt'    
    CHARACTER(LEN=50) :: outfile1                           ! Outfile data, for example 'clustsplitter_iris.txt'
    INTEGER :: features                                     ! Number of features in infile data                                 
    INTEGER :: records                                      ! Number of records in infile data
    INTEGER :: max_cl                                       ! Maximum number of clusters
    INTEGER :: used_method                                  ! Method to use
                                                                ! 4 = Clust-Splitter
                                                                ! 2 = Clust-Splitter without k-clustering problem  
    INTEGER :: k                                            ! Loop variable
    INTEGER :: task                                         ! Number of data sets in the loop
    INTEGER :: n_outlier                                    ! Limit for outlier points
                                                                ! if a cluster contains fewer than n_outlier points, we do not split the cluster
    INTEGER :: delete_outlier                               ! Do you want to delete outlier points from the data?
                                                                ! 0=yes, 1=no
    INTEGER :: opt_startingpoint                            ! Do you want to optimize starting points using the starting point auxiliary problem?
                                                                ! 1=yes, 0=no
    INTEGER :: ncenter1, ncenter2                           ! Number of random points used to calculate the average center in auxiliary problems
    INTEGER :: opt1, opt2, opt3                             ! Number of different starting points when opt_startingpoint=1
    INTEGER :: noopt1, noopt2, noopt3                       ! Number of different starting points when opt_startingpoint=0
    INTEGER :: noopt4, noopt5, noopt6                       ! Number of different starting points when opt_startingpoint=0
    INTEGER :: allstart                                     ! Calculate 2-clustering auxiliry problem from all starting points or just from the best one?
                                                                ! 1=best, 0=all
    INTEGER :: index1, index2                               ! Do you want to make Davies-Bouldin indices (DBI) better (at the cost of the function value)?
                                                                ! 1=yes 0=no
                                                                    ! index1=using 2-clustering auxiliary problem
                                                                    ! index2=using k-clustering problem
    REAL(KIND=dp) :: readingtime                            ! Time for reading the data
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: x         ! Cluster centers
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f           ! Function values for clusters of different sizes
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f_method4   ! Function values for clusters of different sizes for method 4 when we want to delete outliers
                                                                ! f_method4 differs from f only when delete_outlier=0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: nf          ! Number of function evaluations
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: nsub        ! Number of subgradient evaluations
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpu         ! CPU times (without data reading time)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: db          ! Davies-Bouldin indexes
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dn          ! Dunn indexes
    
    !----------------------------
    ! Initialize the parameters
    !----------------------------
    used_method = 4                 ! Method to use
                                        ! 4 = Clust-Splitter (default)
                                        ! 2 = Clust-Splitter without k-clustering problem
    max_cl = 25                     ! Number of clusters to calculate
                                        ! max_cl > 1, 25=default
    n_outlier = 6                   ! A cluster is an outlier-cluster if it contains fewer than n_outlier points
                                        ! n_outlier >= 0, 6=default
    delete_outlier = 1              ! Do you want to delete the data points that belong to clusters with fewer than n_outlier points
                                        ! 0=yes, 1=no (default)
    opt_startingpoint = 1           ! Do you want to optimize starting points using the starting point auxiliary problem
                                        ! 1=yes (default), 0=no
    allstart = 1                    ! Calculate 2-clustering auxiliry problem from all starting points or just from the best
                                        ! 1=best (default), 0=all
    ncenter1 = 10                   ! Number of random points used to calculate the average center in auxiliary problems
                                        ! This is applied when the distance of the average point from the cluster center does not matter
                                        ! ncenter1 > 0, 10=default
    ncenter2 = 7                    ! Number of random points used to calculate the average center in auxiliary problems
                                        ! This is applied when the distance of the average point from the cluster center matters
                                        ! ncenter2 > 0, 7=default
    index1 = 0                      ! Do you want to make DBI indices better (at the cost of the function value) using 2-clustering auxiliary problem
                                        ! 1=yes 0=no (default)
    index2 = 0                      ! Do you want to make DBI indices better (at the cost of the function value) using k-clustering problem
                                        ! 1=yes 0=no (default)
    
    !----------------------------------------------------
    ! Initialize the number of different starting points
    !----------------------------------------------------
    
    ! Starting points for starting point auxiliary problem when we optimize starting point (opt_startingpoint = 1)
    opt1 = 1    ! Average center of random data points (when the distance of the average point from the cluster center does not matter)
                    ! opt1 >= 0, 1=default
    opt2 = 1    ! Average center of random data points (when the distance of the average point from the cluster center matters)
                    ! opt2 >= 0, 1=default
    opt3 = 1    ! The cluster center
                    ! opt3=1, do not change
    
    ! Starting points for 2-clustering problem when we do not optimize starting point (opt_startingpoint = 0)
    noopt1 = 2  ! The cluster center + random data point
                    ! noopt1 >= 0, 2=default
    noopt2 = 2  ! The cluster center + the average of random data points
                    ! noopt2 >= 0, 2=default
    noopt3 = 1  ! The cluster center + The cluster center
                    ! noopt3 = 1, do not change
    noopt4 = 2  ! Random data point + random data point
                    ! noopt4 >= 0, 2=default
    noopt5 = 2  ! Random data point + the average of random data points
                    ! noopt5 >= 0, 2=default
    noopt6 = 2  ! The average of random data points + the average of random data points
                    ! noopt6 >= 0, 2=default
                    
                    
    !----------------------------------------------------
    ! Checking that parameter values are not forbidden
    !
    ! If they are, the default parameters are set
    !----------------------------------------------------
    
    IF (used_method /= 4 .AND. used_method /= 2) THEN
        PRINT*, 'The value of "used_method" was forbidden'
        used_method = 4   ! set the default value
    END IF
    IF (max_cl < 2) THEN
        PRINT*, 'The value of "max_cl" was forbidden'
        max_cl = 25   ! set the default value
    END IF
    IF (n_outlier < 0) THEN
        PRINT*, 'The value of "n_outlier" was forbidden'
        n_outlier = 6   ! set the default value
    END IF
    IF (delete_outlier /= 0 .AND. delete_outlier /= 1) THEN
        PRINT*, 'The value of "delete_outlier" was forbidden'
        delete_outlier = 1   ! set the default value
    END IF
    IF (opt_startingpoint /= 0 .AND. opt_startingpoint /= 1) THEN
        PRINT*, 'The value of "opt_startingpoint" was forbidden'
        opt_startingpoint = 1   ! set the default value
    END IF
    IF (allstart /= 0 .AND. allstart /= 1) THEN
        PRINT*, 'The value of "allstart" was forbidden'
        allstart = 1   ! set the default value
    END IF
    IF (ncenter1 < 1) THEN
        PRINT*, 'The value of "ncenter1" was forbidden'
        ncenter1 = 10   ! set the default value
    END IF
    IF (ncenter2 < 1) THEN
        PRINT*, 'The value of "ncenter2" was forbidden'
        ncenter2 = 7   ! set the default value
    END IF
    IF (index1 /= 0 .AND. index1 /= 1) THEN
        PRINT*, 'The value of "index1" was forbidden'
        index1 = 0   ! set the default value
    END IF
    IF (index2 /= 0 .AND. index2 /= 1) THEN
        PRINT*, 'The value of "index2" was forbidden'
        index2 = 0   ! set the default value
    END IF
    IF (opt1 < 0) THEN
        PRINT*, 'The value of "opt1" was forbidden'
        opt1 = 1   ! set the default value
    END IF
    IF (opt2 < 0) THEN
        PRINT*, 'The value of "opt2" was forbidden'
        opt2 = 1   ! set the default value
    END IF
    IF (opt3 /= 1) THEN
        PRINT*, 'The value of "opt3" was forbidden'
        opt3 = 1   ! set the default value
    END IF
    IF (noopt1 < 0) THEN
        PRINT*, 'The value of "noopt1" was forbidden'
        noopt1 = 2   ! set the default value
    END IF
    IF (noopt2 < 0) THEN
        PRINT*, 'The value of "noopt2" was forbidden'
        noopt2 = 2   ! set the default value
    END IF
    IF (noopt3 /= 1) THEN
        PRINT*, 'The value of "noopt3" was forbidden'
        noopt3 = 1   ! set the default value
    END IF
    IF (noopt4 < 0) THEN
        PRINT*, 'The value of "noopt4" was forbidden'
        noopt4 = 2   ! set the default value
    END IF
    IF (noopt5 < 0) THEN
        PRINT*, 'The value of "noopt5" was forbidden'
        noopt5 = 2   ! set the default value
    END IF
    IF (noopt6 < 0) THEN
        PRINT*, 'The value of "noopt6" was forbidden'
        noopt6 = 2   ! set the default value
    END IF
    
    
    DO task=1,1  ! Give here the number of data/datas you want to run
                    ! For example, task=2,2 if you only want to run data number 2
                    ! For example, task=1,3 if you want to run data numbers 1, 2, and 3
        
        !---------------------------------------------------------
        ! Give here the data, features, records, and the outfile
        !---------------------------------------------------------
        SELECT CASE(task)
            CASE(1)
                infile = 'd15112.txt'  ! Data set: D15112
                features = 2
                records = 15112
                outfile1 = 'clustsplitter_d15112.txt'  ! The name of the result file
            CASE(2)
                infile = 'pla85900.txt'  ! Data set: Pla85900
                features = 2
                records = 85900
                outfile1 = 'clustsplitter_pla85900.txt'  ! The name of the result file
            CASE(3)
                infile = 'Skin_segmentation.txt'  ! Data set: Skin Segmentation
                features = 3
                records = 245057
                outfile1 = 'clustsplitter_skin.txt'  ! The name of the result file
            CASE(4)
                infile = '3D_road.txt'  ! Data set: 3D Road Network
                features = 3
                records = 434874
                outfile1 = 'clustsplitter_3Droad.txt'  ! The name of the result file
            CASE(5)
                infile = 'shuttle.txt'  ! Data set: Shuttle Control
                features = 9
                records = 58000
                outfile1 = 'clustsplitter_shuttle.txt'  ! The name of the result file
            CASE(6)
                infile = 'Sensorless_drive_diagnosis.txt'  ! Data set: Sensorless Drive Diagnosis
                features = 49
                records = 58509
                outfile1 = 'clustsplitter_sensorless.txt'  ! The name of the result file
            CASE(7)
                infile = 'KEGGmetabolic.txt'  ! Data set: KEGG Metabolic
                features = 20
                records = 53413
                outfile1 = 'clustsplitter_kegg.txt'  ! The name of the result file
            CASE(8)
                infile = 'OnlineNews.txt'  ! Data set: Online News Popularity
                features = 58
                records = 39644
                outfile1 = 'clustsplitter_news.txt'  ! The name of the result file
            CASE(9)
                infile = 'gas.txt'  ! Data set: Gas Sensor Array Drift
                features = 128
                records = 13910
                outfile1 = 'clustsplitter_gas.txt'  ! The name of the result file
            CASE(10)
                infile = 'EEG_Eye_State.txt'  ! Data set: EEG Eye State
                features = 14
                records = 14980
                outfile1 = 'clustsplitter_eegeye.txt'  ! The name of the result file
            CASE(11)
                infile = 'isolet.txt'  ! Data set: ISOLET
                features = 616
                records = 7797
                outfile1 = 'clustsplitter_isolet.txt'  ! The name of the result file
            CASE(12)
                infile = 'gisette.txt'  ! Data set: Gisette
                features = 5000
                records = 13500
                outfile1 = 'clustsplitter_gisette.txt'  ! The name of the result file
            CASE(13)
                infile = 'MiniBooNE.txt'  ! Data set: MiniBooNE Particle Identification
                features = 50
                records = 130064
                outfile1 = 'clustsplitter_miniboone.txt'  ! The name of the result file
            CASE(14)
                infile = 'mfcc.txt'  ! Data set: MFCCs for Speech Emotion Recognition
                features = 58
                records = 85134
                outfile1 = 'clustsplitter_mfcc.txt'  ! The name of the result file
            CASE(15)
                infile = 'music.txt'  ! Data set: Music Analysis
                features = 518
                records = 106574
                outfile1 = 'clustsplitter_music.txt'  ! The name of the result file
            CASE(16)
                infile = 'prot_homo.txt'  ! Data set: Protein Homology
                features = 74
                records = 145751
                outfile1 = 'clustsplitter_protein.txt'  ! The name of the result file
            CASE(17)
                infile = 'rq3.txt'  ! Data set: Range Queries Aggregates
                features = 7
                records = 200000
                outfile1 = 'clustsplitter_rangeque.txt'  ! The name of the result file
            CASE(18)
                infile = 'Covertype.txt'  ! Data set: Covertype
                features = 10
                records = 581012
                outfile1 = 'clustsplitter_covertype.txt'  ! The name of the result file
            CASE(19)
                infile = 'MiniBooNE_uusi.txt'  ! Data set: MiniBooNE Particle Identification (values -999 have been deleted)
                features = 50
                records = 129596
                outfile1 = 'clustsplitter_miniboone_new.txt'  ! The name of the result file
            CASE(20)
                infile = 'generoitu_50pros_data10.txt'  ! Data set: Own data set
                features = 2
                records = 360
                outfile1 = 'clustsplitter_generoitu.txt'  ! The name of the result file
            CASE(21)
                infile = 'iris.txt' ! Data set: Iris
                features = 4
                records = 150
                outfile1 = 'clustsplitter_iris.txt'  ! The name of the result file
            CASE(22)
                infile = 'soybean.txt'  ! Data set: Soybean
                features = 35
                records = 47
                outfile1 = 'clustsplitter_soybean.txt'  ! The name of the result file
            CASE(23)
                infile = 'arcane_train.txt'  ! Data set: Arcane (training set only)
                features = 10000
                records = 100
                outfile1 = 'clustsplitter_arcane.txt'  ! The name of the result file
        END SELECT
        
        OPEN(46,file=outfile1) 
        WRITE(46,*)  ' k ', ' f(k) ', ' f_method4(k) ', ' Davies-Bouldin ', ' Dunn ', ' nf(k) ', ' nsub(k) ', ' cpu(k) '
                          
        ! Allocations
        ALLOCATE(x(features*max_cl,max_cl))  
        ALLOCATE(f(max_cl))
        ALLOCATE(f_method4(max_cl))
        ALLOCATE(nf(max_cl))
        ALLOCATE(nsub(max_cl))
        ALLOCATE(cpu(max_cl))
        ALLOCATE(db(max_cl))
        ALLOCATE(dn(max_cl))
        
        ! Call the Clust-Splitter method
        CALL clustsplitter(infile, used_method, max_cl, delete_outlier, n_outlier, opt1, opt2, opt3,&
            noopt1, noopt2, noopt3, noopt4, noopt5, noopt6, opt_startingpoint, ncenter1, ncenter2, allstart,&
            index1, index2, features, records, x, f, f_method4, db, dn, nf, nsub, cpu, readingtime)
        
        DO k=1,max_cl
            WRITE(46,*)  k, f(k), f_method4(k), db(k), dn(k), nf(k), nsub(k), cpu(k)
        END DO
        
        WRITE (46,*) '-----------------------------------------------------'
        WRITE(46,*) 'Time for reading the data ', readingtime
        
        CLOSE(46)
        
        ! Deallocations
        DEALLOCATE(x)
        DEALLOCATE(f)
        DEALLOCATE(f_method4)
        DEALLOCATE(nf)
        DEALLOCATE(nsub)
        DEALLOCATE(cpu)
        DEALLOCATE(db)
        DEALLOCATE(dn)
    
    END DO
    
END PROGRAM clust_splitter