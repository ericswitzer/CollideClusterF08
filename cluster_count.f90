MODULE cluster_count_module
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, error_unit, output_unit, iostat_end, iostat_eor
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: cluster_count
    
CONTAINS
    
    SUBROUTINE cluster_count(npart_total,r_cut,r,num_clust,npart_clust,clust_point,r_cut_sq_change)
	
	IMPLICIT NONE
	
	INTEGER*8,						        INTENT(in)	:: npart_total		! Number of particles
	INTEGER*8,								INTENT(out) :: num_clust		! Total number of clusters
    INTEGER*8								            :: num_clust_max	! Maximum possible total number of clusters
    INTEGER*8								            :: i		        ! Dummy index variable
    INTEGER*8								            :: j		        ! Dummy index variable
    INTEGER*8								            :: k,n		        ! Dummy index variable
    INTEGER*8								            :: icount	        ! Count variable to determine cluster end
	REAL*8,							        INTENT(in)	:: r_cut			! Cut off distance of potential
    REAL*8							                	:: r_cut_sq			! Cut off squared
    REAL*8,	   DIMENSION(3)			                	:: rjk 			    ! Distance between particle j and k
    REAL*8                  		                	:: rjk_sq           ! Squared distance between partickle j and k
	REAL*8,    DIMENSION(:,:), ALLOCATABLE, INTENT(in)	:: r				! Position array
	INTEGER*8, DIMENSION(:),   ALLOCATABLE, INTENT(out)	:: npart_clust		! Number of particles for cluster i
	INTEGER*8, DIMENSION(:),   ALLOCATABLE, INTENT(out) :: clust_point		! Array of which particles belong to which cluster
	INTEGER*8, DIMENSION(:),   ALLOCATABLE          	:: flag_index	 	! Counted flag for each particles
    INTEGER*8                                           :: cur_clust_point_index
    REAL*8, OPTIONAL                                 :: r_cut_sq_change  ! Override
    
    
	ALLOCATE(npart_clust(npart_total),clust_point(npart_total),flag_index(npart_total))

    num_clust_max = npart_total   
    r_cut_sq = r_cut**2.0d0
    IF(PRESENT(r_cut_sq_change)) THEN
        r_cut_sq = r_cut_sq_change
    END IF
	num_clust = 0       ! Equivalent to Dr. B's 'ncl' variable
	npart_clust = 0     ! Equivalent to Dr. B's 'nsc' variable
	clust_point = 0     ! Equivalent to Dr. B's 'nmc' variable    
	flag_index = 0      ! Equivalent to Dr. B's 'index' variable
    

	DO i=1,npart_total
        
        IF(flag_index(i) == 0) THEN
            num_clust = num_clust + 1
            
            IF(num_clust > num_clust_max) THEN
                WRITE(unit=error_unit,fmt='(a)') 'Error in count cluster routine: number of clusters &
                    & exceeds maximum possible value'
                STOP
            END IF
            
            npart_clust(num_clust) = npart_clust(num_clust) + 1
            
            ! START ERIC FIX TO DR. B CODE
            cur_clust_point_index = 0
            
            DO n=1,npart_total    
                cur_clust_point_index = cur_clust_point_index + npart_clust(n)
            END DO
            
            clust_point(cur_clust_point_index) = i
            
            ! clust_point(npart_clust(num_clust)) = i
            ! END ERIC FIX TO DR. B CODE
            
            flag_index(i) = 1
            j = i
            icount = 1
            
30          CONTINUE
            
            DO k=1,npart_total
            
                IF(flag_index(k) == 0) THEN
                    
                    rjk(:) = r(:,j) - r(:,k)
                    rjk_sq = SUM(rjk**2.0d0)
                    
                    IF (rjk_sq < r_cut_sq) THEN
                        flag_index(k) = 1
                        npart_clust(num_clust) = npart_clust(num_clust) + 1
                        ! START ERIC FIX TO DR. B CODE
                        cur_clust_point_index = cur_clust_point_index + 1
                        clust_point(cur_clust_point_index) = k
                        !clust_point(npart_clust(num_clust)) = k
                        ! END ERIC FIX TO DR. B CODE
                    END IF
                    
                END IF
            END DO
                            
            icount = icount + 1
                
            IF(icount > npart_clust(num_clust)) THEN
                WRITE(unit=output_unit,fmt='(a)') 'Cluster found!'
            ELSE
                j = clust_point(icount)
                GOTO 30
            END IF
                
        END IF
    END DO	
	
    
    DEALLOCATE(flag_index)
    
	RETURN
	
	END SUBROUTINE cluster_count
    
END MODULE cluster_count_module