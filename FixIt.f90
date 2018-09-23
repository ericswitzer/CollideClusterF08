MODULE fix_files_module
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, error_unit, output_unit, iostat_end, iostat_eor
    USE            :: cluster_count_module, ONLY : cluster_count
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: fix_files
    
CONTAINS
    
    SUBROUTINE fix_files()
        
        IMPLICIT NONE
    
        INTEGER*8                              :: nc1
        INTEGER*8                              :: i,j,k                     ! Dummy index

        REAL*8,    DIMENSION(:,:), ALLOCATABLE :: r                   ! Postion, Velocity, Force, and Displacement arrays
            
        INTEGER*8                              :: num_clust             ! Number of clusters after collision
        INTEGER*8, DIMENSION(:), ALLOCATABLE   :: npart_clust           ! Total number of particles in each cluster
        INTEGER*8, DIMENSION(:), ALLOCATABLE   :: clust_point           ! Pointer of cluster(i)
        REAL*8,    DIMENSION(:,:), ALLOCATABLE :: r_1                   ! Subset of r - particles belonging to initial cluster
        INTEGER*8                              :: biggest_cluster_index ! Location of biggest cluster in pointer array
        INTEGER*8                              :: npart_adsorb       
        INTEGER*8                              :: biggest_cluster_npart      
        INTEGER*8                              :: biggest_cluster_start       
        INTEGER*8                              :: biggest_cluster_end
        REAL*8                                 :: r_adsorb,r_cut
        REAL*8                                 :: r_adsorb_sq
        INTEGER*8                              :: r_temp
        REAL*8, DIMENSION(3)                   :: r_adsorb_i
        REAL*8, DIMENSION(3)                   :: r_adsorb_ij
        REAL*8                                 :: r_adsorb_ij_sq
        INTEGER*8                              :: file_num = 1
        CHARACTER*2                            :: file_num_str
        CHARACTER*2                            :: throw
        
        OPEN(unit=10,file='.\FixMe\cluster_stats.dat',status='replace')
        WRITE(unit=10,fmt='(a,a,a,a)')'Run #    ','N    ','N_cls    ','N_adh    '
        CLOSE(10)
        
        DO k=1,50
            file_num = k
            
            IF(file_num < 10 ) THEN
                WRITE(file_num_str,'(a1,i1)')'0',file_num
            ELSE
                WRITE(file_num_str,'(i2)')file_num
            END IF
        
            ALLOCATE(r(3,14700))
        
            OPEN(unit=10,file='.\FixMe\collision_cluster01-position001-velocity0'//file_num_str//'.xyz',status='old',action='read')
            DO i=1,2720
                READ(unit=10,fmt=*)
            END DO
            DO i = 1, 300
                READ(unit=10,fmt=*)throw,r(:,i)
            END DO
            CLOSE(10)
        
            OPEN(unit=10,file='.\FixMe\collision_target_cluster01-position001-velocity0'//file_num_str//'.xyz',status='old',action='read')
            DO i=1,129620
                READ(unit=10,fmt=*)
            END DO
            DO i = 1, 14400
                READ(unit=10,fmt=*)throw,r(:,300+i)
            END DO
            CLOSE(10)
        
            nc1 = 300
            r_cut = 3.0d0

        
        
            ALLOCATE (npart_clust(nc1),clust_point(nc1))
            clust_point = 0
            ALLOCATE (r_1(3,nc1))
            r_1 = r(:,:nc1)
            npart_adsorb = 0
            r_adsorb = 1.6
            r_adsorb_sq = r_adsorb**2.0d0
            CALL cluster_count(nc1,r_cut,r_1,num_clust,npart_clust,clust_point,r_adsorb_sq)
            OPEN(unit=10,file='.\FixMe\cluster_stats.dat',status='old',position='append')
            biggest_cluster_index = MAXLOC(npart_clust,DIM=1)
            biggest_cluster_npart = npart_clust(biggest_cluster_index)
            biggest_cluster_end = 0
                
            DO i=1,biggest_cluster_index
                biggest_cluster_end = biggest_cluster_end + npart_clust(i)
            END DO
                
            biggest_cluster_start = biggest_cluster_end - biggest_cluster_npart + 1
            IF (biggest_cluster_start == 0) THEN
                biggest_cluster_start = 1
            END IF
            
            DO i=biggest_cluster_start,biggest_cluster_end,1
                r_temp = clust_point(i)
                r_adsorb_i(:) = r(:,r_temp)
                    
                DO j=1,14400
                    r_adsorb_ij(:) = r_adsorb_i(:) - r(:,nc1+j)
                    r_adsorb_ij_sq = SUM(r_adsorb_ij**2.0d0)
                
                    IF (r_adsorb_ij_sq < r_adsorb_sq) THEN
                        npart_adsorb = npart_adsorb + 1
                        EXIT
                    END IF
                        
                END DO
                        
            END DO
            
            !WRITE(unit=10,fmt='(a)')'!START Saitoh and Hayakawa Check!'
            !WRITE(unit=10,fmt='(a,a2)')'Run Number = ',file_num_str
            !WRITE(unit=10,fmt='(a,i8)')'Start number of particles in cluster = ',nc1
            !WRITE(unit=10,fmt='(a,i8)')'End number of particles in largest cluster = ',MAXVAL(npart_clust)
            !WRITE(unit=10,fmt='(a,i8)')'End number of particles in largest surviving cluster adsorbed = ',npart_adsorb
            !WRITE(unit=10,fmt='(a)')'!END Saitoh and Hayakawa Check!'
            !WRITE(unit=10,fmt=*)
            WRITE(unit=10,fmt='(i2,i8,i8,i8)')file_num,nc1,MAXVAL(npart_clust),npart_adsorb
            CLOSE(unit=10)
            DEALLOCATE (npart_clust,clust_point,r_1,r)
        
        END DO
    END SUBROUTINE fix_files
                
END MODULE fix_files_module