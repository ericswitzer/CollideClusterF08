PROGRAM main
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor
    USE collide_cluster_module, ONLY : collide_cluster
    USE fix_files_module, ONLY : fix_files

    IMPLICIT NONE
  
    INTEGER*8                              :: nsteps                ! Total number of steps
    INTEGER*8                              :: energy_snapshot       ! Snapshot size for energy measurement
    INTEGER*8                              :: position_snapshot     ! Snapshot size for energy measurement
    REAL*8                                 :: i                     ! Dummy real index
    REAL*8                                 :: j                     ! Dummy real index

    REAL*8                                 :: r_cut                 ! Potential cut off distance
    REAL*8                                 :: dt                    ! Time step size
    
    CHARACTER*15                           :: integration_type      ! Type of integration performed
    CHARACTER*9                            :: simulation_type       ! Simulation type: cluster to wall or cluster to cluster

    REAL*8,    DIMENSION(3)                :: r_sep                 ! Start seperation between c.o.m. of clusters
    REAL*8,    DIMENSION(3)                :: v_com                 ! Start velocities between c.o.m. of clusters
    
    INTEGER*8                              :: cohesive_flag_num     ! Flag for cohesive in 1 or 0
    REAL*8                                 :: cohesive              ! Cohesive LJ modification factor
    
    REAL*8                                 :: k_spring              ! Spring constant for wall
    REAL*8                                 :: lambda                ! Viscocity coefficient for wall
    
    CHARACTER*33                           :: clst_str              ! Current cluster name
        
    INTEGER*8                              :: multi_init_clst       ! Multi-initial cluster configuration flag
    INTEGER*8                              :: multi_init_clst_num   ! Number of initial cluster configurations to test over
    INTEGER*8                              :: clst_index = 1        ! Current cluster number
    INTEGER*8                              :: pos_index = 1         ! Current position number
    INTEGER*8                              :: vel_index = 1         ! Current position number
    
    INTEGER*8                              :: multi_pos_flag        ! Multi-position flag
    REAL*8                                 :: multi_pos_start       ! Multi-position start value
    REAL*8                                 :: multi_pos_stop        ! Multi-position stop value
    REAL*8                                 :: multi_pos_step        ! Multi-position step size
    INTEGER*8                              :: multi_vel_flag        ! Multi-velicity flag
    REAL*8                                 :: multi_vel_start       ! Multi-velocity start value
    REAL*8                                 :: multi_vel_stop        ! Multi-velocity stop value
    REAL*8                                 :: multi_vel_step        ! Multi-velocity step size
    
    INTEGER*8                              :: movie_flag            ! Movie flag 1 or 0
    
    INTEGER*8                              :: epsilon_modifier_flag      ! Epsilon modifier for mantle flag
    REAL*8                                 :: epsilon_modifier_strength  ! Epsilon magnitude for the mantle
    
    
    ! Fix to add routines after initial simulations were finished.
    ! CALL fix_files 
    ! Saved almost 34 hours of work by doing this!
    
    OPEN(unit=10,file='.\Input\config.in',status='old',action='read')
    READ(unit=10,fmt=*)
    READ(unit=10,fmt=*)
    READ(unit=10,fmt=*) r_cut
    READ(unit=10,fmt=*) r_sep(1)
    READ(unit=10,fmt=*) r_sep(2)
    READ(unit=10,fmt=*) r_sep(3)
    READ(unit=10,fmt=*) v_com(1)
    READ(unit=10,fmt=*) v_com(2)
    READ(unit=10,fmt=*) v_com(3)
    READ(unit=10,fmt=*) dt
    READ(unit=10,fmt=*) nsteps
    READ(unit=10,fmt=*) integration_type
    READ(unit=10,fmt=*) energy_snapshot
    READ(unit=10,fmt=*) position_snapshot
    READ(unit=10,fmt=*) cohesive_flag_num
    READ(unit=10,fmt=*) cohesive
    READ(unit=10,fmt=*) simulation_type
    READ(unit=10,fmt=*) k_spring
    READ(unit=10,fmt=*) lambda
    READ(unit=10,fmt=*) multi_init_clst
    READ(unit=10,fmt=*) multi_init_clst_num
    READ(unit=10,fmt=*) multi_pos_flag
    READ(unit=10,fmt=*) multi_pos_start
    READ(unit=10,fmt=*) multi_pos_stop
    READ(unit=10,fmt=*) multi_pos_step
    READ(unit=10,fmt=*) multi_vel_flag
    READ(unit=10,fmt=*) multi_vel_start
    READ(unit=10,fmt=*) multi_vel_stop
    READ(unit=10,fmt=*) multi_vel_step
    READ(unit=10,fmt=*) movie_flag
    READ(unit=10,fmt=*) epsilon_modifier_flag
    CLOSE(unit=10)
    
    ! Create a new stats file which will be used for analysis
    OPEN(unit=10,file='.\Output\cluster_stats.dat',status='replace')
    WRITE(unit=10,fmt='(a,f15.8)')'Potential Cut-off Value: ',r_cut
    WRITE(unit=10,fmt='(a,i8)')'Total number of steps: ',nsteps
    WRITE(unit=10,fmt='(a,f15.8)')'Cohesive factor: ',cohesive
    WRITE(unit=10,fmt='(a,f15.8)')'Differential time step: ',dt
    WRITE(unit=10,fmt='(a)')
    CLOSE(unit=10)
    
    ! Create a new check file for debugging purposes
    OPEN(unit=10,file='.\Output\check.dat',status='replace')
    CLOSE(10)
    
    ! If turning off movie, only snapshot the final state
    IF(movie_flag == 0) THEN
        position_snapshot = nsteps
    END IF
    
    ! If a multi-cluster flag is used, measure stats used in Saitoh-Hayakaway paper
    IF(multi_init_clst /= 0 .or. multi_pos_flag /= 0 .or. multi_vel_flag /= 0) THEN
    OPEN(unit=10,file='.\Output\cluster_data_saitoh_hayakawa.dat',status='replace')
    WRITE(unit=10,fmt='(a,a,a,a)')'Absolute Velocity    ','N_Total  ','N_cls    ','N_adh    '
    CLOSE(unit=10)
    
    ! Create a cluster adhere check file
    OPEN(unit=10,file='.\Output\cluster_data_adhere.dat',status='replace')
    WRITE(unit=10,fmt='(a,a,a,a,a,a)')'          v_imp','          E_int','        E_k_cls','        E_k_trg','          E_kin','     Adherence?'
    CLOSE(unit=10)
    END IF
 
    
    !!! BEGIN MULTIPLE COLLIDE ROUTINE !!!
    
    
    
    !!! Multiple clusters, positions, and velocities
    IF(multi_init_clst == 1 .and. multi_pos_flag /= 0 .and. multi_vel_flag /= 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple cluster, position, and velocity simulation'
        WRITE(unit=output_unit,fmt='(a)')
        IF (multi_pos_flag == 1 .or. multi_pos_flag == 2 .or. multi_pos_flag == 3) THEN
            IF (multi_vel_flag == 1 .or. multi_vel_flag == 2 .or. multi_vel_flag == 3) THEN
                DO clst_index=1,multi_init_clst_num,1   ! Loop over cluster numbers
                    pos_index = 1
                    DO i=multi_pos_start, multi_pos_stop, multi_pos_step ! Loop over positions
                        r_sep(multi_pos_flag) = i
                        vel_index = 1
                        DO j=multi_vel_start, multi_vel_stop, multi_vel_step ! Loop over velocities
                            v_com(multi_vel_flag) = j
                            ! Setup proper nomenclature for clst_str
                            
                            ! cluster is [0,10)         pos is [0,10)          vel is [0,10)
                            IF(clst_index < 10 .and. pos_index < 10 .and. vel_index < 10) THEN
                                WRITE(clst_str,'(a8,i1,a11,i1,a11,i1)') 'cluster0', clst_index, '-position00', pos_index, '-velocity00', vel_index
                            ! cluster is [10,100)       pos is [0,10)          vel is [0,10)
                            ELSEIF(clst_index >= 10 .and. pos_index < 10 .and. vel_index < 10) THEN
                                WRITE(clst_str,'(a7,i2,a11,i1,a11,i1)') 'cluster', clst_index, '-position00', pos_index, '-velocity00', vel_index
                            
                            ! cluster is [0,10)         pos is [10,100)        vel is [0,10)
                            ELSEIF(clst_index < 10 .and. pos_index >= 10 .and. pos_index < 100 .and. vel_index < 10) THEN
                                WRITE(clst_str,'(a8,i1,a10,i2,a11,i1)') 'cluster0', clst_index, '-position0', pos_index, '-velocity00', vel_index
                            ! cluster is [10,100)       pos is [10,100)        vel is [0,10)
                            ELSEIF(clst_index >= 10 .and. pos_index >= 10 .and. pos_index < 100 .and. vel_index < 10) THEN
                                WRITE(clst_str,'(a7,i2,a10,i2,a11,i1)') 'cluster', clst_index, '-position0', pos_index, '-velocity00', vel_index
                                
                            ! cluster is [0,10)         pos is [100,1000)      vel is [0,10)
                            ELSEIF(clst_index < 10 .and. pos_index >= 100 .and. vel_index < 10) THEN
                                WRITE(clst_str,'(a8,i1,a9,i3,a11,i1)') 'cluster0', clst_index, '-position', pos_index, '-velocity00', vel_index
                            ! cluster is [10,100)       pos is [100,1000)      vel is [0,10)
                            ELSEIF(clst_index >= 10 .and. pos_index >= 100 .and. vel_index < 10) THEN
                                WRITE(clst_str,'(a7,i2,a9,i3,a11,i1)') 'cluster', clst_index, '-position', pos_index, '-velocity00', vel_index    

                                
                                
                            ! cluster is [0,10)         pos is [0,10)          vel is [10,100)
                            ELSEIF(clst_index < 10 .and. pos_index < 10 .and. vel_index >= 10 .and. vel_index < 100) THEN
                                WRITE(clst_str,'(a8,i1,a11,i1,a10,i2)') 'cluster0', clst_index, '-position00', pos_index, '-velocity0', vel_index
                            ! cluster is [10,100)       pos is [0,10)          vel is [10,100)
                            ELSEIF(clst_index >= 10 .and. pos_index < 10 .and. vel_index >= 10 .and. vel_index < 100) THEN
                                WRITE(clst_str,'(a7,i2,a11,i1,a10,i2)') 'cluster', clst_index, '-position00', pos_index, '-velocity0', vel_index
                            
                            ! cluster is [0,10)         pos is [10,100)        vel is [10,100)
                            ELSEIF(clst_index < 10 .and. pos_index >= 10 .and. pos_index < 100 .and. vel_index >= 10 .and. vel_index < 100) THEN
                                WRITE(clst_str,'(a8,i1,a10,i2,a10,i2)') 'cluster0', clst_index, '-position0', pos_index, '-velocity0', vel_index
                            ! cluster is [10,100)       pos is [10,100)        vel is [10,100)
                            ELSEIF(clst_index >= 10 .and. pos_index >= 10 .and. pos_index < 100 .and. vel_index >= 10 .and. vel_index < 100) THEN
                                WRITE(clst_str,'(a7,i2,a10,i2,a10,i2)') 'cluster', clst_index, '-position0', pos_index, '-velocity0', vel_index
                                
                            ! cluster is [0,10)         pos is [100,1000)      vel is [10,100)
                            ELSEIF(clst_index < 10 .and. pos_index >= 100 .and. vel_index >= 10 .and. vel_index < 100) THEN
                                WRITE(clst_str,'(a8,i1,a9,i3,a10,i2)') 'cluster0', clst_index, '-position', pos_index, '-velocity0', vel_index
                            ! cluster is [10,100)       pos is [100,1000)      vel is [10,100)
                            ELSEIF(clst_index >= 10 .and. pos_index >= 100 .and. vel_index >= 10 .and. vel_index < 100) THEN
                                WRITE(clst_str,'(a7,i2,a9,i3,a10,i2)') 'cluster', clst_index, '-position', pos_index, '-velocity0', vel_index                                  
                                
                                
                                
                            ! cluster is [0,10)         pos is [0,10)          vel is [100,1000)
                            ELSEIF(clst_index < 10 .and. pos_index < 10 .and. vel_index >= 100) THEN
                                WRITE(clst_str,'(a8,i1,a11,i1,a9,i3)') 'cluster0', clst_index, '-position00', pos_index, '-velocity', vel_index
                            ! cluster is [10,100)       pos is [0,10)          vel is [100,1000)
                            ELSEIF(clst_index >= 10 .and. pos_index < 10 .and. vel_index >= 100) THEN
                                WRITE(clst_str,'(a7,i2,a11,i1,a9,i3)') 'cluster', clst_index, '-position00', pos_index, '-velocity', vel_index
                            
                            ! cluster is [0,10)         pos is [10,100)        vel is [100,1000)
                            ELSEIF(clst_index < 10 .and. pos_index >= 10 .and. pos_index < 100 .and. vel_index >= 100) THEN
                                WRITE(clst_str,'(a8,i1,a10,i2,a9,i3)') 'cluster0', clst_index, '-position0', pos_index, '-velocity', vel_index
                            ! cluster is [10,100)       pos is [10,100)        vel is [100,1000)
                            ELSEIF(clst_index >= 10 .and. pos_index >= 10 .and. pos_index < 100 .and. vel_index >= 100) THEN
                                WRITE(clst_str,'(a7,i2,a10,i2,a9,i3)') 'cluster', clst_index, '-position0', pos_index, '-velocity', vel_index
                                
                            ! cluster is [0,10)         pos is [100,1000)      vel is [100,1000)
                            ELSEIF(clst_index < 10 .and. pos_index >= 100 .and. vel_index >= 100) THEN
                                WRITE(clst_str,'(a8,i1,a9,i3,a9,i3)') 'cluster0', clst_index, '-position', pos_index, '-velocity', vel_index
                            ! cluster is [10,100)       pos is [100,1000)      vel is [100,1000)
                            ELSEIF(clst_index >= 10 .and. pos_index >= 100 .and. vel_index >= 100) THEN
                                WRITE(clst_str,'(a7,i2,a9,i3,a9,i3)') 'cluster', clst_index, '-position', pos_index, '-velocity', vel_index  
                                
    
                                
                            ELSE
                                WRITE(unit=output_unit,fmt='(a)') 'Out of bounds - check input values'
                            END IF
                            ! Collide Clusters
                            
                            WRITE(unit=output_unit,fmt='(a,i,f15.8,f15.8)') 'On cluster, position, and velocity: ',clst_index,i,j
                            CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                            & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)                            
                            vel_index = vel_index + 1
                        END DO
                        pos_index = pos_index + 1
                    END DO
                END DO
            ELSE
                WRITE(unit=output_unit,fmt='(a,i)') 'Velocity flag has incorrect value: ',multi_vel_flag      
            END IF
        ELSE
            WRITE(unit=output_unit,fmt='(a,i)') 'Position flag has incorrect value: ',multi_pos_flag
        END IF
        
    !!! Multiple clusters and positions
    ELSEIF(multi_init_clst == 1 .and. multi_pos_flag /= 0 .and. multi_vel_flag == 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple cluster and position simulation'
        WRITE(unit=output_unit,fmt='(a)')        
        IF (multi_pos_flag == 1 .or. multi_pos_flag == 2 .or. multi_pos_flag == 3) THEN         
            DO clst_index=1,multi_init_clst_num,1   ! Loop over cluster numbers
                pos_index = 1
                DO j=multi_pos_start, multi_pos_stop, multi_pos_step ! Loop over positions
                    r_sep(multi_pos_flag) = j 
                    ! Setup proper nomenclature for clst_str
                    IF(clst_index < 10 .and. pos_index < 10) THEN
                        WRITE(clst_str,'(a8,i1,a11,i1,a12)') 'cluster0', clst_index, '-position00', pos_index, '-velocity001'
                    ELSEIF(clst_index >= 10 .and. pos_index < 10) THEN
                        WRITE(clst_str,'(a7,i2,a11,i1,a12)') 'cluster', clst_index, '-position00', pos_index, '-velocity001'
                    ELSEIF(clst_index < 10 .and. pos_index >= 10 .and. pos_index < 100) THEN
                        WRITE(clst_str,'(a8,i1,a10,i2,a12)') 'cluster0', clst_index, '-position0', pos_index, '-velocity001'
                    ELSEIF(clst_index >= 10 .and. pos_index >= 10 .and. pos_index < 100) THEN
                        WRITE(clst_str,'(a7,i2,a10,i2,a12)') 'cluster', clst_index, '-position0', pos_index, '-velocity001'
                    ELSEIF(clst_index < 10 .and. pos_index >= 100) THEN
                        WRITE(clst_str,'(a8,i1,a9,i3,a12)') 'cluster0', clst_index, '-position', pos_index, '-velocity001'
                    ELSEIF(clst_index >= 10 .and. pos_index >= 100) THEN
                        WRITE(clst_str,'(a7,i2,a9,i3,a12)') 'cluster', clst_index, '-position', pos_index, '-velocity001'
                    ELSE
                        WRITE(unit=output_unit,fmt='(a)') 'Out of bounds - check input values'
                    END IF
                    ! Collide Clusters
                    WRITE(unit=output_unit,fmt='(a,i,f15.8)') 'On cluster, position: ',clst_index,j
                    CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                    & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
                    ! Advance position index
                    pos_index = pos_index + 1
                END DO
            END DO
        ELSE
            WRITE(unit=output_unit,fmt='(a,i)') 'Position flag has incorrect value: ',multi_pos_flag
        END IF

    !!! Multiple clusters and velocities
    ELSEIF(multi_init_clst == 1 .and. multi_pos_flag == 0 .and. multi_vel_flag /= 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple cluster and velocity simulation'
        WRITE(unit=output_unit,fmt='(a)')        
        IF (multi_vel_flag == 1 .or. multi_vel_flag == 2 .or. multi_vel_flag == 3) THEN         
            DO clst_index=1,multi_init_clst_num,1   ! Loop over cluster numbers
                vel_index = 1
                DO j=multi_vel_start, multi_vel_stop, multi_vel_step ! Loop over velocities
                    v_com(multi_vel_flag) = j 
                    ! Setup proper nomenclature for clst_str
                    IF(clst_index < 10 .and. vel_index < 10) THEN
                        WRITE(clst_str,'(a8,i1,a23,i1)') 'cluster0', clst_index, '-position001-velocity00', vel_index
                    ELSEIF(clst_index >= 10 .and. vel_index < 10) THEN
                        WRITE(clst_str,'(a7,i2,a23,i1)') 'cluster', clst_index, '-position001-velocity00', vel_index
                    ELSEIF(clst_index < 10 .and. vel_index >= 10 .and. vel_index < 100) THEN
                        WRITE(clst_str,'(a8,i1,a22,i2)') 'cluster0', clst_index, '-position001-velocity0', vel_index
                    ELSEIF(clst_index >= 10 .and. vel_index >= 10 .and. vel_index < 100) THEN
                        WRITE(clst_str,'(a7,i2,a22,i2)') 'cluster', clst_index, '-position001-velocity0', vel_index
                    ELSEIF(clst_index < 10 .and. vel_index >= 100) THEN
                        WRITE(clst_str,'(a8,i1,a21,i3)') 'cluster0', clst_index, '-position001-velocity', vel_index
                    ELSEIF(clst_index >= 10 .and. vel_index >= 100) THEN
                        WRITE(clst_str,'(a7,i2,a21,i3)') 'cluster', clst_index, '-position001-velocity', vel_index
                    ELSE
                        WRITE(unit=output_unit,fmt='(a)') 'Out of bounds - check input values'
                    END IF
                    ! Collide Clusters
                    WRITE(unit=output_unit,fmt='(a,i,f15.8)') 'On cluster, velocities: ',clst_index,j
                    CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                    & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
                    ! Advance velocity index
                    vel_index = vel_index + 1
                END DO
            END DO
        ELSE
            WRITE(unit=output_unit,fmt='(a,i)') 'Velocity flag has incorrect value: ',multi_vel_flag
        END IF    
           
    !!! Multiple positions and velocities
    ELSEIF(multi_init_clst == 0 .and. multi_pos_flag /= 0 .and. multi_vel_flag /= 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple position and velocity simulation'
        WRITE(unit=output_unit,fmt='(a)')                
        IF (multi_pos_flag == 1 .or. multi_pos_flag == 2 .or. multi_pos_flag == 3) THEN
            IF (multi_vel_flag == 1 .or. multi_vel_flag == 2 .or. multi_vel_flag == 3) THEN
                pos_index = 1
                DO i=multi_pos_start, multi_pos_stop, multi_pos_step ! Loop over positions
                    r_sep(multi_pos_flag) = i
                    vel_index = 1
                    DO j=multi_vel_start, multi_vel_stop, multi_vel_step ! Loop over velocities 
                        v_com(multi_vel_flag) = j
                        ! Setup proper nomenclature for clst_str below
                        
                        ! r is [0,10)       v is [0,10)
                        IF(pos_index < 10 .and. vel_index < 10) THEN
                            WRITE(clst_str,'(a20,i1,a11,i1)') 'cluster01-position00', pos_index, '-velocity00', vel_index
                        ! r is [10,100)     v is [0,10)
                        ELSEIF(pos_index >=10 .and. pos_index < 100 .and. vel_index < 10) THEN
                            WRITE(clst_str,'(a19,i2,a11,i1)') 'cluster01-position0', pos_index, '-velocity00', vel_index
                        ! r is [100,1000)   v is [0,10)
                        ELSEIF(pos_index >= 100 .and. vel_index < 10) THEN
                            WRITE(clst_str,'(a18,i3,a11,i1)') 'cluster01-position', pos_index, '-velocity00', vel_index

                        ! r is [0,10)       v is [10,100)
                        ELSEIF(pos_index < 10 .and. vel_index >= 10 .and. vel_index < 100) THEN
                            WRITE(clst_str,'(a20,i1,a10,i2)') 'cluster01-position00', pos_index, '-velocity0', vel_index
                        ! r is [10,100)     v is [10,100)
                        ELSEIF(pos_index >=10 .and. pos_index < 100 .and. vel_index >= 10 .and. vel_index < 100) THEN
                            WRITE(clst_str,'(a19,i2,a10,i2)') 'cluster01-position0', pos_index, '-velocity0', vel_index
                        ! r is [100,1000)   v is [10,100)
                        ELSEIF(pos_index >= 100 .and. vel_index >= 10 .and. vel_index < 100) THEN
                            WRITE(clst_str,'(a18,i3,a10,i2)') 'cluster01-position', pos_index, '-velocity0', vel_index
                        
                        ! r is [0,10)       v is [100,1000)
                        ELSEIF(pos_index < 10 .and. vel_index >= 100) THEN
                            WRITE(clst_str,'(a20,i1,a9,i3)') 'cluster01-position00', pos_index, '-velocity', vel_index
                        ! r is [10,100)     v is [100,1000)
                        ELSEIF(pos_index >=10 .and. pos_index < 100 .and. vel_index >= 100) THEN
                            WRITE(clst_str,'(a19,i2,a9,i3)') 'cluster01-position0', pos_index, '-velocity', vel_index
                        ! r is [100,1000)   v is [100,1000)
                        ELSEIF(pos_index >= 100 .and. vel_index >= 100) THEN
                            WRITE(clst_str,'(a18,i3,a9,i3)') 'cluster01-position', pos_index, '-velocity', vel_index
                            
                        ELSE
                            WRITE(unit=output_unit,fmt='(a)') 'Out of bounds - check input values'
                        END IF
                        ! Collide Clusters
                        WRITE(unit=output_unit,fmt='(a,f15.8,f15.8)') 'On position, velocity: ',i,j
                        CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                        & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
                        ! Advance Velocity Index
                        vel_index = vel_index + 1
                    END DO
                    
                    ! Advance Position Index
                    pos_index = pos_index + 1
                END DO
            ELSE
                WRITE(unit=output_unit,fmt='(a,i)') 'Velocity flag has incorrect value: ',multi_vel_flag  
            END IF            
        ELSE
            WRITE(unit=output_unit,fmt='(a,i)') 'Position flag has incorrect value: ',multi_pos_flag
        END IF       
        
    !!! Multiple clusters  
    ELSEIF(multi_init_clst == 1 .and. multi_pos_flag == 0 .and. multi_vel_flag == 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple cluster simulation'
        WRITE(unit=output_unit,fmt='(a)')        
        DO clst_index = 1, multi_init_clst_num ! Loop over clusters
            IF(clst_index < 10 ) THEN
                WRITE(clst_str,'(a8,i1,a24)') 'cluster0', clst_index, '-position001-velocity001'
            ELSE
                WRITE(clst_str,'(a7,i2,a24)') 'cluster',clst_index, '-position001-velocity001'
            END IF
            WRITE(unit=output_unit,fmt='(a,i)') 'On cluster: ',clst_index
            CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                    & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
        END DO
        
    !!! Multiple positions  
    ELSEIF(multi_init_clst == 0 .and. multi_pos_flag /= 0 .and. multi_vel_flag == 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple position simulation'
        WRITE(unit=output_unit,fmt='(a)')        
        IF (multi_pos_flag == 1 .or. multi_pos_flag == 2 .or. multi_pos_flag == 3) THEN
            pos_index = 1
            DO i=multi_pos_start, multi_pos_stop, multi_pos_step ! Loop over positions
                r_sep(multi_pos_flag) = i
                IF(pos_index < 10 ) THEN
                    WRITE(clst_str,'(a20,i1,a12)') 'cluster01-position00', pos_index, '-velocity001'
                ELSEIF(pos_index >= 10 .and. pos_index < 100) THEN
                    WRITE(clst_str,'(a19,i2,a12)') 'cluster01-position0', pos_index, '-velocity001'
                ELSEIF(pos_index >= 100) THEN
                    WRITE(clst_str,'(a18,i3,a12)') 'cluster01-position', pos_index, '-velocity001'
                ELSE
                    WRITE(unit=output_unit,fmt='(a)') 'Out of bounds - check input values'
                END IF
                WRITE(unit=output_unit,fmt='(a,f15.8)') 'On position: ',i
                CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                        & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
                pos_index = pos_index + 1
            END DO
        ELSE
            WRITE(unit=output_unit,fmt='(a,i)') 'Position flag has incorrect value: ',multi_pos_flag
        END IF
        
    !!! Multiple velocities
    ELSEIF(multi_init_clst == 0 .and. multi_pos_flag == 0 .and. multi_vel_flag /= 0) THEN
        WRITE(unit=output_unit,fmt='(a)') 'Starting multiple velocity simulation'
        WRITE(unit=output_unit,fmt='(a)')        
        IF (multi_vel_flag == 1 .or. multi_vel_flag == 2 .or. multi_vel_flag == 3) THEN
            vel_index = 1
            DO i=multi_vel_start, multi_vel_stop, multi_vel_step ! Loop over positions
                v_com(multi_vel_flag) = i
                IF(vel_index < 10 ) THEN
                    WRITE(clst_str,'(a32,i1)') 'cluster01-position001-velocity00', vel_index
                ELSEIF(vel_index >= 10 .and. vel_index < 100) THEN
                    WRITE(clst_str,'(a31,i2)') 'cluster01-position001-velocity0', vel_index
                ELSEIF(vel_index >= 100) THEN
                    WRITE(clst_str,'(a30,i3)') 'cluster01-position001-velocity', vel_index
                ELSE
                    WRITE(unit=output_unit,fmt='(a)') 'Out of bounds - check input values'
                END IF
                WRITE(unit=output_unit,fmt='(a,f15.8)') 'On velocity: ',i
                CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                        & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
                vel_index = vel_index + 1
            END DO
        ELSE
            WRITE(unit=output_unit,fmt='(a,i)') 'Velocity flag has incorrect value: ',multi_vel_flag
        END IF
        
    !!! Nothing Multiple
    ELSE
        clst_str = 'cluster01-position001-velocity001'
        WRITE(unit=output_unit,fmt='(a)') 'Starting single-run simulation'
        CALL collide_cluster(r_cut,r_sep,v_com,dt,nsteps,integration_type,energy_snapshot,position_snapshot,&
                & cohesive_flag_num,cohesive,simulation_type,k_spring,lambda,clst_str,epsilon_modifier_flag)
    END IF
    
    
    
    
    
    !!! END MULTIPLE COLLIDE ROUTINE !!!
    
    
    
    
    
    
    WRITE(unit=output_unit,fmt='(a)') 'Simulation Complete!'
    WRITE(unit=output_unit,fmt='(a)') 'Press any key to exit'
    READ(unit=input_unit,fmt=*)
    STOP
    
END PROGRAM main