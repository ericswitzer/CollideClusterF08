MODULE config_io_module
    USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, error_unit, output_unit, iostat_end, iostat_eor
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: read_config_atoms, snap_positions_vmd, time_stamp
    
CONTAINS
    
    SUBROUTINE read_config_atoms(filename,n,r,v)
        IMPLICIT NONE
            
        CHARACTER(len=*), INTENT(in) :: filename
        REAL*8, DIMENSION(:,:), OPTIONAL, INTENT(out) :: r, v            
        INTEGER*8, INTENT(inout) :: n
            
        INTEGER :: cnf_unit, ioerr, i
            
        ! Read r and v if they exist
        IF (PRESENT(r) .and. PRESENT (v)) THEN
            ! Open the position file
            OPEN(newunit=cnf_unit,file=filename//'_pos.dat',status='old',action='read',iostat=ioerr)
            IF (ioerr /= 0) THEN
                WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error opening ', filename//'_pos.dat', ioerr
                READ(unit=input_unit,fmt=*)
                STOP 'Error in read_config_atoms'
            END IF
            ! Read the position file
            READ(unit=cnf_unit,fmt=*,iostat=ioerr)
            DO i = 1, n
                READ(unit=cnf_unit,fmt=*,iostat=ioerr) r(:,i)
                IF (ioerr /= 0) THEN
                    WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error reading r from ', filename//'_pos.dat', ioerr
                    READ(unit=input_unit,fmt=*)
                    IF (ioerr == iostat_eor) THEN
                        WRITE(unit=error_unit, fmt='(a)') 'End of record'
                        READ(unit=input_unit,fmt=*)
                    END IF
                    IF (ioerr == iostat_end) THEN
                        WRITE(unit=error_unit, fmt='(a)') 'End of file'
                        READ(unit=input_unit,fmt=*)
                    END IF
                    STOP 'Error in read_config_atoms'
                END IF
            END DO
            ! Close the position file
            CLOSE (unit=cnf_unit)
            ! Open the velocity file
            OPEN(newunit=cnf_unit,file=filename//'_vel.dat',status='old',action='read',iostat=ioerr)
            IF (ioerr /= 0) THEN
                WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error opening ', filename//'_vel.dat', ioerr
                READ(unit=input_unit,fmt=*)
                STOP 'Error in read_config_atoms'
            END IF
            ! Read the velocity file
            READ(unit=cnf_unit,fmt=*,iostat=ioerr)
            DO i = 1, n
                READ (unit=cnf_unit,fmt=*,iostat=ioerr) v(:,i)
                IF (ioerr /= 0) THEN
                    WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error reading v from ', filename//'_vel.dat', ioerr
                    IF (ioerr == iostat_eor) THEN
                        WRITE(unit=error_unit,fmt='(a)') 'End of record'
                        READ(unit=input_unit,fmt=*)
                    END IF
                    IF (ioerr == iostat_end) THEN
                        WRITE(unit=error_unit,fmt='(a)') 'End of file'
                        READ(unit=input_unit,fmt=*)
                    END IF
                    STOP 'Error in read_config_atoms'
                END IF
            END DO
            ! Close the velocity file
            CLOSE (unit=cnf_unit)
        ELSE
            ! Read n
            OPEN(newunit=cnf_unit,file=filename//'_pos.dat',status='old',action='read',iostat=ioerr)
            READ(unit=cnf_unit,fmt=*,iostat=ioerr) n
            IF (ioerr /= 0) THEN
                WRITE(unit=error_unit,fmt='(a,a,i15)') 'Error reading n from ', filename//'_pos.dat', ioerr
                READ(unit=input_unit,fmt=*)
                IF (ioerr == iostat_eor) THEN
                    WRITE(unit=error_unit,fmt='(a)') 'End of record'
                    READ(unit=input_unit,fmt=*)
                END IF
                IF (ioerr == iostat_end) THEN
                    WRITE(unit=error_unit,fmt='(a)') 'End of file'
                    READ(unit=input_unit,fmt=*)
                END IF
                STOP 'Error in read_config_atoms'
            END IF
        END IF
    END SUBROUTINE read_config_atoms    

    SUBROUTINE snap_positions_vmd(nc1,nc2,n,r,filename,flag)
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(in) :: filename
        INTEGER*8, INTENT(in) :: nc1, nc2, n
        REAL*8, DIMENSION(:,:), INTENT(in) :: r
        CHARACTER(len=*), OPTIONAL, INTENT(in) :: flag
        INTEGER :: i,vmd_unit
        vmd_unit=45
        
        OPEN(unit=vmd_unit,file='.\Output\collision_'//filename//'.xyz',status='old',position='append')
        WRITE(unit=vmd_unit,fmt='(i15)') nc1
        WRITE(unit=vmd_unit,fmt=*)
        DO i=1,nc1
            WRITE(unit=vmd_unit,fmt='(a3,3f15.8)')'Ar ',r(:,i)
        END DO
        CLOSE(unit=vmd_unit)
        
        IF(PRESENT(flag)) THEN
            OPEN(unit=vmd_unit,file='.\Output\collision_target_'//flag//'.xyz',status='old',position='append')
            WRITE(unit=vmd_unit,fmt='(i15)') nc2
            WRITE(unit=vmd_unit,fmt=*)
            DO i=nc1+1,n
                WRITE(unit=vmd_unit,fmt='(a3,3f15.8)')'Ar ',r(:,i)
            END DO
            CLOSE(unit=vmd_unit)
        ELSE
            OPEN(unit=vmd_unit,file='.\Output\collision_target_'//flag//'.xyz',status='old',position='append')
            WRITE(unit=vmd_unit,fmt='(i15)') nc2
            WRITE(unit=vmd_unit,fmt=*)
            DO i=nc1+1,n
                WRITE(unit=vmd_unit,fmt='(a3,3f15.8)')'Ar ',r(:,i)
            END DO
            CLOSE(unit=vmd_unit)
        END IF

    END SUBROUTINE snap_positions_vmd
    
    SUBROUTINE time_stamp()
        IMPLICIT NONE

        CHARACTER(len=8)  :: date
        CHARACTER(len=10) :: time
        REAL              :: cpu

        CALL DATE_AND_TIME ( date, time )
        CALL CPU_TIME ( cpu )
        WRITE(unit=output_unit,fmt='(a,t45,a4,a1,a2,a1,a2)') 'Date: ', date(1:4), '/', date(5:6), '/', date(7:8)
        WRITE(unit=output_unit,fmt='(a,t47,a2,a1,a2,a1,a2)') 'Time: ', time(1:2), ':', time(3:4), ':', time(5:6)
        WRITE(unit=output_unit,fmt='(a,t40,f15.6)') 'CPU time: ', cpu
    END SUBROUTINE time_stamp
    
END MODULE config_io_module