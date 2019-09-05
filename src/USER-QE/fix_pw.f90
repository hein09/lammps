! this file reimplements the functionality of run_pwscf.f90

! we keep PWScf running as long as the fix is alive,
! and perform manual stepping tied to lammps
! This allows us to control pw while keeping extrapolation

MODULE fix_pw

  USE kinds,         ONLY: DP
  USE io_global,     ONLY: stdout, ionode, ionode_id
  USE util_param,    ONLY: utilout => stdout
  USE ions_base,     ONLY: nat

  IMPLICIT NONE

  CONTAINS

    PURE FUNCTION energy_pw() BIND(C)
        USE ISO_C_BINDING
        USE ener, ONLY: etot
        IMPLICIT NONE
        REAL(kind=C_DOUBLE) :: energy_pw
        energy_pw = etot
        RETURN
    END FUNCTION

    SUBROUTINE start_pw(comm, nimage, npool, ntask, nband, ndiag, infile, outfile, c_nat) BIND(C)
      USE ISO_C_BINDING
      USE command_line_options, ONLY: set_command_line
      USE mp_global,            ONLY: mp_startup
      USE environment,          ONLY: environment_start
      USE read_input,           ONLY: read_input_file
      USE check_stop,           ONLY: check_stop_init
      USE input_parameters,     ONLY: outdir, calculation, ion_dynamics
      USE mp_pools,          ONLY : intra_pool_comm
      USE mp_bands,          ONLY : intra_bgrp_comm, inter_bgrp_comm
      USE parallel_include
      IMPLICIT NONE
      !
      INTEGER, VALUE :: comm
      INTEGER (kind=C_INT), VALUE :: nimage, npool, ntask, nband, ndiag
      INTEGER :: comm_, nimage_, npool_, ntask_, nband_, ndiag_
      INTEGER (kind=C_INT), INTENT(OUT) :: c_nat
      CHARACTER (kind=C_CHAR), INTENT(IN) :: infile(*), outfile(*)
      CHARACTER(LEN=80) :: infile_, outfile_
      INTEGER :: i
      LOGICAL, SAVE :: first_time=.true.
      !
      ! ... Copying a string from C to Fortran is a bit ugly.
      infile_ = ' '
      DO i=1,80
          IF (infile(i) == C_NULL_CHAR) EXIT
          infile_ = TRIM(infile_) // infile(i)
      END DO
      outfile_ = ' '
      DO i=1,80
          IF (outfile(i) == C_NULL_CHAR) EXIT
          outfile_ = TRIM(outfile_) // outfile(i)
      END DO
      IF (ionode) THEN
        IF ( first_time ) THEN
          OPEN(UNIT=55, FILE=outfile_, STATUS='REPLACE', ACTION='WRITE', IOSTAT=i)
        ELSE
          OPEN(UNIT=55, FILE=outfile_, STATUS='OLD', ACTION='WRITE', IOSTAT=i)
        END IF
        IF ( i /= 0 ) THEN
          c_nat = -1
          RETURN
        END IF
        stdout = 55
        utilout = 55
      END IF
      !
      ! regular initialization
      nimage_ = nimage; npool_ = npool; ntask_ = ntask; nband_ = nband; ndiag_ = ndiag
      CALL set_command_line(nimage=nimage_, npool=npool_, ntg=ntask_, nband=nband_, ndiag=ndiag_)
      comm_ = comm
      CALL mp_startup(my_world_comm=comm_)
      CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm )
      CALL environment_start('PWSCF')
      CALL read_input_file('PW', infile_)
      ! overwrite some settings
      ! if we use 'md', we should benefit most easily from extrapolation
      calculation = 'md'
      ! dynamics not used as positions are integrated by lammps
      ion_dynamics = 'verlet'
      ! reroute secondary output
      outdir = outfile_(5:)
      ! continue initialization
      CALL iosys()
      CALL plugin_initialization()
    ! TODO: can we use the stop-mechanism somehow?
      CALL check_stop_init()
      CALL setup()
      CALL init_run()
    ! TODO: send positions a first time?
    !
      ! mark as executed
      first_time = .false.
      ! send back nat to calling code
      c_nat = nat
    END SUBROUTINE start_pw

    SUBROUTINE calc_pw(f) BIND(C)
      USE ISO_C_BINDING
      USE extrapolation, ONLY: update_file
      USE force_mod,     ONLY: force
      IMPLICIT NONE
      REAL(kind=c_double), INTENT(OUT) :: f(3, nat)
      ! scf calculation
      CALL electrons()
      ! force calculation
      CALL forces()
      ! print coordinates in expected place
      CALL output_tau(.false., .false.)
      ! exchange forces
      f = force
      ! make sure extrapolation works
      CALL update_file()
    END SUBROUTINE calc_pw

    SUBROUTINE update_pw(pos) BIND(c)
      USE ISO_C_BINDING
      USE extrapolation, ONLY: update_pot
      USE constants,     ONLY: bohr_radius_angs
      USE cell_base,     ONLY: alat, at
      USE ions_base,     ONLY: tau
      USE mp_world,      ONLY: world_comm
      USE mp,            ONLY: mp_bcast
      IMPLICIT NONE
      REAL(kind=c_double), INTENT(IN) :: pos(3, nat)
      REAL(DP), DIMENSION(3) :: com
      REAL(DP), DIMENSION(3) :: coc
      INTEGER :: i
      ! exchange and re-center positions
      IF (ionode) THEN
        ! center of cell
        coc = MATMUL(at, (/0.5d0, 0.5d0, 0.5d0/))
        ! center of molecule
        com = SUM(pos, dim=2) / nat
        DO i = 1, nat
            tau(:, i) = pos(:, i) / (alat * bohr_radius_angs) + coc - com
        END DO
      END IF
      CALL mp_bcast(tau, ionode_id, world_comm)
      ! update wavefunctions
      CALL update_pot()
      ! re-initialize atomic position-dependent quantities
      CALL hinit1()
    END SUBROUTINE update_pw

    SUBROUTINE end_pw(retval) BIND(C)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER (kind=C_INT), INTENT(OUT) :: retval
      CALL stop_run(retval)
      IF (ionode) THEN
        CLOSE(UNIT=55)
      END IF
      stdout = 6
      utilout = 6
    END SUBROUTINE end_pw

END MODULE fix_pw
