! this file reimplements the functionality of pwscf.f90 and run_pwscf.f90

! we keep PWScf running as long as the fix is alive,
! and perform manual stepping tied to LAMMPS
! This allows us to control pw while keeping extrapolation

MODULE fix_pw

  USE ISO_C_BINDING
  USE kinds,         ONLY: DP
  USE io_global,     ONLY: stdout, ionode, ionode_id
  USE util_param,    ONLY: utilout => stdout
  USE ions_base,     ONLY: nat

  IMPLICIT NONE

  INTEGER (kind=C_LONG) :: nec = 0

  CONTAINS

    PURE FUNCTION energy_pw() BIND(C)
      USE ener,      ONLY: etot
      USE constants, ONLY: rytoev
      IMPLICIT NONE
      REAL(kind=C_DOUBLE) :: energy_pw
      energy_pw = etot * rytoev
      RETURN
    END FUNCTION

    SUBROUTINE start_pw(comm, partitions, infile, outfile, c_nat, pos_qm, c_nec) BIND(C)
      USE command_line_options, ONLY: set_command_line
      USE cell_base,            ONLY: alat, at
      USE control_flags,        ONLY: gamma_only
      USE constants,            ONLY: bohr_radius_angs
      USE mp_global,            ONLY: mp_startup
      USE environment,          ONLY: environment_start
      USE read_input,           ONLY: read_input_file
      USE check_stop,           ONLY: check_stop_init
      USE input_parameters,     ONLY: outdir, calculation, &
                                      ion_dynamics, restart_mode, &
                                      rd_pos
      USE mp,                   ONLY: mp_bcast
      USE mp_pools,             ONLY: intra_pool_comm
      USE mp_bands,             ONLY: intra_bgrp_comm, inter_bgrp_comm
      USE mp_diag,              ONLY: mp_start_diag
      USE parallel_include
      IMPLICIT NONE
      !
      ! MPI arguments
      INTEGER, VALUE :: comm
      INTEGER (kind=C_INT), DIMENSION(4) :: partitions
      INTEGER :: comm_, npool_, ntask_, nband_, ndiag_
      ! Filenames
      CHARACTER (kind=C_CHAR), INTENT(IN) :: infile(*), outfile(*)
      ! QM atoms
      INTEGER (kind=C_INT), INTENT(INOUT) :: c_nat
      REAL(kind=c_double), INTENT(IN) :: pos_qm(3, c_nat)
      ! Number of EC atoms
      INTEGER (kind=C_LONG), INTENT(IN), VALUE :: c_nec
      ! local variables
      CHARACTER(LEN=80) :: infile_, outfile_
      INTEGER :: i
      LOGICAL, SAVE :: first_time=.true.
      REAL(DP), DIMENSION(3) :: com
      REAL(DP), DIMENSION(3) :: coc
      !
      ! handle filenames, redirect output to outfile
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
      ! mpi initialization
      npool_ = partitions(1); ntask_ = partitions(2)
      nband_ = partitions(3); ndiag_ = partitions(4)
      CALL set_command_line(npool=npool_, ntg=ntask_, nband=nband_, ndiag=ndiag_)
      comm_ = comm
      CALL mp_startup(my_world_comm=comm_, start_images=.true.)
      CALL mp_start_diag(ndiag_, comm_, intra_bgrp_comm, do_distr_diag_inside_bgrp_=.true.)
      CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm )
      ! read input file
      CALL environment_start('PWSCF')
      CALL read_input_file('PW', infile_)
      ! Check if input contains suitable number of atoms
      IF (c_nat /= nat) THEN
        ! send back mismatching nat to calling code
        c_nat = nat
        ! early return, not recoverable
        RETURN
      END IF
      ! initialize positions from lammps
      IF (ionode) THEN
        ! center of cell
        coc = MATMUL(at, (/0.5d0, 0.5d0, 0.5d0/))
        ! center of molecule
        com = SUM(pos_qm, dim=2) / (nat * alat * bohr_radius_angs)
        DO i = 1, nat
            rd_pos(:, i) = pos_qm(:, i) / (alat * bohr_radius_angs) + coc - com
        END DO
      END IF
      CALL mp_bcast(rd_pos, ionode_id, comm_)
      ! overwrite some settings
      ! if we use 'md', we should benefit most easily from extrapolation
      calculation = 'md'
      ! dynamics not used as positions are integrated by lammps
      ion_dynamics = 'verlet'
      ! force no restart as we cannot know what LAMMPS has done between runs
      restart_mode = 'from_scratch'
      ! reroute secondary output
      outdir = outfile_(5:)
      ! finalize initialization
      CALL iosys()
      if ( gamma_only ) write(stdout, *) "gamma-point specific algorithms are used"
    ! TODO: can we use the stop-mechanism somehow?
      CALL check_stop_init()
      CALL setup()
      CALL init_run()
      ! save number of EC atoms
      nec = c_nec
      ! mark as executed
      first_time = .false.
    END SUBROUTINE start_pw

    SUBROUTINE calc_pw(f_qm, f_ec) BIND(C)
      USE extrapolation, ONLY: update_file
      USE force_mod,     ONLY: force
      USE constants,     ONLY: rytoev, bohr_radius_angs
      IMPLICIT NONE
      REAL(kind=c_double), INTENT(OUT) :: f_qm(3, nat)
      REAL(kind=c_double), INTENT(OUT) :: f_ec(3, nec)
      ! scf calculation
      CALL electrons()
      ! force calculation
      CALL forces()
      ! print coordinates in expected place
      CALL output_tau(.false., .false.)
      ! exchange forces
      f_qm = force * rytoev / bohr_radius_angs
      ! TODO: implement forces on ec atoms
      f_ec = 0
      ! make sure extrapolation works
      CALL update_file()
    END SUBROUTINE calc_pw

    SUBROUTINE update_pw(pos_qm, pos_ec) BIND(c)
      USE extrapolation, ONLY: update_pot
      USE constants,     ONLY: bohr_radius_angs
      USE cell_base,     ONLY: alat, at
      USE ions_base,     ONLY: tau
      USE mp_world,      ONLY: world_comm
      USE mp,            ONLY: mp_bcast
      IMPLICIT NONE
      REAL(kind=c_double), INTENT(IN) :: pos_qm(3, nat)
      REAL(kind=c_double), INTENT(IN) :: pos_ec(3, nec)
      REAL(DP), DIMENSION(3) :: com
      REAL(DP), DIMENSION(3) :: coc
      INTEGER :: i
      ! exchange and re-center positions
      IF (ionode) THEN
        ! center of cell
        coc = MATMUL(at, (/0.5d0, 0.5d0, 0.5d0/))
        ! center of molecule
        com = SUM(pos_qm, dim=2) / (nat * alat * bohr_radius_angs)
        DO i = 1, nat
            tau(:, i) = pos_qm(:, i) / (alat * bohr_radius_angs) + coc - com
        END DO
      END IF
      CALL mp_bcast(tau, ionode_id, world_comm)
      ! update wavefunctions
      CALL update_pot()
      ! re-initialize atomic position-dependent quantities
      CALL hinit1()
    END SUBROUTINE update_pw

    SUBROUTINE end_pw(retval) BIND(C)
      IMPLICIT NONE
      INTEGER (kind=C_INT), INTENT(OUT) :: retval
      CALL laxlib_free_ortho_group()
      CALL stop_run(retval)
      IF (ionode) THEN
        CLOSE(UNIT=55)
      END IF
      stdout = 6
      utilout = 6
    END SUBROUTINE end_pw

END MODULE fix_pw
