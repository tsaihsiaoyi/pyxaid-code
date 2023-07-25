! ------------------------------------------------------------------------
! Using new version of VASP_4.6 reads
! mpif90 -assume byterecl -o ovlap_NORM_OS_allk_mpi ovlap_NORM_OS_allk_mpi.f
! wavefunctions as adiabatic basis states.
! Then finds the overlap between this states at adjacent time steps.
! This OVERLAP is multiplied by 2*hbar*dt and "NORMALIZED".
! Output is presented in file "coupling" - normalised overlap
! and billdata - nonnormalised overlap and other intermidiate calculations
!-------------------------------------------------------------------------
! Also
! --------------------------------------------------------------------
! Finds oscillator strength f_nm=C*Mu^2
! where C=2*m*omega_nm/hbar
!      Mu=|<ksi_n(r)|r|ksi_m(r)>|
! which can be writen in momentum representation as:
! Mu=<n|p|m>/omega_nm
! BUT WE do NOT TAKE CARE OF C, working in atomic units (!)
! --------------------------------------------------------------------

PROGRAM MAIN
    USE MPI
    IMPLICIT NONE

    INTEGER :: I, J, K, K1, K2

    !!!CONFIG!!!
    INTEGER :: NBANDMIN, NBANDMAX, NBAND, TOTAL_K, DEBUG
    CHARACTER :: SUFX * 4, TMP_CHAR * 99
    REAL(8) :: POTIM
    !!!CONFIG!!!

    !!!READ LATTICE!!!
    REAL(8) :: A(3, 3), AK(3, 3), B1(3), B2(3), B3(3), B1MAG, B2MAG, B3MAG
    !!!READ LATTICE!!!

    !!!READ WAVECAR!!!
    INTEGER :: IRECLOLD, IRECLNEW, INKPTOLD, INKPTNEW, &
               INBANDOLD, INBANDNEW, ISPINOLD, ISPINNEW
    REAL(8) :: RRECLOLD, RRECLNEW, RNKPTOLD, RNKPTNEW, &
               RNBANDOLD, RNBANDNEW, RSPINOLD, RSPINNEW
    REAL(8) :: EMAX, RTAG
    !!!READ WAVECAR!!!

    !!!CALCULATE MAX G!!!
    REAL(8) :: C = 0.262465831D0
    INTEGER :: NB1MAXA, NB2MAXA, NB3MAXA, &
               NB1MAXB, NB2MAXB, NB3MAXB, &
               NB1MAXC, NB2MAXC, NB3MAXC, &
               NB1MAX, NB2MAX, NB3MAX
    REAL(8) :: PHI12, PHI13, PHI23, SINPHI123, VTMP(3), VMAG
    !!!CALCULATE MAX G!!!

    !!!MPI SPLIT!!!
    INTEGER :: NPROCS, IERR, RANK, COUNT_K, II_INIT, II_FINAL
    !!!MPI SPLIT!!!

    !!!READ K-POINTS!!!
    INTEGER :: INPWK1, INPWK2
    INTEGER, ALLOCATABLE :: XXK1(:), YYK1(:), ZZK1(:), XXK2(:), YYK2(:), ZZK2(:)
    REAL(8) :: RNPWK1, RNPWK2, KPTK1(3), KPTK2(3)
    !!!READ K-POINTS!!!

    !!!READ COEFFICIENTS!!!
    INTEGER :: IC, IBAND
    INTEGER :: IREC1, IREC2
    COMPLEX(4), ALLOCATABLE :: COEF(:)
    COMPLEX(8), ALLOCATABLE :: CK1OLD(:, :)
    COMPLEX(8), ALLOCATABLE :: CK1NEW(:, :)
    COMPLEX(8), ALLOCATABLE :: CK2OLD(:, :)
    COMPLEX(8), ALLOCATABLE :: CK2NEW(:, :)
    !!!READ COEFFICIENTS!!!

    !!!CALCULATE OVERLAP!!!
    INTEGER :: KEQUAL
    REAL(8), ALLOCATABLE :: DDNORM(:, :)
    COMPLEX(8), ALLOCATABLE :: DD(:, :, :)
    !!!CALCULATE OVERLAP!!!

    !!!OUTPUT NAC!!!
    CHARACTER :: K1_CHAR * 40, K2_CHAR * 40, COUP_FILE * 20
    !!!OUTPUT NAC!!!

    !!!PARAMETER!!!
    REAL(8), PARAMETER :: PI = 4.*ATAN(1.)
    REAL(8), PARAMETER :: HBAR = -0.658218 !!!eV·fs!!!
    !!!PARAMETER!!!

    !!!CONFIG!!!
    DEBUG = 0
    CALL GETARG(1, TMP_CHAR)
    READ (TMP_CHAR, *) NBANDMIN
    CALL GETARG(2, TMP_CHAR)
    READ (TMP_CHAR, *) NBANDMAX
    CALL GETARG(3, TMP_CHAR)
    READ (TMP_CHAR, *) POTIM
    CALL GETARG(4, TMP_CHAR)
    READ (TMP_CHAR, *) TOTAL_K
    CALL GETARG(5, SUFX)
    CALL GETARG(6, TMP_CHAR)
    IF (TMP_CHAR .NE. '') THEN
        READ (TMP_CHAR, *) DEBUG
    END IF

    NBAND = NBANDMAX - NBANDMIN + 1
    !!!CONFIG!!!

    !!!READ LATTICE!!!
    OPEN (9999, FILE="OUTCAR")

    DO WHILE (.TRUE.)
        READ (9999, '(A99)') TMP_CHAR
        IF (TMP_CHAR(46:71) .EQ. 'reciprocal lattice vectors') THEN
            DO K = 1, 3
                READ (9999, *) (A(I, K), I=1, 3), (AK(J, K), J=1, 3)
            END DO
            EXIT
        END IF
    END DO
    CLOSE (9999)
    AK = AK * PI * 2.

    DO I = 1, 3
        B1(I) = AK(I, 1)
        B2(I) = AK(I, 2)
        B3(I) = AK(I, 3)
    END DO

    B1MAG = SQRT(B1(1)**2 + B1(2)**2 + B1(3)**2)
    B2MAG = SQRT(B2(1)**2 + B2(2)**2 + B2(3)**2)
    B3MAG = SQRT(B3(1)**2 + B3(2)**2 + B3(3)**2)

    IF (DEBUG .GE. 1) THEN
        WRITE (*, '(A)') "REAL LATTICE = "
        WRITE (*, '(3(3F16.8,/))') A
        WRITE (*, '(A)') "RECIPROCAL LATTICE = "
        WRITE (*, '(3(3F16.8,/))') AK
        WRITE (*, '(A)') 'RECIPROCAL LATTICE VECTOR MAGNITUDES = '
        WRITE (*, '(3F16.8)') B1MAG, B2MAG, B3MAG
        WRITE (*, *)
    END IF
    !!!READ LATTICE!!!

    !!!READ WAVECAR!!!
    OPEN (12, FILE="WAVECAROLD", STATUS="OLD", FORM="UNFORMATTED", &
          ACCESS='DIRECT', RECL=1000000)
    READ (12, REC=1) RRECLOLD, RSPINOLD, RTAG

    IRECLOLD = RRECLOLD

    CLOSE (12)

    OPEN (12, FILE="WAVECAROLD", STATUS="OLD", FORM="UNFORMATTED", &
          ACCESS='DIRECT', RECL=IRECLOLD)
    READ (12, REC=2) RNKPTOLD, RNBANDOLD, EMAX

    ISPINOLD = RSPINOLD
    INKPTOLD = RNKPTOLD
    INBANDOLD = RNBANDOLD

    IF (DEBUG .GE. 1) THEN
        WRITE (*, '(A, I8)') 'IRECLOLD  = ', IRECLOLD
        WRITE (*, '(A, I8)') 'ISPINOLD  = ', ISPINOLD
        WRITE (*, '(A, F8.2)') 'RTAG      = ', RTAG
        WRITE (*, '(A, I8)') 'INKPTOLD  = ', INKPTOLD
        WRITE (*, '(A, I8)') 'INBANDOLD = ', INBANDOLD
        WRITE (*, '(A, F8.2)') 'EMAX      = ', EMAX
        WRITE (*, *)
    END IF

    OPEN (13, FILE="WAVECARNEW", STATUS="OLD", FORM="UNFORMATTED", &
          ACCESS='DIRECT', RECL=1000000)
    READ (13, REC=1) RRECLNEW, RSPINNEW, RTAG

    IRECLNEW = RRECLNEW

    CLOSE (13)

    OPEN (13, FILE="WAVECARNEW", STATUS="OLD", FORM="UNFORMATTED", &
          ACCESS='DIRECT', RECL=IRECLNEW)
    READ (13, REC=2) RNKPTNEW, RNBANDNEW, EMAX

    ISPINNEW = RSPINNEW
    INKPTNEW = RNKPTNEW
    INBANDNEW = RNBANDNEW

    IF (DEBUG .GE. 1) THEN
        WRITE (*, '(A, I8)') 'IRECLNEW  = ', IRECLNEW
        WRITE (*, '(A, I8)') 'ISPINNEW  = ', ISPINNEW
        WRITE (*, '(A, F8.2)') 'RTAG      = ', RTAG
        WRITE (*, '(A, I8)') 'INKPTNEW  = ', INKPTNEW
        WRITE (*, '(A, I8)') 'INBANDNEW = ', INBANDNEW
        WRITE (*, '(A, F8.2)') 'EMAX      = ', EMAX
        WRITE (*, *)
    END IF
    !!!READ WAVECAR!!!

    !!!CALCULATE MAX G!!!
    PHI12 = ACOS((B1(1) * B2(1) + B1(2) * B2(2) + B1(3) * B2(3)) / (B1MAG * B2MAG))
    CALL VCROSS(VTMP, B1, B2)
    VMAG = SQRT(VTMP(1)**2 + VTMP(2)**2 + VTMP(3)**2)
    SINPHI123 = (B3(1) * VTMP(1) + B3(2) * VTMP(2) + B3(3) * VTMP(3)) / (VMAG * B3MAG)
    NB1MAXA = (SQRT(EMAX * C) / (B1MAG * ABS(SIN(PHI12)))) + 1
    NB2MAXA = (SQRT(EMAX * C) / (B2MAG * ABS(SIN(PHI12)))) + 1
    NB3MAXA = (SQRT(EMAX * C) / (B3MAG * ABS(SINPHI123))) + 1

    PHI13 = ACOS((B1(1) * B3(1) + B1(2) * B3(2) + B1(3) * B3(3)) / (B1MAG * B3MAG))
    CALL VCROSS(VTMP, B1, B3)
    VMAG = SQRT(VTMP(1)**2 + VTMP(2)**2 + VTMP(3)**2)
    SINPHI123 = (B2(1) * VTMP(1) + B2(2) * VTMP(2) + B2(3) * VTMP(3)) / (VMAG * B2MAG)
    NB1MAXB = (SQRT(EMAX * C) / (B1MAG * ABS(SIN(PHI13)))) + 1
    NB2MAXB = (SQRT(EMAX * C) / (B2MAG * ABS(SINPHI123))) + 1
    NB3MAXB = (SQRT(EMAX * C) / (B3MAG * ABS(SIN(PHI13)))) + 1

    PHI23 = ACOS((B2(1) * B3(1) + B2(2) * B3(2) + B2(3) * B3(3)) / (B2MAG * B3MAG))
    CALL VCROSS(VTMP, B2, B3)
    VMAG = SQRT(VTMP(1)**2 + VTMP(2)**2 + VTMP(3)**2)
    SINPHI123 = (B1(1) * VTMP(1) + B1(2) * VTMP(2) + B1(3) * VTMP(3)) / (VMAG * B1MAG)
    NB1MAXC = (SQRT(EMAX * C) / (B1MAG * ABS(SINPHI123))) + 1
    NB2MAXC = (SQRT(EMAX * C) / (B2MAG * ABS(SIN(PHI23)))) + 1
    NB3MAXC = (SQRT(EMAX * C) / (B3MAG * ABS(SIN(PHI23)))) + 1

    NB1MAX = MAX(NB1MAXA, NB1MAXB, NB1MAXC)
    NB2MAX = MAX(NB2MAXA, NB2MAXB, NB2MAXC)
    NB3MAX = MAX(NB3MAXA, NB3MAXB, NB3MAXC)

    IF (DEBUG .GE. 1) THEN
        WRITE (*, '(A,3I8)') 'MAX G1, G2, G3 = ', NB1MAX, NB2MAX, NB3MAX
        WRITE (*, *)
    END IF
    !!!CALCULATE MAX G!!!

    !!!MPI SPLIT!!!
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)
    COUNT_K = TOTAL_K / NPROCS
    II_INIT = RANK * COUNT_K + 1
    II_FINAL = II_INIT + COUNT_K - 1
    IF (RANK .EQ. (NPROCS - 1)) II_FINAL = TOTAL_K
    !!!MPI SPLIT!!!

    DO K1 = II_INIT, II_FINAL
        DO K2 = 1, TOTAL_K

            !!!READ K-POINTS!!!
            READ (12, REC=3 + (K1 - 1) * (INBANDOLD + 1)) RNPWK1, KPTK1(1:3)
            INPWK1 = RNPWK1

            IF (DEBUG .GE. 1) THEN
                WRITE (*, "(A35,I8)") 'NUMBER OF PLANE WAVES NPL FOR K1 = ', INPWK1
                WRITE (*, "(A35,3F8.2)") 'WAVECAROLD FOR K1: K POINT = ', KPTK1
                WRITE (*, *)
            END IF

            READ (12, REC=3 + (K2 - 1) * (INBANDOLD + 1)) RNPWK2, KPTK2(1:3)
            INPWK2 = RNPWK2

            ALLOCATE (XXK1(INPWK1))
            ALLOCATE (YYK1(INPWK1))
            ALLOCATE (ZZK1(INPWK1))
            ALLOCATE (XXK2(INPWK2))
            ALLOCATE (YYK2(INPWK2))
            ALLOCATE (ZZK2(INPWK2))

            IF (DEBUG .GE. 1) THEN
                WRITE (*, "(A35,I8)") 'NUMBER OF PLANE WAVES NPL FOR K2 = ', INPWK2
                WRITE (*, "(A35,3F8.2)") 'WAVECAROLD FOR K2: K POINT = ', KPTK2
            END IF
            !!!READ K-POINTS!!!

            !!!FIND G!!!
            CALL FIND_G(B1, B2, B3, KPTK1, NB1MAX, NB2MAX, NB3MAX, &
                        EMAX, XXK1, YYK1, ZZK1)
            CALL FIND_G(B1, B2, B3, KPTK2, NB1MAX, NB2MAX, NB3MAX, &
                        EMAX, XXK2, YYK2, ZZK2)
            !!!FIND G!!!

            ALLOCATE (COEF(MAX(INPWK1, INPWK2)))

            ALLOCATE (CK1OLD(INPWK1, NBAND))
            ALLOCATE (CK2OLD(INPWK2, NBAND))
            ALLOCATE (CK1NEW(INPWK1, NBAND))
            ALLOCATE (CK2NEW(INPWK2, NBAND))

            !!!READ COEFFICIENTS!!!
            IREC1 = 4 + (K1 - 1) * (INBANDOLD + 1)
            IC = 1
            DO IBAND = 1, INBANDOLD
                READ (12, REC=IREC1) (COEF(I), I=1, INPWK1)
                IREC1 = IREC1 + 1
                IF (IBAND .GE. NBANDMIN .AND. IBAND .LE. NBANDMAX) THEN
                    DO I = 1, INPWK1
                        CK1OLD(I, IC) = COEF(I)
                    END DO
                    IC = IC + 1
                END IF
            END DO

            IREC2 = 4 + (K2 - 1) * (INBANDOLD + 1)
            IC = 1
            DO IBAND = 1, INBANDOLD
                READ (12, REC=IREC2) (COEF(I), I=1, INPWK2)
                IREC2 = IREC2 + 1
                IF (IBAND .GE. NBANDMIN .AND. IBAND .LE. NBANDMAX) THEN
                    DO I = 1, INPWK2
                        CK2OLD(I, IC) = COEF(I)
                    END DO
                    IC = IC + 1
                END IF
            END DO

            CLOSE (12)

            IREC1 = 4 + (K1 - 1) * (INBANDNEW + 1)
            IC = 1
            DO IBAND = 1, INBANDNEW
                READ (13, REC=IREC1) (COEF(I), I=1, INPWK1)
                IREC1 = IREC1 + 1
                IF (IBAND .GE. NBANDMIN .AND. IBAND .LE. NBANDMAX) THEN
                    DO I = 1, INPWK1
                        CK1NEW(I, IC) = COEF(I)
                    END DO
                    IC = IC + 1
                END IF
            END DO

            IREC2 = 4 + (K2 - 1) * (INBANDNEW + 1)
            IC = 1
            DO IBAND = 1, INBANDNEW
                READ (13, REC=IREC2) (COEF(I), I=1, INPWK2)
                IREC2 = IREC2 + 1
                IF (IBAND .GE. NBANDMIN .AND. IBAND .LE. NBANDMAX) THEN
                    DO I = 1, INPWK2
                        CK2NEW(I, IC) = COEF(I)
                    END DO
                    IC = IC + 1
                END IF
            END DO

            CLOSE (13)

            DEALLOCATE (COEF)
            !!!READ COEFFICIENTS!!!

            !!!CALCULATE OVERLAP!!!
            ALLOCATE (DD(NBAND, NBAND, 2))
            ALLOCATE (DDNORM(NBAND, 4))

            DD = (0., 0.)
            DDNORM = 0.

            DO K = 1, INPWK1
                KEQUAL = 0
                CALL FIND_KEQUAL(K, KEQUAL, XXK1, YYK1, ZZK1, XXK2, YYK2, ZZK2, INPWK2)
                IF (KEQUAL .NE. 0) THEN
                    DO I = 1, NBAND
                        DO J = 1, NBAND
                            DD(I, J, 1) = DD(I, J, 1) + &
                                          CONJG(CK1OLD(K, I)) * CK2NEW(KEQUAL, J)
                            DD(I, J, 2) = DD(I, J, 2) - &
                                          CONJG(CK1NEW(K, I)) * CK2OLD(KEQUAL, J)
                        END DO
                        DDNORM(I, 1) = DDNORM(I, 1) + &
                                       CONJG(CK1OLD(K, I)) * CK1OLD(K, I)
                        DDNORM(I, 2) = DDNORM(I, 2) + &
                                       CONJG(CK2OLD(KEQUAL, I)) * CK2OLD(KEQUAL, I)
                        DDNORM(I, 3) = DDNORM(I, 3) + &
                                       CONJG(CK1NEW(K, I)) * CK1NEW(K, I)
                        DDNORM(I, 4) = DDNORM(I, 4) + &
                                       CONJG(CK2NEW(KEQUAL, I)) * CK2NEW(KEQUAL, I)
                    END DO
                END IF
            END DO
            DO I = 1, NBAND
                DO j = 1, NBAND
                    DD(I, J, 1) = (HBAR / (2 * POTIM)) * &
                                  (DD(I, J, 1) / SQRT(DDNORM(I, 1) * DDNORM(J, 4)) + &
                                   DD(I, J, 2) / SQRT(DDNORM(I, 3) * DDNORM(J, 2)))
                END DO
            END DO
            !!!CALCULATE OVERLAP!!!

            !!!OUTPUT NAC!!!
            WRITE (K1_CHAR, '(I0)') K1
            WRITE (K2_CHAR, '(I0)') K2

            COUP_FILE = 'real'//SUFX//'_'//TRIM(K1_CHAR)//'_'//TRIM(K2_CHAR)
            OPEN (141, FILE=COUP_FILE)

            DO I = 1, NBAND
                DO J = 1, NBAND
                    WRITE (141, '(2E10.2,1X)', ADVANCE='NO') DD(I, J, 1)
                END DO
                WRITE (141, *)
            END DO

            CLOSE (141)

            DEALLOCATE (CK1OLD, CK2OLD, CK1NEW, CK2NEW, DD, &
                        DDNORM, XXK1, YYK1, ZZK1, XXK2, YYK2, ZZK2)
            !!!OUTPUT NAC!!!

        END DO
    END DO

    CALL MPI_FINALIZE(ierr)

END PROGRAM

SUBROUTINE FIND_G(B1, B2, B3, KPT, NB1MAX, NB2MAX, NB3MAX, EMAX, XX, YY, ZZ)
    IMPLICIT NONE

    REAL(8) :: B1(3), B2(3), B3(3)
    REAL(8) :: ETOT, GTOT, EMAX
    REAL(8) :: VTMP(3), SUMKG(3), KPT(3)

    INTEGER :: XX(*), YY(*), ZZ(*)
    INTEGER :: N5, I, J, K, I1, I2, I3
    INTEGER :: NB1MAX, NB2MAX, NB3MAX
    INTEGER :: IG1P, IG2P, IG3P

    REAL(8) :: C = 0.262465831D0
    !!!C=2m/hbar^2 in Angstrom^2/eV unit!!!

    N5 = 0

    DO I3 = 0, 2 * NB3MAX
        IG3P = I3
        IF (I3 .GT. NB3MAX) IG3P = I3 - 2 * NB3MAX - 1
        DO I2 = 0, 2 * NB2MAX
            IG2P = I2
            IF (I2 .GT. NB2MAX) IG2P = I2 - 2 * NB2MAX - 1
            DO I1 = 0, 2 * NB1MAX
                IG1P = I1
                IF (I1 .GT. NB1MAX) IG1P = I1 - 2 * NB1MAX - 1
                DO J = 1, 3
                    SUMKG(J) = (KPT(1) + IG1P) * B1(J) + &
                               (KPT(2) + IG2P) * B2(J) + &
                               (KPT(3) + IG3P) * B3(J)
                END DO
                GTOT = SQRT(SUMKG(1)**2 + SUMKG(2)**2 + SUMKG(3)**2)
                ETOT = GTOT**2 / C
                IF (ETOT .LT. EMAX) THEN
                    N5 = N5 + 1
                    XX(N5) = IG1P
                    YY(N5) = IG2P
                    ZZ(N5) = IG3P
                END IF
            END DO
        END DO
    END DO

END SUBROUTINE

SUBROUTINE VCROSS(A, B, C)
    REAL(8) :: A(3), B(3), C(3)
    A(1) = B(2) * C(3) - B(3) * C(2)
    A(2) = B(3) * C(1) - B(1) * C(3)
    A(3) = B(1) * C(2) - B(2) * C(1)
END SUBROUTINE VCROSS

SUBROUTINE FIND_KEQUAL(K, KEQUAL, XXK1, YYK1, ZZK1, XXK2, YYK2, ZZK2, NPWK)
    IMPLICIT NONE
    INTEGER I, K, KEQUAL, NPWK
    INTEGER XXK1(*), YYK1(*), ZZK1(*)
    INTEGER XXK2(*), YYK2(*), ZZK2(*)

    KEQUAL = 0
    DO I = 1, NPWK
        IF (XXK2(I) .EQ. XXK1(K) .AND. &
            YYK2(I) .EQ. YYK1(K) .AND. &
            ZZK2(I) .EQ. ZZK1(K)) THEN
            KEQUAL = I
            CYCLE
        END IF
    END DO
END SUBROUTINE FIND_KEQUAL
