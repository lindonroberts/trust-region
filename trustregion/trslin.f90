Subroutine trslin(N,GOPT,SL,SU,DELTA,D,CRVMIN)
!
!     TRSLIN finds a solution D to the problem
!
!     MIN_{D} GOPT*D
!     S.T.    ||D||_2 <= DELTA
!             SL <= D <= SU
!
!     The solver also returns CRVMIN:
!
!     CRVMIN = 0 if ||D||_2 = DELTA
!     CRVMIN = -1 if every step of the search hit the SL/SU bounds
!
!     This routine is based on Algorithm B.1 in
!     L. Roberts (2019), Derivative-Free Algorithms for Nonlinear Optimisation 
!     Problems, PhD Thesis, University of Oxford
!
!     Use double precision type as per this machine
      IMPLICIT REAL(kind=kind(1.d0)) (A-H,O-Z)
      Integer, Parameter :: dp = kind(1.d0)
!     Input/output variables
      Integer, Intent(in) :: N
      Real(kind=dp), Intent(in) :: DELTA
      Real(kind=dp), Intent(in) :: GOPT(N), SL(N), SU(N)
      Real(kind=dp), Intent(out) :: D(N)
      Real(kind=dp), Intent(out) :: CRVMIN
!f2py integer intent(hide),depend(GOPT) :: N=shape(GOPT,0)
!     Local vectors
      Real(kind=dp) :: DIRN(N)
      Integer :: BDRY(N)
!     BDRY(I)=0 if unconstrained, 1 if constrained by SL/SU

!     Zero threshold value
      ZERO=0.0D0
      ONEMIN=-1.0D0
!     We will always hit ||D||=DELTA unless on box bdry,
!     so set CRVMIN=0 and only change to -1 if stuck in box (below)
      CRVMIN=ZERO
      EPS=1.0D-14
      EPSSQ=EPS**2
!     Never search in directions where GOPT(I)=0
      DO I=1,N
          IF (ABS(GOPT(I)) < EPS) THEN
              DIRN(I)=ZERO
              BDRY(I)=1
          ELSE
              DIRN(I)=-GOPT(I)
              BDRY(I)=0
          END IF
      END DO
!
!     Main loop
!
      DO I=1,N
!         Terminate if ||DIRN||=0 (i.e. on box bdry)
          DIRNSQ=ZERO
          DO J=1,N
              DIRNSQ=DIRNSQ+DIRN(J)**2
          END DO
          IF (DIRNSQ<EPSSQ) THEN
              CRVMIN=ONEMIN
              RETURN
          END IF
!         Terminate if on ball boundary
          DSQ=ZERO
          DO J=1,N
              DSQ=DSQ+D(J)**2
          END DO
          IF (DSQ>=DELTA**2-EPSSQ) THEN
              CRVMIN=ZERO
              RETURN
          END IF
!         Unconstrained step length
!         Solve ||D+ALPHA1*DIRN||^2=DELTA^2 s.t. ALPHA1>=0
!         Include some checks to ensure feasibility, SQRT>=0
          DDIRN=ZERO
          DO J=1,N
              DDIRN=DDIRN+DIRN(J)*D(J)
          END DO
          TMP=MAX(DDIRN**2 + DIRNSQ*(DELTA**2-DSQ),ZERO)
          ALPHA1=(SQRT(TMP)-DDIRN)/DIRNSQ
          ALPHA1=MAX(ALPHA1,ZERO)
!         Check for box boundary
          IDXBOX=0
          STEP=ZERO
          DO J=1,N
!             Only check currently unconstrained directions
              IF (BDRY(J)==0) THEN
                  STEPJ=D(J)+ALPHA1*DIRN(J)
                  IF ((STEPJ<=SL(J)+EPS) .OR. (STEPJ>=SU(J)-EPS)) THEN
                      IDXBOX=J
                      STEP=STEPJ
                  END IF
              END IF
          END DO
!         IDXBOX=0 if unconstrained by box, otherwise project
          IF (IDXBOX==0) THEN
              DO J=1,N
                  D(J)=D(J)+ALPHA1*DIRN(J)
              END DO
              CRVMIN=ZERO
              RETURN
          END IF
!         Project into box along coordinate IDXBOX
          BDRY(IDXBOX)=1
          ALPHA2=ZERO
          IF (STEP<=SL(IDXBOX)+EPS) THEN
              ALPHA2=(SL(IDXBOX)-D(IDXBOX))/DIRN(IDXBOX)
          ELSE
              ALPHA2=(SU(IDXBOX)-D(IDXBOX))/DIRN(IDXBOX)
          END IF
!         Take step
          DO J=1,N
              D(J)=D(J)+ALPHA2*DIRN(J)
          END DO
          DIRN(IDXBOX)=ZERO
      END DO
!     If here, every coordinate has been fixed
      CRVMIN=ONEMIN
      RETURN
End Subroutine trslin
