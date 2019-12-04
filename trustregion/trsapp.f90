Subroutine trsapp(N,XOPT,GQ,HQ,DELTA,STEP,CRVMIN)
!
!     TRSAPP tries to find an approximate solution D to the problem
!
!     MIN_{STEP} GQ*(XOPT+STEP) + 0.5*(XOPT+STEP)*HQ*(XOPT+STEP)
!     S.T.       ||STEP||_2 <= DELTA
!
!     The solver also returns CRVMIN:
!
!     CRVMIN = 0 if ||STEP||_2 = DELTA
!     Otherwise, CRVMIN > 0 is the smallest curvature of HQ found
!
!     This routine is taken (without change) from M. J. D. Powell's NEWUOA
!
!     Use double precision type as per this machine
      IMPLICIT REAL(kind=kind(1.d0)) (A-H,O-Z)
      Integer, Parameter :: dp = kind(1.d0)
!     Input/output variables
      Integer, Intent(in) :: N
      Real(kind=dp), Intent(in) :: DELTA
      Real(kind=dp), Intent(in) :: XOPT(N), GQ(N)
      Real(kind=dp), Intent(in) :: HQ(N,N)
      Real(kind=dp), Intent(out) :: STEP(N)
      Real(kind=dp), Intent(out) :: CRVMIN
!f2py integer intent(hide),depend(XOPT) :: N=shape(XOPT,0)
! These variables were inputs to original TRSAPP, but we define as local here
      Real(kind=dp) :: D(N), G(N), HD(N), HS(N)
!     Minor change: initialize some variables zero to avoid compiler warnings
      TEMPB=0.0D0
      TEMPA=0.0D0
      SHS=0.0D0
      SG=0.0D0
      QRED=0.0D0
      GGBEG=0.0D0
      GG=0.0D0
      DD=0.0D0
      BSTEP=0.0D0
!
!     N is the number of variables of a quadratic objective function, Q say.
!     The arguments XOPT, GQ and HQ have their usual meanings,
!       in order to define the current quadratic model Q.
!     DELTA is the trust region radius, and has to be positive.
!     STEP will be set to the calculated trial step.
!     The arrays D, G, HD and HS will be used for working space.
!     CRVMIN will be set to the least curvature of H along the conjugate
!       directions that occur, except that it is set to zero if STEP goes
!       all the way to the trust region boundary.
!
!     The calculation of STEP begins with the truncated conjugate gradient
!     method. If the boundary of the trust region is reached, then further
!     changes to STEP may be made, each one being in the 2D space spanned
!     by the current STEP and the corresponding gradient of Q. Thus STEP
!     should provide a substantial reduction to Q within the trust region.
!
!     Initialization, which includes setting HD to H times XOPT.
!
      HALF=0.5D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(1.0D0)
      DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      DO 10 I=1,N
   10 D(I)=XOPT(I)
      GOTO 170
!
!     Prepare for the first line search.
!
   20 QRED=ZERO
      DD=ZERO
      DO 30 I=1,N
      STEP(I)=ZERO
      HS(I)=ZERO
      G(I)=GQ(I)+HD(I)
      D(I)=-G(I)
   30 DD=DD+D(I)**2
      CRVMIN=ZERO
      IF (DD .EQ. ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
!
!     Calculate the step to the trust region boundary and the product HD.
!
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
      GOTO 170
   50 DHD=ZERO
      DO 60 J=1,N
   60 DHD=DHD+D(J)*HD(J)
!
!     Update CRVMIN and set the step-length ALPHA.
!
      ALPHA=BSTEP
      IF (DHD .GT. ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC .EQ. 1) CRVMIN=TEMP
          CRVMIN=DMIN1(CRVMIN,TEMP)
          ALPHA=DMIN1(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
!
!     Update STEP and HS.
!
      GGSAV=GG
      GG=ZERO
      DO 70 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      HS(I)=HS(I)+ALPHA*HD(I)
   70 GG=GG+(G(I)+HS(I))**2
!
!     Begin another conjugate direction iteration if required.
!
      IF (ALPHA .LT. BSTEP) THEN
          IF (QADD .LE. 1.0D-6*QRED) GOTO 160
!     Note: replacing test with one from TRSBOX, as modified in DFBOLS [H. Zhang, 2010]
          IF (GG .LE. MIN(1.0D-6*GGBEG,1.0D-18)) GOTO 160
!          IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
          IF (ITERC .EQ. ITERMAX) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO 80 I=1,N
          D(I)=TEMP*D(I)-G(I)-HS(I)
          DD=DD+D(I)**2
          DS=DS+D(I)*STEP(I)
   80     SS=SS+STEP(I)**2
          IF (DS .LE. ZERO) GOTO 160
          IF (SS .LT. DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
!
!     Test whether an alternative iteration is required.
!     Note: replacing test with one from TRSBOX, as modified in DFBOLS [H. Zhang, 2010]
!
   90 IF (GG .LE. MIN(1.0D-6*GGBEG,1.0D-18)) GOTO 160
!      IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
      SG=ZERO
      SHS=ZERO
      DO 100 I=1,N
      SG=SG+STEP(I)*G(I)
  100 SHS=SHS+STEP(I)*HS(I)
      SGK=SG+SHS
      ANGTEST=SGK/DSQRT(GG*DELSQ)
      IF (ANGTEST .LE. -0.99D0) GOTO 160
!
!     Begin the alternative iteration by calculating D and HD and some
!     scalar products.
!
      ITERC=ITERC+1
      TEMP=DSQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO 110 I=1,N
  110 D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO 130 I=1,N
      DG=DG+D(I)*G(I)
      DHD=DHD+HD(I)*D(I)
  130 DHS=DHS+HD(I)*STEP(I)
!
!     Seek the value of the angle that minimizes Q.
!
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
      IF (QNEW .LT. QMIN) THEN
          QMIN=QNEW
          ISAVE=I
          TEMPA=QSAV
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=QNEW
      END IF
  140 QSAV=QNEW
      IF (ISAVE .EQ. ZERO) TEMPA=QNEW
      IF (ISAVE .EQ. IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+ANGLE)
!
!     Calculate the new STEP and HS. Then test for convergence.
!
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO 150 I=1,N
      STEP(I)=CTH*STEP(I)+STH*D(I)
      HS(I)=CTH*HS(I)+STH*HD(I)
  150 GG=GG+(G(I)+HS(I))**2
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
!     Note: replacing test with one from TRSBOX, as modified in DFBOLS [H. Zhang, 2010]
!      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01D0) GOTO 90
      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 1.0D-6) GOTO 90
  160 RETURN
!
!     The following instructions act as a subroutine for setting the vector
!     HD to the vector D multiplied by the second derivative matrix of Q.
!     They are called from three different places, which are distinguished
!     by the value of ITERC.
!
  170 DO 180 I=1,N
  180     HD(I)=ZERO
      DO 210 J=1,N
          DO 210 I=1,N
  210         HD(I)=HD(I)+HQ(I,J)*D(J)
      IF (ITERC .EQ. 0) GOTO 20
      IF (ITERC .LE. ITERSW) GOTO 50
      GOTO 120
End Subroutine trsapp
