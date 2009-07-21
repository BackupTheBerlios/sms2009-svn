C
C
C     PROGRAM EISM1.FOR
C
C
C     This program calculates the evolution of an ice sheet on a
C     bed with constant slope.  The surface mass balance is taken
C     constant.
C
C     Boundary conditions are an ice divide at I = 1.  At the other
C     end, the flotation criterion is used to determine the position
C     of the grounding line.  The profile is extended two gridpoints
C     using the equilibrium profile of a free-floating ice shelf.
C
C     For a description of the model, and discussion of the results,
C     see Applied Glaciology, Section 8.3.
C     This version was used for the EISMINT model intercomparison
C     experiments (Sept. 1997).
C
C
C
C
      DIMENSION H(201),DIFF(201),FLUX(201),DHDT(201)
      DIMENSION HB(201),HS(201)
      DIMENSION UX(201),US(201)
C
C
C---  Input parameters
C
      AFL=1.0E-18
      ROI=910.0
      ROW=1028.0
      G=9.81
C      YRS=31556926.0
C      GYR=YRS*YRS*G
      ROIG=ROI*G
      CONST=((2.0/5.0)*AFL)*(ROIG**3.0)
C
C      
      ALPHA=5.0*0.001
      HBNUL=-250.0
      ACC=0.3
      CFLOAT=1.0-(ROI/ROW)
      C1=(ROW-ROI)/(4.0*ROW)
      CSHELF=(2.6*CONST)*(C1**3.0)
C
C
      DX=2000.0
      DT=0.50
      TIME=0.0
      TMAX=15000.0
      TPR=0.0
      DTPR=250.0
      IMAX=51
      CDIFF=CONST/((2.0*DX)**2.0)
      OPEN(2,FILE='SSHAPENB.CVV',STATUS='NEW')
      OPEN(4,FILE='TSHAPENB.CVV',STATUS='NEW')
      OPEN(3,FILE='SSHASTNB.CVV',STATUS='OLD')
C
C
C---  Initial geometry.
C
      DO 10 I=1,IMAX
        READ(3,1000) DIST,H(I),HS(I),UX(I),US(I),FL,EXX
         DIST=(FLOAT(I)-1.0)*DX
         HB(I)=HBNUL-(ALPHA*DIST)
C         H(I)=570.0
C         HS(I)=HB(I)+H(I)
10    CONTINUE
C
C
C---  START TIME INTEGRATION
C
50    CONTINUE
      WRITE(0,*) TIME,IGR
C
C
C---  Calculate diffusivity.
C
      DO 100 I=2,IMAX 
         IF(H(I).EQ.0.0) THEN
            DIFF(I)=0.0
         ELSE
            DHSDX=HS(I+1)-HS(I-1)
            DIFF(I)=(CDIFF*DHSDX*DHSDX)*(H(I)**5.0)
         ENDIF
100   CONTINUE
      DIFF(1)=0.0
C
C
C---  Calculate ice flux.
C
      DO 150 I=1,IMAX-1
         FLUX(I)=((DIFF(I)+DIFF(I+1))/(2.0*DX))*(HS(I+1)-HS(I))
150   CONTINUE
C
C
C---  Calculate change in ice thickness.
C
      DO 200 I=2,IMAX-1
         DHDT(I)=((FLUX(I)-FLUX(I-1))/DX)+ACC
200   CONTINUE
      DHDT(1)=((2.0*FLUX(1))/DX)+ACC
      DHDT(IMAX)=0.0
C
C
C---  Calculate new ice thickness.
C
      DO 250 I=1,IMAX-1
         H(I)=H(I)+(DHDT(I)*DT)
         HFLOAT=(ROI/ROW)*H(I)
         DEPTH=-HB(I)
         IF(HFLOAT.LE.DEPTH.AND.DEPTH.GE.0.0) THEN
             IGR=I-1
            GOTO 260
         ENDIF
         IF(H(I).LE.0.0) H(I)=0.0
         HS(I)=HB(I)+H(I)
         IGR=I
250   CONTINUE
260   CONTINUE
C
C
C---  Extend ice-shelf profile one gridpoint.
C
265   CONTINUE
      DHSDX=(HS(IGR)-HS(IGR-1))/DX
      UGR=((DHSDX*DHSDX*DHSDX)*(H(IGR)**4.0))*CONST
      UGR=-UGR
C
C
      X1=CSHELF/ACC
      X2=((H(IGR)**4.0)*X1)-1.0
      X3=(UGR*UGR*UGR*UGR)*X2
      X4A=(ACC*DX)+(H(IGR)*UGR)
      X4A=X4A**4.0
      X4B=(ACC*2.0*DX)+(H(IGR)*UGR)
      X4B=X4B**4.0
C
C
      H(IGR+1)=(X1-(X3/X4A))**(-0.25)
      IF(H(IGR+1).GT.0.577*H(IGR)) H(IGR+1)=0.577*H(IGR)
      HS(IGR+1)=CFLOAT*H(IGR+1)
C
C
      HFLOAT=(ROI/ROW)*H(IGR+1)
      DEPTH=-HB(IGR+1)
      IF(HFLOAT.GT.DEPTH) THEN
         HS(IGR+1)=H(IGR+1)+HB(IGR+1)
         IGR=IGR+1
         GOTO 265
      ENDIF
C
C
      DO 270 I=IGR+2,IMAX
         H(I)=0.0
         HS(I)=0.0
270   CONTINUE      
C
C
C---  Write results to file.
C
      TIME=TIME+DT
      IF(TIME.GE.TPR) THEN      
         VOL=0.0
         DO 280 I=1,IGR
            VOL=VOL+H(I)
280      CONTINUE
         VOL=VOL*DX
         UX(21)=-(FLUX(21)+FLUX(20))*(0.5/H(21))
         UX(31)=0.0
         WRITE(4,2000) TIME,VOL,IGR,H(1),UX(21),UX(31)
         TPR=TPR+DTPR
      ENDIF
      IF(TIME.LT.TMAX) GOTO 50
C
C
      DO 290 I=2,IMAX
         IF(H(I).EQ.0.0) THEN
            UX(I)=0.0
            US(I)=0.0
         ELSE
            UX(I)=-(FLUX(I)+FLUX(I-1))*(0.5/H(I))
            US(I)=0.0
         ENDIF
290   CONTINUE
C
C
      DO 300 I=1,IMAX
         DIST=(FLOAT(I)-1.0)*(DX/1000.0)
         HBASE=HS(I)-H(I)
         FL=(UX(I)+US(I))*H(I)
         IF(I.EQ.1) THEN
            EXX=(UX(2)+US(2))/DX
         ELSE
            IF(I.EQ.IMAX) THEN
               EXX=0.0
            ELSE
               EXX=((UX(I+1)+US(I+1))-(UX(I-1)+US(I-1)))/(2.0*DX)
            ENDIF
         ENDIF
         WRITE(2,1000) DIST,H(I),HS(I),UX(I),US(I),FL,EXX
300      CONTINUE
310      CONTINUE
C
C
C
C
1000  FORMAT(5(1X,F8.2),2(1XE14.6))
2000  FORMAT(1X,F8.0,1X,E14.6,1X,I2,3(1X,F8.2))
C
C
C
C
      END
C
C
C
C          