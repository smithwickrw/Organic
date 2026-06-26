C     QCP5.FOR   Modified
      PROGRAM MAIN                                                      
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N) 
      COMMON/ARRAYS/A(800,800),AX(800),B(800,800),BX(800)
      COMMON/HCORE/H(800,800)
      COMMON/LAMMA/G(200,200)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SMK/C1,C2,C3,C4,C5,C6,C7,C8,C9,CA,CB,CC,CD,CE,CF,CG,CH,CI,
     1 CJ,CK,CL,CM,CN,CO,CP1,CP2,CP3,CQ1,CQ2,CQ3,CR1,CR2,CS,CHHH
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      COMMON/BOB/UZ
      DIMENSION ANN(20),QQ1(200),NPSAV(200),DSAV(200)
      DIMENSION INN(12),INO(12),Q1(200,2),Q2(200,2)
      INTEGER U,OCCA,OCCB,AN,CHARGE,EL,ORB
      INTEGER*4 IPROD,NUMBIG
      CHARACTER*4 IOP1,IOP2,ANN,FLAG
      CHARACTER*1 STAR
      CHARACTER(8) TODAY
      CHARACTER(11) TIMENOW1,TIMENOW2
      DATA A0/.529167/,IRR/10/,IW/7/
      WRITE(*,189)
  189 FORMAT(1X,'START:')
      OPEN(8,FILE='FORT8')
      CALL DATE_AND_TIME(DATE=TODAY,TIME=TIMENOW1)
      C1= 0.3126 
      C2= 0.0628 
      C3= 0.0631 
      C4= 0.0273 
      C6= 0.0453 
      C7= 0.0193 
      C8= 0.0212 
      CA= 0.0363 
      CB= 0.0560 
      CC= 0.0383 
      CHHH=0.
      CD= -0.436
      CE= -0.311
      CF= -0.230
C   CG IS NOT USED
      CH=  0.304
      CI= -0.426
      CJ= -0.172
C   CK IS NOT USED
      CL=  0.362
      CM= -0.327
      CN=  0.542
C   CO IS NOT USED
      CP1= 0.3997 
      CP2= 0.3967 
      CP3= 0.4549 
      CQ1= 0.2928 
      CQ2= 0.3613 
      CQ3= 0.4593 
      CR1= 0.3505 
      CR2= 0.4052 
      CS=  0.4122 
      SUMIN=0.15E30
      SUMO= SUMIN
      ILM= 1
      ILL= 1        
      SUM=0.
        IPP=0
      IPOS=1
        INEG=1
        IPROD=0
        IF(ILM.EQ.1) DEL=0.00
 191    FORMAT(1X,F7.3)
        NPASS=0.
        SUMO=SUMIN
        ICOUNT=0
        IRW=0
        DELOLD=DEL
C        OPEN(8,FILE='FORT8')
 1      CONTINUE
        UZ=0.
        DEL=0.
        NPT=0
        INEG=1
        IPOS=1
        IPROD=0
22      CONTINUE
        IF(SUMIN.LT.SUMO) SUMO=SUMIN
        IPP=IPP+1      
        OPEN(IRR,FILE='INPUT2')
        NDPC=0
C     READ TITLE CARD
        SUM= 0.
        SUMA=0.
        SUMB=0.
        SUMC=0.
      NCOUNT=0
      ICONT=0
      ICONT1=0
      ICONT2=0
      ICONT3=0
      ICONT4=0
      ICONT5=0
      ICONT6=0
   15 CONTINUE
C  EACH MOLECULE NEEDS TO START AT 0. AND BE SELF-CONSISTENT IN CHARGE:
      DO 234 NN=1,200
      QQ1(NN)=0.
      CHRG(NN)=0.
      C(NN,1)=0.
      C(NN,2)=0.
      C(NN,3)=0.
      STAR(NN)=' '
 234  CONTINUE   
      IDOUBLE=0
      IPASS=0
      NCOUNT=NCOUNT+1     
C      IF(DEL.EQ.0.) OPEN(7,FILE='OUTPUT')
C      IF(DEL.EQ.0.) OPEN(9,FILE='FORTDM')
  991 FORMAT(1X,'NCOUNT= ',I3)
      READ(IRR,20) (ANN(I),I=1,20)                                      
      IF(ANN(1).EQ.'FINI') GO TO 100                                        
      IF(DEL.EQ.0.) WRITE(8,991) NCOUNT
C      IF(DEL.EQ.0.) WRITE(7,30) (ANN(I),I=1,18)                     
      IF(DEL.EQ.0.) WRITE(8,30) (ANN(I),I=1,18)
C     READ OPTIONS  (CNDO,INDO)  (CLSD,OPEN,NMR )                   
      READ(IRR,21) IOP1,IOP2                                        
      READ(IRR,50) NATOMS,CHARGE,MULTIP
C      WRITE(8,60) NATOMS,CHARGE,MULTIP                
C      WRITE(8,887) MULTIP,IOP2
887   FORMAT(1X,'MULTIPLICITY IS ',I2,2X,'(',A4,')')
C      IF(IOP2.NE.'NMR ') GO TO 14                                    
C      KATOM=MULTIP                                                      
C      MULTIP=1                                                        
   14 IF(MULTIP.GT.0) GO TO 16                                        
      MULTIP=1                                                        
   16 CONTINUE                                                       
C     READ IN THE ATOMIC NUMBERS AND COORDINATES                       
      DO 10 I = 1,NATOMS                                              
      READ(IRR,70) AN(I),STAR(I),C(I,1),C(I,2),C(I,3)                 
      IF(IOP1.EQ.'INDO'.AND.AN(I).GT.10) GO TO 100                    
      GO TO 12                                                       
   11 CONTINUE
C   WRITE(IW,80)                                                     
      CLOSE(8)
      STOP
   12 CONTINUE
      C(I,1) = C(I,1)/A0
      C(I,2) = C(I,2)/A0
      C(I,3) = C(I,3)/A0
   10 CONTINUE
      RYE=0.
      DP=0.
      XDE=0.
      XDM=0.
C     READ EXPERIMENTAL DISSOCIATION ENERGY AND DIPOLE MOMENT
      READ(IRR,215) XDE,XDM
      IF(XDM.NE.0..AND.NATOMS.EQ.2) NDPC=NDPC+1
C      IF(DEL.EQ.0.) WRITE(3,334) NCOUNT,(ANN(K),K=1,6),XDE,XDM
 334  FORMAT(1X,I3,1X,6A4,1X,'XDE=',F8.4,1X,'XDM=',F8.4,/)
C      IF(DEL.EQ.0.) WRITE(8,216) XDE,XDM,Q1(1),Q1(1)     
215   FORMAT(2F10.4)
216   FORMAT(1X,'EXPERIMENTAL ATOMIZATION ENERGY=',F9.4,2X,
     1 'EXPERIMENTAL DIPOLE MOMENT=',F8.4)                           
217   FORMAT(11X,'AN(2)= ',I2,3X,'XYZ COORDINATES= ',3F10.4)
C     CALCULATE THE NUMBER OF ELECTRONS, BASIS FUNCTIONS AND OCCUPIED   
C     ORBITALS (OCCA AND OCCB)                                          
      N=0                                                               
      NELECS= -CHARGE                                                   
      DO 17 I=1,NATOMS                                                  
      K=AN(I)                                                           
      N=N+NORB(K)                                                       
      NELECS=NELECS+INT(CZ(K)+.01)                                      
   17 CONTINUE                                                          
      NP=2
   18 CONTINUE
      IDOUBLE=IDOUBLE+1
C   UP TO 20 PASSES ARE MADE FOR MOLECULES TO DETERMINE Q1:
      IF(IDOUBLE.EQ.1.AND.DEL.EQ.0.) IPASS=1
      IF(IDOUBLE.EQ.1.AND.DEL.EQ.0.) DEL=0.0001
      IF(IDOUBLE.GE.NP.AND.IPASS.EQ.1) DEL=0.
      IF(IDOUBLE.GE.NP.AND.IPASS.EQ.1) IPASS=0
      OCCA=(NELECS+MULTIP-1)/2                                          
      OCCB=(NELECS-MULTIP+1)/2                                         
      CALL INTGRL                                                    
      CALL GUESSP(IOP1)                                           
      IF (IOP2.NE.'CLSD') GO TO 91                                      
      CALL SCFC(IOP1)                                           
      GO TO 95                                                         
   91 IF (IOP2.NE.'OPEN') GO TO 95                                   
      CALL SCFO(IOP1)                                               
C      GO TO 95                                                   
C   92 IF (IOP2.NE.'NMR ') GO TO 100                                  
C      CALL SCFN(IOP1,KATOM)                                        
C      GO TO 15                                                       
   95 CONTINUE
      CALL CPRINT(IOP1,IOP2)
      TOTS=0.
      SUMERR=0.
      DO 196 I=1,NATOMS
      IF(NATOMS.EQ.2.AND.AN(1).EQ.AN(2)) CHRG(I)=0.
      IF(CHRG(I).GT.0.) TOTS=TOTS+CHRG(I)
 196  CONTINUE         
      QQ1(IDOUBLE)=TOTS
      IF(IDOUBLE.GE.2) DELTA=QQ1(IDOUBLE)-QQ1(IDOUBLE-1)
      IF(IDOUBLE.GE.2.AND.ABS(DELTA).GT.0.001.AND.NP.GT.1.AND.NP.LT.20)
     .    NP=NP+1
      IF(IDOUBLE.LT.NP) NPSAV(NCOUNT)=NP
C      IF(IDOUBLE.LT.NP) DSAV(NCOUNT)=DELTA
      IF(IDOUBLE.LT.NP) GO TO 18
      IF(DEL.EQ.0.) WRITE(8,226) UZ 
      DO 776 J=1,NATOMS
      IF(DEL.EQ.0.) WRITE(8,777) J,AN(J),STAR(J),CHRG(J)
  776 CONTINUE
  777 FORMAT(1X,'ATOM#',I3,2X,'ELEMENT#=',I2,1X,A1,' CHARGE=',F8.4)  
      NPT=NPT+NP 
      IF(XDE.GT.0.) SUM= SUM+ (RYE-XDE)*(RYE-XDE)
      IF(XDE.GT.0.) SUMA= SUMA+ (RYE-XDE)*(RYE-XDE)
      IF(XDE.GT.0.) SUMERR= (RYE-XDE)*(RYE-XDE)
      IF(XDE.NE.0..AND.NATOMS.EQ.2.AND.CHRG(1).LE.0.) DP=  ABS(DP)
      IF(XDE.NE.0..AND.NATOMS.EQ.2.AND.CHRG(1).GT.0.) DP= -ABS(DP)      
      IF(XDE.GT.0..AND.XDM.NE.0..AND.NATOMS.EQ.2) SUM=SUM+
     1 10.*(DP-XDM)*(DP-XDM) 
      IF(XDE.GT.0..AND.XDM.NE.0..AND.NATOMS.EQ.2) SUMB=SUMB+
     1 10.*(DP-XDM)*(DP-XDM) 
      IF(XDE.GT.0..AND.XDM.NE.0.0.AND.NATOMS.EQ.2) SUMERR=
     1 (RYE-XDE)*(RYE-XDE) +10.*(DP-XDM)*(DP-XDM)  
677    FORMAT(1X,'CARBON HYBRIDIZATION IS ',A4)
      IF(DEL.EQ.0..AND.CHARGE.EQ.0) WRITE(8,104) -EUA(NATOMS), 
     1 -EUA(NATOMS)-RYE
 104  FORMAT(1X,'TOTAL ELECTRON ENERGY OF UNBOUND ATOMS=',F13.3,
     1 ' eV; OF MOLECULE=',F13.3,' eV')
      IF(DEL.EQ.0.) WRITE(8,105) RYE,DP     
 226  FORMAT(1X,'TOTAL REPULSION WORK FUNCTION, UZ= ',F9.4,2X,'eV')   
C      IF(DEL.EQ.0..AND.DP.LT.0..AND.XDM.GT.0.) XDM= -XDM
      IF(DEL.EQ.0.) WRITE(8,216) XDE,XDM
      IF(DEL.EQ.0.) WRITE(8,219) NPSAV(NCOUNT), SUMERR
      IF(DEL.EQ.0..AND.SUMERR.GT.2.25) WRITE(8,277)
      IF(DEL.EQ.0..AND.SUMERR.GT.2.25) GO TO 311
      IF(DEL.EQ.0..AND.SUMERR.GT.1.5625) WRITE(8,278)
      IF(DEL.EQ.0..AND.SUMERR.GT.1.5625) GO TO 311
      IF(DEL.EQ.0..AND.SUMERR.GT.1.) WRITE(8,279)
      IF(DEL.EQ.0..AND.SUMERR.GT.1.) GO TO 311
      IF(DEL.EQ.0..AND.SUMERR.GT.0.5625) WRITE(8,287)
      IF(DEL.EQ.0..AND.SUMERR.GT.0.5625) GO TO 311
      IF(DEL.EQ.0..AND.SUMERR.GT.0.25) WRITE(8,288)
      IF(DEL.EQ.0..AND.SUMERR.GT.0.25) GO TO 311
      IF(DEL.EQ.0..AND.SUMERR.GT.0.0625) WRITE(8,289)
      IF(DEL.EQ.0..AND.SUMERR.GT.0.0625) GO TO 311
      IF(DEL.EQ.0..AND.SUMERR.LE.0.0625) WRITE(8,290)
 311  CONTINUE
      IF(DEL.EQ.0..AND.SUMERR.LT.0.0625) ICONT=
     1 ICONT+1
      IF(DEL.EQ.0..AND.SUMERR.GT.0.0625.AND.SUMERR.LE.0.25) 
     1 ICONT1=ICONT1+1
      IF(DEL.EQ.0..AND.SUMERR.GT.0.25.AND.SUMERR.LE.0.5625) 
     1 ICONT2=ICONT2+1
      IF(DEL.EQ.0..AND.SUMERR.GT.0.5625.AND.SUMERR.LE.1.) 
     1 ICONT3=ICONT3+1
      IF(DEL.EQ.0..AND.SUMERR.GT.1..AND.SUMERR.LE.1.5625) 
     1 ICONT4=ICONT4+1
      IF(DEL.EQ.0..AND.SUMERR.GT.1.5625.AND.SUMERR.LE.2.25) 
     1 ICONT5=ICONT5+1
      IF(DEL.EQ.0..AND.SUMERR.GT.2.25) ICONT6=ICONT6+1      
 277  FORMAT(31X,'************** SS ERROR >2.25 ***********')       
 278  FORMAT(40X,'************** SS ERROR >1.5625 *********')      
 279  FORMAT(40X,'************** SS ERROR >1. <1.5625 ******')
 287  FORMAT(40X,'************** SS ERROR >0.5625 <1. *****')
 288  FORMAT(40X,'************** SS ERROR >0.25 <0.5625 ***')      
 289  FORMAT(40X,'************** SS ERROR >0.0625 <0.25 ***')
 290  FORMAT(40X,'************** SS ERROR <0.0625 ********')
 219  FORMAT(1X,'# OF ITERATIONS=',I3,3X,
     1 'SUM SQUARED ERROR OF MOLECULE=',F9.5)
C      IF(DEL.EQ.0..AND.XDE.NE.0.AND.XDM.NE.0..AND.NATOMS.EQ.2) WRITE(9,
C     1 905) NCOUNT,XDM,DP
      SUMC=SUMA +SUMB
  99  GO TO 15                                                          
C    'FINI' OF INPUT (AT INDO 103) SENDS THE PROGRAM TO STATEMENT 100:
C    'FINI' MUST HAVE CAPITAL LETTERS
  100 CONTINUE
C905   FORMAT(1X,I5,5X,2F8.4)                                                
      CLOSE(6)
C  MODIFIED HERE:
      IF(DEL.EQ.0..AND.SUMC.LE.SUMIN) SUMIN=SUMC
      IF(DEL.EQ.0.) WRITE(8,27)
       IF(DEL.EQ.0.) WRITE(8,725) C1,C2,C3,C4,C6,C7,C8,CA,CB,CC,0.,CP1,
     1 CP2,CP3,CQ1,CQ2,CQ3,CR1,CR2, CS,CHHH,CD,CE,CF,CH,CI,CJ,CL,CM,
     2 CN,SUMA,SUMB, SUMIN,SUMIN,SUMIN
 725  FORMAT(3X,10F7.4,/,3X,10F7.4,/,3X,10F7.3,/,25X,4F8.3,F12.4,/)  
C      IF(DEL.EQ.0.) WRITE(8,888) NDPC
C   2 MOLECULES WERE CALCULATED BUT NOT USED TO CALCULATE ERRORS:
       IF(DEL.EQ.0.) WRITE(8,312) ICONT-2,ICONT1,ICONT2,ICONT3,ICONT4,
     1 ICONT5
C  PERCENTAGES FOR 117 MOLECULES (+4 MOLECULES ARE "OUT"):
       IF(DEL.EQ.0.) WRITE(8,313) FLOAT(ICONT-4)/1.17
       IF(DEL.EQ.0.) WRITE(8,314) FLOAT(ICONT1+ICONT-4)/1.17
       IF(DEL.EQ.0.) WRITE(8,315) FLOAT(ICONT2+ICONT1+ICONT-4)/1.17
       IF(DEL.EQ.0.) WRITE(8,316) FLOAT(ICONT3+ICONT2+ICONT1+ICONT-4)/
     1 1.17
       IF(DEL.EQ.0..AND.ICONT4.GT.0.) WRITE(8,317) FLOAT(ICONT4+ICONT3+
     1 ICONT2+ICONT1+ICONT-4)/1.17
       IF(DEL.EQ.0..AND.ICONT5.GT.0.) WRITE(8,318) FLOAT(ICONT5+ICONT4+
     1 ICONT3+ICONT2+ICONT1+ICONT-4)/1.17
312   FORMAT(2X,'DISTRIBUTION OF SUM-SQUARED ERROR FOR 117 MOLECULES',
     1/,'(OTHER MOLECULES WERE CALCULATED BUT NOT USED IN CALC:',/,2X,
     2 '<0.0625=', I3,2X,'<0.25= ',I3,2X,'<0.5625= ',I3,2X,'<1.=',I2,
     3 2X,'<1.5625=',I2,2X,'>1.5625=',I2,/ )
313   FORMAT('  % OF CALC MOLECULES WITH ENERGY DIFF OF <0.25 eV:',
     1 F7.2,' %')
314   FORMAT('  % OF CALC MOLECULES WITH ENERGY DIFF OF <0.50 eV:',
     1 F7.2,' %')
315   FORMAT('  % OF CALC MOLECULES WITH ENERGY DIFF OF <0.75 eV:',
     1 F7.2,' %')
316   FORMAT('  % OF CALC MOLECULES WITH ENERGY DIFF OF <1.00 eV:',
     1 F7.2,' %')
317   FORMAT('  % OF CALC MOLECULES WITH ENERGY DIFF OF <1.25 eV:',
     1 F7.2,' %')
318   FORMAT('  % OF CALC MOLECULES WITH ENERGY DIFF OF <1.50 eV:',
     1 F7.2,' %')
       IF(DEL.EQ.0.) WRITE(8,94)
C  KK IS THE ATOMIC NUMBER:
 62   FORMAT(4X,'BINDING ENERGY    = ',F10.4,' eV  Positive if bonding')
 94     FORMAT(1X,/)
105   FORMAT(3X,'CALCULATED ATOMIZATION ENERGY=',F9.4,3X,
     1 '0.50 x CALCULATED DIPOLE MOMENT =',F8.4)
27      FORMAT(1X)         
111      CONTINUE
59      FORMAT(1X,'DEL=',F6.4,2X,'SUMO=',E13.6,2X,'SUM= ',E13.6,2X,'ILM=
     1 ',I3,2X,'ILL= ',I3)
26      FORMAT(1X,'NPTOTAL= ',I4)
C        IF(DEL.EQ.0.) OPEN(2,FILE='FORT22')
        IF(DEL.EQ.0.) GO TO 98
        GO TO 1        
98      CONTINUE               
        CLOSE(10)
  190 FORMAT(1X,'DATE(YR,MO,DA)= ',A8,A5,' HR',A5,' MIN',2X,A5,' SEC')
      WRITE(8,189)
      CALL DATE_AND_TIME(DATE=TODAY,TIME=TIMENOW2)
      WRITE(8,190) TODAY,TIMENOW1(1:2),TIMENOW1(3:4),TIMENOW1(5:9)
      WRITE(8,188)
      WRITE(8,190) TODAY,TIMENOW2(1:2),TIMENOW2(3:4),TIMENOW2(5:9)
  188 FORMAT(/,1X,'FINISH:')    
        CLOSE(8)
300     CONTINUE
301     CONTINUE
413     FORMAT(1X,I3,2X,8F9.3)
 192    FORMAT(1X,F7.3)
 194    FORMAT(1X,F7.3)               
  20   FORMAT(20A4)                                                     
   21 FORMAT(A4,1X,A4)                                                  
   30 FORMAT(1X,20A4)
   31 FORMAT(1X,A4,1X,A4)                                               
   50 FORMAT(I4,I4,I4)                                                  
   60 FORMAT(/5X,I4,' ATOMS   CHARGE  =',I4,'   MULTIPLICITY  =',I4)    
   70 FORMAT(I4,A1,2X,F12.7,2(3X,F12.7))                                
   75 FORMAT(1X,I4,3(3X,F12.7))
   80 FORMAT(51H THIS PROGRAM IS NOT PARAMETERIZED FOR 2ND ROW INDO)      
      STOP
      END
      BLOCK DATA                                                    
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON/COEFF/YZ(311),IA(311),IB(311),LCOUNT(38),LLOW(38),
     1 LLLL(38),LPOINT(183)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      DIMENSION YZ1(90),YZ2(90),YZ3(90),YZ4(41),IA1(133),IA2(133),
     1 IA3(45),IB1(133),IB2(133),IB3(45)
      INTEGER EL,AN,U,ORB
      EQUIVALENCE (YZ4,YZ(271))
      EQUIVALENCE (YZ3,YZ(181))
      EQUIVALENCE (YZ2,YZ(91))
      EQUIVALENCE (YZ1,YZ(1))
      EQUIVALENCE (IA3,IA(267))
      EQUIVALENCE (IA2,IA(134))
      EQUIVALENCE (IA1,IA(1))
      EQUIVALENCE (IB3,IB(267))
      EQUIVALENCE (IB2,IB(134))
      EQUIVALENCE (IB1,IB(1))
      DATA EL/4H   H,4H  HE,4H  LI,4H  BE,4H   B,4H   C,4H   N,4H   O,  
     1  4H   F,4H  NE,4H  NA,4H  MG,4H  AL,4H  SI,4H   P,4H   S,4H  CL, 
     2  4H  AR/                                                         
      DATA ORB/4H   S,4H  PX,4H  PY,4H  PZ,4H DZ2,4H DXZ,4H DYZ,4HDXYZ, 
     1  4H DXY/                                                         
      DATA NORB/1,1,4,4,4,4,4,4,4,4,9,9,9,9,9,9,9,9/                    
      DATA CZ/ 1.,2.,1.,2.,3.,4.,5.,6.,7.,8.,1.,2.,3.,4.,5.,6.,7.,8./   
      DATA YZ1/  1.,  -1.,  -1.,   1., -64.,  64.,  64., -64.,  -1.,    
     1           3.,  -3.,   1.,  -1.,   2.,  -2.,   1.,  64., -64.,    
     2         -64., -64.,  64.,  64.,  64., -64., -96.,  32., 128.,    
     3          96., -96.,-128., -32.,  96.,  -1.,   5., -10.,  10.,    
     4          -5.,   1.,  -1.,   4.,  -5.,   5.,  -4.,   1.,  -1.,    
     5           1.,  -1.,   1.,   1.,  -2.,   1., -64., -64.,  64.,    
     6          64.,  64.,  64., -64., -64.,  16., -16., -16.,  16.,    
     7          16., -16., -64.,  64.,  64., -64.,  -1.,   3.,  -3.,    
     8           1.,  -1.,   1.,   2.,  -2.,  -1.,   1.,  64., -64.,    
     9        -128., 128.,  64., -64., -96.,  32., -96., 160.,  96./    
      DATA YZ2/128., -96.,  96.,-128., -96.,-160.,  96., -32.,  96.,    
     1         -16.,  16.,  16., -16.,  16., -16., -16.,  16., -16.,    
     2          16.,  16., -16.,  48., -48., -48.,  48., -48.,  48.,    
     3          48., -48.,  48., -48., -48.,  48., -64.,-128.,  64.,    
     4         128., 128.,  64.,-128., -64.,  64., -64., -64., -64.,    
     5          64.,  64.,  64., -64., -96.,  32.,  32.,  32.,  96.,    
     6          32., -32., -96., -32., -32., -32.,  96.,  -1.,   2.,    
     7           2.,  -6.,   6.,  -2.,  -2.,   1.,  -1.,   3.,  -1.,    
     8          -5.,   5.,   1.,  -3.,   1.,  64.,  64., -64.,-128.,    
     9         -64.,-128., 128.,  64., 128.,  64., -64., -64., -96./    
      DATA YZ3/ 32.,-192., 192., 288., -96., 192.,-192.,  96.,-288.,    
     1        -192., 192., -32.,  96., -16.,  16.,  32., -16., -16.,    
     2         -16., -16.,  32.,  16., -16.,  48., -48.,  48., -96.,    
     3          48., -48., -48.,  96., -48.,  48.,  48., -48.,  96.,    
     4         -48., -48.,  48., -96.,  48., -48.,  48.,  64.,-128.,    
     5         -64.,  64., 128., -64., -96.,  32., -96.,  64.,  32.,    
     6          96.,  32.,  64.,  96., -32.,  32., -96., -64., -32.,    
     7         -96., -32., -64.,  96., -32.,  96.,-144.,  96., -16.,    
     8         144., -48.,  96., -96.,  48.,-144.,  16., -96., 144.,    
     9         144.,-144.,-144., 144., 144.,-144.,-144., 144., -16./    
      DATA YZ4/ 32., -16.,  16., -48.,  32., -32.,  48., -16.,  16.,    
     1         -32.,  16.,  -1.,   5., -10.,  10.,  -5.,   1.,  -1.,    
     2           1.,   4.,  -4.,  -6.,   6.,   4.,  -4.,  -1.,   1.,    
     3           1.,  -3.,   2.,   2.,  -3.,   1.,   1.,  -1.,  -3.,    
     4           3.,   3.,  -3.,  -1.,   1./                            
      DATA LLLL/35,183,42,34,182,179,40,32,36,44,28,170,73,128,21,27,   
     1          169,166,72,69,158,127,124,19,25,157,154,64,61,119,116,  
     2          98,43,22,10,14,33,20/ 
      DATA LLOW/0,4,8,12,16,24,32,38,44,46,48,51,59,65,69,73,79,85,99,                                 
     1          111,123,131,139,151,159,167,179,193,203,223,229,249,261,
     2          269,281,287,297,303/                                    
      DATA IA1/1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 4, 5, 1, 2, 2, 
     1         3, 3, 4, 4, 5, 1, 1, 2, 3, 3, 4, 5, 5, 1, 2, 3, 4, 5, 6, 
     2         1, 2, 3, 5, 6, 7, 1, 3, 1, 2, 1, 3, 5, 1, 2, 2, 3, 3, 4, 
     3         4, 5, 1, 1, 3, 3, 5, 5, 1, 3, 3, 5, 1, 3, 5, 7, 1, 2, 3, 
     4         4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 
     5         5, 5, 6, 6, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 1, 1, 2, 
     6         2, 3, 3, 4, 4, 5, 5, 6, 6, 1, 2, 2, 3, 4, 5, 5, 6, 1, 2/ 
      DATA IA2/3, 3, 4, 4, 5, 6, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 1, 
     1         2, 3, 4, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 2, 3, 
     2         3, 4, 4, 5, 5, 6, 6, 7, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 
     3         6, 7, 7, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 1, 1, 2, 2, 2, 3, 
     4         3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 1, 3, 3, 5, 5, 
     5         7, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 
     6         7, 7, 1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 1, 1, 3, 3, 5/ 
      DATA IA3/5, 7, 7, 1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 1, 3, 5, 7, 
     1         9,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 1, 
     2         2, 3, 4, 5, 6, 7, 8/                                     
      DATA IB1/4, 3, 2, 1, 3, 4, 1, 2, 4, 3, 2, 1, 5, 4, 2, 1, 4, 3, 5, 
     1         2, 4, 1, 3, 2, 3, 5, 4, 1, 5, 2, 1, 3, 6, 5, 4, 3, 2, 1, 
     2         7, 6, 5, 3, 2, 1, 3, 1, 2, 1, 5, 3, 1, 4, 3, 5, 2, 4, 1, 
     3         3, 2, 3, 5, 1, 5, 1, 3, 3, 1, 5, 3, 7, 5, 3, 1, 6, 5, 4, 
     4         3, 2, 1, 5, 6, 3, 4, 1, 2, 4, 6, 3, 5, 2, 4, 6, 1, 3, 5,
     5         2, 4, 1, 3, 4, 6, 3, 5, 2, 6, 1, 5, 2, 4, 1, 3, 3, 5, 4, 
     6         6, 1, 5, 2, 6, 1, 3, 2, 4, 5, 4, 6, 5, 2, 1, 3, 2, 4, 3/ 
      DATA IB2/2, 6, 1, 5, 4, 3, 3, 5, 4, 6, 1, 5, 2, 6, 1, 3, 2, 4, 9, 
     1         8, 7, 6, 4, 3, 2, 1, 8, 7, 6, 5, 4, 3, 2, 1, 6, 5, 7, 4, 
     2         6, 3, 5, 2, 4, 1, 3, 2, 5, 7, 4, 6, 5, 7, 2, 6, 1, 3, 2, 
     3         4, 1, 3, 5, 7, 3, 5, 7, 1, 3, 5, 1, 3, 4, 6, 3, 5, 7, 2, 
     4         4, 6, 1, 3, 5, 7, 2, 4, 6, 1, 3, 5, 2, 4, 5, 3, 7, 1, 5, 
     5         3, 4, 6, 3, 5, 7, 2, 4, 6, 1, 3, 5, 7, 2, 4, 6, 1, 3, 5, 
     6         2, 4, 3, 5, 7, 1, 5, 7, 1, 3, 7, 1, 3, 5, 3, 5, 1, 7, 1/ 
      DATA IB3/7, 3, 5, 3, 5, 7, 1, 5, 7, 1, 3, 7, 1, 3, 5,11, 9, 7, 5, 
     1         3, 1,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 8, 
     2         7, 6, 5, 4, 3, 2, 1/                                     
      DATA LCOUNT/4,4,4,4,8,8,6,6,2,2,3,8,6,4,4,6,6,14,12,12,8,8,12,8,8,
     1          12,14,10,20,6,20,12,8,12,6,10,6,8/                      
      END
      SUBROUTINE INTGRL                                                 
C     INTEGRALS OVER ATOMIC ORBITALS FOR CNDO CALCULATIONS              
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON/COEFF/YZ(311),IA(311),IB(311),LCOUNT(38),LLOW(38)
     1 ,LLLL(38),LPOINT(183)
      COMMON/ARRAYS/S(800,800),FACT(800),B(800,800),BX(800)
      COMMON/HCORE/GPRINT(800,800)
      COMMON/LAMMA/GAMMA(200,200)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      DIMENSION BETAS(18),BETAP(18),DUMMY(800)
      DIMENSION TEMP(9,9),T(9,9),PAIRS(9,9)
      DIMENSION NC(18),LC(9),MC(9),E(3)
      REAL NUM,K1,K2
      INTEGER AN,U,ANL,ANK,OCCA,OCCB,CHARGE,EL,ORB
      CHARACTER*1 STAR
      DATA NC/1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/                      
      DATA LC/0,1,1,1,2,2,2,2,2/                                        
      DATA MC/0,1,-1,0,0,1,-1,2,-2/                                     
      DATA ZERO/0.0/,ONE/1.0/,TWO/2.0/,HALF/0.5/,IW/7/,IRR/6/
      DO 5555 I=1,38                                                    
      L=LLLL(I)                                                         
CC      WRITE(IW,996) I,L
 5555 LPOINT(L)=I                                                       
 996  FORMAT(1X,'I=  ',I2,3X,'LLL(I)=   ',I4)
C     FACT(I) CONTAINS THE FACTORIAL OF I-1                             
      FACT(1)=ONE                                                       
      FI=ZERO                                                           
      DO 15 I=2,75                                                      
      FI=FI+ONE                                                         
      FACT(I)=FACT(I-1)*FI                                              
   15 CONTINUE
C     FILL U ARRAY---U(J) IDENTIFIES THE ATOM TO WHICH ORBITAL J IS     
C     ATTACHED E.G. ORBITAL 32 ATTACHED TO  ATOM 7, ETC.                
      J=0                                                               
      DO 70 K=1,NATOMS                                                  
      LIM=NORB(AN(K))                                                   
C      WRITE(IW,995) LIM
      DO 69 I=1,LIM                                                     
      J=J+1                                                               
      U(J) = K                                                          
   69 CONTINUE
   70 CONTINUE
      DO 95 I=1,N                                                       
      DO 96 J=1,N                                                       
   96 S(I,J)=ZERO                                                   
   95 S(I,I)=ONE                                                   
C     STEP THRU PAIRS OF ATOMS                                    
      ETOTAL=ZERO                                                
      LLK=0                                                         
      DO 320 K=1,NATOMS                                              
      LLL=LLK                                                      
      ANK=AN(K)
      NORBK=NORB(ANK)
      N1=NC(ANK)
      NCANK=2*N1+1
      DO 319 L=K,NATOMS                                           
      ANL=AN(L)                                                    
      NORBL=NORB(ANL)                                             
      N2=NC(ANL)                                                  
      NCANL=2*N2+1                                                 
      IF (K.EQ.L) GO TO 265                                        
      E(1)=C(L,1)-C(K,1)                                            
      E(2)=C(L,2)-C(K,2)                                             
      E(3)=C(L,3)-C(K,3)                                              
      R=SQRT(E(1)*E(1)+E(2)*E(2)+E(3)*E(3))
C     ADD CORE REPULSION BETWEEN K AND L TO ENERGY                  
C     MODIFICATION DOES NOT USE THE VALENCE-CORE REPULSION:
C  REMOVED:
C      ETOTAL=ETOTAL+CZ(ANK)*CZ(ANL)/R
C     CALCULATE UNIT VECTOR ALONG INTERATOM AXIS,E               
      E(1)=E(1)/R                                                   
      E(2)=E(2)/R                                                   
      E(3)=E(3)/R
C     LOOP THRU PAIRS OF BASIS FUNCTIONS, ONE ON EACH ATOM         
C    2S AND 2P ORBITAL HAVE DIFFERENT VALUES OF MU:
      DO 200 I=1,NORBK                                                 
      IF(I.EQ.1.AND.ANK.LE.2) K1=EMU(ANK,K,0,0.,STAR(K))
      IF(I.EQ.1.AND.ANK.GT.2) K1=EMU(ANK,K,1,0.,STAR(K))      
      IF(I.GT.1.AND.ANK.GT.2) K1=EMU(ANK,K,2,0.,STAR(K))
C   2S AND 2P ORBITALS HAVE DIFFERENT VALUES OF MU:
      DO 200 J=1,NORBL
      IF(J.EQ.1.AND.ANL.LE.2) K2=EMU(ANL,L,0,0.,STAR(L))
      IF(J.EQ.1.AND.ANL.GT.2) K2=EMU(ANL,L,1,0.,STAR(L))
      IF(J.GT.1.AND.ANL.GT.2) K2=EMU(ANL,L,2,0.,STAR(L))
  110 IF(MC(I).NE.MC(J)) GO TO 150                                     
  120 IF(MC(I).LT.0) GO TO 140                                         
C  130 PAIRS(I,J)=SQRT((R*K1)**NCANK*(R*K2)**NCANL/(FACT(NCANK)*FACT(NCA 
C     1NL)))*(-ONE)**(LC(J)+MC(J))*SS(N1,LC(I),MC(I),N2,LC(J),K1*R,K2*R)
      RK1=R*K1
      RK2=R*K2
      FACTK=FACT(NCANK)
      FACTL=FACT(NCANL)
      LCJ=LC(J)
      LCI=LC(I)
      MCJ=MC(J)
      MCI=MC(I)
      SSS1=SS(N1,LCI,MCI,N2,LCJ,RK1,RK2)
      PAIRS(I,J)=SQRT(RK1**NCANK*RK2**NCANL/(FACTK*FACTL))*
     1 (-ONE)**(LCJ+MCJ)*SSS1
 130  CONTINUE
      GO TO 190                                                      
  140 PAIRS(I,J)=PAIRS(I-1,J-1)                                 
      GO TO 190                                                  
  150 PAIRS(I,J)=ZERO                                             
  190 CONTINUE                                                    
  200 CONTINUE                                                    
      LCULK=LC(NORBK)                                            
      LCULL=LC(NORBL)                                              
      MAXL=LCULL                                                  
      IF(LCULL.GT.LCULK) GO TO 210                                  
      MAXL=LCULK                                                    
  210 CONTINUE                                                      
C     ROTATE INTEGRALS FROM DIATOMIC BASIS TO MOLECULAR BASIS        
  220 CALL HARMTR(T,MAXL,E)                                          
      DO 230 I=1,NORBK                                             
      DO 230 J=1,NORBL                                               
      TEMP(I,J) = ZERO                                             
      DO 230 KK=1,NORBL                                              
      TEMP(I,J) = TEMP(I,J)+T(J,KK)*PAIRS(I,KK)                    
  230 CONTINUE                                                       
      DO 240 I=1,NORBK                                              
      DO 240 J=1,NORBL                                              
      PAIRS(I,J) = ZERO                                              
      DO 240 KK=1,NORBK                                           
      PAIRS(I,J) = PAIRS(I,J)+T(I,KK)*TEMP(KK,J)                     
  240 CONTINUE                                                       
C     FILL S MATRIX                                                
      DO 260 I=1,NORBK                                             
      LLKP=LLK+I                                                      
      DO 260 J=1,NORBL                                            
      LLLP=LLL+J                                                  
      S(LLKP,LLLP)=PAIRS(I,J)                                     
  260 CONTINUE                                                                                              
  265 CONTINUE                                                                                               
C     K AND L ARE THE ATOM NUMBERS OF THE PAIR OF ATOMS, K AND L
C  GAMMA SET TO ZERO HERE (AND BELOW INDO 309)  :
  310 GAMMA(K,L) = 0.
C  removed: GAMMA(K,L)=((TWO*K2)**(2*N2+1)/FACT(NCANL))*(TERM1-TERM2)   
C      WRITE(IW,997) K,L,GAMMA(K,L)
C997   FORMAT(1X,2I3,3X,'GAMMA(K,L)= ',F8.3)
      LLL=LLL+NORB(ANL)                                              
  319 CONTINUE                                                     
      LLK=LLK+NORB(ANK)                                              
  320 CONTINUE                                                        
C **********************
C     SYMMETRIZATION OF OVERLAP AND COULOMB INTEGRAL MATRICES           
C     N IS THE TOTAL NUMBER OF ORBITALS IN THE PAIR OF ATOMS, K AND L
      DO 330 I=1,N                                                    
      DO 330 J=I,N                                                 
      S(J,I) = S(I,J)                                                
  330 CONTINUE
C     IF(DEL.EQ.0.) WRITE(IW,1001)                                   
      CALL MATOUT(S,DUMMY(N),0)                                       
C     TRANSFER GAMMA TO GPRINT FOR PRINTING PURPOSES ONLY         
      DO 340 I=1,NATOMS                                             
      DO 340 J=I,NATOMS                                                
C   BOTH ONE- AND TWO- CENTER GAMMA SET TO ZERO !
      GAMMA(I,J)=0.
      GPRINT(I,J)=GAMMA(I,J)                                         
      GPRINT(J,I)=GAMMA(I,J)                                         
C      WRITE(IW,99) I,J,GAMMA(I,J)  
      GAMMA(J,I) = GAMMA(I,J)                                         
  340 CONTINUE
C99    FORMAT(1X,'I= ',I2,2X,'J= ',I2,2X,'GAMMA(I,J)= ',F10.3)
C      IF(DEL.EQ.0.) WRITE(IW,1000)                                   
      CALL MATOUT(GPRINT,DUMMY(N),NATOMS)                             
 1000 FORMAT(18H COULOMB INTEGRALS)                              
 1001 FORMAT(15H OVERLAP MATRIX/)                                     
      RETURN                                                           
      END                                                              
      FUNCTION SS(NN1,LL1,MM,NN2,LL2,ALPHA,BETA)                       
C     PROCEDURE FOR CALCULATING REDUCED OVERLAP INTEGRALS               
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON/COEFF/YZ(311),IA(311),IB(311),LCOUNT(38),LLOW(38),
     1 LLLL(38),LPOINT(183)
      COMMON/ARRAYS/S(800,800),FACT(800),A2(640800)
      DIMENSION A(17),B(17)
      REAL ZERO,HALF,ONE,FOUR,EIGHT                                   
      DATA ZERO/0.0/,HALF/0.5/,ONE/1.0/,FOUR/4.0/,EIGHT/8.0/          
      N1=NN1                                                        
      L1=LL1                                                          
      M=MM                                                         
      N2=NN2                                                         
      L2=LL2                                                         
      P =(ALPHA + BETA)*HALF                                         
      PT=(ALPHA - BETA)*HALF                                            
      X = ZERO                                                         
      M=IABS(M)                                                     
C     REVERSE QUANTUM NUMBERS IF NECESSARY                           
      IF((L2.LT.L1).OR.((L2.EQ.L1).AND.(N2.LT.N1))) GO TO 20           
   10 GO TO 30                                                        
   20 K = N1                                                          
      N1= N2                                                        
      N2= K                                                         
      K= L1                                                            
      L1= L2                                                          
      L2= K                                                           
      PT=-PT                                                         
   30 CONTINUE                                                      
C     FIND A AND B INTEGRALS                                          
      CALL ABINT(A,B,P,PT,N1+N2)                                       
      IF((L1.GT.0).OR.(L2.GT.0)) GO TO 60                              
C     BEGIN SECTION USED FOR OVERLAP INTEGRALS INVOLVING S FUNCTIONS   
C     FIND Z TABLE NUMBER L                                         
      L = (90-17*N1+N1**2-2*N2)/2                                   
      FACTOR=HALF
      GO TO 40                                                       
C     BEGIN SECTION USED FOR OVERLAPS INVOLVING NON-S FUNCTIONS         
C     FIND Y TABLE NUMBER L                                          
   60 L=(5-M)*(24-10*M+M**2)*(83-30*M+3*M**2)/120+                   
     1  (30-9*L1+L1**2-2*N1)*(28-9*L1+L1**2-2*N1)/8+                 
     2  (30-9*L2+L2**2-2*N2)/2                                         
      FACTOR=(FACT(M+2)/EIGHT)**2*SQRT( FLOAT(2*L1+1)*FACT(L1-M+1)*    
     1  FLOAT(2*L2+1)*FACT(L2-M+1)/(FOUR*FACT(L1+M+1)*FACT(L2+M+1)))   
   40 IL =LPOINT(L)                                                   
      IJ=LLOW(IL)                                                  
      NCOUNT=LCOUNT(IL)                                            
      DO 50 I=1,NCOUNT                                              
      IJ=IJ+1                                                       
   50 X=X+YZ(IJ)*A(IA(IJ))*B(IB(IJ))                                
      SS=X*FACTOR                                                   
      RETURN                                                      
      END                                                      
      SUBROUTINE HARMTR(T,MAXL,E)                           
      DIMENSION T(81),E(3)
      REAL ZERO,HALF,ONE,TWO,THREE                               
      DATA ZERO/0.0/,HALF/0.5/,ONE/1.0/,TWO/2.0/,THREE/3.0/       
      COST = E(3)                                                
      IF((ONE -COST**2).GT.0.0000000001) GO TO 20              
      SINT = ZERO                                            
      COSP = ONE                                                
      SINP = ZERO                                             
      GO TO 70                                                  
   20 SINT= SQRT(ONE -COST**2)                                  
      COSP = E(1)/SINT                                           
      SINP = E(2)/SINT                                             
   70 CONTINUE                                                      
      DO 80 I=1,81                                          
   80 T(I)=ZERO                                                
      T(1)=ONE                                                   
      IF (MAXL.GT.1) GO TO 100                                   
   90 IF (MAXL.GT.0) GO TO 110                                    
      RETURN                                                    
  100 COS2T = COST**2-SINT**2                                        
      SIN2T = TWO *SINT*COST                                     
      COS2P = COSP**2-SINP**2                                  
      SIN2P = TWO *SINP*COSP                                     
C     TRANSFORMATION MATRIX ELEMENTS FOR D FUNCTIONS          
      SQRT3= SQRT(THREE)                                            
      T(41)  =(THREE*COST**2-ONE )*HALF                          
      T(50)  = -SQRT3   *SIN2T*HALF                              
      T(68)  = SQRT3   *SINT**2*HALF                                 
      T(42)  = SQRT3   *SIN2T*COSP*HALF                           
      T(51)  = COS2T*COSP                                    
      T(60)  = -COST*SINP                                         
      T(69)  =-T(42) /SQRT3                                       
      T(78)  = SINT*SINP                                             
      T(43)  = SQRT3   *SIN2T*SINP*HALF                         
      T(52)  = COS2T*SINP                                           
      T(61)  = COST*COSP                                             
      T(70)  = -T(43) /SQRT3                                         
      T(79)  = -SINT*COSP                                            
      T(44)  = SQRT3   *SINT**2*COS2P*HALF                           
      T(53)  = SIN2T*COS2P*HALF                                      
      T(62)  = -SINT*SIN2P                                           
      T(71)   = (ONE +COST**2)*COS2P*HALF                            
      T(80)  = -COST*SIN2P                                         
      T(45)  = SQRT3   *SINT**2*SIN2P*HALF                         
      T(54)  = SIN2T*SIN2P*HALF                                   
      T(63)  = SINT*COS2P                                            
      T(72)  = (ONE +COST**2)*SIN2P*HALF                               
      T(81)  = COST*COS2P                                              
C     TRANSFORMATION MATRIX ELEMENTS FOR P FUNCTIONS                   
  110 T(11)  = COST*COSP                                       
      T(20)  = -SINP                                            
      T(29)  = SINT*COSP                                      
      T(12)  = COST*SINP                                          
      T(21)  = COSP                                            
      T(30)  = SINT*SINP                                     
      T(13)  = -SINT                                          
      T(31)  = COST                                             
      RETURN                                                   
      END                                                      
      SUBROUTINE ABINT(A,B,P,PT,K)                                 
C     FILLS ARRAY OF B-INTEGRALS. NOTE THAT B(I) IS B(I-1) IN THE       
C     USUAL NOTATION                                                   
C     FOR X.GT.3                    EXPONENTIAL FORMULA IS USED        
C     FOR 2.LT.X.LE.3 AND K.LE.10   EXPONENTIAL FORMULA IS USED        
C     FOR 2.LT.X.LE.3 AND K.GT.10   15 TERM SERIES IS USED          
C     FOR 1.LT.X .E.2 AND K.LE.7    EXPONENTIAL FORMULA IS USED     
C     FOR 1.LT.X.LE.2 AND K.GT.7    12 TERM SERIES IS USED         
C     FOR .5.LT.X.LE.1 AND K.LE.5   EXPONENTIAL FORMULA IS USED     
C     FOR .5.LT.X.LE.1 AND K.GT.5    7 TERM SERIES IS USED          
C     FOR X.LE..5                    6 TERM SERIES IS USED         
C     ************************************************************
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(17),B(17)
      DATA ZERO/0.0/,HALF/0.5/,ONE/1.0/,TWO/2.0/,THREE/3.0/        
      DO 5 I=1,17                                                    
      A(I)=ZERO                                                      
      B(I)=ZERO                                                     
    5 CONTINUE                                                        
      X=PT                                                           
      KP=K+1                                                         
      ABSX= ABS(X)                                                   
      EXPX=EXP(X)                                                    
      EXPMX= ONE/EXPX                                                
      IF(ABSX.GT.THREE) GO TO 120                                    
   10 IF(ABSX.GT.TWO ) GO TO 20                                      
   40 IF(ABSX.GT.ONE ) GO TO 50                                     
   70 IF(ABSX.GT.HALF) GO TO 80                                      
  100 IF(ABSX.GT..000001D0) GO TO 110                                
      GO TO 170                                                      
  110 LAST =7                                                    
      GO TO 140                                                      
   80 IF(K.LE.5) GO TO 120                                         
   90 LAST=8                                                       
      GO TO 140                                                    
   50 IF(K.LE.7) GO TO 120                                          
   60 LAST=13                                                        
      GO TO 140                                                      
   20 IF(K.LE.10) GO TO 120                                         
   30 LAST=16                                                      
      GO TO 140                                                    
  120 CONTINUE                                                     
      B(1)=(EXPX-EXPMX)/X                                         
      FI=ZERO                                                 
      DO 130 I=1,K                                              
      FI=FI+ONE                                                  
  130 B(I+1)=(FI*B(I)+(-ONE)**I*EXPX-EXPMX)/X                       
      GO TO 190                                                      
  140 FI=ZERO                                               
      DO 160 I=1,KP                                        
      Y=ZERO                                              
      FI=FI+ONE                                                
      FACTM=ONE                                               
      FM=ZERO                                                     
      DO 150 M=1,LAST                                           
      FACTMM=FACTM                                                
      FM=FM+ONE                                                    
      FACTM=FACTM*FM                                                 
  150 Y=Y+(-X)**(M-1)*(ONE-(-ONE)**(M+I-1))/(FACTMM*(FM+FI-ONE))     
  160 B(I)=Y                                                         
      GO TO 190                                                    
  170 FI=ZERO                                                  
      DO 180 I=1,KP                                                 
      FI=FI+ONE                                                       
  180 B(I)=(ONE-(-ONE)**I)/FI                                   
  190 CONTINUE                                                  
      X=P                                                      
      EXPMX=ONE/EXP(X)                                          
      A(1)=EXPMX/X                                              
      FI=ZERO                                                    
      DO 200 I=1,K                                               
      FI=FI+ONE                                                
  200 A(I+1)=(A(I)*FI+EXPMX)/X                                     
      RETURN                                                     
      END                                                     
      SUBROUTINE GUESSP(IOP1)                             
C     CALCULATES AN EXTENDED HUCKEL TYPE OF HAMILTONIAN FOR THE        
C     INITIAL GUESS AND THEN ADDS TERMS TO FORM THE CORE HAMILTONIAN   
C     FOR CNDO OR INDO.  AT THE END OF THIS ROUTINE, THE FIRST GUESS    
C     AT THE EIGENVECTORS IS IN A AND THE CORE HAMILTONIAN IS IN H.     
C     OVERLAPS ARE IN MATRIX A, COULOMB INTEGRALS (GAMMA) ARE IN MATRIX 
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON/ARRAYS/A(800,800),EIG(800),B(800,800),BX(800)
      COMMON/HCORE/H(800,800)
      COMMON/LAMMA/G(200,200)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SMK/C1,C2,C3,C4,C5,C6,C7,C8,C9,CA,CB,CC,CD,CE,CF,CG,CH,CI,
     1 CJ,CK,CL,CM,CN,CO,CP1,CP2,CP3,CQ1,CQ2,CQ3,CR1,CR2,CS,CHHH
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      COMMON/BOB/UZ
      CHARACTER*4 IOP1,IOP2
      CHARACTER*1 STAR,STRK,STRL
      DIMENSION BASE(200),AK(800),EGS(800),DD1(200,200),UXX(200,200)
      DIMENSION BETAS(18),BETAP(18),BB(200,2),UU(200,200),Q1(200,2)
      DIMENSION ENEGS(18),ENEGP(18),ENEGD(18),BETA0(18)
      INTEGER OCCA,OCCB,UL,AN,U,ANI,CHARGE,ANK,ANL,EL,ORB
      REAL ZERO,HALF,ONE,TWO,THREE,FOUR,FIVE,SIX,TWOFIV,TWELV,C2P5,P75,
     1 C1P5,IP,NETAX                                                                           
      DATA ZERO/0./,HALF/.5/,ONE/1./,TWO/2./,THREE/3./,FOUR/4./,SIX/6./ 
      DATA TWOFIV/25./,TWELV/12./,C2P5/2.5/,P75/.75/,C1P5/1.5/          
      DATA EVAU/27.2114/,RHO/1.0E-6/,IW/7/,A0/0.529167/
C  DD IS IN ATOMIC UNITS, CTAX AND UZ ARE IN UNITS OF eV (UNTIL INDO 560):
      J=0
      DO 50 I=1,NATOMS                                                  
      ANI=AN(I)                                                         
      NORBI=NORB(ANI)
C     Q1(I,L) IS INDEXED BY ATOM ORDER, I, IN THE MOLECULE AND
C     BY L=1 FOR S ORBITALS OR L=2 FOR P ORBITALS:
C  HETERONUCLEAR DIATOMIC NEUTRAL AND IONIC MOLECULES:  
C  THE CHARGES ARE ASSIGNED TO THE OUTER, OCCUPIED, VALENCE ELECTRONS:
      IF(NATOMS.EQ.2.AND.AN(1).EQ.AN(2)) CHRG(I)=0.
      IF(ANI.EQ.1) Q1(I,1)=CHRG(I)
      IF(ANI.EQ.1) Q1(I,2)=0.
      IF(ANI.GT.4) Q1(I,1)=0.
      IF(ANI.GT.4) Q1(I,2)=CHRG(I)
      CCQ= 0.
      IF(ANI.EQ.1)                    BASE(I)= CHHH +CCQ
      IF(ANI.EQ.6.AND.STAR(I).EQ.'A') BASE(I)= CD   +CCQ
      IF(ANI.EQ.6.AND.STAR(I).EQ.'B') BASE(I)= CE   +CCQ
      IF(ANI.EQ.6.AND.STAR(I).EQ.'C') BASE(I)= CF   +CCQ
      IF(ANI.EQ.7.AND.STAR(I).EQ.'A') BASE(I)= CH   +CCQ
      IF(ANI.EQ.7.AND.STAR(I).EQ.'B') BASE(I)= CI   +CCQ 
      IF(ANI.EQ.7.AND.STAR(I).EQ.'C') BASE(I)= CJ   +CCQ
      IF(ANI.EQ.8.AND.STAR(I).EQ.'A') BASE(I)= CL   +CCQ
      IF(ANI.EQ.8.AND.STAR(I).EQ.'B') BASE(I)= CM   +CCQ
      IF(ANI.EQ.9)                    BASE(I)= CN   +CCQ
C   
      DO 46 L=1,NORBI                                                   
      J=J+1                                                             
C  J IS THE ORBITAL #:
      IF(ANI.EQ.1)                               AK(J)= C1
C  THE MULTIPLIERS FOR 2P ORBITALS (L>1) ARE ASSIGNED TO J ORBITAL #S:
      IF(L.GT.1.AND.ANI.EQ.6.AND.STAR(I).EQ.'A') AK(J)= C2
      IF(L.GT.1.AND.ANI.EQ.6.AND.STAR(I).EQ.'B') AK(J)= C3
      IF(L.GT.1.AND.ANI.EQ.6.AND.STAR(I).EQ.'C') AK(J)= C4
      IF(L.GT.1.AND.ANI.EQ.7.AND.STAR(I).EQ.'A') AK(J)= C6 
      IF(L.GT.1.AND.ANI.EQ.7.AND.STAR(I).EQ.'B') AK(J)= C7 
      IF(L.GT.1.AND.ANI.EQ.7.AND.STAR(I).EQ.'C') AK(J)= C8
      IF(L.GT.1.AND.ANI.EQ.8.AND.STAR(I).EQ.'A') AK(J)= CA
      IF(L.GT.1.AND.ANI.EQ.8.AND.STAR(I).EQ.'B') AK(J)= CB
      IF(L.GT.1.AND.ANI.EQ.9)                    AK(J)= CC
C
C  THE MULTIPLIERS FOR 2S ORBITALS (L=1) ARE ASSIGNED TO THE J ORBITAL #S:
      IF(L.EQ.1.AND.ANI.EQ.6.AND.STAR(I).EQ.'A') AK(J)= CP1 
      IF(L.EQ.1.AND.ANI.EQ.6.AND.STAR(I).EQ.'B') AK(J)= CP2
      IF(L.EQ.1.AND.ANI.EQ.6.AND.STAR(I).EQ.'C') AK(J)= CP3
      IF(L.EQ.1.AND.ANI.EQ.7.AND.STAR(I).EQ.'A') AK(J)= CQ1
      IF(L.EQ.1.AND.ANI.EQ.7.AND.STAR(I).EQ.'B') AK(J)= CQ2  
      IF(L.EQ.1.AND.ANI.EQ.7.AND.STAR(I).EQ.'C') AK(J)= CQ3 
      IF(L.EQ.1.AND.ANI.EQ.8.AND.STAR(I).EQ.'A') AK(J)= CR1      
      IF(L.EQ.1.AND.ANI.EQ.8.AND.STAR(I).EQ.'B') AK(J)= CR2
      IF(L.EQ.1.AND.ANI.EQ.9)                    AK(J)= CS
C
C  1S:     
      IF(L.EQ.1.AND.ANI.LE.2)        EMUS= EMU(ANI,I,0,0.,STAR(I))
C  2S:      
      IF(L.EQ.1.AND.ANI.GT.4)        EMUS= EMU(ANI,I,1,0.,STAR(I))
C  1S OR 2S ENERGY:
      IF(ANI.EQ.1)                    ENEGS(ANI)= 13.59843      
      IF(L.EQ.1.AND.ANI.GT.4)         ENEGS(ANI)= 13.6057*EMUS*EMUS
      IF(L.EQ.1)                      H(J,J)=ENEGS(ANI)
      IF(L.EQ.1.AND.ANI.EQ.1)         EGS(J)= 1.
      IF(L.EQ.1.AND.ANI.GT.4)         EGS(J)= 1.
C  2P: 
      IF(L.GT.1.AND.ANI.LE.4)         EMUP=1.
      IF(L.GT.1.AND.ANI.GT.4)         EMUP= EMU(ANI,I,2,0.,STAR(I))
      IF(L.GT.1.AND.ANI.GT.4)         ENEGP(ANI)= 13.6057*EMUP*EMUP   
      IF(L.GT.1.AND.ANI.GT.4)         H(J,J)=ENEGP(ANI)
      IF(L.GT.1.AND.ANI.GT.4)         EGS(J)= 0.
   46 CONTINUE
   50 CONTINUE
C  DD IS IN ATOMIC UNITS, CTAX AND UZ ARE IN UNITS OF eV (UNTIL INDO 560):
      UZ=0.  
      DO 177 IK=1,NATOMS
      ANK=AN(IK)
      DO 176 IL=IK,NATOMS
      ANL=AN(IL)
      IF(IK.EQ.IL) GO TO 176
      ED1=C(IK,1)-C(IL,1)
      ED2=C(IK,2)-C(IL,2)
      ED3=C(IK,3)-C(IL,3)      
      DD1(IK,IL)=SQRT(ED1*ED1+ED2*ED2+ED3*ED3)
      DD1(IL,IK)=DD1(IK,IL) 
C  CHANGED TO UX6:
C  UX6 REPULSIONS FOR INTERATOMIC DISTANCES GREATER THAN DD1= 5. BOHR 
C  OR ~2.7 ANGSTROMS RETURNED BY FUNCTION UX6 AS ZERO:  
      DD=DD1(IK,IL)
C      IF(DD.LT.1.5) WRITE(8,889) IK,IL,DD*A0,C(IK,1)*A0,C(IK,2)*A0,
C     1 C(IK,3)*A0,C(IL,1)*A0,C(IL,2)*A0,C(IL,3)*A0
C  889 FORMAT(1X,2I4,3X,7F7.3)    
      STRK=STAR(IK)
      STRL=STAR(IL)
      IF(DD.LE.5.AND.ANK.LE.ANL) UZ=UZ+ UX6(ANK,STRK,ANL,STRL,DD)
      IF(DD.LE.5.AND.ANK.GT.ANL) UZ=UZ+ UX6(ANL,STRL,ANK,STRK,DD)
C      REPULSION BY LIKE CHARGES:
      IF(CHRG(IK).GE.0..AND.CHRG(IL).GE.0.) ADJ=1.     
      IF(CHRG(IK).LT.0..AND.CHRG(IL).LT.0.) ADJ=1.
C      ATTRACTION BY OPPOSITE CHARGES:
      IF(CHRG(IK)*CHRG(IL).LT.0.)           ADJ=1./0.50      
      CHIK=CHRG(IK)
      IF(CHIK.LT.-1.) CHIK= -1.
      IF(CHIK.GT.+1.) CHIK= +1.
      CHIL=CHRG(IL)
      IF(CHIL.LT.-1.) CHIL= -1.
      IF(CHIL.GT.+1.) CHIL= +1.
      IF(DD1(IK,IL).GT.5.) UU(IK,IL)=0.     
      IF(DD1(IK,IL).LE.5.) UU(IK,IL)=ADJ*CHIK*CHIL/DD1(IK,IL)
      UU(IL,IK)=UU(IK,IL) 
C  UX6 REPULSION HAS BEEN CONVERTED TO EXPONENTIAL FORM IN UNITS OF eV:    
176   CONTINUE
177   CONTINUE
C     FORM HUCKEL HAMILTONIAN IN A (OFF DIAGONAL TWO CENTER TERMS)      
C      WRITE(IW,998) 1,1,A(1,1),H(1,1),H(1,1)
C      WRITE(IW,998) 2,2,A(2,2),H(2,2),H(2,2) 
      DO 90 I=2,N                                                
      K=U(I)                                                    
      L=AN(K)                                                
      UL=I-1                                                    
      DO 90 J=1,UL                                              
      KK=U(J)                                                    
      LL=AN(KK)                                                
CC      WRITE(IW,986) N,I,J,A(I,J),H(I,I),H(J,J),ZZK
 986  FORMAT(1X,3I3,4F9.4)
C  OFF-DIAGONAL TERMS ARE ASSIGNED HERE (I CANNOT EQUAL J):
      IF(K.EQ.KK) GO TO 90
C
      ABC=(0.6 +2.6*UU(K,KK))*EGS(I)*EGS(J)
      H(I,J)=A(I,J)*((AK(I)*H(I,I)+AK(J)*H(J,J))/2. -ABC)
C         
C  K AND KK ARE THE ATOM#s THAT GO WITH THE MOLECULAR ORBITALS, I AND J           
   90 CONTINUE
      BBB=0.0
      J=0
      DO 55 I=1,NATOMS
      BB(I,1)=0.
      BB(I,2)=0.
      ANI=AN(I)
      NORBI=NORB(ANI)                                           
      DO 47 L=1,NORBI                                                  
      J=J+1 
      IF(ANI.LE.4) BB(I,1)=BASE(I)
      IF(ANI.GT.4) BB(I,1)=ENEGS(ANI) -ENEGP(ANI) +BASE(I)    
      IF(L.EQ.1)   H(J,J)= BB(I,1)
C  ENERGIES OF 2P ORBITALS SET TO BASE(ANI):[ENEGP(ANI)-ENEGP(ANI)]=0 +BASE(ANI):      
      IF(L.EQ.2)   BB(I,2)= BASE(I)
      IF(L.GE.2)   H(J,J)=  BB(I,2)
C  THE FOLLOWING EQUATIONS SUBTRACT AWAY THE BASE VALUES AND COMPENSATE
C  FOR THE ELECTRONS (1S, 2S, and 2P) OF EACH ATOM IN THE MOLECULE:
C  EXPECTED CORRECTIONS:
      IF(L.EQ.1.AND.ANI.EQ.1) BBB=BBB+BB(I,1)*1.
      IF(L.EQ.1.AND.ANI.GT.4) BBB=BBB+BB(I,1)*2.
      IF(L.EQ.2.AND.ANI.GT.4) BBB=BBB+BB(I,2)*(CZ(ANI)-2.)
C  *****************************
 47   CONTINUE 
 55   CONTINUE
C      IF(DEL.EQ.0.) WRITE(8,888) BBB
C 888  FORMAT(1X,'BBB=', F14.4)  
115    FORMAT(1X,I2,' GROUND=',F11.4,2X,'BOUND=',F11.4,2X,'EBASE=',
     1  F11.4,2X,'Q(I,1/2)=',F6.3)
C      IF(DEL.EQ.0.) WRITE(8,222) DD1*A0,UZ
C 222  FORMAT(1X,'DD= ',F7.4,3X,'UZ= ',F10.4,2X,'eV')      
      BBB=BBB +UZ
C     TRANSFER A COPY OF H INTO A FOR USE AS INITIAL GUESS OF F MATRIX  
      DO 105 I=1,N                                                      
      DO 105 J=I,N                                                  
C  ELECTRON VOLTS ARE CONVERTED TO -AU BY DIVIDING BY EVAU =27.2114
      TEMP= -H(J,I)/EVAU                                             
      H(J,I)=TEMP                                                 
      H(I,J)=TEMP                                                  
      A(J,I)=TEMP                                                  
  105 CONTINUE                                                                         
  300 CONTINUE
C      IF(DEL.EQ.0.) WRITE(IW,1000)                                   
 1000 FORMAT(17H CORE HAMILTONIAN)                                   
C  99  FORMAT(1X,'J=',2X,I2,3X,'H(J,J)=',2X,F10.2)
      CALL MATOUT(H,EIG,0)                                          
      CALL EIGN1M(N,RHO,A,EIG)                                      
      RETURN                                                        
      END                                                           
      SUBROUTINE SCFC(IOP1)                                          
C     CNDO/INDO CLOSED SHELL SCF SEGMENT                             
C     CORE HAMILTONIAN IS CONTAINED IN H, EIGENVECTORS IN A         
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)     
      COMMON/ARRAYS/A(800,800),EIG(800),B(800,800),BX(800)
      COMMON/HCORE/H(800,800)
      COMMON/LAMMA/G(200,200)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      INTEGER ORB,U,AN,OCCA,OCCB,ONCE,CHARGE,EL
      CHARACTER*4 IOP1,IOP2                                             
      CHARACTER*1 STAR
C  DATA F2 AND DATA G1:
      DATA IW/7/,NIT/25/,CNVERG/1.0E-6/,RHO/1.0E-6/                     
      DATA ZERO/0.0/,HALF/0.5/,TWO/2.0/,SIX/6.0/,P14/.14/,P22/.22/      
      IT=0                                                              
      ENERGY=ZERO                                                       
   10 IT=IT+1                                                           
      ONCE=ZERO
      OLDENG=ENERGY                                                     
C     FORM THE DENSITY MATRIX IN A FROM THE EIGENVECTORS IN A.         
CC      IF(DEL.EQ.0.) WRITE(IW,918) OCCA
 918  FORMAT(1X,'OCCA (SCFC)= ',I4)
      DO 280 I=1,N                                                  
      II=U(I)
      JJ=0
      IF((I+1).LE.N) JJ=U(I+1)
      DO 260 J=I,N                                                       
      EIG(J)=ZERO                                                      
      DO 250 K=1,OCCA                                                 
  250 EIG(J)=EIG(J)+A(I,K)*A(J,K)                                    
  260 CONTINUE                                                       
      DO 270 J=I,N                                                    
      EIG(J)=TWO*EIG(J)                                              
  270 A(I,J)=EIG(J)                                                   
  280 CONTINUE                                                         
      IF (IT.GT.NIT) GO TO 300                                         
C     FORM THE CNDO PART OF THE F MATRIX                                
C     AT THE BEGINNING OF THIS SECTION, THE ELEMENTS A(I,I) AND EIG(I)  
C     CONTAIN THE DIAGONAL ELEMENTS OF THE DENSITY MATRIX.  AFTER THE   
C     DO LOOP, A(I,I) CONTAINS THE DIAGONAL ELEMENTS OF THE FOCK MATRIX.
      DO 20 I=1,N                                                      
      II=U(I)                                                         
      DO 40 J=1,I                                                    
      JJ=U(J)                                                       
      A(I,J)=H(I,J)-A(J,I)*HALF*G(II,JJ)                              
   40 CONTINUE                                                         
      DO 30 J=1,N                                                     
      JJ=U(J)                                                         
      A(I,I)=A(I,I)+EIG(J)*G(II,JJ)                                   
   30 CONTINUE                                                        
   20 CONTINUE                                                        
      IF(IOP1.EQ.'CNDO') GO TO 101                                      
      KK=U(J)                                                        
      I=1                                                           
      DO 80 II=1,NATOMS                                              
      K=AN(II)                                                        
      IF (K.LE.2) GO TO 90                                                          
   70 CONTINUE                                                        
   90 CONTINUE                                                         
      I=I+NORB(K)                                                       
   80 CONTINUE                                                        
C     COMPUTE THE ELECTRONIC ENERGY                                   
  101 ENERGY=ZERO                                                    
      DO 100 I=1,N                                                   
  100 ENERGY=ENERGY+HALF*EIG(I)*(A(I,I)+H(I,I))                        
      NM=N-1                                                           
      DO 105 I=1,NM                                                    
      IP=I+1                                                            
      DO 105 J=IP,N                                                     
      ENERGY=ENERGY+A(I,J)*(A(J,I)+H(J,I))                             
  105 CONTINUE
CC      WRITE(IW,130)  IT,ENERGY                                       
  130 FORMAT(21H ENERGY ON ITERATION ,I2,4H IS ,E20.10)                
      IF(ABS(ENERGY-OLDENG).GE.CNVERG) GO TO 150                       
C     IF CONVERGENCE, PUT THE ITERATION COUNTER AT ITS LIMIT          
      IT=NIT                                                           
C     COMPLETE HARTREE-FOCK MATRIX MAY BE PRINTED HERE                
      DO 170 I=1,N                                                     
      DO 170 J=1,I                                                     
C      IF(DEL.EQ.ZERO.AND.I.EQ.J) WRITE(IW,99) A(I,I)
99    FORMAT(1X,'J=',2X,I2,3X,'H(J,J)=',2X,F10.2)
      A(J,I)=A(I,J)                                                    
  170 CONTINUE
CC      IF(DEL.EQ.0.) WRITE(IW,1000)                                   
 1000 FORMAT(12H1FOCK MATRIX)                                           
      CALL MATOUT(A,EIG,0)                                            
  150 CALL EIGN1M(N,RHO,A,EIG)                                      
C     IF NO CONVERGENCE, GO BACK UP TOP AGAIN                       
      IF(IT.LT.NIT) GO TO 10                                        
C     EIGENVALUES AND EIGENVECTORS MAY BE PRINTED HERE            
CC      IF(DEL.EQ.0.) WRITE(IW,1010)                               
 1010 FORMAT(29H EIGENVALUES AND EIGENVECTORS)                      
      CALL MATOUT(A,EIG,1)                                          
  190 CONTINUE                                                      
C     IF CONVERGENCE, GO BACK TOP ONCE MORE FOR DENSITY MATRIX         
      GO TO 10                                                        
  300 CONTINUE                                                        
      DO 180 I=1,N                                                    
      DO 180 J=I,N                                                   
      A(J,I)=A(I,J)                                                   
  180 CONTINUE
CC      IF(DEL.EQ.0.) WRITE(IW,1030)                                 
 1030 FORMAT(15H DENSITY MATRIX)                                      
      CALL MATOUT(A,EIG,0)                                             
      RETURN                                                           
      END                                                              
      SUBROUTINE SCFO(IOP1)                                             
C     OPEN-SHELL UNRESTRICTED CALCULATION USING POPLE-NESBET EQUATIONS  
C     CORE HAMILTONIAN IS CONTAINED IN H, EIGENVECTORS IN A             
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)     
      COMMON/ARRAYS/A(800,800),EIGA(800),B(800,800),EIGB(800)
      COMMON/HCORE/H(800,800)
      COMMON/LAMMA/G(200,200)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      INTEGER OCCA,OCCB,U,AN,CHARGE,EL,ORB,ANK,ANL
      REAL ZERO,HALF,TWO,SIX,THREE,P2,P24,TWOFIV                               
      CHARACTER*4 IOP1,IOP2                                             
      CHARACTER*1 STAR
      DATA IW/7/,NIT/5/,CNVERG/1.0E-6/,RHO/1.0E-6/                      
      DATA ZERO/0./,HALF/.5/,TWO/2./,SIX/6./,THREE/3./,P2/.2/,P24/.24/  
      DATA TWOFIV/25./                                                 
      IT=0                                                             
      DO 5 I=1,N                                                        
      DO 5 J=1,N                                                       
      B(I,J)=A(I,J)                                                    
    5 CONTINUE
      ENERGY=ZERO                                                      
   10 IT=IT+1                                                          
      OLDENG=ENERGY                                                     
C     FORM DENSITY MATRIX IN A FROM EIGENVECTORS IN A AND DENSITY MATRIX
C     IN B FROM EIGENVECTORS IN B                                       
C     AT THE END OF THIS SECTION THE ARRAYS EIGA AND EIGB CONTAIN       
C     THE DIAGONAL VALUES OF P(ALPHA) AND P(BETA)                      
      DO 280 I=1,N                                                   
      DO 260 J=I,N                                                     
      EIGA(J)=ZERO                                                     
      DO 250 K=1,OCCA                                                 
  250 EIGA(J)=EIGA(J)+A(I,K)*A(J,K)                                    
      EIGB(J)=ZERO                                                      
      DO 251 K=1,OCCB                                                  
  251 EIGB(J)=EIGB(J)+B(I,K)*B(J,K)                                    
  260 CONTINUE                                                        
      DO 270 J=I,N                                                      
      A(I,J)=EIGA(J)                                                    
      B(I,J)=EIGB(J)                                                    
  270 CONTINUE                                                         
  280 CONTINUE                                                         
      IF (IT.GT.NIT) GO TO 300                                          
C     FORM F(ALPHA) IN THE LOWER TRIANGLE OF A AND F(BETA) IN THE LOWER 
C     TRIANGLE OF B FROM P(ALPHA) IN THE UPPER TRIANGLE OF A (DIAGONAL  
C     ELEMENTS IN EIGA) AND P(BETA) IN THE UPPER TRIANGLE OF B (DIA-    
C     GONAL ELEMENTS IN EIGB)                                           
C     FORM THE CNDO PART OF THE F MATRIX                                
      DO 20 I=1,N                                                       
      II=U(I)                                                           
      DO 40 J=1,I                                                       
      JJ=U(J)                                                        
      A(I,J)=H(I,J)-A(J,I)*G(II,JJ)                                  
      B(I,J)=H(I,J)-B(J,I)*G(II,JJ)                                  
   40 CONTINUE                                                      
      DO 30 J=1,N                                                    
      JJ=U(J)                                                      
      TEMP=(EIGA(J)+EIGB(J))*G(II,JJ)                                
      A(I,I)=A(I,I)+TEMP                                             
      B(I,I)=B(I,I)+TEMP                                             
   30 CONTINUE                                                      
   20 CONTINUE                                                        
      IF(IOP1.EQ.'CNDO') GO TO 101                                   
C     INDO MODIFICATION TO COULOMB AND EXCHANGE INTEGRALS            
      I=1                                                           
      DO 99 II=1,NATOMS                                             
      K=AN(II)                                                       
      IF (K.LE.2) GO TO 90                                            
   90 CONTINUE                                                       
      I=I+NORB(K)                                                    
   99 CONTINUE                                                      
C     COMPUTE THE ELECTRONIC ENERGY                                
  101 ENERGY=ZERO                                                
      DO 100 I=1,N                                                   
C      EGO(I)= EIGA(I)*(A(I,I)+H(I,I))+EIGB(I)*(B(I,I)+H(I,I))
C      WRITE(7,350) I,EGO(I)
      ENERGY=ENERGY+EIGA(I)*(A(I,I)+H(I,I))+EIGB(I)*(B(I,I)+H(I,I))     
  100 CONTINUE                                                     
C 350  FORMAT(10X,I2,F20.4)
      ENERGY=HALF*ENERGY                                           
      NM=N-1                                                          
      DO 105 I=1,NM                                                   
      IP=I+1                                                          
      DO 105 J=IP,N                                                   
      ENERGY=ENERGY+A(I,J)*(A(J,I)+H(J,I))+B(I,J)*(B(J,I)+H(J,I))      
  105 CONTINUE                                                         
CC      IF(DEL.EQ.0.) WRITE(IW,130)  IT,ENERGY                          
  130 FORMAT(21H ENERGY ON ITERATION ,I2,4H IS ,E20.10)                 
      IF(ABS(ENERGY-OLDENG).GE.CNVERG) GO TO 150                        
C     IF CONVERGENCE, PUT THE ITERATION COUNTER AT ITS LIMIT           
      IT=NIT                                                            
C     COMPLETE HARTREE-FOCK MATRIX MAY BE PRINTED HERE                
      DO 170 I=1,N                                                      
      DO 170 J=1,I                                                      
      B(J,I)=B(I,J)                                                    
      A(J,I)=A(I,J)                                                  
  170 CONTINUE
CC      WRITE(IW,1000)                                                  
 1000 FORMAT(18H1ALPHA FOCK MATRIX)                                    
      CALL MATOUT(A,EIGA,0)                                          
CC      WRITE(IW,1005)                                                 
 1005 FORMAT(18H  BETA FOCK MATRIX)                                     
      CALL MATOUT(B,EIGB,0)                                            
  150 CALL EIGN1M(N,RHO,A,EIGA)                                       
      CALL EIGN1M(N,RHO,B,EIGB)                                        
C     THE EXPECTATION VALUE OF S**2 IS COMPUTED HERE BUT ONLY PRINTED  
C     AFTER CONVERGENCE.                                               
      MU=FLOAT(OCCB)                                                   
      S=0.5*FLOAT(OCCA-OCCB)                                         
      S=S*(S+1)+MU                                                  
      DO 181 I=1,OCCA                                                  
      DO 181 J=1,OCCB                                                
      TEMP=0.0                                                         
      DO 182 K=1,N                                                   
      TEMP=TEMP+A(K,I)*B(K,J)                                          
  182 CONTINUE                                                       
      S=S-TEMP*TEMP                                                     
  181 CONTINUE                                                         
C     IF NO CONVERGENCE, GO BACK UP TOP AGAIN                          
      IF(IT.LT.NIT) GO TO 10                                         
CC      WRITE(IW,1040) S                                                
 1040 FORMAT(//42H THE EXPECTATION VALUE OF S(SQUARED) IS   ,F10.5//)   
C     EIGENVALUES AND EIGENVECTORS MAY BE PRINTED HERE                  
CC      IF(DEL.EQ.0.) WRITE(IW,1010)                                   
 1010 FORMAT(35H ALPHA EIGENVALUES AND EIGENVECTORS)                   
      CALL MATOUT(A,EIGA,1)                                            
CC      WRITE(IW,1015)                                                 
 1015 FORMAT(35H  BETA EIGENVALUES AND EIGENVECTORS)                   
      CALL MATOUT(B,EIGB,1)                                           
C     IF CONVERGENCE, GO BACK TOP ONCE MORE FOR DENSITY MATRIX         
      GO TO 10                                                       
  300 CONTINUE                                                        
      DO 180 I=1,N                                                  
      DO 180 J=I,N                                                    
      A(J,I)=A(I,J)                                                  
      B(J,I)=B(I,J)                                                  
  180 CONTINUE
CC      IF(DEL.EQ.0.) WRITE(IW,1030)                                  
 1030 FORMAT(21H ALPHA DENSITY MATRIX)                              
      CALL MATOUT(A,EIGA,0)                                         
CC      IF(DEL.EQ.0.) WRITE(IW,1035)                               
 1035 FORMAT(21H  BETA DENSITY MATRIX)                              
      CALL MATOUT(B,EIGB,0)                                          
      RETURN                                                         
      END            
      SUBROUTINE CONVRG(IT,ACURCY,ICOUNT,N,NTT2,NO)                    
C     THIS ROUTINE IS A MODIFIED VERSION OF THE CONVERGENCE           
C     ROUTINE APPEARING IN GAUSSIAN 70.                               
C      N IS THE SIZE OF THE BASIS SET,NTT2 IS THE NUMBER OF DENSITY   
C     MATRIX ELEMENTS AND NO IS 1 FOR CLOSED SHELLS, 2 FOR OPEN SHELLS
C     IT IS THE ITERATION NUMBER, ICOUNT IS THE NUMBER OF ITERATIONS  
C     SINCE THE LAST EXTRAPOLATION.                                   
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON/ARRAYS/A3(640800),A2(640800)
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      INTEGER CHARGE
      REAL PT99,PT995,ONE,ONEPT9,TWO,FOUR
      DATA IOUT/7/                                                    
      DATA PT99/.99/,PT995/.995/,ONE/1./,ONEPT9/1.9/,TWO/2./,FOUR/4./ 
   10 FORMAT(17H CONVERGENCE DADA,20X,E15.6,2X,E15.6,2X,E15.6)  
   20 FORMAT(22H 3-POINT EXTRAPOLATION)                              
   30 FORMAT(22H 4-POINT EXTRAPOLATION)                               
      ICOUNT=ICOUNT+1                                                 
      LOC1=33                                                   
      IF(ICOUNT-(ICOUNT/2)*2)40,50,40                                  
   40 LOC1=32                                                          
   50 CONTINUE                                                       
      LOC2=65-LOC1                                                    
C     READ IN CURRENT DENSITY MATRIX (SQUEEZED)                                          
C     SKIP TO  END IN FIRST CYCLE OR IMMEDIATELY AFTER EXTRAPOLATION  
      IF(ICOUNT-1)70,250,70                                          
C     P(N) NOW IN A3, READ P(N-1) INTO A2                                                                    
   70 CONTINUE
C     FORM P(N)-P(N-1) IN A2                                       
      DO 80 I=1,NTT2                                                   
   80 A2(I)=A3(I)-A2(I)                                                
C     FIND LENGTH DP1                                                  
      SP11=TRACE(A2,A2,N,NO)                                           
      IF(SP11.GT.0.) DP1=SQRT(SP11/TWO)                                
      IF(SP11.LE.0.) DP1=0.
C     TEST FOR CONVERGENCE BY FINDING IF ROOT MEAN SQUARE DP IS LESS   
C     THAN ACURCY                                                      
      RMSDP=DP1/N                                                     
      IF(RMSDP-ABS(ACURCY))90,100,100                                  
C     SETS EXIT FLAG (IT) FOR SCF                                     
C     IF CONVERGENCE, PUT THE ITERATION COUNTER AT ITS LIMIT         
   90 IT=25                                                            
      GO TO 250                                                        
  100 IF(ICOUNT-2)120,110,120                                          
  110 CONTINUE                                                          
C      WRITE(IOUT,10) RMSDP                                             
      GO TO 245                                                                                               
  120 CONTINUE
      IF(ICOUNT-4) 140,130,130                                                                                  
  130 CONTINUE
      SP23=SP12                                                       
      SP33=SP22                                                        
      SP13=TRACE(A3,A2,N,NO)                                          
C     FIND LENGTH DP3                                                   
      IF(SP33.GT.0.) DP3=SQRT(SP33/TWO)                                
      IF(SP33.LE.0.) DP3=0. 
C     READ P(N-1)-P(N-2) INTO A3                                                                                
  140 CONTINUE
      SP12=TRACE(A3,A2,N,NO)                                           
      SP22=TRACE(A3,A3,N,NO)                                            
C     FIND LENGTH DP2                                                  
      IF(SP22.GT.0.) DP2=SQRT(SP22/TWO)                                
      IF(SP22.LE.0.) DP2=0.
C     FIND COSINE OF ANGLE BETWEEN SUCCESSIVE DISPLACEMENTS             
      COSPHI=SP12/(TWO*DP1*DP2)                                        
      IF(ICOUNT-3)160,150,160                                        
  150 WRITE(IOUT,10) RMSDP,COSPHI                                     
      GO TO 240                                                       
C     FIND COSINE OF ANGLE BETWEEN DP3 AND PLANE OF DP1 AND DP2       
  160 X=(SP13*SP22-SP12*SP23)/(SP11*SP22-SP12*SP12)                     
      Y=(SP23*SP11-SP12*SP13)/(SP11*SP22-SP12*SP12)                     
      COSPSI=SQRT((X*X*SP11+Y*Y*SP22+TWO*X*Y*SP12)/TWO)/DP3          
C      WRITE(IOUT,10) RMSDP,COSPHI,COSPSI                                                                       
C     DO NOT EXTRAPOLATE UNLESS 4 CONSECUTIVE POINTS ARE NEARLY COPLANAR
      IF(COSPSI-PT99)240,240,170                                       
C     EXPRESS VECTOR DP1 AS X*DP3(PROJECTED)+Y*DP2                     
  170 Y=-Y/X                                                          
      X=ONE/X                                                        
C     TEST IF 2*2 MATRIX HAS REAL EIGENVALUES BETWEEN -.95 AND +.95   
      XY=Y*Y+FOUR*X                                                  
      IF(XY)190,180,180                                                 
  180 XY=ABS(Y)+SQRT(XY)                                               
      IF(XY-ONEPT9)220,220,190                                         
C     IF 4-POINT EXTRAPOLATION IS NOT POSSIBLE TRY 3-POINT          
  190 IF(ABS(COSPHI)-PT995)240,240,200                                
  200 X=DP1/(DP2*COSPHI-DP1)                                           
      DO 210 I=1,NTT2                                                 
  210 A3(I)=A3(I)+X*A2(I)                                                                                           
      ICOUNT=0                                                       
      GO TO 240                                                     
  220 XXX=X/(ONE-X-Y)                                                
      YYY=(X+Y)/(ONE-X-Y)                                             
      DO 230 I=1,NTT2                                                 
  230 A3(I)=A3(I)+YYY*A2(I)                                                                                
      DO 231 I=1,NTT2                                                 
  231 A3(I)=A3(I)+XXX*A2(I)                                                      
      ICOUNT=0                                                                                                   
      GO TO 260                                                                                            
  240 CONTINUE
      GO TO 260                                                                                           
  245 CONTINUE                                                       
  250 CONTINUE
  260 CONTINUE                                                       
      RETURN                                                           
      END                                                             
      FUNCTION TRACE(A,B,N,NO)                                        
C     TRACE OF A PRODUCT OF DOUBLE SYMMETRIC MATRICES A AND B STORED   
C     AS LINEAR VECTORS                                                
C     NO=1 FOR CLOSED SHELLS, =2 FOR OPEN SHELLS                      
      DIMENSION A(640800),B(640800)                                   
      REAL ZERO,TWO
      DATA ZERO/0./,TWO/2./                                          
      TRACE=ZERO                                                        
      K=0                                                              
      DO 3 L=1,NO                                                     
      DO 2 J=1,N                                                       
      DO 1 I=1,J                                                       
      K=K+1                                                             
    1 TRACE =TRACE +TWO*A(K)*B(K)                                      
    2 TRACE =TRACE -A(K)*B(K)                                          
    3 CONTINUE                                                       
      RETURN                                                         
      END                                                              
      SUBROUTINE CPRINT(IOP1,IOP2)                                   
C     CALCULATION OF MOLECULAR PROPERTIES                               
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)  
      COMMON/ARRAYS/A(800,800),AX(800),B(800,800),BX(800)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      INTEGER AN,U,EL,OCCA,OCCB,ORB,ANI,CHARGE
      CHARACTER*4 IOP1,IOP2
      CHARACTER*1 STAR
      REAL CP65,C4P95,C10P27,CP325,C2P54,C7P33
      REAL ZERO,ONE,TWO,THREE,FIVE,SEVEN,SQRT3
      DIMENSION DPM(3),DM(3),DMSP(3),DMPD(3),CIS(18),ATENG(18,2),
     1 CHRG1OLD(200),DUMMY(800)
      DATA CIS/539.86,0.,0.,0.,0.,820.10,379.36,888.68,44829.2,       
     1 0., 0., 0., 0., 0., 0., 0., 0., 0./
      DATA IW/7/,AUKCAL/627.088/                                      
      DATA ZERO/0./,ONE/1./,TWO/2./,THREE/3./,FIVE/5./,SEVEN/7./      
      DATA CP65/.65/,C4P95/4.95/,C10P27/10.27175/,CP325/.325/       
      DATA C2P54/2.5416/,C7P33/7.33697/                              
      ETOTAL=ENERGY+ETOTAL                                            
      EBIND=ETOTAL                                                  
C      J=1                                                             
      IF(IOP1.EQ.'CNDO') GO TO 5                                        
C      J=2                                                               
    5 CONTINUE                                                         
C ATENG NOT USED
      DO 10 I=1,NATOMS                                                
      K=AN(I)                                                        
C      ETOTAL=ETOTAL+ATENG(K,J)
C   10 EBIND=EBIND-ATENG(K,J)                                        
   10 EBIND=EBIND
      EBIND= -EBIND*27.2114
CC      IF(DEL.EQ.0.) WRITE(IW,60) ENERGY                              
      RYE=EBIND -BBB
CC      IF(DEL.EQ.0.) WRITE(IW,62) RYE                                
CC      IF(DEL.EQ.0.) WRITE(8,62)  RYE 
CC      IF(DEL.EQ.0.) WRITE(IW,61) ETOTAL                              
C      EKCAL= AUKCAL*ETOTAL                                           
C      WRITE(IW,64) EKCAL                                             
   60 FORMAT(10X,'ELECTRONIC ENERGY = ',F16.10,' HARTREE')           
   61 FORMAT(10X,'TOTAL ENERGY      = ',F16.10,' HARTREE')          
   62 FORMAT(4X,'BINDING ENERGY    = ',F10.4,' eV ')                                                                                 
C   64 FORMAT(10X,22H ENERGY IN KCAL/MOLE = ,F16.7)                   
      IF(IOP2.EQ.'CLSD') GO TO 30                                      
C***** THIS SECTION OF CODING SHOULD BE REMOVED IF B ARRAY IS REMOVED   
C     FOR OPEN-SHELL CALCULATION, COMPUTE P(TOTAL) AND P(SPIN) FROM    
C     P(ALPHA) AND P(BETA)                                            
C     P(TOTAL) WILL BE IN A AND P(SPIN) IN B                          
      DO 20 I=1,N                                                   
      DO 20 J=1,N                                                     
      TEMP=A(I,J)-B(I,J)                                             
      A(I,J)=TEMP+TWO*B(I,J)                                          
      B(I,J)=TEMP                                                    
   20 CONTINUE                                                        
CC      IF(DEL.EQ.0.) WRITE(IW,1000)                                    
 1000 FORMAT(21H TOTAL DENSITY MATRIX)                                
      CALL MATOUT(A,DUMMY(N),0)                                        
CC      WRITE(IW,1001)                                                 
 1001 FORMAT(21H  SPIN DENSITY MATRIX)                              
      CALL MATOUT(B,DUMMY(N),0)                                        
C***** REMOVE THE ABOVE CODING IF THE B ARRAY IS REMOVED             
   30 CONTINUE                                                        
      DO 120 I=1,3                                                  
      DM(I)=ZERO                                                      
      DMSP(I)=ZERO                                                    
  120 DMPD(I)=ZERO                                                      
C     THIS SECTION CALCULATES ELECTRON DENSITIES FOR ATOMS, DIPOLE      
C     MOMENT CONTRIBUTIONS FROM EACH ATOM, AND HYPERFINE COUPLING      
C     CONSTANTS IN OPEN-SHELL MOLECULES.                                
      J=0                                                           
      INDEX=1                                                       
      SQRT3=SQRT(3.)                                                 
C      MLL=1
      DO 200 I=1,NATOMS                                               
      CHRG1OLD(I)=CHRG(I)
      IF(IDOUBLE.EQ.1) CHRG1OLD(I)=0.
      ANI=AN(I)                                                     
CC      WRITE(IW,100) I,EL(ANI)                                     
      TCHG=ZERO                                                     
      NORBI=NORB(ANI)                                                 
      IF(IOP2.EQ.'OPEN') GO TO 115                                   
      DO 90 K=1,NORBI                                                
      J=J+1                                                       
      TCHG=TCHG+A(J,J)                                               
CCC      IF(DEL.EQ.0.) WRITE(IW,101) J,ORB(K),A(J,J)                    
   90 CONTINUE                                                        
      CHRG(I)= CZ(ANI)-TCHG                                               
C  CHRG1 IS THE TOTAL POSITIVE CHARGE SEPARATED WITHIN THE MOLECULE:
      CHRG(I)= (CHRG(I)+FLOAT(2*IDOUBLE-2)*CHRG1OLD(I))/
     1 FLOAT(2*IDOUBLE-1)
CCC      IF(DEL.EQ.0.) WRITE(IW,102) TCHG,CHRG(I)                      
C      IF(DEL.EQ.0.) WRITE(8,102) TCHG,CHRG(I)
      GO TO 110                                                         
  115 CONTINUE                                                          
C***** THIS SECTION OF CODING SHOULD BE REMOVED IF B ARRAY IS REMOVED   
      DO 91 K=1,NORBI                                                  
      J=J+1                                                             
      TCHG=TCHG+A(J,J)                                                     
C     IF(DEL.EQ.0.) WRITE(IW,103) J,ORB(K),A(J,J),B(J,J)              
   91 CONTINUE                                                            
      CHRG(I)= CZ(ANI)-TCHG
      CHRG(I)= (CHRG(I)+FLOAT(2*IDOUBLE-2)*CHRG1OLD(I))/
     1 FLOAT(2*IDOUBLE-1)
C      IF(DEL.EQ.0.) WRITE(IW,102) TCHG,CHRG(I)
C      IF(DEL.EQ.0.) WRITE(8,102) TCHG,CHRG(I)
      K=J+1-NORBI                                                     
CC      HFC=CIS(ANI)*B(K,K)                                             
C  CIS VALUES ARE USED TO CALCULATE ESR PARAMETERS (POPLE AND BEVERIDGE pg 131)
CC      IF(DEL.EQ.0.) WRITE(IW,104) TCHG,CHRG(I),HFC                   
CC      IF(DEL.EQ.0.) WRITE(8,104) TCHG,CHRG(I),HFC
C***** REMOVE THE ABOVE CODING IF THE B ARRAY IS REMOVED               
  100 FORMAT(2X,I3,2X,A4)                                              
  101 FORMAT(8X,I3,2X,A4,2X,F13.10)                                    
  102 FORMAT(2X,24HVALENCE ELECTRON DENSITY,F8.4,8H, CHARGE,F8.4)     
  103 FORMAT(8X,I3,2X,A4,2(2X,F13.10))                                  
  104 FORMAT(2X,24HVALENCE ELECTRON DENSITY,F8.4,8H, CHARGE,F8.4,      
     1 6H   A =,F10.3)                                                 
  110 CONTINUE                                                        
      IF(ANI.EQ.1) GO TO 180                                          
      IF(ANI.EQ.2) GO TO 180
      IF(ANI.GT.2) GO TO 140                                           
C  FOR SECOND-ROW ELEMENTS WITH D ORBITALS (COMMENTED OUT):
C  160 SLTR1=(CP65 * FLOAT(ANI)-C4P95)/THREE                            
C      FACTOR=C2P54  *SEVEN/(SQRT(5.)*SLTR1)                           
C      DO 170 K=1,3                                                      
C  170 DMSP(K)=DMSP(K)-A(INDEX,INDEX+K)*C10P27/SLTR1                    
C      IF(NORBI.EQ.4) GO TO 180                                         
C      DMPD(1)=DMPD(1)-FACTOR*(A(INDEX+2,INDEX+8)+A(INDEX+3,INDEX+5)     
C     1  +A(INDEX+1,INDEX+7)-ONE / SQRT3     *A(INDEX+1,INDEX+4))        
C      DMPD(2)=DMPD(2)-FACTOR*(A(INDEX+1,INDEX+8)+A(INDEX+3,INDEX+6)     
C     1  -A(INDEX+2,INDEX+7)-ONE / SQRT3     *A(INDEX+2,INDEX+4))       
C      DMPD(3)=DMPD(3)-FACTOR*(A(INDEX+1,INDEX+5)+A(INDEX+2,INDEX+6)     
C     1  +2.00/ SQRT3     *A(INDEX+3,INDEX+4))                        
C      GO TO 180                                                       
  140 CONTINUE                                                          
      DO 150 K=1,3                                                     
C CHANGED TO USE MU OF 2P ORBITAL: EMU(ANI,0,2,0.5*FLOAT(CHARGE),STAR(I)) 
      IF(ANI.GE.4) DENOM=EMU(ANI,0,2,0.,STAR(I))
      IF(ANI.GE.4) DMSP(K)=DMSP(K)-A(INDEX,INDEX+K)*C7P33/DENOM
C   THE FACTOR 7.3370 IS FROM POPLE AND SEGAL, J CHEM PHYS, VOL 43
C   (1965), EQN 4.2 ON PAGE S140
C  SLATER EFFECTIVE NUCLEAR CHARGE HAD BEEN: /(CP325*FLOAT(ANI-1)). THIS
C  DENOMINATOR, WHEN DIVIDED BY 2, GIVES THE EXPONENT EMU OF 0.325*(ATOMIC#-1) 
  150  CONTINUE
  180 DO 190 K=1,3                                                   
  190 DM(K)=DM(K)+CHRG(I)*C2P54*C(I,K)
      INDEX=INDEX+NORBI                                              
  200 CONTINUE                                                       
      DO 210 I=1,3                                                   
  210 DPM(I)=DM(I)+DMSP(I)+DMPD(I)                                   
C     PRINT OUT THE DIPOLE MOMENT.                                                                      
  220 FORMAT(20X,16H  DIPOLE MOMENTS)                                                                                  
  230 FORMAT(5X,11H COMPONENTS,3X,2H X,8X,2H Y,8X,2H Z)             
CC      IF(DEL.EQ.0.) WRITE(IW,240) DM(1),DM(2),DM(3)                  
  240 FORMAT(5X,10H DENSITIES,3(1X,F9.5))                           
CC      IF(DEL.EQ.0.) WRITE(IW,250) DMSP(1),DMSP(2),DMSP(3)           
  250 FORMAT(5X,4H S.P,6X,3(1X,F9.5))                                 
CC      IF(DEL.EQ.0.)  WRITE(IW,260) DMPD(1),DMPD(2),DMPD(3)       
  260 FORMAT(5X,4H P.D,6X,3(1X,F9.5))                                
CC      IF(DEL.EQ.0.) WRITE(IW,270) DPM(1),DPM(2),DPM(3)             
  270 FORMAT(5X,6H TOTAL,4X,3(1X,F9.5))                               
      XYQ=DPM(1)*DPM(1)+DPM(2)*DPM(2)+DPM(3)*DPM(3)
C MODIFICATION TO DIPOLE MOMENT CALCULATION !!!!!:
      DP= 0.50*SQRT(XYQ)
C      IF(DEL.EQ.0.) WRITE(8,105) RYE,DP,CHRG1,BBB
  105 FORMAT(5X,'DE= ',F7.4,7X,'DP= ',F9.4,3X,'CHRG1= ',F7.4,3X,
     1 'BBB= ',F8.4)
CC      IF(DEL.EQ.0.) WRITE(IW,280) CHRG1,DP
  280 FORMAT(3X,'CHRG1= ',F8.4,3X,15H DIPOLE MOMENT=,F9.5,7H DEBYES)   
      RETURN                                                         
      END                                                             
      SUBROUTINE EIGN1M(N,RHO,A,EIG)                                   
C     A GIVENS TYPE DIAGONALIZATION ROUTINE USING THE HOUSEHOLDER AND   
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(800,800),EIG(800)
      DIMENSION GAMMA(800),BETA(800),BETASQ(800),W(800)
C     THE FOLLOWING DIMENSIONED VARIABLES ARE EQUIVALENCED
      DIMENSION P(800),Q(800)
      DIMENSION IPOSV(800),IVPOS(800),IORD(800)
      REAL ZERO,ONE,HALF,TWO,DEC12                                   
      EQUIVALENCE (P(1),BETA(1)),(Q(1),BETA(1))                     
      EQUIVALENCE (IPOSV(1),GAMMA(1)),(IVPOS(1),BETA(1)),          
     1(IORD(1),BETASQ(1))                                           
      DATA ZERO/0./,ONE/1./,HALF/.5/,TWO/2./,DEC12/1.E-12/           
      RHOSQ=RHO*RHO                                                 
C PART I***        HOUSEHOLDER TRIDIAGONALIZATION PROCEDURE          
      IF(N.EQ.0) GO TO 600                                           
      N1=N-1                                                         
      N2=N-2                                                          
      GAMMA(1)=A(1,1)                                                
      IF(N2) 200,190,40                                               
   40 DO 180 NR=1,N2                                                    
      B=A(NR+1,NR)                                                     
      S=ZERO                                                         
      DO 50 I=NR,N2                                                     
   50 S=S+A(I+2,NR)**2                                               
C     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION                  
      A(NR+1,NR)=ZERO                                                
      IF (S) 170,170,60                                          
   60 S=S+B*B                                                    
      SGN=ONE                                                  
      IF (B) 70,80,80                                                 
   70 SGN=-ONE                                                    
   80 IF(S.GT.0.) SQRTS= SQRT(S)                                    
      IF(S.LE.0.) SQRTS= 0.
      D=SGN/(SQRTS+SQRTS)                                            
      IF((HALF+B*D).GT.0.) TEMP= SQRT(HALF+B*D)                    
      IF((HALF+B*D).LE.0.) TEMP=0.
      W(NR)=TEMP                                                       
      A(NR+1,NR)=TEMP                                                 
      D=D/TEMP                                                      
      B=-SGN*SQRTS                                                      
C     D IS FACTOR OF PROPORTIONALITY. NOW COMPUTE AND SAVE W VECTOR.   
C     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED.              
      DO 90 I=NR,N2                                                    
      TEMP=D*A(I+2,NR)                                                 
      W(I+1)=TEMP                                                      
   90 A(I+2,NR)=TEMP                                                  
C     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.             
C     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)        
      WTAW=ZERO                                                    
      DO 140 I=NR,N1                                                  
      SUM=ZERO                                                    
      DO 100 J=NR,I                                                 
  100 SUM=SUM+A(I+1,J+1)*W(J)                                         
      I1=I+1                                                      
      IF(N1-I1) 130,110,110                                         
  110 DO 120 J=I1,N1                                             
  120 SUM=SUM+A(J+1,I+1)*W(J)                                        
  130 P(I)=SUM                                                         
  140 WTAW=WTAW+SUM*W(I)                                            
C     P VECTOR AND SCALAR K  NOW STORED. NEXT COMPUTE Q VECTOR      
      DO 150 I=NR,N1                                                
  150 Q(I)=P(I)-WTAW*W(I)                                           
C     NOW FORM PAP MATRIX, REQUIRED PART                            
      DO 160 J=NR,N1                                                 
      QJ=Q(J)                                                        
      WJ=W(J)                                                       
      DO 160 I=J,N1                                                 
      A(I+1,J+1)=A(I+1,J+1)-TWO *(W(I)*QJ+WJ*Q(I))                    
  160 CONTINUE
  170 BETA(NR)=B                                                     
      BETASQ(NR)=B*B                                                 
  180 GAMMA(NR+1)=A(NR+1,NR+1)                                     
  190 B=A(N,N-1)                                                   
      BETA(N-1)=B                                                      
      BETASQ(N-1)=B*B                                                 
      GAMMA(N)=A(N,N)                                               
  200 BETASQ(N)=ZERO                                                  
C     PUT AN IDENTITY MATRIX IN THE UPPER TRIANGLE               
      DO 220 I=1,N                                                    
      DO 210 J=I,N                                                   
  210 A(I,J)=ZERO                                                    
  220 A(I,I)=ONE                                                     
C PART II***       FORM PRODUCT P1*P2*...*PN-2                       
      K=N                                                             
  221 KP=K                                                           
      K=K-1                                                        
      KM=K-1                                                     
      DO 222 I=K,N                                               
  222 W(I)=A(I,KM)                                                   
      AKKM=W(K)                                                        
      C2AKKM=-(AKKM+AKKM)                                           
      DO 225 J=KP,N                                                  
      SUM=ZERO                                                        
      DO 223 L=K,N                                                    
  223 SUM=SUM+W(L)*A(L,J)                                           
      SUM=SUM+SUM                                                 
      DO 224 I=K,N                                                    
      A(I,J)=A(I,J)-SUM*W(I)                                         
  224 CONTINUE
  225 CONTINUE                                                       
      DO 226 I=K,N                                                
  226 A(I,K)=C2AKKM*W(I)                                             
      A(K,K)=ONE+A(K,K)                                          
      IF(K.GT.2) GO TO 221                                            
      DO 227 I=2,N                                                    
  227 A(I,1)=ZERO                                                    
C PART III***      QR ALGORITHM                                       
      M=N                                                             
      SUM=ZERO                                                        
      SHIFT=ZERO                                                     
      GO TO 350                                                      
  230 SUM=SUM+SHIFT                                                
      COSA=ONE                                                     
      G=GAMMA(1)-SHIFT                                           
      PP=G                                                           
      PPBS=PP*PP+BETASQ(1)                                     
      IF(PPBS.GT.0.) PPBR= SQRT(PPBS)                             
      IF(PPBS.LE.0.) PPBR=0.
      DO 320 J=1,M                                              
      COSAP=COSA                                                    
      IF (PPBS.GT. DEC12) GO TO 250                                 
  240 SINA=ZERO                                                     
      SINA2=ZERO                                                    
      COSA=ONE                                                      
      GO TO 290                                                    
  250 SINA=BETA(J)/PPBR                                              
      SINA2=BETASQ(J)/PPBS                                           
      COSA=PP/PPBR                                                    
C      POSTMULTIPLY P1*...*PN-2 BY ROTATION MATRIX  (P TRANSPOSE)     
      DO 260 I=1,N                                                    
      VTEMP=A(I,J)                                                     
      VJTEMP=A(I,J+1)                                                
        A(I,J)=COSA*VTEMP+SINA*VJTEMP                                  
        A(I,J+1)=-SINA*VTEMP+COSA*VJTEMP                            
  260 CONTINUE                                                    
  290 DIA=GAMMA(J+1)-SHIFT                                              
      U=SINA2*(G+DIA)                                                
      GAMMA(J)=G+U                                                   
      G=DIA-U                                                         
      PP=DIA*COSA-SINA*COSAP*BETA(J)                                    
      IF(J .NE. M) GO TO 310                                           
  300 BETA(J)=SINA*PP                                                 
      BETASQ(J)=SINA2*PP*PP                                        
      GO TO 330                                                   
  310 PPBS=PP*PP+BETASQ(J+1)                                      
      IF(PPBS.GT.0.) PPBR= SQRT(PPBS)                               
      IF(PPBS.LE.0.) PPBR=0.
      BETA(J)=SINA*PPBR                                            
  320 BETASQ(J)=SINA2*PPBS                                         
  330 GAMMA(M+1)=G                                                
C     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT                  
  335 IF(BETASQ(M)-RHOSQ)  340,340,370                            
  340 EIG(M+1)=GAMMA(M+1)+SUM                                         
  350 BETA(M)=ZERO                                                    
      BETASQ(M)=ZERO                                                   
      M=M-1                                                          
      IF(M) 400,400,335                                                
C     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE    
C     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT                       
  370 A2=GAMMA(M+1)                                                  
      R2=HALF *A2                                                     
      R1=HALF *GAMMA(M)                                               
      R12=R1+R2                                                       
      DIF=R1-R2                                                       
      IF((DIF*DIF+BETASQ(M)).GT.0.) TEMP= SQRT(DIF*DIF+BETASQ(M))   
      IF((DIF*DIF+BETASQ(M)).LE.0.) TEMP=0. 
      R1=R12+TEMP                                                  
      R2=R12-TEMP                                                   
      DIF= ABS(A2-R1)- ABS(A2-R2)                                    
      IF(DIF .LT. ZERO) GO TO 390                                   
  380 SHIFT=R2                                                       
      GO TO 230                                                      
  390 SHIFT=R1                                                      
      GO TO 230                                                      
  400 EIG(1)=GAMMA(1)+SUM                                          
      IF(N1.EQ.0) GO TO 600                                          
C PART IV***       SORTING AND ORDERING OF EIGENVALUES                
C     INITIALIZE AUXILIARY TABLES REQUIRED FOR REARRANGING THE VECTORS 
      DO 410 J=1,N                                                    
      IPOSV(J)=J                                                  
      IVPOS(J)=J                                                   
  410 IORD(J)=J                                                        
C     USE A TRANSPOSITION SORT TO ORDER THE EIGENVALUES              
      M=N                                                           
      GO TO 450                                                    
  420 DO 440 J=1,M                                               
      IF (EIG(J) .LE. EIG(J+1)) GO TO 440                          
  430 TEMP=EIG(J)                                                      
      EIG(J)=EIG(J+1)                                                
      EIG(J+1)=TEMP                                                  
      ITEMP=IORD(J)                                                  
      IORD(J)=IORD(J+1)                                             
      IORD(J+1)=ITEMP                                               
  440 CONTINUE                                                    
  450 M=M-1                                                         
      IF(M .NE. 0) GO TO 420                                      
C PART V***        ORDERING OF EIGENVECTORS                          
  470 DO 500 L=1,N1                                              
      NV=IORD(L)                                                    
      NP=IPOSV(NV)                                                  
      IF(NP .EQ. L) GO TO 500                                         
  480 LV=IVPOS(L)                                                     
      IVPOS(NP)=LV                                                 
      IPOSV(LV)=NP                                                     
      DO 490 I=1,N                                                      
      TEMP=A(I,L)                                                    
      A(I,L)=A(I,NP)                                                   
  490 A(I,NP)=TEMP                                                   
  500 CONTINUE                                                         
  600 RETURN                                                           
      END                                                              
      SUBROUTINE MATOUT(A,MDIAG,LOP)                                  
C     THIS SUBROUTINE PRINTS OUT SQUARE ARRAYS                          
C     IF LOP = 0 AN NXN ARRAY IS PRINTED WITH ORBITAL LABELS.          
C     IF LOP = 1 THE CONTENTS OF MDIAG (USUALLY EIGENVALUES) ARE       
C                PRINTED ALONG THE TOP.                               
C     IF LOP GT 1 AN LOP X LOP ARRAY IS PRINTED WITH NO LABELS.         
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      COMMON/SCLRS/N,NATOMS,OCCA,OCCB,ENERGY,ETOTAL,BBB,CHARGE,
     1 STAR(200)
      COMMON/SSS/RYE,DP,CHRG1,DEL,XBO,ACHR,CHRG(200),IDOUBLE
      DIMENSION A(800,800),MDIAG(800)
      INTEGER ANII,ORB,EL,AN,U,OCCA,OCCB,CHARGE                    
      CHARACTER*1 STAR
      REAL MDIAG
      DATA IW/7/                                                     
      IF(LOP.GT.1) GO TO 200                                          
      DO 110 M=1,N,11                                                
      K=M+10                                                         
      IF (K.LE.N) GO TO 30                                           
   20 K=N                                                              
   30 IF (LOP.EQ.0) GO TO 50                                          
CC      IF(DEL.EQ.0.) WRITE(IW,10) (MDIAG(I),I=M,K)                   
   50 CONTINUE
CC      IF(DEL.EQ.0.) WRITE(IW,60) (I,I=M,K)                           
      I=0                                                            
      DO 110 II=1,NATOMS                                                                                               
      ANII=AN(II)                                                   
      NORBI=NORB(ANII)                                                
      DO 105 L=1,NORBI                                                 
      I=I+1                                                            
   70 CONTINUE
CC      IF(DEL.EQ.0.) WRITE(IW,80) I,II,EL(ANII),ORB(L),(A(I,J),J=M,K)  
  105 CONTINUE
  110 CONTINUE                                                        
      GO TO 300                                                        
  200 CONTINUE                                                       
      DO 210 M=1,LOP,11                                              
      K=M+10                                                         
      IF(K.LE.LOP) GO TO 230                                            
      K=LOP                                                           
  230 CONTINUE                                                   
CC      IF(DEL.EQ.0.) WRITE(IW,61)(I,I=M,K)                           
      DO 210 I=1,LOP                                                 
CC      IF(DEL.EQ.0.) WRITE(IW,81) I,(A(I,J),J=M,K)                    
  210 CONTINUE                                                       
  300 CONTINUE                                                   
   10 FORMAT(15X,11(F9.3))                                         
   60 FORMAT(13X,10I9)                                               
   61 FORMAT(1X,10I9)                                                  
   80 FORMAT(1X,I2,I3,A4,1X,A4,10(F9.4))                            
   81 FORMAT(1X,I2,10(F9.3))                                         
  100 FORMAT(1X)                                                     
      RETURN                                                        
      END                                                            
      FUNCTION EUA(NATOMS)
C  CALCULATES TOTAL ELECTRON ENRGY OF UNBOUND ATOMS IN MOLECULE
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)                
      COMMON/MOLINF/AN(200),CZ(18),C(200,3),U(800),NORB(18),EL(18),
     1 ORB(9)
      INTEGER AN,U,EL,ORB      
      FH=   13.598435
      FC= 1030.1085
      FN= 1486.058
      FO= 2043.8429
      FF= 2715.891
      EUA=0.
      DO 5 I=1,NATOMS
      IF(AN(I).EQ.1) EUA=EUA+FH
      IF(AN(I).EQ.6) EUA=EUA+FC
      IF(AN(I).EQ.7) EUA=EUA+FN
      IF(AN(I).EQ.8) EUA=EUA+FO
      IF(AN(I).EQ.9) EUA=EUA+FF
   5  CONTINUE
      RETURN
      END
      FUNCTION EMU(II,KK,LZ,Q,STARX)
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 STARX,STAR
      DIMENSION A(10,3),Z(3)
C     1 Bohr unit is 0.529 A 
C     ATOMIC CHARGES QQ ARE ALL ZERO:      
      QQ=0.         
      B1=  0.3125
      B2=  0.030
      B3=  0.006
      B4=  0.0221
      B5=  0.0074
      B6=  0.6989
      B7=  0.3008
      B8=  0.2299
      B9=  0.9270
      B10= 0.3668
      B11= 0.3526
      B12= 0.1618
      A(1,1)= 1. 
      A(1,2)= 0.
      A(1,3)= 0.
      A(2,1)= 2. 
      A(2,2)= 0.
      A(2,3)= 0.
      A(3,1)= 2.
      A(3,2)= 1. 
      A(3,3)= 0.
      A(4,1)= 2.
      A(4,2)= 1.      
      A(4,3)= 1.      
      A(5,1)= 2.
      A(5,2)= 2.
      A(5,3)= 1.
      A(6,1)= 2.
      A(6,2)= 2.
      A(6,3)= 2.
      A(7,1)= 2.
      A(7,2)= 2.
      A(7,3)= 3.
      A(8,1)= 2.
      A(8,2)= 2.
      A(8,3)= 4.
      A(9,1)= 2.
      A(9,2)= 2.
      A(9,3)= 5.           
C
      A(10,1)=2.
      A(10,2)=2.
      A(10,3)=6.
      T1S=0.
      T2S=0.
      T2P=0.
      IF(II.LE.2) T1S=A(II,1) -QQ
      IF(II.GT.2) T1S=A(II,1)
      IF(II.LE.2) T2S=A(II,2)
C      IF(II.EQ.3) T2S=A(II,2) -QQ
C      IF(II.EQ.4) T2S=A(II,2)
C      IF(II.EQ.4) T2P=A(II,3) -QQ
C      IF(II.EQ.4.AND.QQ.GT.1.) T2P= 0.
C      IF(II.EQ.4.AND.QQ.GT.1.) T2S=A(II,2) -(QQ-1.)
C      IF(II.LT.4) T2P=A(II,3)           
      IF(II.GE.4) T2S=A(II,2)       
      IF(II.GE.4) T2P=A(II,3) -QQ
      Z(1)=0.
      Z(2)=0.
      Z(3)=0.
      Z1=FLOAT(II)
      IF(II.LE.2.AND.T1S.LE.1.) Z(1)=Z1
      IF(II.LE.2.AND.T1S.GT.1.) Z(1)=Z1-(B1-B2/Z1-B3/(Z1*Z1))*(T1S-1.D0) 
C      IF(II.EQ.3) Z(1)=Z1-(B1-B2/Z1-B3/(Z1*Z1))*(T1S-1.)-B4*T2S
C      IF(II.EQ.4) Z(1)=Z1-(B1-B2/Z1-B3/(Z1*Z1))*(T1S-1.)-B4*T2S
      IF(II.GE.4) Z(1)=Z1-(B1-B2/Z1-B3/(Z1*Z1))*(T1S-1.)-B4*T2S-B5*T2P
C      IF(II.EQ.3.AND.T2S.LE.1.) Z(2)=Z1-B6*T1S
C      IF(II.EQ.3.AND.T2S.GT.1.) Z(2)=Z1-B6*T1S-B7*(T2S-1.)
C      IF(II.EQ.4.AND.T2S.GT.1.) Z(2)=Z1-B6*T1S-B7*(T2S-1.)
C      IF(II.EQ.4.AND.T2S.LE.1.) Z(2)=Z1-B6*T1S
      IF(II.GE.4.AND.T2S.GT.1.) Z(2)=Z1-B6*T1S-B7*(T2S-1.)-B8*T2P  
      IF(II.GE.4.AND.T2S.LE.1.) Z(2)=Z1-B6*T1S            -B8*T2P
C      IF(II.EQ.3)               Z(3)=Z1-B9*T1S-B10*1.0
C      IF(II.EQ.4)             Z(3)=Z1-B9*T1S-B10*T2S
      IF(II.GE.4.AND.T2P.LE.1.) Z(3)=Z1-B9*T1S-B10*T2S
      IF(II.GE.4.AND.T2P.GT.1.) Z(3)=Z1-B9*T1S-B10*T2S-(B11-B12/Z1)*
     1 (T2P-1.)      
      IF(LZ.EQ.0) EMU=Z(1)
      IF(LZ.EQ.1) EMU=Z(2)/2.
      IF(LZ.EQ.2) EMU=Z(3)/2.
      IF(EMU.LT.0.) EMU=0.
      RETURN
      END
C FUNCTION UX6
      FUNCTION UX6(NA,STR1,NB,STR2,DD)
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*1 STR1,STR2
C  ENERGIES ARE FIRST CALCULATED IN HARTREES (27.2114 eV/HARTREE) AND THEN
C  CONVERTED TO EVs
C  NA MUST NOT BE LARGER THAN NB WHEN THIS FUNCTION IS CALLED       
C  THESE VALUES OF EXP ARE THE EMU VALUES (Zeff/2) FOR 2S ORBITALS
C  OF ATOMs OF C, N, O, AND F. EXPH OF 1S FOR H IS 1.0000
      I=NA*NB
      EXP1=  1.0000
      EXPH=    DD
      EXP6S= 1.9208
      EXP6S5=1.9208*1.9208*1.9208*1.9208*1.9208*DD*DD*DD
      EXP7S= 2.3059
      EXP7S5=2.3059*2.3059*2.3059*2.3059*2.3059*DD*DD*DD
      EXP8S= 2.6909
      EXP8S5=2.6909*2.6909*2.6909*2.6909*2.6909*DD*DD*DD
      EXP9S= 3.0759
      EXP9S5=3.0759*3.0759*3.0759*3.0759*3.0759*DD*DD*DD
      IF(I.EQ.1) PROD=2.*DD*(EXP1 +EXP1)
      IF(I.EQ.1) UX=  4.1*EXPH*EXPH
      IF(I.EQ.6) PROD=2.*DD*(EXP1 +EXP6S)
      IF(I.EQ.6) UX= 14.*EXPH*EXP6S5
      IF(I.EQ.7) PROD=2.*DD*(EXP1 +EXP7S)
      IF(I.EQ.7) UX= 12.*EXPH*EXP7S5
      IF(I.EQ.8) PROD=2.*DD*(EXP1 +EXP8S)
      IF(I.EQ.8) UX= 18.*EXPH*EXP8S5
      IF(I.EQ.9) PROD=2.*DD*(EXP1 +EXP9S)
      IF(I.EQ.9) UX= 28.*EXPH*EXP9S5
      IF(I.EQ.36) PROD=2.*DD*(EXP6S +EXP6S)
      IF(I.EQ.36) B= EXP6S5*EXP6S5
      IF(I.EQ.36) CA= 15.6
      IF(I.EQ.36) CB=  7.5
      IF(I.EQ.36) CC=  4.7
      IF(I.EQ.36) CAA= CA*CA
      IF(I.EQ.36) CAB= CA*CB
      IF(I.EQ.36) CAC= CA*CC
      IF(I.EQ.36) CBB= CB*CB
      IF(I.EQ.36) CBC= CB*CC
      IF(I.EQ.36) CCC= CC*CC
      IF(I.EQ.36.AND.STR1.EQ.'B'.AND.STR2.EQ.'B') UX=B*CBB
      IF(I.EQ.36.AND.STR1.EQ.'C'.AND.STR2.EQ.'C') UX=B*CCC
      IF(I.EQ.36.AND.STR1.EQ.'A'.AND.STR2.EQ.'A') UX=B*CAA
      IF(I.EQ.36.AND.STR1.EQ.'B'.AND.STR2.EQ.'C') UX=B*CBC
      IF(I.EQ.36.AND.STR1.EQ.'C'.AND.STR2.EQ.'B') UX=B*CBC
      IF(I.EQ.36.AND.STR1.EQ.'B'.AND.STR2.EQ.'A') UX=B*CAB
      IF(I.EQ.36.AND.STR1.EQ.'A'.AND.STR2.EQ.'B') UX=B*CAB
      IF(I.EQ.36.AND.STR1.EQ.'C'.AND.STR2.EQ.'A') UX=B*CAC
      IF(I.EQ.36.AND.STR1.EQ.'A'.AND.STR2.EQ.'C') UX=B*CAC
      IF(I.EQ.42) PROD=2.*DD*(EXP6S +EXP7S)
      IF(I.EQ.42) UX= 36.*EXP6S5*EXP7S5
      IF(I.EQ.48) PROD=2.*DD*(EXP6S +EXP8S)
      IF(I.EQ.48) UX=  98.*EXP6S5*EXP8S5
      IF(I.EQ.54) PROD=2.*DD*(EXP6S +EXP9S)
      IF(I.EQ.54) UX= 900.*EXP6S5*EXP9S5
      IF(I.EQ.49) PROD=2.*DD*(EXP7S +EXP7S)
      IF(I.EQ.49) UX= 42.5*EXP7S5*EXP7S5
      IF(I.EQ.56) PROD=2.*DD*(EXP7S +EXP8S)
      IF(I.EQ.56) UX= 115.*EXP7S5*EXP8S5
      IF(I.EQ.63) PROD=2.*DD*(EXP7S +EXP9S)
      IF(I.EQ.63) UX= 3100.*EXP7S5*EXP9S5
      IF(I.EQ.64) PROD=2.*DD*(EXP8S +EXP8S)
      IF(I.EQ.64) UX= 740.*EXP8S5*EXP8S5
      IF(I.EQ.72) PROD=2.*DD*(EXP8S +EXP9S)
      IF(I.EQ.72) UX=11000.*EXP8S5*EXP9S5
      IF(I.EQ.81) PROD=2.*DD*(EXP9S +EXP9S)
      IF(I.EQ.81) UX=69000.*EXP9S5*EXP9S5
C     DD IS IN ATOMIC UNITS
C     UX6 IS IN UNITS OF eV WHEN MULTIPLIED BY 27.2114 eV/hartree
      UX6= 27.2114*UX*EXP(-PROD)
      RETURN
      END
      
