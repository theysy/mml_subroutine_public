CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     THE ORIGINAL UMAT CODE WAS WRITTEN BY J.-W. LEE (2010. 11)              C
C                                                                             C
C     NAME    : SEONG-YONG  YOON                                              C
C     ADVISOR : FREDERIC BARLAT                                               C
C     E-MAIL  : theysy@postech.ac.kr                                          C
C     GIT-HUB : https://github.com/theysy                                     C
C     AFFILIATION:                                                            C
C         GRADUATE INSTITUTE OF FERROUS TECHNOLOGY (GIFT)                     C
C         POHANG UNIVERSITY OF SCIENCE AND TECHNOLOGY (POSTECH)               C
C                                                                             C
C     VERSION INFORMATION                                                     C
C         1. MML_U0: ONLY FOR PLANE STRESS (SHELL ELEMENT)                    C
C         2. MML_U1: SHELL OR SOLID ELEMENT                                   C
C         3. MML_UF: MML_V1 + FRACTURE MODEL                                  C
C         4. MML_U2: HAH20 AND CHABOCHE WERE UPDATED.                         C
C            DEVELOPMENT LOG:                                                 C
C             -TANGENT MODULUS METHOD IS UPDATED.(190328)                     C
C             -LINE SEARCH EULER BACKWARD IS UPDATED.(190702)                 C
C             -HAH20E IS UPDATED. (200318)                                    C
C             -MODIFIED YOSHIDA-UEMORI IS UPDATED. (200318)                   C
C             -PERMANENT SOFTENING OF HAH20 IS ENHANCED. (20200331)           C
C             -HAH14 IS UPDATED. (20200610)                                   C
C             -HAH20 LATENT HARDENING EQUATION IS UPDATED (20200810)          C
C             -RGBV MODEL IS UPDATED (20201102)                               C
C             -BROYDEN METHOD IS UPDATED (20210113)                           C
C             -DIRECT METHOD IS UPDATED (20210623)                            C
C             -HEXAH MODEL IS UPDATED (20220101)                              C
C                                                                             C
C     REFERENCES                                                              C
C         [1] F.BARLAT  ET AL. IJP     (2003): YLD2000-2D                     C
C         [2] F.BARLAT  ET AL. IJP     (2012): HAH11                          C
C         [3] F.BARLAT  ET AL. IJP     (2014): HAH14                          C
C         [4] J.W. LEE  ET AL. IJP     (2012): CUTTING PLANE METHOD           C
C         [5] J.W. LEE  ET AL. CMAM    (2012): CLOSEST POINT PROJECTION       C
C         [6] J.W. YOON ET AL. IJP     (2004): MULTI-STAGE EULER BACKWARD     C
C                                            : THE ANALYTICAL DERIVATIVE OF   C
C                                              YLD2000_2D MODEL               C
C         [7] F. BARLAT ET AL.  IJP     (2005): YLD2004_18P                   C
C         [8] H. ARETZ          ESAFORM (2007): SCALED FDM                    C
C         [9] T.J. PARK ET AT.  IJSS    (2012): CHABOCHE1-SIMPLER MODEL       C
C         [10]J.Y. LEE  ET AL.  IJSS    (2012): CHABOCHE2                     C
C         [11]F. YOSHIDA ET AL. IJP     (2002): YOSHIDA-UEMORI                C
C         [12]M. ORTIZ  ET AL.  IJNME   (1985): TANGENT MODULUS METHOD        C
C         [13]W.M. SCHERZINGER  CMAME   (2017): LINE SEARCH EULER BACKWARD    C
C         [14]F. BARLAT ET AL.  IJSS    (2020): HAH20H & HAH20E               C
C         [15]E. RAUCH ET AL.   MSMSE   (2011): RGBV                          C
C         [16]K. KITAYAMA ET AL.IJP     (2013): CRYSTALLOGRAPHIC RGBV         C
C         [17]S.Y. YOON ET AL.  CMAME   (2020): LINE-SEARCH FOR HARDENING     C
C         [18]C.G. BROYDEN      MOC     (1965): BROYDEN METHOD                C
C         [19]B. REYNE ET AL.   IJP     (2022): HEXAH (HAH22)                 C
C         [20]S. YOON ET AL.    IJMS    (2023): NSPM                          C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C-----------------------------------------------------------------------
C     STATE VARIABLES FOR UMAT
C     STATEV(1)   : EQUIVALENT PLASTIC STRAIN
C     STATEV(2)   : EQUIVALENT STRESS
C     STATEV(3)   : YOUNG'S MODULUS
C     STATEV(4)   : YIELD CRITERION(0: ELASTIC / 1: PLASTIC)
C     STATEV(5)   : ALOGRITHMIC ITERATION
C     STATEV(6)   : EQUIVALENT PLASTIC STRAIN INCREMENT
C-----------------------------------------------------------------------
C     # STATE VARIABLES OF CHABOCHE
C     STATEV(7:12)    : BACK-STRESS: ALPHA
C-----------------------------------------------------------------------
C     # STATE VARIABLES OF YOSHIDA-UEMORI
C     STATEV(7:12)    : BACK-STRESS: ALPHA
C     STATEV(13:18)   : CENTER OF BOUNDING SURFACE: BETA
C     STATEV(19:24)   : CENTER OF NON-IH SURFACE  : Q
C     STATEV(25)      : RADIUS OF BOUNDING SURFACE: R
C-----------------------------------------------------------------------
C     # STATE VARIABLES IN HAH11 MODELS
C     [1] G-VALUES
C     STATEV(7) : STRAIN-PATH CHANGE PARAMETER
C     STATEV(8) : G(1): G1
C     STATEV(9) : G(2): G2
C     STATEV(10): G(3): G3
C     STATEV(11): G(4): G4
C     [2] MICROSTRUCTURE DEVIATOR
C     STATEV(12): H11
C     STATEV(13): H22
C     STATEV(14): H33
C     STATEV(15): H12
C     STATEV(16): H23
C     STATEV(17): H31
C----------------------------------------------------------------------
C     # STATE VARIABLES IN HAH14 MODELS
C     [1] G-VALUES
C     STATEV(7) : STRAIN-PATH CHANGE PARAMETER
C     STATEV(8) : G(1): G1
C     STATEV(9) : G(2): G2
C     STATEV(10): G(3): G3
C     STATEV(11): G(4): G4
C     STATEV(12): G(5): GS
C     STATEV(13): G(6): GL
C     [2] MICROSTRUCTURE DEVIATOR
C     STATEV(14): H11
C     STATEV(15): H22
C     STATEV(16): H33
C     STATEV(17): H12
C     STATEV(18): H23
C     STATEV(19): H31
C-----------------------------------------------------------------------
C     # STATE VARIABLES IN HAH20 & HAH20 MODELS
C     [1] G-VALUES
C     STATEV(7) : STRAIN-PATH CHANGE PARAMETER
C     STATEV(8) : G(1): G-    *
C     STATEV(9) : G(2): G+    *
C     STATEV(10): G(3): G3    *
C     STATEV(11): G(4): G4    *
C     STATEV(12): G(5): GP    *
C     STATEV(13): G(6): GS    *
C     STATEV(14): G(7): GC    *
C     STATEV(15): G(8): GL    *
C     STATEV(16): G(9): G3*   *
C     STATEV(17): G(10): GP*  *
C     [2] MICROSTRUCTURE DEVIATOR
C     STATEV(18): H11        *
C     STATEV(19): H22        *
C     STATEV(20): H33        *
C     STATEV(21): H12        *
C     STATEV(22): H23
C     STATEV(23): H31
C     STATEV(24): HP11       *
C     STATEV(25): HP22       *
C     STATEV(26): HP33       *
C     STATEV(27): HP12       *
C     STATEV(28): HP23
C     STATEV(29): HP31
C     STATEV(30): HS11       *
C     STATEV(31): HS22       *
C     STATEV(32): HS33       *
C     STATEV(33): HS12       *
C     STATEV(34): HS23
C     STATEV(35): HS31
C-----------------------------------------------------------------------
C     # STATE VARIABLES IN HEXAH MODELS
C     [1] GLOBAL STATE VARIABLES
C     STATEV(7) : STRAIN-PATH CHANGE PARAMETER
C     STATEV(8) : GP
C     STATEV(9) : GL
C     STATEV(10): NSS
C     STATEV(11): PLAS(1): E11/E11   *
C     STATEV(12): PLAS(2): E22/E22   *
C     STATEV(13): PLAS(3): E33/E12   *
C     STATEV(14): PLAS(4): E12
C     STATEV(15): PLAS(5): E23 <-> E23
C     STATEV(16): PLAS(6): E13
C     [2] LOCAL STATE VARIABLES AT PREVIOUSLY ACTIVATED SLIP SYSTEM (NSS=1)
C     STATEV(17) : G(1,1): G-
C     STATEV(18) : G(1,2): G+
C     STATEV(19) : G(1,3): GC
C     STATEV(20) : DEV_H(1,1): H11    *
C     STATEV(21) : DEV_H(1,2): H22    *
C     STATEV(22) : DEV_H(1,3): H33    *
C     STATEV(23) : DEV_H(1,4): H12
C     STATEV(24) : DEV_H(1,5): H23
C     STATEV(25) : DEV_H(1,6): H13
C     STATEV(26) : EQPLASI(1): LOCAL EFFECTIVE STRAIN
C     STATEV(27) : COS_X(1) : LOCAL STRAIN PATH CHANGE PARAMETER
C     [3] LOCAL STATE VARIABLES AT CURRENTLY ACTIVATED SLIP SYSTEM (NSS=2)
C     STATEV(28) : G(2,1): G-
C     STATEV(29) : G(2,2): G+
C     STATEV(30) : G(2,3): GC
C     STATEV(31) : DEV_H(2,1): H11    *
C     STATEV(32) : DEV_H(2,2): H22    *
C     STATEV(33) : DEV_H(2,3): H33    *
C     STATEV(34) : DEV_H(2,4): H12
C     STATEV(35) : DEV_H(2,5): H23
C     STATEV(36) : DEV_H(2,6): H13
C     STATEV(37) : EQPLASI(2): LOCAL EFFECTIVE STRAIN
C     STATEV(38) : COS_X(1) : LOCAL STRAIN PATH CHANGE PARAMETER
C     [3] FOR ADDITIONAL SLIP SYSTEM
C     NES=3 -> STATEV(39:49)....*
C     NES=4 -> STATEV(50:60)
C     NES=5 -> STATEV(61:71)
C-----------------------------------------------------------------------
C     # STATE VARIABLES OF RGBV
C     STATEV(50): RHO_F
C     STATEV(51): RHO_R1(+)
C     STATEV(52): RHO_R2(-)
C     STATEV(53): RHO_L
C     STATEV(54): RHO_F01
C     STATEV(55): RHO_F02
C     STATEV(56): RHO_TOT
C     STATEV(57): DRDE
C     STATEV(58): BACK-STRESS(OPTIONAL)
C-----------------------------------------------------------------------
C     # GRADIENT OF YIELD SURFACE
C     STATEV(59:65): DFDS(1:6)
C     NSTAT_VAR= 65
C-----------------------------------------------------------------------
C     #   UMAT PARAMETERS
C     PROPS(1): OPTIONS FOR ANISOTROPIC HARDENING
C                 0- ISOTROPIC HARDENING
C                 1- CHABOCHE [9][10]
C                 2- YOSHIDA-UEMORI [11] | 2.5- MODIFIED YOSHIDA-UEMORI
C                 3- HAH11 [2]
C                 4- HAH14 [3]
C                 5- HAH20H [14] | 5.5-HAH20E [14]
C                 6- HEXAH (HAH22)
C     PROPS(2): OPTIONS FOR YTELD FUNCTION MODEL
C                 1- VON_MISES
C                 2- HILL 1948
C                 3- YLD2000_2D   [1]
C                 4- YLD2004_18P  [7]
C     PROPS(3): OPTIONS FOR HARDENING INPUT: 
C                 1- SWIFT: SIGMA=P1*(P2+EQPLAS)**P3
C                 2- VOCE: SIGMA=P1-P2*EXP(-P3*EQPLAS)
C                 3- MODIFIED VOCE: SIGMA= P1 + P2*EQPLAS + P3*(1-EXP(-P4*EQPLAS))
C                 4- SWIFT+MODIFIED VOCE: P1*(P2+EQPLAS)**P3 + P4 + P5*EQPLAS + P6*(1-EXP(-P7*EQPLAS))
C                 5- HOCKETT-SHERBY + LINEAR: P1-(P1-P2)*EXP(-P3*EQPLAS**P4)+P5
C                 6- COMBINED SWIFT-VOCE: P1*[P2*(P3+EQPLAS)**P4]+(1-P1)*[P5-P6*EXP(-P7*EQPLAS)]
C                 7- DISLOCATION-BASED HARDENING MODEL: P1*[P2+P3*P4*P5*SQRT(RHO(EQPLAS))
C     PROPS(4): OPTIONS FOR E-MODULUS DEGRADATION INPUT:
C                 1- CONSTANT E
C                 2- CHORD E
C     PROPS(5): OPTIONS FOR STRESS-UPDATE ALGORITHM (ITERATIVE)
C                 1- CUTTING PLANE MEHTOD     (EXPLICIT)
C                 2- CLOSEST POINT PROJECTION (FULLY-IMPLICIT)
C                 3- CLOSEST POINT PROJECTION (SEMI-IMPLICIT)
C                 4- TANGENT MODULUS METHOD   (IMPLICIT)
C                 5- BROYDEN METHOD           (FULLY-IMPLICIT)
C                 SUA_PAR=UDT_PAR+SSC_PAR
C                 EX) 1.55 = CPM + UPDATE METHOD + LINE-SEARCH
C     PROPS(5): OPTIONS FOR STRESS-UPDATE METHOD (NON-ITERATIVE)
C                -1- FIRST ORDER EULER FORWARD METHOD
C                -2- NON-ITERATIVE STRESS UPDATE METHOD (CONTINUUM TANGENT MODULUS)
C                -3- NON-ITERATIVE STRESS UPDATE METHOD (CONSISTENT TANGENT MODULUS)
C     PROPS(6): OPTIOSN FOR DIFFERENTIATION METHOD
C                 1- ANAYLTICAL DIFFERENTIATION
C                 2- NUMERICAL DIFFERENTIATION
C-----------------------------------------------------------------------
C     #   MATERIAL PARAMETERS
C     PROPS(7) : YOUNG'S MODULUS
C     PROPS(8) : POISSON'S RATIO
C-----------------------------------------------------------------------
C     #   HARDENING LAW
C     PROPS(9:18)=P(1:10): HARDENING LAW PARAMETERS
C-----------------------------------------------------------------------
C     #   ELASTIC MODULUS DEGRADATION
C     ##  PROPS(4)=1
C     PROPS(7) : EMOD_A
C     0        : EMOD_B
C     ##  PROPS(4)=2
C     PROPS(19): EMOD_A (FOR CHORD-E) EX. EMOD_A=175700.
C     PROPS(20): EMOD_B (FOR CHORD-E) EX. EMOD_B=49.6
C-----------------------------------------------------------------------
C     ##  PROPS(1)=1: CHABOCHE
C     #   PROPS(21:29): COEFFICIENTS FOR CHABOCHE
C     #   CHABOCHE1 [9]: TABLE6
C     PROPS(21):      H1: H1
C     PROPS(22):      H2: H2
C     PROPS(23):      H3: RATIO OF ISO- TO KINEMATIC HARDENING
C     #   CHABOCHE2 [10]
C     PROPS(21):      H1: A
C     PROPS(22):      H2: B
C     PROPS(23:29):   H3:H9
C
C     ##  PROPS(1)=2: YOSHIDA-UEMORI
C     #   PROPS(3)=0: ORIGINAL | >0: MODIFIED
C     #   PROPS(21:31): COEFFICIENTS FOR YOSHIDA-UEMORI MODEL [11]
C     PROPS(21):      H1: B
C     PROPS(22):      H2: Y
C     PROPS(23):      H3: C
C     PROPS(24):      H4: m
C     PROPS(25):      H5: b
C     PROPS(26):      H6: h
C     PROPS(27):      H7: YH1: R_sat
C     PROPS(28):      H8: YH2
C     PROPS(29):      H9: YH3
C     PROPS(30):      H10:YH4
C     PROPS(31):      H11:YH5
C
C     ##  PROPS(1)=3: HAH11 MODEL
C     #   PROPS(21:27): COEFFICIENTS USED FOR HAH11 MODEL [2]
C     PROPS(21):      H1: QM: EXPONENT OF HAH MODEL
C     PROPS(22):      H2: K
C     PROPS(23:27):   H3~H7: K1~K5
C
C     ##  PROPS(1)=4: HAH14 MODEL
C     #   PROPS(21:31): K-VALUES USED FOR HAH14 MODEL [3]
C     PROPS(21):      H1: EXPONENT OF HAH MODEL
C     PROPS(22):      H2: K
C     PROPS(23:27):   H3~H7: K1~K5
C     PROPS(28:31):   H8~H11: CROSS-LOADING / LATENT-HARDENING COEFFICIENTS
C
C     ##  PROPS(1)=5: HAH20 MODEL
C     #   PROPS(21:31): K-VALUES USED FOR HAH20 MODEL [14]
C     PROPS(21:22):   H1, H2: QM, PM: EXPONENT OF HAH MODEL
C     PROPS(23:27):   H3~H7: K1~K5
C     PROPS(28:31):   H8~H11: CROSS-LOADING / LATENT-HARDENING COEFFICIENTS
C
C     ##  PROPS(1)=6: HEXAH(2022) MODEL
C     #   PROPS(21:31): COEFFICIENTS USED FOR HEXAH MODEL [19]
C     PROPS(21:23):   H1~H3: XTB/XRB/XBB : BAUSCHINGER EFFECT
C     PROPS(24:25):   H4~H5: XTP/XRP     : PERMANENT SOFTENING
C     PROPS(26:28):   H6~H8: XTC/XRC/XG  : CROSS-LOADING CONTRACTION
C     PROPS(29:31):   H9~H11:XTL/XDEL/XEL: MACROSCOPIC LATENT HARDENING
C-----------------------------------------------------------------------
C     ##  PROPS(2)=2: GENERALIZED HILL
C     #   PROPS(32:37): COEFFICIENTS FOR HILL YIELD FUNCTION
C     PROPS(32):      YM: E (EXPONENT)
C     PROPS(33):      Y1: F
C     PROPS(34):      Y2: G
C     PROPS(35):      Y3: H
C     PROPS(36):      Y4: L
C     PROPS(37):      Y5: M
C     PROPS(38):      Y6: N
C
C     ##  PROPS(2)=3: YLD2000_2D MODEL [1]
C     #   PROPS(32:40): COEFFICIENTS FOR YLD2000-2D
C     PROPS(32):      YM: AM (EXPONENT)
C     PROPS(33:40):   Y1~Y8: A1~A8
C
C     ##  PROPS(2)=4: YLD2004_18P MODEL
C     #   PROPS(32:50): COEFFICIENTS FOR YLD2004-18P
C     PROPS(32):      YM: CM (EXPONENT)
C     PROPS(33:41):   Y1~Y9: C1~C9
C     PROPS(42:50):   Y10~Y18: C10~C18
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.1   VUMAT                                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE VUMAT(NBLOCK,NDIR,NSHR,NSTATEV,NFIELDV,NPROPS,LANNEAL,
C  READ ONLY (UNMODIFIABLE) VARIABLES -
     1 STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS,
     2 DENSITY, STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, 
     3 FIELDOLD, STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     4 TEMPNEW, STRETCHNEW, DEGRADNEW, FIELDNEW,
C  WRITE ONLY (MODIFIABLE) VARIABLES -
     5 STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW )
C
      INCLUDE 'VABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     1 CHARLENGTH(NBLOCK), STRAININC(NBLOCK,NDIR+NSHR),
     2 RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK),
     3 STRETCHOLD(NBLOCK,NDIR+NSHR),
     4 DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     5 FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR),
     6 STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),
     7 ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     8 STRETCHNEW(NBLOCK,NDIR+NSHR),
     8 DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
     9 FIELDNEW(NBLOCK,NFIELDV),
     1 STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     2 ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
C

      REAL*8, DIMENSION(:), ALLOCATABLE:: STRESS, DFDS, STAT_VAR
      REAL*8, DIMENSION(:), ALLOCATABLE:: DEV_H0, DEV_SIG, DSTRAN
      REAL*8, DIMENSION(:,:), ALLOCATABLE:: DDSDDE, DDFDDS
C
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KYLD1/ YM, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9
      COMMON /KYLD2/ Y10, Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18
      COMMON /KIHARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
      COMMON /KHEXAH/ MSS
C-----------------------------------------------------------------------
C     1.1     UMAT PARAMETERS
C
C     [1]     ANISOTROPIC HARDENING MODEL PARAMETERS
      HARD_PAR=PROPS(1)
C     [2]     YIELD FUNCTION PARAMETER
      YLD_PAR=PROPS(2)
C     [3]     ISOTORPIC HARDENING LAW PARAMETER
      FLOW_PAR=PROPS(3)
C     [4]     ELASTIC MODULUS DEGRADATION
      CHORD_PAR=PROPS(4)
C     [5]     STRESS UPDATE ALGORITHM PARAMETER
      SUA_PAR=PROPS(5)
C     [6]     DIFFERENTIATION METHOD
      GRAD_PAR=2.D0
C-----------------------------------------------------------------------
C     1.2     ARRAY DIMENSION
C     [1]     NDIM3: STRESS TENSOR DIMENSION
      NDIM2=NSHR
      IF(NDIM2 .EQ. 1) THEN ! PLANE STRESS
          NDIM1=NDIR-1
      ELSEIF(NDIM2 .EQ. 3) THEN ! 3-D STRESS
          NDIM1=NDIR
      END IF
      NDIM3=NDIM1+NDIM2
C     [2]     NDIM4: SIZE OF DEVIATORIC STRESS TENSOR
      IF(NDIM3 .EQ. 3) THEN ! PLANE STRESS
          NDIM4=NDIM3+1
      ELSEIF(NDIM3 .EQ. 6) THEN ! 3-D STRESS
          NDIM4=NDIM3
      ELSE
          WRITE(*,*) '# ERROR: CHECK THE INPUT FILE'
          STOP
      END IF
C     [3]     NDIM5: DIMENSION OF G-VALUES
      IF(HARD_PAR .EQ. 3.D0) THEN   ! HAH11
          NDIM5= 4
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN   ! HAH14
          NDIM5= 6
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
          NDIM5= 10
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH (HAH22)
          NDIM5= 3 ! FOR LOCAL VALUES
      ELSE
          NDIM5= 1
      END IF
C     [4]     NDIM6: DIMENSION FOR HARDENING STATE VARIABLES
      IF(HARD_PAR .EQ. 1.D0) THEN ! CHABOCHE
C     #   STAT_VAR(1:6) = ALPHA
          NDIM6= NDIM3
      ELSEIF(HARD_PAR .EQ. 2.D0) THEN ! YOSHIDA-UEMORI
C     #   STAT_VAR(1:6)  = ALPHA
C     #   STAT_VAR(7:12) = BETA
C     #   STAT_VAR(13:18)= Q
C     #   STAT_VAR(19)   = R
          NDIM6= 3*NDIM3+1
      ELSEIF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
C     #   STAT_VAR(1:4)= G
C     #   STAT_VAR(5:10)= DEV_H
          NDIM6= NDIM5+NDIM4
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
C     #   STAT_VAR(1:7)= G
C     #   STAT_VAR(8:13)= DEV_H
          NDIM6= NDIM5+NDIM4
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
C     #   STAT_VAR(1:10)= G
C     #   STAT_VAR(11:16)= DEV_H
C     #   STAT_VAR(17:22)= DEV_HP
C     #   STAT_VAR(23:28)= DEV_HS
          NDIM6= NDIM5+3*NDIM4
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH (HAH22)
C     #   STAT_VAR(1:3)= GP/GL/NSS
C     #   STAT_VAR(4:NDIM)= PLAS ! NDIM=3+NDIM3
C     #   STAT_VAR(NDIM+1:NDIM+NSS*(NDIM5+NDIM4))= G-/G+/GC/DEV_H
C     #   STAT_VAR(NDIM+NSS*(NDIM5+NDIM4+1):NDIM+NSS*(NDIM5+NDIM4+1))= EQPLASI
C     #   STAT_VAR(NDIM+NSS*(NDIM5+NDIM4+2):NDIM+NSS*(NDIM5+NDIM4+2))= COS_XI
          MSS=3                               ! MAXIMUM NUMBER OF SLIP SYSTEM
          NDIM6= 3+NDIM3                      ! GLOBAL VARIABLES
          NDIM6= NDIM6+(MSS)*(NDIM5+NDIM4+2)  ! LOCAL VARIABLES
      ELSE
          NDIM6= 1
      END IF
C     [5]     NDIM7: DIMENSION OF DISLOCATION HARDENING STATE VARIABLES
      IF(FLOW_PAR.EQ.7.D0) THEN ! RGBV (K. KITAYAMA ET AL. [16])
C     #   STAT_VAR(NDIM7+1)=RHO_F
C     #   STAT_VAR(NDIM7+2)=RHO_R1
C     #   STAT_VAR(NDIM7+3)=RHO_R2
C     #   STAT_VAR(NDIM7+4)=RHO_L
C     #   STAT_VAR(NDIM7+5)=RHO_F01
C     #   STAT_VAR(NDIM7+6)=RHO_F02
C     #   STAT_VAR(NDIM7+7)=RHO_RTOT
C     #   STAT_VAR(NDIM7+8)=RHO_TOT
C     #   STAT_VAR(NDIM7+9)=DRDE
          NDIM7=9
      ELSE
          NDIM7=0
      END IF
      NDIM6=NDIM6+NDIM7
C     [6]     ALLOCATE VARIABLE SIZE OF ARRAY
      ALLOCATE(STRESS(NDIM3),DFDS(NDIM3),DSTRAN(NDIM3),DEV_H0(NDIM4),
     1         DEV_SIG(NDIM4),STAT_VAR(NDIM6),DDSDDE(NDIM3,NDIM3),
     2         DDFDDS(NDIM3,NDIM3))
C-----------------------------------------------------------------------
C     1.3     READ USER MATERIAL PROPERTIES
C     1.3.1   MECHANICAL PROPERTIES
C     [1]     YOUNG'S MODULUS
      EMOD0=PROPS(7)
      EMOD=EMOD0
C     [2]     POISSON'S RATIO
      ENU=PROPS(8)
C
C     1.3.2   READ PLASTICITY MODEL COEFFICIENTS
C     [1]     ISOTROPIC HARDENING LAW PARAMETERS
      P1=PROPS(9)
      P2=PROPS(10)
      P3=PROPS(11)
      P4=PROPS(12)
      P5=PROPS(13)
      P6=PROPS(14)
      P7=PROPS(15)
      P8=PROPS(16)
      P9=PROPS(17)
      P10=PROPS(18)
C     [2]     ELASTIC MODULUS DEGRADATION
      IF(CHORD_PAR .EQ. 0.D0) THEN
          EMOD_A=EMOD0
          EMOD_B=0.D0
      ELSE
          EMOD_A=PROPS(19)
          EMOD_B=PROPS(20)
      END IF
C     [3]     ANISOTROPIC HARDENING COEFFICIENTS
      H1 =PROPS(21)
      H2 =PROPS(22)
      H3 =PROPS(23)
      H4 =PROPS(24)
      H5 =PROPS(25)
      H6 =PROPS(26)
      H7 =PROPS(27)
      H8 =PROPS(28)
      H9 =PROPS(29)
      H10=PROPS(30)
      H11=PROPS(31)
C     [4]     YIELD FUNCTION COEFFICIENT
      YM= PROPS(32)
      Y1= PROPS(33)
      Y2= PROPS(34)
      Y3= PROPS(35)
      Y4= PROPS(36)
      Y5= PROPS(37)
      Y6= PROPS(38)
      Y7= PROPS(39)
      Y8= PROPS(40)
      Y9= PROPS(41)
      Y10=PROPS(42)
      Y11=PROPS(43)
      Y12=PROPS(44)
      Y13=PROPS(45)
      Y14=PROPS(46)
      Y15=PROPS(47)
      Y16=PROPS(48)
      Y17=PROPS(49)
      Y18=PROPS(50)
C-----------------------------------------------------------------------
C     1.4     RECOVER STATE
C
C     [1]     ELASTIC TANGENT MODULUS
C     #   INITIALIZATION OF ELASTIC STIFFNESS MATRIX
C         ISOTROPIC ELASTIC STIFFNESS MATRIX FOR PLANE STRESS
      XMU= EMOD/(2.D0*(1.D0+ENU))                 ! #SHEAR MODULUS:   ¥ì=E/[2*(1+¥í)]
      XLM= ENU*EMOD/((1.D0+ENU)*(1.D0-2.D0*ENU))  ! #LAME'S CONSTANT: ¥ë=E¥í/[(1+¥í)(1-2¥í)]
      XVAL=EMOD/(1.D0-ENU**2)
C     [2]     TRIAL STRESS FOR RECOVER STATE
C     #   AT THE INITIAL STEP, VUMAT RETURNS TRIAL STRESS AS AN INITIAL STRESS STATE
      IF(STEPTIME .EQ. 0.) THEN
          DO K=1, NBLOCK
          TRACE= STRAININC(K,1)+STRAININC(K,2)+STRAININC(K,3)
          IF(NDIM3 .EQ. 3) THEN ! PLANE STRESS
              STRESSNEW(K,1)= STRESSOLD(K,1)
     1                       +XVAL*(STRAININC(K,1)+ENU*STRAININC(K,2))
              STRESSNEW(K,2)= STRESSOLD(K,2)
     1                       +XVAL*(STRAININC(K,2)+ENU*STRAININC(K,1))
              STRESSNEW(K,3)= 0.D0
              STRESSNEW(K,4)= STRESSOLD(K,4)+2.D0*XMU*STRAININC(K,4)
          ELSEIF(NDIM3 .EQ. 6) THEN ! 3D STRESS
              DO I=1, NDIM1
                  STRESSNEW(K,I)= STRESSOLD(K,I)
     1                       +2.D0*XMU*STRAININC(K,I)+XLM*TRACE
              END DO
              DO I=NDIM1+1,NDIM3
                  STRESSNEW(K,I)= STRESSOLD(K,I)+2.D0*XMU*STRAININC(K,I)
              END DO
          END IF
          END DO
          RETURN
      END IF
C-----------------------------------------------------------------------
      DO 100 K=1, NBLOCK
C     1.5     ELASTIC MATETRIAL PROPERTIES
      EMOD=STATEOLD(K,3)
      EMOD0=PROPS(7)
      IF(EMOD.EQ.0.D0) THEN
          EMOD=EMOD0
      END IF
      ENU=PROPS(8)
C-----------------------------------------------------------------------
C     1.6     ELASTIC TANGENT MODULUS
      DDSDDE=0.D0
      XMU= EMOD/(2.D0*(1.D0+ENU))                 ! #SHEAR MODULUS:   ¥ì=E/[2*(1+¥í)]
      XLM= ENU*EMOD/((1.D0+ENU)*(1.D0-2.D0*ENU))  ! #LAME'S CONSTANT: ¥ë=E¥í/[(1+¥í)(1-2¥í)]
      XVAL=EMOD/(1.D0-ENU**2)
      IF(NDIM3 .EQ. 3) THEN   !PLANE STRESS
C     #   DDSDDE(1:2, 1:2)
          DO I=1, NDIM1
          DO J=1, NDIM1
              DDSDDE(I,J)= XVAL*ENU
          END DO
              DDSDDE(I,I)= XVAL*1.D0
          END DO
C     #   DDSDDE(3,3)
          DDSDDE(NDIM3,NDIM3)= XMU
      ELSE                    !3D-STRESS
C     #   DDSDDE(1:3, 1:3)
          DO I=1, NDIM1
          DO J=1, NDIM1
              DDSDDE(I,J)=XLM
          END DO
          DDSDDE(I,I)=XLM+2.D0*XMU
          END DO
C     #   DDSDDE(4:6, 4:6)
          DO I=NDIM1+1, NDIM3
              DDSDDE(I,I)= XMU
          END DO
      END IF
C-----------------------------------------------------------------------
C     1.6     TOTAL STRAIN INCREMENT
      IF(NDIM3 .EQ. 3) THEN
          DSTRAN(1)=STRAININC(K,1)
          DSTRAN(2)=STRAININC(K,2)
          DSTRAN(3)=2.D0*STRAININC(K,4)
      ELSE
          DO I=1, NDIM1
              DSTRAN(I)=STRAININC(K,I)
          END DO
          DO I=NDIM1+1,NDIM3
              DSTRAN(I)=2.D0*STRAININC(K,I)
          END DO
      END IF
C-----------------------------------------------------------------------
C     1.7     TRIAL STRESS
C     #   THE UPDATED STRESS IS INITIALLY ASSUMED TO BE PURELY ELASTIC
C         FOR A GIVEN STRAIN INCREMENT, WHILE PLASTIC STATE VARIABLES
C         ARE KEPT CONSTANT [4].
C     C = CE
C     SIGT(N+1) = SIG(N) + C:DSTRAN(N+1)  [EQ.(23)][4]
      IF(NDIM3 .EQ. 3) THEN ! PLANE STRESS
        STRESS(1)= STRESSOLD(K,1)+XVAL*(DSTRAN(1)+ENU*DSTRAN(2))
        STRESS(2)= STRESSOLD(K,2)+XVAL*(DSTRAN(2)+ENU*DSTRAN(1))
        STRESS(3)= STRESSOLD(K,4)+XMU*DSTRAN(3)
      ELSEIF(NDIM3 .EQ. 6) THEN ! 3D STRESS
          TRACE= DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
          DO I=1, NDIM1
              STRESS(I)= STRESSOLD(K,I)
     1                       +2.D0*XMU*DSTRAN(I)+XLM*TRACE
          END DO
          DO I=NDIM1+1,NDIM3
              STRESS(I)= STRESSOLD(K,I)+XMU*DSTRAN(I)
          END DO
      END IF
C-----------------------------------------------------------------------
C     1.6     READ STATE VARIABLES
C     [1] EQUIVALENT PLASTIC STRAIN
      EQPLAS=STATEOLD(K,1)
C     [2] EFFECTIVE PLASTIC STRAIN INCREMENT
      D_EQPLAS=STATEOLD(K,6)
C     [3] STATE VARIABLES OF HARDENING MODELS
C     #   WHEN EQPLAS=0, THE STATE VARIABLES WOULD BE INITIALIZED
      DFDS=0.D0
      STAT_VAR=0.D0
      IF(EQPLAS .EQ. 0.D0) THEN
          IF(HARD_PAR .EQ. 1.D0) THEN ! CHABOCHE
              STAT_VAR= 0.D0
          ELSEIF(HARD_PAR .EQ. 2.D0) THEN ! YOSHIDA-UEMORI
              STAT_VAR= 0.D0
          ELSEIF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
C         #   G-VALUES
              STAT_VAR(1:NDIM5)=1.D0
C         #   DEV_H
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO I=1, NDIM4
                  STAT_VAR(NDIM5+I)=DEV_H0(I)
              END DO
          ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
C         #   G-VALUES
              STAT_VAR(1:NDIM5)=1.D0
C         #   DEV_H
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO I=1, NDIM4
                  STAT_VAR(NDIM5+I)=DEV_H0(I)
              END DO
          ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
C         #   G-VALUES
              STAT_VAR(1:NDIM5)= 1.D0
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO I=1, NDIM4
C         #   DEV_H
                  STAT_VAR(NDIM5+I)=DEV_H0(I)
C         #   DEV_HP
                  STAT_VAR(NDIM5+NDIM4+I)=DEV_H0(I)
C         #   DEV_HS
                  STAT_VAR(NDIM5+2*NDIM4+I)=DEV_H0(I)
              END DO
          ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
C         [1] GLOBAL STATE VARIABLES
              NDIM0=7
C         #   GLOBAL G-VALUES
              STAT_VAR(1:2)=1.D0
C         #   NSS: NUMBER OF ACTIVATED SLIP SYSTEMS
              NSS=1
              STAT_VAR(3)= NSS
C         #   PLASTIC STRAIN: PLAS
              DO I=1,NDIM3
                  STAT_VAR(3+I)=0.D0
              END DO
C         [2] LOCAL STATE VARIABLES
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO N=1, NSS
                  NDIM =3+NDIM3+(N-1)*(NDIM5+NDIM4+2)
                  NDIM0=7+3+6+(N-1)*(NDIM5+6+2) ! 16+X
C         #   LOCAL G-VALUES
                  DO I=1,NDIM5
                      STAT_VAR(NDIM+I)=1.D0
                  END DO
C         #   LOCAL MICROSTRUCTURE DEVIATOR
                  DO I=1,NDIM4
                      STAT_VAR(NDIM+NDIM5+I)=DEV_H0(I)
                  END DO
C         #   LOCAL EQUIVALENT PLASTIC STRAIN
                  STAT_VAR(NDIM+NDIM5+NDIM4+1)=0.D0
C         #   LOCAL STRAIN PATH CHANGE PARAMETER
                  STAT_VAR(NDIM+NDIM5+NDIM4+2)=1.D0
              END DO
          END IF
      ELSE
          IF(HARD_PAR .EQ. 1.D0) THEN     ! CHABOCHE
              NDIM0=6
              DO I=1, NDIM3
                  STAT_VAR(I)=STATEOLD(K,NDIM0+I)
              END DO
          ELSEIF(HARD_PAR .EQ. 2.D0) THEN ! YOSHIDA-UEMORI
              NDIM0=6
              DO I=1, NDIM3
                  STAT_VAR(I)=STATEOLD(K,NDIM0+I)
                  STAT_VAR(NDIM3+I)=STATEOLD(K,NDIM0+6+I)
                  STAT_VAR(2*NDIM3+I)=STATEOLD(K,NDIM0+12+I)
              END DO
              STAT_VAR(3*NDIM3+1)=STATEOLD(K,NDIM0+18+1)
          ELSEIF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
              NDIM0=7
C         #   G-VALUES
              STAT_VAR(1:NDIM5)=STATEOLD(K,NDIM0+1:NDIM0+NDIM5)
C         #   DEV_H
              DO I=1, NDIM4
                  STAT_VAR(NDIM5+I)=STATEOLD(K,NDIM0+NDIM5+I)
              END DO
          ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
              NDIM0=7
C         #   G-VALUES
              STAT_VAR(1:NDIM5)=STATEOLD(K,NDIM0+1:NDIM0+NDIM5)
C         #   DEV_H
              DO I=1, NDIM4
                  STAT_VAR(NDIM5+I)=STATEOLD(K,NDIM0+NDIM5+I)
              END DO
          ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
              NDIM0=7
C         #   G-VALUES
              STAT_VAR(1:NDIM5)= STATEOLD(K,NDIM0+1:NDIM0+NDIM5)
              DO I=1, NDIM4
C         #   DEV_H
                  STAT_VAR(NDIM5+I)=STATEOLD(K,NDIM0+NDIM5+I)
C         #   DEV_HP
                  STAT_VAR(NDIM5+NDIM4+I)=STATEOLD(K,NDIM0+NDIM5+6+I)
C         #   DEV_HS
                  STAT_VAR(NDIM5+2*NDIM4+I)=STATEOLD(K,NDIM0+NDIM5+12+I)
              END DO
          ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
C         [1] GLOBAL STATE VARIABLES
              NDIM0=7
C         #   GLOBAL G-VALUES
              STAT_VAR(1:2)=STATEOLD(K,NDIM0+1:NDIM0+2)
C         #   NSS: NUMBER OF ACTIVATED SLIP SYSTEMS
              STAT_VAR(3)= STATEOLD(K,NDIM0+3)
              NSS=IDINT(STAT_VAR(3))
C         #   PLASTIC STRAIN: PLAS
              DO I=1,NDIM3
                  STAT_VAR(3+I)=STATEOLD(K,NDIM0+3+I)
              END DO
C         [2] LOCAL STATE VARIABLES
              DO N=1, NSS
                  NDIM =3+NDIM3+(N-1)*(NDIM5+NDIM4+2)
                  NDIM0=7+3+6+(N-1)*(NDIM5+6+2) ! 16+X
C         #   LOCAL G-VALUES
                  DO I=1,NDIM5
                      STAT_VAR(NDIM+I)=STATEOLD(K,NDIM0+I)
                  END DO
C         #   LOCAL MICROSTRUCTURE DEVIATOR
                  DO I=1,NDIM4
                      STAT_VAR(NDIM+NDIM5+I)=STATEOLD(K,NDIM0+NDIM5+I)
                  END DO
C         #   LOCAL EQUIVALENT PLASTIC STRAIN
                  STAT_VAR(NDIM+NDIM5+NDIM4+1)=STATEOLD(K,NDIM0+NDIM5+7)
C         #   LOCAL STRAIN PATH CHANGE PARAMETER
                  STAT_VAR(NDIM+NDIM5+NDIM4+2)=STATEOLD(K,NDIM0+NDIM5+8)
              END DO
          END IF
      END IF
C     [4] READ STATE VARIABLES OF DISLOCATION-BASED HARDENING MODEL
      IF(EQPLAS .EQ. 0.D0) THEN
          IF(FLOW_PAR .EQ. 7.D0) THEN ! RGBV
              RHO0=P10
              STAT_VAR(NDIM6-NDIM7+1)=RHO0 !RHO_F
              STAT_VAR(NDIM6-NDIM7+2)=0.D0 !RHO_R1
              STAT_VAR(NDIM6-NDIM7+3)=0.D0 !RHO_R2
              STAT_VAR(NDIM6-NDIM7+4)=0.D0 !RHO_L
              STAT_VAR(NDIM6-NDIM7+5)=RHO0 !RHO_F01
              STAT_VAR(NDIM6-NDIM7+6)=RHO0 !RHO_F02
              STAT_VAR(NDIM6-NDIM7+7)=RHO0 !RHO_TOT
              STAT_VAR(NDIM6-NDIM7+8)=0.D0 !DRDE
              STAT_VAR(NDIM6-NDIM7+9)=0.D0 !RHO_RTOT (TEMPORARY)
          END IF
      ELSE
          IF(FLOW_PAR .EQ. 7.D0) THEN ! RGBV
              DO I=1, NDIM7
                  STAT_VAR(NDIM6-NDIM7+I)=STATEOLD(K,49+I)
              END DO
          END IF
      END IF
C-----------------------------------------------------------------------
C     1.7     YIELD CONDITION
C     [1]     STRESS=STRESS-BACKSTRESS
      IF(HARD_PAR.EQ.1.D0 .OR. HARD_PAR.EQ.2.D0) THEN
          DO I=1, NDIM3
              STRESS(I)= STRESS(I)-STAT_VAR(I)
          END DO
      END IF
C     [2]     RESIDUAL FUNCTION VALUE
      CALL YIELD_CONDITION(STAT_VAR, EQPLAS, STRESS, FVAL)
C-----------------------------------------------------------------------
C     1.8     STRESS INTEGRATION METHODS
      ITER=0
      STATE=0.D0
C     [1]     PLASTIC CORRECTOR
      IF(FVAL .GE. 1.D-6) THEN
          STATE=1.D0
          IF(SUA_PAR .GE. 0.D0) THEN
              CALL STRESS_UPDATE  ! ITERATIVE STRESS UPDATE METHOD
     1                   (STAT_VAR,EQPLAS,STRESS,DDSDDE,ITER)
          ELSE
              CALL STRESS_UPDATE2 ! NON-ITERATIVE STRESS UPDATE METHOD
     1                   (STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,STRESS,DDSDDE)
          END IF
      END IF
C     [2]     EFFFECTIVE STRESS
      CALL DOUBLE_DOT(STRESS, STRESS, NDIM3, XSIG)
      CALL EFF_SIG(STAT_VAR, STRESS, SIG_BAR)
C     [3]     GRADIENT OF YIELD SURFACE
      CALL GRAD(2,1,STAT_VAR,STRESS,DFDS,DDFDDS)
C     [4]     STRESS=STRESS+BACKSTRESS
      IF(HARD_PAR.EQ.1.D0 .OR. HARD_PAR.EQ.2.D0) THEN
          DO I=1, NDIM3
              STRESS(I)=STRESS(I)+STAT_VAR(I)
          END DO
      END IF
C     RETURN:
C             (1) STRESS  : UPDATED STRESS TENSOR
C             (2) EQPLAS  : EQUIVALENT PLASTIC STRAIN
C             (3) STAT_VAR: STATE VARIABLES
C             (3) DFDS    : YIELD SURFACE GRADIENT
C-----------------------------------------------------------------------
C     1.9     UPDATE STRESS AND STRAIN
C
C     [1]     EFFECTIVE STRAIN INCREMENT
      D_EQPLAS=EQPLAS-STATEOLD(K,1)
C     [2]     STRAIN INCREMENT 33 IN PLANE STRESS
      IF(NDIM3 .EQ. 3) THEN
          D_PLAS1=D_EQPLAS*DFDS(1)
          D_PLAS2=D_EQPLAS*DFDS(2)
          DEL_SIG1=STRESS(1)-STRESSOLD(K,1)
          DEL_SIG2=STRESS(2)-STRESSOLD(K,2)
          D_PLAS3=-(D_PLAS1+D_PLAS2)
          STRAININC(K,3)=-(ENU/EMOD)*(DEL_SIG1+DEL_SIG2)+D_PLAS3
      END IF
C     [3]     UPDATED STRESS TENSOR
      IF(NDIM3 .EQ. 3) THEN ! PLANE STRESS
          STRESSNEW(K,1)=STRESS(1)
          STRESSNEW(K,2)=STRESS(2)
          STRESSNEW(K,3)=0.D0
          STRESSNEW(K,4)=STRESS(3)
      ELSEIF(NDIM3 .EQ. 6) THEN
          STRESSNEW(K,:)=STRESS(:)
      END IF
C-----------------------------------------------------------------------
C     1.10    UPDATE SATATE VARIABLES
C     [1]     EFFECTIVE STRAIN
      STATENEW(K,1)=EQPLAS
C     [2]     EFFECTIVE STRESS
      STATENEW(K,2)=SIG_BAR
C     [3]     YOUNG'S MODULUS
      STATENEW(K,3)=EMOD0-(EMOD0-EMOD_A)*(1.D0-EXP(-EMOD_B*EQPLAS))
C     [4]     DEFORMATION STATE
      STATENEW(K,4)=STATE
C     [5]     ALGORITHMIC CONVERGENCE ITERATION
      STATENEW(K,5)=ITER
C     [6]     EFFECTIVE PLASTIC STRAIN INCREMENT
      STATENEW(K,6)=D_EQPLAS
C     [7]     STATE VARIABLES OF ANISOTROPIC HARDENING MODELS
      IF(HARD_PAR .EQ. 1.D0) THEN ! CHABOCHE
          NDIM0=6
C         #STATEV(7:12)    : ALPHA
          DO I=1, NDIM3
              STATENEW(K,NDIM0+I)= STAT_VAR(I)
          END DO
      ELSEIF(HARD_PAR .EQ. 2.D0) THEN ! YOSHIDA-UEMORI
          NDIM0=6
C         #STATEV(7:12)    : ALPHA
C         #STATEV(13:18)   : BETA
C         #STATEV(19:24)   : Q
C         #STATEV(25)      : R
C         #STATEV(26)      : STAG_R
          DO I=1, NDIM3
              STATENEW(K,NDIM0+I)=STAT_VAR(I)
              STATENEW(K,NDIM0+6+I)=STAT_VAR(NDIM3+I)
              STATENEW(K,NDIM0+12+I)=STAT_VAR(2*NDIM3+I)
          END DO
          STATENEW(K,NDIM0+18+1)=STAT_VAR(3*NDIM3+1)
      ELSEIF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
          NDIM0=7
C         #STATEV(7)      : COS_X
          CALL COSX(STAT_VAR,STRESS,COS_X)
          STATENEW(K,NDIM0)= COS_X
C         #STATEV(8:11)    : G-VALUES
          DO I=1, NDIM5
              STATENEW(K,NDIM0+I)=STAT_VAR(I)
          END DO
C         #STATEV(12:17)  : DEV_H
          DO I=1, NDIM4
              STATENEW(K,NDIM0+NDIM5+I)=STAT_VAR(NDIM5+I)
          END DO
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
          NDIM0=7
C         #STATEV(7)      : COS_X
          CALL COSX(STAT_VAR,STRESS,COS_X)
          STATENEW(K,NDIM0)= COS_X
C         #STATEV(8:13)    : G-VALUES
          DO I=1, NDIM5
              STATENEW(K,NDIM0+I)=STAT_VAR(I)
          END DO
C         #STATEV(14:19)  : DEV_H
          DO I=1, NDIM4
              STATENEW(K,NDIM0+NDIM5+I)=STAT_VAR(NDIM5+I)
          END DO
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
          NDIM0=7
C         #STATEV(7)      : COS_X
          CALL COSX(STAT_VAR,STRESS,COS_X)
          STATENEW(K,NDIM0)= COS_X
C         #STATEV(8:17)   : G-VALUES
          DO I=1, NDIM5
              STATENEW(K,NDIM0+I)=STAT_VAR(I)
          END DO
C         #STATEV(18:23)  : DEV_H
C         #STATEV(24:29)  : DEV_HP
C         #STATEV(30:35)  : DEV_HS
          DO I=1, NDIM4
              STATENEW(K,17+I)=STAT_VAR(NDIM5+I)
              STATENEW(K,23+I)=STAT_VAR(NDIM5+NDIM4+I)
              STATENEW(K,29+I)=STAT_VAR(NDIM5+2*NDIM4+I)
          END DO
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
          NDIM0=7
C         #STATEV(7)      : COS_X
          CALL COSX(STAT_VAR,STRESS,COS_X)
          STATENEW(K,NDIM0)= COS_X
C         [1] GLOBAL STATE VARIABLES
C         #   GLOBAL G-VALUES: GP/GL
          STATENEW(K,NDIM0+1:NDIM0+2)=STAT_VAR(1:2)
C         #   NSS: NUMBER OF ACTIVATED SLIP SYSTEMS
          STATENEW(K,NDIM0+3)=STAT_VAR(3)
          NSS=IDINT(STAT_VAR(3))
C         #   PLASTIC STRAIN: PLAS=PLAS+D_EQPLAS*DFDS
          DO I=1,NDIM3
              STATENEW(K,NDIM0+3+I)=STAT_VAR(3+I)
          END DO
C         [2] LOCAL STATE VARIABLES
          DO N=1, NSS
              NDIM =NDIM3+3+(N-1)*(NDIM5+NDIM4+2)
              NDIM0=7+3+6+(N-1)*(NDIM5+6+2)
C         #   LOCAL G-VALUES
              DO I=1,NDIM5
                  STATENEW(K,NDIM0+I)=STAT_VAR(NDIM+I)
              END DO
C         #   LOCAL MICROSTRUCTURE DEVIATOR
              DO I=1,NDIM4
                  STATENEW(K,NDIM0+NDIM5+I)=STAT_VAR(NDIM+NDIM5+I)
              END DO
C         #   LOCAL EQUIVALENT PLASTIC STRAIN
              STATENEW(K,NDIM0+NDIM5+6+1)=STAT_VAR(NDIM+NDIM5+NDIM4+1)
C         #   LOCAL STRAIN PATH CHANGE PARAMETER
              STATENEW(K,NDIM0+NDIM5+6+2)=STAT_VAR(NDIM+NDIM5+NDIM4+2)
          END DO
      ELSEIF(HARD_PAR .EQ. 0.D0) THEN ! ISOTROPIC HARDENING
          STATENEW(K,7:49)=0.D0
      END IF
C     [8]     STATE VARIALBES OF DISLOACTION HARDENING MODELS
      IF(FLOW_PAR .EQ. 7.D0) THEN
C         # STATE VARIABLE OF RGBV #
C         STATEV(50): RHO_F
C         STATEV(51): RHO_R1
C         STATEV(52): RHO_R2
C         STATEV(53): RHO_L
C         STATEV(54): RHO_F01
C         STATEV(55): RHO_F02
C         STATEV(56): RHO_TOT
C         STATEV(57): DRDE
C         STATEV(58): RHO_RTOT
          DO I=1,NDIM7
              STATENEW(K,50+I)=STAT_VAR(NDIM6-NDIM7+I)
          END DO
      END IF
C     [9]     GRADIENT OF YIELD SURFACE
      DO I=1, NDIM3
          STATENEW(K,58+I)=DFDS(I)
      END DO
C-----------------------------------------------------------------------
C     1.11    UPDATE INTERNAL ENERGY PER UNIT MASS
      DO I=1, NDIM1
          STRESSPOWER=STRESSPOWER+(STRESSOLD(K,I)+
     1                STRESSNEW(K,I))*STRAININC(K,I)
      END DO
      DO I=NDIM1+1,NDIM1+NDIM2
          STRESSPOWER=STRESSPOWER+2.D0*(STRESSOLD(K,I)+
     1                STRESSNEW(K,I))*STRAININC(K,I)
      END DO
      ENERINTERNNEW(K)=ENERINTERNOLD(K)+ 0.5*STRESSPOWER/DENSITY(K)

C     1.4.11  UPDATE THE DISSIPATED INELASTIC SPECIFIC ENERGY
      PWORKINC=0.5*(STATENEW(K,2)+STATEOLD(K,2))*D_EQPLAS
      ENERINELASNEW(K)=ENERINELASOLD(K)+PWORKINC/DENSITY(K)

100   CONTINUE ! END OF BLOCK
C-----------------------------------------------------------------------
      DEALLOCATE(STRESS,DFDS,DSTRAN,STAT_VAR,DEV_H0,DEV_SIG,
     1           DDSDDE,DDFDDS)
      RETURN
      END SUBROUTINE VUMAT
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.2   STRESS UPDATE ALGORITHMS                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE STRESS_UPDATE(STAT_VAR,EQPLAS,SIG,CE,ITER)
C     # THE SUBROUTINES NUMERICALLY INTEGRATE THE STRESS TENSOR THROUGH
C       THE INDIRECT METHOD.
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION STAT_VAR(NDIM6),SIG(NDIM3),CE(NDIM3,NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KITER/ MAX_ITER

C     [1] MAXIMUM ITERATION
      MAX_ITER=2000
C     [2] SELECTION OF STRESS UPDATE ALGORITHM
      IF(FLOOR(SUA_PAR) .EQ. 1.D0) THEN
          CALL CPM(STAT_VAR,EQPLAS,SIG,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 2.D0) THEN
          CALL CPPM(STAT_VAR,EQPLAS,SIG,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 3.D0) THEN
          CALL CPPM2(STAT_VAR,EQPLAS,SIG,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 4.D0) THEN
          CALL TMM(STAT_VAR,EQPLAS,SIG,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 5.D0) THEN
C          CALL BROYDEN(STAT_VAR,EQPLAS,SIG,CE,ITER)
      END IF

      RETURN
      END SUBROUTINE STRESS_UPDATE
C-----------------------------------------------------------------------
      SUBROUTINE STRESS_UPDATE2(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
C     # THE SUBROUTINES NUMERICALLY INTEGRATE THE STRESS TENSOR THROUGH
C       THE INDIRECT METHOD.
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION STAT_VAR(NDIM6),DSTRAN(NDIM3),SIG(NDIM3),CE(NDIM3,NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER

C     [1] ORIGINAL STRESS TENSOR AT THE PREVIOUS DEFORMATION STEP
      SIG=SIG-MATMUL(CE,DSTRAN)
C-----------------------------------------------------------------------
C     [2] NON-ITERATIVE STRESS INTEGRATION METHOD
      IF(SUA_PAR .EQ. -1.D0) THEN
          CALL FOEFM(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
      ELSEIF(SUA_PAR .EQ. -2.D0) THEN
          CALL NSPM(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
      ELSEIF(SUA_PAR .EQ. -3.D0) THEN
          CALL NSPM2(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
      END IF

      RETURN
      END SUBROUTINE STRESS_UPDATE2
C-----------------------------------------------------------------------
C
      SUBROUTINE CPM(STAT_VAR,EQPLAS,SIG,CE,ITER)
C     EULER FORWARD METHOD= CUTTING PLANE METHOD (EXPLICIT SCHEME)#
C     REF: [4] J.W. LEE ET AL. IJP  (2012): CUTTING PLANE METHOD  #
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), CE(NDIM3,NDIM3)
      DIMENSION SIG_K(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3),
     1          DEL_SIG(NDIM3), DDFDEDS(NDIM3), CEP(NDIM3,NDIM3),
     2          R1(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      COMMON /KITER/ MAX_ITER
      PARAMETER(TOL=1.D-6)
C-----------------------------------------------------------------------
C     1.0     INITIALIZATION
C     ¥Å_BAR(0, N+1) = ¥Å_BAR(N)
      EQPLAS_K = EQPLAS
C     SIG(0, N+1) = SIG_TRIAL(N+1)
      SIG_K = SIG
      DEL_EQPLAS=0.D0
      DEL_SIG=0.D0
      FLAG=0.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0)        ! STEP-SIZE CONTROL PARAMETER
C-----------------------------------------------------------------------
210   CONTINUE
C     1.1     GRADIENTS OF YIELD SURFACE: DFDS
      CALL GRAD(2,1,STAT_VAR,SIG_K,DFDS,DDFDDS)
C     1.2     EFFECTIVE STRESS & FLOW STRESS
C     [1] EFFECTIVE STRESS
      CALL EFF_SIG(STAT_VAR, SIG_K, SIG_BAR_K)
C     [2] FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, EQPLAS_K, FLOW_SIG_K)
C     1.3    RESIDUAL FUNCTION: YIELD CONDITION
      R2=SIG_BAR_K-FLOW_SIG_K
      IF(ITER.EQ.0) R0=R2
C-----------------------------------------------------------------------
C     1.4     CONVERGENCE CHECK
      IF(DABS(R2) .LT. TOL) THEN
          GOTO 215
      END IF
      IF(ITER .GT. MAX_ITER) THEN
          WRITE(*,*) '#DIVERGENCE', R2, ITER
          FLAG=1.D0
      END IF
      IF(ISNAN(R2)) THEN
          WRITE(*,*) '#ERROR: SIG_BAR_K=NaN'
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          WRITE(*,*) 'INITIAL RESIDUAL', R0
          STOP
      END IF
C-----------------------------------------------------------------------
C     1.5     DENOMIATOR OF TANGENT STIFFNESS TENSOR
C     [1]     DFDS_C_DFDS
      DFDS_CE_DFDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_CE_DFDS= DFDS_CE_DFDS+DFDS(I)*CE(I,J)*DFDS(J)
      END DO
      END DO
C     [2]     KINEMATIC COMPONENTS: DFDE & DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS_K,SIG_K,DFDE,DDFDEDS)
C     [3]     DENOMINATOR
      XDENOM=DFDS_CE_DFDS-DFDE
C-----------------------------------------------------------------------
C     1.6     STEP-SIZE CONTROL
      R1=0.D0
      VAL=R2/XDENOM
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     1.7     STEP SIZE OF SOLUTION VARIABLE INCREMENTS
C     [1]     EQUIVALENT PLASTIC STRAIN
      DEL_EQPLAS=ALPHA*VAL
C     [2]     STRESS
      DEL_SIG= MATMUL(-CE,DFDS*DEL_EQPLAS)
C     1.8     INTEGRATE SOLUTION VARIABLES
C     [1]     EQUIVALENT PLASTIC STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS
C     [2]     STRESS TENSOR
      SIG_K=SIG_K+DEL_SIG
C     1.9     UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS_K, SIG_K, DEL_EQPLAS)
      END IF
      ITER=ITER+1
      GOTO 210
C-----------------------------------------------------------------------
215   CONTINUE
C     1.10    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     1.11    UPDATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
      END IF
C     1.12    CONTINUUM TANGENT MODULUS
      CALL CONTINUUM(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
      CE=CEP

      RETURN
      END SUBROUTINE CPM
C-----------------------------------------------------------------------
C
C     EULER BACKWARD METHOD= CLOSETST POINT PROJECTION METHOD (FULLY IMPLICIT)
      SUBROUTINE CPPM(STAT_VAR,EQPLAS,SIG,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), CE(NDIM3,NDIM3)
      DIMENSION SIG_K(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3),
     1          R1(NDIM3), D_PLAS(NDIM3), SE(NDIM3,NDIM3), 
     2          DR1DS(NDIM3,NDIM3), DR1DE(NDIM3), DDFDEDS(NDIM3),
     4          DEL_SIG(NDIM3), D_SIG(NDIM3), CEP(NDIM3,NDIM3)
      DIMENSION XJACOB(NDIM3+1,NDIM3+1), XJACOBI(NDIM3+1,NDIM3+1),
     1          XLHS(NDIM3+1),XRHS(NDIM3+1)
      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      COMMON /KITER/ MAX_ITER
C-----------------------------------------------------------------------
C     2.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_PLAS      = 0.D0
      DEL_EQPLAS  = 0.D0
      DEL_SIG     = 0.D0
      FLAG        = 0.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0)        ! STEP-SIZE CONTROL PARAMETER
C-----------------------------------------------------------------------
220   CONTINUE
C     2.1     GRADIENTS OF YIELD SURFACE: DFDS
      CALL GRAD(2,2,STAT_VAR,SIG_K,DFDS,DDFDDS)
C     2.2     EFFECTIVE STRESS & FLOW STRESS
C     [1] EFFECTIVE STRESS
      CALL EFF_SIG(STAT_VAR, SIG_K, SIG_BAR_K)
C     [2] FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, EQPLAS_K, FLOW_SIG_K)
C     2.3     RESIDUAL FUNCTIONS
C     [1]     ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_EQPLAS*DFDS(I) !(1,3)
          R1_NORM=R1_NORM+R1(I)**2
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     YIELD CONDITION
      R2=SIG_BAR_K-FLOW_SIG_K
      IF(ITER.EQ.0) R0=R2
C-----------------------------------------------------------------------
C     2.4    CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 225
      END IF
      IF(ITER .GT. MAX_ITER) THEN
          WRITE(*,*) '#DIVERGENCE', R1_NORM, DABS(R2)
          WRITE(*,*) '#TOLERANCE', TOL, TOL
          FLAG=1.D0
      END IF
      IF(ISNAN(R2)) THEN
          WRITE(*,*) '#ERROR: SIG_BAR_K=NaN'
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          WRITE(*,*) 'INITIAL RESIDUAL', R0
          STOP
      END IF
C-----------------------------------------------------------------------
C     2.5     GRADIENT OF RESIDUAL FUNCTION I
C     [1]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, SE)
C     [2]     DR1DS: INVERSE HESSIAN MATRIX
      DO I=1, NDIM3
      DO J=1, NDIM3
          DR1DS(I,J)= SE(I,J)+D_EQPLAS*DDFDDS(I,J)
      END DO
      END DO
C     [3]     KINEMATIC COMPONENTS: DFDE & DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS_K,SIG_K,DFDE,DDFDEDS)
C     [4]     DR1DE: XK
      DR1DE=D_EQPLAS*DDFDEDS+DFDS
C-----------------------------------------------------------------------
C     2.6     JACOBIAN MATRIX
C     NOTE: DR2DS=DFDS & DR2DE=DFDE
      XJACOB(1:NDIM3,1:NDIM3)=DR1DS
      XJACOB(NDIM3+1,NDIM3+1)=DFDE
      XJACOB(1:NDIM3,NDIM3+1)=DR1DE
      XJACOB(NDIM3+1,1:NDIM3)=DFDS
C     #       INVERSE JACOBIAN MATRIX
      CALL INVERSE(-XJACOB, NDIM3+1, XJACOBI)
C     2.7     SOLUTION VARIABLES: LHS & RHS
C     [1] LEFT-HAND SIDE: LHS
      DO I=1, NDIM3+1
          IF(I .LE. NDIM3) THEN
              XLHS(I)=R1(I)
          ELSE
              XLHS(I)=R2
          END IF
      END DO
C     [2] RIGHT-HAND SIDE: RHS
      XRHS=MATMUL(XJACOBI,XLHS)
C-----------------------------------------------------------------------
C     2.8     STEP-SIZE CONTROL
      VAL=XRHS(NDIM3+1)
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     2.9     STEP-SIZE OF SOLUTION VARIABLE INCREMENTS
C     [1]     EQUIVLANET PLASTIC STRAIN
      DEL_EQPLAS=ALPHA*VAL
C     [2]     STRES
      DEL_SIG(1:NDIM3)=ALPHA*XRHS(1:NDIM3)
C     2.10    INCREMENT AT A TIME-STEP
C     [1]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [2]     EQUIVALENT PALSTIC STRAIN INCREMENT AT A TIME-STEP
      D_EQPLAS=D_EQPLAS+DEL_EQPLAS
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-SE,D_SIG)
C     2.11    INTEGRATE SOLUTION VARIALBES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS
C     2.12    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS_K, SIG_K, DEL_EQPLAS)
      END IF
      ITER=ITER+1
      GOTO 220
C-----------------------------------------------------------------------
225   CONTINUE
C     2.13    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     2.14    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
      END IF
C     2.15    CONSISTENT TANGENT MODULUS
      CALL CONSISTENT(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
      CE=CEP

      RETURN
      END SUBROUTINE CPPM
C-----------------------------------------------------------------------
C
C     EULER BACKWARD METHOD= CLOSETST POINT PROJECTION METHOD (SEMI-IMPLICIT)
      SUBROUTINE CPPM2(STAT_VAR,EQPLAS,SIG,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), CE(NDIM3,NDIM3)
      DIMENSION SIG_K(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3),
     1          R1(NDIM3), D_PLAS(NDIM3), SE(NDIM3,NDIM3), 
     2          DR1DS(NDIM3,NDIM3), DR1DE(NDIM3), DDFDEDS(NDIM3),
     4          DEL_SIG(NDIM3), D_SIG(NDIM3), CEP(NDIM3,NDIM3)
      DIMENSION XHI(NDIM3,NDIM3), XH(NDIM3,NDIM3)
      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      COMMON /KITER/ MAX_ITER
C-----------------------------------------------------------------------
C     3.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_PLAS      = 0.D0
      DEL_EQPLAS  = 0.D0
      DEL_SIG     = 0.D0
      FLAG        = 0.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C-----------------------------------------------------------------------
230   CONTINUE
C     3.1     GRADIENTS OF YIELD SURFACE: DFDS
      CALL GRAD(2,2,STAT_VAR,SIG_K,DFDS,DDFDDS)
C     3.2     EFFECTIVE STRESS & FLOW STRESS
C     [1] EFFECTIVE STRESS
      CALL EFF_SIG(STAT_VAR, SIG_K, SIG_BAR_K)
C     [2] FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, EQPLAS_K, FLOW_SIG_K)
C     3.3     RESIDUAL FUNCTIONS
C     [1]     ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_EQPLAS*DFDS(I) !(1,3)
          R1_NORM=R1_NORM+R1(I)**2
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     YIELD CONDITION
      R2=SIG_BAR_K-FLOW_SIG_K
      IF(ITER.EQ.0) R0=R2
C-----------------------------------------------------------------------
C     3.4     CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 235
      END IF
      IF(ITER .GT. MAX_ITER) THEN
          WRITE(*,*) '#DIVERGENCE', R1_NORM, DABS(R2)
          WRITE(*,*) '#TOLERANCE', TOL, TOL
          FLAG=1.D0
      END IF
      IF(ISNAN(R2)) THEN
          WRITE(*,*) '#ERROR: SIG_BAR_K=NaN'
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          WRITE(*,*) 'INITIAL RESIDUAL', R0
          STOP
      END IF
C-----------------------------------------------------------------------
C     3.5     GRADIENT OF RESIDUAL FUNCTION:DR1DS / DR1DE
C     [1]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, SE)
C     [2]     DR1DS=XHI: INVERSE HESSIAN MATRIX
      DO I=1, NDIM3
      DO J=1, NDIM3
          XHI(I,J)= SE(I,J)+D_EQPLAS*DDFDDS(I,J)
      END DO
      END DO
C     [3]     KINEMATIC COMPONENTS: DFDE & DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS_K,SIG_K,DFDE,DDFDEDS)
C     [4]     DR1DE: XK
      DR1DE=D_EQPLAS*DDFDEDS+DFDS
C-----------------------------------------------------------------------
C     3.6     DENOMINATOR AND NUMERATOR
C     [1]     HESSIAN MATRIX: XH=DR1DSI
      CALL INVERSE(XHI, NDIM3, XH)
C     [2]     DFDS_XH_DR1DE
      DFDS_XH_DR1DE=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_XH_DR1DE=DFDS_XH_DR1DE+DFDS(I)*XH(I,J)*DR1DE(J)
      END DO
      END DO
C     [3]     XDENOM
      XDENOM=DFDS_XH_DR1DE-DFDE
C     [4]     DFDS_XH_R1
      DFDS_XH_R1=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_XH_R1= DFDS_XH_R1 + DFDS(I)*XH(I,J)*R1(J)
      END DO
      END DO
C     [5]     XNUMER
      XNUMER= R2 - DFDS_XH_R1
C-----------------------------------------------------------------------
C     3.7     STEP-SIZE CONTROL
      VAL=XNUMER/XDENOM
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     3.8     STEP-SIZE OF SOLUTION VARIABLE INCREMENTS AT A SUB-STEP
C     [1]     EQUIVALENT PLASTIC STRAIN
      DEL_EQPLAS= ALPHA*VAL
C     [2]     STRESS
      DEL_SIG= MATMUL(-XH,R1+DR1DE*DEL_EQPLAS)
C     3.9     INCREMENTS AT A TIME-STEP
C     [1]     EQUIVALENT PALSTIC STRAIN INCREMENT AT A TIME-STEP
      D_EQPLAS=D_EQPLAS+DEL_EQPLAS
C     [2]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-SE,D_SIG)
C     3.10    UPDATE STATE VARIALBES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS
C     3.11    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS_K, SIG_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 230
C-----------------------------------------------------------------------
235   CONTINUE
C     3.12    UPDATE SOLUTION VARIABLES
C     [1]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [2]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [3]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     3.13    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
      END IF
C     3.14    CONSISTENT TANGENT MODULUS
      CALL CONSISTENT(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
      CE=CEP
      
      RETURN
      END SUBROUTINE CPPM2
C-----------------------------------------------------------------------
C
C     TANGENT MODULUS METHOD (SEMI-IMPLICIT)
C     REF: [12] M. ORTIZ  ET AL.  IJNME   (1985)
      SUBROUTINE TMM(STAT_VAR,EQPLAS,SIG,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), CE(NDIM3,NDIM3)
      DIMENSION SIG_K(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3),
     1          R1(NDIM3), D_PLAS(NDIM3), SE(NDIM3,NDIM3), 
     2          DR1DS(NDIM3,NDIM3), DR1DE(NDIM3), DDFDEDS(NDIM3),
     4          DEL_SIG(NDIM3), D_SIG(NDIM3), CEP(NDIM3,NDIM3)
      DIMENSION XHI(NDIM3,NDIM3), XH(NDIM3,NDIM3)
      DIMENSION DEL_SIG0(NDIM3), DEL_SIG1(NDIM3)

      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      COMMON /KITER/ MAX_ITER
C-----------------------------------------------------------------------
C     4.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_PLAS      = 0.D0
      DEL_EQPLAS  = 0.D0
      DEL_SIG     = 0.D0
      DEL_SIG0    = 0.D0
      DEL_SIG1    = 0.D0
      FLAG        = 0.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C     NUMERICAL SCHEME PARAMETER
C     TH=1    : IMPLICIT INTEGRATION
C     TH=0    : EXPLICIT INTEGRATION
C     TH=0.5  : TRAPEZOIDAL RULE
      TH=0.5D0
C-----------------------------------------------------------------------
240   CONTINUE
C     4.1     GRADIENTS OF YIELD SURFACE: DFDS
      CALL GRAD(2,2,STAT_VAR,SIG_K,DFDS,DDFDDS)
C     4.2     EFFECTIVE STRESS & FLOW STRESS
C     [1] EFFECTIVE STRESS
      CALL EFF_SIG(STAT_VAR, SIG_K, SIG_BAR_K)
C     [2] FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, EQPLAS_K, FLOW_SIG_K)
C     4.3     RESIDUAL FUNCTIONS
C     [1]     ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_EQPLAS*DFDS(I) !(1,3)
          R1_NORM=R1_NORM+R1(I)**2
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     YIELD CONDITION
      R2=SIG_BAR_K-FLOW_SIG_K
      IF(ITER.EQ.0) R0=R2
C-----------------------------------------------------------------------
C     4.4    CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 245
      END IF
      IF(ITER .GT. MAX_ITER) THEN
          WRITE(*,*) '#DIVERGENCE', R1_NORM, DABS(R2)
          WRITE(*,*) '#TOLERANCE', TOL, TOL
          FLAG=1.D0
      END IF
      IF(ISNAN(R2)) THEN
          WRITE(*,*) '#ERROR: SIG_BAR_K=NaN'
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          WRITE(*,*) 'INITIAL RESIDUAL', R0
          STOP
      END IF
C-----------------------------------------------------------------------
C     4.5     DENOMINATOR
C     [0]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, SE)
C     [1]     INVERSE HESSIAN MATRIX: XHI
      DO I=1, NDIM3
      DO J=1, NDIM3
          XHI(I,J)= SE(I,J)+D_EQPLAS*DDFDDS(I,J)
      END DO
      END DO
C     [2]     HESSIAN MATRIX: XH
      CALL INVERSE(XHI, NDIM3, XH)
C     [3]     KINEMATIC COMPONENTS: DFDE & DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS_K,SIG_K,DFDE,DDFDEDS)
C     [4]     DR1DE: XK
      DR1DE=D_EQPLAS*DDFDEDS+DFDS
C     [5]     DFDS_XH_DR1DE / DFDS_CE_DFDS
      DFDS_XH_DR1DE=0.D0
      DFDS_CE_DFDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_XH_DR1DE=DFDS_XH_DR1DE+DFDS(I)*XH(I,J)*DR1DE(J)
          DFDS_CE_DFDS=DFDS_CE_DFDS+DFDS(I)*CE(I,J)*DFDS(J)
      END DO
      END DO
C     [6]     DENOMINATOR
      XDENOM=TH*DFDS_XH_DR1DE+(1.D0-TH)*DFDS_CE_DFDS-DFDE
C-----------------------------------------------------------------------
C     4.6     NUMERATOR
C     [1]     DFDS_XH_R1
      DFDS_XH_R1=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_XH_R1= DFDS_XH_R1 + DFDS(I)*XH(I,J)*R1(J)
      END DO
      END DO
C     [2]     NUMERATOR
      XNUMER= R2 - TH*DFDS_XH_R1
C-----------------------------------------------------------------------
C     4.7     STEP-SIZE CONTROL
      VAL=XNUMER/XDENOM
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     4.8     STEP-SIZE OF SOLUTION VARIABLES AT A SUB-STEP
C     [1]     EQUIVALENT PLASTIC STRAIN INCREMENT
C     DEL_EQPLAS(K+1/2)
      DEL_EQPLAS= ALPHA*VAL
C     [2]     PLASTIC CORRECTOR INCREMENT
      DEL_SIG1= MATMUL(-XH,R1+DR1DE*DEL_EQPLAS)
      DEL_SIG0= MATMUL(-CE,DFDS*DEL_EQPLAS)
      DEL_SIG=TH*DEL_SIG1+(1.D0-TH)*DEL_SIG0
C-----------------------------------------------------------------------
C     4.9     SOLUTION VARIABLE INCREMENTS AT A TIME-STEP
C     [1]     EQUIVALENT PALSTIC STRAIN INCREMENT AT A TIME-STEP
      D_EQPLAS=D_EQPLAS+DEL_EQPLAS
C     [2]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-SE,D_SIG)
C-----------------------------------------------------------------------
C     4.10    UPDATE STATE VARIALBES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS
C-----------------------------------------------------------------------
C     4.11    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS_K, SIG_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 240
C-----------------------------------------------------------------------
245   CONTINUE
C     4.12    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     4.13    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
      END IF
C     4.14    CONSISTENT TANGENT MODULUS
      CALL CONSISTENT(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)

      RETURN
      END SUBROUTINE TMM
C-----------------------------------------------------------------------
      SUBROUTINE LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,DEL_EQPLAS)
C     # THIS SUBROUTINE CONTROLS THE STEP SIZE OF THE SOLUTION VARIABLE
C       BASED ON THE LINE-SEARCH TECHNIQUE.
C     REF: [13] W. M. SCHERZINGER  CMAME   (2017)
C          [17] S. Y. YOON ET AL. CMAME (2020)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION R1(NDIM3), CE(NDIM3,NDIM3)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      PARAMETER(TOL=1.D-6)
C     [0] COEFFICIENTS OF LINE-SEARCH
      BETA= 1.D-4
      ETA= 1.D-1
      EPS_MAX=1.D-4
      EM=2.D0
      GAM1= 0.4D0
      GAM2= 9.D0
      GAM3= 0.5D0
C     [1] MERIT FUNCTION: PSI  [13][EQ. (46)]
      CALL DOUBLE_DOT(R1,R1,NDIM3,DOUBLE_R1)
      XMU= CE(NDIM3,NDIM3)
      PSI= (5.D-1)*(DOUBLE_R1**(EM/2.D0)+(DABS(R2)/(2.D0*XMU))**EM)
C     [2] GOLDSTEIN'S CONDITION   [13][EQ. (51)]
      IF(ITER .EQ. 0) THEN
          PSI0= PSI
          ALPHA1= ALPHA
      ELSE
          PSI2= (1.D0-2.D0*BETA*ALPHA)*PSI0
          IF(PSI-PSI2.GT.0.D0) THEN
              ALPHA1= PSI0/(PSI0+PSI)
          ELSE
              ALPHA1= ALPHA
          END IF
      END IF
      VAL2=ALPHA*DEL_EQPLAS
C     [3] ADDITIONAL STEP SIZE CONTROL CONDITION
      IF(VAL2.LT.0.D0 .OR. VAL2.GT.EPS_MAX) THEN
C      IF(DABS(VAL2).GT.EPS_MAX) THEN
          ALPHA1= GAM1/((1.D0+GAM2*(PSI/PSI0)**GAM3))
      ELSE
          ALPHA1=1.D0
      END IF
C-----------------------------------------------------------------------
      PSIK=PSI
      ALPHA= MAX(ETA, ALPHA1)
      RETURN
      END SUBROUTINE LINESEARCH
C-----------------------------------------------------------------------
C
C     THE FIRST-ORDER EULER FORWARD METHOD
C     REF: F. DUNNE ET AL. (2005)
      SUBROUTINE FOEFM(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), DSTRAN(NDIM3), SIG(NDIM3), 
     1          CE(NDIM3,NDIM3)
      DIMENSION DFDS(NDIM3), DDFDDS(NDIM3,NDIM3), D_SIG(NDIM3),
     1          DELAS(NDIM3), CEP(NDIM3,NDIM3), CE_DFDS(NDIM3),
     2          CE_DFDS_CE_DFDS(NDIM3,NDIM3), DDFDEDS(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      COMMON /KITER/ MAX_ITER
      PARAMETER(TOL=1.D-6)
C-----------------------------------------------------------------------
C     1.1     GRADIENT OF YIELD SURFACE: DFDS, DDFDDS
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
C     1.2     EFFECTIVE STRESS & FLOW STRESS
C     [1] EFFECTIVE STRESS
      CALL EFF_SIG(STAT_VAR, SIG, SIG_BAR)
C     [2] FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     1.3     KINEMATIC COMPONENT: DFDE, DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS,SIG,DFDE,DDFDEDS)
C     1.4     THE INCREMENT OF EFFECTIVE PLASTIC STRAIN
C     [1]     NUMERATOR
      DFDS_CE_DSTRAN=0.D0
      DO I=1,NDIM3
      DO J=1,NDIM3
          DFDS_CE_DSTRAN=DFDS_CE_DSTRAN+DFDS(I)*CE(I,J)*DSTRAN(J)
      END DO
      END DO
C     [2]   DENOMINATOR
      DFDS_CE_DFDS=0.D0
      DO I=1,NDIM3
      DO J=1,NDIM3
          DFDS_CE_DFDS=DFDS_CE_DFDS+DFDS(I)*CE(I,J)*DFDS(J)
      END DO
      END DO
      XDENOM=DFDS_CE_DFDS-DFDE
C
      D_EQPLAS=DFDS_CE_DSTRAN/XDENOM
C     1.5     THE INCREMENT OF STRESS TENSOR BASED ON HOOKE'S LAW
      DELAS=DSTRAN-D_EQPLAS*DFDS
      D_SIG=MATMUL(CE,DELAS)

C     1.6     NUMERICAL INTEGRATION
      SIG=SIG+D_SIG
      EQPLAS=EQPLAS+D_EQPLAS
      CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
C-----------------------------------------------------------------------
C     1.7     CONTINUUM TANGENT MODULUS
      CALL CONTINUUM(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
      CE=CEP

      RETURN
      END SUBROUTINE FOEFM
C-----------------------------------------------------------------------
C
C     NON-ITERATIVE STRESS UPDATE METHOD
C     REF: S. YOON AND F. BARLAT, IJMS, 2023
      SUBROUTINE NSPM(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), DSTRAN(NDIM3), SIG(NDIM3),
     1          CE(NDIM3,NDIM3)
      DIMENSION SIG0(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3), 
     1          CEP(NDIM3,NDIM3), D_SIG(NDIM3), SIGT(NDIM3),
     2          SE(NDIM3,NDIM3), DPLAS(NDIM3), DELAS(NDIM3)
      DIMENSION D_SIG1(NDIM3), D_SIG2(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC, KLAYER
      PARAMETER(TOL=1.D-6)
C-----------------------------------------------------------------------
C     2.0     INITIALIZATION
C     [1]     STRESS
      SIG0=SIG
C     [2]     DOUBLE_SIG
      CALL DOUBLE_DOT(SIG, SIG, NDIM3, XSIG)
      IF(XSIG.EQ.0.D0) SIG=SIG0+MATMUL(CE,DSTRAN)
C-----------------------------------------------------------------------
C     2.1     ELASTOPLASTIC TANGENT MODULUS
      CALL CONTINUUM(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
C-----------------------------------------------------------------------
C     2.2     CALCULATION OF THE STRESS TENSOR INCREMENT
      IF(D_EQPLAS .EQ. 0.D0) THEN
C     [1] EFFECTIVE STRESS & FLOW STRESS
          CALL EFF_SIG(STAT_VAR, SIG0, SIG_BAR)
          CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     [2] STRESS UPDATE WITH ELASTIC STIFFNESS
          D_SIG1=MATMUL(CE, DSTRAN)
          SIG=SIG0+D_SIG1
          CALL EFF_SIG(STAT_VAR, SIG, SIG_BARE)
C     [3] STRESS UPDATE WITH ELASTOPLASTIC STIFFNESS
          D_SIG2=MATMUL(CEP,DSTRAN)
          SIG=SIG0+D_SIG2
          CALL EFF_SIG(STAT_VAR, SIG, SIG_BAREP)
C     [4] PORTION OF ELASTICITY
          GAM=DABS(FLOW_SIG-SIG_BAR)/DABS(SIG_BARE-SIG_BAR)
C     [5] INTERMEDIATE STRESS UPDATE
          D_SIG=GAM*D_SIG1+(1.D0-GAM)*D_SIG2
          SIG=SIG0+D_SIG
      ELSE
          D_SIG=MATMUL(CEP,DSTRAN)
          SIG=SIG0+D_SIG
      END IF
C-----------------------------------------------------------------------
C     2.3     GRADIENT OF YIELD SURFACE: DFDS, DDFDDS
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
C-----------------------------------------------------------------------
C     2.4     POSITIVE-DEFINITE EFFECTIVE STRAIN INCREMENT
C     [1]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, SE)
C     [2]     ELASTIC STRAIN INCREMENT
      DELAS=MATMUL(SE,D_SIG)
      DPLAS=DSTRAN-DELAS
C     [3]     INCREMENT OF EFFECTIVE PLASTIC STRAIN
      CALL EFF_STR(DPLAS, DFDS, D_EQPLAS)
C     [4]     EFFECTIVE STRAIN INTEGRATION
      EQPLAS=EQPLAS+D_EQPLAS
C-----------------------------------------------------------------------
C     2.5     UPDATE STATE VARIABLES OF ANISOTROPIC HARDENING MODEL
      CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
      CALL EFF_SIG(STAT_VAR, SIG, SIG_BAR)
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C-----------------------------------------------------------------------
C     2.6     STRESS PROJECTION
      SIGT=SIG/SIG_BAR
      SIG=SIGT*FLOW_SIG
C-----------------------------------------------------------------------
C     2.7     RETURN ELASTOPLASTIC TANGENT STIFFNESS TENSOR
      CE=CEP

      RETURN
      END SUBROUTINE NSPM
C-----------------------------------------------------------------------
C
C     NON-ITERATIVE STRESS UPDATE METHOD
C     REF: S. YOON AND F. BARLAT, IJMS, 2023
      SUBROUTINE NSPM2(STAT_VAR,DSTRAN,D_EQPLAS,EQPLAS,SIG,CE)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), DSTRAN(NDIM3), SIG(NDIM3),
     1          CE(NDIM3,NDIM3)
      DIMENSION SIG0(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3), 
     1          CEP(NDIM3,NDIM3), D_SIG(NDIM3), SIGT(NDIM3),
     2          SE(NDIM3,NDIM3), DPLAS(NDIM3), DELAS(NDIM3)
      DIMENSION D_SIG1(NDIM3), D_SIG2(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      PARAMETER(TOL=1.D-6)
C-----------------------------------------------------------------------
C     3.0     INITIALIZATION
C     [1]     STRESS
      SIG0=SIG
C     [2]     DOUBLE_SIG
      CALL DOUBLE_DOT(SIG, SIG, NDIM3, XSIG)
      IF(XSIG.EQ.0.D0) SIG=SIG0+MATMUL(CE,DSTRAN)
C-----------------------------------------------------------------------
C     3.1     ELASTOPLASTIC TANGENT MODULUS
      CALL CONSISTENT(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
C-----------------------------------------------------------------------
C     3.2     CALCULATION OF THE STRESS TENSOR INCREMENT
      IF(D_EQPLAS .EQ. 0.D0) THEN
C     [1] EFFECTIVE STRESS & FLOW STRESS
          CALL EFF_SIG(STAT_VAR, SIG0, SIG_BAR)
          CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     [2] STRESS UPDATE WITH ELASTIC STIFFNESS
          D_SIG1=MATMUL(CE, DSTRAN)
          SIG=SIG0+D_SIG1
          CALL EFF_SIG(STAT_VAR, SIG, SIG_BARE)
C     [3] STRESS UPDATE WITH ELASTOPLASTIC STIFFNESS
          D_SIG2=MATMUL(CEP,DSTRAN)
          SIG=SIG0+D_SIG2
          CALL EFF_SIG(STAT_VAR, SIG, SIG_BAREP)
C     [4] PORTION OF ELASTICITY
          GAM=DABS(FLOW_SIG-SIG_BAR)/(SIG_BARE-SIG_BAR)
C     [5] INTERMEDIATE STRESS UPDATE
          D_SIG=GAM*D_SIG1+(1.D0-GAM)*D_SIG2
          SIG=SIG0+D_SIG
      ELSE
          D_SIG=MATMUL(CEP,DSTRAN)
          SIG=SIG0+D_SIG
      END IF
C-----------------------------------------------------------------------
C     3.3     GRADIENT OF YIELD SURFACE: DFDS, DDFDDS
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
C-----------------------------------------------------------------------
C     3.4     POSITIVE-DEFINITE EFFECTIVE STRAIN INCREMENT
C     [1]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, SE)
C     [2]     ELASTIC STRAIN INCREMENT
      DELAS=MATMUL(SE,D_SIG)
      DPLAS=DSTRAN-DELAS
C     [3]     INCREMENT OF EFFECTIVE PLASTIC STRAIN
      CALL EFF_STR(DPLAS, DFDS, D_EQPLAS)
C     [4]     EFFECTIVE STRAIN INTEGRATION
      EQPLAS=EQPLAS+D_EQPLAS
C-----------------------------------------------------------------------
C     3.5     UPDATE STATE VARIABLES OF ANISOTROPIC HARDENING MODEL
      CALL UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, D_EQPLAS)
      CALL EFF_SIG(STAT_VAR, SIG, SIG_BAR)
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C-----------------------------------------------------------------------
C     3.6     STRESS PROJECTION
      SIGT=SIG/SIG_BAR
      SIG=SIGT*FLOW_SIG
C-----------------------------------------------------------------------
C     3.7     RETURN ELASTOPLASTIC TANGENT STIFFNESS TENSOR
      CE=CEP

      RETURN
      END SUBROUTINE NSPM2
C-----------------------------------------------------------------------
      SUBROUTINE CONTINUUM(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
C     #   THIS SUBROUTINE RETURNS CONTINUUM TANGENT MODULUS.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), CE(NDIM3,NDIM3),
     1          CEP(NDIM3,NDIM3)
      DIMENSION DFDS(NDIM3), DDFDDS(NDIM3,NDIM3), DDFDEDS(NDIM3),
     1          CE_DFDS(NDIM3), XNUMER(NDIM3,NDIM3)
      
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     1.  PREQUISITES
C     [1]     GRADIENT OF PLASTIC POTENTIAL: DFDS
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
C     [2]     HARDENING RATE: DHDE
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     [3]     KINEMATIC COMPONENT: DFDE, DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS,SIG,DFDE,DDFDEDS)
C     2.  DENOMIATOR OF TANGENT STIFFNESS TENSOR
C     [1]     DFDS_CE_DFDS
      DFDS_CE_DFDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_CE_DFDS= DFDS_CE_DFDS+DFDS(I)*CE(I,J)*DFDS(J)
      END DO
      END DO
C     [2]     DENOMINATOR
      XDENOM=DFDS_CE_DFDS-DFDE
C     3.   NUMERATOR OF TANGENT STIFFNESS TENSOR
C     [1]     C_DFDS
      CE_DFDS=MATMUL(CE,DFDS)
C     [2]     XNUMER: C_DFDS_C_DFDS
      CALL DYAD(CE_DFDS, CE_DFDS, NDIM3, XNUMER)
C     4.  ELASTOPLASTIC TANGENT STIFFNESS
      DO I=1, NDIM3
      DO J=1, NDIM3
          CEP(I,J)=CE(I,J)-XNUMER(I,J)/XDENOM
      END DO
      END DO

      RETURN
      END SUBROUTINE CONTINUUM
C-----------------------------------------------------------------------
      SUBROUTINE CONSISTENT(STAT_VAR,D_EQPLAS,EQPLAS,SIG,CE,CEP)
C     #   THIS SUBROUTINE RETURNS CONSISTENT TANGENT MODULUS.
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), CE(NDIM3,NDIM3),
     1          CEP(NDIM3,NDIM3)
      DIMENSION DFDS(NDIM3), DDFDDS(NDIM3,NDIM3), SE(NDIM3,NDIM3),
     1          XHI(NDIM3,NDIM3), XH(NDIM3,NDIM3), DDFDEDS(NDIM3),
     2          XK(NDIM3), XH_XK(NDIM3), XH_DFDS(NDIM3), 
     3          XNUMER(NDIM3,NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     1.  PREQUISITES
C     [1]     GRADIENT OF PLASTIC POTENTIAL: DFDS & DDFDDS
      CALL GRAD(2,2,STAT_VAR,SIG,DFDS,DDFDDS)
C     [2]     HARDENING RATE: DHDE
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     [3]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, SE)
C     [4]     INVERSE HESSIAN MATRIX: XHI
      DO I=1, NDIM3
      DO J=1, NDIM3
          XHI(I,J)= SE(I,J)+D_EQPLAS*DDFDDS(I,J)
      END DO
      END DO
C     [5]     HESSIAN MATRIX: XH
      CALL INVERSE(XHI, NDIM3, XH)
C     2.  DENOMINATOR
C     [1]     KINEMATIC COMPONENT: DFDE, DDFDEDS
      CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS,SIG,DFDE,DDFDEDS)
C     [2]     XK: DR1DE
      XK=D_EQPLAS*DDFDEDS+DFDS
C     [3]     DFDS_XH_XK
      DFDS_XH_XK=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_XH_XK=DFDS_XH_XK+DFDS(I)*XH(I,J)*XK(J)
      END DO
      END DO
C     [4]     XDENOM
      XDENOM=DFDS_XH_XK-DFDE
C     3.  NUMERATOR
C     [1]     XH_XK
      XH_XK=MATMUL(XH,XK)
C     [2]     XH_DFDS
      XH_DFDS=MATMUL(XH,DFDS)
C     [3]     XNUMER:XH_XK_XH_DFDS
      CALL DYAD(XH_XK, XH_DFDS, NDIM3, XNUMER)
C     4.  ELASTOPLASTIC CONSISTENT TANGENT STIFFNESS
      CEP=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          CEP(I,J)=XH(I,J)-XNUMER(I,J)/XDENOM
      END DO
      END DO

      RETURN
      END SUBROUTINE CONSISTENT
C-----------------------------------------------------------------------
      SUBROUTINE KINEMATIC_COMPONENT(STAT_VAR,EQPLAS,SIG,DFDE,DDFDEDS)
C     #   THIS SUBROUTINE CALCULATES THE KINEMATIC COMPONENTS: DFDE & DDFDEDS
C         VIA THE PSEUDO-UPDATE OF STATE VARIABLE.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3), DDFDEDS(NDIM3)
      DIMENSION DFDS(NDIM3), DDFDDS(NDIM3,NDIM3)
      DIMENSION STAT_VAR1(NDIM6), SIG1(NDIM3), DFDS1(NDIM3), 
     1          DDFDDS1(NDIM3,NDIM3)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     [0] PREQUISITES
      D_EQPLAS0=1.D-9
      EQPLAS1=EQPLAS+D_EQPLAS0
      DO I=1,NDIM6
          STAT_VAR1(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM3
          SIG1(I)=SIG(I)
      END DO
c     [1] PSEUDO-UPDATE OF STATE VARIABLE
      CALL UPDATE_HARD_VAR(STAT_VAR1, EQPLAS1, SIG1, D_EQPLAS0)
C     [2] DFDE: THE FIRST-ORDER KINEMATIC COMPONENT + HARDENING RATE
      CALL YIELD_CONDITION(STAT_VAR, EQPLAS, SIG, FVAL)
      CALL YIELD_CONDITION(STAT_VAR1, EQPLAS1, SIG1, FVAL1)
      DFDE=(FVAL1-FVAL)/D_EQPLAS0
C     [3] DDFDEDS: THE SECOND-ORDER KINEMATIC COMPONENT
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
      CALL GRAD(2,1,STAT_VAR1,SIG1,DFDS1,DDFDDS1)
      DO I=1,NDIM3
          DDFDEDS(I)=(DFDS1(I)-DFDS(I))/D_EQPLAS0
      END DO

      RETURN
      END SUBROUTINE KINEMATIC_COMPONENT
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.3  HARDENING LAWS                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     #   THIS SUBROUTINE RETURNS FLOW STRESS
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION STAT_VAR(NDIM6)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

      IF(HARD_PAR.EQ.1.D0 .OR. HARD_PAR.EQ.2.D0) THEN
          IF(FLOW_PAR .EQ. 7.) THEN
              WRITE(*,*) '#ERROR: RGBV IS NOT AVAILABLE IN KINEMATIC HAR
     1                    DENING'
              STOP
          END IF
      END IF

      IF(HARD_PAR.EQ.1.D0) THEN
          IF(H4.EQ.0.D0) THEN
C     #   REF: T.J. PARK, IJSS, 2012
              CALL ISO_HARD(STAT_VAR, 0.D0, FLOW_SIG0)
              CALL ISO_HARD(STAT_VAR, EQPLAS, FLOW_SIG)
              FLOW_SIG=H3*(FLOW_SIG+FLOW_SIG0)
              DHDE=H3*DHDE
          ELSE
              CALL CB_HARD(EQPLAS, FLOW_SIG)
          END IF
      ELSEIF(HARD_PAR.EQ.2.D0) THEN
          FLOW_SIG= H2
          DHDE= 0.D0
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN
          CALL ISO_HARD(STAT_VAR, EQPLAS, FLOW_SIG)
          GP= STAT_VAR(5)
          DHDE= GP*DHDE
          FLOW_SIG=GP*FLOW_SIG
      ELSE
          CALL ISO_HARD(STAT_VAR, EQPLAS, FLOW_SIG)
      END IF
      DHDE=0.D0
      RETURN
      END SUBROUTINE FLOW_STRESS
C-----------------------------------------------------------------------
      SUBROUTINE ISO_HARD(STAT_VAR, EQPLAS, FLOW_SIG)
C     #   THIS SUBROUTINE RETURNS FLOW STRESS OF ISOTROPIC HARDENING LAW.
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION STAT_VAR(NDIM6)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KIHARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10

      IF(FLOW_PAR.EQ.1.D0) THEN
C     #   SWIFT HARDENING LAW #
          FLOW_SIG=P1*(P2+EQPLAS)**P3             ! H(¥Å_BAR) = P1*(P2 + ¥Å_BAR)^P3
          DHDE    =P1*P3*(P2+EQPLAS)**(P3-1.D0)    ! ¡ÓH/¡Ó(¥Å_BAR) = P1*P3*(P2 + ¥Å_BAR)^(P3-1)
      ELSEIF(FLOW_PAR.EQ.2.D0) THEN
C     #   VOCE HARDENING LAW #
          FLOW_SIG=P1-P2*DEXP(-P3*EQPLAS)    ! H(¥Å_BAR) = P1-P2*EXP(-P3*¥Å_BAR)]
          DHDE    =P2*P3*DEXP(-P3*EQPLAS)    ! ¡ÓH/¡Ó(¥Å_BAR) = P2*P3*EXP(-P3*¥Å_BAR)
      ELSEIF(FLOW_PAR.EQ.3.D0) THEN
C     #   MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=P1+P2*EQPLAS +P3*(1.D0-DEXP(-P4*EQPLAS))
          DHDE    =P2+ P3*P4*DEXP(-P4*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.4.D0) THEN
C     #   SWIFT+MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=P1*(P2+EQPLAS)**P3+P4+P5*EQPLAS
     1             +P6*(1.D0-DEXP(-P7*EQPLAS))
          DHDE    =P5+P1*P3*(P2 + EQPLAS)**(P3-1.D0)
     1             +P6*P7*DEXP(-P7*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.5.D0) THEN
C     #   HOCKETT-SHERBY + LINEAR HARDENING LAW #
          IF(EQPLAS .EQ. 0.D0) THEN
              EQPLAS0= 5.D-10
          ELSE
              EQPLAS0= EQPLAS
          END IF
          FLOW_SIG=P1-(P1-P2)*DEXP(-P3*EQPLAS0**P4) + P5*EQPLAS0
          DHDE=(P1-P2)*P3*P4*EQPLAS0**(P4-1.D0)
     1        *DEXP(-P3*EQPLAS0**P4) + P5
      ELSEIF(FLOW_PAR.EQ.6.D0) THEN
C     #   COMBINED SWIFT-VOCE #
          FLOW_SIG= P1*(P2*(P3+EQPLAS)**P4)
     1             +(1.D0-P1)*(P5-P6*DEXP(-P7*EQPLAS))
          DHDE= P1*P2*P4*(P3+EQPLAS)**(P4-1.D0)
     1         +(1.D0-P1)*P6*P7*DEXP(-P7*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.7.D0) THEN
C     #   DISLOCATION-BASED HARDENING: RGBV [15]#
          IF(EQPLAS .EQ. 0.D0) THEN
              RHO_TOT=P10
              DRDE=0.D0
          ELSE
              RHO_TOT=STAT_VAR(NDIM7-NDIM8+7)
              DRDE=DABS(STAT_VAR(NDIM7-NDIM8+8))
          END IF
C     #   RESOLVED SHEAR STRESS
          XAL=0.5D0 ! DISLOCATION INTERACTION
          TAU=P2+XAL*P3*P4*DSQRT(RHO_TOT)
C     #   WEI WEN BACK-STRESS
C          FLOW_SIG=P1*(TAU*(1.D0-P10*XX))
C          DHDE=0.5D0*P1*P3*P4*P5/DSQRT(RHO_TOT)*DRDE
C     #   ORIGINAL EQUATION
          FLOW_SIG=P1*(TAU)
          DHDE=0.5D0*P1*XAL*P3*P4/DSQRT(RHO_TOT)*DRDE
      END IF

      RETURN
      END SUBROUTINE ISO_HARD
C-----------------------------------------------------------------------
      SUBROUTINE CB_HARD(EQPLAS, FLOW_SIG)
C     #   THIS SUBROUTINE RETURNS STATE VARIABLES OF KINEMATIC HARDENING MODELS.
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   REF: J.Y. LEE, IJSS, 2012
      IF(FLOW_PAR.EQ.1.D0) THEN
C     #   SWIFT HARDENING LAW #
          FLOW_SIG=H3*(H4+EQPLAS)**H5             ! H(¥Å_BAR) = H3*(H4 + ¥Å_BAR)^H5
          DHDE    =H3*H5*(H4+EQPLAS)**(H5-1.D0)    ! DH/D(¥Å_BAR) = H3*H5*(H4 + ¥Å_BAR)^(H5-1)
      ELSEIF(FLOW_PAR.EQ.2.D0) THEN
C     #   VOCE HARDENING LAW #
          FLOW_SIG=H3-H4*DEXP(-H5*EQPLAS)    ! H(¥Å_BAR) = H3-H4*EXP(-H5*¥Å_BAR)]
          DHDE    =H4*H5*DEXP(-H5*EQPLAS)    ! DH/D(¥Å_BAR) = H4*H5*EXP(-H5*¥Å_BAR)
      ELSEIF(FLOW_PAR.EQ.3.D0) THEN
C     #   MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=H3+H4*EQPLAS +H5*(1.D0-DEXP(-H6*EQPLAS))
          DHDE    =H4+ H5*H6*DEXP(-H6*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.4.D0) THEN
C     #   SWIFT+MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=H3*(H4+EQPLAS)**H5+H6+H7*EQPLAS
     1             +H8*(1.D0-DEXP(-H9*EQPLAS))
          DHDE    =H7+H3*H5*(H4 + EQPLAS)**(H5-1.D0)
     1             +H8*H9*DEXP(-H9*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.5.D0) THEN
C     #   HOCKETT-SHERBY + LINEAR HARDENING LAW #
          IF(EQPLAS .EQ. 0.D0) THEN
              EQPLAS0= 1.D-9
          ELSE
              EQPLAS0= EQPLAS
          END IF
          FLOW_SIG=H3-(H3-H4)*DEXP(-H5*EQPLAS0**H6) + H7*EQPLAS0
          DHDE=(H3-H4)*H5*H6*EQPLAS0**(H6-1.D0)
     1            *DEXP(-H5*EQPLAS0**H6) + H7
      ELSEIF(FLOW_PAR.EQ.6.D0) THEN
C     #   COMBINED SWIFT-VOCE #
          FLOW_SIG= H1*(H2*(H3+EQPLAS)**H4)
     1             +(1.D0-H1)*(H5-H6*DEXP(-H7*EQPLAS))
          DHDE= H1*H2*H4*(H3+EQPLAS)**(H4-1.D0)
     1             +(1.D0-H1)*H6*H7*DEXP(-H7*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.7.D0) THEN
C     #   DISLOCATION-BASED HARDENING: RGBV [15]#
          PRINT *, '#MSG: DISLOCATION-BASED HARDENING IS NOT AVAILABLE 
     1              IN KINEMATIC HARDENING'
          STOP
      END IF
      RETURN
      END SUBROUTINE CB_HARD
C-----------------------------------------------------------------------
      SUBROUTINE YU_HARD(EQPLAS, FLOW_SIG)
C     #   THIS SUBROUTINE RETURNS FLOW STRESS OF ISOTROPIC HARDENING TERM OF YU.
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
C     #   ORITINAL HARDENING LAW#
      IF(FLOW_PAR.EQ.0.D0) THEN
          FLOW_SIG= H7*(1.D0-DEXP(-H4*EQPLAS))
          DHDE    = H7*H4*DEXP(-H4*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.1.D0) THEN
C     #   SWIFT HARDENING LAW #
          FLOW_SIG=H7*((H8+EQPLAS)**H9-H8**H9)             ! H(E_BAR) = H7*(H8 + E_BAR)^H9
          DHDE    =H7*H9*(H8+EQPLAS)**(H9-1.D0)    ! DH/D(E_BAR) = H7*H9*(H8 + E_BAR)^(H9-1)
      ELSEIF(FLOW_PAR.EQ.2.D0) THEN
C     #   VOCE HARDENING LAW #
          FLOW_SIG=H7-H8*DEXP(-H9*EQPLAS)    ! H(E_BAR) = H7-H8*EXP(-H9*E_BAR)]
          DHDE    =H8*H9*DEXP(-H9*EQPLAS)    ! DH/D(E_BAR) = H8*H9*EXP(-H9*E_BAR)
      ELSEIF(FLOW_PAR.EQ.3.D0) THEN
C     #   MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=H7+H8*EQPLAS +H9*(1.D0-DEXP(-H10*EQPLAS))
          DHDE    =H8+ H9*H10*DEXP(-H10*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.4.D0) THEN
C     #   SWIFT+MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG= H7*(EQPLAS)**H8- H7**H8
     1             +H9*EQPLAS+H10*(1.D0-DEXP(-H11*EQPLAS))
          DHDE= H7*H8*(EQPLAS)**(H8-1.D0)
     1         +H9+H10*H11*DEXP(-H11*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.5.D0) THEN
C     #   HOCKETT-SHERBY + LINEAR HARDENING LAW #
          IF(EQPLAS .EQ. 0.D0) THEN
              EQPLAS0= 1.D-9
          ELSE
              EQPLAS0= EQPLAS
          END IF
          FLOW_SIG=H7-(H7-H8)*DEXP(-H9*EQPLAS0**H10)+H11*EQPLAS0
          DHDE=(H7-H8)*H9*H10*EQPLAS0**(H10-1.D0)
     1        *DEXP(-H9*EQPLAS0**H10)+H11
C     #   COMBINED HOLLOMON-VOCE #
      ELSEIF(FLOW_PAR.EQ.6.D0) THEN
          FLOW_SIG= H7*(H8*(EQPLAS)**H9- H8**H9)
     1             +(1.D0-H7)*(H10*(1.D0-DEXP(-H11*EQPLAS)))
          DHDE= H7*H8*H9*(EQPLAS)**(H9-1.D0)
     1         +(1.D0-H7)*H10*H11*DEXP(-H11*EQPLAS)
      END IF

      RETURN
      END SUBROUTINE YU_HARD
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.4  EQUIVLANET STRESS                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE YIELD_CONDITION(STAT_VAR, EQPLAS, SIG, FVAL)
C     #   THIS SUBROUTINE RETURNS THE RESIDUAL VALUE OF YIELD CONDITION.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     [1] EFFECTIVE STRESS
      CALL EFF_SIG(STAT_VAR, SIG, SIG_BAR)
C     [2] FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     [3] RESIDUAL VALUE
      FVAL=SIG_BAR-FLOW_SIG
      
      RETURN
      END SUBROUTINE YIELD_CONDITION
C-----------------------------------------------------------------------
      SUBROUTINE EFF_SIG(STAT_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS EQUIVALENT STRESS FROM CONSTITUTIVE MODEL.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM4), STAT_VAR(NDIM6)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     [1] PRE-CONDITION
      CALL DOUBLE_DOT(SIG, SIG, NDIM3, XSIG)
      IF(XSIG .EQ. 0.D0) THEN
          SIG_BAR=0.D0
          RETURN                      ! -> RETURN(END THIS SUBROUTINE) HERE!!
      END IF
C     [2] EFFECTIVE STRESS
      IF(HARD_PAR .LT. 3.D0) THEN ! HAH: OFF
          CALL YLD(SIG, DEV_SIG, SIG_BAR, 1)
      ELSEIF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
          CALL HAH11(STAT_VAR, SIG, SIG_BAR)
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
          CALL HAH14(STAT_VAR, SIG, SIG_BAR)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
          CALL HAH20(STAT_VAR, SIG, SIG_BAR)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
          CALL HEXAH(STAT_VAR, SIG, SIG_BAR)
      END IF

      RETURN
      END SUBROUTINE EFF_SIG
C-----------------------------------------------------------------------
      SUBROUTINE YLD(SIG, DEV_SIG, SIG_BAR, IOPT)
C     #   THIS SUBROUTINE RETURNS THE EFFECITVE STRESS OF FE-IMPLEMENTED YIELD FUNCTION
C     #   IOPT=1: ARGUMENT CORRESPONDS TO CAUCHY STRESS.
C         IOPT=2:                         DEVIATORIC STRESS.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM4), DEV_SIG0(NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     [1] PRE-CONDITION
      IF(IOPT .EQ. 1) THEN
          CALL DOUBLE_DOT(SIG,SIG,NDIM3,XSIG)
      ELSEIF(IOPT .EQ. 2) THEN
          CALL DOUBLE_DOT(DEV_SIG,DEV_SIG,NDIM4,XSIG)
      END IF
      IF(XSIG .EQ. 0.D0) THEN
          SIG_BAR=0.D0
          RETURN                      ! -> RETURN(END THIS SUBROUTINE) HERE!!
      END IF

C     [2] STRESS ARGUMENT
      IF(IOPT .EQ. 1) THEN ! INPUT: NORMAL STRESS
          CALL DEVIATORIC(SIG, DEV_SIG)
      END IF

C     [3] YIELD FUNCTIONS
      IF(YLD_PAR .EQ. 1.D0) THEN        ! VON_MISES
          CALL VON_MISES(DEV_SIG, SIG_BAR)
      ELSEIF(YLD_PAR .EQ. 2.D0) THEN    ! HILL48
          CALL HILL_1948(DEV_SIG, SIG_BAR)
      ELSEIF(YLD_PAR .EQ. 3.D0) THEN    ! YLD2000_2D
          IF(NDIM3 .NE. 3) THEN
              PRINT *, '#ERROR: CHECK THE YIELD FUNCTION OPTION!'
              STOP
          ELSE
              CALL DEVIATORIC_TRANS(DEV_SIG, DEV_SIG0, 2)
              CALL YLD2000_2D(DEV_SIG0, SIG_BAR)
          END IF
      ELSEIF(YLD_PAR .EQ. 4.D0) THEN    ! YLD2004/_18P
          IF(NDIM3 .NE. 6) THEN
              PRINT *, '#ERROR: CHECK THE YIELD FUNCTION OPTION!'
              STOP
          ELSE
              CALL YLD2004_18P(DEV_SIG, SIG_BAR)
          END IF
      END IF

      RETURN
      END SUBROUTINE YLD
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.5   ANISOTROPIC HARDENING MODELS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HAH11(STAT_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EQUIVALENT STRESS OF HAH11 MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), G(NDIM5), DEV_H(NDIM4), DEV_SIG(NDIM4)
      DIMENSION STAT_VAR(NDIM6)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   DEFINE COEFFICIENT
      QM=H1

C     [1] STATEV VARIABLES
      DO I=1, NDIM5
          G(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM5+I)
      END DO
C     [2] PREREQUISITES
C     [2.1]   DEVIATORIC STRESS
      CALL DEVIATORIC(SIG,DEV_SIG)
C     [2.2]   DEV_SH=H^:S
      CALL DOUBLE_DOT(DEV_H, DEV_SIG, NDIM4, DEV_SH)
C-----------------------------------------------------------------------
C     [3] STABLE COMPONENT: EFFECTIVE STRESS OF YIELD_FUNCTION
      CALL YLD(SIG, DEV_SIG, PHI, 1)
C     [3] FLUCTUATION COMPONENT: TRUNCATION COMPONENT
C     [3.1] F-FUNCTION
      F1=DSQRT(3.D0/8.D0)*((1.D0/(G(1)**QM)-1.D0)**(1.D0/QM))
      F2=DSQRT(3.D0/8.D0)*((1.D0/(G(2)**QM)-1.D0)**(1.D0/QM))
C     [3.1] FLUCTUATING COMPONENT: ¥Õ_H
      PHI_H1= F1*DABS(DEV_SH-DABS(DEV_SH))
      PHI_H2= F2*DABS(DEV_SH+DABS(DEV_SH))
      PHI_H = PHI_H1**QM+PHI_H2**QM
C-----------------------------------------------------------------------
C     [4] EQUIVALENT STRESS
      SIG_BAR=(XI**QM+PHI_H)**(1.D0/QM)

      RETURN
      END SUBROUTINE HAH11
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HAH11(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES THE HARDENING STATE VARIABLES OF HAH11
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM4), G(NDIM5), DEL_G(NDIM5),
     1          DEV_H(NDIM4), DEL_DEV_H(NDIM4), DEV_S(NDIM4)
      DIMENSION STAT_VAR(NDIM6)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   DEFINE COEFFICIENT
      XK=H2; XK1=H3; XK2=H4; XK3=H5; XK4=H6; XK5=H7

C     #   STATE VARIABLES FOR HAH11 MODEL
      DO I=1, NDIM5
          G(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM5+I)
      END DO
C     # DEFINE YIELD STRENGHTH AND FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, 0.D0, FLOW_SIG0, DHDE)
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     #   CALCULATE THE INCREMENTS OF INTERNAL PARAMETERS OF HAH MODEL
C     S= DEV_SIG(I+1, N+1)
      CALL DEVIATORIC(SIG,DEV_SIG)
C     S^= DEV_S(I+1, N+1)
      CALL NORM_HAH(DEV_SIG,DEV_S)
C     H^(N):S(I+1, N+1)= DOUBLE
      CALL DOUBLE_DOT(DEV_H, DEV_SIG, NDIM4, DOUBLE)
C     [H^:S^] = COS_X
      CALL DOUBLE_DOT(DEV_S, DEV_H, NDIM4, COS_X)
C     COS_X=(8/3)*(S^:H^)
      COS_X=(8.D0/3.D0)*COS_X
      IF(COS_X .LT. 0.D0) THEN
          XLM= -1.D0
      ELSE
          XLM= 1.D0
      ENDIF
C     #   MICROSTRUCTURE DEVIATOR
      DO I=1, NDIM4
          DEL_DEV_H(I)=XLM*XK*(DEV_S(I)-DEV_H(I)*COS_X)
      END DO
      DEV_H=DEV_H+DEL_DEV_H*DEL_EQPLAS
      CALL NORM_HAH(DEV_H, DEV_H)
C     #   REVERSE LOADING VARIABLES
C     [1] BAUSCHINGER EFFECT
      DEL_G(1)= ((1.D0-XLM)/2.D0)*XK1*(G(4)/G(1)-1.D0)
     1         +((1.D0+XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(1))
      DEL_G(2)= ((1.D0-XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(2))
     1         +((1.D0+XLM)/2.D0)*XK1*(G(3)/G(2)-1.D0)
C     [2] PERMANENT SOFTENING
      DEL_G(3)= ((1.D0-XLM)/2.D0)*XK5*(XK4-G(3))
      DEL_G(4)= ((1.D0+XLM)/2.D0)*XK5*(XK4-G(4))

      G=G+DEL_G*DEL_EQPLAS
      IF(G(1) .GT. 1.D0) G(1)=1.D0
      IF(G(2) .GT. 1.D0) G(2)=1.D0

C     #   UPDATE STATE VARIABLES FOR HAH11 MODEL
      DO I=1, NDIM5
          STAT_VAR(I)=G(I)
      END DO
      DO I=1, NDIM4
          STAT_VAR(NDIM5+I)=DEV_H(I)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HAH11
C-----------------------------------------------------------------------
      SUBROUTINE HAH14(STAT_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EFFECTIVE STRESS GENERATED BY HAH (HAH20H) MODEL.
C     #   BOTH HAH20H & HAH20E ARE IMPLEMENTED IN THIS SUBROUTINE.
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), G(NDIM5), DEV_H(NDIM4), DEV_H0(NDIM3),
     1          DEV_SIG(NDIM4)
      DIMENSION SIG_C(NDIM4), SIG_O(NDIM4), SIG_L(NDIM4), SIG_X(NDIM4)
      DIMENSION STAT_VAR(NDIM6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   DEFINE COEFFICIENT
      QM=H1;

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM5+I)
      END DO
C-----------------------------------------------------------------------
C     [1] PREREQUIRED VARIALBES
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)
      CALL DOUBLE_DOT(DEV_SIG, DEV_H, NDIM4, DEV_SH)
      CALL DOUBLE_DOT(DEV_H0, DEV_H, NDIM4, COS_X)
      HVAL=(8.D0/3.D0)
      COS_X=HVAL*COS_X
C-----------------------------------------------------------------------
C     [2] STRESS DEVIATOR DECOMPOSITION
C     [2.1] COLLINEAR STRESS COMPONENTS
      SIG_C=HVAL*DEV_SH*DEV_H
C     [2.2] ORTHOGONAL STRESS COMPONENTS
      SIG_O=DEV_SIG-SIG_C
C     [2.3] STRESS COMPONENTS FOR LATENT HARDENING
      IF(G(6) .EQ. 0.D0) THEN
          WRITE(*,*) '#ERROR: GL= 0.D0'
          STOP
      END IF
      SIG_L=SIG_C+(1.D0/G(6))*SIG_O
C     [2.4] STRESS COMPONENTS FOR CROSS-LOADING CONTRACTION
      SIG_X=(4.D0*(1.D0-G(5)))*SIG_O
C-----------------------------------------------------------------------
C     [3] STABLE COMPONENT
      CALL YLD(SIG, SIG_L, PHI_L, 2)
      CALL YLD(SIG, SIG_X, PHI_X, 2)
      XI=(PHI_L**2.D0 + PHI_X**2.D0)**(1.D0/2.D0)
C-----------------------------------------------------------------------
C     [4] FLUCTUATION COMPONENT: TRUNCATION COMPONENT
C     [4.1] F-FUNCTION
      F1=DSQRT(3.D0/8.D0)*((1.D0/(G(1)**QM)-1.D0)**(1.D0/QM))
      F2=DSQRT(3.D0/8.D0)*((1.D0/(G(2)**QM)-1.D0)**(1.D0/QM))
C     [4.1] FLUCTUATING COMPONENT: ¥Õ_H
      PHI_H1= F1*DABS(DEV_SH-DABS(DEV_SH))
      PHI_H2= F2*DABS(DEV_SH+DABS(DEV_SH))
      PHI_H = PHI_H1**QM+PHI_H2**QM
C-----------------------------------------------------------------------
C     [5] EQUIVALENT STRESS
      SIG_BAR=(XI**QM+PHI_H)**(1.D0/QM)

      RETURN
      END SUBROUTINE HAH14
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HAH14(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES THE HARDENING STATE VARIABLES OF HAH14
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION G(NDIM5), DEV_H(NDIM4), DEV_HP(NDIM4), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM4), DEV_H0(NDIM4), DEV_HS(NDIM4),
     1          DEL_DEV_H(NDIM4), DEL_DEV_HP(NDIM4), DEL_G(NDIM5)
      DIMENSION STAT_VAR(NDIM6)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   DEFINE COEFFICIENT
      XK=H2; XK1=H3; XK2=H4; XK3=H5; XK4=H6; XK5=H7
      XC=H8; XKC=H9; XL=H10; XKL=H11

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM5+I)
      END DO
C     # DEFINE YIELD STRENGHTH AND FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, 0.D0, FLOW_SIG0, DHDE)
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
C     [1] PREREQUISITE
C     1.1 MICROSCTURE DEVIATOR COEFFICIENTS
      XZ  =5.D0
      XKR =0.1D0
      XKRP=0.02D0
      HVAL=8.D0/3.D0
C      HVAL=1.D0
C     1.2 STRAIN-PATH CHANGE PARAMETER
C     1.2.1
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)

C     1.2.2
      CALL DOUBLE_DOT(DEV_H, DEV_H0, NDIM4, COS_X)
      COS_X=HVAL*COS_X
      IF(DABS(COS_X) .GT. 1.D0) COS_X=COS_X/DABS(COS_X)
      IF(DABS(COS_X) .LT. 1.D-7) COS_X=0.D0
      COS_X2=COS_X*COS_X
      IF(COS_X .LT. 0.D0) THEN
          XLM= -1.D0
      ELSE
          XLM= 1.D0
      ENDIF

C     [2] MICROSTRUCTURE DEVIATOR EVOLUTION
C      VAL=XK*XLM*((DABS(COS_X))**(1.D0/XZ)+G(7))
      VAL=XK*XLM*((DABS(COS_X))**(XKR)+XKRP)
      DO I=1, NDIM4
          DEL_DEV_H(I)=VAL*(DEV_H0(I)-DEV_H(I)*COS_X)
      END DO
      DEV_H=DEV_H+DEL_DEV_H*DEL_EQPLAS
      CALL NORM_HAH(DEV_H, DEV_H)
C-----------------------------------------------------------------------
C     [3] REVERSE LOADING VARIABLE EVOLUTION
      DEL_G=0.D0
      DEL_G(1)= ((1.D0-XLM)/2.D0)*XK1*(G(4)/G(1)-1.D0)
     1         +((1.D0+XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(1))
      DEL_G(2)= ((1.D0-XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(2))
     1         +((1.D0+XLM)/2.D0)*XK1*(G(3)/G(2)-1.D0)
      DEL_G(3)= ((1.D0-XLM)/2.D0)*XK5*(XK4-G(3))
      DEL_G(4)= ((1.D0+XLM)/2.D0)*XK5*(XK4-G(4))
C     [5] CROSS-LOADING CONTRACTION EVOLUTION
      DEL_G(5)=XKC*(1.D0+(XC-1.D0)*COS_X2-G(5))

C     [6] LATENT HARDENING EVOLUTION
      VAL=(FLOW_SIG-FLOW_SIG0)/FLOW_SIG
      VAL2=XL*(1.D0-COS_X2)+COS_X2
      DEL_G(6)= XKL*(VAL*(DSQRT(VAL2)-1.D0)+1.D0-G(6))

C     [8] UPDTAE G-VALUES
      G=G+DEL_G*DEL_EQPLAS
      VAL=XK3*(FLOW_SIG0/FLOW_SIG)
      IF(G(1) .GT. 1.D0) G(1)=1.D0
      IF(G(2) .GT. 1.D0) G(2)=1.D0

C     [9] UPDATE STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          STAT_VAR(I)=G(I)
      END DO
      DO I=1, NDIM4
          STAT_VAR(NDIM5+I)=DEV_H(I)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HAH14
C-----------------------------------------------------------------------
      SUBROUTINE HAH20(STAT_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EFFECTIVE STRESS GENERATED BY HAH (HAH20H) MODEL.
C     #   BOTH HAH20H & HAH20E ARE IMPLEMENTED IN THIS SUBROUTINE.
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), G(NDIM5), DEV_H(NDIM4), DEV_HP(NDIM4), 
     1          DEV_SIG(NDIM4)
      DIMENSION SIG_C(NDIM4), SIG_O(NDIM4), SIG_L(NDIM4), SIG_X(NDIM4)
      DIMENSION DEV_H0(NDIM4), DEV_EH(NDIM4), DEV_EH0(NDIM4),
     1          DFDH(NDIM3), DDFDDH(NDIM3, NDIM3),
     2          T(NDIM3, NDIM3), DEV_SIG0(NDIM4)
      DIMENSION STAT_VAR(NDIM6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   DEFINE COEFFICIENT
      QM=H1; PM=H2;

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM5+I)
          DEV_HP(I)=STAT_VAR(NDIM5+NDIM4+I)
      END DO
C     [1] PREREQUIRED VARIALBES
      XI_L=5.D-1
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL DOUBLE_DOT(DEV_SIG, DEV_H, NDIM4, DEV_SH)
C-----------------------------------------------------------------------
C     [2] STRESS DEVIATOR DECOMPOSITION
C     [2.1] COLLINEAR STRESS COMPONENTS
      SIG_C=DEV_SH*DEV_H
C     [2.2] ORTHOGONAL STRESS COMPONENTS
      SIG_O=DEV_SIG-SIG_C
C     [2.3] STRESS COMPONENTS FOR LATENT HARDENING
      IF(G(8) .EQ. 0.D0) THEN
          WRITE(*,*) '#ERROR: GL= 0.D0'
          STOP
      END IF
      VAL=1.D0/(XI_L*(G(8)-1.D0)+1.D0)
      SIG_L=VAL*SIG_C+(1.D0/G(8))*SIG_O
C     [2.4] STRESS COMPONENTS FOR CROSS-LOADING CONTRACTION
      SIG_X=(4.D0*(1.D0-G(7))/G(8))*SIG_O
C-----------------------------------------------------------------------
C     [3] STABLE COMPONENT
      CALL YLD(SIG, SIG_L, PHI_L, 2)
      CALL YLD(SIG, SIG_X, PHI_X, 2)
      XI=(PHI_L**PM + PHI_X**PM)**(1.D0/PM)
C-----------------------------------------------------------------------
C     [4] HAH20H
      IF(HARD_PAR .EQ. 5.D0) THEN ! HAH20H
C     [4.1] F-FUNCTION
          F1=DSQRT(3.D0/8.D0)*((1.D0/(G(1)**QM)-1.D0)**(1.D0/QM))
          F2=DSQRT(3.D0/8.D0)*((1.D0/(G(2)**QM)-1.D0)**(1.D0/QM))
C     [4.1] FLUCTUATING COMPONENT
          PHI_H1= F1*DABS(DEV_SH-DABS(DEV_SH))
          PHI_H2= F2*DABS(DEV_SH+DABS(DEV_SH))
          PHI_H = PHI_H1**QM+PHI_H2**QM
C-----------------------------------------------------------------------
C     [5] HAH20E
      ELSEIF(HARD_PAR .EQ. 5.5D0) THEN ! HAH20E
C     [5.1] GRADIENT OF YIELD SURFACE AT DEV_H
          CALL DEVIATORIC_TRANS(DEV_H, DEV_H0, 2)
          CALL GRAD(1,1,STAT_VAR,DEV_H0,DFDH,DDFDDH)
          CALL DEVIATORIC_TENS(T)
          DFDH=MATMUL(T, DFDH)
          CALL DEVIATORIC_TRANS(DEV_EH0, DFDH, 1)
          CALL NORM_HAH(DEV_EH0, DEV_EH)

C     [5.2] F-FUNCTION
          F1=DSQRT(3.D0/8.D0)*((1.D0/(G(1)**QM)-1.D0)**(1.D0/QM))
          F2=DSQRT(3.D0/8.D0)*((1.D0/(G(2)**QM)-1.D0)**(1.D0/QM))
C     [5.3] FLUCTUATING COMPONENT
          CALL DOUBLE_DOT(DEV_EH, DEV_SIG, NDIM4, DEV_ES)
          CALL DOUBLE_DOT(DEV_EH, DEV_H, NDIM4, DEV_EHH)
          PHI_H1= F1*DABS((DEV_ES-DABS(DEV_ES))/DEV_EHH)
          PHI_H2= F2*DABS((DEV_ES+DABS(DEV_ES))/DEV_EHH)
          PHI_H = PHI_H1**QM+PHI_H2**QM
      END IF
C-----------------------------------------------------------------------
C     [6] EQUIVALENT STRESS
      SIG_BAR=(XI**QM+PHI_H)**(1.D0/QM)

      RETURN
      END SUBROUTINE HAH20
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HAH20(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES THE HARDENING STATE VARIABLES OF HAH20
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION G(NDIM5), DEV_H(NDIM4), DEV_HP(NDIM4), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM4), DEV_H0(NDIM4), DEV_HS(NDIM4),
     1          DEL_DEV_H(NDIM4), DEL_DEV_HP(NDIM4), DEL_G(NDIM5)
      DIMENSION STAT_VAR(NDIM6)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   DEFINE COEFFICIENT
      XK=H3; XK1=H3; XK2=H4; XK3=H5; XK4=H6; XK5=H7
      XC=H8; XKC=H9; XL=H10; XKL=H11

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=STAT_VAR(I)
      END DO
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM5+I)
          DEV_HP(I)=STAT_VAR(NDIM5+NDIM4+I)
          DEV_HS(I)=STAT_VAR(NDIM5+2*NDIM4+I)
      END DO
C     # DEFINE YIELD STRENGHTH AND FLOW STRESS
      CALL FLOW_STRESS(STAT_VAR, 0.D0, FLOW_SIG0, DHDE)
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)

C     [1] PREREQUIRED VARIABLES
C     1.1 THE OTHER COEFFICIENTS
C     1.1.1   MICROSTRUCTURE DEVIATOR
C      XKP= XK1/2.D0
      XKP= 15.D0
      XI_R= 8.D0
C     1.1.2   BAUSCHINGER RELATED EFFECTS
      XI_B =4.D0
      XI_BP=1.5D0
C     1.1.3   PERMANENT SOFTENING
      XKS =XK1
      XI_S =1.D0
C     1.1.4   CROSS-LOADIN CONTRACTION
      XKCP =XKC
      XI_CP=6.D0
C     1.1.5   LATENT HARDENING
      XI_LP =0.5D0

C     1.2 STRAIN-PATH CHANGE PARAMETER
C     1.2.1
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)

C     1.2.2
      CALL DOUBLE_DOT(DEV_H, DEV_H0, NDIM4, COS_X)
      IF(DABS(COS_X) .GT. 1.D0) COS_X=COS_X/DABS(COS_X)
      COS_X2=COS_X*COS_X
      IF(COS_X .LT. 0.D0) THEN
          XLM= -1.D0
      ELSE
          XLM= 1.D0
      ENDIF

C     1.2.3
      CALL DOUBLE_DOT(DEV_HP, DEV_H0, NDIM4, COS_XP)
      IF(DABS(COS_XP) .GT. 1.D0) COS_XP=COS_XP/DABS(COS_XP)
      COS_XP2=COS_XP*COS_XP
      IF(COS_XP .LT. 0.D0) THEN
          XLMP= -1.D0
      ELSE
          XLMP= 1.D0
      ENDIF

C     [2] MICROSTRUCTURE DEVIATOR EVOLUTION
      DO I=1, NDIM4
C     2.1 DEV_HP
          DEL_DEV_HP(I)= XLMP*XKP*(DEV_H0(I)-COS_XP*DEV_HP(I))

C     2.2 DEV_H
C     #   DEV_H EVOLVES ONLY WHEN ITSELF OR DEV_HP ARE CLOSE ENOUGH FROM
C         THE ACTIVE STRESS STATE DEV_SIG AS CONTROLLED BY THE COEFFICIENT
C         XI_R.
          DEL_DEV_H(I)= XLM*XK*((COS_XP2**XI_R)
     1                 +(COS_X2**XI_R))*(DEV_H0(I)-COS_X*DEV_H(I))
      END DO
      DEV_H=DEV_H+DEL_DEV_H*DEL_EQPLAS
      DEV_HP=DEV_HP+DEL_DEV_HP*DEL_EQPLAS
      CALL NORM_HAH(DEV_H, DEV_H)
      CALL NORM_HAH(DEV_HP, DEV_HP)
C-----------------------------------------------------------------------
C     [3] REVERSE LOADING VARIABLE EVOLUTION
      DEL_G=0.D0
C     3.1 G(1)= G-
      DEL_G(1)= ((1.D0+XLM)/2.D0)*(XK1*(1.D0-G(1))/(G(1)**XI_BP)
     1          *(1.D0-(COS_X2**XI_B))
     2          +XK2*(XK3*(FLOW_SIG0/FLOW_SIG)-G(1))*(COS_X2**XI_B))
     3          +((1.D0-XLM)/2.D0)*(XK1*(1.D0-G(1))/(G(1)**XI_BP))

C     3.2 G(2)= G+
      DEL_G(2)= ((1.D0-XLM)/2.D0)*(XK1*(1.D0-G(2))/(G(2)**XI_BP)
     1          *(1.D0-(COS_X2**XI_B))
     2          +XK2*(XK3*(FLOW_SIG0/FLOW_SIG)-G(2))*(COS_X2**XI_B))
     3          +((1.D0+XLM)/2.D0)*(XK1*(1.D0-G(2))/(G(2)**XI_BP))

C     [4] PERMANENT SOFTENING EVOLUTION
      IF(G(4)*XLM .LT. 0.D0) THEN
          G(4) =XLM
          G(9) =G(3)  ! G3*
          G(10)=G(5)  ! GP*
          DEV_HS=DEV_H! H*
          CALL DOUBLE_DOT(DEV_H0, DEV_HS, NDIM4, DOUBLE2) !DOUBLE2= S:H*
          IF(DABS(DOUBLE2).GT.1.D0) DOUBLE2=DOUBLE2/DABS(DOUBLE2)
          DOUBLE22=DOUBLE2*DOUBLE2
          GPP=G(10)-(G(10)-G(9))*DABS(DOUBLE22)**XI_S   !
          IF((GPP-G(5)) .LT. 0.D0) G(5)=GPP
          G(3)= G(5)
      ELSE
          CALL DOUBLE_DOT(DEV_H0, DEV_HS, NDIM4, DOUBLE2) !DOUBLE2= S:H*
          IF(DABS(DOUBLE2).GT.1.D0) DOUBLE2=DOUBLE2/DABS(DOUBLE2)
          DOUBLE22=DOUBLE2*DOUBLE2
          GPP=G(10)-(G(10)-G(9))*DABS(DOUBLE22)**XI_S   !
          IF((GPP-G(5)) .LT. 0.D0) G(5)=GPP
C     4.1 DEL_G3P= ¥ÄG3'
          DEL_G3P= XK5*G(6)*(XK4-G(3))
          G3P= G(3)+DEL_G3P*DEL_EQPLAS
C     4.2 DEL_G(6)= ¥ÄGS
          DEL_G(6)=-XKS*(1.D0-DABS(DOUBLE2))*G(6)
          IF((G3P-G(5)) .LE. 0.D0) THEN
              G(3)= G3P
          ELSE
              G(3)= G(5)
          END IF
          
      END IF

C     [5] CROSS-LOADING CONTRACTION EVOLUTION
      IF(DABS(COS_X) .GE. 0.996) THEN
          DEL_G(7)=XKC*(XC-G(7))
      ELSE
          DEL_G(7)=XKCP*(1.D0-G(7))/(G(7)**XI_CP)
      END IF

C     [6] LATENT HARDENING EVOLUTION
      VAL1= (XL-G(8))*(1.D0-DABS(COS_XP)**XI_LP)
      VAL2= (1.D0-G(8))*(1.D0+DABS(COS_XP)**XI_LP)
      DEL_G(8)=XKL*(VAL1+VAL2)

C     [7] UPDTAE G-VALUES
      G=G+DEL_G*DEL_EQPLAS
      VAL=XK3*(FLOW_SIG0/FLOW_SIG)
      IF(G(1) .GT. 1.D0) THEN
          G(1)=1.D0
      ELSEIF(G(1) .LT. VAL) THEN
          G(1)= VAL
      END IF
      IF(G(2) .GT. 1.D0) THEN
          G(2)=1.D0
      ELSEIF(G(2) .LT. VAL) THEN
          G(2)= VAL
      END IF
      IF(G(3) .GT. 1.D0) THEN
          G(3)=1.D0
      ELSEIF(G(3) .LT. XK4) THEN
          G(3)=XK4
      END IF

C     [9] UPDATE STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          STAT_VAR(I)=G(I)
      END DO
      DO I=1, NDIM4
          STAT_VAR(NDIM5+I)=DEV_H(I)
          STAT_VAR(NDIM5+NDIM4+I)=DEV_HP(I)
          STAT_VAR(NDIM5+2*NDIM4+I)=DEV_HS(I)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HAH20
C-----------------------------------------------------------------------
      SUBROUTINE HEXAH(STAT_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EFFECTIVE STRESS GENERATED BY HEXAH MODEL.
C     #   REF: B. REYNE, 2022, IJP
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM4), DEV_SH(MSS), DEV_PCI(NDIM4)
C     GLOBAL STATE VARIABLES
      DIMENSION PLAS(NDIM3)
C     LOCAL STATE VARIABLES
      DIMENSION G(MSS,NDIM5), DEV_H(MSS,NDIM4)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
      COMMON /KHEXAH/ MSS

C     #   DEFINE COEFFICIENT
      XQM=2.D0; XPM=2.D0;

C     [1] STATE VARIABLES OF HEXAH MODEL
C     [1.1]   GLOBAL STATE VARIABLES
C     #   GLOBAL G-VALUES
      GP=STAT_VAR(1)
      GL=STAT_VAR(2)
C     #   NUMBER OF ACTIVATED SLIP SYSTEM
      NSS=IDINT(STAT_VAR(3)) ! # OF SLIP SYSTEM
C     #   PLASTIC STRAIN
      DO I=1, NDIM3
          PLAS(I)=STAT_VAR(3+I)
      END DO
C     [1.2]   LOCAL STATE VARIABLES
      DO N=1, NSS
          NDIM=NDIM3+3+(N-1)*(NDIM5+NDIM4+2)
C     #   LOCAL G-VALUES
          DO I=1, NDIM5
              G(N,I)=STAT_VAR(NDIM+I)
          END DO
C     #   LOCAL MICROSTRUCTURE DEVIATOR
          DO I=1, NDIM4
              DEV_H(N,I)=STAT_VAR(NDIM+NDIM5+I)
          END DO
      END DO
C-----------------------------------------------------------------------
C     [2] STABLE COMPONENT
C     [2.1]   EFFECTIVE STRESS OF YIELD FUNCTION
C     #   DEVIATORIC STRESS
      CALL DEVIATORIC(SIG, DEV_SIG)
C     #   EFFECTIVE STRESS
      CALL YLD(SIG, DEV_SIG, PHI, 2)
C     [2.2] CROSS-LOADING CONTRACTION CONTRIBUTION
      PHI_C=0.D0
      DEV_SH=0.D0
      DO N=1, NSS
C     #   ORTHOGONAL STRESS COMPONENTS
          CALL DOUBLE_DOT(DEV_SIG, DEV_H(N,:), NDIM4, DEV_SH(N))
          DEV_PCI=4.D0*(1.D0-G(N,3))*(DEV_SIG-DEV_SH(N)*DEV_H(N,:))
C     #   EFFECTIVE STRESS
          CALL YLD(SIG, DEV_PCI, PHI_CI, 2)
          PHI_C=PHI_C+PHI_CI**XPM
      END DO
      XKSI=(PHI**XPM + PHI_C)**(1.D0/XPM)
C-----------------------------------------------------------------------
C     [3] FLUCTUATION(TRUNCATION) COMPONENT
      XF=0.D0
      DO N=1, NSS
C     [3.1] F-FUNCTION
          FG1=DSQRT(3.D0/2.D0)*((1.D0/(G(N,1)**XQM)-1.D0)**(1.D0/XQM))
     1        *(DABS(DEV_SH(N)-DABS(DEV_SH(N)))/2.D0)
          FG2=DSQRT(3.D0/2.D0)*((1.D0/(G(N,2)**XQM)-1.D0)**(1.D0/XQM))
     1        *(DABS(DEV_SH(N)+DABS(DEV_SH(N)))/2.D0)
C     [3.2] FLUCTUATING COMPONENT
          XF=XF+(FG1**XQM+FG2**XQM)
      END DO
C-----------------------------------------------------------------------
C     [4] EQUIVALENT STRESS
      SIG_BAR=(XKSI**XQM+XF)**(1.D0/XQM)
      SIG_BAR=SIG_BAR/(GP*GL)
      RETURN
      END SUBROUTINE HEXAH
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HEXAH0(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES THE HARDENING STATE VARIABLES OF HEXAH
C     #   A CRITERION FOR RECYCLE OF PRE-ACTIVATED SLIP SYSTEM
C     #   REF: B. REYNE ET AL., 2022
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM4), DEV_S(NDIM4), DFDS(NDIM3), 
     1          DDFDDS(NDIM3,NDIM3), DEL_G(NDIM5), DEL_DEV_H(NDIM4), 
     2          DPLAS(NDIM3), IDSS(MSS+1), CRIT(MSS), GM(MSS+1,NDIM5), 
     3          DEV_HM(MSS+1,NDIM4), EQPLASIM(MSS+1), COS_XM(MSS+1)
      DIMENSION DEV_AUX(NDIM4), DPLAS0(NDIM4), PLAS0(NDIM4)
C     GLOBAL STATE VARIABLES
      DIMENSION PLAS(NDIM3)
C     LOCAL STATE VARIABLES
      DIMENSION G(MSS+1,NDIM5), DEV_H(MSS+1,NDIM4), EQPLASI(MSS+1),
     1          COS_X(MSS+1)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
      COMMON /KHEXAH/ MSS
      PARAMETER(PI=DACOS(-1.D0))

C     #   DEFINE COEFFICIENT
      XTB=H1; XRB=H2; XBB=H3; XTP=H4; XRP=H5;
      XTC=H6; XRC=H7; XGC=H8; XTL=H9; XDEL=H10; XEL=H11;

C     [0] STATE VARIABLES OF HEXAH MODEL
C     [0.1]   GLOBAL STATE VARIABLES
C     #   GLOBAL G-VALUES
      GP=STAT_VAR(1)
      GL=STAT_VAR(2)
C     #   NUMBER OF ACTIVATED SLIP SYSTEM
      NSS=IDINT(STAT_VAR(3))
C     #   PLASTIC STRAIN
      DO I=1, NDIM3
          PLAS(I)=STAT_VAR(3+I)
      END DO
C     [0.2]   LOCAL STATE VARIABLES
      DO N=1, NSS
          NDIM=NDIM3+3+(N-1)*(NDIM5+NDIM4+2)
C     #   LOCAL G-VALUES
          DO I=1, NDIM5
              G(N,I)=STAT_VAR(NDIM+I)
          END DO
C     #   LOCAL MICROSTRUCTURE DEVIATOR
          DO I=1, NDIM4
              DEV_H(N,I)=STAT_VAR(NDIM+NDIM5+I)
          END DO
C     #   LOCAL EQUIVALENT PLASTIC STRAIN
          EQPLASI(N)=STAT_VAR(NDIM+NDIM5+NDIM4+1)
C     #   LOCAL STRAIN PATH CHANGE PARAMETER
          COS_X(N)=STAT_VAR(NDIM+NDIM5+NDIM4+2)
      END DO
C     [1] PREREQUISITES
C      XI_MAX=PI/3.D0
      XI_MAX=PI/6.D0
      XRH=100.D0 ! RATE OF ROTATION OF DEV_H
      COS_XMAX=DCOS(XI_MAX)
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_S)
      CALL DOUBLE_DOT(DEV_H(NSS,:), DEV_S, NDIM4, COS_X(NSS))
      IF(DABS(COS_X(NSS)).GT.1.D0)COS_X(NSS)=COS_X(NSS)/DABS(COS_X(NSS))
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
      DEL_G=0.D0
C-----------------------------------------------------------------------
C     [2.1]   SPAWN OR RECYCLE MICROSTRUCTURE DEVIATOR
C     [2.1.1] ACTIVATION OF NEW SLIP SYSTEM
      IF(DABS(COS_X(NSS)).LT.COS_XMAX) THEN
C     #   SPAWN NEW SLIP SYSTEM
          NSS=NSS+1
C     #   INITIALIZE LOCAL MICROSTRUCTURE DEVIATOR
          DEV_H(NSS,:)=DEV_S
C     #   INITIALIZE LOCAL G-VALUES
          G(NSS,:)=1.D0
C     #   STORE THE EFFECTIVE STRAIN
          EQPLASI(NSS)=EQPLAS
C     #   STRAIN PATH CHANGE PARAMETER
          COS_X(NSS)  =1.D0
      END IF
C-----------------------------------------------------------------------
C     [3] EVOLUTION LAWS OF LOCAL HEXAH VARIALBES
      DO N=1, NSS
      IF(N .EQ. NSS) THEN ! ACTIVE
C     [3.0A]   SIGN OF STRAIN PATH CHANGE PARAMETER: XLM
          IF(COS_X(N) .LT. 0.D0) THEN
              XLM= -1.D0
          ELSE
              XLM= 1.D0
          ENDIF
C     [3.1A]   MICROSTRUCTURE DEVIATOR
          DEL_DEV_H=XLM*XRH*(DEV_S-COS_X(N)*DEV_H(N,:))
C     [3.2A]   BAUSCHINGER EFFECT VARIABLES
C         [3.2A.1] G(1)= G-
          DEL_G(1)= XRB*((1.D0-XLM)/2.D0)*(1.D0-G(N,1))/(G(N,1)**XBB)
     1             +XRB*((1.D0+XLM)/2.D0)*(XTB-G(N,1))/(G(N,1)**XBB)
C         [3.2A.1] G(2)= G+
          DEL_G(2)= XRB*((1.D0+XLM)/2.D0)*(1.D0-G(N,2))/(G(N,2)**XBB)
     1             +XRB*((1.D0-XLM)/2.D0)*(XTB-G(N,2))/(G(N,2)**XBB)
C     [3.3A]   CROSS-LOADING CONTRACTION VARIABLES
          DEL_G(3)= XRC*(XTC-G(N,3))/(G(N,3)**XGC)
      ELSE ! INACTIVE
C     [3.1B]   MICROSTRUCTURE DEVIATOR
          DEL_DEV_H=0.D0
C     [3.2B]   BAUSCHINGER EFFECT VARIABLES
C         [3.2B.1] G(1)= G-
          DEL_G(1)= XRB*(1.D0-G(N,1))/(G(N,1)**XBB)
C         [3.2B.2] G(2)= G+
          DEL_G(2)= XRB*(1.D0-G(N,2))/(G(N,2)**XBB)
C     [3.3B]   CROSS-LOADING CONTRACTION VARIABLES
          DEL_G(3)= XRC*(1.D0-G(N,3))/(G(N,3)**XGC)
      END IF

C     [3.4]   UPDATE LOCAL STATE VARIABLES
      DEV_H(N,:)=DEV_H(N,:)+DEL_DEV_H*DEL_EQPLAS
      G(N,1)=G(N,1)+DEL_G(1)*DEL_EQPLAS ! G+
      G(N,2)=G(N,2)+DEL_G(2)*DEL_EQPLAS ! G-
      G(N,3)=G(N,3)+DEL_G(3)*DEL_EQPLAS ! GC

C     [3.5]   NUMERICAL CONSTRAINT
      IF(G(N,1) .GT. 1.D0) THEN
          G(N,1)=1.D0
      ELSEIF(G(N,1) .LT. XTB) THEN
          G(N,1)=XTB
      END IF
      IF(G(N,2) .GT. 1.D0) THEN
          G(N,2)=1.D0
      ELSEIF(G(N,2) .LT. XTB) THEN
          G(N,2)=XTB
      END IF
      IF(G(N,3) .GT. 1.D0) THEN
          G(N,3)=1.D0
      ELSEIF(G(N,3) .LT. XTC) THEN
          G(N,3)=XTC
      END IF
      END DO
C-----------------------------------------------------------------------
C     [4] EVOLUTION LAWS OF GLOBAL HEXAH VARIALBES
C     [4.1]   PERMANENT SOFTENING VARIABLE
      DPLAS=DEL_EQPLAS*DFDS
      CALL DEVIATORIC_TRANS(DEV_AUX, DPLAS, 1)
      CALL NORM_TENS(DEV_AUX,NDIM4,DPLAS0)
      CALL DEVIATORIC_TRANS(DEV_AUX, PLAS, 1)
      CALL NORM_TENS(DEV_AUX,NDIM4,PLAS0)
      CALL DOUBLE_DOT(-DPLAS0, PLAS0, NDIM4, DEP_EP)
      IF(DEP_EP .GT. 0.D0) THEN
C          DEL_GP_TEM= XRP*(XTP-GP)*DEP_EP
C          DEL_GP_LIM=-(DHDE/FLOW_SIG)*GP
C          DEL_GP=MAX(DEL_GP_TEM,DEL_GP_LIM)
          DEL_GP= XRP*(XTP-GP)*DEP_EP
      ELSE
          DEL_GP=0.D0
      END IF
      GP=GP+DEL_GP*DEL_EQPLAS
C     [4.2]   LATENT HARDENING VARIABLES
      IF(NSS .GT. 1) THEN
C         [2.6.1] LATENT HARDENING RATE
          XRL=DEXP(XEL)/(DEXP(XEL)-1.D0)
C         [2.6.2] AMOUNT OF LATENT HARDENING
C         (1) WEIGHT 1
          W_EI=(EQPLASI(NSS)-EQPLASI(NSS-1))/XDEL
C         (2) WEIGHT 2
          SIN_X=DSQRT(1.D0-COS_X(NSS-1)*COS_X(NSS-1))
          SIN_XMAX=DSIN(XI_MAX)
          VAL=SIN_X-SIN_XMAX
          BRAC_SINX=(VAL+DABS(VAL))/2.D0
          W_XI=BRAC_SINX/(1.D0-SIN_XMAX)
C         (2) L
          XL= W_EI*W_XI*(XTL-1.D0)*(XRL/(XRL-1.D0))**(XRL-1.D0)
C         [2.6.3] STATE VARIABLE UPDATE
          DEQPLAS=EQPLAS-EQPLASI(NSS)
          GL=1.D0+XL*XRL*DEXP(-XRL*DEQPLAS)*(DEXP(DEQPLAS)-1.D0)
      ELSE
          GL=1.D0
      END IF
C-----------------------------------------------------------------------
C     [9] STORE STATE VARIABLES FOR HAH MODEL
C     [9.1]   DEACTIVATION OF SLIP SYSTEM CRITERIA
C     [9.1.1] INEFFECTIVE SLIP SYSTEM
      CRIT=0.D0
      IDSS=1
      DO N=1, NSS
          CRIT(N)=(1.D0-G(N,1))**2+(1.D0-G(N,2))**2+(1.D0-G(N,3))**2
          IF(DSQRT(CRIT(N)).LT.1.D-3.AND.N.NE.NSS.AND.NSS.GT.1) THEN
              G(N,:)=1.D0
          END IF
      END DO
C     [9.1.2] MAXIMUM NUMBER OF SLIP SYSTEM
C     #   NOTE: KILL THE MOST INEFFECTIVE SLIP SYSTEM
      IF(NSS .GT. MSS) THEN
          CRIT0=MINVAL(CRIT(1:MSS))
          DO N=1,MSS
              IF(CRIT(N).EQ.CRIT0) IDSS(N)=0
          END DO
          NSS=MSS
      END IF
C     [9.1.3] PUSH FORWARD TO STATE VARIABLE MATRIX
      INDX=1
      GM=0.D0
      DEV_HM=0.D0
      DO N=1, NSS+1
          IF(IDSS(N) .EQ. 1) THEN ! EFFECTIVE(=ALIVE) (!= ACTIVATED)
C             #   G-VALUES MATRIX
              DO I=1, NDIM5
                  GM(INDX,I)=G(N,I)
              END DO
C             #   MICROSTRUCTURE DEVIATOR MATRIX
              DO I=1, NDIM4
                  DEV_HM(INDX,I)=DEV_H(N,I)
              END DO
C             #   EQUIVALENT STRAIN MATRIX
              EQPLASIM(INDX)=EQPLASI(N)
C             #   STRAIN PATH CHANGE PARAMETER MATRIX
              COS_XM(INDX)=COS_X(N)
              INDX=INDX+1
          END IF
      END DO
C     [9.2]   GLOBAL STATE VARIABLES
C     #   GLOBAL G-VALUES
      STAT_VAR(1)=GP
      STAT_VAR(2)=GL
C     #   NUMBER OF ACTIVATED SLIP SYSTEM
      STAT_VAR(3)=NSS
C     #   PLASTIC STRAIN
      DO I=1, NDIM3
          STAT_VAR(3+I)=PLAS(I)+DPLAS(I)
      END DO
C     [9.3]   LOCAL STATE VARIABLES
      DO N=1, NSS
          NDIM=NDIM3+3+(N-1)*(NDIM5+NDIM4+2)
C     #   LOCAL G-VALUES
          DO I=1, NDIM5
              STAT_VAR(NDIM+I)=GM(N,I)
          END DO
C     #   LOCAL MICROSTRUCTURE DEVIATOR
          DO I=1, NDIM4
              STAT_VAR(NDIM+NDIM5+I)=DEV_HM(N,I)
          END DO
C     #   LOCAL EQUIVALENT PLASTIC STRAIN
          STAT_VAR(NDIM+NDIM5+NDIM4+1)=EQPLASIM(N)
C     #   LOCAL STRAIN PATH CHANGE PARAMETER
          STAT_VAR(NDIM+NDIM5+NDIM4+2)=COS_XM(N)
      END DO
      RETURN
      END SUBROUTINE UPDATE_HEXAH0
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HEXAH1(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES THE HARDENING STATE VARIABLES OF HEXAH (MODIFIED)
C     #   A CRITERION FOR RECYCLE OF PRE-ACTIVATED SLIP SYSTEM
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM4), DEV_S(NDIM4), DFDS(NDIM3), 
     1          DDFDDS(NDIM3,NDIM3), DEL_G(NDIM5), DEL_DEV_H(NDIM4), 
     2          DPLAS(NDIM3), PLAS0(NDIM3), DPLAS0(NDIM3), IDSS(MSS+1),
     3          GM(MSS+1,NDIM5), DEV_HM(MSS+1,NDIM4), EQPLASIM(MSS+1)
C     GLOBAL STATE VARIABLES
      DIMENSION PLAS(NDIM3)
C     LOCAL STATE VARIABLES
      DIMENSION G(MSS+1,NDIM5), DEV_H(MSS+1,NDIM4), EQPLASI(MSS+1)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
      COMMON /KHEXAH/ MSS
      PARAMETER(PI=DACOS(-1.D0))

C     #   DEFINE COEFFICIENT
      XTB=H1; XRB=H2; XBB=H3; XTP=H4; XRP=H5;
      XTC=H6; XRC=H7; XGC=H8; XTL=H9; XDEL=H10; XEL=H11;

C     [0] STATE VARIABLES OF HEXAH MODEL
C     [0.1]   GLOBAL STATE VARIABLES
C     #   GLOBAL G-VALUES
      GP=STAT_VAR(1)
      GL=STAT_VAR(2)
C     #   NUMBER OF ACTIVATED SLIP SYSTEM
      NSS=STAT_VAR(3)
C     #   PLASTIC STRAIN
      PLAS(1:NDIM3)=STAT_VAR(4:NDIM3+3)
C     [0.2]   LOCAL STATE VARIABLES
      DO N=1, NSS
          NDIM=NDIM3+3+(N-1)*(NDIM5+NDIM4+1)
C     #   LOCAL G-VALUES
          DO I=1, NDIM5
              G(N,I)=STAT_VAR(NDIM+I)
          END DO
C     #   LOCAL MICROSTRUCTURE DEVIATOR
          DO I=1, NDIM4
              DEV_H(N,I)=STAT_VAR(NDIM+NDIM5+I)
          END DO
C     #   LOCAL EQUIVALENT PLASTIC STRAIN
          EQPLASI(N)=STAT_VAR(NDIM+NDIM5+NDIM4+I)
      END DO
C     [1] PREREQUISITES
C      XI_MAX=PI/3.D0
      XI_MAX=PI/6.D0
      XB=1.D0
      COS_XMAX=DCOS(XI_MAX)
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_S)
      CALL FLOW_STRESS(STAT_VAR, EQPLAS, FLOW_SIG)
      CALL GRAD(2,1,STAT_VAR,SIG,DFDS,DDFDDS)
      DEL_G=0.D0
C-----------------------------------------------------------------------
C     [2.1]   SPAWN OR RECYCLE MICROSTRUCTURE DEVIATOR
      DO N=1, NSS
          CALL DOUBLE_DOT(DEV_H(N,:), DEV_S, NDIM4, COS_X)
          IF(DABS(COS_X) .GT. 1.D0) COS_X=COS_X/DABS(COS_X)
C         [2.1.1] ACTIVATION OF NEW SLIP SYSTEM
          IF(DABS(COS_X).LT.COS_XMAX) THEN
              NSS=NSS+1
              IASS=NSS ! INDEX FOR ACTIVATED SLIP SYSTEM
              DEV_H(NSS,:)=DEV_S ! SPAWN
              EQPLASI(NSS)=EQPLAS
C         [2.1.2] REACTIVATION OF EXISTING SLIP SYSTEM
          ELSEIF(DABS(COS_X).GE.COS_XMAX) THEN
              IASS=I
          END IF
      END DO
C-----------------------------------------------------------------------
C     [3] EVOLUTION LAWS OF LOCAL HEXAH VARIALBES
      DO N=1, NSS
      IF(N .EQ. IASS) THEN ! ACTIVE
C     [3.1A]   MICROSTRUCTURE DEVIATOR
          CALL DOUBLE_DOT(DEV_H(N,:), DEV_S, NDIM4, COS_X)
          IF(DABS(COS_X) .GT. 1.D0) COS_X=COS_X/DABS(COS_X)
          IF(COS_X .LT. 0.D0) THEN
              XLM= -1.D0
          ELSE
              XLM= 1.D0
          ENDIF
          DEL_DEV_H=XLM*XRH*(DEV_S-COS_X*DEV_H(N,:))
C     [3.2A]   BAUSCHINGER EFFECT VARIABLES
C         [3.2A.1] G(1)= G-
          DEL_G(1)= XRB*((1.D0-XLM)/2.D0)*(1.D0-G(N,1))/(G(N,1)**XB)
     1            +XRB*((1.D0+XLM)/2.D0)*(XTB-G(N,1))/(G(N,1)**XB)
C         [3.2A.1] G(2)= G+
          DEL_G(2)= XRB*((1.D0+XLM)/2.D0)*(1.D0-G(N,2))/(G(N,2)**XB)
     1            +XRB*((1.D0-XLM)/2.D0)*(XTB-G(N,2))/(G(N,2)**XB)
C     [3.3A]   CROSS-LOADING CONTRACTION VARIABLES
          DEL_G(3)=XRC*(XTC-G(N,3))/(G(N,3)**XG)
      ELSE ! INACTIVE
C     [3.1B]   MICROSTRUCTURE DEVIATOR
          DEL_DEV_H=0.D0
C     [3.2B]   BAUSCHINGER EFFECT VARIABLES
C         [3.2B.1] G(1)= G-
          DEL_G(1)= XRB*(1.D0-G(N,1))/(G(N,1)**XB)
C         [3.2B.2] G(2)= G+
          DEL_G(2)= XRB*(1.D0-G(N,2))/(G(N,2)**XB)
C     [3.3B]   CROSS-LOADING CONTRACTION VARIABLES
          DEL_G(3)=XRC*(1.D0-G(N,3))/(G(N,3)**XGC)
      END IF

C     [3.4]   UPDATE LOCAL STATE VARIABLES
      DEV_H(N,:)=DEV_H(N,:)+DEL_DEV_H*DEL_EQPLAS
      G(N,1)=G(N,1)+DEL_G(1)*DEL_EQPLAS ! G+
      G(N,2)=G(N,2)+DEL_G(2)*DEL_EQPLAS ! G-
      G(N,3)=G(N,3)+DEL_G(3)*DEL_EQPLAS ! GC
C     [3.5]   NUMERICAL CONSTRAINT
      CALL NORM_HAH(DEV_H(N,:), DEV_S)
      DEV_H(N,:)=DEV_S
      IF(G(N,1) .GT. 1.D0) THEN
          G(N,1)=1.D0
      END IF
      IF(G(N,2) .GT. 1.D0) THEN
          G(N,2)=1.D0
      END IF
      END DO
C-----------------------------------------------------------------------
C     [4] EVOLUTION LAWS OF GLOBAL HEXAH VARIALBES
C     [4.1]   PERMANENT SOFTENING VARIABLE
      DPLAS=DEL_EQPLAS*DFDS
      CALL NORM_TENS(DPLAS,NDIM3,DPLAS0)
      CALL NORM_TENS(PLAS,NDIM3,PLAS0)
      CALL DOUBLE_DOT(-DPLAS0, PLAS0, NDIM3, DEP_EP)
      IF(DEP_EP .GT. 0.D0) THEN
          DEL_GP_TEM= XRP*(XTP-GP)*DEP_EP
          DEL_GP_LIM=-(DHDE/FLOW_SIG)*GP
          DEL_GP=MAX(DEL_GP_TEM,DEL_GP_LIM)
      ELSE
          DEL_GP=0.D0
      END IF
      GP=GP+DEL_GP*DEL_EQPLAS
C     [4.2]   LATENT HARDENING VARIABLES
      IF(NSS .GT. 1) THEN
C         [2.6.1] LATENT HARDENING RATE
          XRL=DEXP(XEL)/(DEXP(XEL)-1.D0)
C         [2.6.2] AMOUNT OF LATENT HARDENING
C         (1) WEIGHT 1
          W_EI=(EQPLASI(NSS)-EQPLASI(NSS-1))/XDEL
C         (2) WEIGHT 2
          SIN_X=DSQRT(1.D0-COS_X*COS_X)
          SIN_XMAX=DSIN(XI_MAX)
          VAL=SIN_X-SIN_XMAX
          BRAC_SINX=(VAL+DABS(VAL))/2.D0
          W_XI=BRAC_SINX/(1.D0-SIN_XMAX)
C         (2) L
          XL= W_EI*W_XI*(XTL-1.D0)*(XRL/(XRL-1.D0))**(XRL-1.D0)
C         [2.6.3] STATE VARIABLE UPDATE
          DEQPLAS=EQPLAS-EQPLASI(NSS)
          GL=1.D0+XL*XRL*DEXP(-XRL*DEQPLAS)*(DEXP(DEQPLAS)-1.D0)
      ELSE
          GL=1.D0
      END IF
C-----------------------------------------------------------------------
C     [9] STORE STATE VARIABLES FOR HAH MODEL
C     [9.1]   GLOBAL STATE VARIABLES
C     #   GLOBAL G-VALUES
      STAT_VAR(1)=GP
      STAT_VAR(2)=GL
C     #   NUMBER OF ACTIVATED SLIP SYSTEM
      STAT_VAR(3)=NSS
C     #   PLASTIC STRAIN
      STAT_VAR(4:NDIM3+3)=PLAS(1:NDIM3)
C     [9.2]   DEACTIVATION OF SLIP SYSTEM CRITERIA
C     [9.2.1] INEFFECTIVE SLIP SYSTEM
      VAL=0.D0
      IDSS=1
      DO N=1, NSS
          VAL=(1.D0-G(N,1))**2+(1.D0-G(N,2))**2+(1.D0-G(N,3))**2
          IF(DSQRT(VAL) .LT. 1.D-3) THEN
              IDSS(I)=0 ! FLAG FOR DEACTIVATED SLIP SYSTEM
              NSS=NSS-1
          END IF
      END DO
C     [9.2.1] MAXIMUM NUMBER OF SLIP SYSTEM
      IF(NSS .GT. MSS) THEN
          IDSS(1)=0 ! THE OLDEST SLIP SYSTEM IS DELETED.
          NSS=MSS
      END IF
C     [9.2.3] PUSH FORWARD TO STATE VARIABLE MATRIX
      INDX=1
      DO N=1, NSS
          IF(IDSS(N) .EQ. 1) THEN ! EFFECTIVE(=ALIVE) (!= ACTIVATED)
C             #   G-VALUES MATRIX
              DO I=1, NDIM3
                  GM(INDX,J)=G(N,I)
              END DO
C             #   MICROSTRUCTURE DEVIATOR MATRIX
              DO I=1, NDIM4
                  DEV_HM(INDX,J)=DEV_H(N,I)
              END DO
C             #   EQUIVALENT STRAIN MATRIX
              EQPLASIM(INDX)=EQPLASI(N)
              INDX=INDX+1
          END IF
      END DO
C     [9.3]   LOCAL STATE VARIABLES
      DO N=1, NSS
          NDIM=NDIM3+3+(N-1)*(NDIM5+NDIM4+1)
C     #   LOCAL G-VALUES
          DO I=1, NDIM5
              STAT_VAR(NDIM+J)=GM(N,I)
          END DO
C     #   LOCAL MICROSTRUCTURE DEVIATOR
          DO I=1, NDIM4
              STAT_VAR(NDIM+NDIM5+J)=DEV_HM(N,I)
          END DO
C     #   LOCAL EQUIVALENT PLASTIC STRAIN
          STAT_VAR(NDIM+NDIM5+NDIM4+1)=EQPLASIM(N)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HEXAH1
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_CHABOCHE(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES STATE VARIABLES OF KINEMATIC HARDENING.
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      DIMENSION ALPHA(NDIM3), DEL_ALPHA(NDIM3), DDFDEDS(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

      DO I=1, NDIM3
          ALPHA(I)= STAT_VAR(I)
      END DO
      IF(H4 .EQ. 0.D0) THEN   ! CHABOCHE1[9]
          HP1= H1
          HP2= H2
          CALL CB_HARD(EQPLAS, FLOW_SIG)
      ELSE                    ! CHABOCHE2[10]
          CALL CB_HARD(EQPLAS, FLOW_SIG)
          CALL ISO_HARD(STAT_VAR, EQPLAS, FLOW_ISO, DHDE)
          CALL KINEMATIC_COMPONENT(STAT_VAR,EQPLAS,SIG,DFDE,DDFDEDS)
          AP=FLOW_ISO-FLOW_SIG
          HP1= H1*EQPLAS+ H2
          HP2= (HP1-DFDE)/AP
      END IF
      DEL_ALPHA= HP1*SIG/FLOW_SIG-HP2*ALPHA
      DO I=1, NDIM3
          STAT_VAR(I)=ALPHA(I)+DEL_ALPHA(I)*DEL_EQPLAS
      END DO
      SIG= SIG-DEL_ALPHA*DEL_EQPLAS
      RETURN
      END SUBROUTINE UPDATE_CHABOCHE
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_YU(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES BACK-STRESS VALUES OF YOSHIDA-UEMORI.
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      DIMENSION ALPHA(NDIM3), DEL_ALPHA(NDIM3), ALPHAS(NDIM3), 
     1          DEL_ALPHAS(NDIM3), BETA(NDIM3), DEL_BETA(NDIM3), 
     2          Q(NDIM3), DEL_Q(NDIM3)
      DIMENSION DFDS(NDIM3), DDFDDS(NDIM3,NDIM3), DEV_ALPHAS(NDIM4),
     1          DEV_SIG(NDIM4), DEV_BETAQ(NDIM4), BETAQ(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KAHARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   ISOTORPIC-KINEMATIC HARDENING
C     [1] STATE VARIALBES
      DO I=1, NDIM3
          ALPHA(I)= STAT_VAR(I)
          BETA(I) = STAT_VAR(NDIM3+I)
          Q(I)    = STAT_VAR(2*NDIM3+I)
      END DO
      R= STAT_VAR(3*NDIM3+1) !BOUNDING SURFACE

C     [2] STATE VARIALBE INCREMENTS
C     #   PREREQUISITES
      ALPHAS=ALPHA-BETA
      A=DABS(H1+R-H2) ! [EQ. 7b][10]: A=B+R-Y
      CALL YLD(ALPHAS, DEV_ALPHAS, ALPHAS_BAR, 1) ! [EQ. 7b][10]: ALPHAS_BAR
      CALL YLD(SIG+ALPHA, DEV_SIG, SIG_BAR, 1)
      CALL GRAD(1,1,STAT_VAR,SIG+ALPHA,DFDS,DDFDDS)
      IF(ALPHAS_BAR .EQ. 0.D0) THEN
          DEL_ALPHAS=H3*(A/H2*SIG)
      ELSE
          DEL_ALPHAS= H3*((A/H2)*SIG-DSQRT(A/ALPHAS_BAR)
     1                      *ALPHAS)
      END IF
      DEL_BETA=H4*((2.D0/3.D0)*H5*DFDS-BETA)
      DEL_ALPHA= DEL_ALPHAS+DEL_BETA

C     #   ISOTROPIC HARDENING FOR BOUNDING SURFACE
      IF(FLOW_PAR .EQ. 0.D0) THEN
          DEL_R=H4*(H7-R)
      ELSE
          CALL YU_HARD(EQPLAS, R)
          DEL_R=0.D0
      END IF

C     #   WORK-HARDENING STAGNATION
C     [1] STAGNATION SURFACE: VON-MISES YIELD FUNCTION
      BETAQ=BETA-Q
      CALL YLD(BETAQ, DEV_BETAQ, STAG_BAR, 1)
C     [2] STATE VARIALBE INCREMENTS
      CALL DOUBLE_DOT(BETAQ,DEL_BETA*DEL_EQPLAS,NDIM3,STAG_DOUBLE)! [EQ. 23a][10]
      IF(STAG_DOUBLE.GT.0.D0) THEN
          ! #ASSUMPTION: STAG_R= STAG_BAR
          STAG_GAMMA=(3.D0/2.D0)*STAG_DOUBLE/STAG_BAR
          STAG_MU=(1.D0-H6)*STAG_GAMMA/STAG_BAR         ! [EQ. 27][10]
      ELSE
C         DEL_R=0.D0
          STAG_GAMMA=0.D0
          STAG_MU=0.D0
      END IF
      DEL_Q=STAG_MU*(BETAQ)
c     [3] UPDATE STATE VARIALBES
      DO I=1, NDIM3
          STAT_VAR(I)         = ALPHA(I)+DEL_ALPHA(I)*DEL_EQPLAS
          STAT_VAR(NDIM3+I)   = BETA(I)+DEL_BETA(I)*DEL_EQPLAS
          STAT_VAR(2*NDIM3+I) = Q(I)+DEL_Q(I)*DEL_EQPLAS
          SIG(I)= SIG(I)-DEL_ALPHA(I)*DEL_EQPLAS
      END DO
      STAT_VAR(3*NDIM3+1)= R + DEL_R*DEL_EQPLAS

      RETURN
      END SUBROUTINE UPDATE_YU
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_RGBV(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
C     # THIS SUBROUTINE FOLLOWS THE DISLOCATION EVOLUTION LAWS BY K. KITAYAMA'S WORK [16].
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION STAT_VAR(NDIM6),SIG(NDIM3),DEV_SIG(NDIM4),H(3),PSIG(3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KIHARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10

C     # DETERMINE STRAIN-PATH CHANGE PARAMETER
      IF(HARD_PAR.GE.3.) THEN ! HAH MODEL
          CALL COSX(STAT_VAR,SIG,COS_X)
      ELSE
          CALL PRIN_YLD2004_18P(SIG, H, PSIG)
          VAL=SUM(PSIG)
          IF(VAL .GE. 0.D0) THEN
              COS_X=1.D0
          ELSE
              COS_X=-1.D0
          END IF
      END IF
C     # READ STATE VARIABLE
      RHO_F=STAT_VAR(NDIM7-NDIM8+1)
      RHO_R1=STAT_VAR(NDIM7-NDIM8+2)
      RHO_R2=STAT_VAR(NDIM7-NDIM8+3)
      RHO_L=STAT_VAR(NDIM7-NDIM8+4)
      RHO_F01=STAT_VAR(NDIM7-NDIM8+5)
      RHO_F02=STAT_VAR(NDIM7-NDIM8+6)
      RHO_TOT=STAT_VAR(NDIM7-NDIM8+7)
      RHO_RTOT=STAT_VAR(NDIM7-NDIM8+9)
C     # VARIOUS REVERSIBLE PARAMETER
      VAL1=(RHO_MAX+RHO_MIN)/(RHO_MAX-RHO_MIN)
      VAL2=(RHO_MAX+RHO_MIN)/2.D0
      VAL3=(DSQRT(RHO_TOT)-VAL2)/VAL2
      RP=0.5D0-0.5D0*DTANH(3.D0*VAL1*VAL3)
C     # UPDATE DISLOCATION DENSITY
      XLAM=DSQRT(RHO_TOT)/(P4*P5)+1.D0/(P4*P6) ! MEAN FREE-PATH
      XLAM=P1*XLAM
      DEL_RHOF=(1.D0-P8)*XLAM-P7*RHO_F1
      XM=1.D0
      IF(COS_X .GE. 0.D0) THEN
          DEL_RHOR1=P8*XLAM-P7*RHO_R1
          DEL_RHOR2=-XLAM*(DABS(RHO_R2/RHO_F02))**XM
          RHO_F01=RHO_TOT
          XX=(RHO_R1/RHO_TOT)**0.5D0
      ELSE
          DEL_RHOR1=-XLAM*(DABS(RHO_R1/RHO_F01))**XM
          DEL_RHOR2=P8*XLAM-P7*RHO_R2
          RHO_F02=RHO_TOT
          XX=(RHO_R2/RHO_TOT)**0.5D0
      END IF
      RHO_F=RHO_F+DEL_RHOF*DABS(DEL_EQPLAS)
      RHO_R1=RHO_R1+DEL_RHOR1*DABS(DEL_EQPLAS)
      RHO_R2=RHO_R2+DEL_RHOR2*DABS(DEL_EQPLAS)
      RHO_TOT=RHO_F+RHO_R1+RHO_R2
      RHO_RTOT=RHO_R1+RHO_R2
      DRDE=DEL_RHOF+DEL_RHOR1+DEL_RHOR2
C     # UPDATE BACK-STRESS
C      XX=STAT_VAR(NDIM7-NDIM8+9)
C      TAU=P2+P3*P4*P5*DSQRT(RHO_TOT)
C      XXS=P10*(TAU-P2)
C      DEL_XX=P11*(XXS-XX)
C      XX=XX+DEL_XX*DABS(DEL_EQPLAS)
C     # UPDATE STATE VARIABLE
      STAT_VAR(NDIM7-NDIM8+1)=RHO_F
      STAT_VAR(NDIM7-NDIM8+2)=RHO_R1
      STAT_VAR(NDIM7-NDIM8+3)=RHO_R2
      STAT_VAR(NDIM7-NDIM8+4)=RHO_L
      STAT_VAR(NDIM7-NDIM8+5)=RHO_F01
      STAT_VAR(NDIM7-NDIM8+6)=RHO_F02
      STAT_VAR(NDIM7-NDIM8+7)=RHO_TOT
      STAT_VAR(NDIM7-NDIM8+8)=DRDE
      STAT_VAR(NDIM7-NDIM8+9)=RHO_RTOT
      RETURN
      END SUBROUTINE UPDATE_RGBV
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HARD_VAR(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION STAT_VAR(NDIM6), SIG(NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     # UPDATE ANISOTROPIC HARDENING MODEL
      IF(HARD_PAR .EQ. 1.D0) THEN ! CHABOCHE
          CALL UPDATE_CHABOCHE(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      ELSEIF(HARD_PAR .EQ. 2.D0) THEN ! YOSHIDA-UEMORI
          CALL UPDATE_YU(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      ELSEIF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
          CALL UPDATE_HAH11(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
          CALL UPDATE_HAH14(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
          CALL UPDATE_HAH20(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
          CALL UPDATE_HEXAH0(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      END IF
      
C     # UPDATE DISLOCATION HARDENING MODEL
      IF(FLOW_PAR.EQ.7.D0) THEN
          CALL UPDATE_RGBV(STAT_VAR, EQPLAS, SIG, DEL_EQPLAS)
      END IF

      RETURN
      END SUBROUTINE UPDATE_HARD_VAR
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.6 ANISOTROPIC YIELD FUNCTIONS                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE VON_MISES(DEV_SIG, PHI)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM4)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

      CALL DOUBLE_DOT(DEV_SIG, DEV_SIG, NDIM4, PSI)
      PHI= DSQRT((3.D0/2.D0)*PSI)

      RETURN
      END SUBROUTINE VON_MISES
C-----------------------------------------------------------------------
      SUBROUTINE HILL_1948(DEV_SIG, PHI)
C     #   THIS SUBROUTINE CALCULATES THE EQUIVALENT STRESS OF GENERALIZED
C         HILL YIELD FUNCTION UNDER PLANE STRESS ASSUMPTION.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DEV_SIG(NDIM4), S(6), HP(3), PS(3)

      COMMON /KYLD1/ YM, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     #   DEFINE COEFFICIENT
      HE=YM; HF=Y1; HG=Y2; HH=Y3; HN=Y4; HL=Y5; HM=Y6

C     [1] STRESS TENSOR DIMENSION RESET
C     S(6)=[S11 S22 S33 S12 S23 S13]
      S=0.D0
      DO I=1, NDIM4
          S(I)= DEV_SIG(I)
      END DO

C     [2] QUADRATIC HILL'S YIELD FUNCTION
      PSI= HF*(S(2)-S(3))**2+HG*(S(3)-S(1))**2+HH*(S(1)-S(2))**2
     1    +2.D0*HN*S(4)**2+2.D0*HL*S(5)**2+2.D0*HM*S(6)**2
      PHI=DSQRT(PSI)

      RETURN
      END SUBROUTINE HILL_1948
C-----------------------------------------------------------------------
      SUBROUTINE YLD2000_2D(DEV_SIG, PHI)
C     #   THIS SUBROUTINE RETURNS THE EQUIVALENT STRESS OF YLD2000_2D MODEL
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION DEV_SIG(NDIM3), X1P(3), X2P(3), PX1P(2), PX2P(2), ST(3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KYLD1/ YM, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     #   DEFINE COEFFICIENT
      AM=YM

C     [1] LINEAR TRANSFORMATION
      CALL LTM_YLD2000_2D(DEV_SIG, X1P, X2P)
C     [2] PRINCIPAL STRESS
      CALL PRIN_YLD2000_2D(X1P, PX1P)
      CALL PRIN_YLD2000_2D(X2P, PX2P)

C     [3] ANISOTROPIC YILED FUNCTION [EQ. 10 & EQ. 15][1]
      ST(1)=DABS(PX1P(1)-PX1P(2))            !-> |X'_1 - X'_2|
      ST(2)=DABS(2.D0*PX2P(2)+PX2P(1))        !-> |2*X''_2 + X''_1|
      ST(3)=DABS(2.D0*PX2P(1)+PX2P(2))        !-> |2*X''_1 + X''_2|
C
      PHI1P=ST(1)**AM
      PHI2P=ST(2)**AM + ST(3)**AM

      PSI= PHI1P + PHI2P
      PHI=((1.D0/2.D0)*(PSI))**(1.D0/AM)

      RETURN
      END SUBROUTINE YLD2000_2D
C-----------------------------------------------------------------------
      SUBROUTINE LTM_YLD2000_2D(DEV_SIG, X1P, X2P)
C     #   THIS SUBROUTINE RETURNS THE LINEAR TRANSFORMATION MATIRX REQUIRED
C         IN YLD2000_2D MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DEV_SIG(3), C1P(3,3), C2P(3,3), X1P(3), X2P(3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KYLD1/ YM, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9

C     #   DEFINE COEFFICIENT
      AM=YM; A1=Y1; A2=Y2; A3=Y3; A4=Y4; A5=Y5; A6=Y6; A7=Y7; A8=Y8;

      C1P= 0.D0
      C1P(1,1)= A1
      C1P(2,2)= A2
      C1P(3,3)= A7

      C2P= 0.D0
      C2P(1,1)= (4.D0/3.D0)*A5 - (1.D0/3.D0)*A3
      C2P(1,2)= (2.D0/3.D0)*(A6-A4)
      C2P(2,1)= (2.D0/3.D0)*(A3-A5)
      C2P(2,2)= (4.D0/3.D0)*A4 - (1.D0/3.D0)*A6
      C2P(3,3)= A8

C     LINEAR TRANSFORMATION [EQ. 14][1]
      X1P=MATMUL(C1P, DEV_SIG)
      X2P=MATMUL(C2P, DEV_SIG)

      RETURN
      END SUBROUTINE LTM_YLD2000_2D
C-----------------------------------------------------------------------
      SUBROUTINE PRIN_YLD2000_2D(XP, PXP)
C     # THIS SUBROUTINE RETURNS PRINCIPAL VALUES OF STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION XP(3), PXP(2)

C     THE PRINCIPAL VALUES OF X' & X'' [EQ. 17][1]
      XRP=DSQRT(((XP(1)-XP(2))/2.D0)**2.D0 + XP(3)**2.D0)    ! -> SQRT(¥Ä'/4)
C
      PXP(1)= (XP(1)+XP(2))/2.D0+XRP     !-> X'_1
      PXP(2)= (XP(1)+XP(2))/2.D0-XRP     !-> X'_2


      RETURN
      END SUBROUTINE PRIN_YLD2000_2D
C-----------------------------------------------------------------------
      SUBROUTINE YLD2004_18P(DEV_SIG, PHI)
C     # THIS SUBROUTINES RETURNS THE EQUIVALENT STRESS OF YLD2004-18P MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DEV_SIG(NDIM4), S1P(NDIM3), S2P(NDIM3),
     1          H1P(3), H2P(3), PS1P(3), PS2P(3)

      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KYLD1/ YM, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9
      COMMON /KYLD2/ Y10, Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18

C     #   DEFINE COEFFICIENT
      CM=YM

C     [1] LINEAR TRANSFORMATION
      CALL LTM_YLD2004_18P(DEV_SIG, S1P, S2P)

C     [2] PRINCIPAL STRESS
      CALL PRIN_YLD2004_18P(S1P, H1P, PS1P)
      CALL PRIN_YLD2004_18P(S2P, H2P, PS2P)

C     [3] ANISOTROPIC YIELD FUNCTION
      PSI=0.D0
      DO 1100 IP=1, 3
      DO 1100 IQ=1, 3
          PSI=PSI+DABS(PS1P(IP)-PS2P(IQ))**CM
1100  CONTINUE

      PHI= (PSI/4.D0)**(1.D0/CM)
      RETURN
      END SUBROUTINE YLD2004_18P
C-----------------------------------------------------------------------
      SUBROUTINE LTM_YLD2004_18P(DEV_SIG, S1P, S2P)
C     # THIS SUBROUTINE RETURNS LINEAR-TRANSFORMED STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DEV_SIG(NDIM4), S1P(NDIM3), S2P(NDIM3)
      DIMENSION C1P(NDIM3, NDIM3), C2P(NDIM3, NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      COMMON /KYLD1/ YM, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9
      COMMON /KYLD2/ Y10, Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18

C     #   DEFINE COEFFICIENT
      C1=Y1; C2=Y2; C3=Y3; C4=Y4; C5=Y5; C6=Y6; C7=Y7; C8=Y8;
      C9=Y9; C10=Y10; C11=Y11;C12=Y12; C13=Y13; C14=Y14; C15=Y15;
      C16=Y16; C17=Y17; C18=Y18;

      C1P=0.D0
      C1P(1,2)=-C1
      C1P(1,3)=-C2
      C1P(2,1)=-C3
      C1P(2,3)=-C4
      C1P(3,1)=-C5
      C1P(3,2)=-C6
      C1P(4,4)= C9
      C1P(5,5)= C7
      C1P(6,6)= C8

      C2P=0.D0
      C2P(1,2)=-C10
      C2P(1,3)=-C11
      C2P(2,1)=-C12
      C2P(2,3)=-C13
      C2P(3,1)=-C14
      C2P(3,2)=-C15
      C2P(4,4)= C18
      C2P(5,5)= C16
      C2P(6,6)= C17

C     #   LINEAR TRANSFORMATION
      S1P=MATMUL(C1P,DEV_SIG)
      S2P=MATMUL(C2P,DEV_SIG)

      RETURN
      END SUBROUTINE LTM_YLD2004_18P
C-----------------------------------------------------------------------
      SUBROUTINE PRIN_YLD2004_18P(SIG, H, PS)
C     #   THIS SUBROUTINE RETURNS THE PRINCIPAL STRESSES USED IN YLD2004 MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), S(6), PS(3), H(3)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      PARAMETER(PI=DACOS(-1.D0))

C     [0] STRESS TENSOR
C     S(6)=[S11 S22 S33 S12 S23 S13]
      S=0.D0
      IF(NDIM3 .EQ. 3) THEN
          S(1)=SIG(1)
          S(2)=SIG(2)
          S(4)=SIG(3)
      ELSE
          S=SIG
      END IF
C     [1] CALCULATE INVARIANTS
      H(1)=(S(1)+S(2)+S(3))/3.D0
      H(2)=(S(4)**2+S(5)**2+S(6)**2-S(1)*S(2)-S(1)*S(3)-S(2)*S(3))/3.D0
      H(3)=(2.D0*S(4)*S(5)*S(6)+S(1)*S(2)*S(3)-S(1)*S(5)**2-S(2)*S(6)**2
     1     -S(3)*S(4)**2)/2.D0

C     [2] PARAMETERS FROM INVARIANTS
      P=H(1)**2+H(2)
      Q=(2.D0*H(1)**3+3.D0*H(1)*H(2)+2.D0*H(3))/2.D0

      IF(DABS(Q/(P**(3.D0/2.D0))) .GT. 1.D0) THEN
          TH=0.D0
      ELSE
          TH=DACOS(Q/(P**(3.D0/2.D0)))
      END IF

C     [3] PRINCIPAL STRESSES
      PS(1)=2.D0*DSQRT(H(1)**2+H(2))*DCOS(TH/3.D0)+H(1)
      PS(2)=2.D0*DSQRT(H(1)**2+H(2))*DCOS((TH+4.D0*PI)/3.D0)+H(1)
      PS(3)=2.D0*DSQRT(H(1)**2+H(2))*DCOS((TH+2.D0*PI)/3.D0)+H(1)

      RETURN
      END SUBROUTINE PRIN_YLD2004_18P
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.7  NUMERICAL DIFFERENTIATION                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------------------------------------------
      SUBROUTINE GRAD(IOPT1, IOPT2, STAT_VAR, SIG, DFDS, DDFDDS)
C     # THIS SUBROUTINE CALCULATES PARTIAL DERIVATIVE THROUGH NUMERICAL
C       DIFFERENTIATION WHICH IS FINITE DIFFERENCE METHOD.
C       THREE-POINT MIDPOINT RULE WAS USED.
C     # IOPT1=1: YIELD FUNCTION DERIVATIVE
C            =2: EQUIVALENT STRESS DERIVATIVE (HAH OR YLD)
C     # IOPT2=1: DFDS
C             2: DFDS + DDFDDS
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3),SIGP(NDIM3),DFDS(NDIM3),DDFDDS(NDIM3,NDIM3),
     1          ASIG(NDIM3,3,NDIM3),ASIG_BAR(NDIM3,3),HM(NDIM3,3,NDIM3),
     2          H(NDIM3),DGDS(3)
      DIMENSION STAT_VAR(NDIM6), DEV_ASIG(NDIM4)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     #   CHECK STRESS TENSOR
      XSIG=0.D0
      DO I=1,NDIM3
          XSIG=XSIG+SIG(I)*SIG(I)
      END DO
      IF(XSIG.EQ.0.D0) THEN
          DFDS=0.D0
          DDFDDS=0.D0
          RETURN
      END IF

C     SCALING FACTOR [8]
      CALL DOUBLE_DOT(SIG, SIG, NDIM3, XI)
      XI=1.D0/DSQRT(XI)
C     SCALED STRESS COMPONENTS
      SIGP=SIG*XI

C     # DEFINE THE STEP SIZE
      H=1.D-6

      HM= 0.D0
      DFDS=0.D0
      DDFDDS=0.D0
      DO I=1, NDIM3
          DO J=1, 3
C         HM: MULTIPLIEROF STEP-SIZE, H
              HM(I, J, I)= 2.D0-J*1.D0

C         INPUT MATRIX FOR EQUIVALENT STRESS
C         ASIG=   {[SIG(1)+H1     SIG(2)  SIG(3)] [SIG(1) SIG(2)+H2   SIG(3)] [SIG(1) SIG(2)  SIG(3)+H3   ]}
C                 {[SIG(1)        SIG(2)  SIG(3)] [SIG(1) SIG(2)      SIG(3)] [SIG(1) SIG(2)  SIG(3)      ]}
C                 {[SIG(1)-H1     SIG(2)  SIG(3)] [SIG(1) SIG(2)-H2   SIG(3)] [SIG(1) SIG(2)  SIG(3)-H3   ]}
              ASIG(I, J, :) = SIGP
              DO K=1, NDIM3
                  ASIG(I, J, K)= SIGP(K)+HM(I,J,K)*H(K)
              END DO
C         CALCULATE EQUIVALENT STRESS
              IF(IOPT1 .EQ. 1) THEN
                  CALL YLD(ASIG(I,J,:), DEV_ASIG, ASIG_BAR(I,J), 1)
              ELSEIF(IOPT1 .EQ. 2) THEN
                  CALL EFF_SIG(STAT_VAR, ASIG(I,J,:), ASIG_BAR(I,J))
              ELSEIF(IOPT1 .EQ. 3) THEN
                  CALL YIELD_CONDITION
     1                (STAT_VAR,0.D0,ASIG(I,J,:),ASIG_BAR(I,J))
              END IF
          END DO
C         CALCULATE THE FIRST-ORDER PARTIAL DERIVATIVE
          CALL FDM_CAL(H(I), ASIG_BAR(I,:), DFDS(I))
      END DO
C-----------------------------------------------------------------------
C         CALCULATE THE SECOND-ORDER PARTIAL DERIVATIVE
      IF(IOPT2 .EQ. 2) THEN
      INDX=0.D0
      DO I=1, NDIM3
          DO J=1, NDIM3
              IF(I .EQ. J) THEN
C         #   I=J: (1,1), (2,2), (3,3)
                  DDFDDS(I, J)=
     1             (ASIG_BAR(I,1)-2.D0*ASIG_BAR(I,2)
     2            +ASIG_BAR(I,3))/(H(I)*H(I))
              ELSEIF(I .LT. J) THEN
C         #   I<J: (1,2), (1,3), (2,3)
                  INDX=INDX+1.D0
                  CALL FDM_NUMER(STAT_VAR,INDX,I,J,H,HM,SIGP,DGDS)
                  CALL FDM_CAL(H(I), DGDS, DDFDDS(I, J))
                  DDFDDS(J, I)= DDFDDS(I, J)
              END IF
          END DO
      END DO
      END IF
C-----------------------------------------------------------------------
      DDFDDS=XI*DDFDDS

      RETURN
      END SUBROUTINE GRAD
C-----------------------------------------------------------------------
      SUBROUTINE FDM_NUMER(STAT_VAR, INDX, I, J, H, HM, SIG, DGDS)
C     #   THIS SUBROUTINE RETURNS G-FUNCTION REQUIRED TO CALCULATE
C         SECOND-ORDER DERIVATIVE.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3),ASIG(9,NDIM3),ASIG_BAR(9),H(NDIM3),
     1          HM(NDIM3,3,NDIM3),HM2ND(9,NDIM3),DGDS(3)
      DIMENSION STAT_VAR(NDIM6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     #   INITIALIZATION
      HM2ND=0.D0
      ASIG=0.D0
      ASIG_BAR=0.D0
      DGDS=0.D0

C     #   H-VALUE MATRIX
      DO J0=1, 3
          N0=3*J0-2
          N1=3*J0
          HM2ND(N0:N1,I)=HM(I,J0,I)
          HM2ND(N0:N1,J)=HM(J,:,J)
      END DO

      DO J1=1, 9
C     INPUT MATRIX FOR EQUIVALENT STRESS FOR THREE-POINT RULE
C     ASIG=   {[SIG(1)+H1  SIG(2)+H2  SIG(3)] [SIG(1)+H1  SIG(2)  SIG(3)+H3] [SIG(1)  SIG(2)+H2  SIG(3)+H3]} 1
C             {[SIG(1)+H1  SIG(2)     SIG(3)] [SIG(1)+H1  SIG(2)  SIG(3)   ] [SIG(1)  SIG(2)+H2  SIG(3)   ]} 2
C             {[SIG(1)+H1  SIG(2)-H2  SIG(3)] [SIG(1)+H1  SIG(2)  SIG(3)-H3] [SIG(1)  SIG(2)+H2  SIG(3)-H3]} 3
C             {[SIG(1)     SIG(2)+H2  SIG(3)] [SIG(1)     SIG(2)  SIG(3)+H3] [SIG(1)  SIG(2)     SIG(3)+H3]} 4
C             {[SIG(1)     SIG(2)     SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)   ] [SIG(1)  SIG(2)     SIG(3)   ]} 5
C             {[SIG(1)     SIG(2)-H2  SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)-H3] [SIG(1)  SIG(2)     SIG(3)-H3]} 6
C             {[SIG(1)-H1  SIG(2)+H2  SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)+H3] [SIG(1)  SIG(2)-H2  SIG(3)+H3]} 7
C             {[SIG(1)-H1  SIG(2)     SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)   ] [SIG(1)  SIG(2)-H2  SIG(3)   ]} 8
C             {[SIG(1)-H1  SIG(2)-H2  SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)-H3] [SIG(1)  SIG(2)-H2  SIG(3)-H3]} 9
          DO K1=1, NDIM3
              ASIG(J1,K1)= SIG(K1) + HM2ND(J1,K1)*H(K1)
          END DO

C     CALCULATE EQUIVALENT STRESS
          CALL EFF_SIG(STAT_VAR,ASIG(J1,:),ASIG_BAR(J1))
      END DO

C     CALCULATE G-FUNCTION FOR THE SECOND-DERIVATIVE
      DO N2=1, 3
          N3=3*N2-2
          N4=3*N2
          CALL FDM_CAL(H(J), ASIG_BAR(N3:N4), DGDS(N2))
      END DO

      RETURN
      END SUBROUTINE FDM_NUMER
C-----------------------------------------------------------------------
      SUBROUTINE FDM_CAL(H, F, DF)
C     # THIS SUBROUTINE CALCULATES THREE MID-POINT RULE OF FDM.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(3)

      DF=(F(1)-F(3))/(2.D0*H)

      RETURN
      END SUBROUTINE FDM_CAL
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.8   OPERATORS                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------------------------------------------
      SUBROUTINE EFF_STR(DPLAS, DFDS, D_EQPLAS)
C     #   THIS SUBROUTINE CALCULATES DIRECTLY THE EFFECTIVE PLASTIC STRAIN.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DPLAS(NDIM3), DFDS(NDIM3)
      DIMENSION DPLAS1(NDIM4), DFDS1(NDIM4)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
C     [1] EXTENSION OF THE TENSOR DIMENSION (ONLY FOR PLANE STRESS)
      CALL DEVIATORIC_TRANS(DPLAS1, DPLAS, 1) ! RETURN: DPLAS1
      CALL DEVIATORIC_TRANS(DFDS1, DFDS, 1) ! RETURN: DFDS1
C     [2] SCALAR PRODUCT OF THE VARIABLES
      CALL DOUBLE_DOT(DPLAS1, DPLAS1, NDIM4, DOUBLE_DPLAS)
      CALL DOUBLE_DOT(DFDS1, DFDS1, NDIM4, DOUBLE_DFDS)
C     [3] THE EFFECTIVE PLASTIC STRAIN INCREMENT
      D_EQPLAS=DSQRT(DOUBLE_DPLAS/DOUBLE_DFDS)
      RETURN
      END SUBROUTINE EFF_STR
C-----------------------------------------------------------------------
      SUBROUTINE DEVIATORIC(SIG, DEV_SIG)
C     #   THIS SUBROUTINE CALCULATE DEVIATORIC STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM4), DELTA(NDIM4)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

      DEV_SIG=0.D0
      NDI=NDIM4-NDIM2
C     INITIAL DEVIATORIC STRESS TENSOR
      IF(NDIM3 .EQ. 3) THEN ! PLANE STRESS
          DEV_SIG(1)= SIG(1)
          DEV_SIG(2)= SIG(2)
          DEV_SIG(3)= 0.D0
          DEV_SIG(4)= SIG(3)
      ELSE
          DEV_SIG= SIG
      END IF

      SIG_KK=0.D0
      DELTA=0.D0
      DO I=1, NDI
          SIG_KK=SIG_KK+DEV_SIG(I)/3.D0
          DELTA(I)=1.D0
      END DO
      DEV_SIG=DEV_SIG-SIG_KK*DELTA

      RETURN
      END SUBROUTINE DEVIATORIC
C-----------------------------------------------------------------------
      SUBROUTINE DEVIATORIC_TRANS(TENS, TENS0, IOPT)
C     #   IOPT=1: THIS SUBROUTINE REMOVES COMPONENT SIG_33.
C         -> RETURN: TENS0
C     #   IOPT=2: THE RECOVERY TO THE ORIGINAL DIMENSION OF THE DEVIATORIC TENSOR.
c         -> RETURN: TENS
C     #   NOTE: THIS SUBROUTINE IS NOT NECESSARY FOR THE 3D STRESS STATE.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION TENS(NDIM4), TENS0(NDIM3)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7
      
      IF(NDIM3 .EQ. 3) THEN
          IF(IOPT .EQ. 1) THEN
C             A. DIAGONAL COMPONENTS
              DO I=1, NDIM1
                  TENS(I)=TENS0(I)
              END DO
C             B. SHEAR COMPONENTS
              TENS(3)=-(TENS0(1)+TENS0(2))
              TENS(4)= TENS0(3)
          ELSEIF(IOPT .EQ. 2) THEN
C             A. DIAGONAL COMPONENTS
              TENS0(1)= TENS(1)
              TENS0(2)= TENS(2)
C             B. SHEAR COMPONENTS
              TENS0(3)= TENS(4)
          END IF
      ELSEIF(NDIM3 .EQ. 6) THEN
          IF(IOPT .EQ. 1) THEN
              TENS=TENS0
          ELSEIF(IOPT .EQ. 2) THEN
              TENS0=TENS
          END IF
      END IF

      RETURN
      END SUBROUTINE DEVIATORIC_TRANS
C-----------------------------------------------------------------------
      SUBROUTINE DEVIATORIC_TENS(T)
C     #   THIS SUBROUTINE REMOVES DEV_SIG33.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION T(NDIM3, NDIM3)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

      T=0.D0
      IF(NDIM3 .EQ. 3) THEN
          T(1,1)= 2.D0
          T(1,2)=-1.D0
          T(2,1)=-1.D0
          T(2,2)= 2.D0
          T(3,3)= 3.D0
      ELSEIF(NDIM3 .EQ. 6) THEN
          T(1,1)= 2.D0
          T(1,2)=-1.D0
          T(1,3)=-1.D0
          T(2,1)=-1.D0
          T(2,2)= 2.D0
          T(2,3)=-1.D0
          T(3,1)=-1.D0
          T(3,2)=-1.D0
          T(3,3)= 2.D0
          T(4,4)= 3.D0
          T(5,5)= 3.D0
          T(6,6)= 3.D0
      END IF
      T=T/3.D0

      RETURN
      END SUBROUTINE DEVIATORIC_TENS
C-----------------------------------------------------------------------
      SUBROUTINE NORM_TENS(TENS,NDIM,TENS0)
C     #   THIS SUBROUTINE NORMALIZES AN ARBITRARY TENSOR FOR GENERAL PURPOSES.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION TENS(NDIM), TENS0(NDIM)

C     #   DENOMINATOR
C     XDENOM = TENS_KL*TENS_KL
      CALL DOUBLE_DOT(TENS,TENS,NDIM,XDENOM)

      IF (XDENOM .EQ. 0.D0) THEN
          TENS0=0.D0
          RETURN
      END IF

C     #   NORMALIZED TENSOR: TENS0
      TENS0=TENS/DSQRT(XDENOM)

      RETURN
      END SUBROUTINE NORM_TENS
C-----------------------------------------------------------------------
      SUBROUTINE NORM_HAH(H,DEV_H)
C     #   THIS SUBROUTINE RETURNS THE NORMALIZED MICRO STRUCTURE DEVIATOR.
C         [EQ.(2)][3] WAS USED.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION H(NDIM4), DEV_H(NDIM4)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     #   NUMERATOR
C     XDENOM = H_KL*H_KL
      CALL DOUBLE_DOT(H,H,NDIM4,XDENOM)

      IF (XDENOM .EQ. 0.D0) THEN
          DEV_H=0.D0
          RETURN
      END IF
      HVAL=8.D0/3.D0
      IF(HARD_PAR .EQ. 3.D0) THEN !HAH11
C         XDENOM = SQRT[(8/3)*H_KL*H_KL]
          XDENOM= DSQRT(HVAL*XDENOM)
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN !HAH14
          XDENOM= DSQRT(HVAL*XDENOM)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN !HAH20
          XDENOM= DSQRT(XDENOM)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN !HEXAH
          XDENOM= DSQRT(XDENOM)
      END IF

C     #   NORMALIZED DEVIATOR
C     DEV_H = H^_IJ = H_IJ / SQRT[(8/3)*H_KL*H_KL]
      DEV_H=H/XDENOM

      RETURN
      END SUBROUTINE NORM_HAH
C-----------------------------------------------------------------------
      SUBROUTINE COSX(STAT_VAR, STRESS, COS_X)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION STAT_VAR(NDIM6), STRESS(NDIM3), DEV_SIG(NDIM4)
      DIMENSION DEV_S(NDIM4), DEV_H(NDIM4)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

      IF(HARD_PAR .EQ. 3.D0) THEN ! HAH11
          VAL=8.D0/3.D0
      ELSEIF(HARD_PAR .EQ. 4.D0) THEN ! HAH14
          VAL=8.D0/3.D0
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.D0) THEN ! HAH20
          VAL=1.D0
      ELSEIF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
          VAL=1.D0
      ELSE
          RETURN
      END IF

C     [1] NORMALIZED DEVIATORIC STRESS
      CALL DEVIATORIC(STRESS, DEV_SIG)
      CALL NORM_HAH(DEV_SIG,DEV_S)
C     [2] MICROSTRUCTURE DEVIATOR
      IF(FLOOR(HARD_PAR) .EQ. 6.D0) THEN ! HEXAH
          NSS=STAT_VAR(3)
          NDIM=NDIM3+3+(NDIM5+NDIM4+1)*(NSS-1)
      ELSE
          NDIM=0
      END IF
      DO I=1, NDIM4
          DEV_H(I)=STAT_VAR(NDIM+NDIM5+I)
      END DO
C     [3] STRAIN-PATH CHANGE PARAMETER
      CALL DOUBLE_DOT(DEV_S, DEV_H, NDIM4, COS_X)
      COS_X=VAL*COS_X
      IF(DABS(COS_X).GE.1.D0) COS_X=COS_X/DABS(COS_X)

      RETURN
      END SUBROUTINE COSX
C-----------------------------------------------------------------------
      SUBROUTINE DOUBLE_DOT(A,B,NDIM,DOUBLE)
C     #   A:B = A_IJ*B_IJ = TR(A^T*B)
C         A:B = A11*B11 + A12*B12 + A13*B13 + A21*B21 + A22*B22 + A23*B23 + A31*B31 + A32*B32 + A33*B33
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION A(NDIM), B(NDIM)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

C     DOUBLE = ¥Ò DEV_H(I)*DEV_SIG(I)
      DOUBLE = 0.D0
      NDI= NDIM-NDIM2
      DO I=1,NDI
          DOUBLE = DOUBLE + A(I)*B(I)
      END DO
      DO I=NDI+1, NDIM
          DOUBLE= DOUBLE + 2.D0*A(I)*B(I)
      END DO

      RETURN
      END SUBROUTINE DOUBLE_DOT
C-----------------------------------------------------------------------
      SUBROUTINE DOUBLE_DOT21(A,B,DOUBLE)
C     #   THIS SUBROUTINE CALCULATE DOUBLE DOT PRODUCT FOR NON-DEVIATORIC VARAIBLES
C         IN PLANE STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION A(NDIM3,NDIM3), B(NDIM3), DOUBLE(NDIM3)
      COMMON /KSIZE/ NDIM1, NDIM2, NDIM3, NDIM4, NDIM5, NDIM6, NDIM7

      DOUBLE = 0.D0
      NDI= NDIM3-NDIM2
      DO I=1, NDIM3
      DO J=1, NDIM3
          IF(J .LE. NDI) THEN
              DOUBLE(I) = DOUBLE(I) + A(I,J)*B(J)
          ELSE
              DOUBLE(I) = DOUBLE(I) + 2.D0*A(I,J)*B(J)
          END IF
      END DO
      END DO

      RETURN
      END SUBROUTINE DOUBLE_DOT21
C-----------------------------------------------------------------------
      SUBROUTINE DYAD(A, B, NDIM, AB)
c     #   THIS SUBROUTINE CALCULATES THE DYAD PRODUCT OF TWO TENSORS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION A(NDIM), B(NDIM), AB(NDIM,NDIM)

      DO I=1, NDIM
      DO J=1, NDIM
          AB(I,J)=A(I)*B(J)
      END DO
      END DO
      
      RETURN
      END SUBROUTINE DYAD
C-----------------------------------------------------------------------
      SUBROUTINE INVERSE(A, N, INV_A)
C     #   THIS SUBROUTINE RETURNS THE INVERSE MATRIX A(NXN) USING LU-DECOMPOSITION.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N), INDX(N)
      REAL*8 INV_A(N,N), LU_A(N,N)

      LU_A= A
      INV_A= 0.D0
      DO I=1,N
          INV_A(I,I)= 1.D0
      END DO
      CALL LUDCMP(LU_A, N, N, INDX, DDCMP)
      DO J=1,N
          CALL LUBKSB(LU_A, N, N, INDX, INV_A(1,J))
      END DO

      RETURN
      END SUBROUTINE INVERSE
C-----------------------------------------------------------------------
      SUBROUTINE LUDCMP (A, N, NP, INDX, D)
C-----  LU decomposition
C-----  Use single precision on cray
C-----  THIS CODE RETURN THE MIXED LU DECOMPOSED MATRIX
C-----  A = L*U = LU'
C-----  LU' = [U11 U12 U13]
C             [L21 U22 U23]
C             [L31 L32 L33]
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.D0
      DO I=1,N
         AAMAX=0.D0

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.D0) THEN
             WRITE(*, *) '#####################'
             WRITE(*, *) '## SINGULAR MATRIX ##'
             WRITE(*, *) '#####################'

             STOP
         END IF
         VV(I)=1.D0/AAMAX
      END DO
C     
      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO
      RETURN
      END SUBROUTINE LUDCMP
C-----------------------------------------------------------------------
      SUBROUTINE LUBKSB (A, N, NP, INDX, B)
C-----  Linear equation solver based on LU decomposition
C-----  Use single precision on cray
C----   IN THIS CODE, WE CAN GET THE INVERSE MATRIX OF A BY CALCULATING A MATRIX X
C----   A*X = I -> X = A^-1
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)
C     ## A = L*U
C     ## A*X = (L*U)*X = L*(U*X) = B
C     ## FORWARD SUBSTITUTION ##
C     ## L*Y = B: GET Y(I)!!
C     L(I,I) = 1
      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF
C     ## HERE: B(I) = SUM -> Y(I)
         B(I)=SUM 
      END DO
C     ## BACKWARD SUBSTITUTION ##
C     ## U*X = Y: GET X(I)!!
      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF
C     ## HERE: B(I) = SUM/A(I,I) -> X(I)
         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END SUBROUTINE LUBKSB
C-----------------------------------------------------------------------