# Information of developer
    1. Name: Seong-Yong Yoon
    2. Affiliation: Pohang university of science and technology (POSTECH)
    3. Advisor: Frederic Barlat
    4. E-mail: theysy@postech.ac.kr

# Brief description
    - This subroutine is developed to apply various macroscopic plasticty models into sheet metal FE simulations.
    - Plasticity theories are modulized and made easy to be appended within UMAT code.
    - Various stress update algorithms are available.

# Code name
    1. MML_U: UMAT subroutine
    2. MML_V: VUMAT subroutine

# Version information
    1. Version 1.0: Macroscopic plasticity models for plane stress
    2. Version 2.0: Ver.1 is extended to (2) 3D stress state
    3. Version 3.0: Pressure-dependent yield condition
    4. Version 2.F: Fracture criterion are coupled with Ver.2

# Implemented plasticity models
    1. Anisotropic hardening models
        - 0: Isotropic hardening
        - 1: Chaboche-type isotropic-kinematic hardening [9]
        - 2: Yoshida-Uemori kinematic hardening | 2.5: Modified versions [11]
        - 3: HAH11 [2]
        - 4: HAH14 [3]
        - 5: HAH20h | 5.5: HAH20e [14]
    2. Anisotropic yield functions
        - 1: Von-Mises
        - 2: Hill 1948
        - 3: Yld2000-2d [1]
        - 4: Yld2004-18p [7]
    3. Isotropic hardening laws 
        - 1: SWIFT: SIGMA=P1*(P2+EQPLAS)**P3
        - 2: VOCE: SIGMA=P1-P2*EXP(-P3*EQPLAS)
        - 3: MODIFIED VOCE: SIGMA= P1 + P2*EQPLAS + P3*(1-EXP(-P4*EQPLAS))
        - 4: SWIFT+MODIFIED VOCE: P1*(P2+EQPLAS)**P3 + P4 + P5*EQPLAS + P6*(1-EXP(-P7*EQPLAS))
        - 5: HOCKETT-SHERBY: P1-(P1-P2)*EXP(-P3*EQPLAS**P4)
        - 6: COMBINED SWIFT-VOCE: P1*[P2*(P3+EQPLAS)**P4]+(1-P1)*[P5-P6*EXP(-P7*EQPLAS)]
        - 7: MACROSCOPIC RGBV MODEL: P1*[P2+P3*P4*P5*SQRT(RHO(EQPLAS))
     4. Stress update algorithms
        - 1: CUTTING PLANE METHOD (EXPLICIT)
        - 2: CLOSEST POINT PROJECTION METHOD(FULLY-IMPLICIT)
        - 3: CLOSEST POINT PROJECTION METHOD(SEMI-IMPLICIT)
        - 4: TANGENT MODULUS METHOD
        - 5: LINE SEARCH EULER BACKWARD METHOD (FULLY-IMPLICIT)
     5. DIFFERENTIATION METHOD
        - 1: Analytical derivative
        - 2: Numerical derivative
# Reference
    [1] F.BARLAT  ET AL. IJP      (2003): YLD2000-2D
    [2] F.BARLAT  ET AL. IJP      (2012): HAH11
    [3] F.BARLAT  ET AL. IJP      (2014): HAH14
    [4] J.W. LEE  ET AL. IJP      (2012): CUTTING PLANE METHOD
    [5] J.W. LEE  ET AL. CMAM     (2012): CLOSEST POINT PROJECTION
    [6] J.W. YOON ET AL. IJP      (2004): MULTI-STAGE EULER BACKWARD
    [7] F. BARLAT ET AL.  IJP     (2005): YLD2004_18P
    [8] H. ARETZ          ESAFORM (2007): SCALED FDM
    [9] T.J. PARK ET AT.  IJSS    (2012): CHABOCHE1-SIMPLER MODEL
    [10]J.Y. LEE  ET AL.  IJSS    (2012): CHABOCHE2
    [11]F. YOSHIDA ET AL. IJP     (2002): YOSHIDA-UEMORI
    [12]M. ORTIZ  ET AL.  IJNME   (1985): TANGENT MODULUS METHOD
    [13]W.M. SCHERZINGER  CMAME   (2017): LINE SEARCH EULER BACKWARD
    [14]F. BARLAT ET AL.  IJSS    (2020): HAH20H & HAH20E
    [15]E. RAUCH ET AL.   MSMSE   (2011): RGBV
    [16]K. KITAYAMA ET AL.IJP     (2013): CRYSTALLOGRAPHIC RGBV
