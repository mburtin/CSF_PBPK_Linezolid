;; 1. Based on: run16
;; 2. Description: Minimal PBPK from Jusko - defix plasma - fit cerveau - fd fix
;; x1. Author: user




$PROBLEM    PBPK Critically ill patients metronidazole study

$INPUT      ID TIME TAD TIN CONC=DROP DV LNCO=DROP LNCU=DROP FU AMT II ADDL RATE EVID CMT VDVE QDVE 
			AGE SEX HT WT SALB=DROP SPRO=DROP CRCL DIAG=DROP STUD FLAG

$DATA       ../data/metronidazole_QDVEind.csv IGNORE=# 

$SUBROUTINE ADVAN13 TOL=9 

;$PRIOR      NWPRI NTHETA=12 NETA=6 NTHP=3 NETP=0 NPEXP=1

$MODEL  
NCOMP = 7   
COMP=(CENTRAL,DEFDOSE) 	; 1- Plasma
COMP=(RES) 				; 2- Rest
COMP=(BVAS) 			; 3- Brain vasculature
COMP=(ECF) 				; 4- Extracellular fluid
COMP=(MICRO)			; 5- Dialysate
COMP=(CSF)  			; 6- Cerebrospinal fluid
COMP=(DVE)  			; 7- DVE

$PK       

;; --------------------------- Parameters estimated ---------------------------- ;;
TVKRES   	= THETA(1)
KRES     	= TVKRES*EXP(ETA(1))

TVCL10   	= THETA(2)
CL10     	= TVCL10*EXP(ETA(2))

TVfd		= THETA(3)
fd			= TVfd*EXP(ETA(3))		;fraction of CO for the rest compartment

TVVBLO		= THETA(4)
VBLO		= TVVBLO*EXP(ETA(4))

TVCO		= THETA(5)
CO			= TVCO*EXP(ETA(5))

TVPS1     = THETA(6)
PS1       = TVPS1*EXP(ETA(6))

TVPS2     = THETA(7)
PS2       = TVPS2*EXP(ETA(7))


QCSFOUT = 0
IF (STUD.EQ.1) QCSFOUT = 0.024
IF (STUD.EQ.2.AND.QDVE.LT.0.024) QCSFOUT = 0.024-QDVE				;patients LCR



;; ------------------------- Physiological parameters -------------------------- ;;
; Volumes:
VBVAS = 0.0637                            	;; Brain vascular volume L (paper from Ball K)
;VBLO  = 5.85			                   	;; Blood volume L (Brown et al. Physiological parameter values for physiologically based pharmacokinetic models. Toxicol. Ind. Health 13, 407−84)
VECF  = 0.24			                   	;; ECF brain tissue volume L (Westerhout et al)
VVENT = 0.130                            	;; CSF volume L (Gaohua et al, Drug Metabolism and Pharmacokinetics Volume 31, Issue 3, June 2016, Pages 224-233)
VRES  = WT-(VBLO+VBVAS+VECF+VVENT)			;; Rest volume L

; Blood flows:
;CO   = 312			             			;; Cardiac output L/h (Brown et al.)
QBRA = 42				           			;; Brain blood flow L/h (paper from Ball K)
QRES = fd*CO                    			;; Rest blood flow
QECF = 0.0105 	                 			;; Bulk flow : flow rate from ECF to CSF (Westerhout et al)
RBP  = 0.8                         			;; blood/plasma ratio (simcyp predictions)


$DES       
;; ------------------------- Compartment concentrations ----------------------- ;;
C1 = A(1)/VBLO                                       ;; AMT/Blood  volume
Cplasma = C1/RBP									 ;; Conc plasma
C2 = A(2)/VRES                                       ;; AMT/Plasma rest tissue volume
C3 = A(3)/VBVAS                                      ;; AMT/Brain vascular volume
C4 = A(4)/VECF                                       ;; AMT/Brain ECF volume
C5 = A(5)/TIN										 ;; Conc dialysate
C6 = A(6)/VVENT                                      ;; AMT/CSF ventricle volume
C7 = A(7)/(VDVE+0.000001)							 ;; Conc dve
;; -------------------------- Differential equations -------------------------- ;;
DADT(1)= (C2*QRES/KRES) + C3*QBRA - C1*(QRES+CL10+QBRA)
DADT(2)= (C1*QRES) - (C2*QRES/KRES)

DADT(3)= (C1-C3)*QBRA + C4*PS1 + C6*(PS2+QCSFOUT) - C3*(PS1+PS2)

DADT(4)=  C3*PS1 - C4*(PS1+QECF)
DADT(5) = C4

DADT(6)=  C3*PS2 + (C4*QECF) - C6*(QCSFOUT+PS2) - C6*QDVE
DADT(7)=  C6*QDVE


$ERROR       

AA1 = A(1)			;; AMT/Plasma
AA2 = A(2)			;; AMT/Rest
AA3 = A(3)			;; AMT/Brain vasculature
AA4 = A(4)			;; AMT/ECFV
AA5 = A(5)			;; AMT/ECFI
AA6 = A(6)			;; AMT/CSF
AA7 = A(7)			;; AMT/DVE

	IPRED = A(1)/(VBLO*RBP)
		IF(CMT.EQ.5)    IPRED = (A(5)/TIN)
		IF(CMT.EQ.7)    IPRED = (A(7)/(VDVE+0.000001))

	W = SQRT(SIGMA(1)*IPRED**2 + SIGMA(2))
		IF(CMT.EQ.5)    W = SQRT(SIGMA(3)*IPRED**2)
		IF(CMT.EQ.7)    W = SQRT(SIGMA(4)*IPRED**2)
		IF(W.EQ.0) 		W = 1 

	IRES  = DV-IPRED
    IWRES = IRES/W 
      
    RESV = EPS(1)
    IF(CMT.EQ.5) RESV = EPS(3)
    IF(CMT.EQ.7) RESV = EPS(4)
	
	RESV2 = 0
    IF(CMT.EQ.1) RESV2 = EPS(2)
  
            
     Y = IPRED + IPRED*RESV + RESV2      


;; ------------------------- Initial estimates THETA -------------------------- ;;
$THETA
(0, 0.796,1) ; 1 ~KRES
(0, 7.28,100) ; 2 ~CL10
(0, 0.86,1) FIX ; 3 ~fd
(5.85) FIX ; 4 ~VBLO
(312) FIX ; 5 ~CO
(0, 6.4) FIX ; 6 ~PS1
(0, 3.2) FIX ; 7 ~PS2

$OMEGA 0 FIX  ; 1 ~IIV_Kp
$OMEGA 0.124 ; 2 ~IIV_CL10
$OMEGA 0 FIX  ; 3 ~IIV_fd
$OMEGA 0 FIX  ; 4 ~IIV_VBLO
$OMEGA 0 FIX  ; 5 ~IIV_CO
$OMEGA 0 FIX  ; 6 ~IIV_PS1
$OMEGA 0 FIX  ; 7 ~IIV_PS2
;; --------------------------- Initial estimates SIGMA ------------------------- ;;

$SIGMA 0.0207 ; 1 ~Plasma Prop RES_ERR
$SIGMA 1.18 ; 2 ~Plasma Add RES_ERR
$SIGMA 0.0517 ; 3 ~ECF Prop RES_ERR
$SIGMA 0.079 ; 4 ~DVE Prop RES_ERR

$ESTIMATION METHOD=1 INTER PRINT=1 MAX=9999 POSTHOC NOABORT

$COVARIANCE PRINT=E

$TABLE ID TIME TAD TIN EVID CMT FLAG DV PRED IPRED IRES WRES IWRES CWRES NOPRINT ONEHEADER FILE=sdtab22
;$TABLE ID AGE HT WT SALB SPRO CRCL NOPRINT ONEHEADER NOAPPEND FILE=cotab22
;$TABLE ID SEX DIAG STUD NOPRINT ONEHEADER NOAPPEND FILE=catab22
$TABLE ID TIME TAD TIN EVID CMT KRES CL10 fd VDVE NOPRINT ONEHEADER NOAPPEND FILE=patab22
$TABLE ID TIME TAD CMT EVID C1 Cplasma C2 C3 C4 C5 C6 C7 NOPRINT ONEHEADER FILE=mytab22