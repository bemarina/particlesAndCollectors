PROGRAM Patch_Cell_Collector
IMPLICIT NONE

!*** Includes corrections of functions lambda and ji
!***Reads an already generated PATCHY SPHERE, w/o a collector discretization

!ready to go in loop for many parts

INTEGER, PARAMETER :: tsteps = 500000
INTEGER, PARAMETER :: Np = 6489 !6489 !25963  !Np is the number of ELEMENTS on the SPHERE

!INTEGER, PARAMETER :: percentInt = $PERCENT$ !For PYTHON script
!INTEGER, PARAMETER :: percentInt = 17     !This is the percent of heterogeneity ON THE SPHERE

!INTEGER, PARAMETER :: Ninitial = $NINITIAL$ !For PYTHON script
INTEGER, PARAMETER :: Ninitial = 1		     !For compilation script

!INTEGER, PARAMETER :: Nparticles = $NPARTICLES$ !For PYTHON script
INTEGER, PARAMETER :: Nparticles = 1		     !For compilation script

!DOUBLE PRECISION, PARAMETER :: percent = DBLE(percentInt)/100 
!INTEGER, PARAMETER :: Npatches = floor(percent*Np)
!INTEGER, DIMENSION (1: Npatches, 1:(Nparticles-Ninitial+1)) :: patchlocs

INTEGER, DIMENSION(1:Np) :: Stripe, contactSph
!To store 2 time steps of THETA and PHI For EACH PATCH
DOUBLE PRECISION, DIMENSION(1:Np, 1:4):: ThetaPhiPatches
DOUBLE PRECISION, DIMENSION(1:Np, 1:2):: vectorThetaPhi, lastLocs
!For the WHOLE CELL/SPHERE
DOUBLE PRECISION, DIMENSION (1:2*tsteps) :: hor, ver, saveTime, Vz, Vy

!For the COLLECTOR
INTEGER, PARAMETER :: CollecGrid = 227, cols = 9090  !this no. of cols gives a MAX distance of 20 microns. 
DOUBLE PRECISION, PARAMETER :: el_size = 2.2026e-9  !size of the square element in the collector
DOUBLE PRECISION, DIMENSION (1:CollecGrid) :: surfcent !stores CENTERS of square elements on the collector
INTEGER, DIMENSION (1:Np) :: indx, indy
!These are the COLLECTOR surfaces (one is the real one, the other stores the path, and the third stores the time step of first contact)
INTEGER, DIMENSION (1:CollecGrid, 1:cols) :: collecSurf, surfctc, firstT
INTEGER :: disp
INTEGER, PARAMETER :: NOnColl = 3245
DOUBLE PRECISION, DIMENSION (1:NOnColl, 1:3) :: forIJconv
INTEGER, DIMENSION (1:CollecGrid) :: allcols

!GENERAL parameters
DOUBLE PRECISION, PARAMETER :: h_initial = 20e-9, shear = 25 !1/sec
DOUBLE PRECISION, PARAMETER :: a = 250e-9 !250e-9 !500e-9   !a is PARTICLE RADIUS
DOUBLE PRECISION, PARAMETER :: deb = 1e-9, AH = 5.000000e-21 
DOUBLE PRECISION, PARAMETER :: fiHomoSph = -0.025, fiHetSph = 0.05
DOUBLE PRECISION, PARAMETER :: fiHomoColl = -0.025, fiHetColl = -0.025
DOUBLE PRECISION, PARAMETER :: delta = 1e-9
!BE CAREFUL with max_dist ---> so that it won't diverge when it reaches the end
!ESPECIALLY important when there's a PATCHY collector
DOUBLE PRECISION, PARAMETER :: max_dist = 20e-6     !10e-6   !30e-6 
DOUBLE PRECISION, PARAMETER :: mur = 1.3e-4 !friction coefficient	

DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846
DOUBLE PRECISION, DIMENSION (1:6, 1:6):: MobM
DOUBLE PRECISION, DIMENSION (1:6,1:1) :: ForceVec, Vel
DOUBLE PRECISION:: hdist, vdist, dt, h, Fs, Ts, mm, mu
DOUBLE PRECISION:: lambda, ji, Ft, Fr, Tt, Tr, Df
DOUBLE PRECISION:: FcollZZ, dPhiDt, dThetaDt, xx, yy, zz, Rzoi 
DOUBLE PRECISION:: factor11, factor22, sumdA, hnew
INTEGER :: t, i, j, condition, time, nInZoi, jj, total, isumtot, tt2
INTEGER :: opq, sumCtcColl, pl, pl2, t1, k1, k2, iii
INTEGER :: ColsMin, ColsMax, mb1, hj, mbmax

INTEGER :: mb, k, condition1, kk, kl, ij, km, kp, ab, kz, klp 
INTEGER :: mt, gg, jkl
DOUBLE PRECISION :: i1, j1, hdistprev
CHARACTER(4) :: num, nump
CHARACTER(5) :: numcond1, numcond2 
CHARACTER(16):: horfile, verfile, partnum, spherefile, timefile
CHARACTER(16) :: rotvelfile, transvelfile, sumdAfile
CHARACTER(20) ::  firstPatchfile, velratiofile
CHARACTER(31):: condfile, condfileS, condfileNS, condfileTO
CHARACTER(25):: DAngDtFirst, DAngDtMiddle  
	
	!For ROLLING
	DOUBLE PRECISION :: Fyn, fd, td
	
        mu = 0.001 ; mm = mu*a*pi       !mu is the viscosity
        Rzoi = 2*sqrt(deb*a)
		dt = 2e-5
		fd = 6*pi*mu*(a**2)*shear
		td = 8*pi*mu*(a**3)*shear
        
        contactSph = 0 
        
        !****** define the vector SURFCENT ***********
        !Setting the mid-point (the zero) manually and then compute the positive parts and 
        !in the same loop set the negative parts    
        !So, set the mid-point manually
        surfcent((CollecGrid+1)/2) = 0
        mt = 1;
        !the POSITIVE part of the vector
        DO gg = (CollecGrid+3)/2 , CollecGrid
            surfcent(gg) = surfcent(gg-1) + el_size
            surfcent(gg-2*mt) = -surfcent(gg)
            mt = mt + 1
        END DO

        !**************** end SURFCENT definition
        
        
        !***** Create the collector Surface: Either read patches locations (as before) or define
        !**** a uniform collector
        !-----------------------------------------------------------------------------------
        collecSurf = 0 
        
        !********************
        !Initialize the matrix that will store the path
        surfctc = 0
        
        !Initialize the matrix that will store the t steps of first contact
        firstT = 0 
        !*****************************
        
        OPEN(unit = 225, file = "TestingPatchesNo.txt")
        OPEN(unit = 97, file = "LocsOnCollCtc.txt")
        
        
		 !To convert INTEGER to STRING
        write(numcond1,'(i5.5)') Ninitial
        write(numcond2,'(i5.5)') Nparticles
		
		condfile = 'condition'//numcond1//'-'//numcond2//'.txt'
        condfileS = 'condition'//numcond1//'-'//numcond2//'Stuck.txt'
        condfileNS = 'condition'//numcond1//'-'//numcond2//'NoStuck.txt'
        condfileTO = 'condition'//numcond1//'-'//numcond2//'TimeOut.txt'
		
		OPEN (unit = 6, file = condfile)
		OPEN (unit = 7, file = condfileS)
		OPEN (unit = 8, file = condfileNS)
        OPEN (unit = 9, file = condfileTO)
   
		!ThetaPhiPatches stores Theta and Phi for each patch        
        !the are 4 columns: Theta(present) Phi(present) Theta(new) Phi(new)
        !with 0<Theta<pi and 0<Phi<2*pi    
		
		!First need to load the elements as produced by the EQSP algorithm 
		!The 0 as the first Theta gives NaN so I changed it into a very small number
		!(2.22E-16 - the "eps" value in matlab) in the file "PlainElementsUsedInSims.txt"
        !So, instead of using "PlainElements.txt", I use "PlainElementsUsedInSims.txt"   
		

		!For MANY patchy particles, store the locations of patches in a matrix
		!of size patchlocs = Npatches x Nparticles
		!For each particle read the respective column and form the vector "Stripe"
        !So, no need to open this "RandomPatches00x.txt" files
        !OPEN(unit=225, file = "RandomPatches003.txt")
        !READ(225,*) (Stripe(i), i=1,Np)
        !CLOSE(225)
        
		!****************
		!Stripe = 0 !for the case of UNIFORM SPHERE
    !*****************

		!**** create the matrix with the patches locations
		!patchlocs = 0  

		!DO mb = 1, (Nparticles-Ninitial+1)
		!	k = 1		
			!Assign the first one
		!	CALL random_seed
		!	CALL random_number (i1)			
		!	ij = 1 + (Np-1)*(i1)
		!	patchlocs(k,mb) = ij
		!	k = 2
			!Assigning the rest - beeware of repetitions
		!	CALL random_seed
			
		!	DO 
		!		condition1 = 0 
		!		IF (k>Npatches) EXIT
				
		!		CALL random_number (i1)			
		!		ij= 1 + (Np-1)*(i1)
				
				!checking the number is not on the list
		!		DO kk = 1, k-1				
		!			IF (patchlocs(kk,mb)==ij) THEN 
		!				condition1 = 1 
		!				EXIT
		!			END IF
		!		END DO	
				!if it's not a repetition, then put it in the list	
		!		IF (condition1 ==0) THEN
		!			patchlocs(k,mb) = ij
		!			k = k+1
		!		END IF
				
		!	END DO			
			!If want to save the random particle to reproduce later
		!	write(nump,'(i4.4)') mb
		!	partnum = 'sphere'//nump//'.txt'				
			!To save when doing a loop, will save later as a vector of 0's and 1's
		!		OPEN (unit = 85, file = partnum)
		!		DO ab = 1, Npatches
		!			WRITE(85,*) patchlocs(ab,mb)
		!		END DO
		!		CLOSE(85)			
		!END DO		
		!*********************** end make matrix with patch locations
		isumtot = 0 

		!*** To use an existing patchy particle (no need to generate it)
		!The inputed file should NOT be the vector with 0's and 1's, but instead
		!the list of the LOCATIONS of the patches
		!OPEN (unit = 10, file = "sphere0005.txt")
		!	READ(10,*) ((patchlocs(i,j), j=1,1), i=1,Npatches)	
		!CLOSE(10)
		!OPEN (unit = 11, file = "sphere0002.txt")
		!	READ(11,*) ((patchlocs(i,j), j=2,2), i=1,Npatches)	
		!CLOSE(11)
        !If the inputed file is the vector with 0's and 1's then use
            Stripe = 0
			OPEN (unit = 10, file = "sphere0001A250DoNotErase.txt")		
            !OPEN (unit = 10, file = "sphere0001A500DoNotErase.txt")		
			READ(10,*) (Stripe(i), i=1,Np)		
			CLOSE(10)
		!to allocate the 0's and 1's into the stripe variable directly
		!but then, also need to change a few other things down below... careful!
		!************************************
		kz = 1
		!***Start loop over many particles
		DO km = Ninitial, Nparticles
		
			!*****For every particle, I have to build it "mehadash"
			ThetaPhiPatches = 0; vectorThetaPhi = 0 
			condition = 0        
			OPEN(unit=23, file = "PlainElementsUsedInSimsA250Nm.txt") 
            !OPEN(unit=23, file = "PlainElementsUsedInSimsA500Nm.txt")   	            
				READ(23,*) ((ThetaPhiPatches(i,j), j=1,2), i=1,Np)		
			CLOSE(23)
			!*****

      !Stripe = 0 - need to annulate this here b/c it's READING the sphere
			!For each particle create the variable "stripe"
			!DO kp = 1, Npatches
			!	Stripe(patchlocs(kp,kz)) = 1
			!END DO
            
			WRITE(225,*) sum(stripe)
			kz = kz+1
			!For each particle create the respective hor/ver/sphere files			
			write(num,'(i4.4)') km			
			horfile = 'hor'//num//'.txt'
			verfile = 'ver'//num//'.txt'
			spherefile = 'sphere'//num//'.txt'
            timefile = 'time'//num//'.txt'
			rotvelfile = 'rotVel'//num//'.txt'
            transvelfile = 'transVel'//num//'.txt'
            firstPatchfile = 'ThetaOnSph'//num//'.txt'
            velratiofile = 'VelRatios'//num//'.txt'
            !sumdAfile = 'SumDA'//num//'.txt'
            DAngDtFirst = 'DAnglesDtFirst'//num//'.txt'
            DAngDtMiddle = 'DAnglesDtMid'//num//'.txt'
            
			IF((km==1).OR.(km==5).OR.(km==10).OR.(mod(km,100)==0)) THEN				
				OPEN (unit = 45, file = horfile)
				OPEN (unit = 46, file = verfile)
				!OPEN (unit = 47, file = spherefile)
                !OPEN (unit = 48, file = timefile)
                !OPEN (unit = 49, file = rotvelfile)
                !OPEN (unit = 50, file = transvelfile)
                !OPEN (unit = 51, file = firstPatchfile)
                !OPEN (unit = 52, file = velratiofile)
				!DO klp = 1, Np
				!	WRITE(47,*) Stripe(klp)
				!END DO
				!CLOSE(47)			
                !OPEN (unit = 53, file = DAngDtFirst)                
                !OPEN (unit = 54, file = DAngDtMiddle)
            END IF	
        
			!CHECK differences between DSIN and SIN
			hdist = 0.000;          vdist = h_initial
			Vz = 0; Vy = 0        ; Vel = 0 ; lastLocs = 0 
			t=1  !for while 
			!DO t = 1, tsteps  !for a specific number of time steps
			DO  WHILE(condition==0)
		
                hor(t)=hdist;   ver(t)=vdist               				
                saveTime(t) = dt*(t-1)              
                
				IF((km==1).OR.(km==5).OR.(km==10).OR.(mod(km,100)==0)) THEN		
					WRITE(45,*) hdist/a 
					WRITE(46,*) vdist/a
                    !WRITE(48,*) saveTime(t)
                END IF
              
                
                h = vdist
                CALL func_ft(h,a,Ft); CALL func_tt(h,a,Tt)
                CALL func_fr(h,a,Fr); CALL func_tr(h,a,Tr)
                Df = Tt*Fr - Ft*Tr
                CALL func_lambda(h,a,lambda) !lambda=2.1957
                CALL func_ji(h,a,ji) ! ji=1.1384
                !Valid only for small sep distances
                !Fs = 6*mm*(a+h)*shear*(1+ 9*a/(16*(a+h)))
                !Ts = 4*mm*a*a*(1-(3/16)*((a/(a+h))**3))*shear 
                !Valid for larger sep distances - From Ranojoy's paper JCIS 2007
                factor11 = (1.7007337+ 1.0221606*(h/a))/(1+1.0458291*(h/a)-0.0014884706*((h/a)**2)) 
                Fs = 6*mm*(a+h)*shear*factor11
                factor22 = 0.054651334*(18.276952- exp(-1.422943*(h/a)))
                Ts = 4*mm*a*a*shear*factor22                

                !Define the Mobility Matrix
                MobM = 0
                MobM(1,1)= Tr/(6*pi*mu*a*Df); MobM(2,2)=	1/(6*pi*mu*a*lambda)
                MobM(3,3)= Tr/(6*pi*mu*a*Df); MobM(4,3)= Fr/(8*pi*mu*(a**2)*Df)
                MobM(6,1)=(-Fr)/(8*pi*mu*(a**2)*Df); MobM(1,6)= (-Tt)/(6*pi*mu*(a**2)*Df)
                MobM(3,4)= Tt/(6*pi*mu*(a**2)*Df); MobM(4,4)= Ft/(8*pi*mu*(a**3)*Df)
                MobM(5,5)= 1/(8*pi*mu*(a**3)*ji); MobM(6,6)= Ft/(8*pi*mu*(a**3)*Df)
             
                		
                vectorThetaPhi = ThetaPhiPatches(:,1:2)
                
                !****** convert points on sphere to locations on collector to 
                !****** find out what is the potential (or AH) at that point on the collector
                indx = 0 ; indy = 0
                CALL ConvertXYtoIJ (a, Np, CollecGrid, vectorThetaPhi, surfcent, el_size, indx, indy)
                !indx is in fact the column
                !indy is in fact the row
                
                !Adding the horizontal displacement
                disp = floor(hdist*(CollecGrid/2)*(1/a))
                indx = indx + disp
                
                !IF(t==66)THEN
                !    OPEN(unit = 54, file = "DotsOnCollt66First.txt")       
                !END IF
                
                
                mb1 = 1
                allcols = 0 
                !IF(t==2000)THEN
                !    OPEN(unit = 57, file = "SphOnCollt2000.txt")
                IF(t==7000)THEN
                    OPEN(unit = 57, file = "SphOnCollt7000.txt")
                    DO t1 = 1, Np
                        WRITE(57,*) indy(t1), indx(t1)
              
                        IF (indy(t1)==114)THEN
                            allcols(mb1)=indx(t1)
                            mb1 = mb1 + 1
                        END IF
                        
                    END DO
                    
                    DO hj = 1,CollecGrid
                        IF (allcols(hj)==0) THEN
                            mbmax = hj-1
                            EXIT
                        END IF
                    END DO
                    
                    
                    colsMin = MINVAL(allcols(1:mbmax))
                    colsMax = MAXVAL(allcols(1:mbmax))
                    
                    print*, colsMin, colsMax
                    
                    CLOSE(57)
                END IF
            
                
                
                !******* Compute DLVO(colloidal) forces ***********************************
                !Call the subroutine to compute Fcoll = Fvdw + Felec for ALL 
                !the patches - subroutine yields total force for all patches at a 
                !given ditance VDIST - The total force is FcollZZ, because it's in the z-direction
                sumdA = 0 
                forIJconv = 0
                CALL  FcollZ (a, vdist, deb,AH,fiHomoSph, fiHetSph, vectorThetaPhi, Np, Stripe,FcollZZ, t, indx, indy, CollecGrid, cols, collecSurf, fiHetColl, fiHomoColl, contactSph, delta, NOnColl, forIJconv)
                	
                !IF(t==66)THEN
                !    OPEN(unit = 233, file = "DotsOnCollt66Only.txt")
                !    DO iii=1,NOnColl
                !        IF(forIJconv(iii,1)==0)THEN
                !            EXIT
                !        ELSE
                    
                !            WRITE(233,*) forIJconv(iii,1), forIJconv(iii,2), forIJconv(iii,3)
                !        END IF
                    
                !    END DO
                !    CLOSE(233)
                !END IF    
                    
                    
                    
                !To save in a new matrix the amount of times each collector element is contacted by the sphere    
                !hopefully be able to trace the sphere's TRAJECTORY on the collector
                CALL convXYtoIJContct(t, NOnColl, forIJconv, CollecGrid, cols, surfcent, el_size, a, disp, surfctc, firstT)    
                   
                    
                !IF(t==65)THEN
                !    OPEN(unit = 98, file = "DotsOnCollt65.txt")
                !    DO k1 = 1, collecGrid
                !        DO k2 = 1, cols                      
                !            IF (surfctc(k1,k2)>0)THEN
                !                WRITE(98,*) k1, k2, surfctc(k1,k2)
                !            END IF                
                !        END DO
                !    END DO
                !    CLOSE(98)
                !END IF
                    
                IF(t==2000)THEN
                    OPEN(unit = 93, file = "DotsOnCollt2000OnProjecOnly.txt")
                    !Here's the problem: it's reading from the WHOLE surface!
                    !If want for a specific tstep have to either choose the projection only
                    !or find a way to identify which ones were added at a specific t-step
                    DO k1 = 1, collecGrid
                        DO k2 = colsMin, colsMax                      
                            IF (surfctc(k1,k2)>0)THEN
                                WRITE(93,*) k1, k2, surfctc(k1,k2)
                            END IF                
                        END DO
                    END DO
                    CLOSE(93)
                END IF    
                
               
                ForceVec = 0 
                ForceVec(1,1)= Fs 
                ForceVec(2,1)= FCollZZ  
                ForceVec(6,1)= Ts
                
                !IF((km==1).OR.(km==5).OR.(km==10).OR.(mod(km,100)==0)) THEN		
                    !Rotational Velocity: Omega_x is Vel(6,1)
                !    WRITE(49,*) Vel(6,1)*a
                    !Translational Velocity: V_y is Vel (1,1)
                !    WRITE(50,*) Vel(1,1)
                    !Sum of areas of elements ON THE SPHERE
                    !should be = 2*(projected area of sphere) = 2*(pi*(a**2))
                !    WRITE(52,*) (Vel(6,1)*a)/Vel(1,1)
                !END IF
                

                !Compute velocities - matrix multiplication
                DO i=1,6 
                        Vel(i,1)= MobM(i,1)*ForceVec(1,1)+MobM(i,2)*ForceVec(2,1)+MobM(i,3)*ForceVec(3,1)&
		&+MobM(i,4)*ForceVec(4,1)+MobM(i,5)*ForceVec(5,1)+MobM(i,6)*ForceVec(6,1)		
                END DO
                !hdist, vdist (centre of sphere wrt surface) as before
                time = t+1				
                
	
				hdistprev = hdist/a	! to keep the hordist before it's updated
                !----------------------------------------------------------------
                !UPDATE the sphere's position - horizontal and vertical positions
                !----------------------------------------------------------------
                hdist=hdist + (dt)*Vel(1,1)
                vdist=vdist + (dt)*Vel(2,1)

                !Updates for each patch: NEED loop here: for i = 1:Np 
                DO i = 1,Np
                        !for these expressions need the PRESENT values of Theta and Phi
                        dThetaDt = (cos(ThetaPhiPatches(i,2)))*Vel(4,1) - (sin(ThetaPhiPatches(i,2)))*Vel(6,1)	               			
                        dPhiDt = ((-cos(ThetaPhiPatches(i,2)))/(tan(ThetaPhiPatches(i,1))))*Vel(6,1)+&
		&((-sin(ThetaPhiPatches(i,2)))/(tan(ThetaPhiPatches(i,1))))*Vel(4,1)+Vel(5,1)
                        !calculate new Theta
                        ThetaPhiPatches(i,3) = dThetaDt*dt + ThetaPhiPatches(i,1)
                        !calculate new Phi
                        ThetaPhiPatches(i,4) = dPhiDt*dt + ThetaPhiPatches(i,2)                           
                END DO  !end loop of updates of patches positions

                !last locations "USED"
                lastLocs(:,1) = ThetaPhiPatches(:,1)
                lastLocs(:,2) = ThetaPhiPatches(:,2)
                
                !Now, need to move all new values from columns 3 and 4 to 1 and 2
				ThetaPhiPatches(:,1) = ThetaPhiPatches(:,3) ; ThetaPhiPatches(:,3) = 0
				ThetaPhiPatches(:,2) = ThetaPhiPatches(:,4) ; ThetaPhiPatches(:,4) = 0
                						
                            
                IF(vdist<=delta)THEN
                        Fyn = -ForceVec(2,1)
						IF((ForceVec(1,1)+(ForceVec(6,1)/a)-(mur*Fyn))>0)THEN	
						!The particle is rolling, so correct the already updated values hdist and vdist 
						!with those for the rolling conditions										
							vdist = delta	!this can be here or in the same loop as the definition of Fyn					
							!***** THIS *****
                            !Vy(t+1) = shear*a*((mur*(Fyn/fd)-(4./3)*(ForceVec(6,1)/td)-(ForceVec(1,1)/fd))/(Ft+Fr+(4./3)*(Tt+Tr)))
							!rolling velocity
							!Vz(t+1) = 0
							!no velocity in perpendicular direction	
                            !hdist = hdistprev*a + (dt)*Vy(t+1) 
                            !******---------------------------------------------
                            
                            !**OR***
                            !--------------- THIS ------------------------------
                            !------- Hydrodynamic functions for the ROLLING velocity?!
                            hnew = vdist
                            CALL func_ft(hnew,a,Ft); CALL func_tt(hnew,a,Tt)
                            CALL func_fr(hnew,a,Fr); CALL func_tr(hnew,a,Tr)
                            !ROLLING VELOCITY
                            Vel(1,1) = shear*a*((mur*(Fyn/fd)-(4./3)*(ForceVec(6,1)/td)-(ForceVec(1,1)/fd))/(Ft+Fr+(4./3)*(Tt+Tr)))
							!NO velocity in perpendicular direction		
							Vel(2,1) = 0
							!UPDATE HORIZONTAL Position
							hdist = hdistprev*a + (dt)*Vel(1,1)
                            !ALSO correct the rotational velocity which has to be equal to the translational velocity
                            Vel(6,1) = (Vel(1,1))/a
                           
						ELSE					 
							condition=1 !particle is arrested
							WRITE(7,*) condition, hdistprev
							isumtot = isumtot + 1                            
                        END IF
				END IF
				
				
				IF(hdist>= max_dist)THEN
                        condition=2
						WRITE(8,*) condition
                END IF		              
				
                !For WHILE 		
                t = t+1 
                IF(t>tsteps)THEN
                !IF(t==7000)THEN        !When want to stop the sphere at a certain time step                
                        condition=3 !max. number of steps is reached
						WRITE(9,*) condition
                END IF  
		
        END DO         !*** end loop of time steps / WHILE "condition==0" loop 
		
		IF (condition==1)THEN
			WRITE(6,*) condition, hdistprev    
		ELSE
			WRITE(6,*) condition
		END IF
		IF((km==1).OR.(km==5).OR.(km==10).OR.(mod(km,100)==0)) THEN		
			CLOSE(45)
			CLOSE(46)
            !CLOSE(48)
            !CLOSE(49)
            !CLOSE(50)
            !CLOSE(51)
            !CLOSE(52)
            !CLOSE(53)
            !CLOSE(54)
		END IF
		
		END DO !**** end loop for many particles
       
        CLOSE(225)

        
        !DO opq = 1, Np
        !    print*, contactSph(opq)
        !END DO
        !print*, 'sumOnSphere', sum(contactSph)
        
        sumCtcColl = 0 
        DO pl = 1, collecGrid
            DO pl2 = 1, cols
                IF (surfctc(pl,pl2)>0)THEN
                    WRITE(97,*) pl, pl2, surfctc(pl,pl2)
                    sumCtcColl = sumCtcColl + 1
                END IF                
            END DO
        END DO
        print*, 'sumOnCollec', sum(surfctc)
        
        
        !DO pl = 1, collecGrid
        !    DO pl2 = 1, cols
        !        IF (firstT(pl,pl2) ==66)THEN
        !            WRITE(54,*) pl, pl2, firstT(pl,pl2)
        !        END IF                
        !    END DO
        !END DO
        
        
        
		WRITE(6,*) isumtot, 'total number adhered'
		WRITE(7,*) isumtot, 'total number adhered'		
		CLOSE(6) 
		CLOSE(7)
		CLOSE(8)
        CLOSE(9)
        CLOSE(97)
        !CLOSE(54)
     
END PROGRAM
 

 
!SUBROUTINES -------------------------------------------------------------
!******* convert points on sphere to locations on the collector
!-------------------------------------------------------------------------
SUBROUTINE ConvertXYtoIJ (a, Np, CollecGrid, sphere, surfcent, el_size, indx, indy)
IMPLICIT NONE

!INPUTS
DOUBLE PRECISION, DIMENSION(1:Np,1:2):: sphere   !sphere contains the "present" theta and phi angles
INTEGER*4 :: Np
INTEGER :: CollecGrid
DOUBLE PRECISION :: a, el_size
DOUBLE PRECISION, DIMENSION (1:CollecGrid) :: surfcent 

!VARIABLES in the subroutine
DOUBLE PRECISION, DIMENSION(1:Np,1:3):: vectorXYZ
INTEGER :: i, j

!OUTPUTS
INTEGER, DIMENSION (1:Np) :: indx, indy

        
        vectorXYZ(:,1) = a*dsin(sphere(:,1))*dcos(sphere(:,2))
        vectorXYZ(:,2) = a*dsin(sphere(:,1))*dsin(sphere(:,2))
        vectorXYZ(:,3) = a*dcos(sphere(:,1))
        
        !First convert the X values to coordinates on the collector 
        DO i = 1, Np !run for ALL the points ON THE SPHERE 
            DO j = 1, CollecGrid
                IF (((surfcent(j)-(el_size/2))< vectorXYZ(i,1)).AND. &
                    &(vectorXYZ(i,1)<(surfcent(j)+ (el_size/2))))THEN
                    indx(i) = j
                    EXIT
                else
                    continue    !maybe here should use CYCLE instead of continue? check wikipedia page
                END IF
            END DO
        END DO
        
        !Second convert the Y values to coordinates on the collector 
        DO i = 1, Np !run for ALL the points ON THE SPHERE 
            DO j = 1, CollecGrid
                IF (((surfcent(j)-(el_size/2))< vectorXYZ(i,2)).AND. &
                    &(vectorXYZ(i,2)<(surfcent(j)+ (el_size/2))))THEN
                    indy(i) = j
                    EXIT
                else
                    continue    !maybe here should use CYCLE instead of continue? check wikipedia page
                END IF
            END DO
        END DO
   
END SUBROUTINE ConvertXYtoIJ
!-------------------------------------------------------------------------------------------
!*****Convert to IJ locations on the collector the elements of the particle that become into "contact"
!***** with the collector
!-------------------------------------------------------------------------------------------
SUBROUTINE convXYtoIJContct(t, NOnColl, forIJconv, CollecGrid, cols, surfcent, el_size, a, disp, surfctc, firstT)
IMPLICIT NONE

    !** Inputs
    INTEGER :: NOnColl, CollecGrid, cols, disp, t
    DOUBLE PRECISION, DIMENSION(1:NOnColl,1:3) :: forIJconv
    DOUBLE PRECISION :: a, el_size
    DOUBLE PRECISION, DIMENSION (1:CollecGrid) :: surfcent 
    
    !** Outputs
    !** This is the COLLECTOR surface that stores the number of times that each element is contacted
    !** by the sphere
    INTEGER, DIMENSION (1:CollecGrid, 1:cols) :: surfctc  
    !** This is the COLLECTOR surface that stores the time step in which the collector element is 
    !** contacted by the sphere the FIRST TIME
    INTEGER, DIMENSION (1:CollecGrid, 1:cols) :: firstT  
    
    !** In subroutine
    INTEGER :: i23, kl, j, jj, indxx, indyy
    DOUBLE PRECISION :: theta, phi, xx, yy, zz 
    
    !IF(t==66)THEN
    !    OPEN(unit = 54, file = "DotsOnCollt66Only.txt")       
    !END IF
    
    !shorten the vector 'forIJconv' - only a few rows were used, not the total NOnColl   
        i23 = 0 
        kl = 1
        DO WHILE (i23 == 0)
            IF (forIJconv(kl,1)==0) THEN
                i23 = 1
            ELSE
                !!store theta
                theta = forIJconv(kl,2)
                !!store phi
                phi = forIJconv(kl,3)
                kl = kl + 1
                
                xx = a*dsin(theta)*dcos(phi)
                yy = a*dsin(theta)*dsin(phi)
                zz = a*dcos(theta)
                
                !Convert the X coordinate to indx 'j' (column)
                DO j = 1, CollecGrid
                    IF (((surfcent(j)-(el_size/2))< xx).AND.(xx<(surfcent(j)+ (el_size/2))))THEN
                        indxx = j
                        EXIT
                    else
                        continue    !maybe here should use CYCLE instead of continue? check wikipedia page
                    END IF
                END DO

                !Second convert the Y values to coordinates on the collector
                !Convert the Y coordinate to indy 'jj' (row)             
                DO jj = 1, CollecGrid
                    IF (((surfcent(jj)-(el_size/2))<yy).AND.(yy<(surfcent(jj)+ (el_size/2))))THEN
                        indyy = jj
                        EXIT
                    else
                        continue    !maybe here should use CYCLE instead of continue? check wikipedia page
                    END IF
                END DO
                
                !Now, hopefully have the pair(i,j) = (indyy, indxx) so want to put a '1' on the surface at that point
                indxx = j + disp
                
                !Use ANOTHER matrix that stores which is the TIME STEP 't' in which the collector element is FIRST contacted
                IF(surfctc(indyy, indxx)==0)THEN
                    firstT(indyy,indxx)= t
                END IF
                
                surfctc(indyy, indxx) = surfctc(indyy, indxx) + 1 
                
                !IF(t==66)THEN
                !    WRITE(54,*) indyy, indxx, surfctc(indyy, indxx)
                !END IF
        
            END IF
        END DO
                 
                
        !IF(t==66)THEN
        !    CLOSE(54)
        !END IF

 END SUBROUTINE  convXYtoIJContct              
                

                
                
!---------------------------------------------------------------------------------                
!************************* Colloidal Forces **************************************
!---------------------------------------------------------------------------------
SUBROUTINE  FcollZ (a, vdist, deb,AH,fiHomoSph, fiHetSph, sphere, Np, Stripe,FcollZZ, t, indx, indy, CollecGrid, cols, collecSurf, fiHetColl, fiHomoColl, contactSph, delta, NOnColl, forIJconv)
IMPLICIT NONE
	
!INPUTS
INTEGER*4 :: Np
INTEGER :: cols, CollecGrid
INTEGER, DIMENSION(1:Np) :: Stripe, contactSph
DOUBLE PRECISION :: a, deb, AH, vdist, fiHetSph, fiHomoSph, fiHetColl, fiHomoColl
DOUBLE PRECISION :: delta
DOUBLE PRECISION, DIMENSION(1:Np,1:2):: sphere   !sphere contains the "present" theta and phi angles
INTEGER, DIMENSION (1:CollecGrid, 1:cols) :: collecSurf 
INTEGER, DIMENSION (1:Np) :: indx, indy
INTEGER :: t, NOnColl
DOUBLE PRECISION, DIMENSION (1:NOnColl, 1:3) :: forIJconv

!OUTPUT
DOUBLE PRECISION:: FcollZZ, sumdA


!VARIABLES IN THE SUBROUTINE
DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846
DOUBLE PRECISION :: h1, h2, fiCollector, fiSphere, constForce, expr1, expr2, k, D, dpatch
DOUBLE PRECISION :: EnrgElec, EnrgVdw, unitArea, ep, dSi, norm, term18
DOUBLE PRECISION :: Bol, Temp, const2, Felec, FVDW, r, dA, ForceElec, ForceVdw
DOUBLE PRECISION, DIMENSION(1:Np,1:3):: vectorXYZ
DOUBLE PRECISION, DIMENSION(3):: normal
INTEGER :: z, i, i1, j1, numcontact
INTEGER :: nn
	
        IF(t==65)THEN
            OPEN(unit=67,file='FXYconv.txt')
        END IF
    
        ep=7.08e-10; k = 1/deb
        Bol = 1.3806503e-23; Temp = 298.15; const2= Bol*Temp
        numcontact = 0
        
        !These variables will accumulate the sum of forces
        FcollZZ = 0;  FElec = 0;  FVDW = 0
        
        !I think this is faster than doing in a loop - like in MATLAB
        vectorXYZ(:,1) = a*dsin(sphere(:,1))*dcos(sphere(:,2))
        vectorXYZ(:,2) = a*dsin(sphere(:,1))*dsin(sphere(:,2))
        vectorXYZ(:,3) = a*dcos(sphere(:,1))
        
        dpatch = 11e-9 ;        dSi = dpatch**2 !AREA IS A SQUARE OF dpatch*dpatch
        sumdA = 0 
        
        !want to know the locations ON the collector of the sphere elements that touched
        !the collector (went below 5 nm) -- this is the purpose
        forIJconv = 0 !define it as a vector of 3245 - to be on the safe side; 
        nn = 1 
        !may need to modify varying on the minimum h1
        !b/c of the arguments to use the subroutine Np will be a different (smaller) number
        
        
        DO i = 1, Np                                                   
                !compute projected area, for that need first the normal vector
                !compute vector n --> for this need x, y, z of the patch!! 
                !normal(1) = vectorXYZ(i,1)/sqrt(a**2-vectorXYZ(i,1)**2-vectorXYZ(i,2)**2);
                !normal(2) = vectorXYZ(i,2)/sqrt(a**2-vectorXYZ(i,1)**2-vectorXYZ(i,2)**2);
                !normal(3) = 1; 
                !unit normal vector
                !normal = (sqrt(a**2-vectorXYZ(i,1)**2-vectorXYZ(i,2)**2)/a)*normal
                
                normal(1) = vectorXYZ(i,1)/a 
                normal(2) = vectorXYZ(i,2)/a            
                term18 = a**2-vectorXYZ(i,1)**2-vectorXYZ(i,2)**2
                
                IF(term18<0)THEN
                        !PRINT*, t, i, sphere(i,1), sphere(i,2)
                        term18 = 0
                END IF
                
                normal(3) = (sqrt(term18))/a               
                !checking it's the unit vector --> it's norm should be 1
                norm = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)                
                                
                dA = abs(normal(3))*dSi  
                
                sumdA = sumdA + dA
                
                !h1 is the vertical distance b/t the point on the sphere's surface
                !and the collector (top of pillars or collector itself)
                h1 = a+ vdist + vectorXYZ(i,3)  
	                        
                IF(h1<=5e-9)THEN
                    numcontact = contactSph(i) + 1
                    contactSph(i) = numcontact 
                    !print*, t, i, sphere(i,1)
                    !store in for IJconv: i theta phi
                    forIJconv(nn,1) = i
                    forIJconv(nn,2) = sphere(i,1)
                    forIJconv(nn,3) = sphere(i,2)
                    print*, h1, vdist     !to get the distances h1 = local distance from the sphere surface to the collector 
                    IF(t==65)THEN
                        WRITE(67,*) forIJconv(nn,1), forIJconv(nn,2), forIJconv(nn,3) 
                    END IF
                    nn = nn + 1
                    
                END IF
                
                !electrostatic force
                IF(stripe(i)==1)THEN
                        fiSphere = fiHetSph  !of SPHERE heterogeneity					        
                ELSEIF(stripe(i)==0)THEN
                        fiSphere= fiHomoSph   !of SPHERE parts w/o heterogeneity
                END IF  
                
                !*** When the collector is NOT discretized
                !fiCollector = fiHomoSph
                
                !*** When the collector IS discretized
                IF(collecSurf(indy(i),indx(i))==1) THEN
                        fiCollector = fiHetColl  !of COLLECTOR heterogeneity					        
                ELSEIF(collecSurf(indy(i),indx(i))==0) THEN
                        fiCollector = fiHomoColl   !of COLLECTOR parts w/o heterogeneity
                END IF  
                
                
                constForce = ep*(fiSphere**2+fiCollector**2)*k*k/2
                expr1 = -constForce*(1/dsinh(k*h1))*(1/dsinh(k*h1) - ((fiSphere*fiCollector*2)/(fiSphere**2+fiCollector**2))/dtanh(k*h1))
                ForceElec = expr1*dA	
                
                !Vdw force                               
                ForceVdw = (-AH/(6*pi*(h1**3)))*dA
                                
                IF(vectorXYZ(i,3)>0)THEN     
                        FVDW = FVDW - ForceVdw
                        FELEC = FELEC - ForceElec
                        FcollZZ = FcollZZ - (ForceElec + ForceVdw)     
                ELSEIF(vectorXYZ(i,3)<=0)THEN                     
                        FVDW = FVDW + ForceVdw
                        FELEC = FELEC + ForceElec
                        FcollZZ = FcollZZ + (ForceElec + ForceVdw)
                END IF

                        
        END DO
        
        IF(t==65)THEN
            CLOSE(67)
        END IF
    
                
END SUBROUTINE FcollZ	





!************* Hydrodynamic functions ****************************
!----------------------------------------------------------------	
	SUBROUTINE func_ft (h,a,Ft)
		IMPLICIT NONE
		DOUBLE PRECISION:: h, a, Ft
			Ft= (-1)/(0.14116 + (0.5967*((h/a)**0.2984)))
		END SUBROUTINE func_ft


	SUBROUTINE func_tt (h,a,Tt)
		IMPLICIT NONE
		DOUBLE PRECISION:: h, a, Tt
			Tt= (0.04362 - (0.0459*((h/a)**0.557)))/(0.06801 + ((h/a)**0.557))
		END SUBROUTINE func_tt


	SUBROUTINE func_fr (h,a,Fr)
		IMPLICIT NONE
		DOUBLE PRECISION:: h, a, Fr
			Fr = (0.05826 - (0.06126*((h/a)**0.557)))/(0.0681 + ((h/a)**0.557)) 
		END SUBROUTINE func_fr


	SUBROUTINE func_tr (h,a,Tr)
		IMPLICIT NONE
		DOUBLE PRECISION:: h, a, Tr
			Tr = (-0.312373 - (0.739*((h/a)**0.4906)))/(0.0954 + ((h/a)**0.4906)) 
		END SUBROUTINE func_tr

	SUBROUTINE func_lambda (h,a,lambda)
		IMPLICIT NONE
		DOUBLE PRECISION:: h, a, lambda, alpha, beta, gamma, sum1, den, n1
        INTEGER :: n
			alpha = log( ((a+h)/a) + sqrt(((1+(h/a)) **2)-1))
			sum1 = 0.0
		DO n = 1,200
            n1 = DBLE(n)
			den = 4*((dsinh((n1+0.5)*alpha))**2)-(((2*n1+1)**2)*(dsinh(alpha))**2)
			beta = (2*dsinh((2*n1+1)*alpha) + (2*n1+1)*dsinh(2*alpha))/den 
			gamma =   ((n1*(n1+1))/((2*n1-1)*(2*n1+3))) *(beta-1)
			sum1 = sum1 + gamma
		END DO
		lambda = (4./3)*(dsinh(alpha))* sum1
		END SUBROUTINE func_lambda


	SUBROUTINE func_ji (h,a,ji)
		IMPLICIT NONE
		DOUBLE PRECISION:: h, a, ji_1, sum1, alpha, ji,n1
        INTEGER :: n
			alpha = log( ((a+h)/a) + sqrt(((1+(h/a))**2)-1))
			sum1=0.0
			DO n=1,200
                n1 = DBLE(n)
				ji_1 = ((1/(dsinh(n1*alpha)))**3)/((1/(dsinh(alpha)))**3)
				sum1 = sum1 + ji_1
			END DO
			ji=sum1
		END SUBROUTINE func_ji
					




