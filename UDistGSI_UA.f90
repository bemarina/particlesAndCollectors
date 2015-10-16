PROGRAM Force_Energy_Distance
IMPLICIT NONE

!PRINTS TOTAL ENERGIES - PROBABILITIES MODEL (later called GSI_UA) 
!Calculates ELECTROSTATIC forces and energies
!BOTH PARTS OF SPHERE
!--------------------------------------------------------------------------------------
DOUBLE PRECISION, parameter :: rows = 101
INTEGER, parameter :: rows1 = 101
INTEGER, DIMENSION (1:rows1, 1:rows1) :: surface
DOUBLE PRECISION, PARAMETER :: a = 500e-9	 !radius of sphere
DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793238
DOUBLE PRECISION, PARAMETER :: depth = 0 !0.5e-9
DOUBLE PRECISION :: fi_sphere, fi_surface, fi_pillars, fi_patches, D
DOUBLE PRECISION :: percent_pillars, percent_patches, percent_pill_effective
INTEGER :: row, column, i, pp, z
INTEGER*4 :: no_1, no_2, no_3
DOUBLE PRECISION :: UTotal1, UTotal2, sum_pot_elec3, sum_pot_elec4, const2, Bol, T
DOUBLE PRECISION :: Uflat, Upill, UTotal 
DOUBLE PRECISION :: FelecTotal, deb, factor

!OPEN (unit = 76, file = "KappaD.txt")
!OPEN (unit = 76, file = "DOverA.txt")
OPEN (unit = 76, file = "Dnm.txt")
OPEN(unit=78,file = "UkT.txt")

	fi_sphere = -0.054  !54 mV     
	fi_surface = -0.054			 
	fi_pillars = -0.054

	!regular surface - pillars and surface
	percent_pillars = 0      
	factor = 1
	Bol = 1.3806503e-23
	T = 298.15
	const2= Bol*T
	deb = 4e-9
	
    DO z = 1,1000
	
        IF (z==1)THEN
            D = 0.0001e-9 
        ELSE
            D=(z-1)*(0.1e-9)  
        END IF
        
        
        !WRITE(76,*) D/deb
        !WRITE(76,*) D/a
        WRITE(76,*) D/1e-9 !D in nm
        
        !Surface generation
        !Don't need specific surfaces for probabilities approach
        surface = 0 

        !To compute potential elec with a flat surface (no additional height)
        !To compute interactions with FLAT SURFACE ONLY 
        CALL potential1 (UTotal1, a, depth, rows, rows1, deb, fi_sphere, fi_surface, D, surface)
        Uflat = UTotal1
        	
        !To compute potential elec with a surface at height depth
        CALL potential2 (UTotal2, a, depth, rows, rows1, deb, fi_sphere, fi_pillars, D, surface)
        Upill = UTotal2

        UTotal = percent_pillars*factor*Upill + (1-percent_pillars*factor)*Uflat
	
        WRITE(78,*) UTotal/const2

    END DO
   
    CLOSE(76)
    CLOSE(78)
    
END PROGRAM



!----------------------------------------------------------------
SUBROUTINE  potential1 (TotalPotential1, a, depth, rows, rows1, deb, fi_sphere, fi_surface, D, surface)
IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: AH=5e-21 
DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793238
DOUBLE PRECISION :: rows
INTEGER :: rows1
DOUBLE PRECISION, DIMENSION (1:rows1, 1:rows1) :: potential_elec, potential_vdw, force_elec, force_vdw, Uelec_av
INTEGER, DIMENSION (1:rows1, 1:rows1) :: surface
DOUBLE PRECISION :: a, k, fip, fip_av, ep, h1, h2, expr1, expr2, r, deb, percent
DOUBLE PRECISION :: const1, Bol, T , const2, Fxs, expr1_f, expr2_f, const_force, mu, gamma
DOUBLE PRECISION :: const22, expr22, sum_UelecAv
DOUBLE PRECISION :: sum_pot_elec1, sum_pot_vdw, sum_force_elec1, sum_force_vdw, D, TotalPotential1, TotalForce, depth !D0
DOUBLE PRECISION :: fi_sphere, fi_surface, fi_patches, fi_pillars
DOUBLE PRECISION :: percent_pillars, percent_patches, fis1, fis2, unitArea
INTEGER :: i,j, i1, j1, z, kk  

	k=1/deb
	ep=7.08e-10
	Bol = 1.3806503e-23
	T = 298.15
	const2= Bol*T	
	mu = 0.001

	potential_elec(:,:)=0
	potential_vdw(:,:)=0
	force_elec(:,:)=0
	force_vdw(:,:)=0
	unitArea = (a*a)/((rows/2)**2)	
	
	sum_pot_elec1 = 0
	sum_pot_vdw =0
	sum_force_elec1 = 0
	sum_force_vdw =0

	DO i=1,rows1
		DO j=1,rows1
				   
			i1 = i - (rows1+1)/2  
			j1 = j - (rows1+1)/2 

			IF ((i1**2+j1**2)<=(((rows/2))**2)) THEN
			
				r =  sqrt(real(i1**2 + j1**2))*a/((rows/2))				
				h1 = a+ D - a*sqrt(1-(r/a)**2) 
				h2 = a+ D + a*sqrt(1-(r/a)**2) 
								
				!fis1 			potential of sphere
				!fis2			potential of surface - changes!!
				
				IF(surface(i,j)==0)THEN
					!Compute for flat surface
					h1 = h1 + depth
					h2 = h2 + depth
					
					fis1 = fi_sphere
					fis2 = fi_surface					
					!electrostatic energy
					const1 = ep*(fis1**2+fis2**2)*k/2
					expr1 =  const1*(1- (1/dtanh(k*h1)) + ((fis1*fis2*2)/(fis1**2+fis2**2))*(1/dsinh(k*h1)))
					expr2 =  const1*(1- (1/dtanh(k*h2)) + ((fis1*fis2*2)/(fis1**2+fis2**2))*(1/dsinh(k*h2)))
					potential_elec(i,j) =(expr1-expr2)*unitArea
                    !electrostatic force
                    const_force  = ep*(fis1**2+fis2**2)*k*k/2
                    expr1_f = -const_force*(1/dsinh(k*h1))*(1/dsinh(k*h1) - ((fis1*fis2*2)/(fis1**2+fis2**2))/dtanh(k*h1))
                    expr2_f = -const_force*(1/dsinh(k*h2))*(1/dsinh(k*h2) - ((fis1*fis2*2)/(fis1**2+fis2**2))/dtanh(k*h2))
                    force_elec(i,j)  = (expr1_f-expr2_f)*unitArea               
                    !Vdw force and energy
                    potential_vdw(i,j)= ((-AH/(12*pi*(h1**2)))-(-AH/(12*pi*(h2**2))))*unitArea
                    force_vdw(i,j)  = (((-AH)/(6*pi*(h1**3)))-((-AH)/(6*pi*(h2**3))))*unitArea		
                
				END IF							
				sum_pot_elec1 = sum_pot_elec1 + potential_elec(i,j)
				sum_force_elec1 = sum_force_elec1 + force_elec(i,j)
                
                sum_pot_vdw = sum_pot_vdw + potential_vdw(i,j)
                sum_force_vdw = sum_force_vdw + force_vdw(i,j)               
			END IF	 
			
		END DO	!DO of j's
	END DO	 !DO of i's
 
    TotalPotential1 = sum_pot_elec1 + sum_pot_vdw
		!TotalForce=	sum_force_elec + sum_force_vdw
				
	!Fxs = 6*pi*mu*a*(a)*gamma

END SUBROUTINE potential1




SUBROUTINE  potential2 (TotalPotential2, a, depth, rows, rows1, deb, fi_sphere, fi_pillars, D, surface)
IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: AH=5e-21 
DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793238
DOUBLE PRECISION :: rows
INTEGER :: rows1
DOUBLE PRECISION, DIMENSION (1:rows1, 1:rows1) :: potential_elec, potential_vdw, force_elec, force_vdw, Uelec_av
INTEGER, DIMENSION (1:rows1, 1:rows1) :: surface
DOUBLE PRECISION :: a, k, fip, fip_av, ep, h1, h2, expr1, expr2, r, deb, percent
DOUBLE PRECISION :: const1, Bol, T , const2, Fxs, expr1_f, expr2_f, const_force, mu, gamma
DOUBLE PRECISION :: const22, expr22, sum_UelecAv
DOUBLE PRECISION :: sum_pot_elec2, sum_pot_vdw, sum_force_elec2, sum_force_vdw, D, TotalPotential2, TotalForce, depth !D0
DOUBLE PRECISION :: fi_sphere, fi_surface, fi_patches, fi_pillars
DOUBLE PRECISION :: percent_pillars, percent_patches, fis1, fis2, unitArea2
INTEGER :: i,j, i1, j1, z, kk  

	k=1/deb
	ep=7.08e-10
	Bol = 1.3806503e-23
	T = 298.15
	const2= Bol*T
	gamma = 50
	mu = 0.001

	potential_elec(:,:)=0
	potential_vdw(:,:)=0
	force_elec(:,:)=0
	force_vdw(:,:)=0
	unitArea2 = (a*a)/((rows/2)**2)	
	
	sum_pot_elec2 = 0
	sum_pot_vdw =0
	sum_force_elec2 = 0
	sum_force_vdw =0
	sum_UelecAv =0

	DO i=1,rows1
		DO j=1,rows1
				   
			i1 = i - (rows1+1)/2  
			j1 = j - (rows1+1)/2 

			IF ((i1**2+j1**2)<=(((rows/2))**2)) THEN
			
				r =  sqrt(real(i1**2 + j1**2))*a/((rows/2))				
				h1 = a+ D - a*sqrt(1-(r/a)**2) 
				h2 = a+ D + a*sqrt(1-(r/a)**2) 
								
				!fis1 			potential of sphere
				!fis2			potential of surface - changes!!
				
				IF(surface(i,j)==0)THEN
					!Compute for NOT flat surface					
					fis1 = fi_sphere
					fis2 = fi_pillars						
					!electrostatic energy
					const1 = ep*(fis1**2+fis2**2)*k/2
					expr1 =  const1*(1- (1/dtanh(k*h1)) + ((fis1*fis2*2)/(fis1**2+fis2**2))*(1/dsinh(k*h1)))
					expr2 =  const1*(1- (1/dtanh(k*h2)) + ((fis1*fis2*2)/(fis1**2+fis2**2))*(1/dsinh(k*h2)))
					potential_elec(i,j) =(expr1-expr2)*unitArea2
                    !electrostatic force
                    const_force  = ep*(fis1**2+fis2**2)*k*k/2
                    expr1_f = -const_force*(1/dsinh(k*h1))*(1/dsinh(k*h1) - ((fis1*fis2*2)/(fis1**2+fis2**2))/dtanh(k*h1))
                    expr2_f = -const_force*(1/dsinh(k*h2))*(1/dsinh(k*h2) - ((fis1*fis2*2)/(fis1**2+fis2**2))/dtanh(k*h2))
                    force_elec(i,j)  = (expr1_f-expr2_f)*unitArea2  				
                    !Vdw force and energy
                    potential_vdw(i,j)= ((-AH/(12*pi*(h1**2)))-(-AH/(12*pi*(h2**2))))*unitArea2
                    force_vdw(i,j)  = (((-AH)/(6*pi*(h1**3)))-((-AH)/(6*pi*(h2**3))))*unitArea2		                 
                END IF	
								
				sum_pot_elec2 = sum_pot_elec2 + potential_elec(i,j)
				sum_force_elec2 = sum_force_elec2 + force_elec(i,j) 
                    
                sum_pot_vdw = sum_pot_vdw + potential_vdw(i,j)
                sum_force_vdw = sum_force_vdw + force_vdw(i,j)    	
			END IF	 			
		END DO	!DO of j's
	END DO	 !DO of i's
 
	TotalPotential2 = sum_pot_elec2 + sum_pot_vdw
		!TotalForce=	sum_force_elec + sum_force_vdw				
		! Fxs = 6*pi*mu*a*(a)*gamma
END SUBROUTINE potential2





 
   
     			
	           		
         

				
                
		
