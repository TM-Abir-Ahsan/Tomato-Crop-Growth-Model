 Subroutine Type216

! Object: Type216
! Simulation Studio Model: Type216
! 

! Author: T M ABIR AHSAN 
! Editor: T M ABIR AHSAN
! Date:	 February 29, 2024
! last modified: February 29, 2024
! 
! 
! *** 
! *** Model Parameters 
! *** 
!			AREA_GH	m^2 [-Inf;+Inf]
!			VEG_FRAC	- [-Inf;+Inf]
!			LAI_Max	- [-Inf;+Inf]
!			CROP_density	- [-Inf;+Inf]

! *** 
! *** Model Inputs 
! *** 
!			T_GH	C [-Inf;+Inf]
!			C_CO2	any [-Inf;+Inf]
!			Q_SOL_IN	W/m^2 [-Inf;+Inf]

! *** 
! *** Model Outputs 
! *** 
!			LAI	- [-Inf;+Inf]
!			W	g [-Inf;+Inf]
!			Wf	g [-Inf;+Inf]
!			Wm	g [-Inf;+Inf]
!			W_Total	g [-Inf;+Inf]
!			Wf_Total	g [-Inf;+Inf]
!			Wm_Total	g [-Inf;+Inf]

! *** 
! *** Model Derivatives 
! *** 

! (Comments and routine interface generated by TRNSYS Simulation Studio)
!************************************************************************

!-----------------------------------------------------------------------------------------------------------------------
! This TRNSYS component skeleton was generated from the TRNSYS studio based on the user-supplied parameters, inputs, 
! outputs, and derivatives.  The user should check the component formulation carefully and add the content to transform
! the parameters, inputs and derivatives into outputs.  Remember, outputs should be the average value over the timestep
! and not the value at the end of the timestep; although in many models these are exactly the same values.  Refer to 
! existing types for examples of using advanced features inside the model (Formats, Labels etc.)
!-----------------------------------------------------------------------------------------------------------------------


      Use TrnsysConstants
      Use TrnsysFunctions

!-----------------------------------------------------------------------------------------------------------------------

!DEC$Attributes DLLexport :: Type216

!-----------------------------------------------------------------------------------------------------------------------
!Trnsys Declarations
      Implicit None

      Double Precision Timestep,Time
      Integer CurrentUnit,CurrentType


!    PARAMETERS
      DOUBLE PRECISION AREA_GH  ! Area of greenhouse [m2]
      DOUBLE PRECISION VEG_FRAC ! Vegtated frac of floor [-]
      DOUBLE PRECISION LAI_Max  ! Maximum leaf area index [m2/m2]
      DOUBLE PRECISION CROP_density ! no[plants]/m2ground

!    INPUTS
      DOUBLE PRECISION T_GH     ! Greenhouse indoor temperature [C]
      DOUBLE PRECISION C_CO2    ! Greenhouse indoor CO2 concentration [ppm]
      DOUBLE PRECISION Q_SOL_IN ! Solar radiation above canopy [W/m2]
      
!    Outputs
      DOUBLE PRECISION LAI      ! Leaf area index [m2/m2]
      DOUBLE PRECISION W        ! Weight of aboveground biomass [g/m2]
      DOUBLE PRECISION Wf       ! Weight of fruit [g/m2]
      DOUBLE PRECISION Wm       ! Wight of mature fruit [g/m2]
      DOUBLE PRECISION W_Total  ! Total Weight of aboveground biomass [g] 
      DOUBLE PRECISION Wf_Total ! Total Weight of fruit [g]
      DOUBLE PRECISION Wm_Total ! Total Wight of mature fruit [g]
      
!    Variables
      DOUBLE PRECISION Cum_Time 
      DOUBLE PRECISION Cum_PPFD 
      DOUBLE PRECISION Cum_dNdt 
      DOUBLE PRECISION Cum_Tday
      DOUBLE PRECISION Cum_Tdaytime	
      DOUBLE PRECISION Cum_ndaytime
      DOUBLE PRECISION N 
      DOUBLE PRECISION fN_h
      DOUBLE PRECISION fN
      DOUBLE PRECISION dNdt
      DOUBLE PRECISION T_day
      DOUBLE PRECISION PPFD_day
      DOUBLE PRECISION T_daytime
      DOUBLE PRECISION lambda_d
      DOUBLE PRECISION lambdas
      DOUBLE PRECISION dLAIdt_d
      DOUBLE PRECISION dLAIdt
      DOUBLE PRECISION fR_d
      DOUBLE PRECISION fR
      DOUBLE PRECISION LFmax_d  		
      DOUBLE PRECISION LFmax
      DOUBLE PRECISION PGRED_d
      DOUBLE PRECISION PGRED
      DOUBLE PRECISION Pg_d
      DOUBLE PRECISION Pg
      DOUBLE PRECISION Rm_d
      DOUBLE PRECISION Rm
      DOUBLE PRECISION GRnet_d
      DOUBLE PRECISION GRnet
      DOUBLE PRECISION fF_d
      DOUBLE PRECISION fF
      DOUBLE PRECISION g_d
      DOUBLE PRECISION g
      DOUBLE PRECISION dWfdt_d
      DOUBLE PRECISION dWfdt
      DOUBLE PRECISION dWdt_d
      DOUBLE PRECISION dWdt
      DOUBLE PRECISION Df_d 
      DOUBLE PRECISION Df
      DOUBLE PRECISION dWmdt_d
      DOUBLE PRECISION dWmdt
      DOUBLE PRECISION is_day
      DOUBLE PRECISION Tmin
      

!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Get the Global Trnsys Simulation Variables
      Time=getSimulationTime()
      Timestep=getSimulationTimeStep()
      CurrentUnit = getCurrentUnit()
      CurrentType = getCurrentType()
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Set the Version Number for This Type
      If(getIsVersionSigningTime()) Then
		Call SetTypeVersion(17)
		Return
      EndIf
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Do Any Last Call Manipulations Here
      If(getIsLastCallofSimulation()) Then
		Return
      EndIf
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
      If(getIsEndOfTimestep()) Then
		Return
      EndIf
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Do All of the "Very First Call of the Simulation Manipulations" Here
      If(getIsFirstCallofSimulation()) Then

		!Tell the TRNSYS Engine How This Type Works
		Call SetNumberofParameters(4)           !The number of parameters that the the model wants
		Call SetNumberofInputs(3)               !The number of inputs that the the model wants
		Call SetNumberofDerivatives(0)          !The number of derivatives that the the model wants
		Call SetNumberofOutputs(7)              !The number of outputs that the the model produces
		Call SetIterationMode(1)                !An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
		Call SetNumberStoredVariables(0,12)     !The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
		Call SetNumberofDiscreteControls(0)     !The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)

		Return

      EndIf
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Do All of the First Timestep Manipulations Here - There Are No Iterations at the Intial Time
      If (getIsStartTime()) Then
        
        AREA_GH = getParameterValue(1)
        VEG_FRAC = getParameterValue(2)
        LAI_Max = getParameterValue(3)
        CROP_density = getParameterValue(4)


        T_GH = GetInputValue(1)
        C_CO2 = GetInputValue(2)
        Q_SOL_IN = GetInputValue(3)
      

        !Here we validate the parameters somewhat
        IF (AREA_GH < 0.) CALL foundBadParameter(1, 'Fatal', &
        'The greenhouse area must be positive')
        IF ((VEG_FRAC < 0.) .or. (VEG_FRAC > 1.)) CALL foundBadParameter(2, 'Fatal', &
        'The fraction of vegatated floor area must be between 0 and 1')
        IF (LAI_Max < 0.) CALL foundBadParameter(1, 'Fatal', &
        'The maximum leaf area index must be positive')
        IF (CROP_density < 0.) CALL foundBadParameter(1, 'Fatal', &
        'The crop density must be positive')      
        IF (ErrorFound()) RETURN


        !Set the Initial Values of the Outputs (#,Value)
        CALL setOutputValue(1, 0.006d0)                     !LAI
        CALL setOutputValue(2, 0.28d0)                      !W
        CALL setOutputValue(3, 0.0d0)                       !Wf
        CALL setOutputValue(4, 0.0d0)                       !Wm
        CALL setOutputValue(5, 0.28d0 * AREA_GH * VEG_FRAC) !W_Total
        CALL setOutputValue(6, 0.d0 * AREA_GH * VEG_FRAC)   !Wf_Total
        CALL setOutputValue(7, 0.d0 * AREA_GH * VEG_FRAC)   !Wm_Total

        !If Needed, Set the Initial Values of the Dynamic Storage Variables (#,Value)
        !Sample Code: Call SetDynamicArrayValueThisIteration(1,20.d0)
        CALL setDynamicArrayValueThisIteration(1, 0.d0)     !Cum_Time
        CALL setDynamicArrayValueThisIteration(2, 0.d0)     !Cum_PPFD
        CALL setDynamicArrayValueThisIteration(3, 0.d0)     !Cum_dNdt
        CALL setDynamicArrayValueThisIteration(4, 0.d0)     !Cum_Tday
        CALL setDynamicArrayValueThisIteration(5, 0.d0)     !Cum_Tdaytime 
        CALL setDynamicArrayValueThisIteration(6, 0.d0)     !Cum_ndaytime
        CALL setDynamicArrayValueThisIteration(7, 6.0d0)    !N
        CALL setDynamicArrayValueThisIteration(8, 0.006d0)  !LAI
        CALL setDynamicArrayValueThisIteration(9, 0.28d0)   !W
        CALL setDynamicArrayValueThisIteration(10, 0.d0)    !Wf
        CALL setDynamicArrayValueThisIteration(11, 0.d0)    !Wm
        CALL setDynamicArrayValueThisIteration(12, 100.d0)  !Tmin
 

		Return

      EndIf
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!ReRead the Parameters if Another Unit of This Type Has Been Called Last
      If(getIsReReadParameters()) Then
		!Read in the Values of the Parameters from the Input File
        AREA_GH = getParameterValue(1)
        VEG_FRAC = getParameterValue(2)
        LAI_Max = getParameterValue(3)
        CROP_density = getParameterValue(4)		
      EndIf
!-----------------------------------------------------------------------------------------------------------------------

!Read the Inputs to the model for current timestep
      T_GH = GetInputValue(1)
      C_CO2 = GetInputValue(2)
      Q_SOL_IN = GetInputValue(3)

! Read the outpu from prev timestep as the intiial value for this timestep
       Cum_Time = GetDynamicArrayValueLastTimestep(1)
       Cum_PPFD = GetDynamicArrayValueLastTimestep(2)
       Cum_dNdt = GetDynamicArrayValueLastTimestep(3)
       Cum_Tday = GetDynamicArrayValueLastTimestep(4)
       Cum_Tdaytime = GetDynamicArrayValueLastTimestep(5)	
       Cum_ndaytime = GetDynamicArrayValueLastTimestep(6)
       N = GetDynamicArrayValueLastTimestep(7)
       LAI = GetDynamicArrayValueLastTimestep(8)
       W = GetDynamicArrayValueLastTimestep(9)
       Wf = GetDynamicArrayValueLastTimestep(10)
       Wm = GetDynamicArrayValueLastTimestep(11)
       Tmin = GetDynamicArrayValueLastTimestep(12)

!-----------------------------------------------------------------------------------------------------------------------
!Get some derived values

! Calculating the daily integral of PPFD, daily average temperature, and daily average day time temperature
	Cum_Time = Cum_Time + Timestep
	Cum_PPFD = Cum_PPFD + Q_SOL_IN * 2.3d0 ! Considering PPFD = Incoming short wave radiation [W/m2] * frac_PAR [default value 0.5] * J_to_mol_conv_fact [default value 4.6] https://search.r-project.org/CRAN/refmans/bigleaf/html/Rg.to.PPFD.html
	Cum_Tday = Cum_Tday + T_GH
    
    If (Q_SOL_IN>0.0d0)Then
        is_day = 1.0d0
    Else
        is_day = 0.0d0
    EndIf    
    
    Tmin = MIN(T_GH, Tmin)
	Cum_Tdaytime = Cum_Tdaytime + is_day * T_GH ! Total day time temperature combining all time step every 24 hour
	Cum_ndaytime = Cum_ndaytime + is_day        ! Cumulative countdown on the number of timesteps it was found day considering an above 0 value of incoming solar radiation
	fN_h = fN(T_GH)                             ! Calculating hourly fN value
	Cum_dNdt = Cum_dNdt + dNdt(fN_h, Timestep)  ! hourly dNdT is multiplied with hourly fN value and integrated daily to finally come up with a daily dNdT i.e Cum_dNdT 

	IF (Cum_Time > 24.0d0) Then ! Update all principal variables every 24 hour, outer loop
    	T_day = Cum_Tday / (24.0d0/Timestep)
		PPFD_day = Cum_PPFD / (24.0d0/Timestep)
		T_daytime = Cum_Tdaytime / Cum_ndaytime
		
		! dLAI/dt
		lambda_d = lambdas(T_day)
		dLAIdt_d = dLAIdt(LAI, LAI_Max, CROP_density, N, lambda_d, Cum_dNdt)

		! dWfdt
		fR_d = fR(N)
		LFmax_d = LFmax(C_CO2)
		PGRED_d = PGRED(T_day)
    	Pg_d = Pg(LFmax_d, PGRED_d, PPFD_day, LAI)
    	Rm_d = Rm(T_day, W, Wm)
    	GRnet_d = GRnet(Pg_d, Rm_d, fR_d)
    	fF_d = fF(T_day,Tmin)
    	g_d = g(T_daytime)
    	dWfdt_d = dWfdt(GRnet_d, N, g_d)
				
		! dWdt
		dWdt_d = dWdt(LAI,LAI_Max, dWfdt_d, GRnet_d, CROP_density, Cum_dNdt)

		! dWmdt
		Df_d = Df(T_day)
		dWmdt_d = dWmdt(Df_d, Wf, Wm, N)

		! Update state variables every 24 hours
		N = N + Cum_dNdt
		LAI = LAI + dLAIdt_d
		Wf = Wf + dWfdt_d
		W = W + dWdt_d
		Wm = Wm + dWmdt_d


		! Reset Counter clock for another 24 hour time calculation
		Cum_Time = 0.0d0
        Cum_PPFD = 0.0d0
        Cum_dNdt = 0.0d0
        Cum_Tday = 0.0d0
        Cum_Tdaytime = 0.0d0
        Cum_ndaytime = 0.0d0
        Tmin = 100.0d0 ! Reset Tmin to a high value
	ELSE

		! hold state variables as is
    	N = N
		LAI = LAI
		Wf = Wf
		W = W
		Wm = Wm
    ENDIF
!-----------------------------------------------------------------------------------------------------------------------
!Set the Outputs from this Model (#,Value)
      CALL setOutputValue(1,LAI) ! Leaf Area index
      CALL setOutputValue(2,W)   !Weight of aboveground biomass [g/m2]
      CALL setOutputValue(3, Wf)  !Weight of fruit [g/m2]
      CALL setOutputValue(4, Wm)  !Wight of mature fruit [g/m2]
      CALL setOutputValue(5, W * AREA_GH * VEG_FRAC)  ! Total Weight of aboveground biomass [g]
      CALL setOutputValue(6, Wf * AREA_GH * VEG_FRAC) ! Total Weight of fruit [g]
      CALL setOutputValue(7, Wm * AREA_GH * VEG_FRAC)  ! Total Weight of mature fruit [g]

!-----------------------------------------------------------------------------------------------------------------------
!Now need to update the state vector:

      CALL setDynamicArrayValueThisIteration(1, Cum_Time)
      CALL setDynamicArrayValueThisIteration(2, Cum_PPFD)
      CALL setDynamicArrayValueThisIteration(3, Cum_dNdt)
      CALL setDynamicArrayValueThisIteration(4, Cum_Tday)
      CALL setDynamicArrayValueThisIteration(5, Cum_Tdaytime)
      CALL setDynamicArrayValueThisIteration(6, Cum_ndaytime)
      CALL setDynamicArrayValueThisIteration(7, N)
      CALL setDynamicArrayValueThisIteration(8, LAI)
      CALL setDynamicArrayValueThisIteration(9, W)
      CALL setDynamicArrayValueThisIteration(10, Wf)
      CALL setDynamicArrayValueThisIteration(11, Wm)
      CALL setDynamicArrayValueThisIteration(12, Tmin)    !Tmin

!-----------------------------------------------------------------------------------------------------------------------
!Now return control to the TRNSYS kernel: 
      Return
  End
!-----------------------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!Function to modify node development rate as a function of hourly temperature (-), Taken from Heuvelink(1994) & Jones(1991)

      DOUBLE PRECISION FUNCTION fN(T)
        DOUBLE PRECISION T
		IF (T > 12.d0 .and. T <= 28.d0) then
    			fN = 1.0d0 + 0.0281d0 * (T - 28.d0)
		ELSEIF (T > 28.d0 .and. T < 50.d0) then
    			fN = 1.0d0 - 0.0455d0 * (T - 28.d0)
		ELSE
    			fN = 0.d0
		ENDIF

        RETURN 
      END

!Function to calculate the node derivative (dNdt) 	
      DOUBLE PRECISION FUNCTION dNdt(fN_hour, Timestep)
	    DOUBLE PRECISION Nm, fN_hour, Timestep
		Nm = (0.5d0/24.0d0) * Timestep ! Nm value per timestep
		dNdt = Nm * fN_hour ! Calculating dNdt for every timestep
	    RETURN
      END

!Function to calculate the rate of development or aging of leaves at temperature T
      DOUBLE PRECISION FUNCTION lambdas(Td)
	    DOUBLE PRECISION Td
		lambdas = 1.0d0
	    RETURN
      END

! Function to calculate the derivative of leaf area index (dLAIdt)
      DOUBLE PRECISION FUNCTION dLAIdt(LAI, LAImax, dens, N, lambda_d, dNdt_d)
	    DOUBLE PRECISION LAI, dens, N, lambda_d, dNdt_d, sigma, beta, Nb, LAImax, a 
 	        sigma = 0.038d0 ! Maximum leaf area expansion per node, coefficient in expolinear equation; Jones(1999)
    	    beta = 0.169d0 ! Coefficient in expolinear equation; Jones(1999)
    	    Nb = 16.0d0 ! Coefficient in expolinear equation, projection of linear segment of LAI vs N to horizontal axis; Jones(1999)
	    
            IF (LAI > LAImax) Then
    		    dLAIdt = 0.0d0
	        ELSE
    		    a = exp(beta * (N - Nb))
    		    dLAIdt = (dens * sigma * lambda_d * a * dNdt_d) / (1.0d0 + a)
            ENDIF
	    RETURN
      END

!Function to calculate the above ground biomass accumulation (dWdt)
      DOUBLE PRECISION FUNCTION dWdt(LAI,LAImax, dWfdt_d, GRnet_daily, dens, dNdt_d)
	    DOUBLE PRECISION  LAI, dWfdt_d, GRnet_d, LAImax, p1, Vmax, dens, dNdt_d, a, b, GRnet_daily  

	    IF (LAI >= LAImax) then
    		p1 = 2.0d0 ! P Jones(1999)
	    ELSE
    		p1 = 0.0d0
	    ENDIF

	    Vmax = 8.0d0 ! P Jones(1999)

	    a = dWfdt_d + (Vmax - p1) * dens * dNdt_d
	    b = GRnet_daily - p1 * dens * dNdt_d
	    dWdt = min(a, b)
	    RETURN
       END

!Function to calculate The rate of development or aging of fruit at temperature T; Jones(1991)
       DOUBLE PRECISION FUNCTION Df(T)
    	    IF (T > 9.0d0 .and. T <= 28.0d0) then
            	Df = 0.0017d0 * T - 0.015d0
            ELSEIF (T > 28.0d0 .and. T <= 35.0d0) then
        	Df = 0.032d0
            ELSE
        	Df = 0.0d0
    	    ENDIF
	RETURN
       END

!Function to calculate the mature fruit biomass accumulation rate (dWmdt)
	DOUBLE PRECISION FUNCTION dWmdt(Df_d, Wf, Wm, N)
		DOUBLE PRECISION Df_d, Wf, Wm, N, NFF, kF
		NFF = 22.0d0 ! P Jones(1999)
    	kF = 5.0d0 ! P Jones(1999)
		IF (N <= NFF + kF) then
    			dWmdt = 0.0d0
		ELSE
    			dWmdt = Df_d * (Wf - Wm)
		ENDIF
	 RETURN
	END

!root phenology-dependent fraction; Jones(1991)
	DOUBLE PRECISION FUNCTION fR(N)
		DOUBLE PRECISION N
	!IF (N >= 30.0d0) Then
    !		fR = 0.07d0
	!ELSE
    !		fR = -0.0046d0 * N + 0.2034d0
	!ENDIF
        
        fR = max(0.02, 0.18 - 0.0032*N)
        
	 RETURN
	END

!Maximum leaf photosynthetic rate; Jones(1991)
	DOUBLE PRECISION FUNCTION LFmax(CO2)
		DOUBLE PRECISION tau, CO2
		tau = 0.0693d0 ! carbon dioxide use efficiency; Jones(1991)
		LFmax = tau * CO2
	 RETURN
	END

!function to modify Pg under suboptimal daytime temperatures; Jones(1991)
	DOUBLE PRECISION FUNCTION PGRED(T)
		IF (T > 0.0d0 .and. T <= 12.0d0) then
    			PGRED = 1.0d0 / 12.0d0 * T
		ELSEIF (T > 12.0d0 .and. T < 35.0d0) then
    			PGRED = 1.0d0
		ELSE
    			PGRED = 0.0d0
		ENDIF
	 RETURN
	END

! function for calculating photosynthesis rate
	DOUBLE PRECISION FUNCTION Pg(LFmax_d, PGRED_d, PPFD_d, LAI)
		DOUBLE PRECISION LFmax_d, PGRED_d, PPFD_d, LAI, D, K, m, Qe, a, b
		D = 2.593d0   ! coefficient to convert Pg from CO2 to CH2O; Jones(1991)
		K = 0.58d0    ! light extinction coefficient; Jones(1991)
		m = 0.1d0     ! leaf light transmission coefficient; Jones(1991)
		Qe = 0.0645d0 ! leaf quantum efficiency; Jones(1991)

		a = D * LFmax_d * PGRED_d / K
		b = log(((1.0d0-m) * LFmax_d + Qe * K * PPFD_d) / ((1.0d0-m) * LFmax_d + Qe * K * PPFD_d * exp(-1.0d0 * K * LAI)))
		Pg = a * b
	 RETURN
	END

! function for calculating daily integral of maintenance respiration
	DOUBLE PRECISION FUNCTION Rm(T, W, Wm)
		DOUBLE PRECISION T, W, Wm, Q10, rm_c
		Q10 = 1.40d0 ! P Jones(1991)
		rm_c = 0.016d0 ! P Jones(1999)
		Rm = (Q10 ** ((T-20.0d0)/10.0d0)) * rm_c * (W - Wm)
	 RETURN
	END

! function to calculate net above ground growth rate
	DOUBLE PRECISION FUNCTION GRnet(Pg_d, Rm_d, fR_d)
    		DOUBLE PRECISION Pg_d, Rm_d, fR_d, E
   	 	    E = 0.717d0 ! Growth efficiency, ratio of biomass to photosynthate available for growth (0.70 g dry weight g–1CH2O) convert efficiency; Dimokas(2009)
    		GRnet = max(0.0d0, E * (Pg_d - Rm_d) * (1.0d0 - fR_d))
	 RETURN
	END

! FUNCTION TO MODIFY PARTITIONING TO FRUIT VS. AVERAGE DAILY TEMPERATURE
	DOUBLE PRECISION FUNCTION fF(Td,Tm)
    		DOUBLE PRECISION Td, Tm

    		! IF (Td > 8.0d0 .and. Td <= 28.0d0) then
        	!	fF = 0.0017d0 * Td - 0.0147d0
    		! ELSEIF (Td > 28.0d0) then
        		fF = 0.032d0
    		! ELSE
        		fF = 0.0d0
    		!ENDIF
            
            fF = MAX(0.0d0,MAX(1.0d0,0.0625d0 * (Td-Tm)))
            
	 RETURN
	END

! Function to reduce growth due to high daytime temperature
	DOUBLE PRECISION FUNCTION g(T_daytime)
    		DOUBLE PRECISION T_daytime, T_CRIT
    		T_CRIT = 24.4d0 ! mean daytime temperature above which fruits abortion starts; Jones(1999)
    		IF (T_daytime <= T_CRIT) then
        		g = 0.0d0
    		ELSE
        		g = max(0.0d0, 1.0d0 - 0.154d0 * (T_daytime - T_CRIT))
    		ENDIF
	 RETURN
	END

! function to calculate fruit dry weight increase rate
	DOUBLE PRECISION FUNCTION dWfdt(GRnet_day, Node, g_day) ! THIS FUNCTION HAS BEEN CORRECTED
    		DOUBLE PRECISION Node, NFF, alpha_F, v, fF_day, GRnet_day, g_day
    		NFF = 22.0d0 ! Nodes per plant when the first fruit appears; Jones(1999)
    		alpha_F = 0.80d0 ! Maximum partitioning of new growth to fruit; Jones(1999)
    		v = 0.135d0 ! Transition coefficient between vegetative and full fruit growth; Jones(1999)
    		fF_day = 0.5d0 ! ORIGINAL
    		IF (Node <= NFF) then
        		dWfdt = 0.0d0
    		ELSE
        		dWfdt = GRnet_day * alpha_F * fF_day * (1.0d0 - exp(v * (NFF - Node))) * g_day
    		ENDIF
	 RETURN
	END