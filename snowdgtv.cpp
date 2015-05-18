//copied from snowdgtv.f90
//_# 6.6.13
# include "uebpgdecls.h"
//**********************************************************************************************
//
//  Copyright (C) 2012  David Tarboton, Utah State University, dtarb@usu.edu.  http://www.engineering.usu.edu/dtarb
//
//  This file is part of UEB.
//
//    UEB is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    UEB is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    A copy of the GNU General Public License is included in the file gpl.txt.
//    This is also available at: http://www.gnu.org/licenses/.
//
//  ifyou wish to use or incorporate this program (or parts of it) into 
//  other software that does not meet the GNU General Public License 
//  conditions contact the author to request permission.
//  David G. Tarboton  
//  Utah State University 
//  8200 Old Main Hill 
//  Logan, UT 84322-8200 
//  USA 
//  http://www.engineering.usu.edu/dtarb/ 
//  email:  dtarb@usu.edu 
//
//**********************************************************************************************  
// Definitions
//  dt  time step in hours
//  nt number of time steps
// input  -- input forcing
//	  input(1,*) air temperature (C)
//	  input(2,*) precipitation (m/hr)
//	  input(3,*) wind speed  (m/s)
//	  input(4,*) relative humidity (on a 0-1 scale)
//	  input(5,*) incoming short wave  (kJ/m^2/h)
//	  input(6,*) net radiation  (kJ/m^2/h)
//	  input(7,*) Cosine of Zenith angle  
// SITEV -- site variables
//        site variables (1-5)
//        sitev(1)  drift factor  (No detailed information give 1)
//        sitev(2)  air pressure (Pa)
//        sitev(3) ground heat flux  Kj/m^2/hr (3.6 = 1 W/m^2)
//        sitev(4) albedo extinction parameter (m)
// STATEV
//        statev(1)  Snow Energy Content  (KJ/m^2)
//        statev(2)  Snow Water Equivalent (m) relative to T = 0 C solid phase
//        statev(4)  Canopy Snow Water Equivalent (m) relative to T = 0 C solid phase
//        statev(3)  Dimensionless age of snow surface (or albedo - depending on flag 4)
//        statev(5)  Refreezing depth (m) used as refdepth
//        statev(6)  Refreezing depth (m) used as totalrefdepth
//  totalrefdepth is a misnomer.  These are the same quantity - the refreezing depth.  They are repeated because 
//  when refdepth exceeds the diurnal penetration depth it is set to 0 to tell the code to use the regular 
//  surface temperature functions while totalrefdepth is only reset to 0 when there is actual melt generated
//  or energy content becomes negative indicating freezing of all liquid phase.  This ensures that the regular 
//  surface temperature functions persist when there is liquid water present but the refreezing front has penetrated
//  to depth greater than diurnal temperature fluctuations.
//        TsPrevday(1:nstepday)   Surface temperature over the last 24 hours
// 	   TavePrevday(1:nstepday)   Depth average te

// Imperature over the last 24 hours
//
// PARAM  --  snowmelt model parameters (see below)
// iflag  -- flags 
//	   iflag(1) 0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7
//        iflag(2)        no 0 (/yes 1) printing
//        iflag(3)  unit to which to print
//        iflag(4)  how albedo calculations are done (a value 1 means albedo is calculated, otherwise statev(3) is albedo
// 	   iflag(5)  model option for surface temperature approximation 
//              1) the Simple Gradient, almost the same as Original UEB,
//              2) Simple Gradient approach with shallow snow correction. 
//              3) The classical force-restore approach.
//              4) Modified force-restore approach.
// cump,cume,cummr  -- cumulative precipitation (with df), cumulative evaporation, cumulative melt over time step in m
// outv  -- output variables 
//       outv(1)=prain   rain  m/hr
//       outv(2)=ps     snow  m/hr
//       outv(3)=a     albedo
//       outv(4)=qh    sensible heat (kJ/m2/hr) 
//       outv(5)=qe    latent heat (kJ/m2/hr) 
//       outv(6)=e     sublimation m/hr
//       outv(7)=mr    melt outflow m/hr
//       outv(8)=qm    heat advected with melt
//       outv(9)=q     Energy rate of change (kJ/m2/hr) 
//       outv(10)=fm   Mass rate of change (m/hr)
//       outv(11)=tave  Average temperature (C)
//       outv(12)=tsurf  Surface temperature C
//       outv(13)=qnet  Net Radiation (kJ/m2/hr) 
//       outv(14)=smelt   Surface melt  m/hr
//
//yjs Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow 
//yjs drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
//yjs snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes 
//yjs back. //
void SNOWUEB2(float dt, float *input, float *sitev, float *statev, float *tsprevday, float *taveprevday, int &nstepday, float *param, int *iflag,
			            float &cump, float &cumes, float &cumEc, float &cumMr, float &cumGM, float &cumEg,float *outv, float *mtime ,float atff, float cf, float *OutArr)
{
    const int niv = 8;
    int pflag, iradfl;
	//FILE *outUnit;
    //	For canopy variables	
	float  AvailableSnow, RemainingAvailableWater;
    //CHANGES TO ACCOMODATE GLACIER
    float  WGT; // WGT=WATER EQUIVALENT GLACIER THICKNESS
   

    int iTsMethod;       //yjs Add model time initialization 09/19/2000
    int windfl;
    float SWIT,SWISM, SWIR,SWIGM;

	//definition added--didn't exit in the forrtran version 5.6.13
	float refDepth, totalRefDepth, Us_old, refDepth_old;
	float RRHOI, RRHO ,RID ,rhom;
	float Ta, P, V, RH, Qsi, Qli, Qnetob, cosZen, Tave;
	//int snowdgtvariteflag = 0;
	float Ps, Pr, Alb/*previously A*/, Rhofs, S;
	//outputs from Predictor corrector
	float Tc, QHc, QEc, Ec, Qpc, Qmc, Mc, FMc, intc, Inmax, ieff, Ur, Taufb, Taufd, Qsib, Qsid, Taub, Taud, Qsnc, Qsns, Qlnc, Qlns, Rkinc, Rkinsc, Vz, Tac;
	// Just for testing (outputs from predictor corre)
    float QHs, QEs, Es, QPs, Mr, QMs, Q, FM, TSURFs,  Tavepc,  Qnet, smelt, Qlis;

	//####***TGBC_6. 5.13
    
//  Parameters
    float Tr=param[0],     //  Temperature above which all is rain [3 C],
    Ts=param[1],     //  Temperature below which all is snow [-1 C],
    Ems=param[2],    //  emmissivity of snow [nominally 0.99],
    cg =param[3],    //  Ground heat capacity [nominally 2.09 KJ/kg/C],
    z=param[4],      //  Nominal meas. height for air temp. and humidity [2m],
    zo=param[5],     //  Surface aerodynamic roughness [m],
    rho=param[6],    //  Snow Density [Nominally 450 kg/m^3],
    rhog=param[7],   //  Soil Density [nominally 1700 kg/m^3],
    lc=param[8],     //  Liquid holding capacity of snow [0.05],
    ks=param[9],    //  Snow Saturated hydraulic conductivity [20 m/hr],
    de=param[10],    //  Thermally active depth of soil [0.1 m],
    abg=param[11],   //  Bare ground albedo  [0.25],
    avo=param[12],   //  Visual new snow albedo [0.95],
    anir0=param[13], //  NIR new snow albedo [0.65],
	lans= param[14], //  the thermal conductivity of fresh [dry], snow [0.0576 kJ/m/k/hr], [Vinod// 0.36 Ref: Snow and Climate :Richard L Armstrong and Eric Brun ], 
	lang= param[15], //  the thermal conductivity of soil [:9.68 kJ/m/k/hr], [TK of ice or wet soil[2.22~ 3.48W/m/k],:Vinod],
	wlf= param[16],  //  Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr], 
	rd1= param[17],  //  Apmlitude correction coefficient of heat conduction [1],
    fstab=param[18], //  Stability correction control parameter 0 = no corrections, 1 = full corrections
	Tref=param[19],  //  Reference temperature of soil layer in ground heat calculation input
	dNewS=param[20], //  The threshold depth of for new snow [0.001 m],
 	gsurf=param[21], //  The fraction of surface melt that runs off [e.g. from a glacier],

// 7 Parameters added for canopy
	EmC      = param[22],			// Emissivity of canopy 
	alpha    = param[23],			// Scattering coefficient for solar radiation
	alphaL   = param[24],		    // Scattering coefficient for long wave radiation
	G        = param[25],            // leaf orientation with respect to zenith angle
    Uc       = param[26],		    // Unloading rate coefficient [Per hour], [Hedstrom and pomeroy, 1998], 
	as       = param[27],			// Fraction of extraterrestaial radiation on cloudy day,Shuttleworth [1993],  
	bs       = param[28],		  	// [as+bs],:Fraction of extraterrestaial radiation on clear day, Shuttleworth [1993], 
	Lambda   = param[29],		    // Ratio of direct atm radiation to diffuse,worked out from Dingman [1993],
	Rimax    = param[30],            // Maximum value of Richardsion number for stability corretion
	Wcoeff   = param[31];            // Wind decay coefficient for the forest
  
  // Some initializations
  Rkinc=0.0;
  Tac=0.0;
//  State Variables - These serve as initial conditions
    float Us = statev[0],				    	// Snow Energy Content  [KJ/m^2];
    Ws = statev[1],						// Snow Water Equivalent [m]; relative to T = 0 C solid phase
    Wc = statev[3];                      // Added for canopy

	/*if(snowdgtvariteflag == 1)
	    for (int ki=0;ki<4; ki++)
		cout<<statev[ki]<<endl;
	cout<<Ws<<" "<<Wc<<endl;
	*/
    if(Us <= 0.0)
	{
		refDepth = 0.0;
	    totalRefDepth = 0.0;
	}
	else
	{
	  refDepth      = statev[4];
	  totalRefDepth = statev[5];
	}

//  Save Old Value 07/23/01     
	Us_old        =  Us;
	refDepth_old  =  refDepth;

//  Site Variables
    float dF = sitev[0],     //  Drift factor
    APr= sitev[1],     //  Atmospheric Pressure [Pa],
    qg = sitev[2],     //  Ground heat flux [KJ/m^2/hr],  This is more logically an
                        //  input variable,but is put here because it is never known at each
                        //  time step. Usually it will be assigned a value 0.
    Aep= sitev[3],     //  Albedo extinction parameter to smooth
                        //  transition of albedo when snow is shallow. Depends on Veg. height [m],
// 7 Site Variables added for canopy
	Cc    = sitev[4],   // Canopy Coverage
	Hcan  = sitev[5],   // Canopy height  
	LAI   = sitev[6],   // Leaf Area Index  
	Sbar  = sitev[7],  // Maximum snow load held per unit branch area[Kg/m2 for Pine],
	Ycage = sitev[8];  // Parameters for wind speed transformation
                                 // Ycage=1 for young pine   Should be in parameters
                                 // Ycage=2 for Leafed deciduous
                                 // Ycage=3 for Old pine with logn stems (Paw U and Meyers, 1987) 
			        		   // Requires for wind speed transformation		
//  Control flags
   iradfl= iflag[0];
   pflag = iflag[1];
   //ounit = iflag[2];      //iflag(4) albedo caculation
	                 
   iTsMethod = iflag[4];  // the method to approximate the surface temperature
						  // 1 normal old snow melt model
						  // 2 revised direct gradient method (New method for ke) and qe
						  // 3 force restore approach
						  // 4 modified force restore approach
    windfl = iflag[5];

//  Model time step information
	float yy = mtime[0],
	mm = mtime[1],
	dd = mtime[2],
	hr = mtime[3];
//**************************************************
//   Different densities and their ratio   
    RRHOI= Rho_i/Rho_w;
	RRHO = rho/Rho_w;
    RID  = 1.0/RRHO-1.0/RRHOI;    //used to compute melt water flux (Fmelt)
	rhom = lc*rho;

	// for Gracier melting calculation
    if((sitev[9] == 0) || (sitev[9] == 3))
        WGT=0.0;
	else
        WGT=1.0;
    Ws = Ws+WGT;
// loop over nt removed here because UEBGrid makes single call to this function
//#_5.6.13
//   DO 2 i = 1,nt
//   Input variables
        Ta =input[0];			// Air temperature input [Degrees C];
        P  =input[1];			// Precipitation rate input [m/hr];
        V  =input[2];			// Wind Speed [m/s];
        RH =input[3];			// Relative humidity [fraction 0-1];
//   if[iradfl.eq.0];THEN		// input is incoming short and longwave
        Qsi=input[4];			// Incoming shortwave radiation [KJ/m^2/hr];
        Qli=input[5];			// Incoming longwave radiation [KJ/m^2/hr];
//   ELSE
        Qnetob = input[6]; // Net allwave radiation [KJ/m^2/hr];
//   ENDIF
        cosZen = input[7];   // Cos[angle between direct sunlight and surface normal];

//         Representative value over time step used in albedo calculation.  
//         We neglect the difference between direct and diffuse albedo.
//DGT Daily average temperatures handled internally so that multiple time steps will work
	  Tssk_ave = daily_ave(tsprevday, nstepday, -100.0) +    T_k; // (K)  Surface temperature average over last 24 hours
	  Tsavek_ave = daily_ave(taveprevday, nstepday, -100.0) +    T_k; // (K)  Depth averaged temperature average over last 24 hours

	//  if(snowdgtvariteflag == 1)
	  //cout<<Tssk_ave<<" "<<Tsavek_ave<<endl;

      Tssk_old = tsprevday[nstepday-1] +    T_k; // (C) Surface temperature from previous time step
      Tsavek_old = taveprevday[nstepday-1] +    T_k; // (C) Average temperature from previous time step
//   ifany of these variables are out of range due to any problem set them back to freezing
	  if(Tssk_old < 0)
	  {
	    if(snowdgtvariteflag == 1)		
	        cout<<"Invalid previous time step surface temperature set to 273 K"<<Tssk_old<<endl;
		Tssk_old = T_k;  
	  }
	  if(Tsavek_old  < 0)
	  {
	    if(snowdgtvariteflag == 1)		
	        cout<<"Invalid previous time step average temperature  set to 273 K"<<Tsavek_old<<endl;
		Tsavek_old  = T_k;  
	  }
	  if(Tssk_ave  < 0)
	  {
	    if(snowdgtvariteflag == 1)		
	        cout<<"Invalid last 24 hr average surface temperature  set to 273 K"<<Tssk_ave<<endl;
		Tssk_ave  = T_k;  
	  }
	  if(Tsavek_ave  < 0)
	  {
	    if(snowdgtvariteflag == 1)		
	        cout<<"Invalid last 24 hr average temperature set to 273 K"<<Tsavek_ave<<endl;
		Tsavek_ave  = T_k;  
	  }	 
//  Separate rain and snow      
        Ps = PARTSNOW(P,Ta,Tr,Ts);
        Pr = P - Ps;
//  Increase precipitation as snow by drift multiplication factor
        Ps = Ps * dF;
//  Calculate Albedo	
        // for Gracier melting calculation        
        if(iflag[3] == 1)
//            if(SITEV(10) .EQ. 0 .OR. SITEV(10) .EQ. 3)THEN
	          Alb = ALBEDO(statev[2], cosZen, (Ws-WGT)/RRHO, Aep, abg, avo, anir0);      // Snow depth (Ws/RRho)
//	        ELSE                                                         // Use of this albedo throughout time step neglects the
//                A=albedo(statev(3),cosZen,(Ws-WGT)/RRHO,aep,abg,avo,anir0)// changes due to new snow within a time step.o
//            ENDIF
		else
            Alb = statev[2];   

//************  Maximum Interception, Literature  ************************
//  This model assumes that rainfall interception and thoroughfall occurs in a similar way that snowfall interception and unloading does. 
//  However this model is not suitable to model the rainfall interception and throughfall
		
	Rhofs = 67.92 + 51.25 * exp(Ta/2.59);		// Density of fresh snow from air temp, Army Corps 1956.
	S     = Sbar * (0.27 + 46/Rhofs );			// canopy holding capacity depends on density of snow (Headstrom and Pomeroy (1998))
	Inmax = S * LAI;							// kg/m2 =lit/m2 = .001m3/m2=.001 m =1/1000 m ()
	if(Cc > 0)

	   Inmax = Inmax/1000 * Cc; 					// convert to m for per sq m of canopy
	else
	  Inmax = 0;	
	//if(snowdgtvariteflag == 1)
	  //cout<<LAI<<" "<<Cc<<" "<<" "<<Sbar<<" "<<Rhofs<<" "<<S<<" "<<Inmax<<endl;

                   
//   Call predictor corrector subroutine to do all the work for Canopy
	 PREDICORRc(Us,Ws,Wc,Alb,dt,RID,P,Pr,Ps,Ta,V,RH,Qsi,Qli,atff, cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,  
                  Tc,QHc,QEc,Ec,Qpc,Qmc,Mc,FMc,intc,Inmax,ieff,Ur, cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc,
				    Rkinsc,Vz,Tac, QHs, QEs, Es, QPs, Mr, QMs, Q, FM, TSURFs,  Tavepc,  Qnet, refDepth, totalRefDepth, smelt, gsurf,Qlis); 
//  DGT 4/2/05   Despite all checks in predicor It can (and does) occur that  we still have Us so large that it results in tavepc greater than 0, which implies that all the 
//   snow is liquid.  In these cases - just force the snow to disappear and add the energy involved to Qm. 
     
	 Tave = TAVG(Us,Ws,Rho_w,C_s,T_0,rhog,de,cg,H_f); //  this call necessary to keep track of average internal temperature used in some surface energy algorithms.
	 if(snowdgtvariteflag == 1)
	      cout<<Tave<<endl;
	 if(Tave > 0)   //  all is liquid so snow must disappear
	 {
		Mr = Mr+Ws/dt;
		QMs = QMs+Ws/dt*Rho_w*H_f;
		Q  = Q-Ws/dt*Rho_w*H_f;
		Us = Us-Ws/dt*Rho_w*H_f;
		Ws = 0.0;
	 }
//  since Us, Ws was changed need T_0 re-evaluate Tave
    Tave  = TAVG(Us,Ws,Rho_w,C_s,T_0,rhog,de,cg,H_f);   
	//if(snowdgtvariteflag == 1)cout<<Tave<<endl;
// DGT 7/25/05  T_0 guard against unreasonable Us when there is no snow do not allow bulk temperature T_0 go above 10 C
    if(Tave > 10)
	{
         Us=rhog*de*cg*10.0;
         Tave = 10.0;		 
	}      
// Update snow surface age based on snowfall in time step
    if(iflag[3] ==1)
		AGESN(statev[2],dt,Ps,TSURFs,T_k,dNewS);        
	  // if(snowdgtvariteflag == 1)
	     //cout<<statev[2]<<endl;

	//  2013 - Introduced glacier substrate into model
//  Partition melt outflow into the portion that is from rain, snow and glacier
//  Inputs are condensation (C=max(-Es,0)), precipitation as rain (pr) and precipitation as snow (Ps)
//  Outputs are melt (Mr) and sublimation (Sub=Max(Es,0))
//  DW is change in snow water equivalent (positive if decrease)
//  DG is change in glacier (positive if decrease)
//  Mass balance
//  DW+DG=Mr+Sub-Pr-Ps-C
//  In the case of snow left on the ground/glacier - any rain is assumed to have added to the snow
    float SnowLeft = findMax(Ws-WGT,0.0);
    float DW = (statev[1] - SnowLeft)/dt;   // change in seasonal snow rate
    // Reduction in snow storage is snow at the beginning minus the snow at the end in excess of glacier.  Can be negative if snow increases
    float DG = findMax(WGT-Ws,0.0)/dt;   //  Reduction in glacier storage rate
 
	float R, MSX, MRX, Ms, MG, Eglacier;
    if(SnowLeft > 0.)
	{
		 R=0.0; // No outflow due to rain
//  Inequalities that hold
//  0<= Ms <= Ps+Pr+C+DW
//  0<= MG <= DG
//  The slack in these inequalities is due to Sub
//  Compute Ms and MG using proportioning
		 MSX = Ps+Pr+DW+ findMax(-Es,0.0);   //  Maximum melt contrib from snow.  Note that condensation is taken to add to snow as is rain
		 MRX = MSX+DG;     //  Maximum melt overall
		// Mr is less than MRX due to sublimation losses.  Do the rain proportioning accordingly  
		 if(MRX <= 0.0)
		 {
			Ms=0.0;
			MG=0.0;
		 }
		 else 
		 {     
			Ms=Mr*MSX/MRX;
			MG=Mr*DG/MRX;
		 }
	}
    else 
	{// The case where all seasonal snow is gone and precipitation may contribute to total outflow
//  Inequalities that hold
//  0<= Ms <= Ps+C+DW
//  0<= R <= Pr
//  0<= MG <= DG
//  The slack in these inequalities is due to Sub.  Note that here Pr may contribute to R but that C is still put to snow
		  MSX = Ps+DW+ findMax(-Es,0.0);   //  Maximum melt contrib from snow.  Note that condensation is taken to add to snow as is rain
		  MRX = MSX+DG+Pr;     //  Maximum melt overall
		  if(MRX <= 0.0)
		  {
			Ms=0.0;
			R=0.0;
			MG=0.0;
		  } 
		  else
		  {
			Ms=Mr*MSX/MRX;
			R=Mr*Pr/MRX;
			MG=Mr*DG/MRX;
		  }
	}

    Eglacier=DG-MG;
    SWISM=Ms;
    SWIGM=MG;
    SWIR=R;
    SWIT=Mr;
    if(Ws < WGT)                  //  There was glacier melt
       Ws=0.0;                             //  Reset      
	else                                  //  Here ws is greater than the glacier thickness so no glacier melt
       Ws=Ws-WGT; 
            
/*        commented out  10.2.13 relacing with the above
// Partition of melt outflow between rain and snow melt for non glacier case
    if(Ws > 0)
	{
        SWIR = 0;
        SWISM = Mr;
	}
	else   //  This else case is never entered for glaciers because the entire WGT is not melted in one time step
	{
       AvailableSnow = (Ps + statev[1]/dt);
       SWISM = findMin(AvailableSnow,Mr);       // This bounds snow surface melt by the available snow - the balance is due T_0 rain
       SWIR = Mr-SWISM;   
	}   
//  Correction T_0 partition T_0 account for presence of glaciers 
   if(Ws < WGT)                  // There was glacier melt
   {
       AvailableSnow = (Ps+statev[1]/dt);
       SWISM = findMin(AvailableSnow,Mr);       //  This bounds snow surface melt by the available snow - the balance is due T_0 rain
       RemainingAvailableWater = Mr-SWISM;  //  This is the remaining melt T_0 allocate
       SWIGM = findMin(RemainingAvailableWater,(WGT-Ws)/dt);
       SWIR = Mr-SWISM-SWIGM;
       Ws = 0;                              //  Reset
   }
   else                                  //  Here ws is greater than the glacier thickness so no glacier melt
   {
       Ws = Ws - WGT;  
       SWIGM = 0;
   } 
//  SWIT=Mr+SWIGM  // DGT 2/22/13.  THis seems wrong
   SWIT = Mr;   //  The total is always Mr    */


//  accumulate for mass balance
   cump  = cump+(Ps+Pr)*dt;
   cumes  = cumes+Es*dt;                  
   cumEc = cumEc+Ec*dt;                // Evaporation from canopy
   cumMr = cumMr+Mr*dt;                // canopy melt not added
   cumGM = cumGM+SWIGM*dt;             //  Cumulative glacier melt
   cumEg= cumEg+Eglacier*dt;

//  yjs update the total depth of the refreezing depth in the snow pack according the the refreezing depth at time step and the positive energy input. 07/22/01
//  DGT's revised logic  1/13/05
      if(lc > 0 )
	  {
         if(refDepth > 0) 
	        totalRefDepth = refDepth;                           // ifrefreezing depth increases totalRefDepth keeps track
	     else
			 {       // here refDepth has gone T_0 0 somehow
	           if( (Mr > 0) || (Us > Us_old && Us > 0))              //   ifthere is melt or an increase in energy refDepth is reset 
     	            totalRefDepth = 0.0; 
		     }
	  }
	  else if ( (Mr > 0) || (Us > Us_old && Us > 0))            //   Here lc=0.  ifthere is melt or an increase in energy refDepth is reset                                                     
	               totalRefDepth =0.0;	                        //   This is likely redundant because iflc=0 then there is no meltwater T_0 refreeze

	if(totalRefDepth < 0)  
		totalRefDepth = 0.0;
//yjs update tsbackup and tavebackup
	for(int ii = 0; ii< nstepday-1; ++ii)
	{
		tsprevday[ii] = tsprevday[ii+1];
	    taveprevday[ii]= taveprevday[ii+1];
	}
	   tsprevday[nstepday-1] = TSURFs;
	   taveprevday[nstepday-1]= Tave;
        
      OutArr[0] = Us;
      OutArr[1] = Ws;
      OutArr[2] = statev[2];
      OutArr[3] = Pr;
      OutArr[4]=Ps;
      OutArr[5]=Alb;
      OutArr[6]=QHs;
      OutArr[7]=QEs;
      OutArr[8]= Es;
      OutArr[9]=SWIT;
      OutArr[10]=QMs;
      OutArr[11]=Q;
      OutArr[12]= FM;
      OutArr[13]=Tave;
      OutArr[14]=TSURFs;
      OutArr[15]=cump;
      OutArr[16]=cumes;
      OutArr[17]=cumMr;
      OutArr[18]=Qnet;
      OutArr[19]=smelt;
      OutArr[20]=refDepth;
      OutArr[21]=totalRefDepth;
      OutArr[22]=cf;
      OutArr[23]=Taufb;
      OutArr[24]=Taufd;
      OutArr[25]=Qsib;
      OutArr[26]=Qsid;
      OutArr[27]=Taub;
      OutArr[28]=Taud;
      OutArr[29]=Qsns;
      OutArr[30]=Qsnc;
      OutArr[31]=Qlns;
      OutArr[32]=Qlnc ; 
      OutArr[33]=Vz;
      OutArr[34]=Rkinsc;
      OutArr[35]=Rkinc; 
      OutArr[36]=Inmax;
      OutArr[37]=intc;
      OutArr[38]=ieff;
      OutArr[39]=Ur;
      OutArr[40]=Wc;
      OutArr[41]=Tc;
      OutArr[42]=Tac;
      OutArr[43]=QHc;
      OutArr[44]=QEc;
      OutArr[45]=Ec;
      OutArr[46]=Qpc;
      OutArr[47]=Qmc;
      OutArr[48]=Mc;
      OutArr[49]=FMc;
      OutArr[50]=SWIGM;
      OutArr[51]=SWISM;
      OutArr[52]=SWIR;       
  //2    continue
       statev[0]=Us;
       statev[1]=Ws;
	   statev[4]=refDepth;
	   statev[5]=totalRefDepth;
	   statev[3]= Wc           ;           //  Added vinod
       outv[0]=Pr;
       outv[1]=Ps;
       outv[2]=Alb;
       outv[3]=QHs;
       outv[4]=QEs;
       outv[5]=Es;
       outv[6]=Mr;
       outv[7]=QMs;
       outv[8]=Q;
       outv[9]=FM;
       outv[10]=Tave;
       outv[11]=TSURFs;
       outv[12]=Qnet;
	   outv[13]=smelt;   // dgt 5/4/04

       return;
}

//************************  Predictor CORRECTOR ***********************************
void PREDICORRc ( float &Us, float &Ws, float &Wc, float Alb, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float Qli,
				    float atff, float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime, 
          	        // Following variables are output
                   float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qmc, float &Mc, float &FMc,float &intc, float &Inmax, float &ieff, float &Ur, float &Cf, 
				float &Taufb, float &Taufd, float &Qsib, float &Qsid, float &Taub,float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &Tac,
				// Just for testing
			       float &QHs, float &QEs,float &Es, float &QPs, float &MR, float &QMs, float &Q,float &FM, float &TSURFs, float &tave, float &Qnet, float &refDepth, float &totalRefDepth, float &smelt,
				   float &gsurf, float &Qlis )
{
	float   LAI, int1, Mc1, Mr1,Apr, Cc;	
	float smeltC, Wc1, FMc1,Ur1,Ec1,Tc1,Ws1,Us1,Q1,FM1,QHs1,QEs1, Es1, smelt1, TSURFs1, QMs1, Qnet1, Wc2, Ws2, Us2, Es2;
	float Ae, rlf, r;
    float  WGT = 1.0;
	Apr  = sitev[1];
	Cc   = sitev[4];
    LAI  = sitev[6];	

    QFMc(Us,Ws,Wc,Alb,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
		   Tc,QHc,QEc,Ec,Qpc,QPs,Qmc,Mc,FMc,intc,Inmax,ieff,Ur, Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, 
		     Rkinsc,Vz,TSURFs,Tac,  tave,Qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, smelt,smeltC);
	
//     PREDICTOR FOR CANOPY SNOW
	Wc1 = Wc + dt*FMc;
	if(Wc1 < 0) 
	{
		Wc1 = 0.0;
	    FMc = (Wc1- Wc)/dt;
		Ur  = findMax (0.0, (intc-FMc-Ec-Mc));
		Ec  = findMax(0.0,(intc-FMc-Mc-Ur));
		Mc  = (intc-FMc-Ec-Ur);
	    FM  = Pr+Ps-intc+Ur+Mc-MR-Es;         // int,ur,Mc added here
	}
//   Save values so that they can be averaged for output
    FMc1 = FMc;
	int1 = intc;
	Ur1  = Ur;
	Mc1  = Mc;
	Ec1  = Ec;
    Tc1  = Tc;
//  PREDICTOR FOR SUB-CANOPY SNOW
    Ws1 = Ws + dt*FM;       
	if(Ws1 < 0) 
	{
         Ws1=0.0;
         PREHELP(Ws1,Ws,dt,FM,0.0,1.0,Ps,Pr,Es,Rho_w,H_f,Q,QMs,MR,QEs, Hne_u);
	}
    Us1 = Us + dt*Q;
    Q1  = Q;
    FM1 = FM; 
//   Save values so that they can be averaged for output
    QHs1 = QHs;
    QEs1 = QEs;
    Es1	= Es;
    Mr1	= MR;
    smelt1	=smelt;                         //cdgt 5/4/04 surface melt smelt
    QMs1	= QMs;
    TSURFs1 = TSURFs;
    Qnet1	= Qnet;

	if(snowdgtvariteflag3 == 1)
	{
		cout<<"\n Predictor: Us1, Ws1, Ts1, Q1, FM1 "<<endl;
		cout<<"   "<< Us1<<" "<< Ws1<<" "<<TSURFs1<<" "<<Q1<<" "<<FM1<<endl<<endl;
	}

    QFMc(Us1,Ws1,Wc1,Alb,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
           Tc,QHc,QEc,Ec,Qpc,QPs,Qmc,Mc,FMc,intc,Inmax,ieff,Ur, Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc,
             Rkinsc,Vz,TSURFs,Tac, tave,Qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, smelt,smeltC);

//     CORRECTOR FOR CANOPY SNOW
    Wc2 = Wc + dt/2.0*(FMc1 + FMc);
	if(Wc2 < 0) 
	{
		Wc2 = 0.0;
		FMc = (Wc2- Wc)/dt*2-FMc1;
		Ur  = findMax(0.,(intc-FMc-Ec-Mc));
		Ec  = findMax(0.,(intc-FMc-Mc-Ur));
		Mc  = (intc-FMc-Ec-Ur);
        FM  = Pr + Ps-intc+Ur+Mc-MR-Es;         // int,ur,Mc added here  
	}
 	intc = (int1 +intc)/2;
	Ur  = (Ur1+Ur)/2;
	Mc  = (Mc1+Mc)/2;
	Ec  = (Ec1+Ec)/2;                            // But Mc and Ec don't change
	Wc  =  Wc2;
	Tc  = (Tc1+ Tc)/2;
	FMc = (FMc1+FMc)/2;

//     CORRECTOR FOR SUB_CANOPY SNOW
    Ws2 = Ws + dt/2.0*(FM1 + FM);
    if(Ws2 < 0) 
	{
        Ws2 = 0.0;
        PREHELP(Ws2,Ws,dt,FM,FM1,2.0,Ps,Pr,Es,Rho_w,H_f,Q,QMs,MR,QEs,Hne_u);
	}
    Us2 = Us + dt/2.0*(Q1 + Q);

	if(snowdgtvariteflag3 == 1)
	{
		cout<<"\n Corrector: Us2, Ws2, Ts, Q, FM "<<endl;
		cout<<"   "<< Us2<<" "<< Ws2<<" "<<TSURFs<<" "<<Q<<" "<<FM<<endl<<endl;
	}
	//   iterate to convergence to enhance stability
    int niter = 1, imax = 5;
    while ((abs(Ws2-Ws1) > wtol || abs(Us2-Us1) > utol) && (niter < imax))
     {
        Ws1 = Ws2;
        Us1 = Us2;
        Wc1 = Wc2;
        QFMc(Us1,Ws1,Wc1,Alb,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,  
                 Tc,QHc,QEc,Ec,Qpc,QPs,Qmc,Mc,FMc,intc,Inmax,ieff,Ur, Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, 
				   Rkinsc,Vz,TSURFs,Tac,  tave,Qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, smelt,smeltC);
//     CORRECTOR FOR CANOPY SNOW AGAIN
	   Wc2 = Wc + dt/2.0*(FMc1 + FMc);
	   if(Wc2 < 0) 
	   {
			Wc2 = 0.0;
			FMc = (Wc2- Wc)/dt*2-FMc1;
			Ur  = findMax(0.0,(intc-FMc-Ec-Mc));
			Ec  = findMax(0.0,(intc-FMc-Mc-Ur));
			Mc  = (intc-FMc-Ec-Ur);
			FM  = Pr+Ps-intc+Ur+Mc-MR-Es;         // int,ur,Mc added here
//			FM  = FM+ (FMc-Ec)  
	   }
 	    intc = (int1 +intc)/2;
		Ur  = (Ur1+Ur)/2;
		Mc  = (Mc1+Mc)/2;
		Ec  = (Ec1+Ec)/2;                         
		Wc  =  Wc2;
		Tc  = (Tc1+ Tc)/2;
		FMc = (FMc1+FMc)/2;
//dgt 5/4/04 surface melt smelt
//    corrector again
         Ws2 = Ws + dt/2.0*(FM1 + FM);
         if(Ws2 < 0) 
		 {
           Ws2 = 0.0;
           PREHELP(Ws2,Ws,dt,FM,FM1,2.0,Ps,Pr,Es,Rho_w,H_f,Q,QMs,MR,QEs,Hne_u);
		 }
         Us2 = Us + dt/2.0*(Q1 + Q);
         niter++;
         if(niter >= imax)                        // had * steps to converge now hit it. What follows is Alb fix to numerical instability that results from
		 {										//  nonlinearity when the snowpack is shallow and melting rapidly.  If convergence does not occur when the snowpack is not melting (Alb very
												//  rare thing) I just accept the predictor corrector solution. //
												//  Modified by DGT 7/27/05 to add complete meltout condition .The quantities that this changes are w2, ub2, MR and QMs
												//  The fix first checks whether the added energy is sufficient to melt the snow completely.  ifthis is the case then the snow disappears.
												  //  In other cases we assume that the liquid fraction of the snow remains constant. This to some extent bypasses the melt outflow estimates. 
												 //  ae is added energy during the time step.
            Ae = (Q1+Q+QMs1+QMs)*0.5*dt;        //   This fix is only physically sensible under melting conditions and when ae is positive and there is snow
            if((Us > 0) && (Ae > 0) && (Ws > 0))
			{
				Es2=(Es+Es1)*0.5;			// This is the average sublimation
											// Check liquid fraction with added energy.  ifthis is greater than 1 then all snow melts .Otherwise implement Alb solution assuming that the liquid fraction remains constant
	            rlf=(Us+Ae)/(Rho_w*Ws*H_f);
			    if(rlf >= 1.0)
				{
					MR = Ws/dt+(Ps+Pr-Es2);     // Here snow disappears
				    if(MR < 0 )                 //  Force this to not go negative. This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
						MR = 0.0;              //   There is Alb resulting energy discrepancy due to Alb limitation on sublimation and latent heat flux 
					                           //This is ignored because once the snow is gone energy computations are not pertinent.
                    QMs=MR*Rho_w*H_f;
				    Ws2 = 0.0;
				    Us2 = Us + Ae - QMs*dt;
				}
			    else                           //   Determine the w/ub ratio at the beginning of the time step.  Physically liquid fraction = ub/(Rho_w*w*H_f) and since Rho_w and H_f are constants
                {                               //   keeping r=w/ub constant is equivalent to keeping liquid fraction constant.  ifliquid fraction is constant this is constant.
					  r = Ws/Us;                  //   Solve for ub2 and w2 that satisfy the three equations
													//            r=w2/ub2 
													//            ub2=ub+Ae-Rho_w*H_f*MR*dt     Energy balance the last term being the energy advected by melt
													//            w2=w+(Ps+prain-e2-MR)*dt    Mass balance 
													//   The unknowns in the above are ub2, w2 and m and the equations are linear. once the first eqn is multiplied by ub2, the solution is          
					  Us2 = (Rho_w*H_f*(Ws+(Ps+Pr-Es2)*dt)-Ae-Us)/(Rho_w*H_f*r - 1);
					  Ws2=r*Us2;
					  if(Ws2 < 0)                   // Avoid negative W 
						 Ws2 = 0.0;				
					  MR = (Ws-Ws2)/dt -Es2 + Ps +Pr; 
					  if(MR < 0 )                   // Avoid negative MR
					  {
						 MR = 0.0;
						 Ws2 = Ws+ (Ps+Pr-Es2)/dt;
						 if(Ws2 <  0)
							Ws2 = 0.0;        // This can only occur ife2 is large compared to other terms.  Setting w2=0 implicitly reduces e2. There is Alb resulting energy discrepancy due to 
											  // Alb limitation on sublimation and latent heat flux.  This is ignored because once the snow is gone energy computations are not pertinent.			
					  }
					  QMs=MR*Rho_w*H_f;
					  Us2 = Us+Ae-QMs*dt;   // redundant most of the time but recalc due to exceptions
				}
//    Check that nothing went wrong
                if(MR < 0)
				{
	               cout<<"Error - negative melt rate in snow"<<endl;
				   getchar();
				}
	            if(Ws2 <  0)
				{
	               cout<<"Error - negative w2 in snow"<<endl;	
				   getchar();
				}
			    Q = Ae/dt-QMs;
//   Now set first pass values equal to final values to fake the averages below  #$%^****take?
			    QMs1 = QMs;
			    Mr1 = MR;
			    Q1 = Q;
				
			}
		 }
		 if(snowdgtvariteflag3 == 1)
		 {
			//if( niter == 1)
			cout<<endl<<"Iteration in pred-cor: niter, Us2, Ws2, Ts, Q, FM: "<<endl;
			cout<<"  "<<(niter-1)<<" "<<Us2<<" "<<Ws2<<" "<<TSURFs<<" "<<Q<<" "<< FM<<endl;
		 }
              
     }   //while ((abs(Ws2-Ws1) > wtol || abs(Us2-Us1) > utol) && (niter < imax))   //go to 1       
    Ws = Ws2;
    Us = Us2;
    //  average values from two time steps for output.  This is done for MR and e to ensure mass balance and the others for better physical comparisons
    QHs = (QHs+QHs1)* 0.5;
    QEs = (QEs+QEs1)* 0.5;
    Es = (Es+Es1)* 0.5;
    MR = (MR+Mr1)* 0.5;
    QMs = (QMs+QMs1)* 0.5;
    TSURFs = (TSURFs+TSURFs1) * 0.5;
    Qnet = (Qnet + Qnet1) * 0.5;
    Q = (Q+Q1) * 0.5;
	smelt = (smelt+smelt1) * 0.5/(H_f*Rho_w);   //cdgt 5/4/04 surface melt smelt .convert from energy KJ/m^2/hr to water depth m/hr of melt.
	
    return;
}

//*************************************************************************
//				CALCULATE CANOPY MASS FLUX AT ANY INSTANT
void QFMc(float Us, float Ws, float Wc, float A, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float Qli,
		      float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime,
			   // Following variables are output
			     float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qps, float &Qmc, float &Mc, float &FMc, float &intc,float &Inmax, float &ieff, float &Ur, float &Cf, float &Taufb,
				 float &Taufd, float &Qsib, float &Qsid, float &Taub, float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &TSURFs, float &Tac,
				 // Just for testing
				 float &tave, float &qnet, float &QHs, float &QEs, float &Es, float &MR, float &QMs, float &Q, float &FM, float &refDepth, float &totalRefDepth, float &Smelt, float &smeltc )
{
	float   TAK, TAVEK, RHOA, rhom, fKappaS, ds, TherC, Zs, Zm,Fs, Qp, Qc1, Qc2 ;
	float    d, Z0c, Rimax, Rcastar, var_a, var_b, fkappaS, Betab,Betad, Qlis, Ess, Esc, QH, QE, E;
	float Tssk, Tck;
	//not used-but not deleted yet
	//#$%____6.11.13
	float Hcan, Ac;
//	DOUBLE PRECISION Ur
//yjs  Constant data set
//common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

//  Parameters
      Ems=param[2];    //  emmissivity of snow [nominally 0.99],
float cg =param[3],    //  Ground heat capacity [nominally 2.09 KJ/kg/C],
      z=param[4],      //  Nominal meas. height for air temp. and humidity [2m],
      zo=param[5],     //  Surface aerodynamic roughness [m],
      rho=param[6],    //  Snow Density [Nominally 450 kg/m^3],
      rhog=param[7],   //  Soil Density [nominally 1700 kg/m^3],
      lc=param[8],     //  Liquid holding capacity of snow [0.05],
      ks=param[9],    //  Snow Saturated hydraulic conductivity [160 m/hr],
      de=param[10],    //  Thermally active depth of soil [0.4 m],
      abg=param[11],   //  Bare ground albedo  [0.25],
      avo=param[12],   //  Visual new snow albedo [0.95],
      anir0=param[13], // NIR new snow albedo [0.65],
	lans= param[14], // the thermal conductivity of fresh [dry], snow [0.0576 kJ/m/k/hr],
	lang= param[15], // the thermal conductivity of soil [9.68 kJ/m/k/hr],
	wlf= param[16],  // Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr], 
	rd1= param[17],  // Apmlitude correction coefficient of heat conduction [1],
      fstab=param[18], // Stability correction control parameter 0 = no corrections, 1 = full corrections
	Tref=param[19],  //  Reference temperature of soil layer in ground heat calculation input
	dNewS=param[20], //  The threshold depth of for new snow [0.001 m],
	gsurf=param[21], //  The fraction of surface melt that runs off [e.g. from a glacier],
                      //  cdgt gsurf added for modeling surface runoff from a glacier
//   Site variables
    df    =  sitev[0],      //  Drift factor
    Apr   =  sitev[1],      //  Atmospheric Pressure [Pa],
    qg    =  sitev[2],      //  Ground heat flux [KJ/m^2/hr],  This is more logically an
	Cc    =  sitev[4],      //  Canopy Coverage
	LAI   =  sitev[6];	   //  Leaf Area Index

	//the following is passed to INTERCEPT instead of param(27)==param[26], same for cc
	 float Uc   = param[26];		    // Unloading rate coefficient (Per hour) (Hedstrom and pomeroy, 1998) 		   

//  To ensure all temperatures in kelvin
    TAK    = Ta + T_k;
    TAVEK  = tave + T_k;
    RHOA   = Apr / (Ra_g*TAK);     // Density of Air in kg/m3 ifAPr in Pa
	rhom   = lc * rho;

//   Some computations for conductivity 			   
	fKappaS = lans/(rho*C_s);                            // lans= lamda= thermao conductivity
	ds		= sqrt(2 * fKappaS /W1da_y);                   // fkappas =k= diffusivity
	TherC	= lans /(ds * rd1 *rho *C_s);					  // In QFM, related to thermal conductivity , lambda, of snow (Rs in Ori)
	Zs      = Ws * Rho_w/rho;							  // snow depth

//yjs save the old values
	float tave_old = tave;
    float EA = svpw(Ta) *RH;             // The saturation vapour pressure over water is used 
//                                   because most instruments report relative humidity relative to water.
//  Fraction of snow on canopy [Dickinson et. al (1993) eqn (50a)]
	if(Inmax == 0.)             //(When LAI is zero for open area model)
	  Fs= 0.0;
	else
	  Fs   = pow((Wc/Inmax),(2.0/3)); 
	if(Fs > 1)
		Fs=1;    

	INTERCEPT(Ta, LAI, P, Wc, dt, Inmax, Uc, Cc, ieff,Ur,intc);
	/*param[26] = Uc;
	sitev[4] = Cc;*/	                
//  Total heat advected by precipitation is 
	 Qp   = QPF(Pr,Ta, T_0, Ps,Rho_w, H_f, C_w, C_s);
//  Heat advected by precipitaion on canopy snow and sub-canopy snow is 
    if(P > 0)
	{
	    if(Wc > 0)
		{
          Qpc  = intc/P * Qp;          // Energy from precipitaiton that hits the canopy snow
          Qps  = Qp - Qpc;            // Energy advected by precipitaion that passes through canopy gap and thoroughfall from canopy
		}
		else 
		{
	      Qpc = 0.0;
          Qps = ((P - intc)/P) * Qp; 
		}
    }
	else
	{
	    Qpc = 0.0;
      	Qps = 0.0;
	}  
                  
	tave = TAVG (Us,Ws, Rho_w, C_s, T_0, rhog, de, cg, H_f);
	//cout<<tave<<endl;
//yjs as described as below, the original UEB does not consider the refreezing. in this
// change, the model will consider the refreezing effect starting from the top to the bottom.
// Here is the predictor.

//yjs the method is: refreezing can be modeled ifthe Us is greater than 0.0 without any precipiation, 
// meanwhile the total refreezing depth in the snowpack is less than the Refreezing depth times the daily damping of wet snow.
 
//       if(Wc.EQ.0.and.p.EQ.0.) THEN   // When there is no snow in the canopy
//		 Tc  =  Ta                   
//  	 ELSE
//		Tc  = 0.0            // First assume snow is melting and then ifresult is negative - calculate canopy temperature having determined snow is not melting
//	 ENDIF
//**********   To Calculate Refr Depth	    
 	if( (Us > 0) && (Ps <= 0) && (Pr <= 0) && (totalRefDepth <= rd1*ds) && (Ws > 0)) 
	{
		Qc1=  QcEst(Ws,P,T_k,Tck,V,Zm,d,Z0c,Rimax,Rcastar,Cf,Fs, Qli,Hcan,Vz,Ta,RH,Rkinsc,Qps,T_0,Ps,Qsi,atff,cosZen, Apr,TAK, EA, A, Ac, Wc,Inmax, Qnetob,iradfl,param,sitev );
		Qc2=  QcEst(Ws,P,T_k -0.01,Tck,V,Zm,d,Z0c,Rimax,Rcastar,Cf,Fs, Qli,Hcan,Vz,Ta,RH,Rkinsc,Qps,T_0,Ps,Qsi,atff,cosZen, Apr,TAK, EA, A, Ac,Wc,Inmax, Qnetob,iradfl,param,sitev );
        Grad(Qc1,Qc2,0.0,-0.01, var_a, var_b);		
	    refDepth = refDep(lans,var_a, var_b, H_f, rhom, dt, refDepth);       //refreezing depth
	}
	else
        refDepth=0.0;
 // call Temperature
    TSURFs= SRFTMPSC (Tssk, Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Qsi,atff,Cf, Qli,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
		                Qpc,Qps, Inmax, Rkinc,Rkinsc,Vz,Tc,T_k,TAK,EA,RHOA, fkappaS,rho,TherC,Fs, tave,refDepth,Smelt, smeltc);    // This will give me Tsurf, Tc, and smeltc
	//cout<<TSURFs<<endl;
    if(iradfl != 1)                                                            // if iradfl==1 Go to 13   //  To avoid these steps ifthe net radiation is input
	{
	PSOLRAD(Qsi,atff,param,Cf,Taufb,Taufd,Qsib,Qsid);    // Output variables:
	TRANSRADCAN (cosZen,sitev,param, Betab,Betad,Taub,Taud);   // Output variables:
	NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, Qsib,Qsid,param,Fs, Qsns,Qsnc ); //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation)
	NETLONGRAD(RH,Ta,TSURFs,Tc,T_k,Fs,EmC,Ems,SB_c,Cf,sitev, Qli,param, Qlis,Qlns,Qlnc );                     //  Output: Qsns,Qsnc 
	}
/*13*/Ess    = svpi(TSURFs);
	Esc    = svpi(Tc);
	TURBFLUX(Ws,Wc,A,T_k,Tc,Ta,TSURFs,RH,V,EA,P,param,sitev, d,Z0c,Vz,Rkinc,Rkinsc,Tac,Fs,Ess,Esc, QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E);	
	INTERCEPT(Ta,LAI,P,Wc,dt,Inmax,Uc,Cc, ieff,Ur,intc);    // Output variables	
// melt from the canopy			                                    	
    Qmc  =   smeltc;
	Mc   =   Qmc/(Rho_w*H_f); 
//     Ridistribution of intercepted snow subtracting in mass balancep
//     We are applying same amount as unloading is lost from the canopy from the wind redistribution
//      Rid = Ur
//*****************  UPDATING THE MASS BALANCE AT CANOPY
/*12*/	FMc = intc-Ec-Mc-Ur;

//*****************  UPDATING THE MASS and ENERGY BALANCE AT SUB-CANOPY
	MR = FMELT(Us,Rho_w,Ws,H_f,lc, rid,ks,Pr);   //yjs Add a fraction to reduce the evaporation after snow gone // MR in m/hr
    QMs = MR*Rho_w*(H_f+(tave-T_0)*C_w);           //  Include advection of
												//  meltwater/rain that is at tave so that the model does not
												// crash when there is no snow and it rains.
												//   QMs in kj/m2/hr 
//  dgt 5/4/04  Add surface melt that runs off to melt runoff and melt energy so that energy and mass balances
//   are still correct
     MR = MR + Smelt/(H_f*Rho_w)*gsurf;  
	 QMs=QMs+Smelt*gsurf;
     if(iradfl == 0) 
//          QNET = QSI*(1.0-A)+QLNET
          qnet  = Qsns+Qlns;                     // At sub canopy snow
	 else
		 qnet = Qnetob;
//      if(Ws.EQ.0.and. p.EQ.0.) THEN   // vinod added to orgi UEB
//		 Es  =  0.
// 		 MR =  0.	
//	ENDIF
       Q = qnet + Qps + qg + QHs + QEs-QMs;       
       FM =1.0*(Pr+Ps)-intc +Ur +Mc -MR -Es;                                 // int,ur,Mc added here

	return; 
}

//************************* SRFTMP () *********************************
//                COMPUTE THE SURFACE TEMPERATURE OF SNOW 
float SRFTMPSC (float &Tssk,float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Qsi,float atff,float cf,float Qli,float cosZen,
				   float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Qpsin,float Inmax,float Rkinc,float Rkinsc,float Vz,
				     float &Tc,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &Tave,float &refDepth,float &smelt,float &smeltC
			    ) 
{ 		
	float  F1, F2, F1ts,F1tc,F2ts,F2tc,J11, J12, J21, J22, delTs, delTc, SRFTMPSC_v;
	float S1num, S2num, S1denom, S2denom;         //
	float Tssk1 , Tck1 , Tss, Tck , S1, S2, Tclast, Tslast, ERc, ER;
	float Tlb, Tub, Flb, Fub;
	int iterC, ibtowrite = 0;
	float Tssk1p, Tck1p;
	float TSURFs;

		//     assumed small increament to estimate the Jocobian matrix
			   delTs = 0.01;
			   delTc = 0.01;
	//9.22.13
	float f1minf2 = 1;

//  Parameters
	Ems=param[2];    //  emmissivity of snow [nominally 0.99],
float cg =param[3],    //  Ground heat capacity [nominally 2.09 KJ/kg/C],
      z    = param[4],	   //  Nominal meas. height for air temp. and humidity [2m],
      zo   = param[5];	   //  Surface aerodynamic roughness [m],
     RHO  = param[6];   //  Snow Density [Nominally 450 kg/m^3],
float      rhog = param[7],	   //  Soil Density [nominally 1700 kg/m^3],
      lans = param[14],	   //  thermal conductivity of fresh [dry], snow [0.0576 kJ/m/k/hr],
      lang = param[15],	   //  the thermal conductivity of soil [9.68 kJ/m/k/hr],
      wlf  = param[16],	   //  Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr], 
      rd1  = param[17],	   //  Apmlitude correction coefficient of heat conduction [1],
      Apr  = sitev[1],      //  Atmospheric Pressure [Pa],     
      Cc    =   sitev[4],   //  Canopy Coverage
      LAI   = sitev[6];     // Leaf Area Index  

//   This version written on 4/23/97 by Charlie Luce solves directly the energy balance equation using Newtons method - avoiding the linearizations
//   used previously.  The derivative is evaluated numerically over the range ts to fff*ts  where fff = 0.999

      float tol_l = 0.005;                 //local tolerance
      float fff    = (273.0-tol_l)/273.0;                 // Scale the difference used for derivatives by the tolerance
      Tak    = Ta + Tk;
      float Tavek  = Tave + Tk;
	  float Qps = Qpsin;
      float Qpc  =  Qpcin;                         // store input variable locally without changing global value
      
       if((Ws <= 0) && (Qps > 0))
		   Qps = 0.0;
//******** FOR OPEN AREA MODLING ********************************
	   if((LAI == 0) || (Cc == 0))
	   {
			Tssk=Tak;
		                              //goto labl18;
	   }
	   else
	   {
	  
//************  SOLVING NON LINEAR EQUATION IN TWO DIMENSTION (newton's method using Jacobian matrix)*********   
			Tssk1   = Tak;                             // first approximation
			Tck1    = Tak;                             // first approximation  

			int niter = 0;
/*15*/      while (niter < nitermax)
            {
				niter++;
				F1 = SURFEBSC(Tssk1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
				                Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck1,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs, Tave,refDepth);   
                                                        		//yjs add three value to reflect model control changes)) 
			    F2 = SURFEBC(Tck1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
					            Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk1,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);   //yjs add three value to reflect model control changes))
			   Tssk1p = Tssk1 + delTs;
			   Tck1p = Tck1 + delTc;
			   F1ts = SURFEBSC (Tssk1p,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf, Qli,Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob, iTsMethod,mtime,
				                  Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck1, Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);
			   F1tc = SURFEBSC (Tssk1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
				                    Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck1p,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);
			   F2ts = SURFEBC (Tck1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
				                   Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk1p,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);
			   F2tc = SURFEBC (Tck1p,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
				                   Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk1,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);				                                      	 
		//      Jacobian matrix
				J11 = (F1ts-F1)/delTs; 
				J12 = (F1tc-F1)/delTc;
				J21 = (F2ts-F2)/delTs;
				J22 = (F2tc-F2)/delTc;

				//9.21.13
				S1num = F1*J22-F2*J12;
				S1denom = J12*J21-J22*J11;
				S2num = F2*J11-F1*J21;
				S2denom = J12*J21-J22*J11;
				 
				S1 = (F1*J22-F2*J12) / (J12*J21-J22*J11);  //S1num/S1denom;         //
				S2 = (F2*J11-F1*J21) / (J12*J21-J22*J11); //  S2num/S2denom;			//
				
				if(S1 > 5.)
					S1= 5.0;
				if(S1 < -5) 
					S1= -5.0;
				if(S2 > 5)
					S2 = 5.0;
				if(S2 < -5)
					S2= -5.0;
				Tssk1 = Tssk1 + S1;
				Tck1  = Tck1 + S2;
				if(snowdgtvariteflag == 1)
				{
					if (niter == 1)
					{
						cout<<"Iteration solving non-linear eqn. in 2D"<<endl;
						cout<<"F1, F2, F1ts,  F2ts, F1tc, F2tc"<<endl;      //, J11, J12, J21, J22 "<<endl;
						cout<<"S1, S2, Tssk1, Tck1"<<endl;					
						//cout<<"S1num, S2num, S1denom, S2denom"<<endl;
					}
					cout<<F1<<" "<<F2<<" "<<F1ts<<" "<<F2ts<<" "<<F1tc<<" "<<F2tc<<endl; //		cout " "<<J11<<" "<<J12<<" "<<J21<<" "<<J22<<endl;
					cout<< S1<<" "<< S2<<" "<< Tssk1<<" "<< Tck1<<endl;
					//cout<<S1num<<" "<< S2num<<" "<< S1denom<<" "<< S2denom<<endl;
				}
				 //(J12*J21 == 0) && (J22*J11 == 0)  gives NAN or infinity value for S1 and S2,  occues when Ws is zero	
				if(((J12*J21 == 0) && (J22*J11 == 0)) || S1denom == 0 || S2denom == 0)		    //S1denom == 0 || S2denom == 0 added to avoid div by 0; 9.21.13
					goto labl14;   		    
				
				
			    if((abs(S1) > 0.001) || (abs(S2) > 0.001))							
					continue;                                                    //Go TO 15       // loop back and iterate
	// When there is snow and temperature is positive,we need to iterate for again putting temp zero 
			    else if((Ws > 0) && (Tssk1 >Tk))                        // Iteration doesnot fail,
				{
    					Tssk1 = Tk;                                // Just for first assumption, doesnot estimate melt
						Tssk  = Tssk1;
						Tck   = Tck1; 
				        goto labl17;                                      //go to 17
				}
				else if((Wc > 0) && (Tck1 > Tk)) 
				{
						Tck1 = Tk;                                 // Just for first assumption, doesnot estimate melt
						Tck  = Tck1;
						Tc   = Tck-Tk;
						Tssk  = Tssk1; 
						goto labl17;                                     //go to 17      
				}
				else 
				{
					SRFTMPSC_v = Tssk1-Tk; 
					smelt    = 0.0;
					Tc       = Tck1-Tk;
					smeltC   = 0.0;
					goto labl21;                     //return SRFTMPSC_v;                       //go to 21
				}
	  
				
             }     // while niter < nitermax                                                       
	
/*    
			ELSE
	//     We consider the iteration is not complete since the values of S1 and S2 are not satisfied
			GO TO 14  // Use another method :one dimension approach
			ENDIF*/
//%&*@#!!%%******^^^^#_6.7.13 Checking these conditons doesn't seem of any significance
//beecause all lead to labed line 14 which is immediately below

	//moved here 9.17.13
		// Doubt the iteration When the temperature difference between snow surface and air is more thar 20 C  *******************
	  	if((Tssk1 < (Tak-20.)) || (Ws <= 0.0 && Tssk1 > (Tak+20.))||(Tck1 < (Tak-10.0))||(Wc <= 0.0 && Tck1 > (Tak+10.)))        
			goto labl14;

	//************  SOLVING NON LINEAR EQUATION IN ONE DIMENSTION *********************************************
	//      ignore the effect of precip advected energy on the calculation of surface temperature when there is no snow.
	//      Without this ridiculously high temperatures can result as the model tries to balance outgoing radiation with precip advected energy.   

/*14:      if ((Ws > 0) || (Tssk1 >Tk) ||(Wc > 0) || (Tck1 > Tk)){*/
labl14:			
		 Tssk = Tak;                         // first approximation		   
		 Tck  = Tak;
	
   }  // if((LAI == 0) || (Cc == 0))

labl17:
	ERc    = tol_l*2;                           // so that it does not end on first loop
	iterC  =  0;

	// Estimate Tc based on assumed snow surface temperature 
	/*13*/ 
	while ((ERc > tol_l) &&  (iterC < ncitermax))
	{
		if((LAI != 0) && (Cc != 0))
		{
			Tclast = Tck;
			Tc    = CanTemp(Tck, Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Qsi,atff,cf, Qli,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
		                   Qpc, Inmax, Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC,Fs, Tave,refDepth,smeltC);                            // Reduced parametre later 
		}

labl18:
/*18*/  ER    = tol_l*2;                           // so that it does not end on first loop // This is the start when LAI is zero
		int niter = 0;
		while((ER  > tol_l) && (niter < nitermax))
		{  
			Tslast = Tssk;			
			F1 = SURFEBSC (Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
				Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);   
			//yjs add three value to reflect model control changes))       
			F2 = SURFEBSC (fff*Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
				Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);          	//yjs add three value to reflect model control changes)) 
			
			if(F1==F2)              //9.17.13
			{
				/*snowdgtvariteflag = 1;
				snowdgtvariteflag2 = 1;
				snowdgtvariteflag3 = 1;*/
				if(snowdgtvariteflag == 1)
				{					
					//F1 += (0.01*F2); 
					cout<<"  Warning! F1 and F2 have the same value, can lead to v. large surface temp Tssk"<<endl;					
					cout<<"F1, F2, Tssk, Tslast "<<endl;
					cout<<"  "<<F1<<" "<<F2<<" "<<Tssk<<" "<<Tslast<<endl;
					getchar();
			    }				
				goto labl11;                                     //10.30.13 go to bisection if F1 ==F2 as div by zero results in #IND for Tssk				
			}			
			Tssk = Tssk - ((1.0 -fff) * Tssk * F1) / (F1 - F2);
			if(snowdgtvariteflag == 1)
			{
				if (niter == 0)
						cout<<"Surface temp iteration for soln in 1D: F1, F2, Tssk, Tslast "<<endl;
				cout<<"  "<<F1<<" "<<F2<<" "<<Tssk<<" "<<Tslast<<endl;
			}
			if(Tssk < Tak - 50) 
				goto labl11;                      //ifit looks like it is getting unstable go straight to bisection
			ER = abs(Tssk - Tslast);
			niter++;
		}        // loop back and iterate	   

labl11: 
		if(ER > tol_l)                                               //if solution has NOT converged
		{                                    //   ifstill not converged use bisection method
			Tlb = Tak - 20.0;                             // First guess at a reasonable range      
			Tub = Tak + 10.0;
			Flb = SURFEBSC (Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
				Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);
			Fub = SURFEBSC(Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
				Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);

			ibtowrite = 1;    

			if(Flb*Fub > 0)     // these are of the same sign so the range needs to be enlarged // an almost ridiculously large range - solution should be in this ifit exists
			{
				Tlb= Tak - 150.0;         
				Tub= Tak + 100.0;
				Flb = SURFEBSC (Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
					Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth); 
				Fub = SURFEBSC (Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
					Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);    
				ibtowrite = 1;
				if( Flb*Fub > 0)   // these are of the same sign so no bisection solution
				{
					//if(snowdgtvariteflag == 1)
					{
						cout<<"Bisection surface temperature solution failed with large range"<<endl;
						cout<<"datetime, model element:"<< mtime[0]<<mtime[1]<<mtime[2]<<mtime[3]<<"  "<<uebCellY<<" "<<uebCellX<<endl;   //mtime[4]<<endl;					
						cout<<" A canopy temperature of 273 K assumed"<<endl;
					}
					Tssk=Tk;
					goto labl10;
				} 
				else
				{
					if(snowdgtvariteflag == 1)
					{
						cout<<"Bisection surface temperature solution with large range"<<endl;
						cout<<"datetime, model element:"<< mtime[0]<<mtime[1]<<mtime[2]<<mtime[3]<<"  "<<uebCellY<<" "<<uebCellX<<endl;   //mtime[4]<<endl;					
						cout<<" This is not a critical problem unless it happens frequently"<<endl;
						cout<<" and solution below appears incorrect"<<endl;
					}
				}	
			}                        //else      endif
			//     Here do the bisection
			niter = log((Tub-Tlb)/tol_l)/log(2.0);   // Number of iterations needed for temperature tolerance
			//if(snowdgtvariteflag == 1) 	cout<<"Bisection iter: "<<niter<<endl;
			for (int iter =0; iter< niter; ++iter)
			{
				Tssk = 0.5*(Tub+Tlb);
				F1   = SURFEBSC(Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
					Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC,  TSURFs,Tave,refDepth);  
				//yjs add three value to reflect model control changes)) 
				if(F1 > 0 )      // This relies on the assumption (fact) that this is a monotonically decreasing function
					Tlb=Tssk;
				else
					Tub=Tssk;
			}
			if(ibtowrite == 1)
				if(snowdgtvariteflag == 1)
				{
					cout<<"Surface temperature: "<<Tssk<<endl;
				    cout<<"Energy closure: "<<F1<<endl;
				    cout<<"Iterations:"<<niter<<endl;
				}	
		}

labl10:
		Tss = Tssk - Tk;
		if(Ws > 0  && Tss > 0)                  //check if snow melting
		{	//dgt 5/4/04 surface melt smelt
			SRFTMPSC_v = 0.0;
			smelt = SURFEBSC(SRFTMPSC_v+Tk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf, Qli,Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod, 
				mtime, Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,Tave,refDepth);      
													//dgt 5/4/04 surface melt smelt is the energy not accommodated by conduction into the snow
													//  so it results in surface melt which then infiltrates and adds energy to the snowpack
													//  through refreezing
													//  This change added to capture this quantity when the model is used for a glacier where
													//  surface melt may run off rather than infiltrate.  
													//  For modeling of glaciers in Antarctica by Mike Gooseff
		}
		else 
		{
			SRFTMPSC_v = Tss; 
			smelt = 0.0;           //  No surface melt this case
		}

		if((LAI == 0)|| (Cc == 0.0))
			goto labl21;                //return SRFTMPSC_v;                                                 
		// ITERATION FOR THE NEW Tc, FROM ESTIMATED Tss
		Tssk = SRFTMPSC_v + Tk;
		Tck  = Tc+Tk;
		ERc  = abs(Tck -  Tclast);
		iterC++;
		if(snowdgtvariteflag == 1)
			cout<<" iterC: "<<iterC<<" ERc: "<<ERc<<endl;
		//go to 13                  // To estimate the new TC for the estimated Tss. loop back

	}// while ((ERc > tol_l) &&  (iterC < 10))

labl21:
	return SRFTMPSC_v;
}

//************  END OF ITERATION TO SOLVE LINEAR EQUATION IN ONE DIMENSTION ************************
//*************************************************************************************************
//                  COMPUTE THE CANOPY TEMPERATURE
//************************************************************************************************** 
float CanTemp(float& Tck,float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Qsi,float atff,float cf,float Qli,
               float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Inmax,float Rkinsc,float Vz,
			     float& Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &tave,float &refDepth,float &smeltC
			) 
{  						        	          
	float CanTemp_v;
	float Tclast, F1, F2;
	float Rkinc, Qps, TSURFs;
	float Tlb, Tub, Flb, Fub, Tc;
	int ibtowrite;
	float f1minf2 = 1;
	//  Parameters
     Ems=param[2];    //  emmissivity of snow [nominally 0.99],
float cg =param[3],    //  Ground heat capacity [nominally 2.09 KJ/kg/C],
      z    = param[4],	   //  Nominal meas. height for air temp. and humidity [2m],
      zo   = param[5];	   //  Surface aerodynamic roughness [m],
      RHO = param[6];	   //  Snow Density [Nominally 450 kg/m^3],
 float     rhog = param[7],	   //  Soil Density [nominally 1700 kg/m^3],
      lans = param[14],	   //  thermal conductivity of fresh [dry], snow [0.0576 kJ/m/k/hr],
      lang = param[15],     //  the thermal conductivity of soil [9.68 kJ/m/k/hr],
      wlf  = param[16],     //  Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr], 
      rd1  = param[17],     //  Apmlitude correction coefficient of heat conduction [1],
      Apr  = sitev[1];      //  Atmospheric Pressure [Pa],

      //tol,nitermax/0.05,20/   // changed 10 to 20 and .05 to .0005
	 float tol_l = 0.05;                 //local tolerance
      float fff    = (273-tol_l)/273.0;                 // Scale the difference used for derivatives by the tolerance
      Tak    = Ta + Tk;
      float Tavek  = tave + Tk;
      float Qpc  =  Qpcin;                          // store input variable locally without changing global value

      if((Wc <=0) && (Qpc > 0)) Qpc=0.0;//
//      ignore the effect of precip advected energy on the calculation of surface temperature when there is no snow.
//      Without this ridiculously high temperatures can result as the model tries to balance outgoing radiation with precip advected energy.    
////       Tck  = TaK                              // first approximation
       float ER    = tol_l*2;                           // so that it does not end on first loop
       int niter = 0;
       while((ER  > tol_l) && (niter < nitermax))
       {
          Tclast = Tck;
		  F1 = SURFEBC (Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
			               Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);             	//yjs add three value to reflect model control changes))
		  F2 = SURFEBC (fff*Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli,Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
			               Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA,fkappaS,RHO,TherC, TSURFs,tave,refDepth);          	//yjs add three value to reflect model control changes))
		  
		  if(F1==F2)
		  {
			  if(snowdgtvariteflag == 1)
			  {
				  f1minf2 = F1 - F2;
				  cout<<"  Warning! F1 and F2 have the same value, can lead to v. large canopy temp Tck"<<endl;
				  cout<<" F1 - F2: " <<f1minf2<<endl;
			  }
			  goto labl11;                                     //10.30.13 go to bisection if F1 == F2 as div by zero results in #IND for Tssk	
		  }
		  Tck = Tck - ((1.0 - fff) * Tck * F1) / (F1 - F2);
		  if(snowdgtvariteflag == 1)
		  {   
			  if (niter == 0)
						cout<<"Canopy temp iteration for soln in 1D: F1, F2, Tck, Tclast "<<endl;
			  cout<<"  "<<F1<<" "<<F2<<" "<<Tck<<" "<<Tclast<<endl;
		  }

		  if(Tck < Tak - 20) 
			  goto labl11;                                        //ifit looks like it is getting unstable go straight to bisection
	      ER = abs(Tck - Tclast);
          niter++;
	           //continue;                                //go to 1        // loop back and iterate
	    }
          //	if(abs(F1-F2).GT.100000 ) goto 11 
labl11:
	     if(ER > tol_l)                               // If the solution has converged, skip the following
		 {
			   //   ifstill not converged use bisection method
			  Tlb = Tak - 20.0;                             // First guess at a reasonable range                 
			  Tub = Tak + 10.0;
			  Flb =  SURFEBC (Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
								  Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);
			  Fub = SURFEBC (Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
								  Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);
			  ibtowrite = 0;

			  if(Flb*Fub > 0)       // these are of the same sign so the range needs to be enlarged // an almost ridiculously large range - solution should be in this ifit exists
			  {
				  Tlb= Tak - 150.0;          
				  Tub= Tak + 100.0;
 				  Flb = SURFEBC (Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
								   Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);
				  Fub = SURFEBC (Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime, 
									Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);
				  ibtowrite = 1;
				  if(Flb*Fub > 0)   // these are of the same sign so no bisection solution
				  {
					  //if(snowdgtvariteflag == 1)
					  {
						  cout<<"Bisection canopy temperature solution failed with large range"<<endl;
						  cout<<"datetime, model element:"<< mtime[0]<<mtime[1]<<mtime[2]<<mtime[3]<<"  "<<uebCellY<<" "<<uebCellX<<endl;   //mtime[4]<<endl;				
						  cout<<" A canopy temperature of 273 K assumed"<<endl;
					  }
					  Tck=Tk;
					  goto labl10;
				  } 
				  else
				  {
					   if(snowdgtvariteflag == 1)
					  {
						  cout<<"Bisection canopy temperature solution with large range"<<endl;
						  cout<<"datetime, model element:"<< mtime[0]<<mtime[1]<<mtime[2]<<mtime[3]<<"  "<<uebCellY<<" "<<uebCellX<<endl;   //mtime[4]<<endl;			
						  cout<<" This is not a critical problem unless it happens frequently"<<endl;
						  cout<<" and solution below appears incorrect"<<endl;
					  }
				  }	
			  }       // else        endif
	//     Here do the bisection
		   niter = log((Tub-Tlb)/tol)/log(2.0);   // Number of iterations needed for temperature tolerance
		   for (int iter =0; iter< niter; ++iter)
		   {
			   Tck = 0.5*(Tub + Tlb);
			   F1   =  SURFEBC (Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,iTsMethod,mtime,
				                   	Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);  
                                              		//yjs add three value to reflect model control changes)) 
			  if(F1 > 0)   // This relies on the assumption (fact) that this is a monotonically decreasing function
				Tlb=Tck;
			  else
				Tub=Tck;
		   }
		   if(ibtowrite == 1)
			  if(snowdgtvariteflag == 1)
			  {
				  cout<<"Canopy temperature: "<<Tck<<endl;
				  cout<<"Energy closure: "<<F1<<endl;
				  cout<<"Iterations:"<<niter<<endl;
			  }		   

		 }

labl10:
     Tc = Tck - Tk;
	if((Wc > 0) && (Tc > 0))
	{             //dgt 5/4/04 surface melt smelt
          CanTemp_v = 0.0;
		  smeltC = SURFEBC (CanTemp_v+Tk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf, Qli,Qsi,atff,cosZen,EmC,Ems,param,sitev,iradfl,Qnetob,
			                   iTsMethod,mtime,Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz, Tssk,Tk,Tak,EA,RHOA, fkappaS,RHO,TherC, TSURFs,tave,refDepth);
		  CanTemp_v =  Tck*(1-Fs) + (CanTemp_v + Tk)*Fs - Tk;      // Equiv canopy temp during melt
																	//dgt 5/4/04 surface melt smelt is the energy not accommodated by conduction into the snow
																	//  so it results in surface melt which then infiltrates and adds energy to the snowpack
																	//  through refreezing
																	//  This change added to capture this quantity when the model is used for a glacier where
																	//  surface melt may run off rather than infiltrate.  
																	//  For modeling of glaciers in Antarctica by Mike Gooseff
	}
	else
	{
		CanTemp_v = Tc;
	    smeltC  = 0.0;
//  No surface melt this case
     }
      return CanTemp_v;
}

//      FUNCTION TO EVALUATE THE SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR SURFACE TEMPERATURE                      DGT and C Luce 4/23/97
float SURFEBSC(float Tssk, float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Fs,float Cf,float Qli,float Qsi,
                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpc,float Qps,float &Inmax,
				    float &Rkinc,float &Rkinsc,float &Vz,float Tck,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &TSURFs,float &Tave,float &refDepth
				)// Heat and vapor conductance for neutral
{
	float LanE_Ze, LanE_de, LanE_Ze2,LanE_de2;
	float Zs, Tsavek, Qp, qcs,fkappaG, d1,ze,de,dlf,ze2, de2;
	float Ess, Esc,d, Z0c,Tac, QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E, SURFEBSC_v;
	float  Taufb,Taufd,Qsib,Qsid, Betab,Betad,Taub,Taud, Qsns,Qsnc, Qlis,Qlns,Qlnc;

//     Parameters
     float cg   = param[3],				    //  Ground heat capacity [nominally 2.09 KJ/kg/C],
	   z   = param[4],					//  Nominal meas. height for air temp. and humidity [2m],
      zo   = param[5],					//  Surface aerodynamic roughness [m],
      rho  = param[6],				    //  Snow Density [Nominally 450 kg/m^3],
      rhog = param[7],					//  Soil Density [nominally 1700 kg/m^3],
      LanS = param[14],					// the thermal conductivity of fresh [dry], snow [0.0576 kJ/m/k/hr],
	  LanG = param[15],					// the thermal conductivity of soil [9.68 kJ/m/k/hr],
	  wlf  = param[16],					// Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr], 
	  rd1  = param[17],					// Apmlitude correction coefficient of heat conduction [1],
  
      APr   =   sitev[1],     //  Atmospheric Pressure [Pa],
	  Cc    =   sitev[4],     //  Canopy Coverage
      Hcan  =   sitev[5],     //  
      LAI   =   sitev[6];     //  

	  Zs     =  Ws*Rho_w/rho;                   // Snow Depth  
	  Tsavek = Tave + Tk;
      Qp     =  Qps;

	 if(snowdgtvariteflag2 == 1)
	  {
		  cout<<"Inputs to surfebsc"<<endl;
		  cout<<setprecision(15)<<Tssk<<" "<<Us<<" "<< Ws<<" "<< Wc<<" "<< A<<" "<< dt<<" "<< P<<" "<< Pr<<" "<< Ps<<" "<< Ta<<" "<< V<<" "<< RH<<" "<< Fs<<" "<< Cf<<" "<< Qli<<" "<< Qsi<<" ";
          cout<<setprecision(15)<<atff<<" "<< cosZen<<" "<< EmC<<" "<< Ems<<" "<< iradfl<<" "<< Qnetob<<" "<<iTsMethod<<endl;
	  }
//     07/25/02   at Boise.  To make the UEB2 work under Unix/linux system the fancy stuff like "Select case" shall not be used
	   //select case(iTsMethod)				//Four choice of the surface temperature modeling
	                                        //1. Old snow, use calibrated snow surface heat conductance         	   
				                          //2. Revised scheme LanE/Ze of the snow surface
										//3. Force restore approach	
										//4. Modified force restore approach.
	  //the follwing line is common to cases 2, 3, 4 #_6.7.13
	fkappaS = LanS/(rho*C_s);
	fkappaG = LanG/(rhog*cg);
	d1 = sqrt(2*fkappaS/W1da_y);
	if(Zs >= rd1*d1)
	{
		LanE_Ze=LanS/(rd1*d1);
	    ze=rd1*d1;
	}
	else
	{			
		LanE_Ze= LanE(LanS, LanG, Zs, rho, rhog, C_s, cg, rd1, ze, W1da_y);						//Cyjs   call the subroutine to update the heat conductivity. LanE()
		LanE_Ze=LanE_Ze/ze;
	}

	switch (iTsMethod)
	{
	case 1:   	                              //if(iTsMethod .eq. 1) then									
		qcs = rho*C_s*TherC*(Tssk-Tsavek);
		break;
	case 2:                                    //elseif (iTsMethod .eq. 2) then			
      	qcs= LanE_Ze*(Tssk-Tsavek);
		break;
	case 3:                                  //		elseif (iTsMethod .eq. 3) then			
		de=ze/rd1;
		LanE_de=LanE_Ze/de * ze;
	    qcs= LanE_de*(Tssk-Tssk_old)/(W1da_y*dt)+LanE_Ze*(Tssk-Tsavek);
		break;
	default:                                      //case (4)						//change to all default cases. ifnot for model comparison
			dlf = sqrt(2*fkappaG/wlf);
			if(Zs >= rd1*dlf)
			{
				LanE_Ze2=LanS/(rd1*dlf);
				ze2=rd1*dlf;
			}
			else
			{
				LanE_Ze2=LanE(LanS,LanG,Zs,rho,rhog,C_s,cg,rd1,ze2,wlf);  //Cyjs   call the subroutine to update the heat conductivity.                   
			    LanE_Ze2=LanE_Ze2/ze2;  
			}
			de = ze/rd1;
			LanE_de = LanE_Ze/de * ze;
			de2 = ze2/rd1;
			LanE_de2 = LanE_Ze2/de2 * ze2;		
			
			if((Us <= 0.0) || (refDepth <= 0.0))
			{
				qcs = LanE_de*(Tssk-Tssk_old)/(W1da_y*dt) +LanE_Ze*(Tssk-Tssk_ave) + LanE_de2*(Tssk_ave-Tsavek_ave);
				
			}
	        else if(refDepth > rd1*d1) 
			{
			     qcs = LanE_de*(Tssk-Tssk_old)/(W1da_y*dt)+LanE_Ze*(Tssk-Tssk_ave)+ LanE_de2*(Tssk_ave-Tsavek_ave);
				 
			}
            else
                qcs=LanE_Ze*ze*(Tssk-Tk)/refDepth;     
		  break;			
	  }		//End select      endif
	
	Ess    = svp(Tssk-Tk);
	Esc    = svp(Tck-Tk);

	if(snowdgtvariteflag2 == 1)
		cout<<setprecision(15)<<"qcs,Ess,Esc:"<<qcs<<" "<<Ess<<" "<<Esc<<endl;

	TURBFLUX(Ws,Wc,A,Tk,Tck-Tk,Ta,Tssk-Tk,RH,V,EA,P,param,sitev, d,Z0c,Vz,Rkinc,Rkinsc,Tac, Fs,Ess,Esc, QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E);        
	SURFEBSC_v = Qp + QHs + QEs- qcs;

	if(iradfl == 0) 
	{
	  PSOLRAD(Qsi,atff,param,Cf, Taufb,Taufd,Qsib,Qsid);    
	  TRANSRADCAN (cosZen,sitev,param, Betab,Betad,Taub,Taud);   
	  NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, Qsib,Qsid,param,Fs, Qsns,Qsnc ); //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
	  NETLONGRAD(RH,Ta,Tssk-Tk,Tck-Tk,Tk,Fs,EmC,Ems,SB_c,Cf,sitev, Qli,param, Qlis,Qlns,Qlnc);                     //  Output: Qsns,Qsnc
	  SURFEBSC_v = SURFEBSC_v + Qsns + Qlns;
	}
	else 
 		SURFEBSC_v = SURFEBSC_v + Qnetob;

	  if(snowdgtvariteflag2 == 1)
	  {
		 cout<<"Outputs from surfebsc"<<endl;
		 cout<<setprecision(15)<<Qpc<<" "<< Qps<<" "<< Inmax<<" "<<Rkinc<<" "<< Rkinsc<<" "<< Vz<<" "<< Tck<<" "<< Tk<<" ";
	     cout<<setprecision(15)<<Tak<<" "<<EA<<" "<<RHOA<<" "<< fkappaS<<" "<< RHO<<" "<< TherC<<" "<< TSURFs<<" ";
		 std::cout<<std::setprecision(15)<<Tave<<" "<<refDepth<<" "<<Qp<<" "<<QHs<<" "<<QEs<<" "<< qcs<<" "<<Qnetob<<" "<<Qsns<<" "<<Qlns<<" "<<Qsnc<<" "<<Qlnc<<" "<<Qlis<<" "<<SURFEBSC_v<<endl;
	  }
	return SURFEBSC_v;
}
	
//*************  FUNCTION TO EVALUATE THE CANPPY SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR CANOPY TEMPERATURE 
float SURFEBC(float Tck,float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Fs,float Cf,float Qli,float Qsi,
                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float QPc,float &QPs,float &Inmax,
				     float &Rkinc,float &Rkinsc,float &Vz,float Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &FkappaS,float RHO, float &TherC,float &TSURFs,float &tave,float &refDepth
			  )                             // Reduced parametre later                                    // Heat and vapor conductance for neutral
{
	float  SURFEBC_v;
	float  Taufb,Taufd,Qsib,Qsid;
	float Betab,Betad,Taub,Taud, Qsns,Qsnc,Qlis,Qlns,Qlnc,Ess,Esc;
	float d, Z0c,Tac, QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E;
//     Parameters
     float cg   = param[3],				    //  Ground heat capacity [nominally 2.09 KJ/kg/C],
	z    = param[4],					//  Nominal meas. height for air temp. and humidity [2m],
      zo   = param[5],					//  Surface aerodynamic roughness [m],
      rho  = param[6],				    //  Snow Density [Nominally 450 kg/m^3],
      rhog = param[7],					//  Soil Density [nominally 1700 kg/m^3],
      lans = param[14],					// the thermal conductivity of fresh [dry], snow [0.0576 kJ/m/k/hr],
	lang = param[15],					// the thermal conductivity of soil [9.68 kJ/m/k/hr],
	wlf  = param[16],					// Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr], 
	rd1  = param[17],					// Apmlitude correction coefficient of heat conduction [1],
  
      APr  = sitev[1];                 //  Atmospheric Pressure [Pa],

       if(iradfl !=1)
	   {
			PSOLRAD(Qsi,atff,param, Cf, Taufb,Taufd,Qsib,Qsid);       
			TRANSRADCAN (cosZen,sitev,param, Betab,Betad,Taub,Taud);                     // Output variables:
			NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, Qsib,Qsid,param,Fs,Qsns,Qsnc);                       //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
            NETLONGRAD(RH,Ta,Tssk-Tk,Tck-Tk,Tk,Fs,EmC,Ems,SB_c,Cf,sitev, Qli,param, Qlis,Qlns,Qlnc);                     //  Output: Qsns,Qsnc
	   }
   	Ess    = svp(Tssk-Tk);
	Esc    = svp(Tck-Tk);  
    TURBFLUX(Ws,Wc,A,Tk,Tck-Tk,Ta,Tssk-Tk,RH,V,EA,P,param,sitev, d,Z0c,Vz,Rkinc,Rkinsc,Tac, Fs,Ess,Esc, QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E);  
//---> //the following line added 6.7.13
	  SURFEBC_v = QPc+QHc+QEc;                                    //- Qcs

    if(iradfl == 0)
		SURFEBC_v = Qsnc + Qlnc+QPc+QHc+QEc;              
	else
		SURFEBC_v += Qnetob;                              //surfeb = surfeb + Qnetob
    
	return SURFEBC_v;
}

//********** FUNCTION TO ESTIMATE SURFACE HEAT CONDUCTION FOR USE IN SOLVING SURF TEMP  QcEst() *******
float QcEst(float Ws, float P, float Tssk,float &Tck, float V, float &Zm, float &d, float &Z0c, float &Rimax, float &Rcastar, float Cf, float Fs, float Qli, 
               float &Hcan, float Vz,float Ta, float Rh, float RKINsc, float QPs, float To, float Ps,float Qsi, float atff, float cosZen, float APr,float Tak,
			     float &EA, float &A, float &Ac, float Wc, float Inmax, float Qnetob, int Iradfl, float *param, float *sitev )
{
	//*********** Hcan never used here
	//TBC 6.11.13
	float  Tss, Ess, Esc, QcEst_v;
	float QHc, QEc, Ec, QHs, QEs, Es, QH, QE, E;
	float Taufb, Taufd, Qsib, Qsid,
		    Betab, Betad, Taub, Taud,
			  Qsns, Qsnc,
			    Qlis,Qlns,Qlnc;
	float Rkinc, Tac;

	float Ems   =   param[2],   //  emmissivity of snow [nominally 0.99],
	cg    =   param[3],   //  Ground heat capacity [nominally 2.09 KJ/kg/C],
    z     =   param[4],   //  Nominal meas. height for air temp. and humidity [2m],
    zo    =   param[5],   //  Surface aerodynamic roughness [m],
    fstab =   param[18],  //  Stability correction control parameter 0 = no corrections, 1 = full corrections
	EmC   =   param[22],	 //  Emissivity of canopy 
  	Cc    =   sitev[4],   //  Canopy Coverage
	LAI   =   sitev[6];   //  Leaf Area Index  

	if(Wc == 0)                          
	  Tck = Tak;                                     // Tc = Ta assumed here when there is no snow in canopy
	else
		Tck= T_k;                                    // When snow at gound is melting, canopy snow is melting (assumed)
	
	Tss = Tssk- T_k;
	Ess  = svp (Tssk - T_k);
	Esc  = svp(Tck - T_k); 

	TURBFLUX (Ws, Wc, A, T_k, Tck-T_k, Ta, Tssk-T_k, Rh, V, EA, P, param, sitev, d, Z0c, Vz, Rkinc, RKINsc, Tac, Fs, Ess, Esc, 		 
		         QHc, QEc, Ec, QHs, QEs, Es, QH, QE, E);          // Output variables

    QcEst_v = QPs + QHs + QEs; 
                                                      
	if(Iradfl == 0)
	{
		PSOLRAD (Qsi, atff, param, Cf, Taufb, Taufd, Qsib, Qsid);    
		TRANSRADCAN (cosZen,sitev,param, Betab,Betad,Taub,Taud);   // Output variables:
	    NETSOLRAD (Ta, A, Betab, Betad, Wc,Taub, Taud, Qsib, Qsid, param, Fs, Qsns, Qsnc); //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
	    NETLONGRAD (Rh,Ta,Tssk-T_k,Tck-T_k, T_k, Fs, EmC, Ems, SB_c, Cf, sitev,	Qli, param, Qlis,Qlns,Qlnc);                     //  Output: Qlis,Qlns,Qlnc

	    QcEst_v = QcEst_v + Qsns + Qlns;
	}
	else
 	  QcEst_v = QcEst_v + Qnetob;   
	  return QcEst_v;
}

//*********************** FMELT () ****************************************************************
//     Calculates the melt rate and melt outflow
float FMELT(float UB,float Rho_w,float W,float HF,float LC,float RID,float KS,float PRAIN)
{
      float  FMELT_v, UU, SS;
//       write(*,*) 'I am in FMELT////////'
      UU = 0.0;
      if(UB < 0)
         FMELT_v = 0.0;
	  else if (W <= 0)
	  {
         FMELT_v = PRAIN;
         if(PRAIN <= 0) 
            FMELT_v = 0.0;
	  }
	  else
	  {
         UU = UB/(Rho_w*W*HF);     //  liquid fraction
         if(UU > 0.99)
            UU = 0.99;      //     TO ENSURE HIGH MELT RATE WHEN ALL LIQUID   
		 if((UU/(1-UU)) <= LC)
			 SS=0.0;
		 else 
			 SS = (UU/((1- UU)*RID) - LC/RID) / (1-LC/RID);      
         FMELT_v = KS*pow(SS,3);
	  }

      if(FMELT_v < 0)
	  {
          cout<<"FMELT is NEGATIVE"<<endl;
	      getchar();
	  }    
      return FMELT_v;
}

//******************************************** Grad () ********************************************
// Linear simple function to calculate the gradient of Q and T
void Grad(float qc1, float qc2,float t1,float t2, float &a,float &b)
{
	if((t2-t1) != 0) 
	{
	  b = (qc2-qc1)/(t1-t2);
	  a = qc1 + b*t1;
	}
	return;
}

//*************************** refDep() ************************************************************
//     function to calculate the refreezing depth
float  refDep(float flans,float a,float b,float hf,float  rhom,float dt,float x1 )
{
	float temp, refDep_v;

	temp = flans*flans + 2*b*(-a*flans*dt/(hf*rhom) + flans*x1 + 0.5*b*x1*x1);
	if((temp < 0) || (a > 0) || (b == 0)) 
	  refDep_v = 0;
	else 
	  refDep_v = (-flans + sqrt(temp))/b; 	   	
	
	return refDep_v;
}

//******************** CALCULATE THE AVE.TEMP.OF SNOW AND INTERACTING SOIL LAYER ********************
float TAVG( float UB, float W, float RHOW, float CS, float TO, float RHOG, float DE, float CG, float HF)
{
	float SNHC, SHC, CHC, TS, AL, TAVG_v;
    SNHC = RHOW * W * CS;							 // SNHC = Snow heat capacity
    SHC  = RHOG * DE * CG;							 // SHC  = Soil heat capacity //DE   = Therm. active depth of soil (param(11)) 
    CHC  = SNHC + SHC;								 // CHC  = Combined heat capacity

    if(UB <= 0) 
        TS = UB/CHC;
	else
	{
        AL = UB /(RHOW*HF);
        if(W > AL)
          TS=TO;                                   
		else
          TS = (UB - W*RHOW*HF)/CHC;       
	}
    TAVG_v = TS;
	return TAVG_v;    
}

//********************* Caculate the daily average  *********************
//yjs  Calculate the daily average value //yjs  n number of records, a minimum value -100 or some
float daily_ave(float *backup, int n, float a)
{
	float daily_avev = 0, sum = 0.0;
	int count = 0;
	for (int i=0;i<n;i++)
		if(backup[i] > a)
		{
		  sum += backup[i];
		  count++; 
		}
	if(count != 0) 
		daily_avev = sum/count;
	return daily_avev;
}

//***************************** LanE () *************************************
//	Function to get the LanE which is the thermal conductivity by ze
float  LanE( float LanS, float LanG, float Zs,float rho, float rhog,float cs, float cg, float r,float &ze,float w1day)
{
	float fKappaS, fKappaG, d1, d2, LanE_v;	

	 fKappaS = LanS/(rho*cs);
	 fKappaG = LanG/(rhog*cg);
	 d1 = sqrt(2*fKappaS/w1day);
	 d2 = sqrt(2*fKappaG/w1day);
	 LanE_v = findMin(r,Zs/d1)*d1/LanS + findMax(0.0,(r-Zs/d1))*d2/LanG;
	 ze=findMin(r,Zs/d1)*d1 + findMax(0.0,(r-Zs/d1))*d2;	 
	 if(LanE_v != 0) LanE_v = 1/LanE_v * ze;
	 return LanE_v;
}

//******************************** ALBEDO () *****************************************************
//     Function to calculate Albedo.   BATS Albedo Model (Dickinson et. al P.21)
float ALBEDO(float tausn, float cosZen, float d, float aep,float abg, float avo,float airo)     // snow age, zenith angle, snow depth
{                                                           // albedo extinpara,bare Gr albedo, visoble and NIR for fresh snow
      float B  = 2.0;                                             // aep.input file
      float CS = 0.2, CN = 0.5;
	  float FZEN, AVIS, AIRD, ANIR, ALBEDO_v, rr, AVD;
      float FAGE = tausn/(1.0 + tausn);

      if(cosZen < 0.5)
		  //#_Check operator precedence issuse here
		//#$%_6.5.13
         FZEN = 1.0 / B*((B+1.0)/(1.0 + 2.0 * B * cosZen) - 1.0); 
	  else
         FZEN = 0.0; 
      AVD = (1.0 - CS*FAGE )*avo;
      AVIS = AVD + 0.4*FZEN *(1.0 - AVD);
      AIRD = (1.0 - CN*FAGE) * airo;
      ANIR = AIRD + 0.4*FZEN*(1.0 - AIRD);
      ALBEDO_v = (AVIS + ANIR)/2.0;
      if(d < aep)                              // need to transition albedo to a bare ground value
	  {
        rr = (1.0 - d/aep) * exp(-d * 0.5 / aep);
		ALBEDO_v = rr*abg + (1.0 - rr)* ALBEDO_v;
	  }
	  return ALBEDO_v;
}

//******************************** AGESN () *******************************************************
//     Function to calculate Dimensionless age of snow for use in BATS Albedo Model (Dickinson et. al P.21)
void AGESN(float &tausn, float dt, float Ps, float tsurf, float tk, float dNewS) 
{	     
//     Input Variables  : dt, ps,tsurf,tk,dnewS (threshold depth of new snow, param 21)
//     Output Variable  : Agesn
      float tsk, R1, R2, R3;
	  tsk = tsurf + tk;                                         // Express surface temperature in kelvin
      R1 = exp(5000.0*(1.0/tk - 1.0/tsk));
      R2 = pow(R1,10);
      if(R2 > 1)
         R2 = 1.0;
      R3 = 0.03;
      tausn = tausn + 0.0036*(R1+R2+R3)*dt;
	  if(Ps > 0) 
		  if(dNewS > 0) 
			tausn = findMax(tausn*(1-Ps*dt/dNewS),0.0);
		  else
			tausn = 0.0;
      return;
}

//************************** PARTSNOW () **********************************************************
//     Partitioning of precipitation into rain and snow      
float PARTSNOW(float P, float TA, float TR, float TS)
{
	float PARTSNOW_v;
    if(TA < TS)
        PARTSNOW_v = P;
    else if(TA > TR)
        PARTSNOW_v = 0.0;
	else
        PARTSNOW_v = P*(TR-TA)/(TR-TS);
	return PARTSNOW_v;
}
