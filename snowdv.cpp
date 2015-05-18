//copied (and modified) from snowdv.f90  
#include "uebpgdecls.h"

//**********************************************************************************************
void RUNUEB(float* RegArray[13], float statesiteValues[], float paramValues[], float** &OutVarValues, int modelStartDate[],
			     double modelStartHour, int modelEndDate[], double modelEndHour, double modelDT, double UTCOffset)
{   
      const int nsv = 10, npar = 32,nxv = 6,niv = 8, nov = 14;
	  float  Inpt[niv], sitev[nsv], outv[nov], statev[nxv],  Param[npar], dtbar[12], mtime[4];                                     // YJS pass to reflect the change of Snow (Year, month, Date, Hour)
	  int irad, ireadalb, subtype, iradfl, iflag[6]; 
	  float  slope,azi, bca,bcc, lat, ts_last, lon;
	  float *tsprevday, *taveprevday;
	  int stepinaDay, nstepinaDay; 
	  //float referenceHour, referenceTime, CTJD;
	  double dHour, EJD, MHour, sHour, currentModelDateTime, UTCHour, OHour, UTCJulDat;
	  float fHour, fModeldt;
	  int istep,  Year, Month, Day, MYear, MMonth, MDay; 	  
	  
	    // CHANGES TO ACCOMODATE GLACIER
      float  WGT; // WGT=WATER EQUIVALENT GLACIER THICKNESS
	  float Us,Ws, Wc, Apr, cg, rhog, de, tave, Ws1, Wc1, cumP,cumEs,cumEc,cumMr, cumGm, cumEg, Ta, P, V, RH, Tmin = 0.0, Tmax = 0.0,Trange, Qsiobs, Qg, Qli, QLif;
	  float Vp; //#12.18.14 Air vapor pressure 
      float Qnetob,Snowalb,cosZen,atff,cf,atfimplied, Ema,Eacl,dStorage,errMB, HRI0, HRI, as, bs;	   
	  float  OutArr[53];  
      //to get max/min daily temperatures
	  int nb, nf, nbstart, nfend;
	  double dstHour; 

     //  Mapping from parameters read to UEB internal interpretation which follows UEBVeg scheme
        irad = (int) paramValues[0];	
        ireadalb=(int)paramValues[1];
		for(int i=0;i<11;i++)
			Param[i] = paramValues[i+2];
		for(int i=12;i<18;i++)
			Param[i] = paramValues[i+1];      
        Param[18]=-9999;
        Param[19]=-9999;
        Param[20]=paramValues[19];
		for(int i=22;i<32;i++)
			Param[i] = paramValues[i-2];      
        bca = paramValues[30];
        bcc = paramValues[31];  
        
    //copied from paramsiteinitial
		for (int i = 0; i<4;i++)
			statev[i]=statesiteValues[i];
		 sitev[0] = statesiteValues[4];
		 sitev[1]   = statesiteValues[5];
		 for (int i = 3; i<9;i++)
			sitev[i]=statesiteValues[i+3];
		slope=statesiteValues[12];
        azi=statesiteValues[13];
        lat=statesiteValues[14];          
        
        Param[11]=statesiteValues[15];
        //subalb=statesiteValues[15]
        sitev[9]=statesiteValues[16];
        subtype= (int)statesiteValues[16];
        Param[21]=statesiteValues[17];
        //gsurf = statesiteValues[17]
		for(int i=0;i<12;i++)
			dtbar[i] = statesiteValues[i+18];
        ts_last = statesiteValues[30];
        lon = statesiteValues[31];          
       //  FIXME: what if the result is fractional
//  time steps must divide exactly in to a day because we use logic that requires the values from the same time
//  step on the previous day.  Consider in future making the specification of time step as number of time
//  steps in a day, not modeldt to ensure this modeldt is recalculated based on the int timesteps in a day
//  assumption: number of model timesteps in a day must be an int             
      stepinaDay= (int) (24.0/modelDT +0.5);  // closest rounding
      modelDT = 24.0/stepinaDay;
	  nstepinaDay = stepinaDay; 	  
      
	  tsprevday = new float[nstepinaDay];
	  taveprevday = new float[nstepinaDay];   
	 
    // calculating model end date-time in julian date
      dHour = modelEndHour;
      EJD = julian(modelEndDate[0], modelEndDate[1], modelEndDate[2],dHour);     
     
        if (sitev[9] ==0 || sitev[9] ==3)
			   WGT = 0.0;
		else
			WGT = 1.0;      
       
        if(sitev[9] != 3)  //  Only do this work for non accumulation cells where model is run
		{
            //  Initialize Tsbackup and TaveBackup
			for(int i =0;i< nstepinaDay;i++)
			{
				tsprevday[i] = -9999.0;
				taveprevday[i] = -9999.0;
			}
			// Take surface temperature as 0 where it is unknown the previous time step
			// This is for first day of the model to get the force restore going
	//#$#$#$#$#_is this all the time steps or the last time?
			if(ts_last <= -9999)
				//for(int i =0;i< nstepinaDay;i++)
					tsprevday[nstepinaDay-1] = 0;
			else 
				//for(int i =0;i< nstepinaDay;i++)
				tsprevday[nstepinaDay-1] = ts_last;
     // compute Ave.Temp for previous day 
            Us   = statev[0];                                  // Ub in UEB
            Ws   = statev[1];                                  // W in UEB
            Wc   = statev[3];                                  // Canopy SWE
            Apr  = sitev[1];                                   // Atm. Pressure  [PR in UEB]
            cg   = Param[3];                                   // Ground heat capacity [nominally 2.09 KJ/kg/C]
            rhog = Param[7];                                   // Soil Density [nominally 1700 kg/m^3]
            de   = Param[10];                                  // Thermally active depth of soil (0.1 m)

			//this are for coudiness computation
			//6.10.13
			as = Param[27];
			bs = Param[28];
 
			tave = TAVG(Us,Ws+WGT,Rho_w,C_s,T_0, rhog, de, cg, H_f);    
            //for(int i =0;i< nstepinaDay;i++)
			taveprevday[nstepinaDay-1] = tave;      
        //  initialize variables for mass balance
            Ws1   = statev[1];
            Wc1   = statev[3];    
            cumP = 0.0;
            cumEs = 0.0;
            cumEc = 0.0;
            cumMr = 0.0;
            cumGm = 0.0;
			cumEg = 0.0;

		} //  end the skip block done only for accumulation cells

          // Variables to keep track of which time step we are in and which netcdf output file we are in
          istep = 0;  // time step initiated as 0
        
          // map on to old UEB names
          Year = modelStartDate[0];
          Month = modelStartDate[1];
		  Day = modelStartDate[2];          
          sHour = modelStartHour;
          currentModelDateTime = julian(Year, Month, Day, sHour);
		  
          //++++++++++++++++++++This is the start of the main time loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                             
		  while(EJD >= currentModelDateTime)
		  {
 			  //istep++;             
			
            if(sitev[9] != 3)
			{            
			   // UTC to local time conversion
			   calendardate(currentModelDateTime,Year, Month, Day, dHour);
			   UTCHour = dHour - UTCOffset;
			   OHour = UTCHour + lon/15.0;
			   UTCJulDat = julian(Year, Month, Day ,OHour);
			   calendardate(UTCJulDat, MYear, MMonth, MDay, MHour);
			   fHour = (float) MHour;
			   fModeldt = (float) modelDT;

		  //      Map from wrapper input variables to UEB variables     
			   if ( inpDailyorSubdaily == 0)       //inputs are given for each time step (sub-daily time step)--in m/hr units 
			   {	
				P = RegArray[1][istep];
				V = RegArray[2][istep];
				   Ta =RegArray[0][istep];

				//get min max temp
				Tmax = RegArray[0][istep];
				Tmin = RegArray[0][istep];
				Trange = 8;
//get max/min temperature during the day 
				 nb = (dHour - 0)/modelDT;		 //number of time steps before current time within same day
				 nf = (24 - dHour)/modelDT;      //no of time steps after current time within the same day
//#_TBC 9.13.13 look for better method for the following
				 if(dHour > 23)                  //to take care of hour 24, <=>0hr
				 {
					 nb = 0;
					 nf = 24/modelDT;
				 }
				 nbstart = findMax(istep-nb,0);  //to gard against going out of lower bound near start time when the start time isnot 0 hr (istep < nb )
				 nfend = istep + nf;
				
				for (int it = nbstart; it < nfend; ++it)
				{
					if (RegArray[0][it] <= Tmin)					
						Tmin = RegArray[0][it];
				    if (RegArray[0][it] >= Tmax)					
						Tmax = RegArray[0][it];
				}
				Trange = Tmax - Tmin;
				//cout<<Trange<<endl;
				if (Trange <= 0)
				{
					if (snowdgtvariteflag==1)	
					{
						cout<<"Input Diurnal temperature range is less than or equal to 0 which is unrealistic "<<endl;
						cout<< "Diurnal temperature range is assumed as 8 degree celsius on "<<endl;
					    cout<< Year<<" "<< Month<<" "<<Day<<endl;
					}
					Trange = 8.0;
				}

		//		 Flag to control radiation (irad)
  //!     0 is no measurements - radiation estimated from diurnal temperature range
  //!     1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
  //!     2 is incoming shortwave and longwave radiation read from file (measured)
  //!     3 is net radiation read from file (measured)
 				switch(irad)
				{
					case 0:
						Qsiobs = RegArray[5][1];
						Qli = RegArray[6][1];
						Qnetob = RegArray[7][1];						
						break;
					case 1:
					   Qsiobs = RegArray[5][istep];                         // *3.6; // Daymet srad in W/m^2
						Qli = RegArray[6][1];
						Qnetob = RegArray[7][1];
						break;
					case 2:
					   Qsiobs = RegArray[5][istep];                         // *3.6; // Daymet srad in W/m^2
						Qli = RegArray[6][istep];
						Qnetob = RegArray[7][1];
						break;
					case 3:
						Qsiobs = RegArray[5][1];
						Qli = RegArray[6][1];
						Qnetob = RegArray[7][istep];
						break;
					default:
						cout<<" The radiation flag is not the right number; must be between 0 and 3"<<endl;
						getchar();
						break;				
				}
				//atm. pressure from netcdf 10.30.13
				//this needs revision 		//####TBC_6.20.13
				if(RegArray[4][0] == 2)
					sitev[1] = RegArray[4][1];
				else
					sitev[1] = RegArray[4][istep];

				//this needs revision 		//####TBC_6.20.13
				if(RegArray[8][0] == 2)
					Qg = RegArray[8][1];
				else
					Qg = RegArray[8][istep];
				 //!     Flag to control albedo (ireadalb)  
				//!     0 is no measurements - albedo estimated internally
                //!     1 is albedo read from file (provided: measured or obtained from another model)
				//these need revision //####TBC_6.20.13
				if(RegArray[9][0] == 2)
					Snowalb = RegArray[9][1];
				else
					Snowalb = RegArray[9][istep];					
				   //12.18.14 Vapor pressure of air
				   if(RegArray[12][0] == 2)
					   Vp = RegArray[12][1];
				   else
					   Vp = RegArray[12][istep];
				   //relative humidity computed or read from file
				   //#12.18.14 needs revision
				   if (RegArray[3][0] == 2)
				   {
					   RH = RegArray[3][1];					   
				   }
				   else if (RegArray[3][0] == -1)          //RH computed internally 
				   {
					   float eSat = 611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
					   RH = Vp / eSat;
				   }
				   else
					   RH = RegArray[3][istep];
				   if (RH > 1)
				   {
					   //cout<<"relative humidity >= 1 at time step "<<istep<<endl;
					   RH = 0.99;
				   }
			   }   
			   else		//inputs are given as average daily values, precip unit is always m/hr, Tmax and Tmin are manadatory       
			   {
				   //average daily value of precip in units of m/hr 4.22.14
				   P = RegArray[1][istep / nstepinaDay];            // 24000; //#12.19.14 Daymet prcp in mm/day           
				   V = RegArray[2][istep/nstepinaDay];
				    
				   //get min max temp
				   Tmin = RegArray[10][istep/nstepinaDay];				  
				   Tmax = RegArray[11][istep/nstepinaDay];
				   //cout << "Tmin = " << Tmin << "Tmax = " << Tmax << " ";

				   Trange = Tmax - Tmin;
				   //cout<<Trange<<endl;
				   if (Trange <= 0)
				   {
					   if (snowdgtvariteflag==1)	
					   {
						   cout<<"Input Diurnal temperature range is less than or equal to 0 which is unrealistic "<<endl;
						   cout<< "Diurnal temperature range is assumed as 8 degree celsius on "<<endl;
						   cout<< Year<<" "<< Month<<" "<<Day<<endl;
					   }
					   Trange = 8.0;
				   }
				   //sin curve describes the diel temperature fluctuations with max at 15 hrs and min at 3 hrs
				   Ta = Tmin + 0.5*Trange + 0.5*Trange*sin(2*P_i*(fHour + 15.0)/24);
				   //! Flag to control radiation (irad)
				   //!     0 is no measurements - radiation estimated from diurnal temperature range
				   //!     1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
				   //!     2 is incoming shortwave and longwave radiation read from file (measured)
				   //!     3 is net radiation read from file (measured)
				   switch(irad)
				   {
				   case 0:
					   Qsiobs = RegArray[5][1];
					   Qli = RegArray[6][1];
					   Qnetob = RegArray[7][1];						
					   break;
				   case 1:
					   Qsiobs = RegArray[5][istep / nstepinaDay];                 // *3.6; // Daymet srad in W/m^2
					   Qli = RegArray[6][1];
					   Qnetob = RegArray[7][1];
					   break;
				   case 2:
					   Qsiobs = RegArray[5][istep / nstepinaDay];                 // *3.6; // Daymet srad in W/m^2
					   Qli = RegArray[6][istep/nstepinaDay];
					   Qnetob = RegArray[7][1];
					   break;
				   case 3:
					   Qsiobs = RegArray[5][1];
					   Qli = RegArray[6][1];
					   Qnetob = RegArray[7][istep/nstepinaDay];
					   break;
				   default:
					   cout<<" The radiation flag is not the right number; must be between 0 and 3"<<endl;
					   getchar();
					   break;				
				   }
				   //atm. pressure from netcdf 10.30.13
				   //this needs revision 		//####TBC_6.20.13
				   if(RegArray[4][0] == 2)
					   sitev[1] = RegArray[4][1];
				   else
					   sitev[1] = RegArray[4][istep/nstepinaDay];

				   //this needs revision 		//####TBC_6.20.13
				   if(RegArray[8][0] == 2)
					   Qg = RegArray[8][1];
				   else
					   Qg = RegArray[8][istep/nstepinaDay];
				   //!     Flag to control albedo (ireadalb)  
				   //!     0 is no measurements - albedo estimated internally
				   //!     1 is albedo read from file (provided: measured or obtained from another model)
				   //these need revision //####TBC_6.20.13
				   if(RegArray[9][0] == 2)
					   Snowalb = RegArray[9][1];
				   else
					   Snowalb = RegArray[9][istep/nstepinaDay];
				   //12.18.14 Vapor pressure of air
				   if (RegArray[12][0] == 2)
					   Vp = RegArray[12][1];
				   else
					   Vp = RegArray[12][istep/nstepinaDay];
				   //relative humidity computed or read from file
				   //#12.18.14 needs revision
				   if (RegArray[3][0] == 2)
				   {
					   RH = RegArray[3][1];
				   }
				   else if (RegArray[3][0] == -1)          //RH computed internally 
				   {
					   float eSat = 611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
					   RH = Vp / eSat;
				   }
				   else
					   RH = RegArray[3][istep/nstepinaDay];
				   if (RH > 1)
				   {
					   //cout<<"relative humidity >= 1 at time step "<<istep<<endl;
					   RH = 0.99;
				   }
			   }	
        
		//  Below is code from point UEB 
				sitev[2]= Qg;   
				Inpt[0] =Ta;
				 Inpt[1]=P;
				Inpt[2]=V;
				 Inpt[3]=RH;
				 Inpt[6] = Qnetob;         
			
			   //Radiation Input Parameterization  
			   hyri(MYear, MMonth, MDay, fHour, fModeldt,slope, azi, lat, HRI, cosZen); 
			   Inpt[7] = cosZen;                
			   if (irad <= 2)
			   {   
				   atf(atff,Trange, Month,dtbar,bca,bcc);      
		// We found that Model reanalysis and dowscaled data may produce some unreasonably negative solar radiation. this is simply bad data and it is generally better policy to try to give a model good data. 
		// If this is not possible, then the UEB checks will avoid the model doing anything too bad, it handles negative solar radiation in following way:
		// "no data in radiation would be to switch to the temperature method just for time steps when radiation is negative." 

				   if( irad ==0 || Qsiobs < 0)     //  For cases where input is strictly negative we calculate QSI from HRI and air temp range.  This covers the case of missing data being flagged with negative number, i.e. -9999.                 
				   {   
					   Inpt[4] = atff* Io *HRI;
					   cloud(as, bs, atff,cf);   // For cloudiness fraction
				   }
				   else   // Here incoming solar is input
				   {
			//      Need to call HYRI for horizontal surface to perform horizontal measurement adjustment
						hyri(MYear, MMonth, MDay, fHour, fModeldt, 0.0, azi, lat, HRI0, cosZen);
			//      If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
			//      not it indicates a potential measurement problem. i.e. moonshine
						if(HRI0 > 0) 
						{
							//cout<<Qsiobs;
							atfimplied =  findMin(Qsiobs/(HRI0*Io),0.9); // To avoid unreasonably large radiation when HRI0 is small
							Inpt[4] = atfimplied * HRI * Io;
						}
						else
						{
							Inpt[4] = Qsiobs; 
							if(Qsiobs != 0)
							{
								if (radwarnflag < 3)   //leave this warning only three times--enough to alert to non- -ve night time solar rad
							  {
									 cout<<"Warning: you have nonzero nightime incident radiation of "<<Qsiobs<<endl;
									 cout<<"at date "<<Year<<"   "<< Month<<"   "<< Day<<"     "<<dHour<<endl;
									 ++radwarnflag;
							  }
							}
						}
						cloud(as,bs,atff,cf);   // For cloudiness fraction  This is more theoretically correct
					}        
					if(irad < 2)
					{
						qlif(Ta, RH, T_k, SB_c, Ema,Eacl,cf,QLif);
						Inpt[5] = QLif;
					}
					else 
					{
						Ema = -9999;  //  These values are not evaluated but may need to be written out so are assigned for completeness
						Eacl = -9999;
						Inpt[5] = Qli;
					} 
					iradfl = 0;  
			   }   // Long wave or shortwave either measured and calculated
			   else
			   {
					  iradfl = 1;                    // This case is when given IRAD =3 (From Net Radiation)  
					  Inpt[6] = Qnetob;           
			   }        
	
		//      set control flags
				iflag[0] = iradfl;   // radiation [0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7]
				//  In the code above radiation inputs were either computed or read from input files
				iflag[1] = 0;        // no 0 [/yes 1] printing
				//iflag[2] = outFile;        // Output unit to which to print
				if(ireadalb == 0)
				   iflag[3] = 1;        // Albedo Calculation [a value 1 means albedo is calculated, otherwise statev[3] is albedo
				else
				{
				   iflag[3]=0;
				   statev[2]=Snowalb;
				}  
				
			/*	if (istep >=48)
				{
					snowdgtvariteflag = 1;
					snowdgtvariteflag2 = 1;
					snowdgtvariteflag2 = 1;
					getchar();
				}*/

				//added 9.16.13
				iflag[4] = 4;
				mtime[0] = Year;
				mtime[1] = Month;
				mtime[2] = Day;
				mtime[3] = dHour;
				
				SNOWUEB2(fModeldt, Inpt, sitev, statev, tsprevday, taveprevday, nstepinaDay, Param,iflag,
						   cumP, cumEs, cumEc, cumMr, cumGm, cumEg, outv, mtime, atff, cf, OutArr);     
		  
				dStorage = statev[1]-Ws1+ statev[3]-Wc1;
				errMB= cumP-cumMr-cumEs-cumEc -dStorage+cumGm - cumEg; 
				
				OutVarValues[0][istep]= Year;
				OutVarValues[1][istep]=Month;
				OutVarValues[2][istep]=Day;
				OutVarValues[3][istep]=dHour;
				OutVarValues[4][istep]=atff;
				OutVarValues[5][istep]=HRI;
				OutVarValues[6][istep]=Eacl;
				OutVarValues[7][istep]=Ema;
				OutVarValues[8][istep]=Inpt[7]; //cosZen
				OutVarValues[9][istep]=Inpt[0];
				OutVarValues[10][istep]=Inpt[1];
				OutVarValues[11][istep]=Inpt[2];
				OutVarValues[12][istep]=Inpt[3];
				OutVarValues[13][istep]=Inpt[4];
				OutVarValues[14][istep]=Inpt[5];
				OutVarValues[15][istep]=Inpt[6];	
							
				for (int i=16;i<69;i++)
				{		   
					OutVarValues[i][istep] = OutArr[i-16];					
				}
				OutVarValues[69][istep] = errMB;

				if (snowdgt_outflag == 1)             //if debug mode 
				{
					for(int uit = 0; uit<70; uit++)
						cout<<OutVarValues[uit][istep]<<" ";
					cout<<endl;

					cout<<endl;
				}

           }
		   else //sitev[9] ==3 // this block entered only if sitev(10)= 3
		   {
                errMB = 0.0;
				for (int i=0;i<53;i++)
					OutArr[i] = 0.0;
				for (int i= 0;i <70;i++)
					OutVarValues[i][istep] = 0.0;              
		   }          
		
			istep++;				
			UPDATEtime(Year, Month, Day, dHour,modelDT);  			
            //ModHour=DBLE(Hour)
            currentModelDateTime = julian(Year, Month, Day, dHour); 

	   }  //End of the main time loop 
	   
	   delete []tsprevday;
	   delete []taveprevday;  
    
    return; 
}
