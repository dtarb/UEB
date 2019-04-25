//copied (and modified) from snowdv.f90  
#include "gpuuebpgdecls.h"


//this uses class forcing arrays
__host__ __device__ void uebCell::runUEB(int dimLen2)
{
	//dimLen2
	runUEB();
}

__host__ __device__ int uebCell::getforcOffset(int ioffst, int dimLen2)
{
	if (ioffst == 1)
		return uebCellY*dimLen2*numTimeStep + uebCellX*numTimeStep;
	else
		return 0;
}
//this is for gpu; enables passing 'outside' forcing arrays
__host__ __device__ void uebCell::runUEB()
{ 
          // Variables to keep track of which time step we are in and which netcdf output file we are in
          istep = 0;  // time step initiated as 0        
          // map on to old UEB names
          Year = modelStartDate[0];
          Month = modelStartDate[1];
		  Day = modelStartDate[2];          
          sHour = modelStartHour;
          currentModelDateTime = julian(Year, Month, Day, sHour);
		  
		  int indx = uebCellY + uebCellX;     //blockIdx.x*blockDim.x + threadIdx.x;
          //++++++++++++++++++++This is the start of the main time loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                             
		  while (istep < numSimTimeSteps)						// && EJD >= currentModelDateTime)
		  {
			  //cout<<istep<<endl;             
	       //if(indx == 0)
		// printf(" time step: %d    ", istep);		
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
				   P = PrecArr[istep];                               // / 24000;   #12.19.14 --Daymet prcp in mm/day
				   V = WindspArr[istep];
				   Ta = TempArr[istep];

				   //get min max temp
				   Tmax = TempArr[istep];
				   Tmin = TempArr[istep];
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
				   nbstart = findMax(istep - nb, 0);  //to guard against going out of lower bound near start time when the start time is not 0 hr (istep < nb )
				   nfend = findMin(istep + nf, numSimTimeSteps - 1); //don't go out of upper limit
				
				   for (int it = nbstart; it < nfend; ++it)
				   {
					   if (TempArr[it] <= Tmin)					
						   Tmin = TempArr[it];
					   if (TempArr[it] >= Tmax)					
						   Tmax = TempArr[it];
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
					   Qsiobs = infrContArr[8].infdefValue;
					   Qli = infrContArr[9].infdefValue;
					   Qnetob = infrContArr[10].infdefValue;
					   break;
				   case 1:
					   Qsiobs = SradArr[istep];                         // *3.6; // Daymet srad in W/m^2
					   Qli = infrContArr[9].infdefValue;
					   Qnetob = infrContArr[10].infdefValue;
					   break;
				   case 2:
					   Qsiobs = SradArr[istep];                         // *3.6; // Daymet srad in W/m^2
					   Qli = LradArr[istep];
					   Qnetob = infrContArr[10].infdefValue;
					   break;
				   case 3:
					   Qsiobs = infrContArr[8].infdefValue;
					   Qli = infrContArr[9].infdefValue;
					   Qnetob = NradArr[istep];
					   break;
				   default:
					   cout<<" The radiation flag is not the right number; must be between 0 and 3"<<endl;
					   getchar();
					   break;				
				   }
				   //atm. pressure from netcdf 10.30.13    //this needs revision 		//####TBC_6.20.13
				   if (infrContArr[7].infType == 2)
					   sitev[1] = infrContArr[7].infdefValue;
				   else
					   sitev[1] = ApresArr[istep];
				   //this needs revision 		//####TBC_6.20.13
				   if (infrContArr[11].infType == 2)
					   Qg = infrContArr[11].infdefValue;
				   else
					   Qg = QgArr[istep];
				   //!     Flag to control albedo (ireadalb)  				 
				   if (infrContArr[12].infType == 2)
					   Snowalb = infrContArr[12].infdefValue;
				   else
					   Snowalb = SnowalbArr[istep];
				   //12.18.14 Vapor pressure of air
				   if (infrContArr[6].infType == 2)
					   Vp = infrContArr[6].infdefValue;
				   else
					   Vp = VpArr[istep];
				   //relative humidity computed or read from file
				   //#12.18.14 needs revision
				   if (infrContArr[5].infType == 2)
				   {
					   RH = infrContArr[5].infdefValue;
				   }
				   else if (infrContArr[5].infType == -1)          //RH computed internally 
				   {
					   float eSat = 611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
					   RH = Vp / eSat;
				   }
				   else
					   RH = RhArr[istep];
				   if (RH > 1)
				   {
					   //cout<<"relative humidity >= 1 at time step "<<istep<<endl;
					   RH = 0.99;
				   }
			   }   
			   else		//inputs are given as AVERAGE daily values, precip unit is always m/hr, Tmax and Tmin are required
			   {
				   //average daily value of precip in units of m/hr 4.22.14
				   P = PrecArr[istep / nstepinaDay];            // /24000 #12.19.14 Daymet prcp in mm/day           
				   V = WindspArr[istep / nstepinaDay];
				    
				   //get min max temp
				   Tmin = TaminArr[istep / nstepinaDay];
				   Tmax = TamaxArr[istep / nstepinaDay];
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
					   Qsiobs = infrContArr[8].infdefValue;
					   Qli = infrContArr[9].infdefValue;
					   Qnetob = infrContArr[10].infdefValue;
					   break;
				   case 1:
					   Qsiobs = SradArr[istep / nstepinaDay];                 // *3.6; // Daymet srad in W/m^2
					   Qli = infrContArr[9].infdefValue;
					   Qnetob = infrContArr[10].infdefValue;
					   break;
				   case 2:
					   Qsiobs = SradArr[istep / nstepinaDay];                 // *3.6; // Daymet srad in W/m^2
					   Qli = LradArr[istep / nstepinaDay];
					   Qnetob = infrContArr[10].infdefValue;
					   break;
				   case 3:
					   Qsiobs = infrContArr[8].infdefValue;
					   Qli = infrContArr[9].infdefValue;
					   Qnetob = NradArr[istep / nstepinaDay];
					   break;
				   default:
					   cout<<" The radiation flag is not the right number; must be between 0 and 3"<<endl;
					   getchar();
					   break;				
				   }
				   //atm. pressure from netcdf 10.30.13
				   //this needs revision 		//####TBC_6.20.13
				   if (infrContArr[7].infType == 2)
					   sitev[1] = infrContArr[7].infdefValue;
				   else
					   sitev[1] = ApresArr[istep / nstepinaDay];

				   //this needs revision 		//####TBC_6.20.13
				   if (infrContArr[11].infType == 2)
					   Qg = infrContArr[11].infdefValue;
				   else
					   Qg = QgArr[istep / nstepinaDay];
				   //!     Flag to control albedo (ireadalb)  
				   //!     0 is no measurements - albedo estimated internally
				   //!     1 is albedo read from file (provided: measured or obtained from another model)
				   //these need revision //####TBC_6.20.13
				   if (infrContArr[12].infType == 2)
					   Snowalb = infrContArr[12].infdefValue;
				   else
					   Snowalb = SnowalbArr[istep / nstepinaDay];
				   //12.18.14 Vapor pressure of air
				   if (infrContArr[6].infType == 2)
					   Vp = infrContArr[6].infdefValue;
				   else
					   Vp = VpArr[istep / nstepinaDay];
				   //relative humidity computed or read from file 			   //#12.18.14 
				   if (infrContArr[5].infType == 2)
				   {
					   RH = infrContArr[5].infdefValue;
				   }
				   else if (infrContArr[5].infType == -1)          //RH computed internally 
				   {
					   float eSat = 611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
					   RH = Vp / eSat;
				   }
				   else
					   RH = RhArr[istep / nstepinaDay];
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
				
				SNOWUEB2();     
		  
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

				if (snowdgt_outflag == 1 && indx ==0 )        //if debug mode 
				{
	            	printf(" time step: %d\n", istep);
					for (int uit = 0; uit<3; uit++)
						printf(" %d   ", (int) OutVarValues[uit][istep]);
					for(int uit = 3; uit< 70; uit++)
						printf(" %16.4f  ", OutVarValues[uit][istep]);
					printf(" \n");				
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
	   //copy next time step
	   modelStartDate[0] = Year;
	   modelStartDate[1] = Month;
	   modelStartDate[2] = Day;
	   modelStartHour = dHour;
	   //delete []tsprevday;
	   //delete []taveprevday;     
    
    return; 
}