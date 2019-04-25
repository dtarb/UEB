#include "gpuuebpgdecls.h"
 
uebCell::uebCell(const char* inpFile, int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT,
	double  UTCoffset, int inpDailyorSubd, int oStride)
{
	float siteVarInitcondefaults[32] = { 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0,
		0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903,
		4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936,
		9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 };
	float paramArr[32];
	setConstantValues();
	//set parameters	
	readParams(inpFile, paramArr);
	/*cout<<"param read..\n ");
	for(int i=0;i<npar;i++)
	cout<<"%f ",parvalArray[i]);	*/
	setParams(paramArr);	
	//model run settings
	setModelRun_Settings(startDate, endDate, startHour, endHour, modeldT, UTCoffset, inpDailyorSubd, oStride);
	//site vars and intitial conditions 
	setSiteVars_and_Initconds(siteVarInitcondefaults);
        //cout<<"UEBCell initialized"<<endl;	
}

uebCell::uebCell()
{
	float paramDefaults[32] = {0, 0, 3, -1, 0.98, 2.09, 2, 0.01, 337, 1700, 0.05, 20, 0.1, 0.85,
		0.65, 0.278, 1.11, 0.0654, 1, 0.001, 0.98, 0.5, 0, 0.5, 0.004626286,
		0.25, 0.5, 0.857143, 0.16, 0.5, 0.8, 2.4 };
	float siteVarInitcondefaults[32] = { 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0,
		0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903,
		4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936,
		9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 };
	setConstantValues();
	setParams(paramDefaults);	
	//model run settings
	//Default run settings 
	int startDate[3] = { 2009, 10, 1 }, endDate[3] = { 2010, 6, 1 };
	double startHour = 0.0, endHour = 0.0, modeldT = 1.0, UTCoffset = -7;
	int inpDailyorSubd = 0, outStride = 4;     //0 subdaily input
	setModelRun_Settings(startDate, endDate, startHour, endHour, modeldT, UTCoffset, inpDailyorSubd, outStride);
	//site vars and intitial conditions 
	setSiteVars_and_Initconds(siteVarInitcondefaults);
	//cout<<"UEBCell initialized"<<endl;	
}
/*
uebCell::uebCell(uebCell& uCell0)
{
	setConstantValues();
	setParams(uCell0.paramValues);
	//site vars and intitial conditions 
	setSiteVars_and_Initconds(uCell0.statesiteValues);
	//model run settings	
	setModelRun_Settings(uCell0.modelStartDate, uCell0.modelEndDate, uCell0.modelStartHour, uCell0.modelEndHour, uCell0.modelDT, uCell0.UTCOffset, uCell0.inpDailyorSubdaily, uCell0.outtStride);
	nstepinaDay = uCell0.nstepinaDay;
	if (uCell0.tsprevday)
	{
		tsprevday = new float[nstepinaDay];
		//  Initialize Tsbackup and TaveBackup
		for (int i = 0; i < nstepinaDay; i++)
		{
			tsprevday[i] = uCell0.tsprevday[i];
		}
	}
	else tsprevday = NULL;
	if (uCell0.taveprevday)
	{
		taveprevday = new float[nstepinaDay];
		for (int i = 0; i< nstepinaDay; i++)
		{
			taveprevday[i] = uCell0.taveprevday[i];
		}
	}
	else taveprevday = NULL;

}

uebCell& uebCell::operator= (uebCell& uCell0)
{
	if (this != &uCell0)
	{   
		setConstantValues();
		setParams(uCell0.paramValues);
		//site vars and intitial conditions 
		setSiteVars_and_Initconds(uCell0.statesiteValues);
		//model run settings	
		setModelRun_Settings(uCell0.modelStartDate, uCell0.modelEndDate, uCell0.modelStartHour, uCell0.modelEndHour, uCell0.modelDT, uCell0.UTCOffset, uCell0.inpDailyorSubdaily, uCell0.outtStride);
		nstepinaDay = uCell0.nstepinaDay;
		delete[] tsprevday;
		delete[] taveprevday;
		if (uCell0.tsprevday)
		{
			tsprevday = new float[nstepinaDay];
			//  Initialize Tsbackup and TaveBackup
			for (int i = 0; i < nstepinaDay; i++)
			{
				tsprevday[i] = uCell0.tsprevday[i];
			}
		}
		else tsprevday = NULL;
		if (uCell0.taveprevday)
		{
			taveprevday = new float[nstepinaDay];
			for (int i = 0; i< nstepinaDay; i++)
			{
				taveprevday[i] = uCell0.taveprevday[i];
			}
		}
		else taveprevday = NULL;
	}
	
	return *this;
}*/

uebCell::~uebCell()
{
	/*delete []tsprevday;
	delete []taveprevday;
	tsprevday = NULL;
	taveprevday = NULL;*/
	/*for(int i=0; i<72;i++)
		delete[] OutVarValues[i];
	delete []OutVarValues;*/
}

__host__ __device__ void uebCell::setConstantValues()
{
	//defalut is accumulation zone is false
	accumulationZone = false;
	//initialize 
	T_0 = 0.0;				// Temperature of freezing (0 C)
	T_k = 273.15;			// Temperature to convert C to K (273.15)
	SB_c = 2.041334e-7;      // Stefan boltzman constant (2.041334e-7 KJ/m^2-hr-K^4) #corrected 12.23.14
	H_f = 333.5;			// Heat of fusion (333.5 KJ= kg)
	Hne_u = 2834.0;		// Heat of Vaporization (Ice to Vapor, 2834 KJ= kg)
	C_w = 4.18;			// Water Heat Capacity (4.18 KJ/ kg-C)
	C_s = 2.09;			// Ice heat capacity (2.09 KJ= kg= C)
	C_p = 1.005;			// Air Heat Capacity (1.005 KJ= kg= K)
	Ra_g = 287.0;			// Ideal Gas constant for dry air (287 J= kg= K)
	K_vc = 0.4;			// Von Karmans constant (0.4)
	Hs_f = 3600.0;		// Factor to convert = s into = hr (3600)
	Rho_i = 917.0;			// Density of Ice (917 kg= m^3)
	Rho_w = 1000.0;		// Density of Water (1000 kg= m^3)
	Gra_v = 9.81;			// Gravitational acceleration (9.81 m= s^2)
	W1da_y = 0.261799;		// Daily frequency (2pi= 24 hr 0.261799 radians= hr) 
	Io = 4914.0;            //  Solar constant  Kj/m^2/hr
	//pi copied from snowdxv.f90
	P_i = 3.141592653589793238462643383279502884197169399375105820974944592308;		// Pi
	//data for pred-corr
	wtol = 0.025;
	utol = 2000.0;
	//from TURBFLUX()
	tol = 0.001;
	nitermax = 20;
	ncitermax = 21;
	// flag to write warnings,...etc
	snowdgtvariteflag = 0;
	snowdgtvariteflag2 = 0;  // 0;
	snowdgtvariteflag3 = 0;
	snowdgt_outflag = 0;
	radwarnflag = 0;
	inpDailyorSubdaily = 0;       // 0: values given at each (sub-daily time steps); 1: daily values
	uebCellX = 0;
	uebCellY = 0;
	Tsk_save = 273.16, Tssk_old = 273.16, Tsavek_old = 273.16, Tsavek_ave = 273.16, Tssk_ave = 273.16/*added 6.7.13*/;

	//## these were copied from snowdgtv, not clear where they are being used 
	/*fStab = -9999;
	Tref = -9999;
	iTsMethod = 4;
	//#_This is not clear 8.28.13
	windfl = 0;	*/
	return;
}

// functions to read params
void uebCell::readParams(const char* inpFile, float Params[32])
{
	ifstream pinFile(inpFile);
	char headerLine[256];	
	pinFile.getline(headerLine, 256, '\n'); //skip header
	for (int i = 0; i < 32; i++)
	{
		pinFile.getline(headerLine, 256, '\n');
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%f ", &Params[i]);
	}
	pinFile.close();
	return;
}

__host__ __device__ void uebCell::setParams(float Params[32])
{
	//copy params class variables
	for (int i = 0; i < 32; i++)
		paramValues[i] = Params[i];
	/*irad = (int) Params[0];	
    ireadalb=(int)Params[1];
	Tr = Params[2];     //  Temperature above which all is rain [3 C];
    Ts = Params[3];     //  Temperature below which all is snow [-1 C];
    Ems= Params[4];    //  emmissivity of snow [nominally 0.99];
    Cg = Params[5];    //  Ground heat capacity [nominally 2.09 KJ/kg/C];
    z  = Params[6];      //  Nominal meas. height for air temp. and humidity [2m];
    Zo = Params[7];     //  Surface aerodynamic roughness [m];
    Rho = Params[8];    //  Snow Density [Nominally 450 kg/m^3];
    Rhog = Params[9];   //  Soil Density [nominally 1700 kg/m^3];
    Lc = Params[10];     //  Liquid holding capacity of snow [0.05];
    Ks = Params[11];    //  Snow Saturated hydraulic conductivity [20 m/hr];
    De = Params[12];    //  Thermally active depth of soil [0.1 m];
	Avo=Params[13];   //  Visual new snow albedo [0.95];
    Anir0=Params[14]; //  NIR new snow albedo [0.65];
	Lans= Params[15]; //  the thermal conductivity of fresh [dry]; snow [0.0576 kJ/m/k/hr]; [Vinod// 0.36 Ref: Snow and Climate :Richard L Armstrong and Eric Brun ]; 
	Lang= Params[16]; //  the thermal conductivity of soil [:9.68 kJ/m/k/hr]; [TK of ice or wet soil[2.22~ 3.48W/m/k];:Vinod];
	Wlf= Params[17];  //  Low frequency fluctuation in deep snow/soil layer [1/4 w1 = 0.0654 radian/hr]; 
	Rd1= Params[18];  //  Apmlitude correction coefficient of heat conduction [1];
	dNewS=Params[19]; //  The threshold depth of for new snow [0.001 m];
	// 7 Parameters added for canopy
	EmC      = Params[20];			// Emissivity of canopy 
	Alpha    = Params[21];			// Scattering coefficient for solar radiation
	AlphaL   = Params[22];		    // Scattering coefficient for long wave radiation
	Gpar        = Params[23];            // leaf orientation with respect to zenith angle
    Uc       = Params[24];		    // Unloading rate coefficient [Per hour]; [Hedstrom and pomeroy; 1998]; 
	As       = Params[25];			// Fraction of extraterrestaial radiation on cloudy day;Shuttleworth [1993];  
	Bs       = Params[26];		  	// [as+bs];:Fraction of extraterrestaial radiation on clear day; Shuttleworth [1993]; 
	Lambda   = Params[27];		    // Ratio of direct atm radiation to diffuse;worked out from Dingman [1993];
	Rimax    = Params[28];            // Maximum value of Richardsion number for stability corretion
	Wcoeff   = Params[29];            // Wind decay coefficient for the forest
	Bca    = Params[30];              //A in Bristow-Campbell formula for atmospheric transmittance
	Bcc  =  Params[31];*/
//5.2.15 from snowdv
	//  Mapping from parameters read to UEB internal interpretation which follows UEBVeg scheme
	irad = (int)paramValues[0];
	ireadalb = (int)paramValues[1];
	for (int i = 0; i<11; i++)
		Param[i] = paramValues[i + 2];
	for (int i = 12; i<18; i++)
		Param[i] = paramValues[i + 1];
	Param[18] = -9999;
	Param[19] = -9999;
	Param[20] = paramValues[19];
	for (int i = 22; i<32; i++)
		Param[i] = paramValues[i - 2];
	bca = paramValues[30];
	bcc = paramValues[31];
	return;
}

//copy site variables and state intitial conditions at a grid (ueb cell)
__host__ __device__ void uebCell::setSiteVars_and_Initconds(float SiteVars[32])
{
	for (int i = 0; i < 32; i++)
		statesiteValues[i] = SiteVars[i];
	//copy initial conditions
	/*Usic = SiteVars[0];
	Wsic = SiteVars[1];
	dlSageic = SiteVars[2];
	Wcic = SiteVars[3];
	//Tcic = SiteVars[4];
	refDepthic = SiteVars[4];
	totalRefDepthic = SiteVars[5];
	Qg = SiteVars[6];
	//copy site variables at each grid poitn
	dF = SiteVars[7],     //  Drift factor
	APr = SiteVars[8],     //  Atmospheric Pressure [Pa],   
	Aep = SiteVars[9],     //  Albedo extinction parameter to smooth
	//  transition of albedo when snow is shallow. Depends on Veg. height [m],
	// 7 Site Variables added for canopy
	Cc = SiteVars[10],   // Canopy Coverage
	Hcan = SiteVars[11],   // Canopy height  
	LAI = SiteVars[12],   // Leaf Area Index  
	Sbar = SiteVars[13],  // Maximum snow load held per unit branch area[Kg/m2 for Pine],
	Ycage = (int)SiteVars[14];  // Parameters for wind speed transformation
						// Ycage=1 for young pine   Should be in parameters
						// Ycage=2 for Leafed deciduous
						// Ycage=3 for Old pine with logn stems (Paw U and Meyers, 1987)-- Requires for wind speed transformation	
	Slope = SiteVars[15];      //slope
	Azi = SiteVars[16];      //aspect
	lat = SiteVars[17];      //latitude
	Subalb = SiteVars[18];    //Substrate albedo
	Subtype = (int)SiteVars[19]; //substrate type beneath snow
	//0 = Ground/Non Glacier, 1=Clean Ice/glacier,
	//2= Debris covered ice/glacier, 3= Glacier snow accumulation zone
	Gsurf = SiteVars[20];     //The fraction of surface melt that runs off (e.g. from a glacier)
	ts_last = SiteVars[33];     //????
	lon = SiteVars[34];       //Longitude	
	for (int i = 0; i<12; i++)         //Bristow-Campbell B for each month
		dtBcb[i] = SiteVars[i + 21];*/
//==============================================changes for new conf 5.1.15
	//from snowdv
	//copied from paramsiteinitial
	for (int i = 0; i<4; i++)
		statev[i] = statesiteValues[i];
	sitev[0] = statesiteValues[4];
	sitev[1] = statesiteValues[5];
	for (int i = 3; i<9; i++)
		sitev[i] = statesiteValues[i + 3];
	slope = statesiteValues[12];
	azi = statesiteValues[13];
	lat = statesiteValues[14];

	Param[11] = statesiteValues[15];
	//subalb=statesiteValues[15]
	sitev[9] = statesiteValues[16];
	subtype = (int)statesiteValues[16];
	Param[21] = statesiteValues[17];
	//gsurf = statesiteValues[17]
	for (int i = 0; i<12; i++)
		dtbar[i] = statesiteValues[i + 18];
	ts_last = statesiteValues[30];
	lon = statesiteValues[31];
	
	if (subtype == 0 || subtype == 3)
		WGT = 0.0;
	else
		WGT = 1.0;

	if (subtype != 3)  //  Only do this work for non accumulation cells where model is run
	{
		//  Initialize Tsbackup and TaveBackup
		for (int i = 0; i< nstepinaDay; i++)
		{
			tsprevday[i] = -9999.0;
			taveprevday[i] = -9999.0;
		}
		// Take surface temperature as 0 where it is unknown the previous time step
		// This is for first day of the model to get the force restore going
		//#$#$#$#$#_is this all the time steps or the last time?
		if (ts_last <= -9999)
			//for(int i =0;i< nstepinaDay;i++)
			tsprevday[nstepinaDay - 1] = 0;
		else
			//for(int i =0;i< nstepinaDay;i++)
			tsprevday[nstepinaDay - 1] = ts_last;
		// compute Ave.Temp for previous day 
		Us = statev[0];                                  // Ub in UEB
		Ws = statev[1];                                  // W in UEB
		Wc = statev[3];                                  // Canopy SWE
		Apr = sitev[1];                                   // Atm. Pressure  [PR in UEB]
		cg = Param[3];                                   // Ground heat capacity [nominally 2.09 KJ/kg/C]
		rhog = Param[7];                                   // Soil Density [nominally 1700 kg/m^3]
		de = Param[10];                                  // Thermally active depth of soil (0.1 m)

		//this are for coudiness computation
		//6.10.13
		as = Param[27];
		bs = Param[28];

		tave = TAVG(Us, Ws + WGT, Rho_w, C_s, T_0, rhog, de, cg, H_f);
		//for(int i =0;i< nstepinaDay;i++)
		taveprevday[nstepinaDay - 1] = tave;
		//  initialize variables for mass balance
		Ws1 = statev[1];
		Wc1 = statev[3];
		cumP = 0.0;
		cumEs = 0.0;
		cumEc = 0.0;
		cumMr = 0.0;
		cumGm = 0.0;
		cumEg = 0.0;
	} //  end the skip block done only for accumulation cells	
	Tmin = 0.0;
	Tmax = 0.0;

	return;
}

__host__ __device__ void uebCell::setModelRun_Settings(int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT, double  UTCoffset, int inpDailyorSubd, int oStride)
{
	for (int i = 0; i<3; i++)
	{
		modelStartDate[i] = startDate[i];
		modelEndDate[i] = endDate[i];
	}
	modelStartHour = startHour;
	modelEndHour = endHour;
	modelDT = modeldT;	
	UTCOffset = UTCoffset;
	inpDailyorSubdaily = inpDailyorSubd;
	//model time steps
	modelSpan = julian(modelEndDate[0], modelEndDate[1], modelEndDate[2], modelEndHour) - julian(modelStartDate[0], modelStartDate[1], modelStartDate[2], modelStartHour); //no of days in model span
	numTimeStep = (int)ceil(modelSpan*(24 / modelDT));	
//5.2.15 from snowdv
	//  FIXME: what if the result is fractional
	//  time steps must divide exactly in to a day because we use logic that requires the values from the same time
	//  step on the previous day.  Consider in future making the specification of time step as number of time
	//  steps in a day, not modeldt to ensure this modeldt is recalculated based on the int timesteps in a day
	//  assumption: number of model timesteps in a day must be an int             
	/*stepinaDay= (int) (24.0/modelDT +0.5);  // closest rounding
	modelDT = 24.0/stepinaDay;
	nstepinaDay = stepinaDay;
	tsprevday = new float[nstepinaDay];
	taveprevday = new float[nstepinaDay];   */

	//  assumption: number of model timesteps in a day must be an int             
	nstepinaDay = (int)(24.0 / modelDT + 0.5);  // closest rounding
	modelDT = 24.0 / nstepinaDay;
	if (inpDailyorSubdaily == 0)
		numSimTimeSteps = 24;
	else
		numSimTimeSteps = 24 * nstepinaDay;

	for (int i = 0; i < 13; i++)
	{
		startIndex[i] = 0;
		ncReadStart[i] = 0;
	}
	tEnd = 0;
	outtStride = oStride;
	timeSeriesIndex = 0;                  //this changes to 1 when a forcing that is applicable for the whole model is read - --- forcing time series from text file need to be read only once
	//allocate memory for output array	
	//OutVarValues = new float *[70];
	//for (int i = 0; i < 70; i++)
	//OutVarValues = new float[70*numTimeStep];         // *outtStride];
	//cout << "number of t " << numTimeStep << endl;
	/*tsprevday.clear();
	tsprevday.resize(nstepinaDay);
	taveprevday.clear();
	taveprevday.resize(nstepinaDay);
	for (int i = 0; i < nstepinaDay; i++)
	{
		tsprevday[i] = -9999.0;		
		taveprevday[i] = -9999.0;
	}	*/
// 5.2.15 from snowdv
	// calculating model end date-time in julian date
	dHour = modelEndHour;
	EJD = julian(modelEndDate[0], modelEndDate[1], modelEndDate[2], dHour);


	return;
}
// function  to read forcing / weather variables control file
void uebCell::readInputForContr(const char* inputconFile)
{
	ifstream pinFile(inputconFile);
	char headerLine[256];
	//istringstream valueLine;	
	pinFile.getline(headerLine, 256);   //skip header	
	for (int i = 0; i<13; i++)
	{
		pinFile.getline(headerLine, 256, ':');
		sscanf(headerLine, "%s ", &infrContArr[i].infName);
		pinFile.getline(headerLine, 256, '\n');
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%d ", &infrContArr[i].infType);
		//headerLine[0] = 0;
		//fscanf(pinFile,"%d\n",&svArr[i].svType);
		switch (infrContArr[i].infType)
		{
		case -1:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &infrContArr[i].infdefValue);
			break;
		case 0:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%s ", &infrContArr[i].infFile);
			break;
		case 1:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%s %s %s %d", &infrContArr[i].infFile, &infrContArr[i].infvarName, &infrContArr[i].inftimeVar, &infrContArr[i].numNcfiles);
			break;
		case 2:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &infrContArr[i].infdefValue);
			break;
		default:
			cout << "Wrong input/forcing type; has to be -1 (compute by the model), 2 (single value) , 0 (time-series text file) or 1 (3D netcdf)" << endl;
			cout << "Using default value..." << endl;
			break; //exit(1); 
		}
		//i++;
		//headerLine[0] = 0;
		//}
	}
	pinFile.close();
	return;
}

void uebCell::getInpForcArr(int numNc[13], float*** RegArray[13], int ncTotaltimestep[13], MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	for (int it = 0; it < 13; it++)
	{
		if (infrContArr[it].infType == 0)
		{
			/*if (numNc == 0)  // for time series from text file read once ---- outside of this function
			{
				RegArray[it] = new float**[1];
				RegArray[it][0] = new float *[1];
				readTStextFile(infrContArr[it].infFile, RegArray[it][0][0], ncTotaltimestep[it]);   //tsvarArray[it][0] 3.19.15      ntimesteps[0] 12.18.14
			}*/
		}
		else if (infrContArr[it].infType == 2 || infrContArr[it].infType == -1)
		{
			// use default value or compute internally		
		}
		else if (infrContArr[it].infType == 1)     // == 0
		{
			tEnd = ncReadStart[it] + 24; 
			//offSet = 1;    // uebCellY*dimLen2*numTimeStep + uebCellX*numTimeStep;
			int retvalue = 0;			
			//read 3D netcdf (regridded array processed by uebInputs)
			char numtoStr[256];
			sprintf(numtoStr, "%d", numNc[it]);
			char tsInputfile[256];
			strcpy(tsInputfile, infrContArr[it].infFile);
			strcat(tsInputfile, numtoStr);
			strcat(tsInputfile, ".nc");
			//cout<<"%s\n",tsInputfile);
			//clear existing memory RegArray[it] before passing to this function // delete[] RegArray[it];
			readNC_yxSlub(tsInputfile, infrContArr[it].infvarName, infrContArr[it].inftimeVar, ncReadStart[it], tEnd, RegArray[it], ncTotaltimestep[it], numNc[it], inpComm, inpInfo);
			//startIndex[it] = 0;
			//endIndex[it] = ncTotaltimestep[it];
			//cout << "nc time = " << ncNtimestes[it][numNc];			
		}
	}	
}

void uebCell::updateInpForcArr(float*** RegArray[13], int ncTotaltimestep[13])
{
	setInpForcArr(0, RegArray[0], PrecArr, ncTotaltimestep[0]);
	setInpForcArr(1, RegArray[1], TempArr, ncTotaltimestep[1]);
	setInpForcArr(2, RegArray[2], TaminArr, ncTotaltimestep[2]);
	setInpForcArr(3, RegArray[3], TamaxArr, ncTotaltimestep[3]);
	setInpForcArr(4, RegArray[4], WindspArr, ncTotaltimestep[4]);
	setInpForcArr(5, RegArray[5], RhArr, ncTotaltimestep[5]);
	setInpForcArr(6, RegArray[6], VpArr, ncTotaltimestep[6]);
	setInpForcArr(7, RegArray[7], ApresArr, ncTotaltimestep[7]);
	setInpForcArr(8, RegArray[8], SradArr, ncTotaltimestep[8]);
	setInpForcArr(9, RegArray[9], LradArr, ncTotaltimestep[9]);
	setInpForcArr(10, RegArray[10], NradArr, ncTotaltimestep[10]);
	setInpForcArr(11, RegArray[11], QgArr, ncTotaltimestep[11]);
	setInpForcArr(12, RegArray[12], SnowalbArr, ncTotaltimestep[12]);
}

void uebCell::setInpForcArr(int it, float ***inArray, float* forcArr, int ncTotaltimestepit)
{
	//need to call each variable array as each array has to be copied separately to device array in cuda
	int tsLength = 24; //default length	
	if (infrContArr[it].infType == 0)
	{
		if (ncTotaltimestepit - startIndex[it] < tsLength)
			tsLength = ncTotaltimestepit - startIndex[it];  //make sure not to go out of array bound
		if (inpDailyorSubdaily == 0)
		{
			if (numSimTimeSteps > tsLength)
				numSimTimeSteps = tsLength;            // use the smallest number of time steps
		}
		else
		{
			if (numSimTimeSteps > (tsLength * nstepinaDay))
				numSimTimeSteps = tsLength * nstepinaDay;            // use the smallest number of time steps

		}
		for (int i = 0; i < tsLength; i++)
			forcArr[i] = inArray[0][0][startIndex[it] + i];
		startIndex[it] += tsLength;
	}
	else if (infrContArr[it].infType == 2 || infrContArr[it].infType == -1)
	{
		// use default value or compute internally		
	}
	else if (infrContArr[it].infType == 1)     // == 0
	{
		if (ncTotaltimestepit - startIndex[it] < tsLength)
			tsLength = ncTotaltimestepit - startIndex[it];  //make sure not to go out of array bound
		if (inpDailyorSubdaily == 0)
		{
			if (numSimTimeSteps > tsLength)
				numSimTimeSteps = tsLength;            // use the smallest number of time steps
		}
		else
		{
			if (numSimTimeSteps > (tsLength * nstepinaDay))
				numSimTimeSteps = tsLength * nstepinaDay;            // use the smallest number of time steps

		}
		for (int i = 0; i < tsLength; i++)
			forcArr[i] = inArray[uebCellY][uebCellX][startIndex[it] + i];
		startIndex[it] = 0;       // +tsLength;
	}
}

//print all output values at a point
void uebCell::printPointOutputs(const char* outFileName)
{
	FILE* outFile = fopen(outFileName,"a");             //can write multiple times appending at the end
	for (int istep = 0; istep < numSimTimeSteps - 1; istep++)         //-2 to be safe againts ceil( ) in computeModelDateTime()
	{
		fprintf(outFile, " %d %d %d %8.3f ", (int)OutVarValues[0][istep], (int)OutVarValues[1][istep], (int)OutVarValues[2][istep], OutVarValues[3][istep]);
		for (int vnum = 4; vnum <70; vnum++)
			fprintf(outFile, " %16.4f ", OutVarValues[vnum][istep]);
	}
	fclose(outFile);
}
//print values at a point for degugging
void uebCell::printDebugOutputs()
{
	char testPrint[256];
	char ind[256];
	strcpy(testPrint, "ZTest");
	sprintf(ind, "%d", uebCellY);
	strcat(testPrint, ind);
	strcat(testPrint, "_");
	sprintf(ind, "%d", uebCellX);
	strcat(testPrint, ind);
	strcat(testPrint, ".txt");	
	FILE* outFile = fopen(testPrint, "a");
	for (int istep = 0; istep < numSimTimeSteps - 1; istep++)         //-2 to be safe againts ceil( ) in computeModelDateTime()
	{
		fprintf(outFile, " %d %d %d %8.3f ", (int)OutVarValues[0][istep], (int)OutVarValues[1][istep], (int)OutVarValues[2][istep], OutVarValues[3][istep]);
		for (int vnum = 4; vnum < 70; vnum++)
			fprintf(outFile, " %16.4f ", OutVarValues[vnum][istep]);
	}
	fclose(outFile);
}
//print SWE (snow water equivalent), Us (Energy content), P(recipitation) and Ta(Temperature) at a point
void uebCell::printSampleOutputs(const char* outFileName)
{
	FILE* outFile = fopen(outFileName,"a");
	for (int istep = 0; istep < numSimTimeSteps - 1; istep++)         //-2 to be safe againts ceil( ) in computeModelDateTime()	
		fprintf(outFile, "%d %d %d %8.3f %16.4f %16.4f %16.4f %16.4f %16.4f %16.4f\n", (int)OutVarValues[0][istep], (int)OutVarValues[1][istep], (int)OutVarValues[2][istep], OutVarValues[3][istep],
		OutVarValues[12][istep], OutVarValues[13][istep], OutVarValues[16][istep], OutVarValues[17][istep], OutVarValues[18][istep], OutVarValues[19][istep]);
	fclose(outFile);
}

__host__ __device__ int findMax(int a, int  b)
{
	return (a>b)?a:b;	

}

__host__ __device__ int findMin(int a, int b)
{
	return (a<b)?a:b;
}
__host__ __device__ float findMax(float a, float  b)
{
	return (a>b)?a:b;
}

__host__ __device__ float findMin(float a, float b)
{
	return (a<b)?a:b;
}