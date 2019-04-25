#include<iomanip>
#include<iostream>
#include <fstream>
#include <cstring>
#include<sstream>
#include <stdio.h>
#include "mpi.h"
#include <netcdf.h>
#include <netcdf_par.h>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#pragma warning(disable : 4996)
//#include <algorithm>
using namespace std;
//error handling for netcdf
#define  ERR(e) {cout<<"Error: "<< nc_strerror(e)<<endl; return 2; }
//error handling for cuda
//define  cuda_checkERR(err) { if (err != cudaSuccess) std::cout << "Error: "<<cudaGetErrorString(err)<< endl; exit(EXIT_FAILURE); }
__device__ __host__ void cuda_checkERR(cudaError_t err);

struct params {
	float irad, ireadalb, tr, ts, ems, cg, z, zo, rho, rhog, lc, ks, de, avo,
	anir0, lans, lang, wlf, rd1, dnews, emc, alpha, alphal, gpar, uc, as, Bs,
	lambda, rimax, wcoeff, apar, cpar;
};

struct sitevar {
	char svName[256];
	int   svType;
	char svFile[256];
	char svVarName[256];
	float svdefValue;
	float** svArrayValues;
};

struct inpforcvar {
	char infName[256];
	int  infType;
	char infFile[256];
	char infvarName[256];
	char inftimeVar[256];
	float infdefValue;
	int numNcfiles;
};
struct pointOutput {
	char outfName[256];
	int ycoord;
	int xcoord;
};

struct ncOutput {
	char outfName[256];
	char symbol[256];
	char units[256];
};

struct aggOutput {
	//char outfName[256];
	char symbol[256];
	char units[256];
	char aggop[256];
};

struct inptimeseries {
	//CTime dtime;
	float datetime;
	float tsValue;
};

class uebCell {  

	public:
		uebCell();
		uebCell(const char* inpFile, int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT,
			double  UTCoffset, int inpDailyorSubd, int oStride);
		//uebCell(uebCell& uCell0);
		//uebCell& operator= (uebCell& uCell0);
		~uebCell();
		//read parameter values
		__host__ __device__  void  setConstantValues();
		__host__ __device__  void  setDefaultParams();
		void  readParams(const char* inpFile, float Params[32]);
		 __host__ __device__  void  setParams(float Params[32]);
		//copy site variable array
		 __host__ __device__  void  setSiteVars_and_Initconds(float SiteVars[32]);
		 __host__ __device__  void  setModelRun_Settings(int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT, 
			                          double  UTCoffset, int inpDailyorSubd, int outtStride);		
		  //##### snowdv.cpp
		 __host__ __device__ int getforcOffset(int ioffst, int dimLen2);
		 __host__ __device__ void runUEB(int dimLen2);
		 __host__ __device__  void  runUEB();		 
		 void readInputForContr(const char* inputconFile);		
		 void getInpForcArr(int numNc[13], float*** RegArray[13], int ncTotaltimestep[13], MPI::Intracomm inpComm, MPI::Info inpInfo);
		 void updateInpForcArr(float*** RegArray[13], int ncTotaltimestep[13]);
		 void setInpForcArr(int it, float ***inArray, float* forcArr, int ncTotaltimestepit);
         //### 
		 void  printPointOutputs(const char* outFileName);
		 void  printDebugOutputs();
		 void  printSampleOutputs(const char* outFileName);

		//accumulation zone?
		bool accumulationZone;
		int modelStartDate[3], modelEndDate[3];
		double modelStartHour, modelEndHour, modelDT, UTCOffset, modelSpan;
		//int ModelStartDate[3], ModelEndDate[3], TstepsinDay, numTimeSteps; //check this 
		//double ModelStartHour, ModelEndHour, ModeldT, ModelUTCOffset, modelSpan;
		//track current grid cell
		int uebCellX;
		int uebCellY;
		//Global variables ###_TBC 5.6.13
		// flag to write warnings,...etc ###***TBC 9.20.13
		int snowdgtvariteflag;  // 0; print debug information
		int snowdgtvariteflag2;  // more detailed print in Surfebsc
		int snowdgtvariteflag3;  // print in predictor corrector function
		int snowdgt_outflag;        //output print
		int radwarnflag;
		int inpDailyorSubdaily;       // 0: values given at each (sub-daily time steps); 1: daily values
		int outtStride;
		int numTimeStep;
		int nstepinaDay;
		int numSimTimeSteps;
		int timeSeriesIndex;       //indicates whether a time series is read or not --- forcing time series from text file need to be read only once
		//float* tsvarArray[13]; 
		inpforcvar infrContArr[13]; 
		//int ncTotaltimestep[13];
		int startIndex[13], ncReadStart[13], tEnd;
		//outputs
		float OutVarValues[70][96];                     //4.28.15 compute max of 96 time steps at time (for daily input, 4 days with dT = 6h)
		//params and site and inp vars
		float statesiteValues[32], paramValues[32];
		//climate forcing arrays, static arrays of length 3264 (multiples of 32, 64, 24, ---data of about 4.5 months at hourly steps)
		float  PrecArr[24], TempArr[24], TaminArr[24], TamaxArr[24], WindspArr[24], RhArr[24],
			VpArr[24], ApresArr[24], SradArr[24], LradArr[24], NradArr[24], QgArr[24], SnowalbArr[24];

		float  tsprevday[24], taveprevday[24];   //assumes the dt >= 1hr which is generally the case in UEB

		//time related funcs
		//snowxv.cpp
		//********UPDATEtime ()  Update time for each time step
		__host__ __device__  void  UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT);
		//    function to return number of days in February checking for leap years
		__host__ __device__ int lyear(int year);		
		//To convert the real date to julian date
		__host__ __device__  int julian(int yy, int mm, int dd);		
		//these were copied from functions.f90
		//COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  
		__host__ __device__  double julian(int I, int M, int K, double H);
		//COMPUTES CALENDAR DATE AND time, GIVEN JULIAN DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
		__host__ __device__  void  calendardate(double TJD, int &I, int &M, int &K, double &H);

    private:
	    //variables to track previous day temperature profile
		//std::vector <float> tsprevday;
		//std::vector <float> taveprevday;
		//js  Constant    floatset 
		float T_0; 				// Temperature of freezing (0 C)
		float T_k;  // 273.15;			// Temperature to convert C to K (273.15)
		float SB_c;  // Stefan boltzman constant (2.041334e-7 KJ/m^2-hr-K^4) #corrected 12.23.14
		float H_f;  // 333.5;			// Heat of fusion (333.5 KJ;  // kg)
		float Hne_u;  // 2834.0;		// Heat of Vaporization (Ice to Vapor, 2834 KJ;  // kg)
		float C_w;  // 4.18;			// Water Heat Capacity (4.18 KJ;  // kg;  // C)
		float C_s;  // 2.09;			// Ice heat capacity (2.09 KJ;  // kg;  // C)
		float C_p;  // 1.005;			// Air Heat Capacity (1.005 KJ;  // kg;  // K)
		float Ra_g;  // 287.0;			// Ideal Gas constant for dry air (287 J;  // kg;  // K)
		float K_vc;  // 0.4;			// Von Karmans constant (0.4)
		float Hs_f;  // 3600.0;		// Factor to convert ;  // s into ;  // hr (3600)
		float Rho_i;  // 917.0;			// Density of Ice (917 kg;  // m^3)
		float Rho_w;  // 1000.0;		// Density of Water (1000 kg;  // m^3)
		float Gra_v;  // 9.81;			// Gravitational acceleration (9.81 m;  // s^2)
		float W1da_y;  // 0.261799;		// Daily frequency (2pi;  // 24 hr 0.261799 radians;  // hr) 
		float  Io;  // 4914.0;            //  Solar constant  Kj/m^2/hr
		//pi copied from snowdxv.f90
		float P_i;  // 3.141592653589793238462643383279502884197169399375105820974944592308;		// Pi
		//data for pred-corr
		float wtol;  // 0.025;
		float utol;  // 2000.0;
		//from TURBFLUX()
		float tol;  // 0.001;
		int  nitermax;  // 20;  
		int  ncitermax; //21      
		//common (surface and average temp of previous time step ? ### 9.20.13?)
		float Tsk_save;  // 273.16,
		float Tssk_old;  // 273.16,  
		float Tsavek_old;  // 273.16, 
		float Tsavek_ave;  // 273.16, 
		float Tssk_ave;  // 273.16/*added 6.7.13*/;
//==============================================changes for new conf 5.1.15
//from snowdv
		//const int nsv = 10, npar = 32, nxv = 6, niv = 8, nov = 14;
		float  Inpt[8], sitev[10], outv[14], statev[6], Param[32], dtbar[12], mtime[4];                                     // YJS pass to reflect the change of Snow (Year, month, Date, Hour)
		int irad, ireadalb, subtype, iradfl, iflag[6];
		float  slope, azi, bca, bcc, lat, ts_last, lon;
		//float *tsprevday, *taveprevday;
		//int stepinaDay, nstepinaDay; 
		//float referenceHour, referenceTime, CTJD;
		double dHour, EJD, MHour, sHour, currentModelDateTime, UTCHour, OHour, UTCJulDat;
		float fHour, fModeldt;
		int istep, Year, Month, Day, MYear, MMonth, MDay;

		// CHANGES TO ACCOMODATE GLACIER
		float  WGT; // WGT=WATER EQUIVALENT GLACIER THICKNESS
		float Us, Ws, Wc, Apr, cg, rhog, de, tave, Ws1, Wc1, cumP, cumEs, cumEc, cumMr, cumGm, cumEg, Ta, P, V, RH, Tmin, Tmax, Trange, Qsiobs, Qg, Qli, QLif;
		float Vp; //#12.18.14 Air vapor pressure 
		float Qnetob, Snowalb, cosZen, atff, cf, atfimplied, Ema, Eacl, dStorage, errMB, HRI0, HRI, as, bs;
		float  OutArr[53];
		//to get max/min daily temperatures
		int nb, nf, nbstart, nfend;
		double dstHour;

		//flags
		/*int irad, iradfl, ireadalb;       //rad flag: iradfl = 0: short and long wave radiations; 1: net radiation
										 //albedo flag: ireadalb = 0 albedo computed internally from dimensionless snow age (dlSage)
		                                 // 1: albedo copied from dlSage 
		                                // Warning: dlSage can contain two variables 
		int iTsMethod;                 // the method to approximate the surface temperature
										// 1 normal old snow melt model
										// 2 revised direct gradient method (New method for ke) and qe
										// 3 force restore approach
										// 4 modified force restore approach
		//#_Windflag???
	    int windfl;   */ 

		
		
		//parameters
	/*	float Tr,           // Params[2];     //  Temperature above which all is rain [3 C];
			Ts ,           // Params[3];     //  Temperature below which all is snow [-1 C];
			Ems,           // Params[4];    //  emmissivity of snow [nominally 0.99];
			Cg ,           // Params[5];    //  Ground heat capacity [nominally 2.09 KJ/kg/C];
			z  ,           // Params[6];      //  Nominal meas. height for air temp. and humidity [2m];
			Zo ,           // Params[7];     //  Surface aerodynamic roughness [m];
			Rho ,           // Params[8];    //  Snow Density [Nominally 450 kg/m^3];
			Rhog ,           // Params[9];   //  Soil Density [nominally 1700 kg/m^3];
			Lc ,           // Params[10];     //  Liquid holding capacity of snow [0.05];
			Ks ,           // Params[11];    //  Snow Saturated hydraulic conductivity [20 m/hr];
			De ,           // Params[12];    //  Thermally active depth of soil [0.1 m];
			Avo,           //Params[13];   //  Visual new snow albedo [0.95];
			Anir0,           //Params[14]; //  NIR new snow albedo [0.65];
			Lans,           // Params[15]; //  the thermal conductivity of fresh [dry]; snow [0.0576 kJ/m/k/hr]; [Vinod// 0.36 Ref: Snow and Climate :Richard L Armstrong and Eric Brun ]; 
			Lang,           // Params[16]; //  the thermal conductivity of soil [:9.68 kJ/m/k/hr]; [TK of ice or wet soil[2.22~ 3.48W/m/k];:Vinod];
			Wlf,           // Params[17];  //  Low frequency fluctuation in deep snow/soil layer [1/4 w1 ,           // 0.0654 radian/hr]; 
			Rd1,           // Params[18];  //  Apmlitude correction coefficient of heat conduction [1];
			dNewS,           //Params[19]; //  The threshold depth of for new snow [0.001 m];
			// 7 Parameters added for canopy
			EmC ,           // Params[20];	 // Emissivity of canopy 
			Alpha,          // Params[21];   // Scattering coefficient for solar radiation
			AlphaL,         // Params[22];   // Scattering coefficient for long wave radiation
			Gpar,           // Params[23];   // leaf orientation with respect to zenith angle
			Uc  ,           // Params[24];   // Unloading rate coefficient [Per hour]; [Hedstrom and pomeroy; 1998]; 
			As  ,           // Params[25];	 // Fraction of extraterrestaial radiation on cloudy day;Shuttleworth [1993];  
			Bs  ,           // Params[26];   // [as+bs];:Fraction of extraterrestaial radiation on clear day; Shuttleworth [1993]; 
			Lambda,         // Params[27];   // Ratio of direct atm radiation to diffuse;worked out from Dingman [1993];
			Rimax,          // Params[28];   // Maximum value of Richardsion number for stability corretion
			Wcoeff,         // Params[29];   // Wind decay coefficient for the forest
			Bca ,          // Params[30];   //A in Bristow-Campbell formula for atmospheric transmittance
			Bcc ;          // Params[31];   //C in Bristow-Campbell formula for atmospheric transmittance

		//## these were copied from snowdgtv, not clear where they are being used 
		int fStab;			//  Stability correction control parameter 0 = no corrections, 1 = full corrections
	    float Tref;			 //  Reference temperature of soil layer in ground heat calculation input

		//site variables
		float dF ,            // SiteVars[4],     //  Drift factor
			 Qg,        //ground heat flux....taken as site var as it is difficult to find temporally varying Qg (often 0)
			APr,            // SiteVars[5],     //  Atmospheric Pressure [Pa],   
			Aep,            // SiteVars[6],     //  Albedo extinction parameter to smooth transition of albedo when snow is shallow. Depends on Veg. height [m],
		// 7 Site Variables added for canopy
			Cc ,            // SiteVars[7],   // Canopy Coverage
			Hcan ,            // SiteVars[8],   // Canopy height  
			LAI ,            // SiteVars[9],   // Leaf Area Index  
			Sbar ;            // SiteVars[10],  // Maximum snow load held per unit branch area[Kg/m2 for Pine],
		int Ycage;           // (int) SiteVars[11];  // Parameters for wind speed transformation  Ycage:
			                              //1 for young pine   Should be in parameters
										  //2 for Leafed deciduous
										  //3 for Old pine with logn stems (Paw U and Meyers, 1987)--Requires for wind speed transformation
		float Slope,             // SiteVars[12];      //slope
			Azi ,               // SiteVars[13];      //aspect
			lat,               // SiteVars[14];      //latitude
			Subalb,               // SiteVars[15];    //Substrate albedo			
			Gsurf,               // SiteVars[17];     //The fraction of surface melt that runs off (e.g. from a glacier)
			ts_last,               // SiteVars[30];     //????
			lon;              // SiteVars[31];       //Longitude	
		int Subtype;          // (int)SiteVars[16]; //substrate type beneath snow
													//0  Ground/Non Glacier, 1=Clean Ice/glacier,
													//2  Debris covered ice/glacier, 3= Glacier snow accumulation zone			
		float dtBcb[12];     // Bristow-Campbell B for each month
	    float cosZen;         //Cosine of Zenith angle (Not read from file--computed at each time step)

		//state variables
		//float Us, Ws, Tss, Tsave, Wc, Tc, refDepth, totalRefDepth, dlSage, WGT;    //, Us1, Ws1, Wc1;       dlSage = dimenstionless age of snow surface    // WGT=WATER EQUIVALENT GLACIER THICKNESS
		//initial conditions (State variables)
		float Usic, Wsic, dlSageic, Wcic, refDepthic, totalRefDepthic; // Tcic,
		//float Us1, Ws1, Tss1, Tsave1, Wc1, Tc1, refDepth1, totalRefDepth1, dlSage1;   //next time step states
		 
		
		//atmospheric attenuation factor and cloud fraction
		//float Atff, Cf;

		//forcing (time-variable input) variables
		//float Ta, P, V, RH, Qg, Qsi, Qli, Qnet, Snowalb;*/

		//###check operator precedencde issue in albedo()and other functions _6.5.13

		//snowdv.cpp
		/* __host__ __device__  void  runUEB();  float* RegArray[13], float statesiteValues[], float paramValues[], float** &OutVarValues, int modelStartDate[],
					     double modelStartHour, int modelEndDate[], double modelEndHour, double modelDT, double UTCOffset);*/

		//snowxv.cpp		
		//    to get the atmospheric transmissivity using the Bristow and Campbell  (1984) approach
		 __host__ __device__  void  atf(float &atff,float trange,int month, float *dtbar, float a, float c);
		// To get hourly radiation index
		 __host__ __device__  void  hyri(int YEAR, int MONTH, int DAY, float HOUR, float DT, float SLOPE, float AZI, float LAT, float &HRI, float &COSZEN);		 
		//    Computes the incoming longwave radiation using satterlund Formula
		 __host__ __device__  void  cloud(float as, float bs, float atff, float &cf);
		//???? long wave radiation from temperatrue and other weather variables??
		//TBC_6.5.13
		 __host__ __device__  void  qlif(float TA, float RH, float TK, float SBC, float &Ema, float &Eacl, float cf, float &qliff );	

		//##### canopy.cpp       
		//     Partition the incoming solar radiation into direct and diffuse components
		 __host__ __device__  void  PSOLRAD( float Qsi, float atff, float *param, float cf,     
		                float &Taufb, float &Taufd, float &Qsib, float &Qsid);     // Output variables
		//     Estimates the direct and diffuse solar radiation fractions transmitted through the canopy      
		 __host__ __device__  void  TRANSRADCAN (float COSZEN, float *sitev, float *param,       
		     	            float &Betab, float &Betad, float &Taub, float &Taud);                     // Output variables: Transmission and Reflection fractions
		//      Computes the net canopy and sub-canopy solar radiation
		//      considering the multiple scatterings of radiation by the canopy and multile reflections 
		//      between the canopy and the snow surface
		 __host__ __device__  void  NETSOLRAD(float Ta,float A,float Betab,float Betad,float Wc,float Taub,float Taud, 
					     float Qsib,float Qsid,float *param, float	Fs,float &Qsns,float &Qsnc); //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
		//     Computes the net canopy and beneath canopy longwave radiation
		 __host__ __device__  void  NETLONGRAD(float RH,float Ta,float Tss,float Tc,float Tk,float Fs,float EmC,float EmS,float SBC,float cf,float *sitev,float Qli,float *param,
						    float &Qlis, float &Qlns,float &Qlnc );
		//****************  BULK AERODYNAMIC AND CANOPY BOUNDARY LAYER RESISTANCES  *******************
		//      Calculates resistances for the above and beneath canopy turbulent heat fluxes 
		 __host__ __device__  void  AeroRes(float P, float Wc, float V, float Ta, float Tss, float Tc, float Fs, float *param, float *sitev, float Tk,
					    float &d, float &Z0c, float &Vz, float &Rc, float &Ra, float &Rbc, float &Rl, float &RKINc, float &RKINa, float &RKINbc, float &RKINl);         // Output variables
		//     Calculates wind at the 2 m above the surface and at the sink using the input of measured wind at 2 m above the canopy
		 __host__ __device__  void  WINDTRANab(float V, float Zm, float *param, float *sitev,float &Vstar,float &Vh,float &d,float &Z0c,float &Vz,float &Vc);               // Out put wind speed within canopy at height 2 m
		//     Calculates the turbulent heat fluxes (sensible and latent heat fluxes) and condensation/sublimation.
		 __host__ __device__  void  TURBFLUX (float Ws, float Wc, float A, float Tk,float Tc,float Ta,float Tss, float RH, float V,float Ea,float P,float *param,float *sitev,
		                  float &d, float &Z0c, float &Vz, float &Rkinc, float &Rkinbc, float &Tac, float &Fs, float &Ess, float &Esc,  // Output variables				 
		                  float &QHc, float &QEc, float &Ec, float &QHs, float &QEs, float &Es, float &QH, float &QE, float &E
					  ); 
		//     Calculates amount of snow intercepted by the canopy
		 __host__ __device__  void  INTERCEPT (float Ta, float LAI,float P, float Wc, float dt, float Inmax, float Uc, float Cc,
		     	         // Output variables
						 float &ieff, float &Ur, float &intc
					   );
		//     Calculates the heat advected to the snowpack due to rain
		 __host__ __device__  float QPF(float PR, float TA, float TO, float PS, float RHOW, float HF, float CW, float CS);
		//      Routine to correct energy and mass fluxes when numerical overshoots dictate that W was changed in the calling routine - either because W went negative
		//      or due to the liquid fraction being held constant.
		 __host__ __device__  void  PREHELP( float W1,float W,float DT,float &FM,float FM1,float fac,float PS,float PRAIN,float &E,float RHOW,float HF,float &Q,float &QM,float &MR, float &QE, float HSF);
		//     Computes the exponential integral function for the given value      
		 __host__ __device__ float EXPINT(float LAI);
		//     Calculates the vapour pressure at a specified temperature over water or ice depending upon temperature.  Temperature is celsius here.
		 __host__ __device__ float svp(float T);
		//     Calculates the vapour pressure at a specified temperature over water using polynomial from Lowe (1977).
		 __host__ __device__ float svpw(float T);
		//     Calculates the vapour pressure at a specified temperature over ice using polynomial from Lowe (1977).
		 __host__ __device__ float svpi(float T);
		//     Estimates reflection and scattering coefficient
		 __host__ __device__ float Tau1(float Rho, float G, float h, float COSZEN, float kk);
		 __host__ __device__ float Tau2(float Rho, float G, float h, float kk, float EXPI);

		//#######  snowdgtv.cpp
		//ueb 'driving' function
		//yjs Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow 
		//yjs drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
		//yjs snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes 
		//yjs back. //
		 __host__ __device__  void  SNOWUEB2(); 
		    //float dt, float *input, float *sitev, float *statev, float* tsprevday, float* taveprevday, int &nstepday, float *param, int *iflag,
			//float &cump, float &cumes, float &cumEc, float &cumMr, float &cumGM, float &cumEg, float *outv, float *mtime, float atff, float cf, float *OutArr);
		//************************  Predictor CORRECTOR ***********************************
		 __host__ __device__  void  PREDICORRc ( float &Us, float &Ws, float &Wc, float Alb, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float Qli,
						    float atff, float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime, 
		          	        // Following variables are output
		                   float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qmc, float &Mc, float &FMc,float &intc, float &Inmax, float &ieff, float &Ur, float &Cf, 
						float &Taufb, float &Taufd, float &Qsib, float &Qsid, float &Taub,float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &Tac,
						// Just for testing
					       float &QHs, float &QEs,float &Es, float &QPs, float &MR, float &QMs, float &Q,float &FM, float &TSURFs, float &tave, float &Qnet, float &refDepth, float &totalRefDepth, float &smelt,
						   float &gsurf, float &Qlis
					   );
		//				CALCULATE CANOPY MASS FLUX AT ANY INSTANT
		 __host__ __device__  void  QFMc(float Us, float Ws, float Wc, float A, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float Qli,
				      float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime,
					   // Following variables are output
					     float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qps, float &Qmc, float &Mc, float &FMc, float &intc,float &Inmax, float &ieff, float &Ur, float &Cf, float &Taufb,
						 float &Taufd, float &Qsib, float &Qsid, float &Taub, float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &TSURFs, float &Tac,
						 // Just for testing
						 float &tave, float &qnet, float &QHs, float &QEs, float &Es, float &MR, float &QMs, float &Q, float &FM, float &refDepth, float &totalRefDepth, float &Smelt, float &smeltc
				 );
		//  COMPUTE THE SURFACE TEMPERATURE OF SNOW 
		 __host__ __device__ float SRFTMPSC(float &Tssk, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float cf, float Qli, float cosZen,
						   float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Qpsin,float Inmax,float Rkinc,float Rkinsc,float Vz,
						     float &Tc,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &Tave,float &refDepth,float &smelt,float &smeltC
					    );
		// COMPUTE THE CANOPY TEMPERATURE
		 __host__ __device__ float CanTemp(float& Tck, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float cf, float Qli,
		               float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Inmax,float Rkinsc,float Vz,
					     float& Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &tave,float &refDepth,float &smeltC
					);
		//      FUNCTION TO EVALUATE THE SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR SURFACE TEMPERATURE                      DGT and C Luce 4/23/97
		 __host__ __device__ float SURFEBSC(float Tssk, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Fs, float Cf, float Qli, float Qsi,
		                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpc,float Qps,float &Inmax,
						    float &Rkinc,float &Rkinsc,float &Vz,float Tck,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &TSURFs,float &Tave,float &refDepth
						);// Heat and vapor conductance for neutral
		//*************  FUNCTION TO EVALUATE THE CANPPY SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR CANOPY TEMPERATURE 
		 __host__ __device__ float SURFEBC(float Tck, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Fs, float Cf, float Qli, float Qsi,
		                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float QPc,float &QPs,float &Inmax,
						     float &Rkinc,float &Rkinsc,float &Vz,float Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &FkappaS,float RHO, float &TherC,float &TSURFs,float &tave,float &refDepth
					  );                             // Reduced parametre later                                    // Heat and vapor conductance for neutral
		//                           ********  QcEst() *******
		 __host__ __device__ float QcEst(float Ws, float P, float Tssk, float &Tck, float V, float &Zm, float &d, float &Z0c, float &Rimax, float &Rcastar, float Cf, float Fs, float Qli,
		               float &Hcan, float Vz,float Ta, float Rh, float RKINsc, float QPs, float To, float Ps,float Qsi, float atff, float cosZen, float APr,float Tak,
					     float &EA, float &A, float &Ac, float Wc, float Inmax, float Qnetob, int Iradfl, float *param, float *sitev
					);
		//     Calculates the melt rate and melt outflow
		 __host__ __device__ float FMELT(float UB, float RHOW, float W, float HF, float LC, float RID, float KS, float PRAIN);
		// Linear simple function to calculate the gradient of Q and T
		 __host__ __device__  void  Grad(float qc1, float qc2,float t1,float t2, float &a,float &b);
		//     function to calculate the refreezing depth
		 __host__ __device__ float  refDep(float flans, float a, float b, float hf, float  rhom, float dt, float x1);
		 __host__ __device__ float TAVG(float UB, float W, float RHOW, float CS, float TO, float RHOG, float DE, float CG, float HF);
		//yjs  Calculate the daily average value //yjs  n number of records, a minimum value -100 or some
		 __host__ __device__ float daily_ave(float* backup, int n, float a);
		//	Function to get the LanE which is the thermal conductivity by ze
		 __host__ __device__ float  LanE(float LanS, float LanG, float Zs, float rho, float rhog, float cs, float cg, float r, float &ze, float w1day);
		//     Function to calculate Albedo  BATS Albedo Model (Dickinson et. al P.21)
		 __host__ __device__ float ALBEDO(float tausn, float coszen, float d, float aep, float abg, float avo, float airo);
		//     Function to calculate Dimensionless age of snow for use in BATS Albedo Model (Dickinson et. al P.21)
		 __host__ __device__  void  AGESN(float &tausn, float dt, float Ps, float tsurf, float tk, float dNewS);
		//     Partitioning of precipitation into rain and snow      
		 __host__ __device__ float PARTSNOW(float P, float TA, float TR, float TS);

};
//read multiple slubs (y,x arrays) along the time dim (for interval between tstart --- tend) 
__host__ __device__ int readNC_yxSlub(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, int tEnd,float*** &pvar_in, int &nrecords, int &numNc, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read multiple slubs (y,x arrays) along the time dim for data with t, y, x config---time as slowely varying array; and pvar_in already allocated
__host__ __device__ int readNC_yxSlub_givenT(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME, int &tStart, int tEnd,
	int Nydim, int Nxdim, float*** pvar_in, int &nrecords, int &numNc, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read whole nc array for a given variable to contiguous 1d array; the y,x,time order is irrelevant inside the function, but it is assumed known by caller,
__host__ __device__ int readNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &Ntdim, float* &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read whole nc array for a given variable to contiguous 3d array; the y,x,time order is irrelevant inside the function, but it is assumed known by caller,
__host__ __device__ int read3DNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &Ntdim, float*** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read nc to contiguous 1d array;  y,x,time coordinate names are same as index names
__host__ __device__ int readNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* tcorvar, float* ycorvar, float* xcorvar, size_t &Ntdim, size_t &Nydim, size_t &Nxdim, float* &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read 3d nc to contiguous array;  y,x,time coordinate names are same as index names
__host__ __device__ int read3DNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* tcorvar, float* ycorvar, float* xcorvar, size_t &Ntdim, size_t &Nydim, size_t &Nxdim, float*** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//writes the 1D array (TS values) at specified location in the netcdf
__host__ __device__ int WriteTSto3DNC(const char* FileName, const char* VarName, int dimOrd, int y_dim, int x_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//writes multiple 1D arrays (TS values) at specified locations in the netcdf
__host__ __device__ int WriteTSto3DNC_Block(const char* FileName, const char* VarName, int dimOrd, int *YindArr, int *XindArr, int bSize, int Nt_dim, float** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//aggregate outputs at specified zone in the netcdf
__host__ __device__ int Write_uebaggTS_toNC(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp);
//aggregate outputs at specified zone in the netcdf in parallel 
__host__ __device__ int Write_uebaggTS_toNC_par(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 2D netcdf and stores dimension variables for UEB aggregated outputs; some attributes are copied from the watershed netCDF
__host__ __device__ int create3DNC_uebAggregatedOutputs(const char* FileName, aggOutput *aggOut, int naggOut, const char* tName, const char* tUnits, const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd,
	float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, int nZones, const char * zName, float* y_inp, float* x_inp);
//creates 3D netcdf in PARALLEL mode, and stores dimension variables for UEB aggregated outputs; some attributes are copied from the watershed netCDF
__host__ __device__ int create3DNC_uebAggregatedOutputs_par(const char* FileName, aggOutput *aggOut, int naggOut, const char* tName, const char* tUnits, const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd,
	float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, int nZones, const char * zName, float* y_inp, float* x_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf and stores dimension variables for a given UEB output; called once for a given output netcdf 
//attributes are copied from the 2D (watershed) netCDF which the 3D netCDF spatial dimension follow
__host__ __device__ int create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnits, const char* tlong_name,
	const char* tcalendar, int Nt_dim, int dimOrd, float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
__host__ __device__ int Create3DNC(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* yName, const char* xName, const char* tUnits,
	const char* yUnits, const char* xUnits, int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float* fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo);
//This one creates and stores array passed to it at once
__host__ __device__ int Write3DNC(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* yName, const char* xName,
	const char* tUnits, const char* yUnits, const char* xUnits,
	int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float*** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
__host__ __device__ int read3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME,
	const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//function to read multiple blocks of single column/rod along time dimension from 3D netcdf file, for given y , x coordinate arrays
__host__ __device__ int readNC_TS_Block(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float** &pvar_in, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo, int *YindArr, int *XindArr, int bSize);
//function to read single column/rod along time dimension from 3D netcdf file, for given y , x coordinates
__host__ __device__ int readNC_TS(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* &pvar_in, float* &tcorvar, int ydim, int xdim, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo);
//function to read single column/rod variable array (no time array) for given y , x coordinates
__host__ __device__ int readNC_TS(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* &pvar_in, int ydim, int xdim, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo);
__host__ __device__ int read2DNC(const char* FILE_NAME, const char* VAR_NAME, float** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//the following is to read watershed file
__host__ __device__ int readwsncFile(const char* FILE_NAME, const char* VAR_NAME, const char* ycor_NAME,
	const char* xcor_NAME, float* &ycorvar, float* &xcorvar, int** &pvar_in, int &ydim, int &xdim, int &fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo);
__host__ __device__ int Write3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME,
	const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_out, MPI::Intracomm inpComm, MPI::Info inpInfo);

//ueb inputs 
 void  readParams(const char* inpFile, float* &parArray, const int nParams);
 void  readParams(const char* inpFile, params strParamValues);
 void  readSiteVars(const char* inpFile, sitevar svArr[], int indx);
//another overload because the code fails at run time for the above_#7.4.13
 void  readSiteVars(const char* inpFile, sitevar *&svArr);
 void  readSiteVars(const char* inpFile, float svSValue[], char* svFile[], char* svVarName[], int svType[] );
 void  readInputForcVars(const char* inputconFile, inpforcvar *frArray);
 void  readOutputControl(const char* outputconFile, pointOutput* &pOut, ncOutput* &ncOut, aggOutput* &aggOut, int &npout, int &nncout, int &naggOut);
// __host__ __device__  void  readTextData(const char* inforcFile, inptimeseries* &strinpts, int &nrecords);
 void  readTextData(const char* inforcFile, float *&tcor_var, float *&tvar_in, int &nrecords);  
 void  readTextData(const char* inforcFile, float *&tvar_in, int &nrecords);
 void  readTStextFile(const char* inforcFile, float *&tvar_in, int &nrecords);

//arrays and matrices
 __host__ __device__ float*** Create3DArray(int nt, int nr, int nc);
 __host__ __device__  void  Delete3DArray(float ***A, int nt, int nr, int nc);
//create 2D array and allocate contiguous memory block this enbales a full block read of netcdf
 __host__ __device__ float** create2DArray_Contiguous(int nr, int nc);
//delets a 2D array (frees memory allocated) contiguously
 __host__ __device__  void  delete2DArray_Contiguous(float** myMatrix);
//create 3D array and allocate contiguous memory block this enbales a full block read of netcdf
 __host__ __device__ float*** create3DArrayblock_Contiguous(int nt, int nr, int nc);      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
//delets a 3D array (frees memory) allocated contiguously
 __host__ __device__  void  delete3DArrayblock_Contiguous(float*** myMatrix); // int nr, int nc)// input: 3D array

//utility functions
//findMin findMax functions
 __host__ __device__ float findMax(float a, float b);
 __host__ __device__ float findMin(float a, float b);
 __host__ __device__ int findMax(int a, int b);
 __host__ __device__ int findMin(int a, int b);
//print array values (for inspection during debugging)
 /*__host__ __device__  void  printArrayvalues(int arrLength, float* arrValues);
 __host__ __device__  void  printArrayvalues(int yLength, int xLength, float **arrValues);*/