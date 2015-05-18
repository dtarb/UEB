#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <stdio.h>
#include "mpi.h"
#include <netcdf.h>
#include <netcdf_par.h>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
//#include <algorithm>
using namespace std;

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

/* Handle errors by printing an error message and exiting with a * non-zero status. */
#define  ERR(e) {cout<<"Error: "<< nc_strerror(e)<<endl; return 2; }

//creates 3D netcdf and stores dimension variables for a given UEB output; called once for a given output netcdf 
//attributes are copied from the 2D (watershed) netCDF which the 3D netCDF spatial dimension follow
int create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnits, const char* tlong_name,
	const char* tcalendar, int Nt_dim, int dimOrd, float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf and stores dimension variables for UEB aggregated outputs; some attributes are copied from the watershed netCDF
int create3DNC_uebAggregatedOutputs(const char* FileName, aggOutput *aggOut, int naggOut, const char* tName, const char* tUnits, const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd,
	float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, int nZones, const char * zName, float* y_inp, float* x_inp);
//creates 3D netcdf in parallel and stores dimension variables for UEB aggregated outputs; some attributes are copied from the watershed netCDF
int create3DNC_uebAggregatedOutputs_par(const char* FileName, aggOutput *aggOut, int naggOut, const char* tName, const char* tUnits, const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd,
	float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, int nZones, const char * zName, float* y_inp, float* x_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//aggregate outputs at specified zone in the netcdf
int Write_uebaggTS_toNC(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp);
int Write_uebaggTS_toNC_par(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
int Create3DNC(const char* FileName, const char* VarName, const char *varUnits,  const char* tName,  const char* yName, const char* xName, const char* tUnits,
	const char* yUnits, const char* xUnits, int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float *fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo);
//writes the 1D array (TS values) at specified location in the netcdf
int WriteTSto3DNC(const char* FileName, const char* VarName, int dimOrd, int y_dim, int x_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);  //ydim, xdim =the coordinate point the data to be written; Nt_dim =the length of the TS
//writes multiple 1D arrays (TS values) at specified locations in the netcdf
int WriteTSto3DNC_Block(const char* FileName, const char* VarName, int dimOrd, int *YindArr, int *XindArr, int bSize, int Nt_dim, float** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
int Write3DNC(const char* FileName, const char* VarName, const char *varUnits,  const char* tName,  const char* yName, const char* xName, 
	const char* tUnits, const char* yUnits, const char* xUnits, int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float*** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo); // int &t_dimid, int &y_dimid, int &x_dimid) //time, y, x
int read3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME, const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
int read2DNC(const char* FILE_NAME, const char* VAR_NAME, float** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//the following is to read watershed file
int readwsncFile(const char* FILE_NAME, const char* VAR_NAME, const char* ycor_NAME,  
	const char* xcor_NAME, float* &ycorvar, float* &xcorvar, int** &pvar_in, int &ydim, int &xdim, int &fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo);  // , std::set<int> zValues, float * z_ycor, float *z_xcor);
int Write3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME, 
	const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_out, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read 3d nc to contiguous array;  y,x,time coordinate names are same as index names
int read3DNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* tcorvar, float* ycorvar, float* xcorvar, size_t &Ntdim, size_t &Nydim, size_t &Nxdim, float*** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//function to read single column/rod along time dimension from 3D netcdf file, for given y , x coordinates
int readNC_TS(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME,const char* xcor_NAME,
	float* &pvar_in, float* &tcorvar, int ydim, int xdim, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo);    //, int xstride); //double* tcorvar, double*** pvar_in)
//function to read multiple blocks of single column/rod along time dimension from 3D netcdf file, for given y , x coordinate arrays
int readNC_TS_Block(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float** &pvar_in, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo, int *YindArr, int *XindArr, int bSize); // /*float* &tcorvar, int ydim, int xdim, */ int tIndx) 

//arrays and matrices
float*** Create3DArray(int nt, int nr, int nc);
void Delete3DArray(float ***A, int nt, int nr, int nc);
//create 2D array and allocate contiguous memory block this enbales a full block read of netcdf
float** create2DArray_Contiguous(int nr, int nc);
//delets a 2D array (frees memory allocated) contiguously
void delete2DArray_Contiguous(float** myMatrix);
//create 3D array and allocate contiguous memory block this enbales a full block read of netcdf
float*** create3DArrayblock_Contiguous(int nt, int nr, int nc);      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
//delets a 3D array (frees memory) allocated contiguously
void delete3DArrayblock_Contiguous(float*** myMatrix); // int nr, int nc)// input: 3D array

//ueb inputs
void readParams(const char* inpFile, float* &parArray, const int nParams);
void readParams(const char* inpFile, params strParamValues);
void readSiteVars(const char* inpFile, sitevar svArr[], int indx);
//another overload because the code fails at run time for the above_#7.4.13
void readSiteVars(const char* inpFile, sitevar *&svArr);
void readSiteVars(const char* inpFile, float svSValue[], char* svFile[], char* svVarName[], int svType[] );
void readInputForcVars(const char* inputconFile, inpforcvar *frArray);
void readOutputControl(const char* outputconFile, const char* aggoutputconFile, pointOutput* &pOut, ncOutput* &ncOut, aggOutput* &aggOut, int &npout, int &nncout, int &naggOut);
//void readTextData(const char* inforcFile, inptimeseries* &strinpts, int &nrecords);
void readTextData(const char* inforcFile, float *&tcor_var, float *&tvar_in, int &nrecords);  
void readTextData(const char* inforcFile, float *&tvar_in, int &nrecords);

//*******************************************************************************

//js  Constant extern float set 
extern float T_0; 				// Temperature of freezing (0 C)
extern float T_k ;  // 273.15;			// Temperature to convert C to K (273.15)
extern float SB_c ;  // Stefan boltzman constant (2.041334e-7 KJ/m^2-hr-K^4) #corrected 12.23.14
extern float H_f ;  // 333.5;			// Heat of fusion (333.5 KJ;  // kg)
extern float Hne_u ;  // 2834.0;		// Heat of Vaporization (Ice to Vapor, 2834 KJ;  // kg)
extern float C_w ;  // 4.18;			// Water Heat Capacity (4.18 KJ;  // kg;  // C)
extern float C_s ;  // 2.09;			// Ice heat capacity (2.09 KJ;  // kg;  // C)
extern float C_p ;  // 1.005;			// Air Heat Capacity (1.005 KJ;  // kg;  // K)
extern float Ra_g ;  // 287.0;			// Ideal Gas constant for dry air (287 J;  // kg;  // K)
extern float K_vc  ;  // 0.4;			// Von Karmans constant (0.4)
extern float Hs_f ;  // 3600.0;		// Factor to convert ;  // s into ;  // hr (3600)
extern float Rho_i ;  // 917.0;			// Density of Ice (917 kg;  // m^3)
extern float Rho_w ;  // 1000.0;		// Density of Water (1000 kg;  // m^3)
extern float Gra_v    ;  // 9.81;			// Gravitational acceleration (9.81 m;  // s^2)
extern float W1da_y ;  // 0.261799;		// Daily frequency (2pi;  // 24 hr 0.261799 radians;  // hr) 
extern float  Io ;  // 4914.0;            //  Solar constant  Kj/m^2/hr
//pi copied from snowdxv.f90
extern float P_i ;  // 3.141592653589793238462643383279502884197169399375105820974944592308;		// Pi

//data for pred-corr
extern float wtol ;  // 0.025;
extern float utol ;  // 2000.0;
//from TURBFLUX()
extern float tol ;  // 0.001;
extern int  nitermax ;  // 20;  
extern  int  ncitermax; //21

// #***^*-flag to write warnings,...etc ###***TBC 9.20.13
//Global variables ###_TBC 5.6.13
extern int snowdgtvariteflag ;  // 0; print debug information
extern int snowdgtvariteflag2 ;  // more detailed print in Surfebsc
extern int snowdgtvariteflag3 ;  // print in predictor corrector function
extern int snowdgt_outflag ;        //output print
extern int radwarnflag;
extern int inpDailyorSubdaily;       // 0: values given at each (sub-daily time steps); 1: daily values
//track current grid cell
extern int uebCellX;
extern int uebCellY;

//common 
extern float Tsk_save ;  // 273.16,
extern float Tssk_old ;  // 273.16,  
extern float Tsavek_old ;  // 273.16, 
extern float Tsavek_ave ;  // 273.16, 
extern float Tssk_ave ;  // 273.16/*added 6.7.13*/;

//findMin findMax functions
float findMax(float a, float b);
float findMin (float a, float b);
int findMax(int a, int b);
int findMin (int a, int b);

//###check operator precedencde issue in albedo()and other functions _6.5.13

//snowdv.cpp
void RUNUEB(float* RegArray[13], float statesiteValues[], float paramValues[], float** &OutVarValues, int modelStartDate[],
			     double modelStartHour, int modelEndDate[], double modelEndHour, double modelDT, double UTCOffset);

//snowxv.cpp
//********UPDATEtime ()  Update time for each time step
void UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT);
//    function to return number of days in February checking for leap years
int lyear(int year);
//    to get the atmospheric transmissivity using the Bristow and Campbell  (1984) approach
void atf(float &atff,float trange,int month, float *dtbar, float a, float c);
// To get hourly radiation index
void hyri(int YEAR, int MONTH, int DAY, float HOUR, float DT, float SLOPE, float AZI, float LAT, float &HRI, float &COSZEN);
  //To convert the real date to julian date
int julian(int yy, int mm, int dd);
//    Computes the incoming longwave radiation using satterlund Formula
void cloud(float as, float bs, float atff, float &cf);
//???? long wave radiation from temperatrue and other weather variables??
//TBC_6.5.13
void qlif(float TA, float RH, float TK, float SBC, float &Ema, float &Eacl, float cf, float &qliff );
//these were copied from functions.f90
//COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  
double julian ( int I,int M, int K,double H);
//COMPUTES CALENDAR DATE AND time, GIVEN JULIAN DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
void calendardate (double TJD,int &I,int &M,int &K, double &H);

//##### canopy.cpp       
//     Partition the incoming solar radiation into direct and diffuse components
void PSOLRAD( float Qsi, float atff, float *param, float cf,     
                float &Taufb, float &Taufd, float &Qsib, float &Qsid);     // Output variables
//     Estimates the direct and diffuse solar radiation fractions transmitted through the canopy      
void TRANSRADCAN (float COSZEN, float *sitev, float *param,       
     	            float &Betab, float &Betad, float &Taub, float &Taud);                     // Output variables: Transmission and Reflection fractions
//      Computes the net canopy and sub-canopy solar radiation
//      considering the multiple scatterings of radiation by the canopy and multile reflections 
//      between the canopy and the snow surface
void NETSOLRAD(float Ta,float A,float Betab,float Betad,float Wc,float Taub,float Taud, 
			     float Qsib,float Qsid,float *param, float	Fs,float &Qsns,float &Qsnc); //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
//     Computes the net canopy and beneath canopy longwave radiation
void NETLONGRAD(float RH,float Ta,float Tss,float Tc,float Tk,float Fs,float EmC,float EmS,float SBC,float cf,float *sitev,float Qli,float *param,
				    float &Qlis, float &Qlns,float &Qlnc );
//****************  BULK AERODYNAMIC AND CANOPY BOUNDARY LAYER RESISTANCES  *******************
//      Calculates resistances for the above and beneath canopy turbulent heat fluxes 
void AeroRes(float P, float Wc, float V, float Ta, float Tss, float Tc, float Fs, float *param, float *sitev, float Tk,
			    float &d, float &Z0c, float &Vz, float &Rc, float &Ra, float &Rbc, float &Rl, float &RKINc, float &RKINa, float &RKINbc, float &RKINl);         // Output variables
//     Calculates wind at the 2 m above the surface and at the sink using the input of measured wind at 2 m above the canopy
void WINDTRANab(float V, float Zm, float *param, float *sitev,float &Vstar,float &Vh,float &d,float &Z0c,float &Vz,float &Vc);               // Out put wind speed within canopy at height 2 m
//     Calculates the turbulent heat fluxes (sensible and latent heat fluxes) and condensation/sublimation.
void TURBFLUX (float Ws, float Wc, float A, float Tk,float Tc,float Ta,float Tss, float RH, float V,float Ea,float P,float *param,float *sitev,
                  float &d, float &Z0c, float &Vz, float &Rkinc, float &Rkinbc, float &Tac, float &Fs, float &Ess, float &Esc,  // Output variables				 
                  float &QHc, float &QEc, float &Ec, float &QHs, float &QEs, float &Es, float &QH, float &QE, float &E
			  ); 
//     Calculates amount of snow intercepted by the canopy
void INTERCEPT (float Ta, float LAI,float P, float Wc, float dt, float Inmax, float Uc, float Cc,
     	         // Output variables
				 float &ieff, float &Ur, float &intc
			   );
//     Calculates the heat advected to the snowpack due to rain
float QPF(float PR, float TA, float TO, float PS, float RHOW, float HF, float CW, float CS);
//      Routine to correct energy and mass fluxes when numerical overshoots dictate that W was changed in the calling routine - either because W went negative
//      or due to the liquid fraction being held constant.
void PREHELP( float W1,float W,float DT,float &FM,float FM1,float fac,float PS,float PRAIN,float &E,float RHOW,float HF,float &Q,float &QM,float &MR, float &QE, float HSF);
//     Computes the exponential integral function for the given value      
float EXPINT (float LAI);
//     Calculates the vapour pressure at a specified temperature over water or ice depending upon temperature.  Temperature is celsius here.
float svp(float T);
//     Calculates the vapour pressure at a specified temperature over water using polynomial from Lowe (1977).
float svpw(float T);
//     Calculates the vapour pressure at a specified temperature over ice using polynomial from Lowe (1977).
float svpi(float T);
//     Estimates reflection and scattering coefficient
float Tau1(float Rho, float G, float h, float COSZEN, float kk);
float Tau2(float Rho, float G, float h, float kk, float EXPI);

//#######  snowdgtv.cpp
//ueb 'driving' function
//yjs Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow 
//yjs drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
//yjs snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes 
//yjs back. //
void SNOWUEB2(float dt, float *input, float *sitev, float *statev, float *tsprevday, float *taveprevday, int &nstepday, float *param, int *iflag,
			            float &cump, float &cumes, float &cumEc, float &cumMr, float &cumGM,float &cumEg,float *outv, float *mtime ,float atff, float cf, float *OutArr);
//************************  Predictor CORRECTOR ***********************************
void PREDICORRc ( float &Us, float &Ws, float &Wc, float Alb, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float Qli,
				    float atff, float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime, 
          	        // Following variables are output
                   float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qmc, float &Mc, float &FMc,float &intc, float &Inmax, float &ieff, float &Ur, float &Cf, 
				float &Taufb, float &Taufd, float &Qsib, float &Qsid, float &Taub,float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &Tac,
				// Just for testing
			       float &QHs, float &QEs,float &Es, float &QPs, float &MR, float &QMs, float &Q,float &FM, float &TSURFs, float &tave, float &Qnet, float &refDepth, float &totalRefDepth, float &smelt,
				   float &gsurf, float &Qlis
			   );
//				CALCULATE CANOPY MASS FLUX AT ANY INSTANT
void QFMc(float Us, float Ws, float Wc, float A, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float Qli,
		      float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime,
			   // Following variables are output
			     float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qps, float &Qmc, float &Mc, float &FMc, float &intc,float &Inmax, float &ieff, float &Ur, float &Cf, float &Taufb,
				 float &Taufd, float &Qsib, float &Qsid, float &Taub, float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &TSURFs, float &Tac,
				 // Just for testing
				 float &tave, float &qnet, float &QHs, float &QEs, float &Es, float &MR, float &QMs, float &Q, float &FM, float &refDepth, float &totalRefDepth, float &Smelt, float &smeltc
		 );
//  COMPUTE THE SURFACE TEMPERATURE OF SNOW 
float SRFTMPSC (float &Tssk,float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Qsi,float atff,float cf,float Qli,float cosZen,
				   float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Qpsin,float Inmax,float Rkinc,float Rkinsc,float Vz,
				     float &Tc,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &Tave,float &refDepth,float &smelt,float &smeltC
			    );
// COMPUTE THE CANOPY TEMPERATURE
float CanTemp(float& Tck,float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Qsi,float atff,float cf,float Qli,
               float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Inmax,float Rkinsc,float Vz,
			     float& Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &tave,float &refDepth,float &smeltC
			);
//      FUNCTION TO EVALUATE THE SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR SURFACE TEMPERATURE                      DGT and C Luce 4/23/97
float SURFEBSC(float Tssk, float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Fs,float Cf,float Qli,float Qsi,
                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpc,float Qps,float &Inmax,
				    float &Rkinc,float &Rkinsc,float &Vz,float Tck,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &TSURFs,float &Tave,float &refDepth
				);// Heat and vapor conductance for neutral
//*************  FUNCTION TO EVALUATE THE CANPPY SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR CANOPY TEMPERATURE 
float SURFEBC(float Tck,float Us,float Ws,float Wc,float A,float dt,float P,float Pr,float Ps,float Ta,float V,float RH,float Fs,float Cf,float Qli,float Qsi,
                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float QPc,float &QPs,float &Inmax,
				     float &Rkinc,float &Rkinsc,float &Vz,float Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &FkappaS,float RHO, float &TherC,float &TSURFs,float &tave,float &refDepth
			  );                             // Reduced parametre later                                    // Heat and vapor conductance for neutral
//                           ********  QcEst() *******
float QcEst(float Ws, float P, float Tssk,float &Tck, float V, float &Zm, float &d, float &Z0c, float &Rimax, float &Rcastar, float Cf, float Fs, float Qli, 
               float &Hcan, float Vz,float Ta, float Rh, float RKINsc, float QPs, float To, float Ps,float Qsi, float atff, float cosZen, float APr,float Tak,
			     float &EA, float &A, float &Ac, float Wc, float Inmax, float Qnetob, int Iradfl, float *param, float *sitev
			);
//     Calculates the melt rate and melt outflow
float FMELT(float UB,float RHOW,float W,float HF,float LC,float RID,float KS,float PRAIN);
// Linear simple function to calculate the gradient of Q and T
void Grad(float qc1, float qc2,float t1,float t2, float &a,float &b);
//     function to calculate the refreezing depth
float  refDep(float flans,float a,float b,float hf,float  rhom,float dt,float x1 );
float TAVG( float UB, float W, float RHOW, float CS, float TO, float RHOG, float DE, float CG, float HF);
//yjs  Calculate the daily average value //yjs  n number of records, a minimum value -100 or some
float daily_ave(float *backup, int n, float a);
//	Function to get the LanE which is the thermal conductivity by ze
float  LanE( float LanS, float LanG, float Zs,float rho, float rhog,float cs, float cg, float r,float &ze,float w1day);
//     Function to calculate Albedo  BATS Albedo Model (Dickinson et. al P.21)
float ALBEDO(float tausn, float coszen, float d, float aep,float abg, float avo,float airo);
//     Function to calculate Dimensionless age of snow for use in BATS Albedo Model (Dickinson et. al P.21)
void AGESN(float &tausn, float dt, float Ps, float tsurf, float tk, float dNewS);
//     Partitioning of precipitation into rain and snow      
float PARTSNOW(float P, float TA, float TR, float TS);
