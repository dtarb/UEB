//copied from File  snowx.f  
//  Subroutines and function subprograms for the Utah Energy Balance
//  Snow Accumulation and Melt Model.
//  David G. Tarboton, Utah Water Research Laboratory, Utah State University
//  
//  Last Change 9/9/12 to accommodate glacier melt.
//
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
//  If you wish to use or incorporate this program (or parts of it) into 
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
//common declarations
#include "gpuuebpgdecls.h"
//********UPDATEtime ()  Update time for each time step
 __host__ __device__  void    uebCell::UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT)
{
	int DM;				 // 30/03/2004 ITB 
						  // 30/03/2004 ITB  
	//real  hour, dt  // DGT Dec 10, 2004.  Fixing ITB errors 

	int DMON[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	HOUR = HOUR + DT;
	DM = DMON[MONTH-1];
	//  check for leap years 
	if(MONTH == 2)
		DM = lyear (YEAR);
	while(HOUR >= 24.0)
	{
		HOUR = HOUR - 24.0;
		DAY++;
	}

	while ( DAY > DM)
	{
		DAY = DAY - DM;
		MONTH++;
		if(MONTH>12){
			MONTH = 1;
			YEAR++;
		}
		//modified from the original by separating the above two lines in the if (month>12)
		//#_6.27.13
			DM = DMON[MONTH-1];
			if(MONTH  ==  2)
				DM= lyear(YEAR);
		//}
	}
	return;
}

// ************************** lyear () ***************************
//    function to return number of days in February checking for leap years
 __host__ __device__ int uebCell::lyear(int year)
{
	int lyear; // Leap years are every 4 years 
			   // - except for years that are multiples of centuries (e.g. 1800, 1900)
			   // - except again that when the century is divisible by 4 (e.g. 1600, 2000)
	if( (year % 4  > 0) || ((year % 100  == 0) && (year %400 != 0)))
		lyear=28;
	else
		lyear=29;
	return lyear;
}
//**************************** atf () ****************************
//    to get the atmospheric transmissivity using the Bristow and Campbell  (1984) approach
 __host__ __device__  void    uebCell::atf(float &atff,float trange,int month, float *dtbar, float a, float c)
{
	//DIMENSION dtbar(12)
	float b = 0.036* expf(-0.154*dtbar[month-1]);
	atff = a*(1-expf(-b * powf(trange,c)));
	//     write(6,*)trange,month,a,c,dtbar(month),atf
	return;
}

//************************** hourlyRI () To get hourly radiation index
 __host__ __device__  void    uebCell::hyri(int YEAR, int MONTH, int DAY, float HOUR, float DT, float SLOPE, float AZI, float LAT, float &HRI, float &COSZEN)
{
	float  LP,LAT1;   
	// lp= latitude of equivalent plane in radians
	// lat1 = latitude in radians
	// lat = latitude in degrees
	// a number that speaks for itself - every kissable digit
	float PI=3.141592653589793238462643383279502884197169399375105820974944592308;
	float CRAD = PI/180.0;
	//  crad = degree to radian conversion factor
	//    CONVERT timeS TO RADIANS FROM NOON
	float T = (HOUR-12.0)*PI/12.0;
	float DELT1= DT*PI/12.0;
	//    CONVERT angles TO RADIANS
	float SLOPE1=SLOPE*CRAD;
	float AZI1=AZI*CRAD;
	LAT1=LAT*CRAD;
	float   FJULIAN = (float) julian(YEAR,MONTH,DAY);
	float D = CRAD*23.5* sin((FJULIAN-82.0)*0.017214206321); 
	// 0.017214206321 is 2 pi / 365  
	// D is solar declination
	LP=asin(sin(SLOPE1)*cos(AZI1)*cos(LAT1) + cos(SLOPE1)*sin(LAT1));
	// LP is latitude of equivalent plane
	//     TD=ACOS(-TAN(LAT1)*TAN(D))  This formula abandoned 1/8/04 
	//     to make the code work for polar conditions
	// TD is half day length, i.e. the time from noon to sunset.  Sunrise is at -TD
	float tanprod = tan(LAT1)* tan(D);
	float td;
	if(tanprod > 1.0)
		td=PI;  // This is the condition for perpetual light
	else if(tanprod < -1.)
		td=0;   // The condition for perpetual night
	else
		td=acos(-tanprod);  // The condition where there is a sunrise and set

	//  Equivalent longitude offset.  Modified on 1/8/04
	//  so that it correctly accounts for shift in longitude if equivalent 
	//  plane slope goes over a pole.  Achieved using atan2.
	//     DDT=ATAN(sin(AZI1)*sin(SLOPE1)/(cos(SLOPE1)*cos(LAT1)
	//    *    -cos(AZI1)*sin(SLOPE1)*sin(LAT1)))
	float ddt= atan2(sin(AZI1)*sin(SLOPE1), (cos(SLOPE1)*cos(LAT1) - cos(AZI1)*sin(SLOPE1)*sin(LAT1)));  

	//  Now similar logic as before needs to be repeated for equivalent plane
	//  but with times reflecting
	float tpeqp = tan(LP)*tan(D);
	// Keep track of beginning and end of exposure of equiv plane to sunlight
	float tpbeg, tpend;
	if(tpeqp > 1.0)
	{
		tpbeg = -PI;   // perpetual light
		tpend= PI;
	}
	else if (tpeqp < -1.)
	{
		tpbeg=0.0;  // perpetual dark
		tpend=0.0 ;
	}
	else
	{
		tpbeg = -acos(-tpeqp) - ddt;
		tpend = acos(-tpeqp) - ddt;
	}


	//  Start and end times for integration of radiation exposure
	//  need to account for both horizon, slope and time step
	float T1, T2;
	T1 = findMax(T,tpbeg);
	T1 = findMax(T1,-td);
	T2 = findMin(T+DELT1,td);
	T2 = findMin(T2,tpend);
	//     write(6,*)t1,t2
	if(T2 <= T1) 
		HRI=0.0;
	else
		HRI = (sin(D)*sin(LP)*(T2-T1) + cos(D)*cos(LP)*(sin(T2+ddt) - sin(T1+ddt)) ) / (cos(SLOPE1)*DELT1);
	//  In the above the divide by cos slope normalizes illumination to per unit horizontal area

	//  There is a special case if tpbeg is less than -pi that occurs in polar regions
	//  where a poleward facing slope may be illuminated at night more than the day.
	//  Add this in
	if(tpbeg < -PI)
	{
		T1 = findMax(T, 2*PI-tpbeg);
		T1 = findMax(T1,-td);
		T2 = findMin(T+DELT1,td);
		if(T2 > T1)
		{
			HRI = HRI + (sin(D)*sin(LP)*(T2-T1) + cos(D)*cos(LP)*(sin(T2+ddt) - sin(T1+ddt))) / (cos(SLOPE1)*DELT1);
		}
	}
	// for the purposes of calculating albedo we need a cosine of the
	// illumination angle.  This does not have slope correction so back
	// this out again.  This is an average over the time step
	COSZEN = HRI*cos(SLOPE1);
	//     write(6,*)hri,coszen
	return;
}
//***************************** JULIAN () ****************************
//             To convert the real date to julian date
// YJS The Julian are change to a new version to take the Leap Yean into consideration
//    in the old version, there are 365 days each year.
//     FUNCTION JULIAN(MONTH,DAY)
__host__ __device__ int uebCell::julian(int yy, int mm, int dd)
{
	int julian;
	int mmstrt[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
	int jday = mmstrt[mm-1] + dd;
	int ileap = yy - ((int)(yy/4)) * 4 ;
	if((ileap == 0) && (mm >=3))
		jday = jday + 1;
	julian = jday;
	return julian; 
}
//******************** For cloudiness fraction cf *********************
//    Computes the incoming longwave radiation using satterlund Formula
//    Modified 10/13/94 to account for cloudiness.  Emissivity of cloud cover fraction is assumed to be 1.
 __host__ __device__  void    uebCell::cloud(float as, float bs, float atff, float &cf)
{
	//as     = param(28)                        // Fraction of extraterrestaial radiation on cloudy day,Shuttleworth (1993)  
	//bs     = param(29)                      // (as+bs):Fraction of extraterrestaial radiation on clear day, Shuttleworth (1993) 
	if (atff >= (as+bs))
	  cf=0;                                                             // Cloudiness fraction
	else if(atff <= as) 
	  cf=1;
	else 
	 cf = 1.0 - (atff - as)/bs;
	return;
}

//************************************ QLIF ()*********************************
//???? long wave radiation from temperatrue and other weather variables??
//TBC_6.5.13
 __host__ __device__  void    uebCell::qlif(float TA, float RH, float TK, float SBC, float &Ema, float &Eacl, float cf, float &qliff )
{
	float TAK  = TA + TK;
	float EA   = RH * svpw(TA);
	//******************************************************  old option
	//    
	Eacl   =  1.08 * (1.0 - expf(-1*powf(EA/100.0, TAK/2016.0)));   // Clear sky emissivity
	Ema   =  (cf + (1.0 - cf)*Eacl);                              // Emissivity for cloudy sky
	qliff =  Ema * SBC * powf(TAK, 4.0);                                 // Incoming longwave 
	return;
}
//The following were copied from functions.f90
//# 6.8.13

//THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT time VALUE
//CAN BE IN ANY UT-LIKE time SCALE (UTC, UT1, TT, ETC.) - OUTPUT. //JULIAN DATE WILL HAVE SAME BASIS.  
//ALGORITHM BY FLIEGEL AND //VAN FLANDERN. //SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
//I = YEAR (IN) //M = MONTH NUMBER (IN) //K = DAY OF MONTH (IN) //H = UT HOURS (IN) //TJD = JULIAN DATE (OUT)
 __host__ __device__ double uebCell::julian(int I, int M, int K, double H)
{
	double TJD,JD;
	//JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
	JD = K-32075 + 1461*(I+4800 + (M-14)/12) / 4 + 367*(M-2-(M-14)/12*12)/12-3*((I+4900+(M-14)/12)/100)/4;
	TJD = JD - 0.5 + H/24.0;
	//##%^_TBC 6.8.13 //powf(10,0) in place of D0
	return TJD;
}

//THIS SUBROUTINE COMPUTES CALENDAR DATE AND time, GIVEN JULIAN DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
//(UTC, UT1, TT, ETC.) - OUTPUT time VALUE WILL HAVE SAME BASIS. OUTPUT CALENDAR DATE WILL BE GREGORIAN.  
//ALGORITHM BY FLIEGEL AND VAN FLANDERN. //SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
//TJD = JULIAN DATE (IN) //I = YEAR (OUT) //M = MONTH NUMBER (OUT) //K = DAY OF MONTH (OUT) //H = UT HOURS (OUT)
 __host__ __device__  void    uebCell::calendardate (double TJD,int &I,int &M,int &K, double &H)
{
	double DJD, JD;
	int L, N;
	DJD = TJD + 0.5; 
    JD = DJD;
	H = fmod(DJD,1.0)*24;    // 24.D0
	//JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
	L = JD + 68569;
	N = 4*L/146097;
	L = L - (146097*N+3)/4;
	//I=YEAR, M=MONTH, K=DAY
	I = 4000*(L+1)/1461001;
	L = L - 1461*I/4 + 31;
	M = 80*L/2447;
	K = L - 2447*M/80;
	L = M / 11;
	M = M + 2 - 12*L;
	I = 100*(N-49) + I + L;
	return;
}


