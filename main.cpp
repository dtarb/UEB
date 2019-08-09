
#include "uebpgdecls.h"
//#include <queue>
#pragma warning(disable : 4996)
using namespace std;

int main(int argc, char* argv[])
{
	float timeControl = 0.0, timeWS = 0.0, timeSitestate = 0.0, timeTSArrays = 0.0, timeParam = 0.0, timeParamSiteInptcontrol = 0.0, timeModelRun = 0.0;
	float** outvarArray; //= new float*[70]; //[70];
	float*** aggoutvarArray;
	float ***ncoutArray;
	char conFile[256], paramFile1[256], sitevarFile1[256], inputconFile1[256], outputconFile1[256], watershedFile1[256], aggoutputconFile1[256], aggoutputFile1[256];
	char wsvarName1[256], wsycorName1[256], wsxcorName1[256];
	int **wsArray = NULL;
	int dimlen1 = 0, dimlen2 = 0, totalgrid = 0;
	int wsfillVal = -9999;
	//int  *wsmissVal = new int;
	//*wsmissVal = -9999;
	float *parvalArray = NULL;
	float SiteState[32];
	float *wsxcorArray = NULL, *wsycorArray = NULL;
	//params ParamVAlues;
	sitevar *strsvArray = new sitevar[32];
	char * svFile[32];
	char * svVarName[32];
	for (int i = 0; i < 32; i++){
		svFile[i] = new char[256];
		svVarName[i] = new char[256];
	}
	int svType[32];
	inpforcvar strinpforcArray[13];
	//outputs
	pointOutput *pOut = NULL;
	ncOutput *ncOut = NULL;
	aggOutput *aggOut = NULL;
	int npout = 0, nncout = 0, naggout = 0, nZones = 0;
	const char * zName = "Outletlocations"; //12.24.14 watershed zonning for aggregation--	
	float *z_ycor = NULL;
	float *z_xcor = NULL;
	int zoneid = 0;
	int *ZonesArr = NULL;
	//inptimeseries *strintsArray[11];
	float *tcorvar[13], *tsvarArray[13], *tsvarArrayTemp[5]; // [xstride]; //assuming max nc files for a variable =5
	int ncTotaltimestep = 0;                                         //6.24.14
	int ntimesteps[5];
	int tinitTime = 0;
	int npar = 32;
	/*int NUMtimeSTEP,NREFYR,NREFMO,NREFDAY,NumOP;*/
	int ModelStartDate[3], ModelEndDate[3]; //check this 
	double ModelStartHour, ModelEndHour, ModelDt, ModelUTCOffset;
	double modelSpan;
	int numTimeStep;
	int numOut = 70;
	char headerLine[256];
	int retvalue = 0;
	int numgrid = 0;
	
	const char* tNameout = "time";
	int outtSteps = 0;
	int outtStride = 1, outyStep = 1, outxStep = 1;
	float* t_out;
	float out_fillVal = -9999.0;
	int outDimord = 0, aggoutDimord = 1;
	int *yIndxArr = NULL, *xIndxArr = NULL;
	const char* tlong_name = "time";
	const char* tcalendar = "standard";
	char* uebVars[70] = { "Year", "Month", "Day", "dHour", "atff", "HRI", "Eacl", "Ema", "conZen", "Ta", "P", "V", "RH", "Qsi", "Qli", "Qnet",
		"Us", "SWE", "tausn", "Pr", "Ps", "Alb", "QHs", "QEs", "Es", "SWIT", "QMs", "Q", "FM", "Tave", "TSURFs", "cump", "cumes",
		"cumMr", "NetRads", "smelt", "refDepth", "totalRefDepth", "cf", "Taufb", "Taufd", "Qsib", "Qsid", "Taub", "Taud",
		"Qsns", "Qsnc", "Qlns", "Qlnc", "Vz", "Rkinsc", "Rkinc", "Inmax", "intc", "ieff", "Ur", "Wc", "Tc", "Tac", "QHc",
		"QEc", "Ec", "Qpc", "Qmc", "Mc", "FMc", "SWIGM", "SWISM", "SWIR", "errMB" };
	int outvarindx = 17, aggoutvarindx = 17;
	int size, rank, irank, jrank;
	double intermStart_Time = 0.0, startTimeT = 0.0, TotalTime = 0.0, paramSite_Time = 0.0, inputTS_Time = 0.0, computeRun_Time = 0.0, outputWrite_Time = 0.0;
	double TsReadTime = 0.0, TSStartTime, ComputeStartTime, ComputeTime = 0.0, OutWriteTime;

	MPI::Init(argc, argv);
	//how many processes
	size = MPI::COMM_WORLD.Get_size(); //	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//which rank is yours? 
	rank = MPI::COMM_WORLD.Get_rank(); //_Comm_rank(MPI_COMM_WORLD,&rank);
	//cout << "\n rank "<< rank << " of "<< size << " processes has started\n" << endl;	
	MPI::Intracomm worldComm = MPI::COMM_WORLD;
	MPI::Info worldInfo = MPI::INFO_NULL;
	if (rank == 0)
	{
		//microsecond wall time: to time block of work
		startTimeT = MPI::Wtime();
		intermStart_Time = MPI::Wtime();
		TsReadTime = 0.0;
		ComputeTime = 0.0;
	}
	//  Input Arguments		
	if (argc > 1)
	{
		//conFile = new char[sizeof(argv[0])];
		strcpy(conFile, argv[1]);
	}
	else
	{
		if (rank == 0)
			cout << "file not found exiting" << endl;
		MPI::Finalize();
		return 1;
		//cin >> conFile;
	}

	FILE* pconFile = fopen(conFile, "rt");
	fgets(headerLine, 256, pconFile);
	fscanf(pconFile, "%s\n %s\n %s\n %s\n %s\n %s\n", paramFile1, sitevarFile1, inputconFile1, outputconFile1, aggoutputFile1, watershedFile1);
	fscanf(pconFile, "%s %s %s\n", wsvarName1, wsycorName1, wsxcorName1);
	
	//new vs2012 appears to have issues with passing char[256] for const char*
	const char* paramFile = paramFile1, * sitevarFile = sitevarFile1, * inputconFile = inputconFile1, * outputconFile = outputconFile1, * aggoutputFile = aggoutputFile1,
		* watershedFile = watershedFile1, * wsvarName = wsvarName1, * wsycorName = wsycorName1, * wsxcorName = wsxcorName1;

	//read simulation related parameters including start and end datetimes, and model time step dt
	fscanf(pconFile, "%d %d %d %lf\n", &ModelStartDate[0], &ModelStartDate[1], &ModelStartDate[2], &ModelStartHour);
	fscanf(pconFile, "%d %d %d %lf\n", &ModelEndDate[0], &ModelEndDate[1], &ModelEndDate[2], &ModelEndHour);
	fscanf(pconFile, "%lf\n %lf\n %d\n %d %d %d\n %d %d\n", &ModelDt, &ModelUTCOffset, &inpDailyorSubdaily, &outtStride, &outyStep, &outxStep, &outDimord, &aggoutDimord);
	//close control file
	fclose(pconFile);
	//time units
	char tunits[256];
	int hhMod = (int)floor(ModelStartHour);
	int mmMod = (int)(remainder(ModelStartHour, 1.0) * 60);
	sprintf(tunits, "hours since %d-%d-%d %d:%d:00 UTC", ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], hhMod, mmMod);
	const char* tUnitsout = tunits;

//7.25.19 for param, site variables in nc root process reads and broadcasts
		
	if (rank == 0)
	{
		//read watershed (model domain) netcdf file
		//retvalue = readwsncFile(watershedFile, wsvarName, wsycorName, wsxcorName, wsycorArray, wsxcorArray, wsArray, dimlen1, dimlen2, wsfillVal, worldComm, worldInfo);
		retvalue = readwsncFile(watershedFile, wsvarName, wsycorName, wsxcorName, wsycorArray, wsxcorArray, wsArray, dimlen1, dimlen2, wsfillVal); // , worldComm, worldInfo);

		//cout<<"dim1 = "<<dimlen1<<" dim2 = "<< dimlen2<<endl;
		/*printf("fillvalue= %d ",wsfillVal);
		for(int i=0;i<dimlen1;i++){
		for(int j=0;j<dimlen2;j++)
		cout<<wsArray[i][j];
		cout<<"\n";
		}*/
	}
//broadcast ws values
	//wsxcorArray = create2DArray_Contiguous_int(Nydim, Nxdim);
	MPI::COMM_WORLD.Bcast(&dimlen1, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&dimlen2, 1, MPI::INT, 0);
	MPI::COMM_WORLD.Bcast(&wsfillVal, 1, MPI::INT, 0);
	if (rank != 0)
	{
		wsycorArray = new float[dimlen1];
		wsxcorArray = new float[dimlen2];
		wsArray = create2DArray_Contiguous_int(dimlen1, dimlen2);
		/*new int* [dimlen1];
		for (size_t j = 0; j < dimlen1; j++)
			wsArray[j] = new int[dimlen2];
		*/
	}
	MPI::COMM_WORLD.Bcast(&wsycorArray[0], dimlen1, MPI::FLOAT, 0);
	MPI::COMM_WORLD.Bcast(&wsxcorArray[0], dimlen2, MPI::FLOAT, 0);
	MPI::COMM_WORLD.Bcast(&wsArray[0][0], dimlen1* dimlen2, MPI::INT, 0);

	//aggregation zone info
	float * wsArray1D = new float[dimlen1*dimlen2];
	for (int i = 0; i < dimlen1; i++)
		for (int j = 0; j < dimlen2; j++)
			wsArray1D[i*dimlen2 + j] = wsArray[i][j];
	std::set<int> zValues(wsArray1D, wsArray1D + (dimlen1*dimlen2));
	//cout << zValues.size() << endl;
	//std::remove_if(zValues.begin(), zValues.end(), [&wsfillVal](int a){ return a == wsfillVal; });
	std::set<int> fillSet;
    fillSet.insert (wsfillVal);
	//cout << "fill: " << fillSet.size() << " value: " << *(fillSet.begin())<<endl;
	std::vector<int> zVal(zValues.size());
	std::vector<int>::iterator it = std::set_difference(zValues.begin(), zValues.end(), fillSet.begin(), fillSet.end(), zVal.begin());  // exclude _FillValue
	zVal.resize(it - zVal.begin());
	//cout << zVal.size()<<endl;
	z_ycor = new float[zVal.size()];
	z_xcor = new float[zVal.size()];
	//cout << zValues.size() << endl;
	nZones = zVal.size();
	for (int iz = 0; iz < zVal.size(); iz++)
	{
		//#_12.24.14 change these with actual outlet locations coordinates
		z_ycor[iz] = 0.0;
		z_xcor[iz] = 0.0;
		//cout << zValues[iz];		
	}
	//read parameters
	//readParams(paramFile,paramValues);
	readParams(paramFile, parvalArray, npar);
	/*cout<<"param read..\n ");
	for(int i=0;i<npar;i++)
	cout<<"%f ",parvalArray[i]);	*/
	//read site vars 
	//cout<<"Reading site variable ");
	readSiteVars(sitevarFile, strsvArray); //svDefaults,svFile,svVarName,svType);
	/*cout<<"\n site variables read \n");
	for(int i=0;i<32;i++)
	cout<<"%f ",strsvArray[i].svdefValue);
	cout<<"\n");*/
	for (int i = 0; i < 32; i++)
	{
		if (strsvArray[i].svType == 1)
		{
			//cout<<"%d %s %s\n",i, strsvArray[i].svFile,strsvArray[i].svVarName);
			if (rank == 0)
			{
				//retvalue = read2DNC(strsvArray[i].svFile, strsvArray[i].svVarName, strsvArray[i].svArrayValues, worldComm, worldInfo);
				retvalue = read2DNC(strsvArray[i].svFile, strsvArray[i].svVarName, strsvArray[i].svArrayValues);
				/*for(int ih=0;ih<13;ih++)
				{
				for(int jv=0;jv<16;jv++)
				cout<<"%f ",strsvArray[i].svArrayValues[ih][jv]);
				cout<<"\n");
				}*/
			}
			else 
			{
				strsvArray[i].svArrayValues = create2DArray_Contiguous(dimlen1, dimlen2);
				/*strsvArray[i].svArrayValues = new float* [dimlen1];
				for (size_t j = 0; j < dimlen1; j++)
					strsvArray[i].svArrayValues[j] = new float[dimlen2];
				*/
			}
			MPI::COMM_WORLD.Bcast(&strsvArray[i].svArrayValues[0][0], dimlen1* dimlen2, MPI::FLOAT, 0);
		}
	}
	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	//read input /forcing control file--all possible entries of input control have to be provided
	readInputForcVars(inputconFile, strinpforcArray);
	modelSpan = julian(ModelEndDate[0], ModelEndDate[1], ModelEndDate[2], ModelEndHour) - julian(ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], ModelStartHour);
	//model time steps
	numTimeStep = (int)ceil(modelSpan*(24 / ModelDt));
	if (rank == 0)
		cout << "number of time steps: " << " " << numTimeStep << endl;
	//read time series forcing data only once outside of the main loop
	for (int it = 0; it < 13; it++)
	{
		if (strinpforcArray[it].infType == 0)
			readTextData(strinpforcArray[it].infFile, tsvarArray[it], ntimesteps[0]);   //ntimesteps[0] 12.18.14
		else if (strinpforcArray[it].infType == 2 || strinpforcArray[it].infType == -1)
		{
			//######TBC 6.20.13 better way to handle this is needed
			tsvarArray[it] = new float[2];
			ntimesteps[0] = 2;
			//just copy the default value if a single value is the option				
			tsvarArray[it][0] = strinpforcArray[it].infType;
			tsvarArray[it][1] = strinpforcArray[it].infdefValue;
		}
	}
	inputTS_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	//allocate memory for output array	
	outvarArray = new float*[numOut];
	for (int i = 0; i<numOut; i++)
		outvarArray[i] = new float[numTimeStep];
	//total grid size to compute progress
	totalgrid = dimlen1*dimlen2;
	//output control
	readOutputControl(outputconFile, pOut, ncOut, aggOut, npout, nncout, naggout);
	//create output netcdf
	outtSteps = numTimeStep / outtStride;             //save SWE every outstrid'th t-step
	t_out = new float[outtSteps];
	for (int it = 0; it < outtSteps; ++it)
		t_out[it] = it*outtStride*ModelDt;      //in hours since model start time 
	//vector of active cells
	std::vector<std::pair<int, int>> activeCells;
	for (int iy = 0; iy < dimlen1; iy++)
		for (int jx = 0; jx < dimlen2; jx++)
			if (wsArray[iy][jx] != wsfillVal && strsvArray[16].svType != 3)  //compute cell && no accumulation zone //***tbc what happens if it is accumulation zone?
				activeCells.push_back(std::make_pair(iy, jx));
	//output array written to netcdf files
	ncoutArray = new float**[nncout];
	for (int inc = 0; inc < nncout; inc++)
	{
		ncoutArray[inc] = new float*[activeCells.size() / size + 1];
		for (int nindx = 0; nindx < (activeCells.size() / size + 1); nindx++)
			ncoutArray[inc][nindx] = new float[outtSteps];		
	}
	yIndxArr = new int[activeCells.size() / size + 1];
	xIndxArr = new int[activeCells.size() / size + 1];
	//aggregated output arrays 
	aggoutvarArray = new float**[nZones];
	float * totalAgg = new float[outtSteps];
	ZonesArr = new int[nZones];
	for (int j = 0; j < nZones; j++)
	{
		ZonesArr[j] = 0;
		aggoutvarArray[j] = new float*[naggout];
		for (int i = 0; i < naggout; i++)
		{
			aggoutvarArray[j][i] = new float[outtSteps];
			for (int it = 0; it < outtSteps; it++)
				aggoutvarArray[j][i][it] = 0.0;
		}
	}
	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	
	for (int icout = 0; icout < nncout; icout++)
		retvalue = create3DNC_uebOutputs(ncOut[icout].outfName, (const char*)ncOut[icout].symbol, (const char*)ncOut[icout].units, tNameout, tUnitsout,
			tlong_name, tcalendar, outtSteps, outDimord, t_out, &out_fillVal, watershedFile, wsvarName, wsycorName, wsxcorName, worldComm, worldInfo);
	//create aggregate ouput file
	retvalue = create3DNC_uebAggregatedOutputs(aggoutputFile, aggOut, naggout, tNameout, tUnitsout, tlong_name, tcalendar, outtSteps, aggoutDimord, t_out, &out_fillVal,
		watershedFile, wsvarName, wsycorName, wsxcorName, nZones, zName, z_ycor, z_xcor, worldComm, worldInfo);

	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	//MPI::COMM_WORLD.Barrier();	
	int remLength = activeCells.size() % size;	
	int cellIndx = -1;
	for (irank = rank; irank < activeCells.size() - remLength; irank += size)
	{		
		//track grid cell
		cellIndx++;
		uebCellY = activeCells[irank].first;
		uebCellX = activeCells[irank].second;
		yIndxArr[cellIndx] = uebCellY;
		xIndxArr[cellIndx] = uebCellX;
		for (int is = 0; is < 32; is++)
		{
			if (strsvArray[is].svType == 1)
				SiteState[is] = strsvArray[is].svArrayValues[uebCellY][uebCellX];
			else
				SiteState[is] = strsvArray[is].svdefValue;
		}
		paramSite_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		for (int it = 0; it < 13; it++)
		{
			if (strinpforcArray[it].infType == 1)
			{
				ncTotaltimestep = 0;
				for (int numNc = 0; numNc < strinpforcArray[it].numNcfiles; numNc++) //if multiple netcdf for a single variable, they are read one by one and copied to single array
				{
					//read 3D netcdf (regridded array processed by uebInputs)
					char numtoStr[256];
					sprintf(numtoStr, "%d", numNc);
					char tsInputfile[256];
					strcpy(tsInputfile, strinpforcArray[it].infFile);
					strcat(tsInputfile, numtoStr);
					strcat(tsInputfile, ".nc");
					//cout<<"%s\n",tsInputfile);
					retvalue = readNC_TS(tsInputfile, strinpforcArray[it].infvarName, strinpforcArray[it].inftimeVar,
						wsycorName, wsxcorName, tsvarArrayTemp[numNc], tcorvar[it], uebCellY, uebCellX, ntimesteps[numNc],worldComm, worldInfo);
					ncTotaltimestep += ntimesteps[numNc];
					/*for(int tps=0;tps<ncTotaltimestep;tps++)
						cout << "  " << tsvarArrayTemp[numNc][xstrt][tps];
						cout<<"  "<<ncTotaltimestep<<endl;*/
				}
				tsvarArray[it] = new float[ncTotaltimestep];
				tinitTime = 0;
				for (int numNc = 0; numNc < strinpforcArray[it].numNcfiles; numNc++)
				{
					for (int tts = 0; tts < ntimesteps[numNc]; tts++)
						tsvarArray[it][tts + tinitTime] = tsvarArrayTemp[numNc][tts];
					tinitTime += ntimesteps[numNc];
				}
			}
		}
		inputTS_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		/*cout << endl << "Tmin " << endl;
		for (int tps = 0; tps < 100; tps++)
			cout << "  " << tsvarArray[10][tps] << " ";*/
		RUNUEB(tsvarArray, SiteState, parvalArray, outvarArray, ModelStartDate, ModelStartHour, ModelEndDate, ModelEndHour, ModelDt, ModelUTCOffset);
		computeRun_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		//write nc outputs
		for (int icout = 0; icout < nncout; icout++)
		{
			for (int vindx = 0; vindx < 70; vindx++)
			{
				if (strcmp(ncOut[icout].symbol, uebVars[vindx]) == 0)
				{
					outvarindx = vindx;
					break;
				}
			}
			for (int it = 0; it < outtSteps; ++it)
				ncoutArray[icout][cellIndx][it] = outvarArray[outvarindx][outtStride*it];        //t_out[it]3.20.15  //use timeStiride to sample outputs if it is dense (e.g hourly data for a year may be too big to save in one nc file)
			//write var values
			//retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, uebCellY, uebCellX, outtSteps, t_out, worldComm, worldInfo);
		}
		//point outputs
		for (int ipout = 0; ipout < npout; ipout++)
		{
			if (uebCellY == pOut[ipout].ycoord && uebCellX == pOut[ipout].xcoord)
			{
				FILE* pointoutFile = fopen((const char*)pOut[ipout].outfName, "w");
				for (int istep = 0; istep < numTimeStep; istep++)
				{
					fprintf(pointoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
					for (int vnum = 4; vnum < 70; vnum++)
						fprintf(pointoutFile, " %16.6f ", outvarArray[vnum][istep]);
				}
				fclose(pointoutFile);
			}
		}
		//#_??aggregated outputs 12.24.14
		zoneid = wsArray[uebCellY][uebCellX] - 1;
		ZonesArr[zoneid] += 1;
		for (int iagout = 0; iagout < naggout; iagout++)
		{
			for (int vindx = 0; vindx < 70; vindx++)
			{
				if (strcmp(aggOut[iagout].symbol, uebVars[vindx]) == 0)
				{
					aggoutvarindx = vindx;
					break;
				}
			}
			for (int it = 0; it < outtSteps; it++)
				aggoutvarArray[zoneid][iagout][it] += outvarArray[aggoutvarindx][outtStride*it];
		}
		outputWrite_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		//debug outputs
		/*if (irank % outyStep == 0 && (jrank + xstrt) % outxStep == 0)
		{
			char testPrint[256];
			char ind[256];
			strcpy(testPrint, "ZTest");
			sprintf(ind, "%d", irank);
			strcat(testPrint, ind);
			strcat(testPrint, "_");
			sprintf(ind, "%d", jrank + xstrt);
			strcat(testPrint, ind);
			strcat(testPrint, ".txt");
			FILE* testoutFile = fopen(testPrint, "w");
			for (int istep = 0; istep < numTimeStep; istep++)
			{
				fprintf(testoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
				for (int vnum = 4; vnum < 70; vnum++)
				fprintf(testoutFile, " %16.6f ", outvarArray[vnum][istep]);
			}
			fclose(testoutFile);
		}*/
		//grid count progress is calculated and written here
		numgrid += size;
		if (rank == 0 && numgrid % dimlen1 == 0 )
			cout << "   percent completed: " << ((float)numgrid / activeCells.size())*100.0 << " %" << endl;
		fflush(stdout);
	} // 
	cellIndx++;
	//cout << "number of active cells = " << activeCells.size() << " number cellIndx = " << cellIndx << endl;
	for (int icout = 0; icout < nncout; icout++)
	{
		//write var values
		retvalue = WriteTSto3DNC_Block((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, yIndxArr, xIndxArr, cellIndx, outtSteps, ncoutArray[icout], worldComm, worldInfo);               
	}
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
	//cout << "Process " << rank << " before final loop with rem Length " << remLength << endl;
	//MPI::COMM_WORLD.Barrier();
	if (remLength > 0)  //if there are remaining compute cells after even distribution 
	{
		int *remRanks = new int[remLength];
		int leftBorder = activeCells.size() - remLength;
		for (int ir = 0; ir < remLength; ir++)           // = leftBorder; ir < activeCells.size(); ir++)
			remRanks[ir] = ir;                           // [ir - leftBorder] = ir;
		MPI::Group worldGroup = MPI::COMM_WORLD.Get_group();
		MPI::Group remGroup = worldGroup.Incl(remLength, remRanks);
		MPI::Intracomm remComm = MPI::COMM_WORLD.Create(remGroup);
		int newSize = -1;  //remComm.Get_size();
		int newRank = 0;  //remComm.Get_rank();
		//cout << " rank "<< rank << " of "<< size << " processes has started\n" << endl;	
		if (rank < remLength)         //only for processes in the new comm group
		{
			newSize = remComm.Get_size();
			newRank = remComm.Get_rank();
			//cout << " new rank of old process "<< rank << " is:  "<< newRank << endl;	
			//track grid cell
			uebCellY = activeCells[leftBorder + newRank].first;
			uebCellX = activeCells[leftBorder + newRank].second;
			for (int is = 0; is < 32; is++)
			{
				if (strsvArray[is].svType == 1)
					SiteState[is] = strsvArray[is].svArrayValues[uebCellY][uebCellX];
				else
					SiteState[is] = strsvArray[is].svdefValue;
			}
			paramSite_Time += (MPI::Wtime() - intermStart_Time);
		    intermStart_Time = MPI::Wtime();
			for (int it = 0; it < 13; it++)            //it < 13,          12.18.14
		    {
			    if (strinpforcArray[it].infType == 1)     // == 0
			    {
					ncTotaltimestep = 0;
					for (int numNc = 0; numNc < strinpforcArray[it].numNcfiles; numNc++) //if multiple netcdf for a single variable, they are read one by one and copied to single array
					{
						//read 3D netcdf (regridded array processed by uebInputs)
						char numtoStr[256];
						sprintf(numtoStr, "%d", numNc);
						char tsInputfile[256];
						strcpy(tsInputfile, strinpforcArray[it].infFile);
						strcat(tsInputfile, numtoStr);
						strcat(tsInputfile, ".nc");
						//cout<<"%s\n",tsInputfile);
						retvalue = readNC_TS(tsInputfile, strinpforcArray[it].infvarName, strinpforcArray[it].inftimeVar,
							wsycorName, wsxcorName, tsvarArrayTemp[numNc], tcorvar[it], uebCellY, uebCellX, ntimesteps[numNc], remComm, worldInfo);
						ncTotaltimestep += ntimesteps[numNc];
						/*for(int tps=0;tps<ncTotaltimestep;tps++)
							cout << "  " << tsvarArrayTemp[numNc][xstrt][tps];
						cout<<"  "<<ncTotaltimestep<<endl;*/
					}
					tsvarArray[it] = new float[ncTotaltimestep];
					tinitTime = 0;
					for (int numNc = 0; numNc < strinpforcArray[it].numNcfiles; numNc++)
					{
						for (int tts = 0; tts < ntimesteps[numNc]; tts++)
							tsvarArray[it][tts + tinitTime] = tsvarArrayTemp[numNc][tts];
						tinitTime += ntimesteps[numNc];
					}
				}
			}
			inputTS_Time += (MPI::Wtime() - intermStart_Time);
			intermStart_Time = MPI::Wtime();
			/*cout << endl << "Tmin " << endl;
			for (int tps = 0; tps < 100; tps++)
			cout << "  " << tsvarArray[10][tps] << " ";*/
			RUNUEB(tsvarArray, SiteState, parvalArray, outvarArray, ModelStartDate, ModelStartHour, ModelEndDate, ModelEndHour, ModelDt, ModelUTCOffset);
			computeRun_Time += (MPI::Wtime() - intermStart_Time);
			intermStart_Time = MPI::Wtime();
			//write nc outputs
			for (int icout = 0; icout < nncout; icout++)
			{
				for (int vindx = 0; vindx < 70; vindx++)
				{
					if (strcmp(ncOut[icout].symbol, uebVars[vindx]) == 0)
					{
						outvarindx = vindx;
						break;
					}
				}
				for (int it = 0; it < outtSteps; ++it)
					t_out[it] = outvarArray[outvarindx][outtStride*it];         //use timeStiride to sample outputs if it is dense (e.g hourly data for a year may be too big to save in one nc file)
				//write var values
				retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, uebCellY, uebCellX, outtSteps, t_out, remComm, worldInfo);
			}
			//point outputs
			for (int ipout = 0; ipout < npout; ipout++)
			{
				if (uebCellY == pOut[ipout].ycoord && uebCellX == pOut[ipout].xcoord)
				{
					FILE* pointoutFile = fopen((const char*)pOut[ipout].outfName, "w");
					for (int istep = 0; istep < numTimeStep; istep++)
					{
						fprintf(pointoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
						for (int vnum = 4; vnum < 70; vnum++)
							fprintf(pointoutFile, " %16.6f ", outvarArray[vnum][istep]);
					}
					fclose(pointoutFile);
				}
			}
			//#_??aggregated outputs 12.24.14
			zoneid = wsArray[uebCellY][uebCellX] - 1;
			ZonesArr[zoneid] += 1;
			for (int iagout = 0; iagout < naggout; iagout++)
			{
				for (int vindx = 0; vindx < 70; vindx++)
				{
					if (strcmp(aggOut[iagout].symbol, uebVars[vindx]) == 0)
					{
						aggoutvarindx = vindx;
						break;
					}
				}
				for (int it = 0; it < outtSteps; it++)
					aggoutvarArray[zoneid][iagout][it] += outvarArray[aggoutvarindx][outtStride*it];
			}
			outputWrite_Time += (MPI::Wtime() - intermStart_Time);
			intermStart_Time = MPI::Wtime();
			//debug outputs
			/*if (irank % outyStep == 0 && (jrank + xstrt) % outxStep == 0)
			{
				char testPrint[256];
				char ind[256];
				strcpy(testPrint, "ZTest");
				sprintf(ind, "%d", irank);
				strcat(testPrint, ind);
				strcat(testPrint, "_");
				sprintf(ind, "%d", jrank + xstrt);
				strcat(testPrint, ind);
				strcat(testPrint, ".txt");
				FILE* testoutFile = fopen(testPrint, "w");
				for (int istep = 0; istep < numTimeStep; istep++)
				{
					fprintf(testoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
					for (int vnum = 4; vnum < 70; vnum++)
					fprintf(testoutFile, " %16.6f ", outvarArray[vnum][istep]);
				}
				fclose(testoutFile);
			}*/
			//grid count progress is calculated and written here
			numgrid += newSize;
			if (rank == 0)
				cout << "\r   percent completed: " << ((float)numgrid / activeCells.size())*100.0 << " %" << endl;
			fflush(stdout);
		} // if(rank < remLength)
	} //if (remLength > 0)		
	//aggregation/ reduction 
	//cout << "process " << rank << " completed computation" << endl;
	for (int it = 0; it < outtSteps; it++)
		totalAgg[it] = 0.0;
	MPI::COMM_WORLD.Barrier();
	int rankrec = 0;               //receiver rank 
	int totalZonecells = 1, zonValue = 0;
	for (int izone = 0; izone < nZones; izone++)
	{
		rankrec = izone*size / nZones;
		//cout << "process " << rank << " before first reduce to rank: " << rankrec << endl;
		zonValue = ZonesArr[izone];
		MPI::COMM_WORLD.Reduce(&zonValue, &totalZonecells, 1, MPI::INT, MPI::SUM, rankrec);
		//cout<<"process "<<rank<<" total zone cells "<<totalZonecells<<endl;
		if (totalZonecells < 1)
			totalZonecells = 1;
		for (int iagout = 0; iagout < naggout; iagout++)
		{
			//cout << "process " << rank << " before reduce of output " << iagout << endl;
			//if (rank == rankrec) MPI::COMM_WORLD.Reduce(MPI::IN_PLACE, aggoutvarArray[izone][iagout],outtSteps, MPI::FLOAT, MPI::SUM, rankrec); else 
			MPI::COMM_WORLD.Reduce(aggoutvarArray[izone][iagout], totalAgg, outtSteps, MPI::FLOAT, MPI::SUM, rankrec);
			//cout << "process " << rank << " waiting for writing" << endl;
			//#_12.28.14 aggregation operation needs defining
			if (rank == rankrec)
			{
				if (strcmp(aggOut[iagout].aggop, "AVE") == 0)
					for (int it = 0; it < outtSteps; it++)
						totalAgg[it] = totalAgg[it] / totalZonecells;					
				//cout << "process " << rank << " before write of output " << iagout << " for zone: " << izone << endl;
				retvalue = Write_uebaggTS_toNC(aggoutputFile, aggOut[iagout].symbol, aggoutDimord, izone, outtSteps, totalAgg); 
				//cout << "process: " << rank << " done writing output: " << iagout << " for zone " << izone << endl;
			}
		}
	}
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	//MPI::COMM_WORLD.Barrier();
	//cout<<"Process "<<rank<<" starting deallocating memory"<<endl;	

	//deallocate memory ====#_*_#______Needs revisiting; some of the arrays are not deleted 6.23.13
	
	delete2DArray_Contiguous_int(wsArray);
	/*for (int i = 0; i < dimlen1; i++)
		delete[] wsArray[i];
	delete[] wsArray;
	*/

	delete[] parvalArray;
	for (int i = 0; i < 32; i++)
	{
		if (strsvArray[i].svType == 1)
		{
			delete2DArray_Contiguous(strsvArray[i].svArrayValues);
			/*for (int j = 0; j < dimlen1; j++)
				delete[] strsvArray[i].svArrayValues[j];
			delete[] strsvArray[i].svArrayValues;
			*/
		}
	}
	delete[] strsvArray;
	//delete[] tsvarArray[kx];
	for (int it = 0; it < 13; it++)       //10-->13   6.26.14
	{
		delete[] tsvarArray[it];
	}
	//delete[] tsvarArray;
	/*for (int it = 0; it < 5; it++)       //  6.26.14
	{
		delete[] tsvarArrayTemp[it];
	}*/
	//delete[] tsvarArrayTemp;
	/*if (rank < remLength)
	{
		for (int it = 0; it < 13; it++)
		{
			if (strinpforcArray[it].infType == 1)
				delete[] tcorvar[it];
		}
		//delete[] tcorvar;
	}*/
	for (int zk = 0; zk < nZones; zk++)
	{
		for (int ig = 0; ig < naggout; ig++)
			delete[] aggoutvarArray[zk][ig];
		delete[] aggoutvarArray[zk];
	}
	delete[] aggoutvarArray;
	/*for(int k=0 ;k<numOut; k++)
		delete[] outvarArray[k];
	delete []outvarArray; */
	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	cout << "Process " << rank << " finished" << endl;
	fflush(stdout);	
	MPI::COMM_WORLD.Barrier();
	if (rank == 0)
	{
		TotalTime = MPI::Wtime() - startTimeT;                   //(float)1000*(endTimeT - startTimeT)/CLOCKS_PER_SEC;
		cout << "Time in seconds" << endl;
		cout << "Reading param  site state input control:  " << paramSite_Time << endl;
		cout << "Reading input TS txt arrays:  " << TsReadTime << endl;
		cout << "Reading input total TS arrays:  " << inputTS_Time << endl;
		cout << "Model simulation run time:  " << computeRun_Time << endl;
		cout << "Outputs write time: "<<outputWrite_Time<<endl;
		cout << "Total time of including overhead :  " << TotalTime << endl;
		cout << "Done! return value: " << retvalue << endl;
		fflush(stdout);
	}
exitlab:
	MPI::Finalize();
	//getchar();
	return 0;
}
