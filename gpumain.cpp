//#include "uebpg.h"
//#include "nctest.h"
//#include"mpi.h"
#include "gpuuebpgdecls.h"
#include <time.h>
//#include <queue>
#pragma warning(disable : 4996)
using namespace std;

__global__ void callUEBRun(uebCell *uebCellArray, int nCells)
{
	int indx = blockIdx.x*blockDim.x + threadIdx.x;
	if (indx < nCells)
	{
		//float *dev_loOutArr = NULL;
		//cudaMalloc(&dev_loOutArr, 700000 * sizeof(float));		
		uebCellArray[indx].runUEB();
	}
}

__device__ __host__ void cuda_checkERR(cudaError_t err)
{
	if (err != cudaSuccess){
	   std::cout << "Error: " << cudaGetErrorString(err) << endl;
	   exit(EXIT_FAILURE);
    }
}

__host__ __device__ void checkDeviceMemory()
{
	//mem check on device
	size_t freeM, totalM;
	float freeMB, totalMB, allocMB;
	cudaMemGetInfo((size_t*)&freeM, (size_t*)&totalM);
	freeMB = (size_t)freeM / (1024*1024);
	totalMB = (size_t)totalM / (1024*1024);
	allocMB = totalMB - freeMB;
	printf(" %f  MB of   %f   MB total available device memory allocated. Remaining memory =   %f MB\n", allocMB, totalMB, freeMB);
}

__host__ __device__ void estimateThroughput(size_t dataSize, clock_t beginTime, clock_t endTime)
{
	double bandWidth = dataSize * 2.0;
	double GFLOPs = (double)(bandWidth * CLOCKS_PER_SEC) / (double)(endTime - beginTime);
	printf(" Estimated throughput =  %lf GFLOPs\n", GFLOPs);
}

int main(int argc, char* argv[])
{
	//gpu control	
	int threadsPerBlock = 255, blocksPerGrid = 1;
	//time
	float timeControl = 0.0, timeWS = 0.0, timeSitestate = 0.0, timeTSArrays = 0.0, timeParam = 0.0, timeParamSiteInptcontrol = 0.0, timeModelRun = 0.0;
	float* OutVarValues; //= new float*[70]; //[70];
	float*** aggoutvarArray = NULL;
	float ***ncoutArray = NULL;
	char conFile[256], paramFile1[256], sitevarFile1[256], inputconFile1[256], outputconFile1[256], watershedFile1[256], aggoutputconFile1[256], aggoutputFile1[256];
	char wsvarName1[256], wsycorName1[256], wsxcorName1[256];
	int **wsArray = NULL;
	int dimlen1 = 0, dimlen2 = 0;
	int wsfillVal = -9999;	
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
	//inpforcvar strinpforcArray[13];
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
	float ***RegArray[13] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };  //3.19.15  [5]; // [xstride]; //assuming max nc files for a variable =5
	//float *tcorvar[13], *tsvarArray[13],
	float *ycorArr = NULL, *xcorArr = NULL;
	int ncTotaltimestep[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int numNc[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };				//6.24.14	
	size_t ntimesteps = 0, nysteps = 0, nxsteps = 0;
	int tinitTime = 0;
	int npar = 32;
	/*int NUMtimeSTEP,NREFYR,NREFMO,NREFDAY,NumOP;*/
	int ModelStartDate[3], ModelEndDate[3]; //check this 
	double ModelStartHour, ModelEndHour, ModelDt, ModelUTCOffset;
	double modelSpan;	
	int inpDailyorSubdaily;
	//int numTimeStep;
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
		"cumMr", "Qnet", "smelt", "refDepth", "totalRefDepth", "cf", "Taufb", "Taufd", "Qsib", "Qsid", "Taub", "Taud",
		"Qsns", "Qsnc", "Qlns", "Qlnc", "Vz", "Rkinsc", "Rkinc", "Inmax", "intc", "ieff", "Ur", "Wc", "Tc", "Tac", "QHc",
		"QEc", "Ec", "Qpc", "Qmc", "Mc", "FMc", "SWIGM", "SWISM", "SWIR", "errMB" };
	int outvarindx = 17, aggoutvarindx = 17;
	int size =1, rank = 0, irank = 0, jrank;
	double intermStart_Time = 0.0, startTimeT = 0.0, TotalTime = 0.0, paramSite_Time = 0.0, inputTS_Time = 0.0, computeRun_Time = 0.0, outputWrite_Time = 0.0, dataCopy_Time = 0.0;
	double TsReadTime = 0.0, TSStartTime, ComputeStartTime, ComputeTime = 0.0, OutWriteTime;
	clock_t beginTime, endTime;
	//beginTime = clock();
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
	const char *paramFile = paramFile1, *sitevarFile = sitevarFile1, *inputconFile = inputconFile1, *outputconFile = outputconFile1, *aggoutputFile = aggoutputFile1,
		*watershedFile = watershedFile1, *wsvarName = wsvarName1, *wsycorName = wsycorName1, *wsxcorName = wsxcorName1;
	//read simulation related parameters including start and end datetimes, and model time step dt
	fscanf(pconFile, "%d %d %d %lf\n", &ModelStartDate[0], &ModelStartDate[1], &ModelStartDate[2], &ModelStartHour);
	fscanf(pconFile, "%d %d %d %lf\n", &ModelEndDate[0], &ModelEndDate[1], &ModelEndDate[2], &ModelEndHour);
	fscanf(pconFile, "%lf\n %lf\n %d\n %d %d %d\n %d %d\n %d\n", &ModelDt, &ModelUTCOffset, &inpDailyorSubdaily, &outtStride, &outyStep, &outxStep, &outDimord, &aggoutDimord, &threadsPerBlock);
	//close control file
	fclose(pconFile);
	//time units
	char tunits[256];
	int hhMod = (int)floor(ModelStartHour);
	int mmMod = (int)(remainder(ModelStartHour, 1.0) * 60);
	sprintf(tunits, "hours since %d-%d-%d %d:%d:00 UTC", ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], hhMod, mmMod);
	const char* tUnitsout = tunits;
	//read watershed (model domain) netcdf file	
	retvalue = readwsncFile(watershedFile, wsvarName, wsycorName, wsxcorName, wsycorArray, wsxcorArray, wsArray, dimlen1, dimlen2, wsfillVal, worldComm, worldInfo);
	//cout<<"dim1 = "<<dimlen1<<" dim2 = "<< dimlen2<<endl;
	/*printf("fillvalue= %d ",wsfillVal);
	for(int i=0;i<dimlen1;i++){
	for(int j=0;j<dimlen2;j++)
	cout<<wsArray[i][j];
	cout<<"\n";
	}*/	
	//aggregation zone info
	float * wsArray1D = new float[dimlen1*dimlen2];
	for (int i = 0; i < dimlen1; i++)
		for (int j = 0; j < dimlen2; j++)
			wsArray1D[i*dimlen2 + j] = wsArray[i][j];
	//set contains unique id values
	std::set<int> zValues(wsArray1D, wsArray1D + (dimlen1*dimlen2));
	//cout << zValues.size() << endl;
	//std::remove_if(zValues.begin(), zValues.end(), [&wsfillVal](int a){ return a == wsfillVal; });
	std::set<int> fillSet;
    fillSet.insert (wsfillVal);
	//cout << "fill: " << fillSet.size() << " value: " << *(fillSet.begin())<<endl;
	std::vector<int> zVal(zValues.size());
	std::vector<int>::iterator it = std::set_difference(zValues.begin(), zValues.end(), fillSet.begin(), fillSet.end(), zVal.begin());  // exclude _FillValue
	zVal.resize(it - zVal.begin());  //now zVal contains unique watershed ids excluding fill value
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
	//read site vars 
	//cout<<"Reading site variable ");
	readSiteVars(sitevarFile, strsvArray); //svDefaults,svFile,svVarName,svType);
	/*cout<<"\n site variables read \n");
	for(int i=0;i<32;i++)
	cout<<"%f ",strsvArray[i].svdefValue);
	cout<<"\n");*/
	for (int i = 0; i < 32; i++)
		if (strsvArray[i].svType == 1)
		{
			//cout<<"%d %s %s\n",i, strsvArray[i].svFile,strsvArray[i].svVarName);
			retvalue = read2DNC(strsvArray[i].svFile, strsvArray[i].svVarName, strsvArray[i].svArrayValues, worldComm, worldInfo);
			/*for(int ih=0;ih<13;ih++)
			{
				for(int jv=0;jv<16;jv++)
					cout<<"%f ",strsvArray[i].svArrayValues[ih][jv]);
				cout<<"\n");
			}*/
	    }
	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
	//vector of active cells
	std::vector<std::pair<int, int> > activeCells;
	for (int iy = 0; iy < dimlen1; iy++)
		for (int jx = 0; jx < dimlen2; jx++)
			if (wsArray[iy][jx] != wsfillVal && strsvArray[16].svType != 3)  //compute cell && no accumulation zone //***tbc what happens if it is accumulation zone?
				activeCells.push_back(std::make_pair(iy, jx));	
	// create ueb model gridcell instance and copy to arrays of grid cells
	uebCell objCell0 (paramFile, ModelStartDate, ModelEndDate, ModelStartHour, ModelEndHour, ModelDt, ModelUTCOffset, inpDailyorSubdaily,outtStride);	

	//output control
	readOutputControl(outputconFile, pOut, ncOut, aggOut, npout, nncout, naggout);
	//create output netcdf
	outtSteps = objCell0.numTimeStep / outtStride;             //save SWE every outstrid'th t-step
	if (rank == 0)
		cout << "number of time steps: " << " " << objCell0.numTimeStep << endl;
	t_out = new float[outtSteps];
	for (int it = 0; it < outtSteps; ++it)
		t_out[it] = it*outtStride*ModelDt;      //in hours since model start time 

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
	//aggregated output arrays # There should be better way than this
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
	//netcdf output files
	for (int icout = 0; icout < nncout; icout++)
		retvalue = create3DNC_uebOutputs(ncOut[icout].outfName, (const char*)ncOut[icout].symbol, (const char*)ncOut[icout].units, tNameout, tUnitsout,
		tlong_name, tcalendar, outtSteps, outDimord, t_out, &out_fillVal, watershedFile, wsvarName, wsycorName, wsxcorName, worldComm, worldInfo);
	//create aggregate ouput file
	retvalue = create3DNC_uebAggregatedOutputs_par(aggoutputFile, aggOut, naggout, tNameout, tUnitsout, tlong_name, tcalendar, outtSteps, aggoutDimord, t_out, &out_fillVal,
		watershedFile, wsvarName, wsycorName, wsxcorName, nZones, zName, z_ycor, z_xcor, worldComm, worldInfo);
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	
	//read input /forcing control file--all possible entries of input control have to be provided
	objCell0.readInputForContr(inputconFile);	
	if (rank == 0){
		//cout << " start date = " << currentModelDateTime << " end date = " << EJD;
		cout << "number of uebcells = " << activeCells.size() << endl;
		cout << " size of uebCell class =  " << sizeof(uebCell) / ((double)(1024 * 1024)) << " MB" << endl;
	}
	for (int it = 0; it < 13; it++) {
		if (objCell0.infrContArr[it].infType == 0) {
			RegArray[it] = new float**[1];
			RegArray[it][0] = new float *[1];
			readTStextFile(objCell0.infrContArr[it].infFile, RegArray[it][0][0], ncTotaltimestep[it]);   //tsvarArray[it][0] 3.19.15      ntimesteps[0] 12.18.14
		}
	}
	TsReadTime += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
    cout<<"proc "<<rank<<" after reading time series"<<endl;
	int remLength = activeCells.size() % size;
	int extraCell = 0;
	if (remLength != 0) {
		if (rank < remLength)              //       (rank/remLength) == 0
		  extraCell = 1;
    }
	cout<<"proc "<<rank<<" creating uebCell arrays"<<endl;
	int numCells = (activeCells.size() / size) + extraCell;
	cout<<"proc "<<rank<<" num cells = "<<numCells<<endl;
	uebCell *uebCellArr = new uebCell[numCells];
	cout<<"proc "<<" before setting site vars"<<endl;   
	//for (irank = rank; irank < activeCells.size() - remLength; irank += size)
	int cellIndx = 0;
	for (irank = rank; irank < activeCells.size(); irank +=size)
	{
		//track grid cell		
		yIndxArr[cellIndx] = activeCells[irank].first;
		xIndxArr[cellIndx] = activeCells[irank].second;
		uebCellArr[cellIndx] = objCell0;
		uebCellArr[cellIndx].uebCellY = activeCells[irank].first;
		uebCellArr[cellIndx].uebCellX = activeCells[irank].second;
		for (int is = 0; is < 32; is++)
		{
			if (strsvArray[is].svType == 1)
				SiteState[is] = strsvArray[is].svArrayValues[uebCellArr[cellIndx].uebCellY][uebCellArr[cellIndx].uebCellX];
			else
				SiteState[is] = strsvArray[is].svdefValue;
		}
		uebCellArr[cellIndx].setSiteVars_and_Initconds(SiteState);
		cellIndx++;
	}		
	paramSite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	int curDv, curStr;    
	//comment the following out when running without gpu
	//memory on device checkDeviceMemory()
	//cuda err check
	cudaError_t err = cudaSuccess;		
	err = cudaGetDevice(&curDv);
	cuda_checkERR(err);
	cout <<"proc "<<rank<< " current device  = " << curDv << endl;
    cudaStream_t oStream;
    //cudaSetDevice(cudIndx);
	err = cudaStreamCreate(&oStream);
	cuda_checkERR(err);
	cout << "proc " << rank << " current stream  = " << oStream << endl;
	uebCell *dev_uebCellArr = NULL;
	err = cudaMalloc(&dev_uebCellArr, numCells*sizeof(uebCell));
	cuda_checkERR(err);
    cout<<"proc "<<rank<<" device memory alloc"<<endl;
	dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
	double EJD = uebCellArr[0].julian(uebCellArr[0].modelEndDate[0], uebCellArr[0].modelEndDate[1], uebCellArr[0].modelEndDate[2], uebCellArr[0].modelEndHour);
	double currentModelDateTime = uebCellArr[0].julian(uebCellArr[0].modelStartDate[0], uebCellArr[0].modelStartDate[1], uebCellArr[0].modelStartDate[2], uebCellArr[0].modelStartHour);
	if (rank ==0)
		printf(" start date = %lf  end date = %lf\n", currentModelDateTime, EJD);
	if (uebCellArr[0].inpDailyorSubdaily == 0)                     //the last 24 steps of forcing are read at once; adjust the EJD so that last time datetime is EndDate - 23*DT
		EJD -= (22 * uebCellArr[0].modelDT);
	else
		EJD -= 22.0;

	int timeOffset = 0, outOffset = 0;
	while (EJD > currentModelDateTime)
	{
		//cudaSetDevice(cudIndx);		
        printf("proc %d  current date =  %lf\n",rank, currentModelDateTime);		
		objCell0.getInpForcArr(numNc, RegArray, ncTotaltimestep, worldComm, worldInfo);
        cout<<"proc "<<rank<<" copied forcing arrays"<<endl;
		for (irank =  0; irank < numCells; irank++)
		{
			uebCellArr[irank].updateInpForcArr(RegArray, ncTotaltimestep);
		}
	    cout<<"proc "<<rank<< " forcing array updated"<<endl;
		cout<<"proc "<<rank<<" Sim numTimeSteps = "<<uebCellArr[0].numSimTimeSteps<<endl;
		inputTS_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();		
		//memory on device checkDeviceMemory()
		//if(numNc==0){
		err = cudaMemcpyAsync(dev_uebCellArr, uebCellArr, numCells*sizeof(uebCell), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
		cuda_checkERR(err);
		err = cudaStreamSynchronize(oStream);
		cuda_checkERR(err);
		//}
		dataCopy_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		cout << "process " << rank << " data copied to device" << endl; 	//cout<<err<<endl;
		//memory on device    checkDeviceMemory()
		// Launch Kernel	
		blocksPerGrid = (numCells + threadsPerBlock - 1) / threadsPerBlock;
		//call device run function	
		callUEBRun << < blocksPerGrid, threadsPerBlock, 0, oStream >> >(dev_uebCellArr, numCells);
		//synchronization
		err = cudaStreamSynchronize(oStream);   /// cudaDeviceSynchronize();
		cuda_checkERR(err);
		computeRun_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		cout << "process " << rank << " finished device compute tasks" << endl;
		//copy data back
		err = cudaMemcpyAsync(uebCellArr, dev_uebCellArr, numCells*sizeof(uebCell), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
		cuda_checkERR(err);
		err = cudaStreamSynchronize(oStream);/// cudaDeviceSynchronize();
		cuda_checkERR(err);
		dataCopy_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		cout << "process " << rank << " data copied to host" << endl;	//cout << err << endl;	
		//end of gpu call
		/*int remLength = activeCells.size() % size;
		for (irank = rank; irank < activeCells.size() - remLength; irank += size)*/
		/*for (int it = 0; it < outtSteps; ++it)
			cout << " " << OutVarValues[70 *it + 17];*/		
		for (irank = 0; irank < numCells; irank++)
		{
			//run without gpu; comment out when running with gpu
			/*uebCellArr[irank].runUEB();
			computeRun_Time += (MPI::Wtime() - intermStart_Time);
			intermStart_Time = MPI::Wtime();*/
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
				for (int it = 0; it < uebCellArr[irank].numSimTimeSteps / outtStride; it++)       // outtSteps; ++it)
					ncoutArray[icout][irank][it + timeOffset] = uebCellArr[irank].OutVarValues[outvarindx][outtStride*it + outOffset];        //t_out[it]3.20.15  //use timeStiride to sample outputs if it is dense (e.g hourly data for a year may be too big to save in one nc file)
				//write var values
				//retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, uebCellY, uebCellX, outtSteps, t_out);                //, worldComm, worldInfo);
			}
			//#_??aggregated outputs 12.24.14
			zoneid = wsArray[uebCellArr[irank].uebCellY][uebCellArr[irank].uebCellX] - 1;
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
				for (int it = 0; it < uebCellArr[irank].numSimTimeSteps / outtStride; it++)          //outtSteps; it++)
					aggoutvarArray[zoneid][iagout][it + timeOffset] += uebCellArr[irank].OutVarValues[aggoutvarindx][outtStride*it + outOffset];
			}
			//point outputs
			for (int ipout = 0; ipout < npout; ipout++)
			{
				if (uebCellArr[irank].uebCellY == pOut[ipout].ycoord && uebCellArr[irank].uebCellX == pOut[ipout].xcoord)
					uebCellArr[irank].printPointOutputs((const char*)pOut[ipout].outfName);
			}
			//debug outputs
			if (irank % (outyStep*dimlen2 + outxStep) == 0)
				uebCellArr[irank].printDebugOutputs();
		}
		outputWrite_Time += (MPI::Wtime() - intermStart_Time);
		intermStart_Time = MPI::Wtime();
		//progress is calculated and written here
		numgrid += uebCellArr[0].numSimTimeSteps;
		if (rank == 0)                               /// && numgrid % (dimlen1*uebCellArr[0].numSimTimeSteps) == 0)
			cout << "   percent completed: " << ((float)numgrid / uebCellArr[0].numTimeStep)*100.0 << " %" << endl;
		fflush(stdout);
			
		timeOffset += (uebCellArr[0].numSimTimeSteps / outtStride);
		//this takes care of the boundary between two arrays, i.e. if there are fewer than "outtSride" points at the end of the output array, they would be "padded" to the start of the next array
		outOffset = 0;           //edit later 5.12.15:  can keep outOffset = 0 as long as numSimTimeSteps is divisible by outtStride --- outtStride - (uebCellArr[0].numSimTimeSteps % outtStride);
	    //cout<<"proc "<<rank<<" before currDT compute"<<endl;
		//udate time
		currentModelDateTime = uebCellArr[0].julian(uebCellArr[0].modelStartDate[0], uebCellArr[0].modelStartDate[1], uebCellArr[0].modelStartDate[2], uebCellArr[0].modelStartHour);
        //cout<<"proc "<<rank<<" at the end of time loop"<<endl;
	}
	//cout << "number of active cells = " << activeCells.size() << " number cellIndx = " << cellIndx << endl;
	cout << "process " << rank << " completed computation" << endl;
	MPI::COMM_WORLD.Barrier();
	//cellIndx++;
	for (int icout = 0; icout < nncout; icout++)
	{
		//write var values
		retvalue = WriteTSto3DNC_Block((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, yIndxArr, xIndxArr, numCells-extraCell, outtSteps, ncoutArray[icout], worldComm, worldInfo);              //, worldComm, worldInfo);
	}
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	cout << "process " << rank << " wrote block outputs" << endl;
	//MPI::COMM_WORLD.Barrier();	
	if (remLength > 0)  //if there are remaining compute cells after even distribution 
	{
		int *remRanks = new int[remLength];	
		for (int ir = 0; ir < remLength; ir++)           // = leftBorder; ir < activeCells.size(); ir++)
			remRanks[ir] = ir;                           // [ir - leftBorder] = ir;
		MPI::Group worldGroup = MPI::COMM_WORLD.Get_group();
		MPI::Group remGroup = worldGroup.Incl(remLength, remRanks);
		MPI::Intracomm remComm = MPI::COMM_WORLD.Create(remGroup);
		int newSize = -1;  //remComm.Get_size();
		int newRank = 0;  //remComm.Get_rank();		
		if (extraCell == 1)                                    //(rank < remLength)         //only for processes in the new comm group
		{
			newSize = remComm.Get_size();
			newRank = remComm.Get_rank();
			cout << " rank " << rank << " of WorldComm has " << newRank << " of new comm group of size " << newSize << " remaining cells "<<remLength<<endl;
			//write var values
			for (int icout = 0; icout < nncout; icout++)
				retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, yIndxArr[numCells - 1], xIndxArr[numCells - 1], outtSteps, ncoutArray[icout][numCells - 1], remComm, worldInfo);
		} // if(rank < remLength)
	}

	if(extraCell==1)
		cout << "process " << rank << " wrote single outputs" << endl;
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();	
	//MPI::COMM_WORLD.Barrier();
	//aggregation/ reduction 
	for (int it = 0; it < outtSteps; it++)
		totalAgg[it] = 0.0;
	//MPI::COMM_WORLD.Barrier();
	int rankrec = 0;               //receiver rank 
	int totalZonecells = 1, zonValue = 0;
	for (int izone = 0; izone < nZones; izone++)
	{
		rankrec =  izone*size / nZones;
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
				{
					for (int it = 0; it < outtSteps; it++)
						totalAgg[it] =  totalAgg[it] / totalZonecells;				//aggoutvarArray[izone][iagout][it] / totalZonecells;	//
				}
				/*else
				{
					for (int it = 0; it < outtSteps; it++)
						totalAgg[it] = aggoutvarArray[izone][iagout][it];	// totalAgg[it] / totalZonecells;
				}*/
				//cout << "process " << rank << " before write of output " << iagout << " for zone: " << izone << endl;
				retvalue = Write_uebaggTS_toNC(aggoutputFile, aggOut[iagout].symbol, aggoutDimord, izone, outtSteps, totalAgg);
				//cout << "process: " << rank << " done writing output: " << iagout << " for zone " << izone << endl;
			}
		}
	}
	outputWrite_Time += (MPI::Wtime() - intermStart_Time);
	intermStart_Time = MPI::Wtime();
	cout << "process " << rank << " wrote aggregated outputs" << endl;
	//MPI::COMM_WORLD.Barrier();
	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);
	//cout<<"Process "<<rank<<" starting deallocating memory"<<endl;	
	//deallocate memory ====#_*_#______Needs revisiting; some of the arrays are not deleted 6.23.13	
	//free device memory
	err = cudaFree(dev_uebCellArr); // != cudaSuccess)	
	cuda_checkERR(err);
	cout << "process " << rank << " device memory freed" << endl;
	// Free host memory
	delete[] uebCellArr;
	cout << "process " << rank << " uebCell arrary freed" << endl;
	for (int i = 0; i < dimlen1; i++)
		delete[] wsArray[i];
	delete[] wsArray;
	cout<<"proc "<<rank<<" freed watershed arrays"<<endl;
	for (int i = 0; i < 32; i++)
	{
		if (strsvArray[i].svType == 1)
		{
			for (int j = 0; j < dimlen1; j++)
				delete[] strsvArray[i].svArrayValues[j];
			delete[] strsvArray[i].svArrayValues;
		}
	}
	delete[] strsvArray;
	cout<<"proc "<<rank<<" freed site vars arrays"<<endl;
	//delete[] tsvarArray[kx];
	/*for (int it = 0; it < 13; it++)       //10-->12   6.26.14
	{
		delete[] tsvarArray[it];
	}*/
	//delete[] tsvarArray;
	/*if (rank < remLength)
	{
		for (int it = 0; it < 13; it++)
		{
			if (strinpforcArray[it].infType == 1)
				delete[] tcorvar[it];
		}
		//delete[] tcorvar;
	}*/
/*	for (int it = 0; it <13; it++)       //  6.26.14
	{
		if (tsvarArrayTemp[it] != NULL)
		delete3DArrayblock_Contiguous(tsvarArrayTemp[it]);
	}
	cout<<"process "<<rank<<"freeed tsvartemp"<<endl;
*/
	for (int inc = 0; inc < nncout; inc++)
	{
		for (int nindx = 0; nindx < (activeCells.size() / size + 1); nindx++)
			delete[] ncoutArray[inc][nindx];
		delete[] ncoutArray[inc];
	}
	delete[] ncoutArray;
	cout << "Process " << rank << " freed ncoutArray" << endl;
	for (int zk = 0; zk < nZones; zk++)
	{
		for (int ig = 0; ig < naggout; ig++)
			delete[] aggoutvarArray[zk][ig];
		delete[] aggoutvarArray[zk];
	}
	delete[] aggoutvarArray;
	cout<<"process "<<rank<<" freed aggregated output array"<<endl;
	/*for(int k=0 ;k<numOut; k++)
		delete[] OutVarValues[k];
	delete []OutVarValues; */
	//paramSite_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	cout << "Process " << rank << " finished" << endl;
	fflush(stdout);
	//MPI::COMM_WORLD.Barrier();
	if (rank == 0)
	{  
		//endTime = clock();
		TotalTime = MPI::Wtime() - startTimeT;                   //
		cout << "Time in seconds" << endl;
		cout << "Reading param  site state input control:  " << paramSite_Time << endl;
		cout << "Reading input TS txt arrays:  " << TsReadTime << endl;
		cout << "Reading input total TS arrays:  " << inputTS_Time << endl;
		cout << "Model simulation run time:  " << computeRun_Time << endl;
		cout << "Host<-->Device data copy time: " << dataCopy_Time << endl;
		cout << "Outptus write time: " << outputWrite_Time << endl;
		cout << "Total time of including overhead :  " << TotalTime << endl;
		cout << "Done! return value: " << retvalue << endl;
		//fflush(stdout);
	}
	//cout << "Done! return value: " << retvalue << endl;
	//fflush(stdout);
exitlab:
	MPI::Finalize();
	//getchar();
	return 0;
}
