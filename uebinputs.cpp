//#include "uebpg.h"
#include "uebpgdecls.h"
#pragma warning(disable : 4996)

//overload functions to read params
//this one would read params as array of floats
void readParams(const char* inpFile, float* &parArray, const int nParams)
{
	parArray = new float[nParams];
	ifstream pinFile(inpFile);		
	char headerLine[256];
	pinFile.getline(headerLine,256,'\n'); //skip header
	for (int i=0;i<nParams; i++)
	{	
		pinFile.getline(headerLine,256,'\n');
		pinFile.getline(headerLine,256,'\n');
		 sscanf(headerLine,"%f ",&parArray[i]);	
	}

	pinFile.close();
	return;
}

//function to read parameters into array of param struct with param name value pair
//Brute force! have to get better ways
void readParams(const char* inpFile, params strParamValues)
{
	FILE* pinFile = fopen(inpFile,"rt");	
	char headerLine[256];
//	char* parID;
	fgets(headerLine,256,pinFile); //skip header line
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.irad);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.ireadalb);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.cg);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.z);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.z);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.rhog);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.lc);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.ks);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.de);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.avo);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.anir0);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.lans);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.lang);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.wlf);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.rd1);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.dnews);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.emc);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.alpha);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.alpha);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.gpar);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.uc);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.as);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.Bs);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.lambda);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.rimax);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.wcoeff);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.apar);
	fgets(headerLine,256,pinFile); //skip id line
	fscanf(pinFile,"%f\n",&strParamValues.cpar);	     
	fclose(pinFile);
	return;
}
void readSiteVars(const char* inpFile, sitevar *&svArr)
{
	ifstream pinFile(inpFile);	
	char headerLine[256];
	//istringstream valueLine;	
	pinFile.getline(headerLine,256);   //skip header	
	float vardefaults[32] = { 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0, 
           0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903, 
           4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936,
		   9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 };

	for(int i=0; i<32;i++)
	{
		pinFile.getline(headerLine,256,':');
		sscanf(headerLine,"%s ",&svArr[i].svName);   
		pinFile.getline(headerLine,256,'\n'); 
		pinFile.getline(headerLine,256,'\n');
		sscanf(headerLine,"%d ",&svArr[i].svType);  		
			//headerLine[0] = 0;
			//fscanf(pinFile,"%d\n",&svArr[i].svType);
			switch (svArr[i].svType)
			{
			case -1: 
				svArr[i].svdefValue = vardefaults[i];
				break;
			case 0:
				pinFile.getline(headerLine,256,'\n');
		        sscanf(headerLine,"%f ",&svArr[i].svdefValue);
				break;
			case 1:
				pinFile.getline(headerLine,256,'\n');
		        sscanf(headerLine,"%s %s ",&svArr[i].svFile, &svArr[i].svVarName);
				break;			
			default:
				cout<<"Wrong site variable type; has to be -1 (Use default values), 0 (single value) or 1 (2D netcdf)"<<endl;
				cout<<"Using default value..."<<endl;
				svArr[i].svdefValue = vardefaults[i]; 
			}
			//i++;
			//headerLine[0] = 0;
		//}
	}
	pinFile.close();

	return;	
}


//function to read site variables (and initial contitions??)
//#_15 ics should come separate?
void readSiteVars(const char* inpFile, float svSValue[], char* svFile[], char* svVarName[], int svType[] )
{
	FILE* pinFile = fopen(inpFile,"rt");	
	char headerLine[256];
	fgets(headerLine,256,pinFile); //skip header line
	int i = 0; //, svtest;
	while (!feof(pinFile)) 
	{
		fgets(headerLine,256,pinFile); 
		if(headerLine[0] != ' ')
		{
			sscanf(headerLine,"%[^':'] %s",&svVarName[i]);
			fgets(headerLine,256,pinFile); 
			sscanf(headerLine,"%d ",&svType[i]);
			//fscanf(pinFile,"%d\n",&svArr[i].svType);
			switch (svType[i])
			{
			case -1: 
				//svArr[i].svdefValue = vardefaults[i];
				break;
			case 0:
				fgets(headerLine,256,pinFile);
				sscanf(headerLine,"%f ",&svSValue[i]);
				break;
			case 1:
				fgets(headerLine,256,pinFile);
				sscanf(headerLine,"%s %s",&svFile[i],&svVarName[i]);
				break;			
			default:
				cout<<"Wrong site variable type; has to be -1 (Use default values), 0 (single value) or 1 (2D netcdf)"<<endl;
				cout<<"Using default value..."<<endl;
				break;       //exit 
			}
			i++; //increment i
			headerLine[0] = ' ';
		}
	}
	fclose(pinFile);
	return;
	
}
void readSiteVars(const char* inpFile, sitevar svArr[], int indx)
{
	FILE* pinFile = fopen(inpFile,"r");	
	char headerLine[256];
	fgets(headerLine,256,pinFile); //skip header line
	
	float vardefaults[32] = { 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0, 
           0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903, 
           4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936,
		   9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 };

	for(int i=0; i<32;i++)
	{
		fgets(headerLine,256,pinFile); 
		if (headerLine[0] != ' ')                                  //while(fgets(headerLine,256,pinFile) != NULL) //
		{
			sscanf(headerLine,"%[^':'] %s",&svArr[i].svName);
			fgets(headerLine,256,pinFile); 
			sscanf(headerLine,"%d ",&svArr[i].svType);
			headerLine[0] = 0;
			//fscanf(pinFile,"%d\n",&svArr[i].svType);
			switch (svArr[i].svType)
			{
			case -1: 
				svArr[i].svdefValue = vardefaults[i];
				break;
			case 0:
				fgets(headerLine,256,pinFile);
				sscanf(headerLine,"%f ",&svArr[i].svdefValue);
				break;
			case 1:
				fgets(headerLine,256,pinFile);
				sscanf(headerLine,"%s %s",&svArr[i].svFile,&svArr[i].svVarName);
				break;			
			default:
				cout<<"Wrong site variable type; has to be -1 (Use default values), 0 (single value) or 1 (2D netcdf)\n"<<endl;
				cout<<"Using default value..."<<endl;
				svArr[i].svdefValue = vardefaults[i]; 
				//exit;
			}
			//i++;
			headerLine[0] = 0;
		}
	}
	fclose(pinFile);
	return;
	
}

//function  to read forcing/weather variables
void readInputForcVars(const char* inputconFile, inpforcvar *frArr)
{
	ifstream pinFile(inputconFile);	
	char headerLine[256];
	//istringstream valueLine;	
	pinFile.getline(headerLine,256);   //skip header	
	for(int i=0; i<13;i++)
	{
		pinFile.getline(headerLine,256,':');
		sscanf(headerLine,"%s ",&frArr[i].infName);   
		pinFile.getline(headerLine,256,'\n'); 
		pinFile.getline(headerLine,256,'\n');
		sscanf(headerLine,"%d ",&frArr[i].infType);  		
			//headerLine[0] = 0;
			//fscanf(pinFile,"%d\n",&svArr[i].svType);
			switch (frArr[i].infType)
			{
			case -1: 				
				pinFile.getline(headerLine, 256, '\n');
				sscanf(headerLine, "%f ", &frArr[i].infdefValue);
				break;
			case 0:
				pinFile.getline(headerLine,256,'\n');
		        sscanf(headerLine,"%s ",&frArr[i].infFile);
				break;
			case 1:
				pinFile.getline(headerLine,256,'\n');
				sscanf(headerLine,"%s %s %s %d",&frArr[i].infFile,&frArr[i].infvarName,&frArr[i].inftimeVar,&frArr[i].numNcfiles);
				break;	
			case 2:
				pinFile.getline(headerLine,256,'\n');
		        sscanf(headerLine,"%f ",&frArr[i].infdefValue);
				break;
			default:
				cout<<"Wrong input/forcing type; has to be -1 (compute by the model), 2 (single value) , 0 (time-series text file) or 1 (3D netcdf)"<<endl;
				cout<<"Using default value..."<<endl;
				break; //exit(1); 
			}
			//i++;
			//headerLine[0] = 0;
		//}
	}
	pinFile.close();	
			
	return;
}

//output control file: details of point, distributed netcdf, aggregated netcdf outputs 
void readOutputControl(const char* outputconFile, const char* aggoutputconFile, pointOutput* &pOut, ncOutput* &ncOut, aggOutput* &aggOut, int &npout, int &nncout, int &naggOut)
{
	ifstream poutFile(outputconFile);
	char headerLine[256];
	int nout;
	//istringstream valueLine;	
	poutFile.getline(headerLine, 256);   //skip header	
	//point outputs
	poutFile.getline(headerLine, 256);
	sscanf(headerLine, "%d ", &nout);
	npout = nout;
	pOut = new pointOutput[nout];
	for (int i = 0; i < nout; i++)
	{
		poutFile.getline(headerLine, 256);
		sscanf(headerLine, "%d %d %s ", &pOut[i].ycoord, &pOut[i].xcoord, &pOut[i].outfName);
	}
	//distributed netcdf output
	poutFile.getline(headerLine, 256);
	sscanf(headerLine, "%d ", &nout);
	nncout = nout;
	ncOut = new ncOutput[nout];
	for (int i = 0; i < nout; i++)
	{
		poutFile.getline(headerLine, 256);
		sscanf(headerLine, "%s %s %s ", &ncOut[i].symbol, &ncOut[i].outfName, &ncOut[i].units);
	}	
	poutFile.close();
	//aggregated outputs
	ifstream paoutFile(aggoutputconFile);
	paoutFile.getline(headerLine, 256);   //skip header	
	paoutFile.getline(headerLine, 256);
	sscanf(headerLine, "%d ", &nout);
	naggOut = nout;
	aggOut = new aggOutput[nout];
	for (int i = 0; i < nout; i++)
	{
		paoutFile.getline(headerLine, 256);
		sscanf(headerLine, "%s %s %s ", &aggOut[i].symbol, &aggOut[i].units, &aggOut[i].aggop);
	}	
	paoutFile.close();

	return;
}

//function to read input forcing time series text file
void readTextData(const char* inforcFile, float *&tcor_var, float *&tvar_in, int &nrecords)
{
	FILE* inputFile = fopen(inforcFile,"r");
	nrecords = 0;
	char commentLine[256];                    //string to read header line

	while(!feof(inputFile))
	{
		commentLine[0] = ' ';
		fgets(commentLine,256,inputFile);
		if (commentLine[0] != ' ')  //condition to make sure empty line is not read; 
			++nrecords;                                        
	}//while(!feof(inputfile)) 

	//strinpts = new inptimeseries[nrecords];                //assign memory to store data records
	tcor_var = new float[nrecords];
	tvar_in = new float[nrecords];

	rewind(inputFile);

	fgets(commentLine,256,inputFile);                     //skip the comment lines              
	for (int i=0; i<nrecords; i++)   
		//#_15 check format of dtime here
			fscanf(inputFile,"%f %f \n",&tcor_var[i],&tvar_in[i]);		

	fclose(inputFile);

}//void readTextData

// overload function with only the variable reading (input forcing time series text file)
void readTextData(const char* inforcFile, float *&tvar_in, int &nrecords)
{
	ifstream inputFile(inforcFile,ios::in);
	nrecords = 0;
	char commentLine[256];                    //string to read header line
	if(!inputFile)
	{
		cout<<"Error opening file: "<<inforcFile<<endl;
		return;
	}
	inputFile.getline(commentLine,256,'\n');  //skip first header line 

	//inputFile.getline(commentLine,256,'\n'); //first data line
	while(!inputFile.eof())
	{
		commentLine[0] = ' ';
		inputFile.getline(commentLine,256,'\n');
		if (commentLine[0] != ' ')  //condition to make sure empty line is not read; 
			++nrecords;                                        
	}//while

	//cout<<"number of records in file: "<<inforcFile<<" "<<nrecords<<endl;
	tvar_in = new float[nrecords];

	inputFile.seekg(0L,ios::beg);  

	inputFile.getline(commentLine,256,'\n');  //skip first header line 
	for (int i=0; i<nrecords-1; i++) 
	{
		inputFile.getline(commentLine,256,'\n');
		sscanf(commentLine,"%*d %*d %*d %*d %f \n",&tvar_in[i]); //&tcor_var[i],&tvar_in[i]);	
	}

	inputFile.close();

}//void readTextData





