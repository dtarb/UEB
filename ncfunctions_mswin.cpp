//#include"nctest.h"
#include "uebpgdecls.h"
//writes the 1D array (TS values) at specified location in the netcdf
int WriteTSto3DNC(const char* FileName, const char* VarName, int dimOrd, int y_dim, int x_dim,int Nt_dim,  float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo)  //ydim, xdim =the coordinate point the data to be written; Nt_dim =the length of the TS
{  
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, y_dimid = 0, x_dimid = 0, t_dimid = 0; 
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0; 
	const char* vunits = "UNITS";
	const int  NDIMS = 3;
	int dimids[NDIMS]; 	
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];  		  
	// Error handling.  
	int retncval = 0;  
		
	// Open netcdf file.  
	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))  
		ERR(retncval);   
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VarName, &v_varid))) 
		ERR(retncval);
	
	/* The dimids array is used to pass the dimids of the dimensions of  
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */  
	 switch(dimOrd)
	 {
	    case 0:        //t,y,x
			 start[0] = 0;  
			 start[1] = y_dim; 
			 start[2] = x_dim;
			 count[0] = Nt_dim; 	
			 count[1] = 1; 
			 count[2] = 1;			 
			 break;
		case 1:   //y,x,t
			 start[0] = y_dim; 
			 start[1] = x_dim;
			 start[2] = 0; 
			 count[0] = 1; 
			 count[1] = 1;
			 count[2] = Nt_dim; 	
			 break;
		case 2:           //x,t,y
			 
			start[0] = x_dim;
			start[0] = 0;  
			start[2] = y_dim; 	
			count[0] = 1; 
			count[1] = Nt_dim; 				
			count[2] = 1;			 
			 break;		
		default:
			cout<<"the dim order has to be between 0 and 2"<<endl;
			getchar();
			break;
	 }	 
	//put variable values
	if ((retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0]))) 
				 ERR(retncval); 		 	
	//close file
	if ((retncval = nc_close(ncid)))    
		ERR(retncval);  
	//delte 3D array
	
	
	//fflush(stdout); 
	return 0;
}
//writes multiple 1D arrays (TS values) at specified locations in the netcdf
int WriteTSto3DNC_Block(const char* FileName, const char* VarName, int dimOrd, int *YindArr, int *XindArr, int bSize, int Nt_dim, float** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo)  //ydim, xdim =the coordinate point the data to be written; Nt_dim =the length of the TS
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, y_dimid = 0, x_dimid = 0, t_dimid = 0;
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0;
	const char* vunits = "UNITS";
	const int  NDIMS = 3;
	int dimids[NDIMS];
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Error handling.  
	int retncval = 0;

	// Open netcdf file.  
	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))             //| NC_MPIIO, inpComm, inpInfo
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VarName, &v_varid)))
		ERR(retncval);
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
	for (int j = 0; j < bSize; j++)
	{
		switch (dimOrd)
		{
		case 0:        //t,y,x
			start[0] = 0;
			start[1] = YindArr[j];   // y_dim;
			start[2] = XindArr[j];   // x_dim;
			count[0] = Nt_dim;
			count[1] = 1;
			count[2] = 1;
			break;
		case 1:   //y,x,t
			start[0] = YindArr[j];   // y_dim;
			start[1] = XindArr[j];	 // x_dim;
			start[2] = 0;
			count[0] = 1;
			count[1] = 1;
			count[2] = Nt_dim;
			break;
		case 2:           //x,t,y

			start[0] = XindArr[j];   // x_dim;
			start[0] = 0;
			start[2] = YindArr[j];   // y_dim;
			count[0] = 1;
			count[1] = Nt_dim;
			count[2] = 1;
			break;
		default:
			cout << "the dim order has to be between 0 and 2" << endl;
			getchar();
			break;
		}		
		//put variable values
		if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[j][0]))
			ERR(retncval);
	}
	//close file
	if ((retncval = nc_close(ncid)))
		ERR(retncval);
	//delte 3D array
	//fflush(stdout); 
	return 0;
}
//aggregate outputs at specified zone in the netcdf
int Write_uebaggTS_toNC(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp)  
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, y_dimid = 0, x_dimid = 0, t_dimid = 0;
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0;
	const char* vunits = "UNITS";
	const int  NDIMS = 3;
	int dimids[NDIMS];
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Error handling.  
	int retncval = 0;

	// Open netcdf file.  
	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VarName, &v_varid)))
		ERR(retncval);

	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
	switch (dimOrd)
	{
	case 0:        //t,y,x
		start[0] = 0;
		start[1] = z_dim;		
		count[0] = Nt_dim;
		count[1] = 1;		
		break;
	case 1:   //y,x,t
		start[0] = z_dim;
		start[1] = 0;		
		count[0] = 1;
		count[1] = Nt_dim;
		break;	
	default:
		cout << "the dim order has to be 0 or 1" << endl;
		getchar();
		break;
	}
	//put variable values
	if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0]))
		ERR(retncval);
	//close file
	if ((retncval = nc_close(ncid)))
		ERR(retncval);	
	//fflush(stdout); 
	return 0;
}
//aggregate outputs at specified zone in the netcdf in parallel 
int Write_uebaggTS_toNC_par(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, y_dimid = 0, x_dimid = 0, t_dimid = 0;
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0;
	const char* vunits = "UNITS";
	const int  NDIMS = 3;
	int dimids[NDIMS];
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Error handling.  
	int retncval = 0;
	// Open netcdf file.  
	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))           //|NC_MPIIO, inpComm, inpInfo,, &ncid)))
		ERR(retncval);
	// get variable id
	if (retncval = nc_inq_varid(ncid, VarName, &v_varid))
		ERR(retncval);
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
	switch (dimOrd)
	{
	case 0:        //t,y,x
		start[0] = 0;
		start[1] = z_dim;
		count[0] = Nt_dim;
		count[1] = 1;
		break;
	case 1:   //y,x,t
		start[0] = z_dim;
		start[1] = 0;
		count[0] = 1;
		count[1] = Nt_dim;
		break;
	default:
		cout << "the dim order has to be 0 or 1" << endl;
		getchar();
		break;
	}

	//put variable values
	if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0]))
		ERR(retncval);
	//close file
	if ((retncval = nc_close(ncid)))
		ERR(retncval);
	//fflush(stdout); 
	return 0;
}
//creates 3D netcdf in PARALLEL mode, and stores dimension variables for UEB aggregated outputs; some attributes are copied from the watershed netCDF
int create3DNC_uebAggregatedOutputs(const char* FileName, aggOutput *aggOut, int naggOut, const char* tName, const char* tUnits, const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd,
	float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, int nZones, const char * zName, float* y_inp, float* x_inp,  MPI::Intracomm inpComm, MPI::Info inpInfo)
{	
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, pxid = 0, pyid = 0;
	//info from watershed nc
	if ((retncval = nc_open(ws_FileName, NC_NOWRITE, &ncid)))
		ERR(retncval);	
	if ((retncval = nc_inq_varid(ncid, ws_VarName, &pvarid)))
		ERR(retncval);
	// Get the varids of the coordinate variables 	
	if ((retncval = nc_inq_varid(ncid, yName, &pyid)))
		ERR(retncval);
	if ((retncval = nc_inq_varid(ncid, xName, &pxid)))
		ERR(retncval);
	//dimensions lengths
	size_t Nzdim = nZones;
	//3D outputs nc	
	int ncid_out = 0, pzid_out = 0, ptid_out = 0;
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0, gridmap_id = 0;
	//attributes
	int natts = 0;
	size_t att_len;
	void * attValue;
	nc_type attType;
	char attName[256]; // = {}; // "grid_mapping_or_any_other_long_enough_name";
	char grid_mappingValue[256] = {NULL}; // "grid_mapping_or_any_other_long_enough_name";	

	const char* fillValueName = "_FillValue";
	//float missVal = -9999;
	const int  NDIMS = 2;
	int dimids[NDIMS];
	int oldFill = 0;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid_out)))
		ERR(retncval);
	//?? set fill on
	if ((retncval = nc_set_fill(ncid_out, NC_FILL, &oldFill)))
		ERR(retncval);
	//copy global attributes from the ws nc
	if (retncval = nc_inq_natts(ncid, &natts))
		ERR(retncval);
	for (int atti = 0; atti < natts; atti++)
	{
		if (retncval = nc_inq_attname(ncid, NC_GLOBAL, atti, attName))
			ERR(retncval);		
		if (retncval = nc_copy_att(ncid, NC_GLOBAL, (const char*)attName, ncid_out, NC_GLOBAL))
			ERR(retncval);
	}
	/* Define the dimensions. record dim can be unlimited*/
	if (retncval = nc_def_dim(ncid_out, tName, Nt_dim, &ptid_out))
		ERR(retncval);
	if (retncval = nc_def_dim(ncid_out, zName, Nzdim, &pzid_out))
		ERR(retncval);
	switch (dimOrd)
	{
	case 0:
		dimids[0] = ptid_out;
		dimids[1] = pzid_out;		
		break;
	case 1:
		dimids[0] = pzid_out;
		dimids[1] = ptid_out;		
		break;	
	default:
		cout << "the dim order has to be either 0 or 1" << endl;
		getchar();
		break;
	}	
	// Define variables	
	for (int ivar = 0; ivar < naggOut; ivar++)
	{		
		if (retncval = nc_def_var(ncid_out, (const char*) aggOut[ivar].symbol, NC_FLOAT, NDIMS, dimids, &v_varid))
			ERR(retncval);
		//assign fill value
		if (retncval = nc_put_att_float(ncid_out, v_varid, fillValueName, NC_FLOAT, 1, fillVal))
			ERR(retncval);
		// Assign units attributes to the netCDF variables.  
		if (retncval = nc_put_att_text(ncid_out, v_varid, "units", strlen((const char*)aggOut[ivar].units), (const char*)aggOut[ivar].units))
			ERR(retncval);	
		//grid mapping (projection) attribute
		if (retncval = nc_copy_att(ncid, pvarid, "grid_mapping", ncid_out, v_varid))
			ERR(retncval);
	}
	//variable containing info about projection / grid mapping	
	if (retncval = nc_get_att_text(ncid, pvarid, "grid_mapping", grid_mappingValue))
		ERR(retncval);
	if (retncval = nc_inq_varid(ncid, (const char*)grid_mappingValue, &gridmap_id))
		ERR(retncval);
	if (retncval = nc_copy_var(ncid, gridmap_id, ncid_out))
		ERR(retncval);
	//time variable
	if (retncval = nc_def_var(ncid_out, tName, NC_FLOAT, 1, &ptid_out, &t_varid))
		ERR(retncval); 
	if (retncval = nc_put_att_text(ncid_out, t_varid, "units", strlen(tUnits), tUnits))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid_out, t_varid, "long_name", strlen(tlong_name), tlong_name))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid_out, t_varid, "calendar", strlen(tcalendar), tcalendar))
		ERR(retncval);
	// y and x vars
	if (retncval = nc_def_var(ncid_out, yName, NC_FLOAT, 1, &pzid_out, &y_varid))
		ERR(retncval);
	if (retncval = nc_def_var(ncid_out, xName, NC_FLOAT, 1, &pzid_out, &x_varid))
		ERR(retncval);
	if (retncval = nc_inq_varnatts(ncid, pyid, &natts))
		ERR(retncval);
	for (int atti = 0; atti < natts; atti++)
	{
		if (retncval = nc_inq_attname(ncid, pyid, atti, attName))
			ERR(retncval);
		if (retncval = nc_copy_att(ncid, pyid, (const char*)attName, ncid_out, y_varid))
			ERR(retncval);
	}
	if (retncval = nc_inq_varnatts(ncid, pxid, &natts))
		ERR(retncval);
	for (int atti = 0; atti < natts; atti++)
	{
		if (retncval = nc_inq_attname(ncid, pxid, atti, attName))
			ERR(retncval);
		if (retncval = nc_copy_att(ncid, pxid, (const char*)attName, ncid_out, x_varid))
			ERR(retncval);
	}
	//put values to dim variables
	if (retncval = nc_put_var_float(ncid_out, t_varid, &t_inp[0]))
		ERR(retncval);
	if (retncval = nc_put_var_float(ncid_out, y_varid, &y_inp[0]))
		ERR(retncval);
	if (retncval = nc_put_var_float(ncid_out, x_varid, &x_inp[0]))
		ERR(retncval);
	//close files
	if (retncval = nc_close(ncid))
		ERR(retncval);
	if (retncval = nc_close(ncid_out))
		ERR(retncval);
	cout << "Sucess creating and storing dimension vars in: " << FileName << endl;
	//fflush(stdout); 
	return 0;
}
//creates 3D netcdf and stores dimension variables for a given UEB output; called once for a given output netcdf 
//attributes are copied from the 2D (watershed) netCDF which the 3D netCDF spatial dimension follow
int create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnits, const char* tlong_name,
	const char* tcalendar, int Nt_dim, int dimOrd, float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, MPI::Intracomm inpComm, MPI::Info inpInfo)
{	
	//const char* yUnits, const char* xUnits;
	float* y_inp;
	float* x_inp;
		//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, pxid = 0, pyid = 0, ndims = 0;
	//dimensions lengths
	size_t Nxdim = 0, Nydim = 0;		
	//array of dimensions
	int pdimids[2]; 
	if ((retncval = nc_open(ws_FileName, NC_NOWRITE, &ncid)))
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, ws_VarName, &pvarid)))
		ERR(retncval);
	// Get the varids of the coordinate variables 	
	if ((retncval = nc_inq_varid(ncid, yName, &pyid)))
		ERR(retncval);
	if ((retncval = nc_inq_varid(ncid, xName, &pxid)))
		ERR(retncval);
	//var information, checking the dimension array
	if (retncval = nc_inq_vardimid(ncid, pvarid, pdimids))
		ERR(retncval);	
	//check dimension sizes
	if ((retncval = nc_inq_dim(ncid, pdimids[0], NULL, &Nydim)))
		ERR(retncval);
	if ((retncval = nc_inq_dim(ncid, pdimids[1], NULL, &Nxdim)))
		ERR(retncval);

	//create the output arrays	
	y_inp = new float[Nydim];
	x_inp = new float[Nxdim];
	// Read the coordinate (dimensions) variable data.  
	if ((retncval = nc_get_var_float(ncid, pyid, &y_inp[0])))
		ERR(retncval);
	if ((retncval = nc_get_var_float(ncid, pxid, &x_inp[0])))
		ERR(retncval);	

	//3D outputs nc	
	int ncid_out = 0, pyid_out = 0, pxid_out = 0, ptid_out = 0;
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0, gridmap_id = 0;
	//attributes
	int natts = 0;
	size_t att_len;
	void * attValue;
	nc_type attType;	
	char attName[256]; // = {}; // "grid_mapping_or_any_other_long_enough_name";
	char grid_mappingValue[256] = { NULL }; // "grid_mapping_or_any_other_long_enough_name";	
	

	const char* fillValueName = "_FillValue";
	//float missVal = -9999;
	const int  NDIMS = 3;
	int dimids[NDIMS];
	int oldFill = 0;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid_out)))
		ERR(retncval);
	//?? set fill on
	if ((retncval = nc_set_fill(ncid_out, NC_FILL, &oldFill)))
		ERR(retncval);
	//copy global attributes from the ws nc
	if (retncval = nc_inq_natts(ncid, &natts))
		ERR(retncval);
	for (int atti = 0; atti < natts; atti++)
	{
		if (retncval = nc_inq_attname(ncid,NC_GLOBAL,atti, attName))
			ERR(retncval);
		/*if (retncval = nc_inq_atttype(ncid, NC_GLOBAL,(const char*)attName,&attType))
			ERR(retncval);
		if (retncval = nc_inq_attlen(ncid, NC_GLOBAL, (const char*)attName, &att_len))
			ERR(retncval);
		attValue = ::operator new (att_len*sizeof(attType));
		if (retncval = nc_get_att(ncid, NC_GLOBAL, (const char*)attName,&attValue))
			ERR(retncval);*/
		if (retncval = nc_copy_att(ncid, NC_GLOBAL, (const char*)attName,ncid_out,NC_GLOBAL))
			ERR(retncval);
	}
	/* Define the dimensions. record dim can be unlimited*/
	if (retncval = nc_def_dim(ncid_out, tName, Nt_dim, &ptid_out))
		ERR(retncval);
	if (retncval = nc_def_dim(ncid_out, yName, Nydim, &pyid_out))
		ERR(retncval);
	if (retncval = nc_def_dim(ncid_out, xName, Nxdim, &pxid_out))
		ERR(retncval);
	
	switch (dimOrd)
	{
	case 0:
		dimids[0] = ptid_out;
		dimids[1] = pyid_out;
		dimids[2] = pxid_out;
		break;
	case 1:
		dimids[0] = pyid_out;
		dimids[1] = pxid_out;
		dimids[2] = ptid_out;
		break;
	case 2:
		dimids[0] = pxid_out;
		dimids[1] = ptid_out;
		dimids[2] = pyid_out;
		break;
	default:
		cout << "the dim order has to be between 0 and 2" << endl;
		getchar();
		break;
	}
	// Define the netCDF variables 
	if (retncval = nc_def_var(ncid_out, VarName, NC_FLOAT, NDIMS, dimids, &v_varid))
		ERR(retncval);
	//assign fill value
	if (retncval = nc_put_att_float(ncid_out, v_varid, fillValueName, NC_FLOAT, 1, fillVal))
		ERR(retncval);
	// Assign units attributes to the netCDF variables.  
	if (retncval = nc_put_att_text(ncid_out, v_varid, "units", strlen(varUnits), varUnits))
		ERR(retncval);	
	//grid mapping
	if (retncval = nc_get_att_text(ncid, pvarid, "grid_mapping",grid_mappingValue))
		ERR(retncval);
	//cout << grid_mappingValue << endl;
	if (retncval = nc_copy_att(ncid, pvarid, "grid_mapping", ncid_out, v_varid))
		ERR(retncval);
	/*if (retncval = nc_put_att_text(ncid_out, v_varid, "grid_mapping", , (const char*)grid_mappingValue))
		ERR(retncval);*/
	if (retncval = nc_inq_varid(ncid,(const char*)grid_mappingValue, &gridmap_id))
		ERR(retncval);
	if (retncval = nc_copy_var(ncid, gridmap_id, ncid_out))
		ERR(retncval);
	//time variable
	if (retncval = nc_def_var(ncid_out, tName, NC_FLOAT, 1, &ptid_out, &t_varid))
		ERR(retncval);
	// units attributes  
	if (retncval = nc_put_att_text(ncid_out, t_varid, "units", strlen(tUnits), tUnits))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid_out, t_varid, "long_name", strlen(tlong_name), tlong_name))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid_out, t_varid, "calendar", strlen(tcalendar), tcalendar))
		ERR(retncval);
	// y and x vars
	if (retncval = nc_def_var(ncid_out, yName, NC_FLOAT, 1, &pyid_out, &y_varid))
		ERR(retncval);	
	if (retncval = nc_def_var(ncid_out, xName, NC_FLOAT, 1, &pxid_out, &x_varid))
		ERR(retncval);	
	if (retncval = nc_inq_varnatts(ncid, pyid, &natts))
		ERR(retncval);
	for (int atti = 0; atti < natts; atti++)
	{
		if (retncval = nc_inq_attname(ncid, pyid, atti, attName))
			ERR(retncval);
		if (retncval = nc_copy_att(ncid, pyid, (const char*)attName, ncid_out, y_varid))
			ERR(retncval);
	}
	if (retncval = nc_inq_varnatts(ncid, pxid, &natts))
		ERR(retncval);
	for (int atti = 0; atti < natts; atti++)
	{
		if (retncval = nc_inq_attname(ncid, pxid, atti, attName))
			ERR(retncval);
		if (retncval = nc_copy_att(ncid, pxid, (const char*)attName, ncid_out, x_varid))
			ERR(retncval);
	}
	/* Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, y_varid, vunits, strlen(yUnits), yUnits)))
		ERR(retncval);
	if ((retncval = nc_put_att_text(ncid, x_varid, vunits, strlen(xUnits), xUnits)))
		ERR(retncval);
	//y and x variables copy
	if (retncval = nc_copy_var(ncid, pyid, ncid_out))
		ERR(retncval);
	if (retncval = nc_copy_var(ncid, pxid, ncid_out))
		ERR(retncval);	*/

	//put values to dim variables
	if ((retncval = nc_put_var_float(ncid_out, t_varid, &t_inp[0])))
		ERR(retncval);
	if ((retncval = nc_put_var_float(ncid_out, y_varid, &y_inp[0])))
		ERR(retncval);
	if ((retncval = nc_put_var_float(ncid_out, x_varid, &x_inp[0])))
		ERR(retncval);
	//close files
	if (retncval = nc_close(ncid))
		ERR(retncval);
	if (retncval = nc_close(ncid_out))
		ERR(retncval);

	cout << "Sucess creating and storing dimension vars in: " << FileName << endl;
	//fflush(stdout); 
	return 0;
}

//creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
int Create3DNC(const char* FileName, const char* VarName, const char *varUnits,  const char* tName,  const char* yName, const char* xName, const char* tUnits,
	const char* yUnits, const char* xUnits, int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float* fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo)
{  
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, y_dimid = 0, x_dimid = 0, t_dimid = 0; 
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0; 
	const char* vunits = "UNITS";
	const char* fillValue = "_FillValue";
	//float missVal = -9999;
	const int  NDIMS = 3;
	int dimids[NDIMS]; 	
	int oldFill = 0;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];  		  
	// Error handling.  
	int retncval = 0;  
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4|NC_CLOBBER, &ncid)))  
		ERR(retncval);   
	//set fill on
	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill))) 
		ERR(retncval); 
	/* Define the dimensions. record dim can be unlimited*/
	if ((retncval = nc_def_dim(ncid, tName, Nt_dim, &t_dimid))) 
		ERR(retncval); 
	if ((retncval = nc_def_dim(ncid, yName, Ny_dim, &y_dimid)))   
		ERR(retncval);  
	 if ((retncval = nc_def_dim(ncid, xName, Nx_dim, &x_dimid)))  
		ERR(retncval);  
	/* The dimids array is used to pass the dimids of the dimensions of  
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */  
	 switch(dimOrd)
	 {
	    case 0:
			 dimids[0] = t_dimid; 
			 dimids[1] = y_dimid; 
			 dimids[2] = x_dimid;
			 break;
		case 1:
			 dimids[0] = y_dimid; 
			 dimids[1] = x_dimid;
			 dimids[2] = t_dimid; 				
			 break;
		case 2:
			 dimids[0] = x_dimid;
			 dimids[1] = t_dimid; 
			 dimids[2] = y_dimid; 			 
			 break;		
		default:
			cout<<"the dim order has to be between 0 and 2"<<endl;
			getchar();
			break;
	 }
	 
	// Define the netCDF variables 
	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS, dimids, &v_varid)))   
		ERR(retncval);  
	//assign missing
	if ((retncval = nc_put_att_float(ncid,v_varid, fillValue,NC_FLOAT,1, fillVal)))    
		ERR(retncval);
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))    
		ERR(retncval); 
	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnits), tUnits)))    
		ERR(retncval);   
	if ((retncval = nc_def_var(ncid, yName, NC_FLOAT, 1, &y_dimid, &y_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, y_varid, vunits, strlen(yUnits), yUnits)))    
		ERR(retncval); 
	if ((retncval = nc_def_var(ncid, xName, NC_FLOAT, 1, &x_dimid, &x_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, x_varid, vunits,  strlen(xUnits), xUnits)))    
		ERR(retncval); 	
	// End define mode. not needed for nc4
	/*if ((retncval = nc_enddef(ncid)))   
		ERR(retncval); 	*/

	//put values to dim variables
	if ((retncval = nc_put_var_float(ncid, t_varid, &t_inp[0]))) 
				 ERR(retncval); 
	if ((retncval = nc_put_var_float(ncid, y_varid, &y_inp[0]))) 
				 ERR(retncval); 
	if ((retncval = nc_put_var_float(ncid, x_varid, &x_inp[0]))) 
				 ERR(retncval); 
	
	//put variable values
	/*if ((retncval = nc_put_var_float(ncid, v_varid, &var_inp[0][0][0]))) 
				 ERR(retncval); 		*/ 	
	//close file
	if ((retncval = nc_close(ncid)))    
		ERR(retncval);  
	//delte 3D array	
	cout<<"Sucess creating and storing dimension vars in: "<<FileName<<endl; 
	//fflush(stdout); 
	return 0;
}

//This one creates and stores array passed to it at once
int Write3DNC(const char* FileName, const char* VarName, const char *varUnits,  const char* tName,  const char* yName, const char* xName, const char* tUnits,  const char* yUnits, const char* xUnits, 
				int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float*** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo)// int &t_dimid, int &y_dimid, int &x_dimid) //time, y, x
{  
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, y_dimid = 0, x_dimid = 0, t_dimid = 0; 
	int v_varid = 0, t_varid = 0, y_varid = 0, x_varid = 0; 
	const char* vunits = "UNITS";
	const int  NDIMS = 3;
	int dimids[NDIMS]; 	
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];  		  
	// Error handling.  
	int retncval = 0;  
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid)))
		ERR(retncval);   
	/* Define the dimensions. record dim can be unlimited*/ 	 
	if ((retncval = nc_def_dim(ncid, tName, Nt_dim, &t_dimid))) 
		ERR(retncval); 
	if ((retncval = nc_def_dim(ncid, yName, Ny_dim, &y_dimid)))   
		ERR(retncval);  
	 if ((retncval = nc_def_dim(ncid, xName, Nx_dim, &x_dimid)))  
		ERR(retncval);  
	/* The dimids array is used to pass the dimids of the dimensions of  
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */  
	 switch(dimOrd)
	 {
	    case 0:
			 dimids[0] = t_dimid; 
			 dimids[1] = y_dimid; 
			 dimids[2] = x_dimid;
			 break;
		case 1:
			 dimids[0] = y_dimid; 
			 dimids[1] = x_dimid;
			 dimids[2] = t_dimid; 				
			 break;
		case 2:
			 dimids[0] = x_dimid;
			 dimids[1] = t_dimid; 
			 dimids[2] = y_dimid; 			 
			 break;		
		default:
			cout<<"the dim order has to be between 0 and 2"<<endl;
			getchar();
			break;
	 }
	 
	// Define the netCDF variables 
	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS,   
		dimids, &v_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, v_varid, vunits,              
		strlen(varUnits), varUnits)))    
		ERR(retncval); 
	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1,   
		&t_dimid, &t_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, t_varid, vunits,              
		strlen(tUnits), tUnits)))    
		ERR(retncval);   
	if ((retncval = nc_def_var(ncid, yName, NC_FLOAT, 1,   
		&y_dimid, &y_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, y_varid, vunits,              
		strlen(yUnits), yUnits)))    
		ERR(retncval); 
	if ((retncval = nc_def_var(ncid, xName, NC_FLOAT, 1,   
		&x_dimid, &x_varid)))   
		ERR(retncval);  
		// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, x_varid, vunits,              
		strlen(xUnits), xUnits)))    
		ERR(retncval); 	
	// End define mode.  
	if ((retncval = nc_enddef(ncid)))   
		ERR(retncval); 	
	//put values to dim variables
	if ((retncval = nc_put_var_float(ncid, t_varid, &t_inp[0]))) 
				 ERR(retncval); 
	if ((retncval = nc_put_var_float(ncid, y_varid, &y_inp[0]))) 
				 ERR(retncval); 
	if ((retncval = nc_put_var_float(ncid, x_varid, &x_inp[0]))) 
				 ERR(retncval); 
	
	//put variable values
	if ((retncval = nc_put_var_float(ncid, v_varid, &var_inp[0][0][0]))) 
				 ERR(retncval); 		 	
	//close file
	if ((retncval = nc_close(ncid)))    
		ERR(retncval);  
	//delte 3D array	
	cout<<"SUCCESS writing file: "<<FileName<<endl; 
	//fflush(stdout); 
	return 0;
}
//read 3d nc to contiguous array;  y,x,time coordinate names are same as index names
int read3DNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* tcorvar, float* ycorvar, float* xcorvar, size_t &Ntdim, size_t &Nydim, size_t &Nxdim, float*** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, pxid = 0, pyid = 0, ptid = 0, ndims = 0; 
	//dimensions lengths
	size_t pdim_sizes;
	//variable data type
	nc_type varType;
	//array of dimensions
	int pdimids[3]; //NC_MAX_DIMS]; 3D file only being read here; expected to get error message otherwise
	//dimension names 
	char pdim_Names[80];
	size_t start[3], count[3];
	//dimensions lengths
	//size_t Ntdim = 0, Nxdim = 0, Nydim = 0;
	//Open the file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))         //NC_MPIIO, inpComm, inpInfo, &ncid)))
		ERR(retncval); 
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid))) 
		ERR(retncval); 
	// Get the varids of the coordinate variables 
	if ((retncval = nc_inq_varid(ncid, tcor_NAME, &ptid))) 
		ERR(retncval);
	if ((retncval = nc_inq_varid(ncid, ycor_NAME, &pyid))) 
		ERR(retncval); 
	if ((retncval = nc_inq_varid(ncid, xcor_NAME, &pxid))) 
		ERR(retncval); 	
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid,pvarid,NULL,&varType,NULL,pdimids,NULL))) 
		ERR(retncval); 
	//check dimension info and set start and count arrays; 
	for (int i = 0; i < 3; i++)
	{
		if ((retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes)))
			ERR(retncval);
		if (strcmp(pdim_Names, tcor_NAME) == 0)
		{
			start[i] = 0;
			Ntdim = pdim_sizes;
			count[i] = Ntdim;
			tcorvar = new float[Ntdim];
		}
		else if (strcmp(pdim_Names, ycor_NAME) == 0)
		{
			start[i] = 0;
			Nydim = pdim_sizes;
			count[i] = Nydim;
			ycorvar = new float[Nydim];
		}
		else
		{
			start[i] = 0;
			Nxdim = pdim_sizes;
			count[i] = Nxdim;
			xcorvar = new float[Nxdim];
		}
	}
	//create the output arrays
	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);          // (Ntdim, Nydim, Nxdim);	

	/* Read the coordinate (dimensions) variable data. */ 
      if (retncval = nc_get_var_float(ncid, ptid, &tcorvar[0]))
		ERR(retncval);
	if (retncval = nc_get_var_float(ncid, pyid, &ycorvar[0]))
		ERR(retncval);
	if (retncval = nc_get_var_float(ncid, pxid, &xcorvar[0]))
		ERR(retncval);  
	 //read variable (input data)	
	if (retncval = nc_get_var_float(ncid, pvarid,&pvar_in[0][0][0]))    
		ERR(retncval); 
	//close netcdf file			
	if ((retncval = nc_close(ncid))) 
		ERR(retncval); 
    //cout<<"SUCCESS reading input file: "<< FILE_NAME<<endl;  
	return 0;
}

int read3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME, const char* ycor_NAME, const char* tcor_NAME, 
                    float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, pxid = 0, pyid = 0, ptid = 0, ndims = 0; 
	//dimensions lengths
	size_t Ntdim = 0, Nxdim = 0, Nydim = 0;
	//variable data type
	nc_type varType;
	//array of dimensions
	int pdimids[3]; //NC_MAX_DIMS]; 3D file only being read here; expected to get error messege otherwise
	size_t start[3], count[3];
	//Open the file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))		
		ERR(retncval); 
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid))) 
		ERR(retncval); 
	// Get the varids of the coordinate variables 
	if ((retncval = nc_inq_varid(ncid, tcor_NAME, &ptid))) 
		ERR(retncval);
	if ((retncval = nc_inq_varid(ncid, ycor_NAME, &pyid))) 
		ERR(retncval); 
	if ((retncval = nc_inq_varid(ncid, xcor_NAME, &pxid))) 
		ERR(retncval); 	
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid,pvarid,NULL,&varType,NULL,pdimids,NULL))) 
		ERR(retncval); 
	//check dimension sizes
	if ((retncval = nc_inq_dim(ncid,pdimids[0],NULL,&Ntdim))) 
		ERR(retncval); 
	if ((retncval = nc_inq_dim(ncid,pdimids[1],NULL,&Nydim))) 
		ERR(retncval); 
	if ((retncval = nc_inq_dim(ncid,pdimids[2],NULL,&Nxdim))) 
		ERR(retncval); 	
	//create the output arrays
	pvar_in = Create3DArray(Ntdim,Nydim,Nxdim);
	tcorvar = new float[Ntdim];
	ycorvar = new float[Nydim];
	xcorvar = new float[Nxdim];
	//set start and count arrays; time, if exists, is expected to to come first
	//although in this code the effect may not be noticed 
	//#_1 (need a check 5/19/13)
	start[0] = 0;   
	start[1] = 0;  
	start[2] = 0;  
	count[0] = 1;  
	count[1] = Nydim; 
	count[2] = Nxdim;
	//temporary 2D array to help in data read
	//#_2 (need check 5/19/13) should read 3d @ once 
	float ** pvar_inp = new float*[Nydim];
	for(size_t j=0; j< Nydim; j++)			
		pvar_inp[j] = new float [Nydim];		
	/* Read the coordinate (dimensions) variable data. */ 
	if ((retncval = nc_get_var_float(ncid, ptid, &tcorvar[0])))
		ERR(retncval);
	if ((retncval = nc_get_var_float(ncid, pyid, &ycorvar[0])))
		ERR(retncval);
	if ((retncval = nc_get_var_float(ncid, pxid, &xcorvar[0]))) 
		ERR(retncval);  
	 //read variable (input data)
	//the readings below are done as Double data type; expected type conversion 
	//to float if input is float (no data loss than the vice versa)
	//#_3 tbc 3.19.13
	for (size_t kt = 0; kt < Ntdim; kt++)  
	{
		start[0] = kt;
		if ((retncval = nc_get_vara_float(ncid, pvarid,start,count, &pvar_inp[0][0])))      
			ERR(retncval); 
		for(size_t i=0; i< Nydim; i++)
			for(size_t j=0;j<Nxdim;j++)
				pvar_in[kt][i][j] = pvar_inp[i][j];
		//cout<<"step no %d\n",k);
	}	/* next record */ 
	//free temporary matrix pvar_inp
	for(size_t i=0; i< Nydim; i++)
		delete []pvar_inp[i];
	delete []pvar_inp;
	//close netcdf file			
	if ((retncval = nc_close(ncid))) 
		ERR(retncval); 
	cout<<"SUCCESS reading input file: "<<FILE_NAME<<endl;  
	return 0;
}
//function to read multiple blocks of single column/rod along time dimension from 3D netcdf file, for given y , x coordinate arrays
int readNC_TS_Block(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float** &pvar_in,  int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo, int *YindArr, int *XindArr, int bSize) // /*float* &tcorvar, int ydim, int xdim, */ int tIndx) 
{
	//float* pvarin_temp = NULL;
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, ptid = 0; // pxid = 0, pyid = 0, ndims = 0; 
	//dimensions lengths
	size_t pdim_sizes; //[3]; //Ntdim = 0, Nxdim = 0, Nydim = 0;
	//variable data type
	nc_type varType;
	//array of dimensions
	int pdimids[3]; //NC_MAX_DIMS]; 3D file only being read here; expected to get error message otherwise
	//dimension names 
	char pdim_Names[80];
	size_t start[3], count[3];
	//Open the netcdf file.  
	if ((retncval = nc_open_par(FILE_NAME, NC_NOWRITE, &ncid)))         // NC_MPIIO, inpComm, inpInfo,
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);
	// Get the varids of the coordinate variables 
	if ((retncval = nc_inq_varid(ncid, tcor_NAME, &ptid)))
		ERR(retncval);
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
		ERR(retncval);	
	for (int j = 0; j < bSize; j++)
	{ 
		//check dimension info and set start and count arrays; 
		for (int i = 0; i < 3; i++)
		{
			if ((retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes)))
				ERR(retncval);
			if (strcmp(pdim_Names, tcor_NAME) == 0)
			{
				start[i] = 0;
				nrecords = pdim_sizes;
				count[i] = nrecords;
				//allocate memory block for the 'data' varialbe 		
				pvar_in[j] = new float[nrecords];
				//pvarin_temp = new float[nrecords];
				//tcorvar = new float[nrecords];
			}
			else if (strcmp(pdim_Names, ycor_NAME) == 0)
			{
				start[i] = YindArr[j];
				count[i] = 1;
			}
			else
			{
				start[i] = XindArr[j];      // xdim;
				count[i] = 1;				
			}
		}		
		//read var data
		if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[j][0]))
			ERR(retncval);
		/* Read the coordinate (dimensions) variable data. */
	}
	/*if (retncval = nc_get_var_float(ncid, ptid, &tcorvar[0]))
		ERR(retncval);*/
	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval);
	return 0;
}

//function to read single column/rod along time dimension from 3D netcdf file, for given y , x coordinates
int readNC_TS(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME,const char* xcor_NAME,
	float* &pvar_in, float* &tcorvar, int ydim, int xdim, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo) //, int xstride) //double* tcorvar, double*** pvar_in)
{
	//6.24.14
	int xdimindx = 1;
	float* pvarin_temp = NULL;
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, ptid = 0; // pxid = 0, pyid = 0, ndims = 0; 
	//dimensions lengths
	size_t pdim_sizes; //[3]; //Ntdim = 0, Nxdim = 0, Nydim = 0;
	//variable data type
	nc_type varType;
	//array of dimensions
	int pdimids[3]; //NC_MAX_DIMS]; 3D file only being read here; expected to get error message otherwise
	//dimension names 
	char pdim_Names[80];
	size_t start[3], count[3];
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))		
		ERR(retncval); 
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid))) 
		ERR(retncval); 	
	// Get the varids of the coordinate variables 
	if ((retncval = nc_inq_varid(ncid, tcor_NAME, &ptid))) 
		ERR(retncval);
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid,pvarid,NULL,&varType,NULL,pdimids,NULL))) 
		ERR(retncval); 
	 //check dimension info and set start and count arrays; 
	for (int i =0; i < 3; i++)	
	{
		if ((retncval = nc_inq_dim(ncid,pdimids[i], pdim_Names, &pdim_sizes))) 
			ERR(retncval); 
	    if (strcmp(pdim_Names, tcor_NAME) == 0 )
		{
			start[i] = 0;		    
			nrecords = pdim_sizes;
			count[i] = nrecords;
			//allocate memory block for the 'data' varialbe 
			pvar_in = new float[nrecords];                    //6.23.14
			pvarin_temp = new float[nrecords];
			tcorvar = new float[nrecords];
		}
		else if (strcmp(pdim_Names, ycor_NAME) == 0 )
		{
			start[i] = ydim;
			count[i] = 1;
		}
		else 
		{
			start[i] = xdim;
			count[i] = 1;     // xstride;                   //6.23.14
			xdimindx= i;
		}
	}	
	//read var data	
	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0]))
		ERR(retncval);	
		/* Read the coordinate (dimensions) variable data. */ 
	if (retncval = nc_get_var_float(ncid, ptid, &tcorvar[0]))
		ERR(retncval);
	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval); 

	return 0;
}

int read2DNC(const char* FILE_NAME, const char* VAR_NAME, float** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 
	//dimensions lengths
	size_t Nxdim = 0, Nydim = 0;
	//variable data type
	nc_type varType;
	//array of dimensions
	int pdimids[2]; //NC_MAX_DIMS]; 2D file only being read here; expected to get error message otherwise
	//size_t start[3], count[3];
	//Open the file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval); 
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid))) 
		ERR(retncval); 
	// Get the varids of the coordinate variables 	
	/*if ((retncval = nc_inq_varid(ncid, ycor_NAME, &pyid))) 
		ERR(retncval); 
	if ((retncval = nc_inq_varid(ncid, xcor_NAME, &pxid))) 
		ERR(retncval); 	*/
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid,pvarid,NULL,&varType,NULL,pdimids,NULL))) 
		ERR(retncval); 
	//check dimension sizes
	if ((retncval = nc_inq_dim(ncid,pdimids[0],NULL,&Nydim))) 
		ERR(retncval); 
	if ((retncval = nc_inq_dim(ncid,pdimids[1],NULL,&Nxdim))) 
		ERR(retncval); 	
	//create the output arrays
	pvar_in = new float*[Nydim];
	for(size_t j=0; j< Nydim; j++)			
		pvar_in[j] = new float [Nxdim];
	float* ta_inp = new float[Nxdim];	  //8.14.13 from ...[Nydim]
	size_t start[2], count[2];
	start[0] = 0;   
	start[1] = 0; 
	count[0] = 1; 
	count[1] = Nxdim;	
	for (int nj = 0; nj < Nydim; nj++)  
	{
		start[0] = nj;
		if ((retncval = nc_get_vara_float(ncid, pvarid,start,count, &ta_inp[0])))      
			ERR(retncval); 
		for(int ix=0; ix< Nxdim; ix++)			
				pvar_in[nj][ix] = ta_inp[ix];
		//cout<<"step no %d\n",nj);
	}	/* next record */ 
	/*ycorvar = new float[Nydim];
	xcorvar = new float[Nxdim];
	// Read the coordinate (dimensions) variable data.  
	if ((retncval = nc_get_var_float(ncid, pyid, &ycorvar[0])))
		ERR(retncval);
	if ((retncval = nc_get_var_float(ncid, pxid, &xcorvar[0]))) 
		ERR(retncval); */ 
	 //read variable (input data)	
	/*if ((retncval = nc_get_var_float(ncid, pvarid, &pvar_in[0][0])))      
			ERR(retncval); 	*/
	//close netcdf file			
	if ((retncval = nc_close(ncid))) 
		ERR(retncval); 

	return 0;
}
//the following is to read watershed file
int readwsncFile(const char* FILE_NAME, const char* VAR_NAME, const char* ycor_NAME,  
	const char* xcor_NAME, float* &ycorvar, float* &xcorvar, int** &pvar_in, int &ydim, int &xdim, int &fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo)       //, std::set<int> zValues,  float * z_ycor, float *z_xcor)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0, pxid = 0, pyid = 0, ndims = 0; 
	//dimensions lengths
	size_t Nxdim = 0, Nydim = 0;
	int fillSet = 0;
	//variable data type
	nc_type varType, missingType;
	//array of dimensions
	int pdimids[2]; //NC_MAX_DIMS]; 2D file only being read here; expected to get error message otherwise
	//size_t start[3], count[3];
	//Open the file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval); 
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid))) 
		ERR(retncval); 
	// Get the varids of the coordinate variables 	
	if ((retncval = nc_inq_varid(ncid, ycor_NAME, &pyid))) 
		ERR(retncval); 
	if ((retncval = nc_inq_varid(ncid, xcor_NAME, &pxid))) 
		ERR(retncval); 	
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid,pvarid,NULL,&varType,NULL,pdimids,NULL))) 
		ERR(retncval); 	  
	//Arcgis uses "missing_vaue TBC # 12.18.14 
	//CF Convension use _FillValue
	if ((retncval = nc_get_att(ncid,pvarid,"_FillValue",&fillVal)))    
		ERR(retncval); 		
	//cout<<"_FillValue: "<<fillVal<<endl;
	/*int iMiss;
	 if ((retncval = nc_inq_var_fill(ncid,pvarid, &fillSet,&iMiss)))     
		 	ERR(retncval); 		
	cout<<" Fill value: "<<iMiss<<endl;
	cout<<" Fill set? "<<fillSet<<endl;*/
	
	//check dimension sizes
	if ((retncval = nc_inq_dim(ncid,pdimids[0],NULL,&Nydim))) 
		ERR(retncval); 
	if ((retncval = nc_inq_dim(ncid,pdimids[1],NULL,&Nxdim))) 
		ERR(retncval); 
	//get dim values
	ydim = (int)Nydim;
	xdim = (int)Nxdim;
	
	//create the output arrays
	pvar_in = new int*[Nydim];
	for(size_t j=0; j< Nydim; j++)			
		pvar_in[j] = new int [Nxdim];
	ycorvar = new float[Nydim];
	xcorvar = new float[Nxdim];
	//
	// Read the coordinate (dimensions) variable data.  
	if ((retncval = nc_get_var_float(ncid, pyid, &ycorvar[0])))
		ERR(retncval);
	if ((retncval = nc_get_var_float(ncid, pxid, &xcorvar[0]))) 
		ERR(retncval);  
	 //read variable (input data)	
	int* ta_inp = new int[Nxdim];	 //8913 frm size of nydim
	size_t start[2], count[2];
	start[0] = 0;   
	start[1] = 0; 
	count[0] = 1; 
	count[1] = Nxdim;
	
	for (int nj = 0; nj < Nydim; nj++)  
	{
		start[0] = nj;
		if ((retncval = nc_get_vara_int(ncid, pvarid,start,count, &ta_inp[0])))      
			ERR(retncval); 
		for(int ix=0; ix< Nxdim; ix++)			
				pvar_in[nj][ix] = ta_inp[ix];
		//cout<<"step no %d\n",nj);
	}	/* next record */ 	
    //cout<<endl;
	/*if ((retncval = nc_get_vara_int(ncid, pvarid,start,count, &pvar_in[0][0])))      
			ERR(retncval); 	*/
	//close netcdf file	
	if ((retncval = nc_close(ncid))) 
		ERR(retncval); 

	return 0;
}

int Write3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME, 
	const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_out, MPI::Intracomm inpComm, MPI::Info inpInfo)
{   
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0 , pxvarid = 0, pyvarid = 0, ptvarid = 0, pxid = 0, pyid = 0, ptid = 0, ndims = 0; 
	//variable data type
	//nc_type varType;
	//array of dimensions
	int pdimids[3]; //NC_MAX_DIMS]; 3D file only being read here; expected to get error messege otherwise
	size_t start[3], count[3];
	//dimensions lengths
	size_t Ntdim = 0, Nxdim = 0, Nydim = 0;
	//computing the coordinate arrays sizes
	//#_5 tbc 5.19.13
	Ntdim = sizeof(tcorvar)/sizeof(tcorvar[0]);
	Nydim = sizeof(ycorvar)/sizeof(ycorvar[0]);
	Nxdim = sizeof(xcorvar)/sizeof(xcorvar[0]); 
	if ((retncval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))  
		ERR(retncval); 
    // Define the dimensions. The time dim is assumed the record dimension and has unlimited length
	if ((retncval = nc_def_dim(ncid, tcor_NAME, NC_UNLIMITED, &ptid))) 
		ERR(retncval); 
	if ((retncval = nc_def_dim(ncid, ycor_NAME, Nydim, &pyid)))  
		ERR(retncval);   
	if ((retncval = nc_def_dim(ncid, xcor_NAME, Nxdim, &pxid)))   
		ERR(retncval);  
	// Define the coordinate variables 
	//the variables here are given same name as the corresponding coordinates
	//#_6 tbc 5.19.13
	if ((retncval = nc_def_var(ncid, tcor_NAME, NC_DOUBLE, 1, &ptid,   
		&ptvarid)))     
		ERR(retncval);  
	if ((retncval = nc_def_var(ncid, ycor_NAME, NC_DOUBLE, 1, &pyid,   
		&pyvarid)))     
		ERR(retncval);  
	if ((retncval = nc_def_var(ncid, xcor_NAME, NC_DOUBLE, 1, &pxid,   
		&pxvarid)))     
		ERR(retncval);
	/* The dimids array is used to pass the dimids of the dimensions of  
	the netCDF variables.  In C, the unlimited dimension must come first on the list of dimids.
	The number of dimensions in this function is 3*/ 
	pdimids[0] = ptid; 
	pdimids[1] = pyid; 
	pdimids[2] = pxid; 
	// Define the output variable  
	if ((retncval = nc_def_var(ncid, VAR_NAME, NC_DOUBLE, 3, pdimids, &pvarid)))   
		ERR(retncval);  
	//##_7 MtoAdd (5.19.13)
	/*------need to write more code lines (possibly call from outside)
	to asign attributes to each coordinate and output variables
	*//*	 
	if ((retncval = nc_put_att_text(ncid, ptvarid, UNITS,   
		strlen(ptvar_UNITS), ptvar_UNITS))) 
		ERR(retncval);  
     */ 
	//End define mode. 
	if ((retncval = nc_enddef(ncid)))   
		ERR(retncval); 	
	//set start and count arrays; time, if exists, is expected to to come first
	//#_8 tbc 5.19.13
	start[0] = 0;   
	start[1] = 0;  
	start[2] = 0;  
	count[0] = 1;  
	count[1] = Nydim; 
	count[2] = Nxdim;
	//temporary 2D array to help in data write 	//#_9 tbc 5.19.13 should write 3d @ once 
	float ** pvar_outp = new float*[Nydim];
	for(size_t j=0; j< Nydim; j++)			
		pvar_outp[j] = new float [Nydim];	
	// Write the coordinate variable data. 
	if ((retncval = nc_put_var_float(ncid, ptvarid, &tcorvar[0]))) 
		ERR(retncval); 
	if ((retncval = nc_put_var_float(ncid, pyvarid, &ycorvar[0]))) 
		ERR(retncval); 
	if ((retncval = nc_put_var_float(ncid, pxvarid, &xcorvar[0]))) 
		ERR(retncval); 
	//write output variable data
	for (size_t kt = 0; kt < Ntdim; kt++)  
	{
		for(size_t i=0; i< Nydim; i++)
			for(size_t j=0;j<Nxdim;j++)
				pvar_outp[i][j] = pvar_out[kt][i][j];
		start[0] = kt;
		if ((retncval = nc_put_vara_float(ncid, pvarid,start,count, &pvar_outp[0][0])))      
			ERR(retncval); 		
		//cout<<"step no %d\n",k);
	}	/* next record */ 
	//free temporary matrix pvar_outp
	for(size_t i=0; i< Nydim; i++)
		delete []pvar_outp[i];
	delete []pvar_outp;
	//close netcdf file			
	if ((retncval = nc_close(ncid))) 
		ERR(retncval); 
	cout<<"SUCCESS writing output file: "<<FILE_NAME<<endl;  
	return 0;	
}
