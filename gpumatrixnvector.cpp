//#include "nctest.h"
#include "gpuuebpgdecls.h"
///******Warning the following code seems to have memo leak------ 5.7.15 
//create 3D array and allocate contiguous memory block this enbales a full block read of netcdf
__host__ __device__ float*** create3DArrayblock_Contiguous(int nt, int nr, int nc)      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
{
	float*** myMatrix = new float**[nt];
	float* ptrMemory = new float[nt*nr*nc];
	for (int t = 0; t<nt; t++)
	{
		//this looks suspicious 
		//#??_TBC 6.18.13
		myMatrix[t] = new float*[nr];
		for (int r = 0; r<nr; r++)
		{
			myMatrix[t][r] = ptrMemory; //new float*[nr];  
			ptrMemory += nc;
		}
	}
	return myMatrix;

}//float*** create3DArrayblock

//delets a 3D array (frees memory) allocated contiguously
__host__ __device__ void delete3DArrayblock_Contiguous(float*** myMatrix) // int nr, int nc)// input: 3D array
{
	/*for (int t = 0; t< nt; t++)
	{
	/*for (int r=0; r<nr; r++)
	delete [] myMatrix[t][r];

	}*/
	delete[] myMatrix[0];
	delete[] myMatrix;
	return;

}//double** CreateMatrix
//******Warning the following code seems to have memo leak------ 5.7.15 
//create 2D array and allocate contiguous memory block this enbales a full block read of netcdf
__host__ __device__ float** create2DArray_Contiguous(int nr, int nc)      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
{
	float** myMatrix = new float*[nr];
	float* ptrMemory = new float[nr*nc];
	for (int r = 0; r<nr; r++)
	{
		myMatrix[r] = ptrMemory; //new float*[nr];  
		ptrMemory += nc;
	}
	return myMatrix;
}//float*** create3DArrayblock

//delets a 2D array (frees memory allocated) contiguously
__host__ __device__ void delete2DArray_Contiguous(float** myMatrix) // int nr, int nc)// input: 2D array
{
	delete[] myMatrix[0];
	delete[] myMatrix;

	return;
}//double** 

//Creates a matrix. The inputs are matrix dimensions. It allocates a memory block of size nrows*ncols* (size of float)
//and returns an array of pointers to the allocated memory block
__host__ __device__ float*** Create3DArray(int nt, int nr, int nc)       //inputs: no. of rows and no. of colos (dimensions) of matrix
{
	float*** myMatrix = new float**[nt];
	for (int i = 0; i<nt; i++)
	{
		myMatrix[i] = new float*[nr];
		for (int j = 0; j< nr; j++)
			myMatrix[i][j] = new float[nc];
	}
	return myMatrix;

}//float** CreateMatrix

//The following program deletes a matrix passed to it (it frees up the memory block allocated to the matrix)
__host__ __device__ void Delete3DArray(float ***A, int nt, int nr, int nc)            //input: A matrix
{
	for (int i = 0; i< nt; i++)
	{
		for (int j = 0; j< nr; j++)
			delete[] A[i][j];
		delete[] A[i];
	}
	delete[] A;

	return;
}//void DeleteMatrix