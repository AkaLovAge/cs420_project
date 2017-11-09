/*
* a parallel huge matrix txt writer.
* M is the row number 
* N is the column number 
*/
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
int main(int argc, char* argv[])
{
	int length;
	int error;
	int m=1000; // row number
	int n=1000; // column number
	int total_size = m*n*6 + m;
	
	int size,rank;
	MPI_Status status;
	MPI_File fh;
	MPI_Offset filesize;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int line_size = n*6+1;
	int fsize = total_size+size;
	int end,start,row;

	if (rank == size -1)
	{
		row = m / size + m % size;
		start = fsize - (row*n*6 + row +1);
		end = total_size;
	}else{
		row = m/size;
		start = rank*(row*n*6 + row +1);
		end = start + row*n*6 + row +1;
	}

	char * buffer = malloc(sizeof(char)*(end-start));
	int i,j;
	int index = 0;
	printf("rank:%d",rank);
	printf("generating buffer in size: %d\n",(end-start));
	for (i=0;i<row;i++)
	{
		for(j=0;j<n-1;j++)
		{
			char s[7];
			double kk = (float)rand()/((float)RAND_MAX/10);
			sprintf(s , "%.3lf,",kk);
			buffer[index] = s[0];
			buffer[index+1] = s[1];
			buffer[index+2] = s[2];
			buffer[index+3] = s[3];
			buffer[index+4] = s[4];
			buffer[index+5] = s[5];
			index+=6;
		}
	//	printf("insert last element in row\n");
		char ss[8];
		double kk = (float)rand()/((float)RAND_MAX/10);
		sprintf(ss,"%.3lf\n",kk);
		buffer[index] = ss[0];
		buffer[index+1] = ss[1];
		buffer[index+2] = ss[2];
		buffer[index+3] = ss[3];
		buffer[index+4] = ss[4];
		buffer[index+5] = ss[5];
		buffer[index+6] = ss[6];
		index += 7;
		
		//printf("%s",buffer);
		//printf("%d,%d,%d\n",rank,i,index);

	}
	printf("%d\n%s\n",rank,buffer);
	error = MPI_File_open(MPI_COMM_WORLD,"test.txt", MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

	if (error != MPI_SUCCESS) {MPI_Finalize();exit(-1);}

	error = MPI_File_write_at(fh, start, buffer, end-start, MPI_DOUBLE, &status);
	if(error != MPI_SUCCESS) {MPI_Finalize();exit(-1);}
	
	MPI_File_close(&fh);

	free(buffer);
	return 0;
}
