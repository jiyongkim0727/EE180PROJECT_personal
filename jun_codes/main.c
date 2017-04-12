#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fann.h"

//[GLOBAL VARIABLES]
//format: <input_csv> <output_csv>
int i, idx;
int rv;

//file reader
char *ifile_name, *ofile_name;
FILE *fp;
char *line = NULL;
size_t len = 0;
ssize_t read;
int N_SAMPLES;
	
//data variables
double *t, *t1,*t2; 
float *x, *y, *z,*gx,*gy,*gz;
float *S_i;
int n_S;



//[FUNCTIONS]


void get_input_file()
{
	printf("Attempting to read from file %s\n", ifile_name);
	
	//opening input file
	fp = fopen(ifile_name, "r");
	if(fp==NULL)
	{
		fprintf(stderr,"Error: cannot open the input file.\n");
		exit(1);
	}

	//count the # of lines
	read = getline(&line,&len,fp);
	N_SAMPLES=0;
	while((read=getline(&line,&len,fp)) != -1)
		N_SAMPLES++;

	//[READING the DATA]
	//initiating
	rewind(fp);
	read = getline(&line, &len, fp);
	i = 0; 
	t = (double *) malloc(sizeof(double) * N_SAMPLES);
	t1 = (double *) malloc(sizeof(double) * N_SAMPLES);
	t2 = (double *) malloc(sizeof(double) * N_SAMPLES);
	x = (float *) malloc(sizeof(float) * N_SAMPLES);
	y = (float *) malloc(sizeof(float) * N_SAMPLES);
	z = (float *) malloc(sizeof(float) * N_SAMPLES);
	gx = (float *) malloc(sizeof(float) * N_SAMPLES);
	gy = (float *) malloc(sizeof(float) * N_SAMPLES);
	gz = (float *) malloc(sizeof(float) * N_SAMPLES);
	
	//format <time_before> <time_after> <x> <y> <z> <gx> <gy> <gz>
	while((read=getline(&line,&len,fp)) != -1)
	{
		rv = sscanf(line, "%lf,%lf,%f,%f,%f,%f,%f,%f\n", &t1[i], &t2[i],
				&x[i],&y[i], &z[i], &gx[i], &gy[i], &gz[i]);
		if(rv != 8)
		{
			fprintf(stderr, "Error: incorrect csv file format %d\n"
					,rv);
			exit(1);
		}

		i++;
	}
	//make into a sample time
	for(i=0; i<N_SAMPLES; i++)
	{
		if(i==0)
			t[0] = t2[0] - t1[0];
		else
		{
			t[i] = t[i-1] + t2[i] - t1[i];
		}
	}
	fclose(fp);

	return;
}
int detect_stair_descend(float* gz_arr)
{
}

void write_to_file()
{
	printf("Attempting to write to file %s\n", ofile_name);
	fp = fopen(ofile_name, "w");
	if(fp==NULL)
	{
		fprintf(stderr, "Error: cannot open the file.\n");
		exit(1);
	}

	//ouput format <stride_index> <stride_time>
	//right now, just outputing gz
	//
	n_S=0;
	fprintf(fp, "S_i,S_t,S_x\n");
	for(i=0; i<N_SAMPLES; i++)
	{
		fprintf(fp, "%d,%20.10lf,%lf\n",i,t[i],gz[i]); 
	}
	fclose(fp);

	return;
}

void train_data(char* input_file, const unsigned int n_input, 
		const unsigned int n_output)
{
	//THESE VARIABLES ASSUME UNCHANGED
	const unsigned int num_layers = 3;
	const unsigned int num_neurons_hidden = 9;
	const float desired_error = (const float) 0.001;
	const unsigned int max_epochs = 10000;
	const unsigned int epochs_between_reports = 100;

	struct fann *ann = fann_create_standard(num_layers, n_input,
		       num_neurons_hidden, n_output);

	fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
	fann_set_activation_function_output(ann, FANN_SIGMOID_SYMMETRIC);
	
	/*fann_train_on_file(ann, input_file, max_epochs, 
			epochs_between_reports, desired_error);*/ // <- to run
	//fann_save(ann, "TEST.net"); <- do we want this?
	
	fann_destroy(ann);

	return;
}

int main(int argc, char **argv)
{
	/*[ARGUMENT CHECKS]*/
	if(argc!=3)
	{
		fprintf(stderr,"Error: incorrect number of arguments.\n");
		exit(1);
	}
	ifile_name = argv[1];
	ofile_name = argv[2];


	/*[READING]*/
	get_input_file();
	

	/*[DETECTION ALGORITHM]*/
	detect_stair_descend(gz);


	/*[WRITING]*/
	write_to_file();


	/*[TRAIN DATA]*/
	train_data("test_data.txt",3,4);


	exit(0);
}
