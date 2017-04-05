#include <unistd.h>
//#include <mraa/aio.h>
#include <stdio.h>
//#include "floatfann.h"
#include <string.h>
#include <stdlib.h>
#define BUFF_SIZE 1024

const char * walk_modes[] = {
        "walk_speed_1_50sec_32m",
        "walk_speed_2_35sec_32m",
        "walk_speed_3_25sec_32m",
        "walk_speed_4_15sec_32m",
                "failure"
};


const unsigned int num_input = 19;
const unsigned int num_output = 4;

int main()
{
//    uint16_t value0, value1, value2, value3;


    int i, j, k, rv;
    int walk_speed;
    int *output;
    int N_SAMPLES;
    float **all_data;
    char *ifile_name;
    FILE *fp;
    char *line = NULL;
    char buffer[1024];
    char *record, *csv_line;
    size_t len = 0;
    ssize_t read;

//    float max;
//    fann_type *calc_out;
//    fann_type input[num_input];
//    struct fann *ann;
//    ann = fann_create_from_file("Walk_Test.net");


    ifile_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
    memset(ifile_name, 0, BUFF_SIZE);
    snprintf(ifile_name,
             BUFF_SIZE,
             "/Users/amaaelantonini/Documents/Winter 2017/EE 180DA/File_Example/test_file.csv"
    );
    fp = fopen(ifile_name, "r");
    if (fp == NULL) {
        fprintf(stderr,
                "Failed to read from file \'%s\'.\n",
                ifile_name
        );
        exit(EXIT_FAILURE);
    }

    /* count the number of lines in the file */
    read = getline(&line, &len, fp); //discard header of file
    N_SAMPLES = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        N_SAMPLES++;
    }

    /* go back to the start of the file so that the data can be read */
    rewind(fp);
    output = (int *) malloc(sizeof(int*) * num_output);
    all_data = (float **) malloc(sizeof(float **) * N_SAMPLES);
    for(i = 0; i < N_SAMPLES; i++)
    {
        all_data[i] = (float *) malloc(sizeof(float) * num_input);
        for(j = 0; j < num_input; j++)
            printf("%f \n", all_data[i][j]);
        printf("next: %d\n", i);
    }
    read = getline(&line, &len, fp); //discard header of file
    j = 0;
    i = 0;
    while ((read = getline(&line, &len, fp)) != -1)
    {
        j = 0;
//        printf("%s\n", line);
        record = strtok(line, ",");
//        printf("%s, %d\n", record, i);
        while (record != NULL) {
//            printf("record : %s\n", record);
//            printf("converted : %f\n", atof(record));
            all_data[i][j++] = (float) atof(record);
            record = strtok(NULL, ",");
        }
        ++i;
    }
    fclose(fp);

    ann = fann_create_from_file("Walk_Test.net");
    int reader;
    walk_speed = 0;
    k = 0;
    for(reader = 0; reader < N_SAMPLES; reader ++) {
        max = -100;
        for(j = 0; j < num_input; j++)
            input[j] = all_data[reader][j];
        calc_out = fann_run(ann, input);
        for (i = 0; i < 5; i++) {
            output[i] = 0;
            if (calc_out[i] > max) {
                max = calc_out[i];
                k = i;
            }
        }
        if(max < 0)
            walk_speed = 4;
        else
            output[k] = 1;

        // print 0 for all, but 1 for the max
        printf("%4.5f,%4.5f,%4.5f,%4.5f, %s"
                , all_data[reader][0], all_data[reader][6], all_data[reader][7], all_data[reader][13]);
        printf("%d\t%d\t%d\t%d\n", output[0], output[1], output[2], output[3]);
    }
    fann_destroy(ann);
    return 0;
}
