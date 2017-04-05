#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Functions.h"
#include <sys/time.h>
#include "classifiers.h"

#define BUFF_SIZE 1024
#define STRIDES 200
#define FILES 12
#define DATA 7
#define NULL 0
#define SLOTS 4
const char * names[] = {"AA2", "DA2"};
const char * file_names[] = {
        "walk_speed_1_50sec_32m",
        "walk_speed_2_35sec_32m",
        "walk_speed_3_25sec_32m",
        "walk_speed_4_15sec_32m",
        "slow_run",
        "medium_run",
        "fast_run",
        "stairs_up",
        "stairs_down",
        "low_jump",
        "medium_jump",
        "high_jump"};


const float times_AA[] = { 0, 30, 0, 30, 0, 30, 2, 22, 0, 15, 0, 12, 2, 12,0, 22, 16, 30,18, 30,0, 20, 0, 15};
const float times_DA[] = {0, 30, 0, 30, 2, 20, 0, 16, 2, 18, 2, 15, 2, 12, 0, 12, 7, 28, 0, 12, 0, 10, 0, 12};
const char * main_path = "/Users/amaaelantonini/Documents/Winter 2017/EE 180DA/File_Example/";
const char * classifiers = "Classifiers/";
const char * class_files[] = {"Classifiers/AA_file/", "Classifiers/DA_file/"};

const char * paths[] = {"Data_sets_3/AA2_data/", "Data_sets_3/DA2_data/"};
const char * paths_out[] = {"Data_sets_out/AA2_data/", "Data_sets_out/DA2_data/"};
const char * stride_path = "Data_strides_out/";

void process() {
    int rv;
    int i, j, k;
    float t, t2, ax, ay, az, gx, gy, gz;
    char **input_files;
    char **output_files;
    char **strides_out;
    char **parameters;
    float ***data;
    float ***sigma;
    float ***mean;
    int ** strides_t1, *count_t1;
    int ** strides_t2, *count_t2;
    float *data2[12][7];
//    int *output;
    int N_SAMPLES;
    float **all_data;
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char buffer[1024];
    char *record, *csv_line;
    int *data_lengths;
    double time;
    char *store_last;


    count_t1 = (int*) malloc(sizeof(int) * FILES);
    count_t2 = (int*) malloc(sizeof(int) * FILES);
    strides_t1 = (int**) malloc(sizeof(int*) * FILES);
    strides_t2 = (int**) malloc(sizeof(int*) * FILES);
    for(i = 0; i < FILES; i++)
    {
        strides_t1[i] = (int*) malloc(sizeof(int) * STRIDES);
        strides_t2[i] = (int*) malloc(sizeof(int) * STRIDES);
    }
    data_lengths = (int*) malloc(sizeof(int) * FILES);
    data = (float ***) malloc(FILES * sizeof(float **));
    mean = (float ***) malloc(FILES * sizeof(float **));
    sigma = (float ***) malloc(FILES * sizeof(float **));
    input_files = (char **) malloc(FILES * sizeof(char *));
    for (i = 0; i < FILES; i++) {
        input_files[i] = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(input_files[i], 0, BUFF_SIZE);
        snprintf(input_files[i],
                 BUFF_SIZE,
                 "%s%s%s_%s.csv",main_path, paths[1], file_names[i], names[1]);
        data[i] = (float **) malloc(sizeof(float*) * DATA);
        mean[i] = (float **) malloc(sizeof(float*) * DATA);
        sigma[i] = (float **) malloc(sizeof(float*) * DATA);
    }


    for (i = 0; i < FILES; i++)
    {
//        printf("Starting\n");
        fp = fopen(input_files[i], "r");
        if (fp == NULL) {
            fprintf(stderr,
                    "Failed to read from file \'%s\'.\n",
                    input_files[i]
            );
            exit(EXIT_FAILURE);
        }
//        printf("File opened\n");
        /* count the number of lines in the file */
        read = getline(&line, &len, fp); //discard header of file
        N_SAMPLES = 0;
        while ((read = getline(&line, &len, fp)) != -1) {
//            if(N_SAMPLES<100)
//                printf("%s\n", line);
            N_SAMPLES++;
        }
//        *lengths[i] = N_SAMPLES;
//        printf("samples counted: %d\n", N_SAMPLES);
        data_lengths[i] = N_SAMPLES;
        for(j = 0; j < DATA; j++)
        {
            data[i][j] = (float *) malloc(sizeof(float) * N_SAMPLES);
            mean[i][j] = (float *) malloc(sizeof(float) * N_SAMPLES);
            sigma[i][j] = (float *) malloc(sizeof(float) * N_SAMPLES);
        }
//        printf("data initialized\n");
        rewind(fp);
        read = getline(&line, &len, fp); //discard header of file
//        printf("file rewined\n");
        j = 0;
        N_SAMPLES = 0;
        k = 0;
        while ((read = getline(&line, &len, fp)) != -1) {
            k = 0;
            N_SAMPLES++;
//            printf("%s\n", line);
            record = strtok(line, ",");
            while (record != NULL) {
//                if(k == 1)
//                    record = strtok(NULL, ",");
                data[i][k++][j] = (float) atof(record);
                record = strtok(NULL, ",");

            }
            j++;
        }
        fclose(fp);
//        *lengths[i] = N_SAMPLES;
//        printf("samples counted: %d\n", N_SAMPLES);
        rewind(fp);
//        printf("%d\n", i);

    }
    for(i = 0; i < 0*FILES; i++) {
        for (j = 0; j < data_lengths[i]; j+=3000) {
            printf("%20.10lf,%20.10lf,%10.10lf,%10.10lf,%10.10lf,%10.10lf,%10.10lf\n",
                    data[i][0][j],
                    data[i][1][j],
                    data[i][2][j],
                    data[i][3][j],
                    data[i][4][j],
                    data[i][5][j],
                    data[i][6][j]

            );
        }
    }
//    printf("offset\n");
    float off_set;
    for(i = 0; i < 0; i++)
    {
        off_set = data[i][0][0];
//        printf("offset: %f\n", off_set);
        for(j = 0; j < data_lengths[i]; j++)
        {
//            printf("time: %f\n", data[i][0][j]);
            data[i][0][j] = data[i][0][j]- off_set;
        }
    }

    printf("mean and sigma\n");
    for(i = 0; i < FILES; i++)
    {
//        printf("i: %d\n", i);
        for(j = 0; j < DATA-1; j++)
        {

            sample_mean(data[i][j+1], mean[i][j], 3000, data_lengths[i]);
//            printf("mean done: %d\n", j);
            sample_sigma(data[i][j+1], mean[i][j], sigma[i][j], 3000, data_lengths[i]);
//            printf("sigma done: %d\n", j);
            center_data(data[i][j+1], mean[i][j], sigma[i][j], 3000, data_lengths[i]);
//            printf("center done: %d\n", j);
        }
//        printf("next_file\n");
        for (j = 100000; j < data_lengths[i]; j+=2900)
        {
//            printf("mean:   %15.10lf,%15.10lf,%10.10lf,%10.10lf,%10.10lf,%10.10lf\n",
//                   mean[i][0][j],
//                   mean[i][1][j],
//                   mean[i][2][j],
//                   mean[i][3][j],
//                   mean[i][4][j],
//                   mean[i][5][j]
//            );
            printf("sigma:   %15.10lf,%15.10lf,%10.10lf,%10.10lf,%10.10lf,%10.10lf\n",
                   sigma[i][0][j],
                   sigma[i][1][j],
                   sigma[i][2][j],
                   sigma[i][3][j],
                   sigma[i][4][j],
                   sigma[i][5][j]
            );
        }
    }


    for (i = 0; i < FILES; i++) {
        count_t1[i] = temp_strides(data[i], strides_t1[i], sigma[i][1], sigma[i][5], 4, 500, 3, 0.2, data_lengths[i]);
//        printf("%d\n", count_t1[i]);
//        printf("[\n");
//        for(j = 0; j < count_t1[i]; j++)
//        {
//            printf("%d, ", strides_t1[i][j]);
//
//        }
//        printf("],\n");
        printf("\n");
        count_t2[i] = select_strides(data[i][0], strides_t1[i], strides_t2[i], 0.4, count_t1[i]);
//        printf("%d\n", count_t2[i]);
        printf("[\n");
        for(j = 0; j < STRIDES; j++)
            strides_t1[i][j] = 0;
        for(j = 0; j < count_t2[i]; j++)
        {
            printf("%d, ", strides_t2[i][j]);
        }
        printf("],\n");
        printf("\n");
        count_t1[i] = peak_strides(data[i][0], data[i][6], strides_t2[i], strides_t1[i], sigma[i][5], 1.4, count_t2[i], data_lengths[i]);
//        printf("strides temp: %d, selected: %d\n", count_t1[i], count_t2[i]);
//        count_t2[i] = select_strides(data[i][0], strides_t2[i], strides_t1[i], 0.4, count_t1[i]);
//        count_t1[i] = count_t2[i];
        count_t2[i] = select_strides(data[i][0], strides_t1[i], strides_t2[i], 0.01, count_t1[i]);
//        printf("before outlier\n");
        count_t1[i] = outlier(data[i][0], strides_t2[i], strides_t1[i], 0.65, count_t2[i]);
//        printf("second\n");
        count_t2[i] = outlier(data[i][0], strides_t1[i], strides_t2[i], 0.65, count_t1[i]);
//        for(j = 0; j < count_t1[i]; j++)
//        {
//            strides_t2[i][j] = strides_t1[i][j];
//        }
//        count_t2[i] = count_t1[i];
    }
    printf("last attempt\n");
    for (i = 0; i < FILES; i++)
    {
//        printf("%d\n", count_t1[i]);
        printf("[\n");
        for(j = 0; j < count_t2[i]; j++)
        {
            printf("%d, ", strides_t2[i][j]);
        }
        printf("\n");
        printf("],\n");
    }
    for(i = 0; i < FILES; i++)
        move_back(strides_t2[i], 100, count_t2[i]);
//        move_back(strides_t1[i], 100, count_t1[i]);


    float ** sigma_high_xa, ** sigma_low_xa;
    float ** sigma_high_ya, ** sigma_low_ya;
    float ** sigma_high_zg, ** sigma_low_zg;

    float *** troughs_xa, *** peaks_xa;
    float *** troughs_ya, *** peaks_ya;
    float *** troughs_zg, *** peaks_zg;


    float *** mean_low_xa, *** mean_high_xa;
    float *** mean_low_ya, *** mean_high_ya;
    float *** mean_low_zg, *** mean_high_zg;

    int *** stride_all;
    float *** stride_float;
    int *num_strides;

    sigma_high_xa = (float **) malloc(FILES * sizeof(float *));
    sigma_low_xa = (float **) malloc(FILES * sizeof(float *));
    sigma_high_ya = (float **) malloc(FILES * sizeof(float *));
    sigma_low_ya = (float **) malloc(FILES * sizeof(float *));
    sigma_high_zg = (float **) malloc(FILES * sizeof(float *));
    sigma_low_zg = (float **) malloc(FILES * sizeof(float *));

    mean_high_xa = (float ***) malloc(FILES * sizeof(float **));
    mean_low_xa = (float ***) malloc(FILES * sizeof(float **));
    mean_high_ya = (float ***) malloc(FILES * sizeof(float **));
    mean_low_ya = (float ***) malloc(FILES * sizeof(float **));
    mean_high_zg = (float ***) malloc(FILES * sizeof(float **));
    mean_low_zg = (float ***) malloc(FILES * sizeof(float **));

    num_strides = (int *) malloc(FILES * sizeof(int ));
    stride_all = (int ***) malloc(FILES * sizeof(int **));
    stride_float = (float ***) malloc(FILES * sizeof(float **));

    peaks_xa = (float ***) malloc(FILES * sizeof(float **));
    peaks_ya = (float ***) malloc(FILES * sizeof(float **));
    peaks_zg = (float ***) malloc(FILES * sizeof(float **));

    troughs_xa = (float ***) malloc(FILES * sizeof(float **));
    troughs_ya = (float ***) malloc(FILES * sizeof(float **));
    troughs_zg = (float ***) malloc(FILES * sizeof(float **));


    for(i = 0; i < FILES; i++) {
        num_strides[i] = 0;
    }


    for(i = 0; i < FILES; i++)
    {
        sigma_high_xa[i] = (float *) malloc(sizeof(float) * STRIDES);
        sigma_low_xa[i] = (float *) malloc(sizeof(float) * STRIDES);
        sigma_high_ya[i] = (float *) malloc(sizeof(float) * STRIDES);
        sigma_low_ya[i] = (float *) malloc(sizeof(float) * STRIDES);
        sigma_high_zg[i] = (float *) malloc(sizeof(float) * STRIDES);
        sigma_low_zg[i] = (float *) malloc(sizeof(float) * STRIDES);

        mean_high_xa[i] = (float **) malloc(STRIDES * sizeof(float *));
        mean_low_xa[i] = (float **) malloc(STRIDES * sizeof(float *));
        mean_high_ya[i] = (float **) malloc(STRIDES * sizeof(float *));
        mean_low_ya[i] = (float **) malloc(STRIDES * sizeof(float *));
        mean_high_zg[i] = (float **) malloc(STRIDES * sizeof(float *));
        mean_low_zg[i] = (float **) malloc(STRIDES * sizeof(float *));

        peaks_xa[i] = (float **) malloc(STRIDES * sizeof(float *));
        troughs_xa[i] = (float **) malloc(STRIDES * sizeof(float *));
        peaks_ya[i] = (float **) malloc(STRIDES * sizeof(float *));
        troughs_ya[i] = (float **) malloc(STRIDES * sizeof(float *));
        peaks_zg[i] = (float **) malloc(STRIDES * sizeof(float *));
        troughs_zg[i] = (float **) malloc(STRIDES * sizeof(float *));

        stride_all[i] = (int **) malloc(STRIDES * sizeof(int *));
        stride_float[i] = (float **) malloc(STRIDES * sizeof(float *));

    }

    for(i = 0; i < FILES; i++)
    {
        for(j = 0; j < STRIDES; j++)
        {
            mean_high_xa[i][j] = (float *) malloc(SLOTS * sizeof(float ));
            mean_low_xa[i][j] = (float *) malloc(SLOTS * sizeof(float ));
            mean_high_ya[i][j] = (float *) malloc(SLOTS * sizeof(float ));
            mean_low_ya[i][j] = (float *) malloc(SLOTS * sizeof(float ));
            mean_high_zg[i][j] = (float *) malloc(SLOTS * sizeof(float ));
            mean_low_zg[i][j] = (float *) malloc(SLOTS * sizeof(float ));

            peaks_xa[i][j] = (float *) malloc(4 * sizeof(float ));
            troughs_xa[i][j] = (float *) malloc(4 * sizeof(float ));
            peaks_ya[i][j] = (float *) malloc(4 * sizeof(float ));
            troughs_ya[i][j] = (float *) malloc(4 * sizeof(float ));
            peaks_zg[i][j] = (float *) malloc(4 * sizeof(float ));
            troughs_zg[i][j] = (float *) malloc(4 * sizeof(float ));

            stride_all[i][j] = (int *) malloc(4 * sizeof(int ));
            stride_float[i][j] = (float *) malloc(STRIDES * sizeof(float));

        }
    }
    int cnt;
    cnt = 0;
    for(i = 0; i < FILES; i++)
    {
        cnt = 0;
        for(j = 1; j < count_t2[i]; j++)
        {
            if(data[i][0][strides_t2[i][j]] - data[i][0][strides_t2[i][j-1]] < 3)
            {
                stride_all[i][cnt][0] = strides_t2[i][j-1];
                stride_all[i][cnt][1] = strides_t2[i][j];
                stride_all[i][cnt][2] = strides_t2[i][j]-strides_t2[i][j-1];
                num_strides[i]++;
                cnt++;
            }
        }
        for(j = 0; j < num_strides[i]-1; j++)
        {
            if(stride_all[i][j][1] == stride_all[i][j+1][0])
            {
                stride_all[i][j][3] = 1;
            }
            else
            {
                stride_all[i][j][3] = -1;
            }
        }
        for(j = 0; j < num_strides[i]; j++)
        {
            printf("stride: %d, end: %d, period: %d, next: %d\n", stride_all[i][j][0],
                   stride_all[i][j][1], stride_all[i][j][2], stride_all[i][j][3]);
        }

    }

    for(i = 0; i < FILES; i++)
    {
        sigma_stride(data[i][1], stride_all[i], sigma[i][0], sigma_low_xa[i], sigma_high_xa[i], num_strides[i]);
        sigma_stride(data[i][2], stride_all[i], sigma[i][1], sigma_low_ya[i], sigma_high_ya[i], num_strides[i]);
        sigma_stride(data[i][6], stride_all[i], sigma[i][5], sigma_low_zg[i], sigma_high_zg[i], num_strides[i]);

        printf("before peaks\n");
        peaks(data[i][1], stride_all[i], peaks_xa[i], troughs_xa[i], num_strides[i]);
        peaks(data[i][2], stride_all[i], peaks_ya[i], troughs_ya[i], num_strides[i]);
        peaks(data[i][6], stride_all[i], peaks_zg[i], troughs_zg[i], num_strides[i]);
        for(j = 0; j < num_strides[i]-1; j++)
        {
            printf("%f\t%f\t%f\t%f\n", peaks_xa[i][j][0], peaks_xa[i][j][1], peaks_xa[i][j][2], peaks_xa[i][j][3]);
            printf("%f\t%f\t%f\t%f\n\n", troughs_xa[i][j][0], troughs_xa[i][j][1], troughs_xa[i][j][2], troughs_xa[i][j][3]);
        }
        printf("next_file\n");
//        sigma_high_xa, ** sigma_low_xa;
//        sigma_high_ya, ** sigma_low_ya;
//        sigma_high_zg, ** sigma_low_zg;
//
//        troughs_xa, *** peaks_xa;
//        troughs_ya, *** peaks_ya;
//        troughs_zg, *** peaks_zg;
//
//
//        mean_low_xa, *** mean_high_xa;
//        mean_low_ya, *** mean_high_ya;
//        mean_low_zg, *** mean_high_zg;
//         sigma_stride()
//        void sigma_stride(float * data, int ** stride, float * level, float * sigma_low, float * sigma_high, int n);
////void center_data(float * data, float * center, int n);
//
//
//
//        void mean_slot(float * data, int ** strides, float * level, float ** mean_low, float ** mean_high, float slots, int n);

    }



    printf("done strides\n");
    strides_out = (char **) malloc(FILES * sizeof(char *));
    output_files = (char **) malloc(FILES * sizeof(char *));
    parameters = (char **) malloc(FILES * sizeof(char *));
    for (i = 0; i < FILES; i++) {
        output_files[i] = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(output_files[i], 0, BUFF_SIZE);
        snprintf(output_files[i],
                 BUFF_SIZE,
                 "%s%s%s_%s_2.csv",main_path, paths_out[1], file_names[i], names[1]);
        strides_out[i] = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(strides_out[i], 0, BUFF_SIZE);
        snprintf(strides_out[i],
                 BUFF_SIZE,
                 "%s%s%s_%s_strides.csv",main_path, stride_path, file_names[i], names[1]);
        parameters[i] = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(parameters[i], 0, BUFF_SIZE);
        snprintf(parameters[i],
                 BUFF_SIZE,
                 "%s%s%s_%s.csv",main_path, class_files[1], file_names[i], names[1]);
//        printf("mem_set_out\n");

    }
    printf("files initialized\n");

    for(i = 0; i < FILES; i++) {

        printf("Attempting to write to file \'%s\'.\n", output_files[i]);
        fp = fopen(output_files[i], "w");
        if (fp == NULL) {
            fprintf(stderr,
                    "Failed to write to file \'%s\'.\n",
                    output_files[i]
            );
            exit(EXIT_FAILURE);
        }
        printf("writing data\n");
        fprintf(fp, "time, acceleration x, acceleration y, acceleration z, gyro x gyro y, gyro z\n");
        for (j = 0; j < data_lengths[i]; j++) {
            fprintf(fp, "%20.10lf,%20.10lf,%10.10lf,%10.10lf,%10.10lf,%10.10lf,%10.10lf\n",
                   data[i][0][j],
                   data[i][1][j],
                   data[i][2][j],
                   data[i][3][j],
                   data[i][4][j],
                   data[i][5][j],
                   data[i][6][j]

            );
        }
    }

    for(i = 0; i < FILES; i++) {

        printf("Attempting to write to file \'%s\'.\n", output_files[i]);
        fp = fopen(strides_out[i], "w");
        if (fp == NULL) {
            fprintf(stderr,
                    "Failed to write to file \'%s\'.\n",
                    output_files[i]
            );
            exit(EXIT_FAILURE);
        }
        printf("Writing strides\n");

        fprintf(fp, "Period,accel_x_peak1_magnitude,accel_x_peak1_location,accel_x_peak2_magnitude,accel_x_peak2_location,"
                "accel_x_trough1_magnitude,accel_x_trough1_location,accel_x_trough2_magnitude,accel_x_trough2_location,"
                "accel_y_peak1_magnitude, accel_y_peak1_location, accel_y_peak2_magnitude, accel_y_peak2_location,"
                "accel_y_trough1_magnitude,accel_y_trough1_location,accel_y_trough2_magnitude,accel_y_trough2_location,"
                "gyro_z_peak1_magnitude,gyro_z_peak1_location,gyro_z_peak2_magnitude,gyro_z_peak2_location,"
                "gyro_z_trough1_magnitude,gyro_z_trough1_location,gyro_z_trough2_magnitude,gyro_z_trough2_location,"
                "standard_deviation_peak,standard_deviation_trough\n");

        for (j = 0; j < num_strides[i]-1; j++) {
            fprintf(fp,     "%10d,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,"
                            "%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,"
                            "%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf\n",
                    stride_all[i][j][2],
                    peaks_xa[i][j][0],
                    peaks_xa[i][j][1],
                    peaks_xa[i][j][2],
                    peaks_xa[i][j][3],
                    troughs_xa[i][j][0],
                    troughs_xa[i][j][1],
                    troughs_xa[i][j][2],
                    troughs_xa[i][j][3],
                    peaks_ya[i][j][0],
                    peaks_ya[i][j][1],
                    peaks_ya[i][j][2],
                    peaks_ya[i][j][3],
                    troughs_zg[i][j][0],
                    troughs_zg[i][j][1],
                    troughs_zg[i][j][2],
                    troughs_zg[i][j][3],
                    peaks_zg[i][j][0],
                    peaks_zg[i][j][1],
                    peaks_zg[i][j][2],
                    peaks_zg[i][j][3],
                    troughs_zg[i][j][0],
                    troughs_zg[i][j][1],
                    troughs_zg[i][j][2],
                    troughs_zg[i][j][3],
                    sigma_high_xa[i][j],
                    sigma_low_xa[i][j]

            );
        }
    }
    for(i = 0; i < FILES; i++) {

        printf("Attempting to write to file \'%s\'.\n", parameters[i]);
        fp = fopen(parameters[i], "w");
        if (fp == NULL) {
            fprintf(stderr,
                    "Failed to write to file \'%s\'.\n",
                    parameters[i]
            );
            exit(EXIT_FAILURE);
        }
        printf("Writing parameters\n");

        fprintf(fp, "Period,accel_x_peak1_magnitude,accel_x_peak1_location,accel_x_peak2_magnitude,accel_x_peak2_location,"
                "accel_x_trough1_magnitude,accel_x_trough1_location,accel_x_trough2_magnitude,accel_x_trough2_location,"
                "accel_y_peak1_magnitude, accel_y_peak1_location, accel_y_peak2_magnitude, accel_y_peak2_location,"
                "accel_y_trough1_magnitude,accel_y_trough1_location,accel_y_trough2_magnitude,accel_y_trough2_location,"
                "gyro_z_peak1_magnitude,gyro_z_peak1_location,gyro_z_peak2_magnitude,gyro_z_peak2_location,"
                "gyro_z_trough1_magnitude,gyro_z_trough1_location,gyro_z_trough2_magnitude,gyro_z_trough2_location,"
                "standard_deviation_peak,standard_deviation_trough\n");

        for (j = 0; j < num_strides[i]-1; j++) {
            fprintf(fp,     "%10d,%20.10lf,%10d,%20.10lf,%10d,%20.10lf,%10d,%20.10lf,%10d,"
                            "%20.10lf,%10d,%20.10lf,%10d,%20.10lf,%10d,%20.10lf,%10d,%20.10lf,"
                            "%10d,%20.10lf,%10d,%20.10lf,%10d,%20.10lf,%10d,%20.10lf,%20.10lf\n",
                    stride_all[i][j][0],
                    peaks_xa[i][j][0],
                    (int) peaks_xa[i][j][1],
                    peaks_xa[i][j][2],
                    (int)  peaks_xa[i][j][3],
                    troughs_xa[i][j][0],
                    (int) troughs_xa[i][j][1],
                    troughs_xa[i][j][2],
                    (int) troughs_xa[i][j][3],
                    peaks_ya[i][j][0],
                    (int) peaks_ya[i][j][1],
                    peaks_ya[i][j][2],
                    (int) peaks_ya[i][j][3],
                    troughs_zg[i][j][0],
                    (int) troughs_zg[i][j][1],
                    troughs_zg[i][j][2],
                    (int) troughs_zg[i][j][3],
                    peaks_zg[i][j][0],
                    (int) peaks_zg[i][j][1],
                    peaks_zg[i][j][2],
                    (int) peaks_zg[i][j][3],
                    troughs_zg[i][j][0],
                    (int) troughs_zg[i][j][1],
                    troughs_zg[i][j][2],
                    (int) troughs_zg[i][j][3],
                    sigma_high_xa[i][j],
                    sigma_low_xa[i][j]

            );
        }
    }
}