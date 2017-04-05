#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Functions.h"
#include <sys/time.h>

#define WAIT 30
#define RUN 30
#define FILES 12
#define IDLE 8
#define OFFSET 1.7

const char * walk_modes[] = {
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
                            "high_jump"
};

//void Functions::Read_file() {
//
//    int i, j, rv;
//    int walk_speed;
//    int *output;
//    int N_SAMPLES;
//    float **all_data;
//    char *ifile_name;
//    FILE *fp;
//    char *line = NULL;
//    char buffer[1024];
//    char *record, *csv_line;
//    size_t len = 0;
//    ssize_t read;
//
//    ifile_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
//    memset(ifile_name, 0, BUFF_SIZE);
//    snprintf(ifile_name,
//             BUFF_SIZE,
//             "/Users/amaaelantonini/Documents/Winter 2017/EE 180DA/File_Example/testing_file_girl_2"
//    );
//    fp = fopen(ifile_name, "r");
//    if (fp == NULL) {
//        fprintf(stderr,
//                "Failed to read from file \'%s\'.\n",
//                ifile_name
//        );
//        exit(EXIT_FAILURE);
//    }
//
//    read = getline(&line, &len, fp); //discard header of file
//    N_SAMPLES = 0;
//    while ((read = getline(&line, &len, fp)) != -1) {
//        N_SAMPLES++;
//    }
//    rewind(fp);
//
//    all_data = (float **) malloc(sizeof(float **) * N_SAMPLES);
//    for (i = 0; i < N_SAMPLES; i++) {
//        all_data[i] = (float *) malloc(sizeof(float) * num_input);
//        for (j = 0; j < num_input; j++)
//            printf("%f ", all_data[i][j]);
//        printf("next: %d\n", i);
//    }
//    read = getline(&line, &len, fp); //discard header of file
//    j = 0;
//    i = 0;
//    while ((read = getline(&line, &len, fp)) != -1) {
//        j = 0;
//        printf("%s\n", line);
//        record = strtok(line, ",");
//        printf("%s, %d\n", record, i);
//        while (record != NULL) {
//            printf("record : %s\n", record);
//            printf("converted : %f\n", atof(record));
//            all_data[i][j++] = (float) atof(record);
//            record = strtok(NULL, ",");
//        }
//        ++i;
//    }
//    fclose(fp);
//}

void fake_data(const char **names) {
    struct timeval current;
    double start, time_wait, over;
    int files, T1, T2;
    over = 0;
    time_wait = 0;
    files = 0;
    gettimeofday(&current, 0);
    start = current.tv_sec;
//    printf("start: %d, elapsed: %d\n", start, time_wait);
    while(time_wait < IDLE)
    {
        gettimeofday(&current, 0);
        time_wait =  current.tv_sec;
        time_wait -= start;
        usleep(1000000);
        printf("idling: %d\n", IDLE - (int) time_wait);
    }
    T1 = 0;
    T2 = 0;
    gettimeofday(&current, 0);
    start = current.tv_sec;
    gettimeofday(&current, 0);
    time_wait = current.tv_sec - start;
    while(files < FILES)
    {
        T1 += WAIT;
        T2 += RUN + WAIT;
        while(time_wait < T1)
        {
            printf("next %s\tin:\t%d\n", names[files], (T1 - (int) time_wait));
            gettimeofday(&current, 0);
            time_wait = current.tv_sec - start;
            usleep(1000000);

        }
        while(time_wait < T2)
        {
            printf("collecting\t%s\t%d\n",names[files],  (int) (T2 -  time_wait));
            gettimeofday(&current, 0);
            time_wait = current.tv_sec;
            time_wait -= start;
            usleep(1000000);
        }
        files++;
        T1 += RUN;
//        T2 += WAIT + RUN;

    }
}


