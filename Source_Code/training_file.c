//
// Created by Amaael Antonini on 4/11/17.
//

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "classifiers.h"
#define BUFFER 1024

void training_file(float ***classifiers, const char * file_name, int *row_counts, int inputs, int n)
{
    char * direct = "/Users/amaaelantonini/Google Drive/Local/Spring_2017/EE 180DB/Real_Time_trial/";
    int total, i, j, k;
    char * train_file;
    train_file = (char *) malloc(sizeof(char) * BUFFER);
    memset(train_file, 0, BUFFER);
    sprintf(train_file, "%s%s", direct, file_name);
    total = 0;
    for(i = 0; i< n; i++)
        total+= row_counts[i];
    FILE *fp;

    fp = fopen(train_file, "w");
    if (fp == NULL) {
        fprintf(stderr,
                "Failed to write to file \'%s\'.\n",
                train_file
        );
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%d\t%d\t%d\n", total, inputs, n);
    printf("header\t");
    for(i = 0; i < n; i++) {
        for (j = 0; j < row_counts[i]; j++) {
//            printf("%d\t %d\t %d\n", n, dimen[i], counts[i][j]);
            for (k = 0; k < inputs; k++) {
                fprintf(fp, "%f\t", (classifiers[i][j][k]< 1024 && -1024< classifiers[i][j][k]?  classifiers[i][j][k]/1024: 0.5));
            }
            fprintf(fp, "\n");
            for(k = 0; k < n; k++){
                fprintf(fp, "%d\t", (k == i? 1: -1));
            }
            fprintf(fp, "\n");

        }
    }
    printf("end, %s, \n", train_file);
    fclose(fp);
}

