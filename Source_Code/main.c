#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "classifiers.h"
#define INPUT 27;
#define CYCLES 10
#define AUX 1

#define BUFF_SIZE 1024
const char * names_1[] = {
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
// const char *train_name = "testing.txt";
const char * test_name = "/test";
const char * test_name2 = "/Users/amaaelantonini/Google Drive/Local/Spring_2017/EE 180DB/Real_Time_trial/test";
const char * data_path = "/Users/amaaelantonini/Google Drive/Local/Spring_2017/EE 180DB/Real_Time_trial/Data_sets/AA2_data/";
const char * name = "AA2";
int main(int argc, char **argv) {

    char * fake_test;
    fake_test = (char*)malloc(sizeof(char)*BUFF_SIZE);
    memset(fake_test, 0, sizeof(char)*BUFF_SIZE);
    sprintf(fake_test, "%s%s_%s.csv", data_path, names_1[0], name);

    //char * test_file_name;
//    test_file_name = (char*)malloc(sizeof(char)*BUFF_SIZE);
//    memset(test_file_name, 0, BUFF_SIZE);

    float ***classifiers = NULL;
    // pointer to 3 d array with first D1 is activities,
    // second D2 is rows per activity,
    // third D3 is number of classifies (27 for now)
    int **counts = NULL;
    // pointer to 2 d array with first D1 is holds number of activities (12 total for now)
    // second D number of rows
    // each element is number of classifiers in each row, (constant of 27, not necessary)
    int *dimen = NULL;
    // dime is a pointer to 1 d array with first D1 number of activities (12 for now)
    // second D2 is rows per activity
    int n;
    // n is integer with number of activities (12 for now)
    int input;
    input = INPUT;
    // when training, there will be sum o n rows if dim[i] for i from 0 to n-1 rows
    // dim[i] is the truth value or element 1 in row of -1 -1 -1 1(ith) -1

    // example follows

    int i, j, k;
    float **test_array;
    int num_of_test_classifers;
    int truth_value;
//    if(access("training.txt", F_OK) == -1 || (argv[1] == "y" || argv[1] == "Y")){

    if(access("/Users/amaaelantonini/Google Drive/Local/Spring_2017/EE 180DB/Real_Time_trial/testing.txt", F_OK) == -1 ){
        printf("training\n");
        collect_training_data();
        train_parameters(&classifiers, &counts, &dimen, &n);
        training_file(classifiers, train_name, dimen, input, n);
        //train neural_network
        train_network(train_name);
    }

//    test_data(test_name2);
//     Check if file exists; -1 means file does not exists
    int loops;
    int check;
    check = 0;
    loops = 0;

    const char *test_file_name;
    test_file_name = (char*)malloc(sizeof(char)*BUFF_SIZE);
    //memset(test_file_name, 0, sizeof(char)*BUFF_SIZE);

//    while (loops < CYCLES && run_flag){
    while (loops < CYCLES ) {

        sprintf(test_file_name, "%s_%d.csv", test_name2, loops);
        // delete next line
        test_data(test_name2);

        if(access(fake_test, F_OK) == -1){
            // File does not exist
            sleep(1);
            printf("%s, sleeping\n", fake_test);
        }

        else {
            // File exists, so populate our data structure

            truth_value = test_classifiers(&test_array, fake_test, &num_of_test_classifers);
            printf("%s\t", test_file_name);
            check = remove(test_file_name);
            printf("%d\n", check);
            test_neural_network(test_array, num_of_test_classifers);
        }
        loops++;
    }


}
