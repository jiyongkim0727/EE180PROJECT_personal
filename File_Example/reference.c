//
// Created by Amaael Antonini on 2/10/17.
//

//317 lines (285 sloc)  7.11 KB
/* for file and terminal I/O */
#include <stdio.h>
/* for string manip */
#include <string.h>
/* for exit() */
#include <stdlib.h>
/* for fabsf() */
//#include <math.h>
#include "Functions.h"

#define BUFF_SIZE 1024

void reintegrate(float *time, float *integrand, float *dest, int n)
{
// integrate again the data with respect to time
// force zeros when integrand is zero
    dest[0] = 0;
    int i;
    for(i = 1; i < n; i++)
    {
        if( integrand[i] > 0)
        {
            dest[i] = dest[i-1] + integrand[i]*(time[i] - time[i-1]);
        }
        else
        {
            dest[i] = 0;
        }
    }
}

int get_strides(float *time, float *data, float *destination, int *sample, float level, int n)
{
// find dthe number of strides by
// counting the peaks higher than level
// print the points of max height
    int test = 1;
    int strides = 0;
    int i;
    int j;
    j = 0;
    for(i = 0; i < n; i++)
    {
        if(data[i] > level)
        {
            if(i == n || (data[i+1] < data[i] && test)) {
                strides++;
                test = 0;
                destination[j] = time[i];
                sample[j] = i;
                // printf("sample: %d, time: %f\n", sample[j], destination[j]);
                j++;
            }
        }
        else if(!test)
        {
            test = 1;
        }
    }
    return strides;
}
//    test = True
//    strides = 0
//    for i in range(len(data) - 1):
//        if data[i] > average:
//            if data[i + 1] < data[i] and test:
//                strides += 1
//                test = False
//                print(time[i], data[i])
//        elif not test:
//            test = True
//    print(strides)
//    return strides

int count_strides(float *time, float *data, float level, int n)
{
    // count the strides to generate right size of data type
    int test = 1;
    int strides = 0;
    int i;
    for(i = 0; i < n; i++)
    {
        if(data[i] > level)
        {
            if(i == n || (data[i+1] < data[i] && test)) {
                strides++;
                test = 0;
                printf("%f, %f\n", time[i], data[i]);
            }
        }
        else if(!test)
        {
            test = 1;
        }
    }
    printf("strides: %d\n", strides);
    return strides;
}
//
//void old_integrate(float *time, float *integrand, float *dest, float l_a, float l_v, int max, int n)
//{
//    // integrate data from integrand
//    // only use values higher than acceleration leve or lower than ints
//    // if integrand is not used and
//    // sum is within range +- velocity level
//    // force a zero in the sum
//    // if nothing added for max interval
//    // for all max sequence to zero
//    int iter(0);
//    dest[0] = 0;
//    if(l_a < 0)
//        l_a = -l_a;
//    if(l_v < 0)
//        l_v = -l_v;
//    int i;
//    float current;
//    for(i = 1; i < n; i++)
//    {
//        if(integrand[i] > l_a || integrand[i] < -l_a)
//        {
//            current = integrand[i] * (time[i] - time[i-1]);
//            dest[i] = dest[i-1] + current;
//        }
//        else if(-l_v < dest[i-1] && l_v > dest[i-1])
//        {
//            dest[i] = 0;
//            iter = 0;
//        }
//        else
//        {
//            dest[i] = dest[i-1];
//            iter++;
//        }
//        if(iter > max)
//        {
//            int j;
//            for(j = 0; j < max; j++)
//            {
//                dest[i-j] = 0;
//            }
//        }
//    }
//}

void integrate(float *time, float *integrand, float *dest, float l_a, int n)
{
    // integrate data from integrand
    // only use values higher than acceleration leve or lower than ints
    dest[0] = 0;
    if(l_a < 0)
        l_a = -l_a;
    int i;
    float current;
    for(i = 1; i < n; i++)
    {
        if(integrand[i] > l_a || integrand[i] < -l_a)
        {
            current = integrand[i] * (time[i] - time[i-1]);
            dest[i] = dest[i-1] + current;
        }
        else
        {
            dest[i] = dest[i-1];
        }
    }
}
//def integrate(time, integrand, acceleration_average):
//    out = np.empty(shape=len(time), dtype="float")
//    out[0] = 0
//    if acceleration_average < 0:
//        acceleration_average = - acceleration_average
//    for i in range(1, len(integrand)):
//        if acceleration_average < integrand[i] or integrand[i] < -acceleration_average:
//            current = integrand[i] * (time[i] - time[i - 1])
//            out[i] = out[i - 1] + current
//        else:
//            out[i] = out[i - 1]
//    return out

//void simple_integrate(float *time, float *data, float *dest, int n)
//{
//    dest[0] = 0;
//    int i;
//    for(i = 1; i < n; i++)
//    {
//        dest[i] = dest[i-1] * (time[i] - time[i-1]);
//    }
//}
//
////def simple_integrate(time, data):
////    out = np.empty(shape = len(time), dtype ="float")
////    out[0] = 0
////    for i in range(len(data)):
////        out[i] = out[i-1] + data[i] * (time[i] - time[i-1])
////    return out

void average_out(float *data, float average, int n)
{
    // subtract the average of the data from it to get net zero sum
    int i;
    for(i = 0; i < n; i++)
        data[i] -= average;
}

float get_average(float *time, float *data, int n)
{
    int i;
    float sum;
    sum = 0;
    for(i = 1; i < n; i++)
        sum += data[i] * (time[i] - time[i-1]);
    return sum / (time[n - 1] - time[0]);
}
//def average_out(data, average):
//    out = np.copy(data)
//    thd = average
//    for i in range(len(data)):
//        out[i] -= thd
//    return out
//
//void clean(float *dest, float *reference, float thd, int n)
//{
//    // copy only data outside range threshold
//    if(thd < 0)
//        thd = -thd;
//    int i;
//    for(i = 0; i < n; i++) {
//        if (-thd < reference[i] && thd > reference[i])
//            dest[i] = 0;
//        else
//            dest[i] = reference[i];
//    }
//}



//void get_half(float *data, int n)
//{
//    // eliminate the smaller half of the data
//    // make the remaining half positive
//    int plus, minus;
//    plus = 0; minus = 0;
//    int i;
//    for(i = 0; i < n; i++)
//    {
//        if(data[i] < 0)
//            minus ++;
//        else if(data[i] > 0)
//            plus ++;
//    }
//    if(plus < minus)
//    {
//        for(i = 0; i < n; i++)
//        {
//            if(data[i] > 0)
//                data[i] = 0;
//            else
//                data[i] = -data[i];
//        }
//    }
//    else
//    {
//        for(i = 0; i < n; i++)
//        {
//            if(data[i] < 0)
//                data[i] = 0;
//        }
//    }
//}

float noise_level(float *data, int n, float step, float limit)
{
    // find the noise level by getting the value whre the greatest number
    float *steps, bound;
    int *count;
    int length;
    bound = step/2;
    length = (int) ((2 * limit + 2) / (step));
    count = (int *) malloc(sizeof(int) * length);
    steps = (float *) malloc(sizeof(float) * length);
    int i, j, k;
    for(i = 0; i < length; i++)
    {
        steps[i] = -limit + i * step;
        count[i] = 0;
    }
    for(i = 0; i < n; i++)
    {
        for(j = 1; j < length; j++)
        {
            if(steps[j-1] + bound < data[i] && steps[j] + bound > data[i])
            {
                count[j]++;
                break;
            }
        }
    }
    k = 0;
    j = 0;
    for(i = 1; i < length; i++)
    {
        if(count[i]> k)
        {
            k = count[i];
            j = i;
        }
    }
    return steps[j];
}
//def noise_level(data, step, limit):
//    bound = step/2
//    temp = np.arange(-limit, limit + step, step)
//    cnt = np.zeros(shape=(len(temp)), dtype = "int")
//    for i in range(len(data)):
//        for j in range(1, len(temp)):
//            if (temp[j-1] + bound) < data[i] <= (temp[j] + bound):
//                cnt[j] += 1
//                break
//    out = 0
//    ind = 0
//    for i in range(1, len(cnt)):
//        if cnt[i] > out:
//            out = cnt[i]
//            ind = i
//    return temp[ind]

void smooth(float *data, float *arr, int length, int n)
{
    int i;
//    arr = (float *) malloc(sizeof(float) * length);
    arr[0] = data[0];
    if(n < length)
        return;
    for(i = 1; i < length; i++)
        arr[i] = (arr[i-1]*(i-1) + data[i])/((float) i);
    for(i = length; i < n; i++)
        arr[i] = arr[i-1] + (data[i] - data[i-length])/length;
}
//def smooth(data, threshold):
//    out = np.empty(shape=(len(data)), dtype = float)
//    if len(data) < threshold:
//        return
//    out[0] = 0
//    for i in range(1,threshold):
//        out[i] = out[i-1]*(i-1)/i + data[i]/i
//    for i in range(threshold, len(data)):
//        out[i] = out[i-1] +(-data[i-threshold] + data[i])/threshold
//    return out

void tail(float *data, float bar, int n)
{
    // remove the tail of the data if it does not come down
    int i, j, b;
    j = 0;
    b = 0;
    for(i = 0; i < n; i++)
    {
        if(data[i] > bar && !b)
        {
            b = 1;
            j = i;
        }
        else if(b && data[i] < bar)
        {
            b = 0;
        }
    }
    if(b)
    {
        for(i = j; i < n; i++)
        {
            data[i] = 0;
        }
    }


}
//def tail(data, a):
//    b = False
//    n = 0
//    for i in range(len(data)):
//        if data[i] > a and not b:
//            b = True
//            n = i
//        elif b and data[i] < a:
//            b = False
//    if b:
//        for i in range(n, len(data)):
//            data[i] = 0
//        print(n)
//    return data

// turn around the data if the area of interest is negative
// find wich area starts first, positive or negative
void turn(float *data, int n, float ref)
{
    int i;
    for(i = 0; i < n; i++)
    {
        if(data[i] > ref)
            return;
        else if (data[i] < - ref)
            break;
    }
    for(i = 0; i < n; i++)
    {
        data[i] = -data[i];
    }
}
//def turn(data, stop):
//    for i in range(len(data)):
//        if data[i] > stop:
//            return data
//        elif data[i] < - stop:
//            break
//    cnt = np.empty(shape=(len(data)), dtype="float")
//    for i in range(len(data)):
//        cnt[i] = - data[i]
//    return cnt

void level(float *data, float *dest,  int n, float step)
{
    // golden funtion that allows appropriate integration
    // gets an oscillating integral and fixes the lows to zero
    // or some space bellow zero if desired
    // assumes the area of interest is the positive one
    float dist;
    dist = step;
    int i;
    int up;
    up = 1;
    dest[0] = data[0] - dist;
    for(i = 1; i < n; i++)
    {
        dest[i] = data[i] - dist;
        if (up)
        {
            if(data[i] >= data[i-1])
                continue;
            else
            {
                up = 0;
            }
        }
        else
        {
            if(data[i] <= data[i-1])
                continue;
            else
            {
                dest[i - 1] = data[i - 1] - dist;
                dest[i] = data[i] - dist;
                dist = data[i-1] + step;
                up = 1;
            }
        }
    }
}
//def min_max(data, space):
//    arr = np.empty(shape=(len(data)), dtype="float")
//    maxim = 0
//    min = 0
//    distance = space
//    arr[0] = data[0]
//    up = True
//    down = False
//    for i in range(1, len(data)):
//        arr[i] = data[i] - distance
//        if up:
//            if data[i] - data[i-1] > 0:
//                continue
//            else:
//                up = False
//                down = True
//                maxim = data[i]
//        else:
//            if data[i-1] - data[i] > 0:
//                continue
//            else:
//                distance = data[i-1] + space
//                min = data[i]
//                up = True
//                down = False
//    return arr

void copy(float *data, float *dest, int n)
{
    int i;
    for(i = 0; i < n; i++)
        dest[i] = data[i];
}

void upper_half(float *data, int n, float offset)
{
    int i;
    for(i = 0; i < n; i++)
    {
        if(data[i] < 0)
        {
            data[i] = - offset;
        }
    }
}
//def upper_half(data):
//    out = np.copy(data)
//    for i in range(len(data)):
//        if data[i] < 0:
//            out[i] = 0
//    return out

void last(float *data, int cut, int n)
{
    int i;
    for(i = n - cut; i < n; i++)
    {
        data[i] = 0;
    }
}

void skip(float *time, float *data, int n, float a, float t, float offset)
{
    // set to zero or -offset the first values t time after each peak
    // Find the first value greater than a, get rise mode,
    // find when it starts falling and get fall mode,
    // when the the value falls below a start cut mode and
    // cut data t time afterwards and start over
    float window;
    int rise, cut;
    window = 0;
    rise = 0; cut = 0;
    int i;
    for(i = 0; i < n; i++)
    {
        if(cut == 1)
        {
            if(window < t)
            {
                data[i] = - offset;
                window += (time[i] - time[i-1]);
//                printf("cutting %f, %f\n, %f", time[i], data[i], window);
                continue;

            }
            else
                window = 0;
            cut = 0;
        }
        else if(data[i] >= a)
        {
            if(rise == 0)
                rise = 1;
            else
                continue;

        }
        else if(data[i] < a)
        {
            if(rise == 1)
            {
                rise = 0;
                cut = 1;
                data[i] = -offset;
            }

        }
    }

}
//def skip(time, data, a, t = 0.2):
//    t2 = 0
//    b = False
//    eliminate = False
//    for i in range(1, len(data)):
//        if eliminate:
//            if t2 < t:
//                data[i] = 0
//                t2 += time[i] - time[i-1]
//                continue
//            else:
//                eliminate = False
//                t2 = 0
//        elif data[i] >= a:
//            if not b:
//                b = True
//            else:
//                continue
//        elif data[i] < a:
//            if b:
//                b = False
//                eliminate = True
//                data[i] = 0
//                print("cut: %f" % time[i])

int find_max(int x, int y, int z)
{
    int out;
    out = x;
    if(y > x)
        out = y;
    if(z > out)
        out = z;
    return out;
}

int main(int argc, char **argv)
{
    /* Generic variables */
    int i, j, idx;
    int rv;
    /* Variables for reading file line by line */
    char *ifile_name, *strides_file_name, *process_file_x, *process_file_y, *process_file_z;
    char *vflie, *dfile, *rfile;
    FILE *fp;
    char *line = NULL;
    char buffer[1024];
    size_t len = 0;
    ssize_t read;
    int N_SAMPLES;
//    float cutoff, step;
//    cutoff = 1.5;
//    step = 0.2;

    /* Variables for storing the data and storing the return values */
    float *t, *x, *y, *z; 	// variables for data collected from input file
    float pk_threshold;	// pk-threshold value
    float **all_data;
    char *record;
    char *lines;
    /* Variables for peak-trough detection */
    int n_d; 	// number of strides from displacement
    int n_r; 	// number of strides from integrate displacement
    /* arrays for x_averaged, velocity aveaged
     * processed velocity displacement and reintegrated */
    float *x_ave, *v_ave;

    // noise, velocity, corrected, displacement, r


    // strides index, time and value both from displacement and integrate displacement
    int vx_strides, dx_strides, rx_strides;
    int vy_strides, dy_strides, ry_strides;
    int vz_strides, dz_strides, rz_strides;
    int max_x_strides, max_y_strides, max_z_strides;
    float *v_tx, *d_tx, *r_tx;
    float *v_ty, *d_ty, *r_ty;
    float *v_tz, *d_tz, *r_tz;
    int *v_nx, *d_nx, *r_nx;
    int *v_ny, *d_ny, *r_ny;
    int *v_nz, *d_nz, *r_nz;

    // data structures for x, y, z with nise corrected
    float *x_new, *y_new, *z_new;

    // data structures for velocity after integrating acceleration
    float *v_x, *v_y, *v_z;

    // data structures for velocity after correcting zeros
    float *vc_x, *vc_y, *vc_z;

    // data structures for velocity after getting upper half
    float *u_x, *u_y, *u_z;

    // temporary holders for velocity after getting upper half
    float *uc_x, *uc_y, *uc_z;

    // data structures for displacements
    float *d_x, *d_y, *d_z;

    // data structures for integrals of displacement
    float *r_x, *r_y, *r_z;

    /*
     * Check if the user entered the correct command line arguments
     * Usage:
     * ./extract_stride_data <ifile_name> <output_peaks> <output_strides>
     * 				<threshold_value_float>
     * Or
     * ./extract_stride_data
     */
    if (argc != 5)
    {
        ifile_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(ifile_name, 0, BUFF_SIZE);
        snprintf(ifile_name,
                 BUFF_SIZE,
                 "Acceleration_Walk_Dataset.csv"
        );
        strides_file_name = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(strides_file_name, 0, BUFF_SIZE);
        snprintf(strides_file_name, BUFF_SIZE, "strides.csv");
        process_file_x = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(process_file_x, 0, BUFF_SIZE);
        snprintf(process_file_x, BUFF_SIZE, "X_DATA.csv");
        process_file_y = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(process_file_y, 0, BUFF_SIZE);
        snprintf(process_file_y, BUFF_SIZE, "Y_DATA.csv");
        process_file_z = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(process_file_z, 0, BUFF_SIZE);
        snprintf(process_file_z, BUFF_SIZE, "Z_DATA.csv");
        dfile = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(dfile, 0, BUFF_SIZE);
        snprintf(dfile, BUFF_SIZE, "acceleration_distance.csv");
        rfile = (char *) malloc(sizeof(char) * BUFF_SIZE);
        memset(rfile, 0, BUFF_SIZE);
        snprintf(rfile, BUFF_SIZE, "acceleration_reintegrate.csv");
        pk_threshold = 6.7;
    }
    else {
        ifile_name = argv[1];
        strides_file_name = argv[2];
        process_file_x = argv[3];
        process_file_y = argv[4];
        process_file_z = argv[5];

        pk_threshold = atof(argv[6]);
    }

    printf("Arguments used:\n\t%s=%s\n\t%s=%s\n\t%s=%s\n\t%s=%f\n",
           "ifile_name", ifile_name,
           "ofile_peak_trough_name", strides_file_name,
           "ofile_stride_name", process_file_x,
           "peak_threshold", pk_threshold
    );

    /* open the input file */
    printf("Attempting to read from file \'%s\'.\n", ifile_name);
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
    read = getline(&line, &len, fp); //discard header of file

    /* start reading the data from the file into the data structures */
    all_data = (float**) malloc(sizeof(float**) * N_SAMPLES);
    i = 0;
    t = (float *) malloc(sizeof(float) * N_SAMPLES);
    x = (float *) malloc(sizeof(float) * N_SAMPLES);
    y = (float *) malloc(sizeof(float) * N_SAMPLES);
    z = (float *) malloc(sizeof(float) * N_SAMPLES);

    /* initialize data structures */
    x_ave = (float *) malloc(sizeof(float) * N_SAMPLES);
    v_ave = (float *) malloc(sizeof(float) * N_SAMPLES);

    // initialize all variables
    x_new = (float *) malloc(sizeof(float) * N_SAMPLES);
    y_new = (float *) malloc(sizeof(float) * N_SAMPLES);
    z_new = (float *) malloc(sizeof(float) * N_SAMPLES);

    v_x = (float *) malloc(sizeof(float) * N_SAMPLES);
    v_y = (float *) malloc(sizeof(float) * N_SAMPLES);
    v_z = (float *) malloc(sizeof(float) * N_SAMPLES);

    vc_x = (float *) malloc(sizeof(float) * N_SAMPLES);
    vc_y = (float *) malloc(sizeof(float) * N_SAMPLES);
    vc_z = (float *) malloc(sizeof(float) * N_SAMPLES);

    uc_x = (float *) malloc(sizeof(float) * N_SAMPLES);
    uc_y = (float *) malloc(sizeof(float) * N_SAMPLES);
    uc_z = (float *) malloc(sizeof(float) * N_SAMPLES);

    u_x = (float *) malloc(sizeof(float) * N_SAMPLES);
    u_y = (float *) malloc(sizeof(float) * N_SAMPLES);
    u_z = (float *) malloc(sizeof(float) * N_SAMPLES);

    d_x = (float *) malloc(sizeof(float) * N_SAMPLES);
    d_y = (float *) malloc(sizeof(float) * N_SAMPLES);
    d_z = (float *) malloc(sizeof(float) * N_SAMPLES);

    r_x = (float *) malloc(sizeof(float) * N_SAMPLES);
    r_y = (float *) malloc(sizeof(float) * N_SAMPLES);
    r_z = (float *) malloc(sizeof(float) * N_SAMPLES);

    // paramenters to find average noise
    // cut is the step between searches: from -bound, -bound + cut, - bound + 2*cut, ... , + bound
    // bound is range, search from -cut to +cut

    float cut, bound, step, offset;
    int dampen;
    step = 0.02; bound = 1.3; cut = 0.01; offset = 0.05; dampen = 100;


    // average noise levels
    float noise_x, noise_y, noise_z;

    // averages for x, y, z in velocity, distance and distance integral
    float v_bx, v_by, v_bz;
    float d_bx, d_by, d_bz;
    float r_bx, r_by, r_bz;

    while ((read = getline(&line, &len, fp)) != -1) {
        /* parse the data */
        rv = sscanf(line, "%f,%f,%f,%f\n", &t[i], &x[i], &y[i], &z[i]);
        if (rv != 4) {
            fprintf(stderr,
                    "%s %d \'%s\'. %s.\n",
                    "Failed to read line",
                    i,
                    line,
                    "Exiting"
            );
            exit(EXIT_FAILURE);
        }
        i++;
    }
    fclose(fp);
    i = 0;
    j = 0;
    fp = fopen("training_file_girl.csv","r");
    if(fp == NULL)
    {
        printf("\n file opening failed ");
        return -1 ;
    }
    printf("file_open\n");
    while((lines=fgets(buffer,sizeof(buffer),fp))!=NULL)
    {
        record = strtok(line,";");
        while(record != NULL)
        {
            printf("record : %s",record) ;    //here you can put the record into the array as per your requirement.
            all_data[i][j++] = atoi(record) ;
            record = strtok(NULL,";");
        }
        ++i ;
    }
    fclose(fp);

    /*
     * From selected thresholds,
     * find indexes of peaks
     * find indexes of troughs
     *   int dx_strides, rx_strides;
     *   float *d_n, *d_tx, *d_x, *r_n, *r_tx, *r_x;
     */
    //cut = 0.02
    //bound = 1.5
    // noise, noise attenuated,  velocity, corrected, displacement, r
    //n_x_a = noise_level(z_a, cut, bound)

    noise_x = noise_level(x, N_SAMPLES, step, bound);
    noise_y = noise_level(y, N_SAMPLES, step, bound);
    noise_z = noise_level(z, N_SAMPLES, step, bound);

    printf("noise: %f, %f, %f\n" ,noise_x, noise_y, noise_z);

    // copy data into new arrays to hold corrected noise level
    // x_new = x; y_new = y; z_new = z;

    copy(x, x_new, N_SAMPLES);
    copy(y, y_new, N_SAMPLES);
    copy(z, z_new, N_SAMPLES);

    //x_c = average_out(x_c, n_x_a)
    average_out(x_new, noise_x, N_SAMPLES);
    average_out(y_new, noise_y, N_SAMPLES);
    average_out(z_new, noise_z, N_SAMPLES);


    // zero out the last samples to avoid confusing the program
    // prevents drift from being integrated as a step
    last(x_new, dampen, N_SAMPLES);
    last(y_new, dampen, N_SAMPLES);
    last(z_new, dampen, N_SAMPLES);


    // integrate acceleration to get velocity
    //v_3 = integrate(t_a, x_c, cut)

    integrate(t, x_new, v_x, step, N_SAMPLES);
    integrate(t, y_new, v_y, step, N_SAMPLES);
    integrate(t, z_new, v_z, step, N_SAMPLES);

    // noise: x_new, attenuated: i_att,  velocity: v_i
    // if the area of interest is bellow zero flip the velocity

    //v_c = turn(v_c, 0.01)
    turn(v_x, N_SAMPLES, offset);
    turn(v_y, N_SAMPLES, offset);
    turn(v_z, N_SAMPLES, offset);


//    int q;
//    for(q = 2000; q < N_SAMPLES; q++)
//        printf("n: %d, z: %f\n", q, v_z[q]);

    // use the minimum points to force a zero on integral
    // since no net zero integral due to drift
    // corrected: vc_i, disp: d_i, integrate dist: r_i
    // the troughs will be at zero minus some offset
    //vc_i = leve(v_c, 0.005)
    level(v_x, vc_x, N_SAMPLES, cut);
    level(v_y, vc_y, N_SAMPLES, cut);
    level(v_z, vc_z, N_SAMPLES, cut);

    // since y axis requires differet process
    // This portion helps avoid work of considering other factors
    // the y axis has 3 peaks together and it woud take more process
    // if instead of even it out something else was done
    smooth(vc_y, uc_y, dampen, N_SAMPLES);
    level(uc_y, vc_y, N_SAMPLES, cut);
    smooth(vc_y, uc_y, dampen, N_SAMPLES);
    level(uc_y, vc_y, N_SAMPLES, cut);

    // get the upper half of the velocity since the rest is not relevant
    // troughs have been already force an offset bellow zero
    // corrected: vc_i, upper half velocity: u_i, disp: d_i, integrate dist: r_i
    //v_cx = upper_half(v_c)


    copy(vc_x, uc_x, N_SAMPLES);
    copy(vc_y, uc_y, N_SAMPLES);
    copy(vc_z, uc_z, N_SAMPLES);

    upper_half(uc_x, N_SAMPLES, 0.0);
    upper_half(uc_y, N_SAMPLES, 0.05);
    upper_half(uc_z, N_SAMPLES, 0.0);

    // get the average of the velocity for 2 purpuses
    // use it as reference to find velocity strides
    // use a ratio of it to eliminate tail peaks that have no troughs

    //v_ax = get_average(t_a, v_cx)

    // even out noise by subtracting to each sample the average of the last i samples


    smooth(uc_x, u_x, dampen, N_SAMPLES);
    smooth(uc_z, u_z, dampen, N_SAMPLES);

    level(u_x, uc_x, N_SAMPLES, cut);
    level(u_z, uc_z, N_SAMPLES, cut);

//    int q;
//    for(q = 2000; q < N_SAMPLES; q++)
//        printf("n: %d, z: %f\n", q, uc_z[q]);

    copy(uc_x, u_x, N_SAMPLES);
    copy(uc_y, u_y, N_SAMPLES);
    copy(uc_z, u_z, N_SAMPLES);

    // again get only the area of interest (upper)
    upper_half(u_x, N_SAMPLES, 0.0);
    upper_half(u_z, N_SAMPLES, 0.0);

//    int q;
//    for(q = 2000; q < N_SAMPLES; q++)
//        printf("n: %d, z: %f\n", q, v_z[q]);


    v_bx = get_average(t, u_x, N_SAMPLES);
    v_by = get_average(t, u_y, N_SAMPLES);
    v_bz = get_average(t, u_z, N_SAMPLES);

    // cut tail if it it does not come down
    // raising tails are assumed to come from drift

    //v_cx = tail(v_cx, v_ax / 4)
    tail(u_x, v_bx/4, N_SAMPLES);
    tail(u_y, v_by/4, N_SAMPLES);
    tail(u_z, v_bz/4, N_SAMPLES);

    skip(t, u_y, N_SAMPLES, 2*v_by, 0.2, 0);


    // get displacement by integrating velocity
    // reintegrate function is not cumulative
    // it restarts evert time the integrand is zero
    // integrate each interval separately

    // corrected: vc_i, upper half velocity: u_i, disp: d_i, integrate dist: r_i
    //d_i = reintegrate(t_a, v_cx)
    reintegrate(t, u_x, d_x, N_SAMPLES);
    reintegrate(t, u_y, d_y, N_SAMPLES);
    reintegrate(t, u_z, d_z, N_SAMPLES);

//    reintegrate(t, d_x, u_x, N_SAMPLES);
//    reintegrate(t, d_y, u_y, N_SAMPLES);
//    reintegrate(t, d_z, u_z, N_SAMPLES);

    // integrate displacement to find result with noise attenuated
    // noise effect is assumed to contribute with less area
    // will tend to disappear when area is found

    //r_x = reintegrate(t_a, d_x)
    reintegrate(t, d_x, r_x, N_SAMPLES);
    reintegrate(t, d_y, r_y, N_SAMPLES);
    reintegrate(t, d_z, r_z, N_SAMPLES);

    // find the average to use as reference to count strides

    //d_ax = get_average(t_a, d_x)
    //r_ax = get_average(t_a, r_x)
    d_bx = get_average(t, d_x, N_SAMPLES);
    d_by = get_average(t, d_y, N_SAMPLES);
    d_bz = get_average(t, d_z, N_SAMPLES);

    r_bx = get_average(t, r_x, N_SAMPLES);
    r_by = get_average(t, r_y, N_SAMPLES);
    r_bz = get_average(t, r_z, N_SAMPLES);


    printf("PROCESSING X_AXIS\n");
    printf("velocity strides\n");
    vx_strides = count_strides(t, u_x, v_bx, N_SAMPLES);
    printf("displacement strides\n");
    dx_strides = count_strides(t, d_x, d_bx, N_SAMPLES);
    printf("reintegrate strides\n");
    rx_strides = count_strides(t, r_x, r_bx, N_SAMPLES);

    printf("PROCESSING Y_AXIS\n");
    printf("velocity strides\n");
    vy_strides = count_strides(t, u_y, v_by, N_SAMPLES);
    printf("displacement strides\n");
    dy_strides = count_strides(t, d_y, d_by, N_SAMPLES);
    printf("reintegrate strides\n");
    ry_strides = count_strides(t, r_y, r_by, N_SAMPLES);

    printf("PROCESSING Z_AXIS\n");
    printf("velocity strides\n");
    vz_strides = count_strides(t, u_z, v_bz, N_SAMPLES);
    printf("displacement strides\n");
    dz_strides = count_strides(t, d_z, d_bz, N_SAMPLES);
    printf("reintegrate strides\n");
    rz_strides = count_strides(t, r_z, r_bz, N_SAMPLES);

    max_x_strides = find_max(vx_strides, dx_strides, rz_strides);
    max_y_strides = find_max(vy_strides, dy_strides, ry_strides);
    max_z_strides = find_max(vz_strides, dz_strides, rz_strides);

    v_tx = (float *) malloc(sizeof(float) * vx_strides);
    d_tx = (float *) malloc(sizeof(float) * dx_strides);
    r_tx = (float *) malloc(sizeof(float) * rx_strides);

    v_ty = (float *) malloc(sizeof(float) * vy_strides);
    d_ty = (float *) malloc(sizeof(float) * dy_strides);
    r_ty = (float *) malloc(sizeof(float) * ry_strides);

    v_tz = (float *) malloc(sizeof(float) * vz_strides);
    d_tz = (float *) malloc(sizeof(float) * dz_strides);
    r_tz = (float *) malloc(sizeof(float) * rz_strides);

    v_nx = (int *)  malloc(sizeof(int) * max_x_strides);
    d_nx = (int *) malloc(sizeof(int) * max_x_strides);
    r_nx = (int *) malloc(sizeof(int) * max_x_strides);

    v_ny = (int *)  malloc(sizeof(int) * max_y_strides);
    d_ny = (int *) malloc(sizeof(int) * max_y_strides);
    r_ny = (int *) malloc(sizeof(int) * max_y_strides);

    v_nz = (int *)  malloc(sizeof(int) * max_z_strides);
    d_nz = (int *) malloc(sizeof(int) * max_z_strides);
    r_nz = (int *) malloc(sizeof(int) * max_z_strides);


    get_strides(t, u_x, v_tx, v_nx, v_bx, N_SAMPLES);
    get_strides(t, d_x, d_tx, d_nx, d_bx, N_SAMPLES);
    get_strides(t, r_x, r_tx, r_nx, r_bx, N_SAMPLES);

    get_strides(t, u_y, v_ty, v_ny, v_by, N_SAMPLES);
    get_strides(t, d_y, d_ty, d_ny, d_by, N_SAMPLES);
    get_strides(t, r_y, r_ty, r_ny, r_by, N_SAMPLES);

    get_strides(t, u_z, v_tz, v_nz, v_bz, N_SAMPLES);
    get_strides(t, d_z, d_tz, d_nz, d_bz, N_SAMPLES);
    get_strides(t, r_z, r_tz, r_nz, r_bz, N_SAMPLES);

    printf("Averages:\naceleration: %f\nvelocity: %f\ndisplacement: %f\nreintegrate: %f\n",
           noise_x, v_bx, d_bx, r_by);

//    for(i = 0; i < dx_strides; i++)
//        d_tx[i] = x[d_nx[i]];
//    for(i = 0; i < rx_strides; i++)
//        r_tx[i] = x[r_nx[i]];

    if (dx_strides < 0 || rx_strides < 0) {
        fprintf(stderr, "find_peaks_and_troughs failed\n");
        exit(EXIT_FAILURE);
    }

    /* DO NOT MODIFY ANYTHING BEFORE THIS LINE */

    /*
     * Insert your algorithm to convert from a series of peak-trough
     * indicies, to a series of indicies that indicate the start
     * of a stride.
     */



    /* DO NOT MODIFY ANYTHING AFTER THIS LINE */

    /* open the output file to write stride data from distance and integrate distance */
    printf("Attempting to write to file \'%s\'.\n", strides_file_name);
    fp = fopen(strides_file_name, "w");
    if (fp == NULL) {
        fprintf(stderr,
                "Failed to write to file \'%s\'.\n",
                strides_file_name
        );
        exit(EXIT_FAILURE);
    }
    printf("X data\n");
    printf("Velocity: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < vx_strides; i++) {
        printf("%d,%20.10lf\n",
               v_nx[i],
               x[v_nx[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                v_nx[i],
                x[v_nx[i]]);
    }
    printf("Distance: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < dx_strides; i++) {
        printf("%d,%20.10lf\n",
               d_nx[i],
               x[d_nx[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                d_nx[i],
                x[d_nx[i]]);
    }
    printf("R Values: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < rx_strides; i++) {
        printf("%d,%20.10lf\n",
               r_nx[i],
               x[r_nx[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                r_nx[i],
                x[r_nx[i]]);
    }
    printf("Y DATA\n");
    printf("Velocity: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < vy_strides; i++) {
        printf("%d,%20.10lf\n",
               v_ny[i],
               y[v_ny[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                v_ny[i],
                y[v_ny[i]]);
    }
    printf("Distance: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < dy_strides; i++) {
        printf("%d,%20.10lf\n",
               d_ny[i],
               y[d_ny[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                d_ny[i],
                y[d_ny[i]]);
    }
    printf("R Values: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < ry_strides; i++) {
        printf("%d,%20.10lf\n",
               r_ny[i],
               y[r_ny[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                r_ny[i],
                y[r_ny[i]]);
    }
    printf("Z_DATA");
    printf("Velocity: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < vz_strides; i++) {
        printf("%d,%20.10lf\n",
               v_nz[i],
               z[v_nz[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                v_nz[i],
                z[v_nz[i]]);
    }
    printf("Distance: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < dz_strides; i++) {
        printf("%d,%20.10lf\n",
               d_nz[i],
               z[d_nz[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                d_nz[i],
                z[d_nz[i]]);
    }
    printf("R Values: sample index, sample time\n");
    fprintf(fp, "sample index, sample time\n");
    for (i = 0; i < rz_strides; i++) {
        printf("%d,%20.10lf\n",
               r_nz[i],
               z[r_nz[i]]);

        fprintf(fp, "%d,%20.10lf\n",
                r_nz[i],
                z[r_nz[i]]);
    }
    fclose(fp);

    /* open the output file to write velocity, displacement and integrate displacement */
    printf("Attempting to write to file \'%s\'.\n", process_file_x);
    fp = fopen(process_file_x, "w");
    if (fp == NULL) {
        fprintf(stderr,
                "Failed to write to file \'%s\'.\n",
                process_file_x
        );
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "time, acceleration, velocity, velocity corrected, velocity turned,"
            " upper half,  displacement, displacement reintegrated\n");
    for (i = 0; i < N_SAMPLES; i++)
    {
        fprintf(fp, "%lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf\n",
                t[i],
                x_new[i],
                v_x[i],
                vc_x[i],
                uc_x[i],
                u_x[i],
                d_x[i],
                r_x[i]
        );
    }
    fclose(fp);

    // Y axis
    /* open the output file to write velocity, displacement and integrate displacement */
    printf("Attempting to write to file \'%s\'.\n", process_file_y);
    fp = fopen(process_file_y, "w");
    if (fp == NULL) {
        fprintf(stderr,
                "Failed to write to file \'%s\'.\n",
                process_file_y
        );
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "time, acceleration, velocity, velocity corrected, velocity turned,"
            " upper half,  displacement, displacement reintegrated\n");
    for (i = 0; i < N_SAMPLES; i++)
    {
        fprintf(fp, "%lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf\n",
                t[i],
                y_new[i],
                v_y[i],
                vc_y[i],
                uc_y[i],
                u_y[i],
                d_y[i],
                r_y[i]
        );
    }
    fclose(fp);

    // Z axis data
    /* open the output file to write velocity, displacement and integrate displacement */
    printf("Attempting to write to file \'%s\'.\n", process_file_z);
    fp = fopen(process_file_z, "w");
    if (fp == NULL) {
        fprintf(stderr,
                "Failed to write to file \'%s\'.\n",
                process_file_z
        );
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "time, acceleration, velocity, velocity corrected, velocity turned,"
            " upper half,  displacement, displacement reintegrated\n");
    for (i = 0; i < N_SAMPLES; i++)
    {
        fprintf(fp, "%lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf,%20.10lf\n",
                t[i],
                z_new[i],
                v_z[i],
                vc_z[i],
                uc_z[i],
                u_z[i],
                d_z[i],
                r_z[i]
        );
    }
    fclose(fp);



    printf("extract_stride_data completed successfuly. Exiting.\n");

    return 0;
}


//
// Created by Amaael Antonini on 3/22/17.
//


//            k = 0;
//            rv = sscanf(line, "%f,%f,%f,%f,%f,%f,%f,%f\n", &t, &t2, &ax, &ay, &az, &gx, &gy, &gz);
//            data[i][k++][j] = t;
//            data[i][k++][j] = ax;
//            data[i][k++][j] = ay;
//            data[i][k++][j] = az;
//            data[i][k++][j] = gx;
//            data[i][k++][j] = gy;
//            data[i][k++][j] = gz;
//            j++;
//            printf("%f,%f,%f,%f,%f,%f,%f,%f\n", t, t2, ax, ay, az, gx, gy, gz);
//            printf("time in: %20.10f\n", data[i][0][j]);
////            printf("line: %s\n", line);

