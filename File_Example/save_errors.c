#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Functions.h"
#include <sys/time.h>
#include <math.h>

#define BUFF_SIZE 1024
#define FILES 12


int temp_strides(float **data, int *out, float *sigma_y, float *sigma_z, float max_accel,
                 float max_gyro, float max_time, float min_time, int n)
{
    float sum_ay, sum_gz, dif_time, dif_false, last, last_false;
    int i, j, cnt;
    sum_ay = 0;
    sum_gz = 0;
    cnt = 0;
    last = data[0][0];
    last_false = data[0][0];
//    printf("Strides start\n");
    for(i = 0; i < n; i++) {
        if (data[2][i] > sigma_y[i]) {
            sum_ay += data[2][i] - sigma_y[i];
//            printf("y");
        }
        if (data[6][i] < - 1.4 * sigma_z[i]) {
            sum_gz -= (data[6][i] + sigma_z[i]);
//            printf("g");
        }
        dif_time = data[0][i] - last;
        dif_false = data[0][i] - last_false;
        if (sum_gz > max_gyro || sum_ay > max_accel) {
            sum_ay = 0;
            sum_gz = 0;
            if (dif_time > max_time)
            {
//                printf("too long ");
                last = data[0][i];
                last_false = data[0][i];
            }
            else if(dif_false < min_time || sigma_y[i] < 0.05)
            {
//                printf("too short ");
                last_false = data[0][i];
            }
            else
            {
                last = data[0][i];
                last_false = data[0][i];
                out[cnt] = i;
                cnt++;
//                printf("Stride detected\n");
            }
        }
    }
    return cnt;
}
//def temp_strides(data, sigmas, max_accel, max_gyro, max_time, min_time = 0.22):
//    sum_acceleration = 0
//    sum_gyro = 0
//    last = data[0][0]
//    last_false = data[0][0]
//    out = []
//    temp = []
//    N = len(data[0])
//    # print(N)
//    for i in range(0, N):
//        # print(len(data[2][i]), len(sigmas[1][0]))
//        if data[2][i] > sigmas[1][i]:
//            sum_acceleration += data[2][i] - sigmas[1][i]
//        if data[6][i] < -sigmas[5][i]:
//            sum_gyro -= data[6][i] + sigmas[5][i]
//
//        dif_time = data[0][i] - last
//        dif_false = data[0][i] - last_false
//        if sum_gyro > max_gyro or sum_acceleration > max_accel:
//            tg = sum_gyro
//            sum_acceleration = 0
//            sum_gyro = 0
//
//            if dif_time > max_time:
//                last = data[0][i]
//                last_false = data[0][i]
//            elif dif_false < min_time or sigmas[1][i] < 0.05:
//                last_false = data[0][i]
//            elif data[6][i] > 2* sigmas[5][i]:
//                last_false = data[0][i]
//            else:
//                # print(tg)
//                out.append(i)
//                last = data[0][i]
//                last_false = data[0][i]
//
//    return out

void sample_mean(float * data, float * out, int chunk, int n)
{
    int step, rest, i, j, start, end;
    float temp_mean, weight;
    step = (int) n/chunk;
    rest = n - chunk * step;
    for( j = 0; j < step +1; j++)
    {
        temp_mean = 0;
        start = j * chunk;
        end = start + chunk;
        weight = chunk;
        if(j == step)
        {
            end = start + rest;
            weight = rest;
        }
        temp_mean = 0;
        for(i = start; i < end; i ++)
        {
            temp_mean += data[i];
        }
//        printf("chunk: %d\n", start);
        temp_mean = temp_mean / weight;
        for(i = start; i < end; i++)
        {
            out[i] = temp_mean;
//            printf("mean: %f\n", out[i]);
        }
//        printf("temp_mean: %f\n", temp_mean);
    }
}


int select_strides(float * times, int * strides, int * out, float min_time, int n)
{
    int i, count, j;
    count = 0;
    float prev_time, current_time;
    if(n < 1)
        return 0;
    out[count] = strides[0];
    count++;

    for(i = 1; i < n; i++)
    {
        current_time = times[strides[i]];
        prev_time = times[strides[i-1]];
        if(current_time - prev_time > min_time)
        {
            out[count] = strides[i];
            count++;
        }
    }
    return count;
}


void sample_sigma(float * data, float * mean, float * out, int chunk, int n)
{
    int step, rest, i, j, start, end;
    float temp_square, weight;
    step = (int) n/chunk;
    rest = n - chunk * step;
    for( j = 0; j < step +1; j++)
    {
        start = j * chunk;
        end = start + chunk;
        weight = chunk;
        if(j == step)
        {
            end = start + rest;
            weight = rest;
        }
        temp_square = 0;
        for(i = start; i < end; i ++)
        {
            temp_square += (float) pow(data[i]-mean[i], 2);
        }
        temp_square = temp_square / weight;
        temp_square = (float) sqrt(temp_square);
        for(i = start; i < end; i++)
        {
            out[i] = temp_square;
        }
    }
}
//def sample_sigma_lh(data, mean, chunk = 3000):
//    N = len(data)
//    rest = N % chunk
//    quot = N // chunk
//    sigma = np.zeros(shape=N, dtype='float')
//    for k in range(quot + 1):
//        upper_sigma = 0
//        start = k * chunk
//        end = start + chunk
//        total = chunk
//        if k == quot:
//            end = start + rest
//            total = rest
//        for i in range(start, end):
//            temp = data[i] - mean[i]
//            upper_sigma += temp * temp
//        sigma[start: end] = upper_sigma / total
//    sigma = np.sqrt(sigma)
//    return sigma


void center_data(float * data, float * mean, float * sigma, int chunk, int n)
{
    int step, rest, i, j, start, end, ratio;
    float temp_center, weight;
    step = n / chunk;
    rest = n - chunk * step;
    for( j = 0; j < step +1; j++)
    {
        temp_center = 0;
        start = j * chunk;
        end = start + chunk;
        weight = chunk;
        if(j == step)
        {
            end = start + rest;
            weight = rest;
        }
        ratio = 1;
        for(i = start; i < end; i ++)
        {
            if(mean[i] - sigma[i] < data[i] && mean[i] + sigma[i] > data[i]){
                ratio ++;
                temp_center += data[i];
            }


        }
        temp_center = temp_center / ratio;
        for(i = start; i < end; i++)
        {
            data[i] -= temp_center;
        }
    }

}

//
//def find_center(data, sigma, mean, chunk = 3000, weight=1.5):
//    N = len(data)
//    rest = N % chunk
//    quot = N // chunk
//    out = np.zeros(shape=N, dtype='float')
//    for k in range(quot):
//        temp_center = 0
//        for i in range(k * chunk, (k + 1) * chunk):
//            if mean[i] - weight*sigma[1][i] < data[i] < mean[i] + weight * sigma[0][i]:
//                temp_center += data[i]
//        out[k * chunk: (k + 1) * chunk] = temp_center / chunk
//    temp_center = 0
//    for i in range(quot * chunk, quot * chunk + rest):
//        # print(len(mean), len(sigma), len(data))
//        if mean[i] - sigma[1][i] < data[i] < mean[i] + sigma[0][i]:
//            temp_center += data[i]
//    out[(quot) * chunk: (quot) * chunk + rest] = temp_center / rest
//    return out

//void center_data(float * data, float * center, int n)
//{
//    int i;
//    for(i = 0; i < n; i++)
//    {
//        data[i] -= center[i];
//    }
//}

int peak_strides(float * data, int * strides, int *holder, float * sigma, float weight, int sn, int n)
{
    int start, end, i , j, count, count2;
    float d1, d2, d3;
    count = 0;
    count2 = 0;
//    printf("\npeaks:\n");
    for(j = 0; j < sn; j++)
    {
        start = strides[j] - 130;
        end = start + 130;
        if(start < 0)
        {
            start = 0;
        }
        if(end > n)
        {
            end = n;
        }
//        printf("start: %d\t stride: %d\t end: %d\n", start, strides[j], end);
        for(i = start+1; i < end; i++)
        {

            if(-1 * weight * sigma[i] > data[i] && data[i] > data[i-1])
            {
                holder[count] = i;
                count ++;
                break;
            }
        }
    }
//    printf("\nholder: %d\n", count);
    for(i = 0; i < sn; i++)
        strides[i] = 0;
    d1 = 0;
    d2 = 0;
    i = 1;
    strides [count2] = holder[0];
    count2 += 1;
    if(count > 1){
        d2 = data[holder[1]] - data[holder[0]];
    }
    printf("%d\n", count);
    while(i < count-1)
    {
        d1 = data[holder[i+1]] - data[holder[i]];
        if(d2 < 0.2)
        {
            i ++;
            continue;
        }
        strides[count2] = holder[i];
        count2 += 1;
        if (d1 / d2 < 0.65 && d2 < 3 && d1 < 3)
        {
//            if(i < count -2 && (data[holder[i+2]] - data[holder[i+1]])/d2 < 0.65)
//                i += 2;
//            else
//                i += 1;
            i+=2;
        }

        else
            i++;
        d2 = d1;
    }
    if (d1 / d2 < 0.65)
        strides[count2++] = holder[i];
    printf("%d\n", count2);
    return count2;
}

//
// Created by Amaael Antonini on 3/23/17.
//

