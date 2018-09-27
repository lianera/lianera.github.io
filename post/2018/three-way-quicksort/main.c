#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

int partition2(float a[], int p, int q)
{
    int i = p - 1;
    float x = a[q];
    for (int j = p; j < q; j++) {
        if (a[j] < x) {
            i++;
            float t = a[j];
            a[j] = a[i];
            a[i] = t;
        }
    }
    i++;
    a[q] = a[i];
    a[i] = x;
    return i;
}

void quick(float a[], int p, int q)
{
    if (p >= q)
        return;
    int r = partition2(a, p, q);
    quick(a, p, r - 1);
    quick(a, r + 1, q);
}

void partition3(float a[], int p, int q, int* less, int* greater)
{
    int lt = p - 1;
    int gt = q + 1;
    float x = a[p];
    int i = p;
    while(i < gt){
        float t = a[i];
        if (a[i] > x) { // swap to right
            gt--;
            a[i] = a[gt];
            a[gt] = t;
            continue;
        }
        if (a[i] < x) { // swap to left
            lt++;
            a[i] = a[lt];
            a[lt] = t;
        }
        i++;
    }
    *less = lt;
    *greater = gt;
}

void quick3way(float a[], int p, int q)
{
    if (p >= q)
        return;
    // random choose pivot
    int k = rand() % (q - p + 1) + p;
    // swap to left 
    float x = a[k];
    a[k] = a[p];
    a[p] = x;
    int lt, gt;
    partition3(a, p, q, &lt, &gt);
    quick3way(a, p, lt);
    quick3way(a, gt, q);
}

void sort_compare(int N, float duprate)
{
    float* a = malloc(N * sizeof(float));
    float* b = malloc(N * sizeof(float));
    // generate random numbers
    for (int i = 0; i < N; i++) {
        a[i] = (float)rand() / RAND_MAX;
        b[i] = a[i];
    }
    // put same numbers
    int sn = (int)(N*duprate);
    float x = (float)rand() / RAND_MAX;
    for (int i = 0; i < sn; i++) {
        int p = rand() % N;
        a[p] = x;
        b[p] = x;
    }

    clock_t t1 = clock();
    quick(a, 0, N - 1);
    clock_t t2 = clock();
    quick3way(b, 0, N - 1);
    clock_t t3 = clock();

    int ms1 = 1000 * (t2 - t1) / CLOCKS_PER_SEC;
    int ms2 = 1000 * (t3 - t2) / CLOCKS_PER_SEC;
    //printf("quick sort: %dms\n", t1);
    //printf("3 way quick sort: %dms\n", t2);
   // printf("%d, %d;\n", ms1, ms2);

    for (int i = 0; i < N; i++) {
        assert(a[i] == b[i]);
    }
    free(a);
    free(b);
}

int main()
{
    int N = 1000;
    for (int i = 0; i < 100; i++) {
        float duprate = (float)i / 99;
        sort_compare(N, duprate);
    }
    return 0;
}