#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "pth_msort.h"

#define RESTRICT /* C89 doesn't have restrict */

/* ========= Shared Helper Functions & Types ========= */

typedef struct {
    int* source;
    int* dest;
    int* temp_a; 
    int* temp_b; 
    size_t start, end, start2, end2, dest_offset;
    size_t len; // Added for clarity in mode2, though start is used
} ThreadArgs;

void serial_merge(const int* RESTRICT source, int* RESTRICT dest, size_t s1, size_t e1, size_t s2, size_t e2, size_t doff) {
    size_t i = s1, j = s2, k = doff;
    while (i <= e1 && j <= e2) {
        if (source[i] <= source[j]) dest[k++] = source[i++];
        else dest[k++] = source[j++];
    }
    while (i <= e1) dest[k++] = source[i++];
    while (j <= e2) dest[k++] = source[j++];
}

void serial_msort_recursive(int* RESTRICT source, int* RESTRICT dest, size_t left, size_t right) {
    size_t mid;
    if (right <= left) return;
    mid = left + (right - left) / 2;
    serial_msort_recursive(dest, source, left, mid);
    serial_msort_recursive(dest, source, mid + 1, right);
    serial_merge(source, dest, left, mid, mid + 1, right, left);
}

void optimized_serial_sort(int* RESTRICT arr_in, int* RESTRICT arr_out, size_t size) {
    int depth = 0;
    if (size > 1) {
        depth = floor(log2(size - 1)) + 1;
    }
    if (depth % 2 == 1) {
        serial_msort_recursive(arr_in, arr_out, 0, size - 1);
        memcpy(arr_in, arr_out, size * sizeof(int));
    } else {
        serial_msort_recursive(arr_in, arr_in, 0, size - 1);
    }
}

size_t binary_search(int val, int* arr, size_t len) {
    size_t low = 0, high = len, mid;
    while (low < high) {
        mid = low + (high - low) / 2;
        if (arr[mid] < val) low = mid + 1;
        else high = mid;
    }
    return low;
}

void* sort_worker_mode1(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    size_t size = args->end - args->start + 1;
    optimized_serial_sort(args->temp_a + args->start, args->temp_b + args->start, size);
    return NULL;
}

void* sort_worker_mode2(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    /* FIX: Use 'start' which holds the length in this mode */
    optimized_serial_sort(args->temp_a, args->temp_b, args->start);
    return NULL;
}

void* merge_worker(void* arg) {
    ThreadArgs* m = (ThreadArgs*)arg;
    /* FIX: Use correct member names from the struct */
    serial_merge(m->source, m->dest, m->start, m->end, m->start2, m->end2, m->dest_offset);
    return NULL;
}

/* ========= The Main Function with Two Modes ========= */
void mergeSortParallel(const int* values, unsigned int N, int* sorted) {
    // Determine which mode to use. 2^29 is 536870912.
    if (N > 536870912) {
        /*
         * MODE 2: M>=30 Memory-Safe Strategy ("2-then-2")
         */
        int* workspace;
        int* scratch_half;
        pthread_t threads[4];
        ThreadArgs sort_args[4];
        ThreadArgs merge_args[4];
        size_t qsz, half_scratch_sz, half_size, chunk_size;
        size_t split_points[3];
        int i, pivot;

        qsz = N / 4;
        half_scratch_sz = qsz * 2;
        scratch_half = (int*)malloc(half_scratch_sz * sizeof(int));
        if (!scratch_half) { exit(1); }

        workspace = (int*)values;
        memcpy(sorted, values, (size_t)N * sizeof(int));

        sort_args[0] = (ThreadArgs){NULL, NULL, sorted,         scratch_half,       qsz, 0, 0, 0, 0, 0};
        sort_args[1] = (ThreadArgs){NULL, NULL, sorted + qsz,   scratch_half + qsz, qsz, 0, 0, 0, 0, 0};
        pthread_create(&threads[0], NULL, sort_worker_mode2, &sort_args[0]);
        pthread_create(&threads[1], NULL, sort_worker_mode2, &sort_args[1]);
        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);

        sort_args[2] = (ThreadArgs){NULL, NULL, sorted + 2 * qsz, scratch_half,       qsz, 0, 0, 0, 0, 0};
        sort_args[3] = (ThreadArgs){NULL, NULL, sorted + 3 * qsz, scratch_half + qsz, N - (3 * qsz), 0, 0, 0, 0, 0};
        pthread_create(&threads[2], NULL, sort_worker_mode2, &sort_args[2]);
        pthread_create(&threads[3], NULL, sort_worker_mode2, &sort_args[3]);
        pthread_join(threads[2], NULL);
        pthread_join(threads[3], NULL);
        
        free(scratch_half);

        half_size = N / 2;
        merge_args[0] = (ThreadArgs){sorted, workspace, NULL, NULL, 0, qsz - 1, qsz, half_size - 1, 0, 0};
        merge_args[1] = (ThreadArgs){sorted, workspace, NULL, NULL, half_size, half_size + qsz - 1, half_size + qsz, N - 1, half_size, 0};
        pthread_create(&threads[0], NULL, merge_worker, &merge_args[0]);
        pthread_create(&threads[1], NULL, merge_worker, &merge_args[1]);
        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);

        chunk_size = half_size / 4;
        for (i = 0; i < 3; i++) {
            pivot = workspace[(i + 1) * chunk_size];
            split_points[i] = binary_search(pivot, workspace + half_size, N - half_size);
        }
        merge_args[0] = (ThreadArgs){workspace, sorted, NULL, NULL, 0, chunk_size - 1, half_size, half_size + split_points[0] - 1, 0, 0};
        merge_args[1] = (ThreadArgs){workspace, sorted, NULL, NULL, chunk_size, 2 * chunk_size - 1, half_size + split_points[0], half_size + split_points[1] - 1, chunk_size + split_points[0], 0};
        merge_args[2] = (ThreadArgs){workspace, sorted, NULL, NULL, 2 * chunk_size, 3 * chunk_size - 1, half_size + split_points[1], half_size + split_points[2] - 1, 2 * chunk_size + split_points[1], 0};
        merge_args[3] = (ThreadArgs){workspace, sorted, NULL, NULL, 3 * chunk_size, half_size - 1, half_size + split_points[2], N - 1, 3 * chunk_size + split_points[2], 0};

        for (i = 0; i < 4; i++) { pthread_create(&threads[i], NULL, merge_worker, &merge_args[i]); }
        for (i = 0; i < 4; i++) { pthread_join(threads[i], NULL); }

    } else {
        /*
         * MODE 1: M < 30 Fast Strategy (4-way parallel sort)
         */
        int* temp_a;
        int* temp_b;
        pthread_t threads[4];
        ThreadArgs args[4];
        size_t quarter_size, half_size, chunk_size;
        size_t split_points[3];
        int i, pivot;

        temp_a = (int*)malloc(N * sizeof(int));
        temp_b = (int*)malloc(N * sizeof(int));
        if (!temp_a || !temp_b) { exit(1); }
        memcpy(temp_a, values, N * sizeof(int));
        
        quarter_size = N / 4;

        for (i = 0; i < 4; i++) {
            args[i].temp_a = temp_a;
            args[i].temp_b = temp_b;
            args[i].start = i * quarter_size;
            args[i].end = (i == 3) ? (N - 1) : ((i + 1) * quarter_size - 1);
            pthread_create(&threads[i], NULL, sort_worker_mode1, &args[i]);
        }
        for (i = 0; i < 4; i++) { pthread_join(threads[i], NULL); }

        half_size = N / 2;
        args[0] = (ThreadArgs){temp_a, temp_b, NULL, NULL, 0, quarter_size - 1, quarter_size, half_size - 1, 0, 0};
        args[1] = (ThreadArgs){temp_a, temp_b, NULL, NULL, half_size, half_size + quarter_size - 1, half_size + quarter_size, N - 1, half_size, 0};
        pthread_create(&threads[0], NULL, merge_worker, &args[0]);
        pthread_create(&threads[1], NULL, merge_worker, &args[1]);
        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);

        chunk_size = half_size / 4;
        for (i = 0; i < 3; i++) {
            pivot = temp_b[(i + 1) * chunk_size];
            split_points[i] = binary_search(pivot, temp_b + half_size, N - half_size);
        }

        args[0] = (ThreadArgs){temp_b, sorted, NULL, NULL, 0, chunk_size - 1, half_size, half_size + split_points[0] - 1, 0, 0};
        args[1] = (ThreadArgs){temp_b, sorted, NULL, NULL, chunk_size, 2 * chunk_size - 1, half_size + split_points[0], half_size + split_points[1] - 1, chunk_size + split_points[0], 0};
        args[2] = (ThreadArgs){temp_b, sorted, NULL, NULL, 2 * chunk_size, 3 * chunk_size - 1, half_size + split_points[1], half_size + split_points[2] - 1, 2 * chunk_size + split_points[1], 0};
        args[3] = (ThreadArgs){temp_b, sorted, NULL, NULL, 3 * chunk_size, half_size - 1, half_size + split_points[2], N - 1, 3 * chunk_size + split_points[2], 0};
        
        for (i = 0; i < 4; i++) { pthread_create(&threads[i], NULL, merge_worker, &args[i]); }
        for (i = 0; i < 4; i++) { pthread_join(threads[i], NULL); }

        free(temp_a);
        free(temp_b);
    }
}