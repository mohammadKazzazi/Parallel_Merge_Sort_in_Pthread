#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "pth_msort.h"

#if defined(__GNUC__)
#  define RESTRICT __restrict__
#  define PREFETCH_R(p) __builtin_prefetch((p), 0, 3)
#else
#  define RESTRICT
#  define PREFETCH_R(p) ((void)0)
#endif

/* Tunables */
#define SMALL_RANGE_THRESHOLD       16    /* switch to simpler method on tiny spans */
#define TAIL_COPY_THRESHOLD         16    /* decide between loop copy vs memcpy */

typedef struct {
    int* source;
    int* dest;
    int* temp_a;
    int* temp_b;
    size_t start, end, start2, end2, dest_offset;
    size_t len;
} ThreadArgs;

/* Small manual copy for short tails */
static inline void tiny_copy_int(int* RESTRICT dst, const int* RESTRICT src, size_t n) {
    size_t k;
    for (k = 0; k < n; ++k) dst[k] = src[k];
}

/* Lightweight sorter for very small ranges */
static inline void small_range_sort(int* arr, size_t lo, size_t hi) {
    size_t i;
    for (i = lo + 1; i <= hi; ++i) {
        int key = arr[i];
        size_t j = i;
        while (j > lo && arr[j - 1] > key) {
            arr[j] = arr[j - 1];
            --j;
        }
        arr[j] = key;
    }
}

/* Merge two sorted runs */
void serial_merge(const int* RESTRICT source, int* RESTRICT dest,
                  size_t s1, size_t e1, size_t s2, size_t e2, size_t doff)
{
    const int* a = source + s1;
    const int* a_end = source + e1 + 1;
    const int* b = source + s2;
    const int* b_end = source + e2 + 1;
    int* out = dest + doff;

    while (a < a_end && b < b_end) {
        PREFETCH_R(a + 32);
        PREFETCH_R(b + 32);
        int av = *a, bv = *b;
        if (av <= bv) { *out++ = av; ++a; }
        else          { *out++ = bv; ++b; }
    }

    if (a < a_end) {
        size_t r = (size_t)(a_end - a);
        if (r <= TAIL_COPY_THRESHOLD) tiny_copy_int(out, a, r);
        else memcpy(out, a, r * sizeof(int));
    } else if (b < b_end) {
        size_t r = (size_t)(b_end - b);
        if (r <= TAIL_COPY_THRESHOLD) tiny_copy_int(out, b, r);
        else memcpy(out, b, r * sizeof(int));
    }
}

/* Recursive mergesort with a shortcut for very small spans */
static void msort_core(int* RESTRICT src, int* RESTRICT dst, size_t L, size_t R)
{
    if (R <= L) { dst[L] = src[L]; return; }

    size_t len = R - L + 1;
    if (len <= SMALL_RANGE_THRESHOLD) {
        memcpy(dst + L, src + L, len * sizeof(int));
        small_range_sort(dst, L, R);
        return;
    }

    size_t M = L + ((R - L) >> 1);
    msort_core(dst, src, L, M);
    msort_core(dst, src, M + 1, R);
    serial_merge(src, dst, L, M, M + 1, R, L);
}

/* Determine odd/even tree depth parity */
static int parity_of_levels(size_t size)
{
    if (size <= 1) return 0;
    size_t n = size - 1;
    int levels = 0;
    while (n) { ++levels; n >>= 1; }
    return levels & 1;
}

/* Entry point for serial sort */
void serial_msort_recursive(int* RESTRICT source, int* RESTRICT dest, size_t left, size_t right)
{
    msort_core(source, dest, left, right);
}

void optimized_serial_sort(int* RESTRICT arr_in, int* RESTRICT arr_out, size_t size)
{
    if (size <= 1) return;
    int odd_depth = parity_of_levels(size);
    if (odd_depth) {
        serial_msort_recursive(arr_in, arr_out, 0, size - 1);
        memcpy(arr_in, arr_out, size * sizeof(int));
    } else {
        serial_msort_recursive(arr_in, arr_in, 0, size - 1);
    }
}

/* Binary search for partitioning */
size_t binary_search(int val, int* arr, size_t len)
{
    size_t lo = 0, hi = len;
    while (lo < hi) {
        size_t mid = lo + ((hi - lo) >> 1);
        if (arr[mid] < val) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/* Thread workers */
void* sort_worker_mode1(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    size_t size = args->end - args->start + 1;
    optimized_serial_sort(args->temp_a + args->start, args->temp_b + args->start, size);
    return NULL;
}
void* sort_worker_mode2(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    optimized_serial_sort(args->temp_a, args->temp_b, args->start);
    return NULL;
}
void* merge_worker(void* arg) {
    ThreadArgs* m = (ThreadArgs*)arg;
    serial_merge(m->source, m->dest, m->start, m->end, m->start2, m->end2, m->dest_offset);
    return NULL;
}

/* ========= The Main Function with Two Modes (kept) ========= */
void mergeSortParallel(const int* values, unsigned int N, int* sorted)
{
    if (N > 536870912u) {
        /* ========== MODE 2: memory-safe for very large N (kept) ========== */
        int* workspace;
        int* scratch_half;
        pthread_t threads[4];
        ThreadArgs sort_args[4];
        ThreadArgs merge_args[4];
        size_t qsz, half_scratch_sz, half_size, chunk_size;
        size_t split_points[3];
        int i, pivot;

        qsz = N / 4u;
        half_scratch_sz = qsz * 2u;
        scratch_half = (int*)malloc(half_scratch_sz * sizeof(int));
        if (!scratch_half) { exit(1); }

        workspace = (int*)values; /* use input as staging like your original */
        memcpy(sorted, values, (size_t)N * sizeof(int));

        /* Sort first half in two pieces with limited scratch */
        sort_args[0] = (ThreadArgs){NULL, NULL, sorted,         scratch_half,       qsz, 0, 0, 0, 0, 0};
        sort_args[1] = (ThreadArgs){NULL, NULL, sorted + qsz,   scratch_half + qsz, qsz, 0, 0, 0, 0, 0};
        pthread_create(&threads[0], NULL, sort_worker_mode2, &sort_args[0]);
        pthread_create(&threads[1], NULL, sort_worker_mode2, &sort_args[1]);
        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);

        /* Sort second half in two pieces with the same scratch (reused) */
        sort_args[2] = (ThreadArgs){NULL, NULL, sorted + 2u * qsz, scratch_half,       qsz, 0, 0, 0, 0, 0};
        sort_args[3] = (ThreadArgs){NULL, NULL, sorted + 3u * qsz, scratch_half + qsz, N - (3u * qsz), 0, 0, 0, 0, 0};
        pthread_create(&threads[2], NULL, sort_worker_mode2, &sort_args[2]);
        pthread_create(&threads[3], NULL, sort_worker_mode2, &sort_args[3]);
        pthread_join(threads[2], NULL);
        pthread_join(threads[3], NULL);

        free(scratch_half);

        /* Two middle merges into workspace (which aliases values) */
        half_size = N / 2u;
        merge_args[0] = (ThreadArgs){sorted, workspace, NULL, NULL, 0, qsz - 1u, qsz, half_size - 1u, 0, 0};
        merge_args[1] = (ThreadArgs){sorted, workspace, NULL, NULL, half_size, half_size + qsz - 1u, half_size + qsz, N - 1u, half_size, 0};
        pthread_create(&threads[0], NULL, merge_worker, &merge_args[0]);
        pthread_create(&threads[1], NULL, merge_worker, &merge_args[1]);
        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);

        /* Final 4-way parallel merge using binary partitions */
        chunk_size = half_size / 4u;
        for (i = 0; i < 3; ++i) {
            pivot = workspace[(size_t)(i + 1) * chunk_size];
            split_points[i] = binary_search(pivot, workspace + half_size, N - half_size);
        }

        merge_args[0] = (ThreadArgs){workspace, sorted, NULL, NULL, 0, chunk_size - 1u, half_size, half_size + split_points[0] - 1u, 0, 0};
        merge_args[1] = (ThreadArgs){workspace, sorted, NULL, NULL, chunk_size, 2u * chunk_size - 1u, half_size + split_points[0], half_size + split_points[1] - 1u, chunk_size + split_points[0], 0};
        merge_args[2] = (ThreadArgs){workspace, sorted, NULL, NULL, 2u * chunk_size, 3u * chunk_size - 1u, half_size + split_points[1], half_size + split_points[2] - 1u, 2u * chunk_size + split_points[1], 0};
        merge_args[3] = (ThreadArgs){workspace, sorted, NULL, NULL, 3u * chunk_size, half_size - 1u, half_size + split_points[2], N - 1u, 3u * chunk_size + split_points[2], 0};

        for (i = 0; i < 4; ++i) pthread_create(&threads[i], NULL, merge_worker, &merge_args[i]);
        for (i = 0; i < 4; ++i) pthread_join(threads[i], NULL);

    } else {
        /* ========== MODE 1: fast path for moderate N (kept) ========== */
        int* temp_a = (int*)malloc((size_t)N * sizeof(int));
        int* temp_b = (int*)malloc((size_t)N * sizeof(int));
        pthread_t threads[4];
        ThreadArgs args[4];
        size_t quarter_size, half_size, chunk_size;
        size_t split_points[3];
        int i, pivot;

        if (!temp_a || !temp_b) { free(temp_a); free(temp_b); exit(1); }

        memcpy(temp_a, values, (size_t)N * sizeof(int));

        /* 4 parallel sorts over quarters of temp_a into temp_b (local buffers per quarter) */
        quarter_size = N / 4u;
        for (i = 0; i < 4; ++i) {
            args[i].temp_a = temp_a;
            args[i].temp_b = temp_b;
            args[i].start  = (size_t)i * quarter_size;
            args[i].end    = (i == 3) ? ((size_t)N - 1u) : ((size_t)(i + 1) * quarter_size - 1u);
            pthread_create(&threads[i], NULL, sort_worker_mode1, &args[i]);
        }
        for (i = 0; i < 4; ++i) pthread_join(threads[i], NULL);

        /* Two middle merges into temp_b -> temp_a (keep your order) */
        half_size = N / 2u;
        args[0] = (ThreadArgs){temp_a, temp_b, NULL, NULL, 0, quarter_size - 1u, quarter_size, half_size - 1u, 0, 0};
        args[1] = (ThreadArgs){temp_a, temp_b, NULL, NULL, half_size, half_size + quarter_size - 1u, half_size + quarter_size, N - 1u, half_size, 0};
        pthread_create(&threads[0], NULL, merge_worker, &args[0]);
        pthread_create(&threads[1], NULL, merge_worker, &args[1]);
        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);

        /* Final 4-way parallel merge using binary partitions, writing to 'sorted' */
        chunk_size = half_size / 4u;
        for (i = 0; i < 3; ++i) {
            pivot = temp_b[(size_t)(i + 1) * chunk_size];
            split_points[i] = binary_search(pivot, temp_b + half_size, N - half_size);
        }

        args[0] = (ThreadArgs){temp_b, sorted, NULL, NULL, 0, chunk_size - 1u, half_size, half_size + split_points[0] - 1u, 0, 0};
        args[1] = (ThreadArgs){temp_b, sorted, NULL, NULL, chunk_size, 2u * chunk_size - 1u, half_size + split_points[0], half_size + split_points[1] - 1u, chunk_size + split_points[0], 0};
        args[2] = (ThreadArgs){temp_b, sorted, NULL, NULL, 2u * chunk_size, 3u * chunk_size - 1u, half_size + split_points[1], half_size + split_points[2] - 1u, 2u * chunk_size + split_points[1], 0};
        args[3] = (ThreadArgs){temp_b, sorted, NULL, NULL, 3u * chunk_size, half_size - 1u, half_size + split_points[2], N - 1u, 3u * chunk_size + split_points[2], 0};

        for (i = 0; i < 4; ++i) pthread_create(&threads[i], NULL, merge_worker, &args[i]);
        for (i = 0; i < 4; ++i) pthread_join(threads[i], NULL);

        free(temp_a);
        free(temp_b);
    }
}
