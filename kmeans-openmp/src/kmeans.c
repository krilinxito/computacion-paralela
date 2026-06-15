/*
 * kmeans.c — K-Means Clustering paralelizado con OpenMP
 * Dataset: sintetico gaussiano escalable (N hasta 1M, D=16, K=10)
 *
 * Compilar: gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm
 * Ejecutar: ./kmeans
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* ------------------------------------------------------------------ */
/*  Constantes del problema                                             */
/* ------------------------------------------------------------------ */
#define N_DIMS    16    /* dimensiones por punto                        */
#define K         10    /* numero de clusters                           */
#define MAX_ITER  100   /* iteraciones maximas antes de forzar parada   */
#define THRESHOLD 1e-4  /* convergencia: desplazamiento maximo          */

/* ================================================================== */
/*  FUNCION: generate_dataset                                           */
/*  Genera n puntos en K clusters gaussianos usando Box-Muller.        */
/*  Cluster c tiene su centro en (c*10, c*10, ...) con ruido N(0,1).  */
/*  Semilla fija (42) para reproducibilidad entre corridas.            */
/* ================================================================== */
static void generate_dataset(int n, double *data) {
    srand(42);
    for (int i = 0; i < n; i++) {
        int cluster = i % K;
        for (int d = 0; d < N_DIMS; d++) {
            double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
            double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
            /* Box-Muller: transforma dos uniformes en una gaussiana    */
            double g  = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
            data[i * N_DIMS + d] = cluster * 10.0 + g;
        }
    }
}

/* ================================================================== */
/*  FUNCION: run_kmeans                                                 */
/*  Ejecuta K-Means sobre 'data' (n puntos, N_DIMS dims) con          */
/*  'num_threads' hilos OpenMP. Retorna el tiempo en segundos.         */
/*                                                                      */
/*  data[i*N_DIMS + d]       : punto i, dimension d (solo lectura)    */
/*  assignments[i]           : cluster asignado al punto i            */
/*  centroids[c*N_DIMS + d]  : centroide c, dimension d               */
/*  init_c                   : copia inicial de centroides (restaurar) */
/* ================================================================== */
static double run_kmeans(int n, int num_threads,
                         const double *data, int *assignments,
                         double *centroids, const double *init_c) {
    /* Restaurar centroides iniciales para comparacion justa entre runs */
    memcpy(centroids, init_c, (size_t)K * N_DIMS * sizeof(double));

    omp_set_num_threads(num_threads);
    double t_start = omp_get_wtime();

    for (int iter = 0; iter < MAX_ITER; iter++) {

        /* ---------------------------------------------------------- */
        /*  FASE 1: ASIGNACION — Paralelizada con OpenMP               */
        /*                                                             */
        /*  Por que se paraleliza:                                     */
        /*  - Es la fase dominante: O(N*K*D) operaciones.             */
        /*    Con N=1M, K=10, D=16 => 160M ops por iteracion.         */
        /*  - Cada iteracion i es COMPLETAMENTE INDEPENDIENTE:         */
        /*    * Lee data[i*N_DIMS..] y centroids (compartidos, solo    */
        /*      lectura: sin conflicto entre hilos)                    */
        /*    * Escribe solo en assignments[i] (indice disjunto:       */
        /*      sin race condition)                                    */
        /*  - Es "embarrassingly parallel": cero comunicacion.        */
        /* ---------------------------------------------------------- */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            double min_dist = 1e300;
            int    best     = 0;
            for (int c = 0; c < K; c++) {
                double dist = 0.0;
                for (int d = 0; d < N_DIMS; d++) {
                    double diff = data[i * N_DIMS + d]
                                - centroids[c * N_DIMS + d];
                    dist += diff * diff;
                }
                if (dist < min_dist) { min_dist = dist; best = c; }
            }
            assignments[i] = best;
        }
        /* Barrera implicita: todos los hilos sincronizan aqui         */

        /* ---------------------------------------------------------- */
        /*  FASE 2: RECALCULO DE CENTROIDES — Serial                   */
        /*                                                             */
        /*  Por que NO se paraleliza:                                  */
        /*  - Es O(N*D), mucho mas barata que la asignacion O(N*K*D). */
        /*  - Requiere acumular sumas por cluster: dos hilos sobre el  */
        /*    mismo cluster causarian race condition en new_c y counts */
        /*    sin usar reduction, lo que complica el codigo sin        */
        /*    beneficio apreciable (contribucion < 10% del tiempo).    */
        /* ---------------------------------------------------------- */
        double new_c[K * N_DIMS];
        int    counts[K];
        memset(new_c,  0, sizeof(new_c));
        memset(counts, 0, sizeof(counts));

        for (int i = 0; i < n; i++) {
            int c = assignments[i];
            counts[c]++;
            for (int d = 0; d < N_DIMS; d++)
                new_c[c * N_DIMS + d] += data[i * N_DIMS + d];
        }

        double max_shift = 0.0;
        for (int c = 0; c < K; c++) {
            if (counts[c] == 0) continue; /* cluster vacio: no dividir */
            for (int d = 0; d < N_DIMS; d++) {
                double new_val = new_c[c * N_DIMS + d] / counts[c];
                double shift   = fabs(new_val - centroids[c * N_DIMS + d]);
                if (shift > max_shift) max_shift = shift;
                centroids[c * N_DIMS + d] = new_val;
            }
        }

        if (max_shift < THRESHOLD) break;
    }

    return omp_get_wtime() - t_start;
}

/* ================================================================== */
/*  MAIN                                                                */
/* ================================================================== */
int main(void) {
    int sizes[]   = {10000, 100000, 500000, 1000000};
    int threads[] = {1, 2, 4};
    int n_sizes   = 4;
    int n_threads = 3;

    FILE *fout = fopen("resultados/benchmark.dat", "w");
    if (!fout) {
        fprintf(stderr, "Error: no se pudo crear resultados/benchmark.dat\n");
        return 1;
    }
    fprintf(fout, "# N  hilos  tiempo_s  speedup\n");

    printf("K-Means Paralelo con OpenMP\n");
    printf("Dataset sintetico gaussiano: D=%d dimensiones, K=%d clusters\n\n",
           N_DIMS, K);
    printf("%-10s %-8s %-14s %-8s\n", "N", "Hilos", "Tiempo (s)", "Speedup");
    printf("%-10s %-8s %-14s %-8s\n",
           "----------", "--------", "--------------", "--------");

    for (int si = 0; si < n_sizes; si++) {
        int n = sizes[si];
        printf("  [generando N=%d puntos...]\n", n);

        double *data        = malloc((size_t)n * N_DIMS * sizeof(double));
        int    *assignments = malloc((size_t)n * sizeof(int));
        double *centroids   = malloc((size_t)K * N_DIMS * sizeof(double));
        double *init_c      = malloc((size_t)K * N_DIMS * sizeof(double));

        if (!data || !assignments || !centroids || !init_c) {
            fprintf(stderr, "Error: sin memoria suficiente para N=%d\n", n);
            return 1;
        }

        /* Generar datos y guardar los primeros K puntos como centroides */
        generate_dataset(n, data);
        memcpy(init_c, data, (size_t)K * N_DIMS * sizeof(double));

        double t1 = 0.0;
        for (int ti = 0; ti < n_threads; ti++) {
            double t = run_kmeans(n, threads[ti], data, assignments,
                                  centroids, init_c);
            if (ti == 0) t1 = t;
            double speedup = (t > 0.0) ? t1 / t : 1.0;

            printf("%-10d %-8d %-14.6f %-8.3f\n",
                   n, threads[ti], t, speedup);
            fprintf(fout, "%d %d %.6f %.4f\n",
                    n, threads[ti], t, speedup);
        }
        printf("\n");

        free(data);
        free(assignments);
        free(centroids);
        free(init_c);
    }

    fclose(fout);
    printf("Benchmark exportado a: resultados/benchmark.dat\n");
    printf("Ejecutar 'gnuplot graficar.gp' para generar las graficas.\n");
    return 0;
}
