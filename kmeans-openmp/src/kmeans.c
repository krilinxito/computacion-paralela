/*
 * kmeans.c — K-Means Clustering paralelizado con OpenMP
 * Dataset: Iris (150 puntos, 4 dimensiones, K=3 clusters)
 *
 * Compilar: gcc -O2 -fopenmp -o kmeans kmeans.c -lm
 * Ejecutar: ./kmeans
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/* ------------------------------------------------------------------ */
/*  Constantes del problema                                             */
/* ------------------------------------------------------------------ */
#define N_POINTS   150   /* filas del dataset Iris                     */
#define N_DIMS     4     /* sepal_length, sepal_width, petal_length,   */
                         /* petal_width                                */
#define K          3     /* numero de clusters (= 3 especies de Iris)  */
#define MAX_ITER   100   /* iteraciones maximas antes de forzar parada */
#define THRESHOLD  0.001 /* convergencia: desplazamiento maximo de     */
                         /* centroides menor a este valor              */
#define LABEL_LEN  64    /* longitud maxima del nombre de especie      */

/* ------------------------------------------------------------------ */
/*  Datos globales                                                      */
/*  Se declaran globales para evitar paso de muchos punteros y para    */
/*  que el compilador los coloque en memoria estatica (no en el stack) */
/* ------------------------------------------------------------------ */
double data[N_POINTS][N_DIMS];          /* matriz 150 x 4 de features */
int    assignments[N_POINTS];           /* cluster asignado a cada pto */
double centroids[K][N_DIMS];            /* posicion actual de centroids*/
double init_centroids[K][N_DIMS];       /* copia inicial para repetir  */
                                        /* el experimento de forma justa*/

/* ================================================================== */
/*  FUNCION: load_iris                                                  */
/*  Lee el archivo CSV y llena el arreglo 'data'.                      */
/*  Retorna el numero de filas leidas (debe ser N_POINTS).             */
/* ================================================================== */
int load_iris(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: no se pudo abrir '%s'\n", filename);
        exit(1);
    }

    int n = 0;
    char label[LABEL_LEN];

    /* fscanf con formato "4 doubles + string separados por comas"     */
    /* El label (nombre de especie) se lee pero se descarta            */
    while (n < N_POINTS &&
           fscanf(f, "%lf,%lf,%lf,%lf,%63s",
                  &data[n][0], &data[n][1],
                  &data[n][2], &data[n][3],
                  label) == 5) {
        n++;
    }

    fclose(f);

    if (n != N_POINTS) {
        fprintf(stderr, "Advertencia: se esperaban %d filas, se leyeron %d\n",
                N_POINTS, n);
    }
    return n;
}

/* ================================================================== */
/*  FUNCION: init_centroids_from_data                                  */
/*  Inicializa centroides tomando filas especificas del dataset.       */
/*  seed_indices[c] = indice de la fila que sera el centroide c.       */
/*  Usamos {0, 50, 100}: primera flor de cada especie, garantizando   */
/*  buena separacion inicial y resultados reproducibles.               */
/*  Tambien guarda una copia en init_centroids para los re-runs.       */
/* ================================================================== */
void init_centroids_from_data(int seed_indices[K]) {
    for (int c = 0; c < K; c++) {
        for (int d = 0; d < N_DIMS; d++) {
            centroids[c][d]      = data[seed_indices[c]][d];
            init_centroids[c][d] = data[seed_indices[c]][d];
        }
    }
}

/* ================================================================== */
/*  FUNCION: run_kmeans                                                 */
/*  Ejecuta el algoritmo K-Means completo con 'num_threads' hilos.    */
/*  Retorna el tiempo transcurrido en segundos.                        */
/* ================================================================== */
double run_kmeans(int num_threads) {

    /* Restaurar centroides iniciales para comparacion justa          */
    /* Todos los experimentos (1, 2, 4 hilos) parten del mismo punto  */
    for (int c = 0; c < K; c++)
        for (int d = 0; d < N_DIMS; d++)
            centroids[c][d] = init_centroids[c][d];

    /* Configurar el numero de hilos OpenMP para esta corrida          */
    omp_set_num_threads(num_threads);

    /* Iniciar medicion de tiempo con el reloj de OpenMP               */
    /* omp_get_wtime() retorna segundos en punto flotante              */
    double t_start = omp_get_wtime();

    /* -------------------------------------------------------------- */
    /*  BUCLE PRINCIPAL DE K-MEANS                                     */
    /* -------------------------------------------------------------- */
    for (int iter = 0; iter < MAX_ITER; iter++) {

        /* ---------------------------------------------------------- */
        /*  FASE 1: ASIGNACION — Paralelizada con OpenMP               */
        /*                                                             */
        /*  Por que se paraleliza esta fase:                           */
        /*  - Es la fase mas costosa: O(N * K * D) operaciones         */
        /*    Con N=150, K=3, D=4 => 1800 multiplicaciones/sumas       */
        /*    Para N grande (millones) esto domina el tiempo total     */
        /*  - Cada iteracion 'i' es COMPLETAMENTE INDEPENDIENTE:       */
        /*    - Solo lee 'data[i]' y 'centroids' (lectura compartida,  */
        /*      sin conflicto)                                         */
        /*    - Solo escribe en 'assignments[i]' (indice unico,        */
        /*      sin race condition)                                    */
        /*  - Esto se llama "embarrassingly parallel" (paralelismo     */
        /*    trivial): ningun hilo interfiere con otro                */
        /* ---------------------------------------------------------- */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N_POINTS; i++) {
            double min_dist = 1e18;
            int    best     = 0;

            /* Calcular distancia euclidiana al cuadrado a cada centroide */
            /* (no usamos sqrt porque solo necesitamos comparar distancias,*/
            /* y x^2 < y^2 iff x < y para x,y >= 0)                      */
            for (int c = 0; c < K; c++) {
                double dist = 0.0;
                for (int d = 0; d < N_DIMS; d++) {
                    double diff = data[i][d] - centroids[c][d];
                    dist += diff * diff;
                }
                if (dist < min_dist) {
                    min_dist = dist;
                    best     = c;
                }
            }

            /* Guardar el cluster ganador para este punto              */
            /* assignments[i] es escrito por un solo hilo -> sin races */
            assignments[i] = best;
        }
        /* FIN zona paralela — aqui todos los hilos se sincronizan     */

        /* ---------------------------------------------------------- */
        /*  FASE 2: RECALCULO DE CENTROIDES — Serial (diseno Opcion A) */
        /*                                                             */
        /*  Por que NO se paraleliza esta fase:                        */
        /*  - Es O(N * D): solo 150*4=600 sumas, mucho mas barata      */
        /*    que la fase de asignacion (1800 ops)                     */
        /*  - Requiere acumular sumas por cluster, lo que introduce    */
        /*    dependencias de datos (race conditions si fuera paralela  */
        /*    sin reduction). Mantenerla serial es mas simple y seguro */
        /*  - En la practica no vale la pena paralelizarla:            */
        /*    su contribucion al tiempo total es < 25%                 */
        /* ---------------------------------------------------------- */
        double new_c[K][N_DIMS];
        int    counts[K];

        /* Inicializar acumuladores a cero                             */
        for (int c = 0; c < K; c++) {
            counts[c] = 0;
            for (int d = 0; d < N_DIMS; d++)
                new_c[c][d] = 0.0;
        }

        /* Acumular suma de coordenadas por cluster                    */
        for (int i = 0; i < N_POINTS; i++) {
            int c = assignments[i];
            counts[c]++;
            for (int d = 0; d < N_DIMS; d++)
                new_c[c][d] += data[i][d];
        }

        /* Calcular nueva posicion = promedio, y medir desplazamiento  */
        double max_shift = 0.0;
        for (int c = 0; c < K; c++) {
            if (counts[c] == 0) continue; /* guarda contra cluster vacio */
            for (int d = 0; d < N_DIMS; d++) {
                double new_val = new_c[c][d] / counts[c];
                double shift   = fabs(new_val - centroids[c][d]);
                if (shift > max_shift) max_shift = shift;
                centroids[c][d] = new_val;
            }
        }

        /* ---------------------------------------------------------- */
        /*  CRITERIO DE CONVERGENCIA                                   */
        /*  Si ningun centroide se movio mas de THRESHOLD, el          */
        /*  algoritmo convergio y podemos terminar antes de MAX_ITER   */
        /* ---------------------------------------------------------- */
        if (max_shift < THRESHOLD) break;

    } /* fin bucle K-Means */

    double t_end = omp_get_wtime();
    return t_end - t_start;
}

/* ================================================================== */
/*  FUNCION: export_results                                             */
/*  Escribe resultados.dat: una linea por punto con formato            */
/*  "petal_length petal_width cluster_id"                              */
/*  Se usan petal_length (dim 2) y petal_width (dim 3) como ejes      */
/*  porque visualmente separan mejor las 3 especies                   */
/* ================================================================== */
void export_results(const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Error: no se pudo crear '%s'\n", filename);
        return;
    }
    for (int i = 0; i < N_POINTS; i++) {
        /* cluster_id se imprime como 1,2,3 (base 1) para gnuplot      */
        fprintf(f, "%.2f %.2f %d\n",
                data[i][2],            /* petal_length                 */
                data[i][3],            /* petal_width                  */
                assignments[i] + 1);   /* cluster 1-indexed            */
    }
    fclose(f);
}

/* ================================================================== */
/*  FUNCION: export_speedup                                             */
/*  Escribe speedup.dat: una linea por experimento con formato         */
/*  "num_threads speedup"                                              */
/* ================================================================== */
void export_speedup(const char *filename,
                    int threads[], double times[], int n) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Error: no se pudo crear '%s'\n", filename);
        return;
    }
    double t1 = times[0]; /* tiempo base con 1 hilo                   */
    for (int i = 0; i < n; i++) {
        double speedup = (times[i] > 0.0) ? (t1 / times[i]) : 1.0;
        fprintf(f, "%d %.4f\n", threads[i], speedup);
    }
    fclose(f);
}

/* ================================================================== */
/*  MAIN                                                                */
/* ================================================================== */
int main(void) {

    /* -------------------------------------------------------------- */
    /*  PASO 1: Cargar dataset                                         */
    /* -------------------------------------------------------------- */
    printf("Cargando iris.csv...\n");
    int n = load_iris("datos/iris.csv");
    printf("  %d puntos cargados, %d dimensiones\n\n", n, N_DIMS);

    /* -------------------------------------------------------------- */
    /*  PASO 2: Inicializar centroides                                  */
    /*  Semillas: filas 0, 50, 100 (una por especie)                   */
    /* -------------------------------------------------------------- */
    int seeds[K] = {0, 50, 100};
    init_centroids_from_data(seeds);

    /* -------------------------------------------------------------- */
    /*  PASO 3: Ejecutar K-Means con 1, 2 y 4 hilos                   */
    /* -------------------------------------------------------------- */
    int    thread_counts[3] = {1, 2, 4};
    double times[3];

    printf("Ejecutando K-Means (K=%d, max_iter=%d)...\n\n", K, MAX_ITER);

    for (int i = 0; i < 3; i++) {
        times[i] = run_kmeans(thread_counts[i]);
    }

    /* -------------------------------------------------------------- */
    /*  PASO 4: Mostrar tabla de resultados                            */
    /* -------------------------------------------------------------- */
    printf("+---------+-------------+---------+\n");
    printf("| Hilos   | Tiempo (s)  | Speedup |\n");
    printf("+---------+-------------+---------+\n");
    for (int i = 0; i < 3; i++) {
        double speedup = times[0] / times[i];
        printf("| %7d | %11.6f | %7.3f |\n",
               thread_counts[i], times[i], speedup);
    }
    printf("+---------+-------------+---------+\n");
    printf("\nNota: con N=150 el overhead de hilos puede superar la ganancia.\n");
    printf("Para N > 100000 el speedup se acerca al numero de hilos.\n\n");

    /* -------------------------------------------------------------- */
    /*  PASO 5: Exportar resultados para gnuplot                       */
    /*  La ultima corrida (4 hilos) deja el estado final en            */
    /*  assignments[], que es lo que exportamos                        */
    /* -------------------------------------------------------------- */
    export_results("resultados/resultados.dat");
    printf("Resultados exportados a: resultados/resultados.dat\n");

    export_speedup("resultados/speedup.dat", thread_counts, times, 3);
    printf("Speedup exportado a:     resultados/speedup.dat\n");

    printf("\nEjecutar 'gnuplot graficar.gp' para generar las graficas.\n");

    return 0;
}
