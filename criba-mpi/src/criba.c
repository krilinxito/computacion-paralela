/*
 * criba.c — Criba de Eratostenes distribuida con MPI.
 * Compila:  mpicc -O2 -Wall -o criba criba.c -lm
 * Ejecuta:  mpirun -np 4 ./criba
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Limite superior de la busqueda */
#define N 10000000

int main(int argc, char *argv[])
{
    int rank, size;

    /* PASO 1: iniciar MPI (rank = id del proceso, size = total de procesos) */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double t_inicio = MPI_Wtime();

    /* PASO 2: rank 0 calcula los primos base hasta sqrt(N) (los unicos que
     * pueden marcar compuestos en cualquier segmento) */
    int raiz_n = (int)sqrt((double)N);
    int base_count = 0;
    int *base_primes = NULL;

    if (rank == 0) {
        char *marcado = calloc(raiz_n + 1, sizeof(char)); /* 0=primo, 1=compuesto */
        if (!marcado) { MPI_Abort(MPI_COMM_WORLD, 1); }

        marcado[0] = 1; /* 0 y 1 no son primos */
        marcado[1] = 1;

        for (int i = 2; (long long)i * i <= raiz_n; i++) {
            if (marcado[i] == 0) {
                for (int j = i * i; j <= raiz_n; j += i)
                    marcado[j] = 1;
            }
        }

        for (int i = 2; i <= raiz_n; i++)
            if (marcado[i] == 0) base_count++;

        base_primes = malloc(base_count * sizeof(int));
        if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }

        int idx = 0;
        for (int i = 2; i <= raiz_n; i++)
            if (marcado[i] == 0) base_primes[idx++] = i;

        free(marcado);
    }

    /* PASO 3: rank 0 envia los primos base a cada proceso (punto a punto).
     * tag 0 = cuantos primos hay; tag 1 = el arreglo de primos */
    if (rank == 0) {
        for (int dest = 1; dest < size; dest++) {
            MPI_Send(&base_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
            MPI_Send(base_primes, base_count, MPI_INT, dest, 1, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&base_count, 1, MPI_INT, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        base_primes = malloc(base_count * sizeof(int));
        if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }

        MPI_Recv(base_primes, base_count, MPI_INT, 0, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /* PASO 4: cada proceso toma su segmento [low, high] del rango [2, N] */
    long long low = 2 + (long long)rank * (N - 1) / size;
    long long high;
    if (rank == size - 1)
        high = N;
    else
        high = 2 + (long long)(rank + 1) * (N - 1) / size - 1;
    long long seg_size = high - low + 1;

    char *segmento = calloc(seg_size, sizeof(char)); /* 0=primo, 1=compuesto */
    if (!segmento) { MPI_Abort(MPI_COMM_WORLD, 1); }

    /* PASO 5: cribar el segmento local marcando los multiplos de cada primo base */
    for (int i = 0; i < base_count; i++) {
        long long p = base_primes[i];

        /* primer multiplo de p >= low */
        long long start = low;
        if (low % p != 0)
            start = low + (p - low % p);
        if (start == p) start += p; /* no marcar p mismo como compuesto */

        for (long long j = start; j <= high; j += p)
            segmento[j - low] = 1;
    }

    /* PASO 6: contar los primos del segmento (celdas con valor 0) */
    long long conteo_local = 0;
    for (long long i = 0; i < seg_size; i++)
        if (segmento[i] == 0) conteo_local++;

    /* PASO 7: rank 0 recolecta los conteos de cada proceso (punto a punto, tag 2) */
    long long *conteos = NULL;
    if (rank == 0) {
        conteos = malloc(size * sizeof(long long));
        conteos[0] = conteo_local;
        for (int src = 1; src < size; src++)
            MPI_Recv(&conteos[src], 1, MPI_LONG_LONG, src, 2,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 2, MPI_COMM_WORLD);
    }

    double t_fin = MPI_Wtime();
    double tiempo = t_fin - t_inicio;

    /* PASO 8: rank 0 muestra el resumen y escribe resultados.dat */
    if (rank == 0) {
        long long total = 0;
        for (int i = 0; i < size; i++) total += conteos[i];

        printf("\n=== Criba de Eratostenes Distribuida (N=%lld) ===\n", (long long)N);
        printf("Procesos MPI        : %d\n", size);
        printf("Primos encontrados  : %lld\n", total);
        printf("Tiempo de ejecucion : %.6f segundos\n\n", tiempo);

        printf("%-8s %-12s %-12s %-10s\n",
               "Proceso", "Inicio", "Fin", "Primos");
        printf("%-8s %-12s %-12s %-10s\n",
               "-------", "------", "---", "------");
        for (int i = 0; i < size; i++) {
            long long lo = 2 + (long long)i * (N - 1) / size;
            long long hi;
            if (i == size - 1)
                hi = N;
            else
                hi = 2 + (long long)(i + 1) * (N - 1) / size - 1;
            printf("%-8d %-12lld %-12lld %-10lld\n", i, lo, hi, conteos[i]);
        }

        printf("\nTIEMPO: %.6f\n", tiempo);

        FILE *fp = fopen("resultados/resultados.dat", "w");
        if (fp) {
            fprintf(fp, "# Criba de Eratostenes distribuida - N=%lld - Procesos=%d\n", (long long)N, size);
            fprintf(fp, "# Proceso  Inicio  Fin  Primos\n");
            for (int i = 0; i < size; i++) {
                long long lo = 2 + (long long)i * (N - 1) / size;
                long long hi;
                if (i == size - 1)
                    hi = N;
                else
                    hi = 2 + (long long)(i + 1) * (N - 1) / size - 1;
                fprintf(fp, "%d %lld %lld %lld\n", i, lo, hi, conteos[i]);
            }
            fclose(fp);
            printf("Resultados guardados en resultados.dat\n");
        }

        free(conteos);
    }

    /* PASO 9: rank 0 recolecta la lista de primos y la escribe en primos.dat.
     * Cada proceso envia sus primos (tag 3 = cuantos, tag 4 = el arreglo) y
     * rank 0 los escribe en orden de rank, asi el archivo queda ordenado.
     * OJO: para N grande este archivo es enorme (N=10^9 -> ~50 millones de lineas). */
    if (rank == 0) {
        FILE *fp_p = fopen("resultados/primos.dat", "w");
        if (fp_p) {
            fprintf(fp_p, "# Lista de primos hasta N=%lld\n", (long long)N);

            for (long long i = 0; i < seg_size; i++)
                if (segmento[i] == 0) fprintf(fp_p, "%lld\n", low + i);

            for (int src = 1; src < size; src++) {
                long long n_src;
                MPI_Recv(&n_src, 1, MPI_LONG_LONG, src, 3,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                long long *buf = malloc(n_src * sizeof(long long));
                if (!buf) { MPI_Abort(MPI_COMM_WORLD, 1); }

                MPI_Recv(buf, n_src, MPI_LONG_LONG, src, 4,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (long long i = 0; i < n_src; i++)
                    fprintf(fp_p, "%lld\n", buf[i]);

                free(buf);
            }

            fclose(fp_p);
            printf("Lista de primos guardada en resultados/primos.dat\n");
        }
    } else {
        long long *buf = malloc(conteo_local * sizeof(long long));
        if (!buf) { MPI_Abort(MPI_COMM_WORLD, 1); }

        long long k = 0;
        for (long long i = 0; i < seg_size; i++)
            if (segmento[i] == 0) buf[k++] = low + i;

        MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 3, MPI_COMM_WORLD);
        MPI_Send(buf, conteo_local, MPI_LONG_LONG, 0, 4, MPI_COMM_WORLD);

        free(buf);
    }

    free(base_primes);
    free(segmento);
    MPI_Finalize();
    return 0;
}
