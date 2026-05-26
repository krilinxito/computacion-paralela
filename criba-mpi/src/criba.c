/*
 * criba.c — Criba de Eratóstenes Distribuida con MPI
 *
 * Compila:  mpicc -O2 -Wall -o criba criba.c -lm
 * Ejecuta:  mpirun -np 4 ./criba
 *
 * Cada proceso trabaja sobre un segmento del rango [2, N].
 * Los primos base (hasta sqrt(N)) se calculan en rank 0 y se
 * distribuyen a todos via MPI_Bcast. Luego cada proceso tamiza
 * su segmento de forma independiente. Al final rank 0 recolecta
 * los conteos y escribe resultados.dat.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Límite superior de la búsqueda */
#define N 10000000

int main(int argc, char *argv[])
{
    int rank, size;

    /* === PASO 1: Inicializar MPI ===
     * MPI_Init arranca el entorno de comunicación.
     * rank  = identificador de este proceso (0, 1, 2, ...)
     * size  = cantidad total de procesos lanzados con mpirun -np
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Iniciamos el cronómetro en todos los procesos al mismo tiempo */
    double t_inicio = MPI_Wtime();

    /* === PASO 2: Calcular primos base hasta sqrt(N) ===
     * Para tamizan cualquier segmento [low, high] necesitamos
     * todos los primos p tales que p <= sqrt(N).
     * Esos primos son los únicos que pueden "alcanzar" compuestos
     * dentro de nuestro segmento.
     *
     * Todos los procesos podrían calcularlo solos (son solo ~446
     * primos para N=10M), pero usamos MPI_Bcast intencionalmente
     * para demostrar el patrón de comunicación colectiva.
     */
    int raiz_n = (int)sqrt((double)N);
    int base_count = 0;
    int *base_primes = NULL;

    if (rank == 0) {
        /* Criba secuencial simple sobre [2, raiz_n] */
        char *marcado = calloc(raiz_n + 1, sizeof(char)); /* 0=primo, 1=compuesto */
        if (!marcado) { MPI_Abort(MPI_COMM_WORLD, 1); }

        marcado[0] = marcado[1] = 1; /* 0 y 1 no son primos */

        for (int i = 2; (long long)i * i <= raiz_n; i++) {
            if (!marcado[i]) {
                for (int j = i * i; j <= raiz_n; j += i)
                    marcado[j] = 1;
            }
        }

        /* Contar y copiar primos al arreglo base_primes */
        for (int i = 2; i <= raiz_n; i++)
            if (!marcado[i]) base_count++;

        base_primes = malloc(base_count * sizeof(int));
        if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }

        int idx = 0;
        for (int i = 2; i <= raiz_n; i++)
            if (!marcado[i]) base_primes[idx++] = i;

        free(marcado);
    }

    /* === PASO 3: Difundir los primos base a todos los procesos ===
     * MPI_Bcast envía un mensaje desde rank 0 a TODOS los demás.
     * Primero enviamos cuántos primos hay (base_count) para que
     * los otros procesos sepan cuánta memoria reservar.
     * Luego enviamos el arreglo completo.
     *
     * Sin esto, los procesos 1, 2, ... no sabrían qué números usar
     * para tamizan su segmento.
     */
    MPI_Bcast(&base_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        base_primes = malloc(base_count * sizeof(int));
        if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }
    }

    MPI_Bcast(base_primes, base_count, MPI_INT, 0, MPI_COMM_WORLD);

    /* === PASO 4: Determinar el segmento de este proceso ===
     * Dividimos el rango [2, N] en 'size' partes iguales.
     * Proceso i maneja los números desde 'low' hasta 'high'.
     *
     * Fórmula de distribución equitativa:
     *   low  = 2 + rank * (N-1) / size
     *   high = low del siguiente proceso - 1   (último proceso llega hasta N)
     *
     * Distribuimos el RANGO y no replicamos todo porque:
     *   - Con N=10M cada proceso solo necesita ~10MB/size de memoria
     *   - El trabajo de marcar compuestos también se divide
     */
    long long low  = 2 + (long long)rank * (N - 1) / size;
    long long high = (rank == size - 1)
                     ? N
                     : 2 + (long long)(rank + 1) * (N - 1) / size - 1;
    long long seg_size = high - low + 1;

    /* arreglo local del segmento: 0=primo, 1=compuesto */
    char *segmento = calloc(seg_size, sizeof(char));
    if (!segmento) { MPI_Abort(MPI_COMM_WORLD, 1); }

    /* === PASO 5: Tamizan el segmento local ===
     * Para cada primo base p, marcamos todos sus múltiplos
     * que caigan dentro de nuestro segmento [low, high].
     *
     * El primer múltiplo de p que está en nuestro segmento es:
     *   start = ceil(low / p) * p
     * Si ese número ES p (es decir, p está dentro del segmento),
     * avanzamos un paso más para no marcar al primo mismo.
     */
    for (int i = 0; i < base_count; i++) {
        long long p = base_primes[i];

        /* primer múltiplo de p >= low */
        long long start = ((low + p - 1) / p) * p;
        if (start == p) start += p; /* no marcar p mismo como compuesto */

        for (long long j = start; j <= high; j += p)
            segmento[j - low] = 1;
    }

    /* === PASO 6: Contar primos en este segmento ===
     * Recorremos el arreglo local; las celdas con valor 0 son primos.
     */
    long long conteo_local = 0;
    for (long long i = 0; i < seg_size; i++)
        if (!segmento[i]) conteo_local++;

    /* === PASO 7: Recolectar conteos en rank 0 ===
     * MPI_Gather envía el valor de cada proceso hacia rank 0,
     * quien los recibe en el arreglo 'conteos[]'.
     * Con esto rank 0 sabe cuántos primos encontró cada proceso.
     */
    long long *conteos = NULL;
    if (rank == 0)
        conteos = malloc(size * sizeof(long long));

    MPI_Gather(&conteo_local, 1, MPI_LONG_LONG,
               conteos, 1, MPI_LONG_LONG,
               0, MPI_COMM_WORLD);

    /* Detenemos el cronómetro */
    double t_fin  = MPI_Wtime();
    double tiempo = t_fin - t_inicio;

    /* === PASO 8: Rank 0 muestra resultados y escribe archivos ===
     * Solo rank 0 tiene todos los conteos, así que solo él imprime
     * el resumen global y escribe los archivos de salida.
     */
    if (rank == 0) {
        long long total = 0;
        for (int i = 0; i < size; i++) total += conteos[i];

        printf("\n=== Criba de Eratostenes Distribuida (N=%d) ===\n", N);
        printf("Procesos MPI        : %d\n", size);
        printf("Primos encontrados  : %lld\n", total);
        printf("Tiempo de ejecucion : %.6f segundos\n\n", tiempo);

        /* Detalle por proceso */
        printf("%-8s %-12s %-12s %-10s\n",
               "Proceso", "Inicio", "Fin", "Primos");
        printf("%-8s %-12s %-12s %-10s\n",
               "-------", "------", "---", "------");
        for (int i = 0; i < size; i++) {
            long long lo = 2 + (long long)i * (N - 1) / size;
            long long hi = (i == size - 1)
                           ? N
                           : 2 + (long long)(i + 1) * (N - 1) / size - 1;
            printf("%-8d %-12lld %-12lld %-10lld\n", i, lo, hi, conteos[i]);
        }

        /* Línea machine-readable para el benchmark del Makefile */
        printf("\nTIEMPO: %.6f\n", tiempo);

        /* --- Escribir resultados.dat ---
         * Formato: proceso  inicio  fin  cantidad_primos
         * Usado por graficar.gp para distribucion.png
         */
        FILE *fp = fopen("resultados/resultados.dat", "w");
        if (fp) {
            fprintf(fp, "# Criba de Eratostenes distribuida - N=%d - Procesos=%d\n", N, size);
            fprintf(fp, "# Proceso  Inicio  Fin  Primos\n");
            for (int i = 0; i < size; i++) {
                long long lo = 2 + (long long)i * (N - 1) / size;
                long long hi = (i == size - 1)
                               ? N
                               : 2 + (long long)(i + 1) * (N - 1) / size - 1;
                fprintf(fp, "%d %lld %lld %lld\n", i, lo, hi, conteos[i]);
            }
            fclose(fp);
            printf("Resultados guardados en resultados.dat\n");
        }

        free(conteos);
    }

    free(base_primes);
    free(segmento);
    MPI_Finalize();
    return 0;
}
