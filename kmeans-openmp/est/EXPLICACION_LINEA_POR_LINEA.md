# Explicación línea por línea — K-Means Paralelo con OpenMP

Esta es una **referencia exhaustiva, línea a línea**, de todo el código del proyecto. A diferencia de `GUIA_ESTUDIO.md` (que es un tutorial conceptual "desde cero"), aquí recorremos **cada línea** de los tres archivos fuente sin saltarnos ninguna, explicando:

- **Qué hace** la línea.
- **Por qué está ahí** (la intención).
- **Sus implicaciones** (rendimiento, memoria, paralelismo, riesgos).
- **En qué contexto se ejecuta**: si está dentro de un loop, de qué loop, o dentro de una región paralela.

---

## Tabla de Contenidos

- [Leyenda de etiquetas de contexto](#leyenda-de-etiquetas-de-contexto)
- [Mapa del proyecto](#mapa-del-proyecto)
- [1. `src/kmeans.c`](#1-srckmeansc)
  - [1.1 Cabecera del archivo (líneas 1-6)](#11-cabecera-del-archivo-líneas-1-6)
  - [1.2 Includes (líneas 8-12)](#12-includes-líneas-8-12)
  - [1.3 Constantes con #define (líneas 14-17)](#13-constantes-con-define-líneas-14-17)
  - [1.4 Función generate_dataset (líneas 19-32)](#14-función-generate_dataset-líneas-19-32)
  - [1.5 Función run_kmeans (líneas 34-96)](#15-función-run_kmeans-líneas-34-96)
  - [1.6 Función main (líneas 98-160)](#16-función-main-líneas-98-160)
- [2. `Makefile`](#2-makefile)
- [3. `graficar.gp`](#3-graficargp)
- [Resumen de paralelismo](#resumen-de-paralelismo)

---

## Leyenda de etiquetas de contexto

A lo largo del documento verás estas etiquetas antes de explicar una línea, para que siempre sepas **dónde** se está ejecutando:

| Etiqueta | Significado |
|----------|-------------|
| `[GLOBAL]` | Fuera de toda función (preprocesador, declaraciones). |
| `[SERIAL]` | Se ejecuta una sola vez, en un único hilo. |
| `[LOOP iteraciones]` | Dentro del bucle externo `for (iter...)` de K-Means. |
| `[PARALELO]` | Dentro de una región `#pragma omp parallel for`: **se ejecuta en varios hilos a la vez**. |
| `[LOOP i]` / `[LOOP c]` / `[LOOP d]` | Dentro de un bucle sobre puntos (`i`), clusters (`c`) o dimensiones (`d`). |
| `[LOOP si]` / `[LOOP ti]` | En `main`: bucle sobre tamaños (`si`) o sobre número de hilos (`ti`). |

Cuando una línea está dentro de varios niveles a la vez, se listan todos, p. ej. `[PARALELO · LOOP i · LOOP c · LOOP d]` = está en el loop más interno de dimensiones, que a su vez está dentro del loop de clusters, dentro del loop de puntos que se reparte entre hilos.

---

## Mapa del proyecto

| Archivo | Qué hace |
|---------|----------|
| `src/kmeans.c` | El programa: genera datos sintéticos, corre K-Means con 1/2/4 hilos, mide tiempos y exporta un benchmark. |
| `Makefile` | Receta de compilación y atajos (`make run`, `make plots`, `make clean`, ...). |
| `graficar.gp` | Script de Gnuplot que dibuja dos gráficas (speedup y escalabilidad) a partir del benchmark. |

El flujo completo es: `make` compila → `./kmeans` genera `resultados/benchmark.dat` → `gnuplot graficar.gp` produce `speedup.png` y `scaling.png`.

---

## 1. `src/kmeans.c`

### 1.1 Cabecera del archivo (líneas 1-6)

```c
/*
 * kmeans.c — K-Means Clustering paralelizado con OpenMP
 *
 * Compilar: gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm
 * Ejecutar: ./kmeans
 */
```

`[GLOBAL]`

- **Línea 1** `/*` — abre un comentario de bloque. En C, todo lo que esté entre `/*` y `*/` lo ignora el compilador; sirve solo para el lector humano.
- **Línea 2** `* kmeans.c — ...` — nombre y propósito del archivo. Una sola frase que dice qué es: una implementación de K-Means paralelizada con OpenMP.
- **Línea 3** `*` — línea en blanco decorativa dentro del comentario (el `*` al inicio es convención de estilo, no obligatorio).
- **Línea 4** `* Compilar: gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm` — documenta el comando exacto de compilación. Cada flag importa:
  - `-O2` activa optimizaciones del compilador (nivel 2): vectorización, inline, etc. **Implicación de rendimiento**: sin esto el código serial sería mucho más lento y el speedup medido estaría distorsionado.
  - `-fopenmp` activa OpenMP: hace que el compilador entienda los `#pragma omp` y enlace la librería de runtime de hilos. **Sin este flag, los `#pragma` se ignoran y el programa corre en un solo hilo.**
  - `-o kmeans` nombra el ejecutable de salida `kmeans`.
  - `-lm` enlaza la librería matemática (necesaria por `sqrt`, `log`, `cos`, `fabs`). **Debe ir al final** (ver el Makefile, línea 6).
- **Línea 5** `* Ejecutar: ./kmeans` — cómo lanzar el binario una vez compilado.
- **Línea 6** `*/` — cierra el comentario de bloque.

### 1.2 Includes (líneas 8-12)

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
```

`[GLOBAL]` — Las directivas `#include` las procesa el **preprocesador** antes de compilar: pegan el contenido del archivo de cabecera indicado. Aquí traen las declaraciones de las funciones de librería que usamos.

- **Línea 8** `#include <stdio.h>` — entrada/salida estándar. Aporta `printf`, `fprintf`, `fopen`, `fclose`, `FILE`, `stderr`. Se usa para imprimir la tabla de resultados y escribir el `.dat`.
- **Línea 9** `#include <stdlib.h>` — utilidades generales. Aporta `malloc`, `free` (memoria dinámica), `srand`, `rand` (números aleatorios) y `RAND_MAX`. Indispensable para generar el dataset y reservar los buffers grandes.
- **Línea 10** `#include <string.h>` — manejo de memoria/cadenas. Aporta `memcpy` (copiar bloques de bytes) y `memset` (rellenar bloques con un valor). Se usan para copiar centroides y poner acumuladores a cero.
- **Línea 11** `#include <math.h>` — funciones matemáticas. Aporta `sqrt`, `log`, `cos`, `fabs` y la macro `M_PI`. Necesario para Box-Muller (generación gaussiana) y para el cálculo del desplazamiento de centroides. **Por esto hace falta `-lm` al enlazar.**
- **Línea 12** `#include <omp.h>` — la API de OpenMP. Aporta `omp_set_num_threads` y `omp_get_wtime`. **Nota**: los `#pragma omp` no necesitan este header (los entiende el compilador con `-fopenmp`), pero las **funciones** de runtime sí.

### 1.3 Constantes con #define (líneas 14-17)

```c
#define N_DIMS    16    /* dimensiones por punto    */
#define K         10    /* numero de clusters       */
#define MAX_ITER  100   /* iteraciones maximas      */
#define THRESHOLD 1e-4  /* umbral de convergencia   */
```

`[GLOBAL]` — Cada `#define` define una **macro**: el preprocesador reemplaza textualmente cada aparición del nombre por su valor antes de compilar. No ocupan memoria en ejecución; son constantes de tiempo de compilación.

- **Línea 14** `#define N_DIMS 16` — cada punto tiene 16 coordenadas (dimensiones). **Implicación**: el coste de calcular una distancia es proporcional a `N_DIMS`, así que este valor multiplica directamente el trabajo del bucle más interno.
- **Línea 15** `#define K 10` — el algoritmo busca 10 clusters (grupos). También es el número de centroides. **Implicación**: el bucle de clusters en la fase de asignación se ejecuta `K` veces por cada punto.
- **Línea 16** `#define MAX_ITER 100` — tope de iteraciones de K-Means. Aunque no converja, el algoritmo se detiene a las 100 vueltas. Es una **red de seguridad** contra bucles infinitos.
- **Línea 17** `#define THRESHOLD 1e-4` — umbral de convergencia (`0.0001`). Si en una iteración ningún centroide se mueve más que esto, consideramos que el algoritmo "ya no cambia" y paramos.
- **Por qué macros y no variables**: porque se usan como **tamaños de arreglos en el stack**, p. ej. `double new_c[K * N_DIMS];` (línea 66). El compilador necesita conocer `K` y `N_DIMS` como constantes para reservar ese espacio. Si fueran variables normales, eso sería un VLA (arreglo de longitud variable), menos portable.

### 1.4 Función generate_dataset (líneas 19-32)

```c
/* Genera n puntos en K clusters gaussianos. Semilla fija para reproducibilidad. */
static void generate_dataset(int n, double *data) {
    srand(42);
    for (int i = 0; i < n; i++) {
        int cluster = i % K;
        for (int d = 0; d < N_DIMS; d++) {
            double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
            double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
            /* Box-Muller: dos uniformes -> una gaussiana */
            double g  = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
            data[i * N_DIMS + d] = cluster * 10.0 + g;
        }
    }
}
```

- **Línea 19** `/* Genera n puntos... */` `[GLOBAL]` — comentario que resume la función: crea `n` puntos repartidos en `K` clusters con ruido gaussiano, con semilla fija para que el resultado sea siempre el mismo.
- **Línea 20** `static void generate_dataset(int n, double *data) {` `[SERIAL]` — firma de la función.
  - `static` significa que la función es **privada a este archivo**: no se puede llamar desde otros `.c`. Buena práctica para funciones auxiliares; evita choques de nombres.
  - `void` = no devuelve nada (escribe el resultado en el buffer que recibe).
  - `int n` = cuántos puntos generar.
  - `double *data` = puntero al bloque de memoria donde se escriben los puntos. El llamador (`main`) ya reservó este espacio; aquí solo lo rellenamos.
- **Línea 21** `srand(42);` `[SERIAL]` — inicializa el generador pseudoaleatorio con la **semilla fija 42**. **Implicación clave**: con la misma semilla, `rand()` produce siempre la misma secuencia, así que **el dataset es idéntico en cada ejecución**. Esto es esencial para un benchmark justo: comparamos tiempos sobre exactamente los mismos datos.
- **Línea 22** `for (int i = 0; i < n; i++) {` `[LOOP i]` — bucle externo sobre los `n` puntos. `i` es el índice del punto actual. Se ejecuta `n` veces (hasta 1.000.000).
- **Línea 23** `int cluster = i % K;` `[LOOP i]` — asigna a cada punto un cluster "verdadero" según `i módulo K`. Así los puntos se reparten cíclicamente: punto 0→cluster 0, 1→1, ..., 9→9, 10→0, etc. **Nota didáctica**: K-Means *no usará* esta etiqueta; sirve solo para colocar cada punto cerca de un centro distinto y que existan grupos reales que descubrir.
- **Línea 24** `for (int d = 0; d < N_DIMS; d++) {` `[LOOP i · LOOP d]` — bucle interno sobre las 16 dimensiones del punto. Genera una coordenada por vuelta.
- **Línea 25** `double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);` `[LOOP i · LOOP d]` — genera un número uniforme en el intervalo **(0, 1]**. `rand()` da un entero entre `0` y `RAND_MAX`. Al sumar `1.0` arriba y abajo, el rango va de `1/(RAND_MAX+1)` hasta `1`, evitando el valor exacto `0`. **Por qué importa**: en la línea 28 se hace `log(u1)`, y `log(0)` es `-infinito` → reventaría el cálculo. Este `+1.0` es una protección.
- **Línea 26** `double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);` `[LOOP i · LOOP d]` — segundo número uniforme independiente, también en (0, 1]. Box-Muller necesita **dos** uniformes para producir un valor gaussiano.
- **Línea 27** `/* Box-Muller: dos uniformes -> una gaussiana */` `[LOOP i · LOOP d]` — comentario que nombra el método. Box-Muller es una transformación matemática que convierte dos números uniformes en uno con **distribución normal (campana de Gauss)**.
- **Línea 28** `double g = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);` `[LOOP i · LOOP d]` — la fórmula de Box-Muller. Desglose:
  - `log(u1)` es negativo (porque `u1 ≤ 1`), así que `-2.0 * log(u1)` es positivo → `sqrt` de un positivo, todo válido.
  - `cos(2.0 * M_PI * u2)` toma un ángulo aleatorio en una vuelta completa.
  - El producto da `g`, un valor con **media 0 y desviación 1** (gaussiano estándar). Esto es el "ruido" que dispersa el punto alrededor del centro de su cluster.
- **Línea 29** `data[i * N_DIMS + d] = cluster * 10.0 + g;` `[LOOP i · LOOP d]` — escribe la coordenada final.
  - `cluster * 10.0` es el **centro** de ese cluster en esta dimensión: cluster 0 centrado en 0, cluster 1 en 10, ..., cluster 9 en 90. La separación de 10 hace que los grupos estén bien diferenciados.
  - `+ g` añade el ruido gaussiano, dispersando el punto alrededor de ese centro.
  - `data[i * N_DIMS + d]` es el **layout row-major (plano)**: aunque conceptualmente `data` es una matriz de `n × N_DIMS`, se almacena como un único arreglo 1D. La fórmula `i * N_DIMS + d` calcula la posición del elemento `(i, d)`. **Por qué un arreglo plano**: es más eficiente en memoria (un solo `malloc`, sin punteros a punteros) y mejora la **localidad de caché** (las dimensiones de un punto están contiguas), lo que acelera mucho los recorridos.
- **Línea 30** `}` `[LOOP i]` — cierra el bucle de dimensiones (`LOOP d`).
- **Línea 31** `}` `[SERIAL]` — cierra el bucle de puntos (`LOOP i`).
- **Línea 32** `}` `[GLOBAL]` — cierra la función `generate_dataset`.

### 1.5 Función run_kmeans (líneas 34-96)

Esta es la función central: ejecuta el algoritmo K-Means completo con un número dado de hilos y devuelve cuánto tardó.

```c
/* Ejecuta K-Means con num_threads hilos. Retorna el tiempo en segundos. */
static double run_kmeans(int n, int num_threads,
                         const double *data, int *assignments,
                         double *centroids, const double *init_c) {
```

- **Línea 34** `/* Ejecuta K-Means... */` `[GLOBAL]` — comentario: ejecuta K-Means con `num_threads` hilos y retorna el tiempo en segundos.
- **Líneas 35-37** `static double run_kmeans(...)` `[SERIAL]` — firma. Devuelve un `double` (el tiempo medido). Parámetros:
  - `int n` — número de puntos.
  - `int num_threads` — cuántos hilos usar en esta corrida.
  - `const double *data` — los puntos. **`const`** promete que la función **no modifica** los datos: es de solo lectura. Esto es importante porque `data` se comparte entre todos los hilos en la región paralela, y compartir datos de solo lectura es seguro (no hay race conditions).
  - `int *assignments` — buffer de salida: a qué cluster pertenece cada punto. **No** es `const` porque sí se escribe.
  - `double *centroids` — los centroides actuales (se leen y se actualizan). No es `const`.
  - `const double *init_c` — los centroides **iniciales**, de solo lectura. Se usan para reiniciar `centroids` al principio de cada corrida.
  - **Por qué `init_c` separado de `centroids`**: K-Means modifica los centroides en cada iteración. Para comparar 1 vs 2 vs 4 hilos de forma justa, todas las corridas deben **empezar desde los mismos centroides**. `init_c` guarda esa copia original intacta; `centroids` es la copia de trabajo que se va modificando.

```c
    /* Restaurar centroides iniciales para una comparacion justa entre corridas */
    memcpy(centroids, init_c, (size_t)K * N_DIMS * sizeof(double));
```

- **Línea 38** `/* Restaurar centroides... */` `[SERIAL]` — comentario que explica el porqué de la línea siguiente.
- **Línea 39** `memcpy(centroids, init_c, (size_t)K * N_DIMS * sizeof(double));` `[SERIAL]` — copia los `K * N_DIMS` valores de `init_c` (originales) sobre `centroids` (trabajo). Así cada corrida arranca desde el mismo punto.
  - `memcpy(destino, origen, bytes)` copia bloques de bytes crudos: rapidísimo.
  - `(size_t)K * N_DIMS * sizeof(double)` calcula los bytes a copiar: `10 * 16 * 8 = 1280` bytes. El cast `(size_t)` asegura que la multiplicación se haga en aritmética de enteros sin signo del ancho del puntero, evitando desbordamientos (aquí el número es pequeño, pero es una buena costumbre, crítica en los `malloc` grandes de `main`).

```c
    omp_set_num_threads(num_threads);
    double t_start = omp_get_wtime();
```

- **Línea 41** `omp_set_num_threads(num_threads);` `[SERIAL]` — le dice al runtime de OpenMP cuántos hilos usar en las siguientes regiones paralelas. **Implicación**: con `num_threads = 1` la región paralela corre en un solo hilo (la corrida base para medir speedup); con `4`, en cuatro.
- **Línea 42** `double t_start = omp_get_wtime();` `[SERIAL]` — toma una marca de tiempo "de reloj de pared" (wall-clock) **antes** del algoritmo. `omp_get_wtime` devuelve segundos con alta resolución. **Por qué wall-clock y no tiempo de CPU**: con varios hilos, el tiempo de CPU sumaría el de todos los núcleos; lo que queremos medir es el **tiempo real transcurrido**, que es lo que el usuario percibe y lo que define el speedup.

```c
    for (int iter = 0; iter < MAX_ITER; iter++) {
```

- **Línea 44** `for (int iter = 0; iter < MAX_ITER; iter++) {` `[LOOP iteraciones]` — **bucle externo de K-Means**. Cada vuelta es una iteración completa del algoritmo: (1) asignar puntos a centroides, (2) recalcular centroides. Se repite hasta `MAX_ITER` (100) o hasta converger (línea 92). Todo lo que sigue, hasta la línea 93, está dentro de este bucle.

#### Fase 1: asignación (paralela) — líneas 46-63

```c
        /* Fase 1: asignacion (paralela) — cada punto al centroide mas cercano */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            double min_dist = 1e300;
            int    best     = 0;
            for (int c = 0; c < K; c++) {
                double dist = 0.0;
                for (int d = 0; d < N_DIMS; d++) {
                    double diff = data[i * N_DIMS + d] - centroids[c * N_DIMS + d];
                    dist += diff * diff;
                }
                if (dist < min_dist) {
                    min_dist = dist;
                    best = c;
                }
            }
            assignments[i] = best;
        }
```

- **Línea 46** `/* Fase 1: asignacion (paralela)... */` `[LOOP iteraciones]` — comentario: esta fase asigna cada punto al centroide más cercano, y es la parte que paralelizamos.
- **Línea 47** `#pragma omp parallel for schedule(static)` `[LOOP iteraciones]` — **la directiva clave de paralelización**. Le dice a OpenMP que reparta el siguiente bucle `for` entre los hilos. Palabra por palabra:
  - `#pragma omp` — directiva de OpenMP.
  - `parallel` — crea un **equipo de hilos** (tantos como `num_threads`).
  - `for` — reparte las iteraciones del bucle siguiente entre esos hilos (cada hilo procesa un subconjunto de los puntos).
  - `schedule(static)` — divide los `n` puntos en bloques contiguos del mismo tamaño, uno por hilo (p. ej. con 4 hilos y 1M puntos, cada uno hace 250.000 puntos consecutivos). **Por qué `static`**: el trabajo por punto es **uniforme** (siempre `K * N_DIMS` operaciones), así que un reparto fijo es óptimo y tiene **overhead casi nulo** (no hay coordinación dinámica entre hilos). Un `dynamic` añadiría coste de sincronización sin beneficio.
  - **Variables privadas vs compartidas** (regla de OpenMP): `i` es privada automáticamente (es la variable del `for` paralelo). Todo lo declarado **dentro** del bucle (`min_dist`, `best`, `c`, `dist`, `d`, `diff`) es **privado** de cada hilo: cada uno tiene su propia copia, sin interferencias. En cambio `data`, `centroids` y `assignments` son **compartidos** (declarados fuera). `data` y `centroids` son de solo lectura aquí → seguro compartirlos. `assignments` se escribe, pero **sin race condition** (ver línea 62).
- **Línea 48** `for (int i = 0; i < n; i++) {` `[PARALELO · LOOP i]` — bucle sobre los puntos, **repartido entre hilos**. Distintos valores de `i` los procesan distintos hilos simultáneamente.
- **Línea 49** `double min_dist = 1e300;` `[PARALELO · LOOP i]` — inicializa la menor distancia encontrada hasta ahora con un número gigantesco (`1e300`, casi el máximo de un `double`). Funciona como "infinito": cualquier distancia real será menor, así que el primer centroide siempre la actualizará. **Privada por hilo.**
- **Línea 50** `int best = 0;` `[PARALELO · LOOP i]` — guarda el índice del centroide más cercano hasta ahora; arranca en 0. **Privada por hilo.**
- **Línea 51** `for (int c = 0; c < K; c++) {` `[PARALELO · LOOP i · LOOP c]` — recorre los `K` centroides para encontrar el más cercano al punto `i`. Se ejecuta `K` (=10) veces por punto.
- **Línea 52** `double dist = 0.0;` `[PARALELO · LOOP i · LOOP c]` — acumulador de la distancia (al cuadrado) entre el punto `i` y el centroide `c`. Se reinicia a 0 para cada centroide. **Privada por hilo.**
- **Línea 53** `for (int d = 0; d < N_DIMS; d++) {` `[PARALELO · LOOP i · LOOP c · LOOP d]` — **bucle más interno**: recorre las 16 dimensiones para sumar las diferencias al cuadrado. Este es el bucle que más veces se ejecuta de todo el programa: `n × K × N_DIMS` veces (¡hasta 160 millones por iteración de K-Means!). Por eso es donde más importa la optimización y la localidad de caché.
- **Línea 54** `double diff = data[i * N_DIMS + d] - centroids[c * N_DIMS + d];` `[PARALELO · LOOP i · LOOP c · LOOP d]` — diferencia entre la coordenada `d` del punto `i` y la del centroide `c`. Usa el mismo indexado row-major plano explicado en la línea 29.
- **Línea 55** `dist += diff * diff;` `[PARALELO · LOOP i · LOOP c · LOOP d]` — eleva al cuadrado la diferencia y la acumula. Al terminar el bucle de dimensiones, `dist` es la **distancia euclídea al cuadrado**.
  - **Por qué al cuadrado y sin `sqrt`**: la raíz cuadrada es cara y **no cambia cuál es el menor**. Si `a² < b²` (con a,b ≥ 0) entonces `a < b`. Como solo comparamos para encontrar el mínimo, nos ahorramos `n × K` raíces cuadradas por iteración. **Optimización importante.**
- **Línea 56** `}` `[PARALELO · LOOP i · LOOP c]` — cierra el bucle de dimensiones.
- **Línea 57** `if (dist < min_dist) {` `[PARALELO · LOOP i · LOOP c]` — ¿este centroide está más cerca que el mejor hasta ahora?
- **Línea 58** `min_dist = dist;` `[PARALELO · LOOP i · LOOP c]` — sí: actualiza la menor distancia.
- **Línea 59** `best = c;` `[PARALELO · LOOP i · LOOP c]` — y recuerda qué centroide es el mejor.
- **Línea 60** `}` `[PARALELO · LOOP i · LOOP c]` — cierra el `if`.
- **Línea 61** `}` `[PARALELO · LOOP i]` — cierra el bucle de centroides. Al salir, `best` contiene el centroide más cercano al punto `i`.
- **Línea 62** `assignments[i] = best;` `[PARALELO · LOOP i]` — guarda el resultado: el punto `i` pertenece al cluster `best`.
  - **Por qué NO hay race condition aquí** (a pesar de que `assignments` es compartido y se escribe en paralelo): cada hilo procesa **valores distintos de `i`**, y cada `i` escribe en **una posición distinta** `assignments[i]`. Dos hilos nunca tocan la misma celda → no hay conflicto. Esta es la razón por la que esta fase es "vergonzosamente paralela" (embarrassingly parallel) y escala bien.
- **Línea 63** `}` `[LOOP iteraciones]` — cierra el bucle paralelo. **Aquí hay una barrera implícita**: OpenMP espera a que **todos** los hilos terminen la fase de asignación antes de continuar. Necesario, porque la fase 2 usa todos los `assignments`.

#### Fase 2: recálculo de centroides (serial) — líneas 65-93

```c
        /* Fase 2: recalculo de centroides (serial) — promedio de cada cluster */
        double new_c[K * N_DIMS];
        int    counts[K];
        memset(new_c,  0, sizeof(new_c));
        memset(counts, 0, sizeof(counts));
```

- **Línea 65** `/* Fase 2: recalculo... */` `[LOOP iteraciones]` — comentario: ahora recalculamos cada centroide como el promedio de los puntos que se le asignaron. Es serial.
- **Línea 66** `double new_c[K * N_DIMS];` `[LOOP iteraciones]` — arreglo temporal en el **stack** para acumular las sumas de coordenadas de cada cluster (`10 × 16 = 160` doubles). Que `K` y `N_DIMS` sean macros constantes permite reservarlo así. **No está inicializado todavía** (contiene basura) — por eso la línea 68.
- **Línea 67** `int counts[K];` `[LOOP iteraciones]` — arreglo en el stack para contar cuántos puntos cayeron en cada cluster (`10` enteros). Necesario para dividir la suma y obtener el promedio.
- **Línea 68** `memset(new_c, 0, sizeof(new_c));` `[LOOP iteraciones]` — pone a **cero** todos los bytes de `new_c`. `sizeof(new_c)` da el tamaño total en bytes del arreglo (`160 × 8 = 1280`). Imprescindible antes de acumular sumas.
- **Línea 69** `memset(counts, 0, sizeof(counts));` `[LOOP iteraciones]` — pone a cero los contadores. `sizeof(counts)` = `10 × 4 = 40` bytes.

```c
        for (int i = 0; i < n; i++) {
            int c = assignments[i];
            counts[c]++;
            for (int d = 0; d < N_DIMS; d++) {
                new_c[c * N_DIMS + d] += data[i * N_DIMS + d];
            }
        }
```

- **Línea 71** `for (int i = 0; i < n; i++) {` `[LOOP iteraciones · LOOP i]` — recorre **todos** los puntos para sumarlos en su cluster. **Serial** (no hay `#pragma`): se ejecuta en un solo hilo.
- **Línea 72** `int c = assignments[i];` `[LOOP iteraciones · LOOP i]` — lee a qué cluster pertenece el punto `i` (calculado en la fase 1).
- **Línea 73** `counts[c]++;` `[LOOP iteraciones · LOOP i]` — incrementa el contador de ese cluster.
  - **Por qué esta fase es serial / por qué tendría race si se paralelizara ingenuamente**: aquí distintos puntos `i` pueden pertenecer al **mismo** cluster `c`, así que `counts[c]` y `new_c[c*...]` los actualizarían **varios hilos a la vez** → race condition (lecturas/escrituras pisándose, resultados corruptos). Paralelizar esto bien exige reducciones o acumuladores por hilo, que añaden complejidad. Como además esta fase es **mucho más barata** que la fase 1 (un solo recorrido sumando, sin el bucle de `K` centroides), se deja serial: paralelizarla aportaría poco al speedup total (ver Ley de Amdahl en la guía).
- **Línea 74** `for (int d = 0; d < N_DIMS; d++) {` `[LOOP iteraciones · LOOP i · LOOP d]` — recorre las dimensiones del punto para sumarlas.
- **Línea 75** `new_c[c * N_DIMS + d] += data[i * N_DIMS + d];` `[LOOP iteraciones · LOOP i · LOOP d]` — acumula la coordenada `d` del punto `i` en la suma de su cluster `c`. Al terminar el bucle externo, `new_c[c*...]` contiene la **suma** de todas las coordenadas de los puntos del cluster `c`.
- **Línea 76** `}` `[LOOP iteraciones · LOOP i]` — cierra el bucle de dimensiones.
- **Línea 77** `}` `[LOOP iteraciones]` — cierra el bucle de acumulación.

```c
        double max_shift = 0.0;
        for (int c = 0; c < K; c++) {
            if (counts[c] == 0) continue; /* cluster vacio: no dividir */
            for (int d = 0; d < N_DIMS; d++) {
                double new_val = new_c[c * N_DIMS + d] / counts[c];
                double shift   = fabs(new_val - centroids[c * N_DIMS + d]);
                if (shift > max_shift) {
                    max_shift = shift;
                }
                centroids[c * N_DIMS + d] = new_val;
            }
        }
```

- **Línea 79** `double max_shift = 0.0;` `[LOOP iteraciones]` — guardará **cuánto se movió el centroide que más se desplazó** en esta iteración. Sirve para el criterio de convergencia. Arranca en 0.
- **Línea 80** `for (int c = 0; c < K; c++) {` `[LOOP iteraciones · LOOP c]` — recorre los `K` clusters para convertir sumas en promedios y medir el desplazamiento. **Serial.**
- **Línea 81** `if (counts[c] == 0) continue;` `[LOOP iteraciones · LOOP c]` — **guarda contra cluster vacío**. Si ningún punto cayó en el cluster `c`, `counts[c]` es 0, y dividir por 0 daría `NaN`/`inf`. El `continue` salta este cluster, dejando su centroide donde estaba. Robustez numérica.
- **Línea 82** `for (int d = 0; d < N_DIMS; d++) {` `[LOOP iteraciones · LOOP c · LOOP d]` — recorre las dimensiones del centroide `c`.
- **Línea 83** `double new_val = new_c[c * N_DIMS + d] / counts[c];` `[LOOP iteraciones · LOOP c · LOOP d]` — **el promedio**: suma de la coordenada dividida por el número de puntos = nueva posición del centroide en esa dimensión. Esto es literalmente "el centro de masa" del cluster.
- **Línea 84** `double shift = fabs(new_val - centroids[c * N_DIMS + d]);` `[LOOP iteraciones · LOOP c · LOOP d]` — cuánto se movió el centroide en esta dimensión: valor absoluto (`fabs`) de la diferencia entre la posición nueva y la vieja. `fabs` porque solo nos importa la magnitud, no la dirección.
- **Línea 85** `if (shift > max_shift) {` `[LOOP iteraciones · LOOP c · LOOP d]` — ¿este desplazamiento supera al máximo registrado?
- **Línea 86** `max_shift = shift;` `[LOOP iteraciones · LOOP c · LOOP d]` — sí: actualiza el máximo. Al final del doble bucle, `max_shift` es el mayor movimiento de cualquier coordenada de cualquier centroide.
- **Línea 87** `}` `[LOOP iteraciones · LOOP c · LOOP d]` — cierra el `if`.
- **Línea 88** `centroids[c * N_DIMS + d] = new_val;` `[LOOP iteraciones · LOOP c · LOOP d]` — **actualiza el centroide** con su nueva posición. A partir de la siguiente iteración, las distancias se medirán contra estos centroides movidos.
- **Línea 89** `}` `[LOOP iteraciones · LOOP c]` — cierra el bucle de dimensiones.
- **Línea 90** `}` `[LOOP iteraciones]` — cierra el bucle de clusters.

```c
        if (max_shift < THRESHOLD) break;
    }

    return omp_get_wtime() - t_start;
}
```

- **Línea 92** `if (max_shift < THRESHOLD) break;` `[LOOP iteraciones]` — **criterio de convergencia**. Si el centroide que más se movió se desplazó menos que `THRESHOLD` (0.0001), los centros prácticamente ya no cambian: el algoritmo ha convergido. `break` sale del bucle de iteraciones antes de agotar las 100. **Implicación**: en datos bien separados (como estos), K-Means converge en pocas iteraciones, ahorrando trabajo.
- **Línea 93** `}` `[SERIAL]` — cierra el bucle de iteraciones (`LOOP iteraciones`).
- **Línea 94** (en blanco).
- **Línea 95** `return omp_get_wtime() - t_start;` `[SERIAL]` — toma una nueva marca de tiempo y le resta la de inicio (línea 42): devuelve el **tiempo total** que tardó el algoritmo, en segundos. Este valor es lo que `main` usa para calcular el speedup.
- **Línea 96** `}` `[GLOBAL]` — cierra la función `run_kmeans`.

### 1.6 Función main (líneas 98-160)

```c
int main(void) {
    int sizes[]   = {10000, 100000, 500000, 1000000};
    int threads[] = {1, 2, 4};
    int n_sizes   = 4;
    int n_threads = 3;
```

- **Línea 98** `int main(void) {` `[SERIAL]` — punto de entrada del programa. `void` = no recibe argumentos de línea de comandos. Devuelve un `int` (0 = éxito).
- **Línea 99** `int sizes[] = {10000, 100000, 500000, 1000000};` `[SERIAL]` — los **tamaños de dataset** a probar: 10 mil, 100 mil, 500 mil y 1 millón de puntos. Probar varios tamaños permite ver cómo el speedup mejora al crecer `N` (con poco `N` el overhead de los hilos domina; con mucho `N`, el speedup se acerca al ideal).
- **Línea 100** `int threads[] = {1, 2, 4};` `[SERIAL]` — los **números de hilos** a probar. El caso de 1 hilo es la **línea base** para medir speedup.
- **Línea 101** `int n_sizes = 4;` `[SERIAL]` — cuántos tamaños hay en `sizes` (4). Se usa como límite de bucle.
- **Línea 102** `int n_threads = 3;` `[SERIAL]` — cuántos valores hay en `threads` (3).

```c
    FILE *fout = fopen("resultados/benchmark.dat", "w");
    if (!fout) {
        fprintf(stderr, "Error: no se pudo crear resultados/benchmark.dat\n");
        return 1;
    }
    fprintf(fout, "# N  hilos  tiempo_s  speedup\n");
```

- **Línea 104** `FILE *fout = fopen("resultados/benchmark.dat", "w");` `[SERIAL]` — abre (o crea) el archivo de salida en modo **escritura** (`"w"`, que **trunca** el archivo si ya existía). `fout` es el manejador del archivo. **Nota**: la carpeta `resultados/` debe existir previamente, o `fopen` falla.
- **Línea 105** `if (!fout) {` `[SERIAL]` — comprueba si `fopen` falló. Si devolvió `NULL` (carpeta inexistente, sin permisos...), `!fout` es verdadero.
- **Línea 106** `fprintf(stderr, "Error: no se pudo crear...\n");` `[SERIAL]` — imprime el error en `stderr` (flujo de errores, separado de la salida normal). Buena práctica: los errores van a `stderr`, no a `stdout`.
- **Línea 107** `return 1;` `[SERIAL]` — termina el programa con código 1 (error). No tiene sentido seguir sin dónde escribir resultados.
- **Línea 108** `}` `[SERIAL]` — cierra el `if`.
- **Línea 109** `fprintf(fout, "# N  hilos  tiempo_s  speedup\n");` `[SERIAL]` — escribe la **cabecera** del archivo de datos. El `#` al inicio es una convención: Gnuplot trata las líneas que empiezan con `#` como comentarios y las ignora al graficar. Documenta las 4 columnas.

```c
    printf("K-Means Paralelo con OpenMP\n");
    printf("Dataset sintetico gaussiano: D=%d dimensiones, K=%d clusters\n\n",
           N_DIMS, K);
    printf("%-10s %-8s %-14s %-8s\n", "N", "Hilos", "Tiempo (s)", "Speedup");
    printf("%-10s %-8s %-14s %-8s\n",
           "----------", "--------", "--------------", "--------");
```

- **Línea 111** `printf("K-Means Paralelo con OpenMP\n");` `[SERIAL]` — imprime el título en pantalla (`stdout`). `\n` es salto de línea.
- **Líneas 112-113** `printf("Dataset sintetico gaussiano: D=%d ..., K=%d ...\n\n", N_DIMS, K);` `[SERIAL]` — imprime los parámetros del experimento. `%d` se sustituye por los enteros `N_DIMS` (16) y `K` (10). El `\n\n` final deja una línea en blanco.
- **Línea 114** `printf("%-10s %-8s %-14s %-8s\n", "N", "Hilos", "Tiempo (s)", "Speedup");` `[SERIAL]` — imprime la **cabecera de la tabla**. Los formatos `%-10s`, `%-8s`, etc. significan: cadena (`s`) alineada a la **izquierda** (`-`) en un ancho fijo (10, 8, 14, 8 caracteres). Esto crea columnas alineadas.
- **Líneas 115-116** `printf("%-10s ...\n", "----------", ...);` `[SERIAL]` — imprime una **línea separadora** de guiones bajo la cabecera, usando los mismos anchos para que quede alineada.

```c
    for (int si = 0; si < n_sizes; si++) {
        int n = sizes[si];
        printf("  [generando N=%d puntos...]\n", n);
```

- **Línea 118** `for (int si = 0; si < n_sizes; si++) {` `[LOOP si]` — **bucle externo sobre los tamaños** de dataset. `si` va de 0 a 3. Todo lo que sigue, hasta la línea 154, se repite una vez por cada tamaño.
- **Línea 119** `int n = sizes[si];` `[LOOP si]` — toma el tamaño actual (p. ej. 10000 en la primera vuelta).
- **Línea 120** `printf("  [generando N=%d puntos...]\n", n);` `[LOOP si]` — avisa por pantalla que se está generando este dataset (puede tardar para N grande, así el usuario sabe que el programa no se colgó).

```c
        double *data        = malloc((size_t)n * N_DIMS * sizeof(double));
        int    *assignments = malloc((size_t)n * sizeof(int));
        double *centroids   = malloc((size_t)K * N_DIMS * sizeof(double));
        double *init_c      = malloc((size_t)K * N_DIMS * sizeof(double));
```

- **Línea 122** `double *data = malloc((size_t)n * N_DIMS * sizeof(double));` `[LOOP si]` — reserva memoria dinámica para los puntos: `n × 16` doubles. Para N=1M son `16M × 8 bytes = 128 MB`. **Por qué `malloc` y no un arreglo en el stack**: el stack es pequeño (típicamente ~8 MB); 128 MB **solo caben en el heap**.
  - **El cast `(size_t)` es crítico aquí**: sin él, `n * N_DIMS` se calcularía en `int` (32 bits). Para N grande eso podría **desbordar** (overflow) y dar un tamaño negativo o diminuto → corrupción de memoria. Al castear `n` a `size_t` (entero sin signo, 64 bits en esta plataforma), toda la multiplicación se hace en 64 bits, sin riesgo.
- **Línea 123** `int *assignments = malloc((size_t)n * sizeof(int));` `[LOOP si]` — reserva `n` enteros para guardar el cluster de cada punto. Mismo cuidado con `(size_t)`.
- **Línea 124** `double *centroids = malloc((size_t)K * N_DIMS * sizeof(double));` `[LOOP si]` — reserva los `K × N_DIMS` centroides de trabajo (pequeño: 1280 bytes).
- **Línea 125** `double *init_c = malloc((size_t)K * N_DIMS * sizeof(double));` `[LOOP si]` — reserva la copia de centroides iniciales (la que se preserva intacta entre corridas).

```c
        if (!data || !assignments || !centroids || !init_c) {
            fprintf(stderr, "Error: sin memoria suficiente para N=%d\n", n);
            return 1;
        }
```

- **Línea 127** `if (!data || !assignments || !centroids || !init_c) {` `[LOOP si]` — comprueba que **los cuatro** `malloc` tuvieron éxito. Si alguno devolvió `NULL` (sin memoria), la condición es verdadera. Verificar siempre el retorno de `malloc` es buena práctica defensiva.
- **Línea 128** `fprintf(stderr, "Error: sin memoria suficiente para N=%d\n", n);` `[LOOP si]` — informa del fallo de memoria por `stderr`.
- **Línea 129** `return 1;` `[LOOP si]` — aborta el programa con código de error.
- **Línea 130** `}` `[LOOP si]` — cierra el `if`.

```c
        /* Generar datos y usar los primeros K puntos como centroides iniciales */
        generate_dataset(n, data);
        memcpy(init_c, data, (size_t)K * N_DIMS * sizeof(double));
```

- **Línea 132** `/* Generar datos... */` `[LOOP si]` — comentario: genera los datos y usa los primeros K puntos como centroides de arranque.
- **Línea 133** `generate_dataset(n, data);` `[LOOP si]` — llama a la función de la sección 1.4 para rellenar `data` con los `n` puntos gaussianos.
- **Línea 134** `memcpy(init_c, data, (size_t)K * N_DIMS * sizeof(double));` `[LOOP si]` — copia los **primeros K puntos** del dataset como centroides iniciales. Es una estrategia de inicialización simple ("Forgy": elegir puntos reales como centros de arranque). Como `data` está ordenado por cluster (punto 0→cluster 0, 1→1, ...), los primeros 10 puntos caen casualmente en 10 clusters distintos, una buena semilla. Estos centros se guardan en `init_c` para que las tres corridas (1/2/4 hilos) arranquen idénticas.

```c
        double t1 = 0.0;
        for (int ti = 0; ti < n_threads; ti++) {
            double t = run_kmeans(n, threads[ti], data, assignments,
                                  centroids, init_c);
            if (ti == 0) t1 = t;
            double speedup = (t > 0.0) ? t1 / t : 1.0;
```

- **Línea 136** `double t1 = 0.0;` `[LOOP si]` — variable que guardará el tiempo con **1 hilo** (la base para el speedup). Se inicializa en 0; se fija en la primera vuelta del bucle siguiente.
- **Línea 137** `for (int ti = 0; ti < n_threads; ti++) {` `[LOOP si · LOOP ti]` — **bucle interno sobre el número de hilos** (1, 2, 4). Para el mismo dataset, corre K-Means tres veces con distinto paralelismo y compara.
- **Líneas 138-139** `double t = run_kmeans(n, threads[ti], data, assignments, centroids, init_c);` `[LOOP si · LOOP ti]` — ejecuta una corrida completa de K-Means con `threads[ti]` hilos y captura el tiempo `t`. Recuerda que `run_kmeans` **reinicia los centroides** desde `init_c` al empezar (línea 39), garantizando una comparación justa.
- **Línea 140** `if (ti == 0) t1 = t;` `[LOOP si · LOOP ti]` — en la **primera** corrida (1 hilo), guarda su tiempo como base `t1`. Las demás corridas se compararán contra este.
- **Línea 141** `double speedup = (t > 0.0) ? t1 / t : 1.0;` `[LOOP si · LOOP ti]` — calcula el **speedup** = tiempo_base / tiempo_actual. Un speedup de 2.0 significa "el doble de rápido".
  - El operador ternario `(t > 0.0) ? t1/t : 1.0` es una **guarda contra división por cero**: si `t` fuera 0 (corrida instantánea, posible con N pequeño y reloj de baja resolución), usa 1.0 en lugar de dividir por cero (que daría `inf`/`NaN`).

```c
            printf("%-10d %-8d %-14.6f %-8.3f\n",
                   n, threads[ti], t, speedup);
            fprintf(fout, "%d %d %.6f %.4f\n",
                    n, threads[ti], t, speedup);
        }
        printf("\n");
```

- **Líneas 143-144** `printf("%-10d %-8d %-14.6f %-8.3f\n", n, threads[ti], t, speedup);` `[LOOP si · LOOP ti]` — imprime una fila de la tabla en pantalla, alineada con la cabecera. Formatos: `%-10d` entero alineado a la izquierda en 10; `%-14.6f` flotante con 6 decimales en ancho 14 (el tiempo); `%-8.3f` flotante con 3 decimales (el speedup).
- **Líneas 145-146** `fprintf(fout, "%d %d %.6f %.4f\n", n, threads[ti], t, speedup);` `[LOOP si · LOOP ti]` — escribe **la misma fila** en el archivo `benchmark.dat`, pero en formato simple separado por espacios (sin anchos fijos), que es lo que Gnuplot espera leer. Las 4 columnas coinciden con la cabecera de la línea 109.
- **Línea 147** `}` `[LOOP si]` — cierra el bucle de hilos (`LOOP ti`).
- **Línea 148** `printf("\n");` `[LOOP si]` — imprime una línea en blanco para separar visualmente los bloques de cada tamaño en la salida.

```c
        free(data);
        free(assignments);
        free(centroids);
        free(init_c);
    }
```

- **Líneas 150-153** `free(data); free(assignments); free(centroids); free(init_c);` `[LOOP si]` — **libera la memoria** de los cuatro buffers. Imprescindible: estamos dentro del bucle de tamaños, así que en la siguiente vuelta se reservará memoria para el siguiente N. Sin estos `free`, la memoria de cada tamaño quedaría "fugada" (memory leak) y se acumularía (128 MB + ... → el programa se quedaría sin RAM). Liberar al final de cada iteración mantiene el uso de memoria acotado.
- **Línea 154** `}` `[SERIAL]` — cierra el bucle de tamaños (`LOOP si`).

```c
    fclose(fout);
    printf("Benchmark exportado a: resultados/benchmark.dat\n");
    printf("Ejecutar 'gnuplot graficar.gp' para generar las graficas.\n");
    return 0;
}
```

- **Línea 156** `fclose(fout);` `[SERIAL]` — cierra el archivo de salida. Esto **vacía el buffer** (flush) asegurando que todo lo escrito quede en disco, y libera el manejador. Olvidar `fclose` puede dejar datos sin escribir.
- **Línea 157** `printf("Benchmark exportado a: resultados/benchmark.dat\n");` `[SERIAL]` — confirma al usuario dónde quedaron los resultados.
- **Línea 158** `printf("Ejecutar 'gnuplot graficar.gp' para generar las graficas.\n");` `[SERIAL]` — indica el siguiente paso para producir las gráficas.
- **Línea 159** `return 0;` `[SERIAL]` — termina el programa con código 0 = **éxito**. Por convención, 0 significa "todo bien".
- **Línea 160** `}` `[GLOBAL]` — cierra `main` y el archivo.

---

## 2. `Makefile`

El Makefile automatiza la compilación y ofrece atajos. `make <target>` ejecuta la receta de ese objetivo. **Importante**: en un Makefile, las líneas de comando de cada receta deben ir indentadas con un **TAB real**, no espacios.

```makefile
CC     = gcc
CFLAGS = -O2 -fopenmp
LIBS   = -lm
TARGET = kmeans
```

- **Línea 1** `CC = gcc` — define la variable `CC` (compilador) como `gcc`. `CC` es el nombre convencional del compilador de C en Make.
- **Línea 2** `CFLAGS = -O2 -fopenmp` — flags de compilación: `-O2` (optimización nivel 2) y `-fopenmp` (activa OpenMP). Reúne en una variable lo explicado en la cabecera de `kmeans.c`.
- **Línea 3** `LIBS = -lm` — librerías a enlazar: `-lm` (librería matemática). Separada de `CFLAGS` a propósito, porque debe ir **al final** del comando (ver línea 6).
- **Línea 4** `TARGET = kmeans` — nombre del ejecutable a producir.

```makefile
# -lm debe ir DESPUES del fuente: GCC resuelve librerias de izq a der
all: $(TARGET)
```

- **Línea 5** `# -lm debe ir DESPUES del fuente...` — comentario (en Make, `#` inicia comentario). Explica una sutileza del enlazador: GCC resuelve símbolos **de izquierda a derecha**. Si pusieras `-lm` antes del fuente, cuando el enlazador procesa `-lm` aún no sabe que `sqrt`/`log` se necesitan, y descartaría la librería → error "undefined reference". Por eso `LIBS` va al final en la línea 10.
- **Línea 6** (continúa el comentario, ya contado arriba como parte del bloque).
- **Línea 7** `all: $(TARGET)` — el target **`all`** (el objetivo por defecto, el que corre `make` sin argumentos) **depende** de `$(TARGET)` (=`kmeans`). No tiene comandos propios; solo dice "para hacer `all`, primero haz `kmeans`". `$(TARGET)` expande la variable a `kmeans`.

```makefile
$(TARGET): src/kmeans.c
	$(CC) $(CFLAGS) -o $(TARGET) src/kmeans.c $(LIBS)
```

- **Línea 9** `$(TARGET): src/kmeans.c` — regla para construir `kmeans`. **Depende** de `src/kmeans.c`: si el fuente cambió desde la última compilación, Make rehace el target; si no, lo deja como está (compilación incremental).
- **Línea 10** `	$(CC) $(CFLAGS) -o $(TARGET) src/kmeans.c $(LIBS)` — **el comando de compilación** (indentado con TAB). Tras expandir las variables queda: `gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm`. Nota que `$(LIBS)` (`-lm`) va **al final**, después del fuente, por lo explicado en la línea 5.

```makefile
run: $(TARGET)
	./$(TARGET)
```

- **Línea 12** `run: $(TARGET)` — target `run`, que depende de `kmeans` (así `make run` compila primero si hace falta).
- **Línea 13** `	./$(TARGET)` — ejecuta el binario (`./kmeans`). El `./` indica "en el directorio actual".

```makefile
# Genera las graficas (requiere gnuplot y que ./kmeans ya se haya ejecutado)
plots: resultados/benchmark.dat
	gnuplot graficar.gp
```

- **Línea 15** `# Genera las graficas...` — comentario explicando el target `plots`.
- **Línea 16** `plots: resultados/benchmark.dat` — target `plots`, que depende del archivo de datos. Si `benchmark.dat` no existe, Make intentará generarlo con la regla de la línea 22.
- **Línea 17** `	gnuplot graficar.gp` — ejecuta Gnuplot con el script, produciendo los PNG.

```makefile
# Ejecuta el programa y luego genera las graficas en un solo paso
all-plots: run plots
```

- **Línea 19** `# Ejecuta el programa...` — comentario.
- **Línea 20** `all-plots: run plots` — target de conveniencia que **encadena** `run` (ejecutar) y luego `plots` (graficar). `make all-plots` hace todo el flujo de un tirón. No tiene comandos propios; solo dependencias en orden.

```makefile
resultados/benchmark.dat: $(TARGET)
	./$(TARGET)
```

- **Línea 22** `resultados/benchmark.dat: $(TARGET)` — regla que dice cómo **producir** el archivo de datos: depende del ejecutable.
- **Línea 23** `	./$(TARGET)` — para generar `benchmark.dat`, ejecuta `./kmeans` (que lo escribe). Así, si pides `plots` sin haber corrido el programa, Make sabe ejecutarlo automáticamente primero.

```makefile
# ── Compilar informe IEEE ────────────────────────────────────
informe:
	cd informe && pdflatex informe_ieee.tex && pdflatex informe_ieee.tex
```

- **Línea 25** `# ── Compilar informe IEEE ──...` — comentario separador de sección.
- **Línea 26** `informe:` — target `informe`, **sin dependencias** (siempre se ejecuta cuando lo invocas).
- **Línea 27** `	cd informe && pdflatex informe_ieee.tex && pdflatex informe_ieee.tex` — entra a la carpeta `informe/` y compila el LaTeX **dos veces**. **Por qué dos veces**: la primera pasada genera las referencias cruzadas, índices y citas; la segunda las **resuelve** correctamente (LaTeX necesita dos pasadas para que las referencias internas queden bien). El `&&` encadena: cada comando solo corre si el anterior tuvo éxito.

```makefile
referencias:
	cd informe && pdflatex referencias.tex
```

- **Línea 29** `referencias:` — target para compilar un documento de referencias aparte.
- **Línea 30** `	cd informe && pdflatex referencias.tex` — compila `referencias.tex` (una sola pasada).

```makefile
clean:
	rm -f $(TARGET)
	rm -f resultados/benchmark.dat
	rm -f resultados/speedup.png resultados/scaling.png
	rm -f informe/*.aux informe/*.log informe/*.out
```

- **Línea 32** `clean:` — target `clean`, convención para **borrar archivos generados** y dejar el proyecto limpio.
- **Línea 33** `	rm -f $(TARGET)` — borra el ejecutable. `-f` (force) evita errores si el archivo no existe.
- **Línea 34** `	rm -f resultados/benchmark.dat` — borra el archivo de datos del benchmark.
- **Línea 35** `	rm -f resultados/speedup.png resultados/scaling.png` — borra las dos gráficas generadas.
- **Línea 36** `	rm -f informe/*.aux informe/*.log informe/*.out` — borra los archivos auxiliares que LaTeX genera al compilar (`.aux`, `.log`, `.out`). El `*` es comodín.

```makefile
.PHONY: all run plots all-plots clean informe referencias
```

- **Línea 38** `.PHONY: all run plots all-plots clean informe referencias` — declara estos targets como **"phony" (falsos)**: nombres de acciones, no archivos reales. **Qué problema resuelve**: si existiera por casualidad un archivo llamado `clean` o `run` en la carpeta, Make pensaría que el target ya está "hecho" y no ejecutaría sus comandos. `.PHONY` le dice "estos no son archivos, ejecútalos siempre". (Nota: `$(TARGET)` y `resultados/benchmark.dat` **no** se listan aquí porque sí son archivos reales que Make debe construir según su fecha.)
- **Línea 39** (fin del archivo).

---

## 3. `graficar.gp`

Script de Gnuplot. Genera dos gráficas PNG a partir de `resultados/benchmark.dat`. En Gnuplot, `#` inicia comentario y `set` configura opciones del gráfico.

```gnuplot
# graficar.gp — Genera graficas del benchmark K-Means con OpenMP
# Uso: gnuplot graficar.gp
# Requiere: resultados/benchmark.dat (generado por ./kmeans)
#
# benchmark.dat tiene 4 columnas:
#   $1 = N (numero de puntos)
#   $2 = hilos
#   $3 = tiempo en segundos
#   $4 = speedup (T_1 / T_N)
```

- **Líneas 1-9** — bloque de comentarios de cabecera. Documentan el propósito, el uso (`gnuplot graficar.gp`), el prerrequisito (que exista `benchmark.dat`) y el **significado de las 4 columnas** del archivo. En Gnuplot, las columnas se referencian como `$1`, `$2`, etc. Esta leyenda es clave para entender los `using` más abajo.

### Gráfica 1: Speedup vs Número de hilos (líneas 11-60)

```gnuplot
# ===================================================================
# GRAFICA 1: Speedup vs Numero de Hilos (una linea por N)
# Muestra como escala el speedup segun el tamano del dataset.
# Para N pequeno el overhead domina; para N grande el speedup
# se acerca al numero de hilos.
# ===================================================================
```

- **Líneas 11-16** — comentario de sección que describe la primera gráfica: speedup en el eje Y, número de hilos en el X, una línea por cada tamaño N. Resume la conclusión esperada: N pequeño → el overhead de los hilos domina → speedup pobre; N grande → speedup cercano al ideal.

```gnuplot
set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "resultados/speedup.png"
```

- **Línea 18** `set terminal pngcairo size 800,600 enhanced font "Arial,12"` — configura el **formato de salida**: imagen PNG usando el motor `pngcairo` (suaviza líneas), de 800×600 píxeles, con texto mejorado (`enhanced`, permite subíndices/superíndices) y fuente Arial 12.
- **Línea 19** `set output "resultados/speedup.png"` — define el **archivo de salida** de esta primera gráfica. Todo lo que se grafique a partir de aquí irá a `speedup.png`.

```gnuplot
set title "Speedup de K-Means con OpenMP\n(Dataset sintetico, D=16, K=10)" \
          font "Arial,14"
set xlabel "Numero de Hilos"       font "Arial,12"
set ylabel "Speedup (T_1 / T_N)"   font "Arial,12"
```

- **Líneas 21-22** `set title "..." font "Arial,14"` — título de la gráfica, en dos líneas (`\n`). El `\` al final de la línea 21 es **continuación**: une la 21 y la 22 como una sola instrucción.
- **Línea 23** `set xlabel "Numero de Hilos" font "Arial,12"` — etiqueta del eje X.
- **Línea 24** `set ylabel "Speedup (T_1 / T_N)" font "Arial,12"` — etiqueta del eje Y. `T_1` y `T_N` aparecen con subíndice gracias al modo `enhanced`.

```gnuplot
set xrange [0.5:5]
set yrange [0:4.5]
set xtics (1, 2, 4)
set ytics 0.5
set key top left
set grid lc rgb "#DDDDDD"
set border lw 1.5
```

- **Línea 26** `set xrange [0.5:5]` — rango del eje X de 0.5 a 5, dejando margen alrededor de los valores de hilos (1, 2, 4).
- **Línea 27** `set yrange [0:4.5]` — rango del eje Y de 0 a 4.5 (un poco más que el speedup máximo posible con 4 hilos).
- **Línea 28** `set xtics (1, 2, 4)` — marcas del eje X **solo** en 1, 2 y 4 (los valores realmente probados), en vez de marcas automáticas.
- **Línea 29** `set ytics 0.5` — marcas del eje Y cada 0.5 unidades.
- **Línea 30** `set key top left` — coloca la **leyenda** (la caja que identifica cada línea) arriba a la izquierda.
- **Línea 31** `set grid lc rgb "#DDDDDD"` — activa una **cuadrícula** de fondo en gris claro (`#DDDDDD`), para leer mejor los valores.
- **Línea 32** `set border lw 1.5` — dibuja el borde del gráfico con grosor de línea 1.5.

```gnuplot
# Colores: N=10k azul, 100k verde, 500k naranja, 1M rojo
# Filtrado: ($1==N ? $2 : 1/0) retorna col2 si N coincide, NaN si no.
# Gnuplot dibuja NaN como punto ausente, logrando una linea por N.

ideal(x) = x
```

- **Líneas 34-36** — comentarios que anticipan dos cosas: el código de colores por tamaño, y **el truco de filtrado** que se usará en el `plot`. Lo explican: la expresión `($1==N ? $2 : 1/0)` devuelve la columna 2 si la fila es del N buscado, y `1/0` (división por cero = NaN) si no. Gnuplot **no dibuja** los puntos NaN, así que de un archivo con todas las filas mezcladas, esto "extrae" solo las de un N concreto → una línea limpia por cada N.
- **Línea 37** (en blanco).
- **Línea 38** `ideal(x) = x` — define una **función** `ideal` que devuelve su argumento. Representa el **speedup ideal lineal** (con 2 hilos, speedup 2; con 4, speedup 4). Se dibuja como referencia para comparar el speedup real contra el perfecto.

```gnuplot
plot \
  ideal(x) \
      title "Ideal (lineal)" lw 2 lc rgb "#AAAAAA" dt 2, \
  'resultados/benchmark.dat' \
      using ($1==10000   ? $2 : 1/0):($1==10000   ? $4 : 1/0) \
      title "N = 10,000"    with linespoints pt 7 ps 1.4 lw 2 \
      lc rgb "#3498DB", \
  'resultados/benchmark.dat' \
      using ($1==100000  ? $2 : 1/0):($1==100000  ? $4 : 1/0) \
      title "N = 100,000"   with linespoints pt 7 ps 1.4 lw 2 \
      lc rgb "#2ECC71", \
  'resultados/benchmark.dat' \
      using ($1==500000  ? $2 : 1/0):($1==500000  ? $4 : 1/0) \
      title "N = 500,000"   with linespoints pt 7 ps 1.4 lw 2 \
      lc rgb "#E67E22", \
  'resultados/benchmark.dat' \
      using ($1==1000000 ? $2 : 1/0):($1==1000000 ? $4 : 1/0) \
      title "N = 1,000,000" with linespoints pt 7 ps 1.4 lw 2 \
      lc rgb "#E74C3C"
```

- **Línea 40** `plot \` — inicia el comando `plot`, que dibuja todas las curvas. El `\` continúa en la línea siguiente; las curvas se separan por comas. Es **una sola instrucción** repartida en muchas líneas.
- **Líneas 41-42** `ideal(x) ... title "Ideal (lineal)" lw 2 lc rgb "#AAAAAA" dt 2,` — primera curva: la función `ideal`. `title` le pone nombre en la leyenda, `lw 2` grosor 2, `lc rgb "#AAAAAA"` color gris, `dt 2` tipo de trazo discontinuo (dashed). La coma final separa de la siguiente curva.
- **Líneas 43-46** — segunda curva, **N = 10.000**:
  - `'resultados/benchmark.dat'` — el archivo de datos.
  - `using ($1==10000 ? $2 : 1/0):($1==10000 ? $4 : 1/0)` — define X e Y. X = columna 2 (hilos) **solo** si la fila tiene `$1==10000`, si no NaN; Y = columna 4 (speedup) con el mismo filtro. Esto es el **truco** explicado en la línea 35: se queda solo con las filas de N=10000.
  - `title "N = 10,000"` — nombre en la leyenda.
  - `with linespoints` — dibuja línea **y** puntos.
  - `pt 7` tipo de punto (círculo relleno), `ps 1.4` tamaño de punto, `lw 2` grosor de línea, `lc rgb "#3498DB"` color azul.
- **Líneas 47-50** — tercera curva, **N = 100.000**, idéntica estructura pero filtrando `$1==100000` y en verde (`#2ECC71`).
- **Líneas 51-54** — cuarta curva, **N = 500.000**, filtro `$1==500000`, naranja (`#E67E22`).
- **Líneas 55-58** — quinta curva, **N = 1.000.000**, filtro `$1==1000000`, rojo (`#E74C3C`). Sin coma al final: es la última curva del `plot`.

```gnuplot
print "Generado: speedup.png"
```

- **Línea 60** `print "Generado: speedup.png"` — imprime un mensaje de progreso en la terminal, confirmando que la primera gráfica se generó.

### Gráfica 2: Tiempo vs N en escala log-log (líneas 62-103)

```gnuplot
# ===================================================================
# GRAFICA 2: Tiempo de Ejecucion vs N (escala log-log)
# Muestra escalabilidad: como crece el tiempo al aumentar N,
# y cuanto ayuda cada numero de hilos.
# En escala log-log, relacion lineal = crecimiento O(N).
# ===================================================================
```

- **Líneas 62-67** — comentario de la segunda gráfica: tiempo (Y) vs N (X) en escala logarítmica en ambos ejes. Apunta una propiedad útil: en log-log, una **recta** indica crecimiento lineal `O(N)` (el tiempo crece proporcional al número de puntos, como se espera de este algoritmo).

```gnuplot
set output "resultados/scaling.png"
```

- **Línea 69** `set output "resultados/scaling.png"` — **cambia el archivo de salida** a `scaling.png`. A partir de aquí lo que se grafique va al segundo PNG. El `set terminal` de la línea 18 sigue vigente (mismo formato/tamaño).

```gnuplot
set title "Tiempo de Ejecucion vs Tamano del Dataset\n(Dataset sintetico, D=16, K=10)" \
          font "Arial,14"
set xlabel "N (numero de puntos)"  font "Arial,12"
set ylabel "Tiempo (segundos)"     font "Arial,12"
```

- **Líneas 71-72** `set title "..." font "Arial,14"` — nuevo título para esta gráfica (sobrescribe el anterior), en dos líneas.
- **Línea 73** `set xlabel "N (numero de puntos)" font "Arial,12"` — etiqueta del eje X (ahora es N).
- **Línea 74** `set ylabel "Tiempo (segundos)" font "Arial,12"` — etiqueta del eje Y (ahora es tiempo).

```gnuplot
set logscale xy
set format x "10^{%L}"
set format y "%.3f"
set xrange [5000:2000000]
set yrange [*:*]
unset ytics
set ytics auto
set xtics auto
set key top left
set grid lc rgb "#DDDDDD"
```

- **Línea 76** `set logscale xy` — activa **escala logarítmica** en ambos ejes (`x` e `y`). Necesario porque N varía de 10 mil a 1 millón (dos órdenes de magnitud); en escala lineal los puntos pequeños se amontonarían.
- **Línea 77** `set format x "10^{%L}"` — formato de las etiquetas del eje X como potencias de 10 (`10^4`, `10^5`...). `%L` es el exponente logarítmico; el modo `enhanced` lo muestra como superíndice.
- **Línea 78** `set format y "%.3f"` — etiquetas del eje Y con 3 decimales.
- **Línea 79** `set xrange [5000:2000000]` — rango del eje X, con margen alrededor de 10.000–1.000.000.
- **Línea 80** `set yrange [*:*]` — rango Y **automático** (`*` = "que Gnuplot lo decida según los datos").
- **Línea 81** `unset ytics` — borra cualquier configuración previa de marcas del eje Y (heredada de la gráfica 1, donde se había puesto `set ytics 0.5`). Limpieza antes de reconfigurar.
- **Línea 82** `set ytics auto` — vuelve a marcas automáticas en Y, apropiadas para la escala logarítmica.
- **Línea 83** `set xtics auto` — marcas automáticas en X (a diferencia de la gráfica 1 que las fijaba en 1,2,4).
- **Línea 84** `set key top left` — leyenda arriba a la izquierda.
- **Línea 85** `set grid lc rgb "#DDDDDD"` — cuadrícula gris claro.

```gnuplot
# Filtrado por numero de hilos con awk via pipe. Las filas de un mismo
# numero de hilos NO son contiguas en el archivo (estan intercaladas por
# N), lo que rompe las lineas si se filtra con el truco '1/0'. Con awk
# extraemos solo las filas del hilo deseado, quedando contiguas y la
# linea conecta limpio.
```

- **Líneas 87-91** — comentario importante que explica **por qué aquí NO se usa el truco `1/0`** de la gráfica 1, sino `awk`. El motivo: en `benchmark.dat` las filas se escriben agrupadas por N (para cada N, las 3 corridas de hilos seguidas). Por tanto las filas de "1 hilo" están **dispersas** por todo el archivo (una por cada N), no contiguas. Si filtráramos con `1/0`, Gnuplot dejaría huecos NaN entre los puntos válidos y **no conectaría la línea** correctamente. La solución: usar `awk` para extraer **solo** las filas de un número de hilos, que así quedan **contiguas**, y la línea se dibuja limpia.

```gnuplot
plot \
  "< awk '$2==1' resultados/benchmark.dat" using 1:3 \
      title "1 hilo"  with linespoints pt 7 ps 1.4 lw 2 lc rgb "#E74C3C", \
  "< awk '$2==2' resultados/benchmark.dat" using 1:3 \
      title "2 hilos" with linespoints pt 9 ps 1.4 lw 2 lc rgb "#3498DB", \
  "< awk '$2==4' resultados/benchmark.dat" using 1:3 \
      title "4 hilos" with linespoints pt 5 ps 1.4 lw 2 lc rgb "#2ECC71"
```

- **Línea 92** `plot \` — inicia el segundo `plot` (continúa con `\`).
- **Líneas 93-94** — primera curva, **1 hilo**:
  - `"< awk '$2==1' resultados/benchmark.dat"` — el `<` indica que Gnuplot debe **ejecutar un comando** y leer su salida como si fuera un archivo. `awk '$2==1'` filtra y deja solo las filas cuya **columna 2** (hilos) es 1.
  - `using 1:3` — usa columna 1 (N) como X y columna 3 (tiempo) como Y.
  - `title "1 hilo"`, `with linespoints`, `pt 7` (círculo), `ps 1.4`, `lw 2`, `lc rgb "#E74C3C"` rojo.
- **Líneas 95-96** — segunda curva, **2 hilos**: `awk '$2==2'`, `pt 9` (triángulo), color azul (`#3498DB`).
- **Líneas 97-98** — tercera curva, **4 hilos**: `awk '$2==4'`, `pt 5` (cuadrado), color verde (`#2ECC71`). Última, sin coma final.

```gnuplot
print "Generado: scaling.png"
print ""
print "Listo. Abre speedup.png y scaling.png para ver las graficas."
```

- **Línea 100** `print "Generado: scaling.png"` — confirma la segunda gráfica.
- **Línea 101** `print ""` — imprime una línea en blanco (separación visual en la terminal).
- **Línea 102** `print "Listo. Abre speedup.png y scaling.png..."` — mensaje final indicando que todo terminó y qué archivos mirar.
- **Línea 103** (fin del archivo).

---

## Resumen de paralelismo

Para fijar lo esencial sobre **dónde** y **por qué** hay paralelismo:

| Parte | Líneas | ¿Paralela? | Razón |
|-------|--------|-----------|-------|
| Generación de datos | 22-31 | No | Se hace una vez; no es el cuello de botella medido. |
| **Fase 1: asignación** | 47-63 | **Sí** (`#pragma omp parallel for`) | Cada punto es independiente; escribe en `assignments[i]` distinto → **sin race conditions**. Es la parte cara (`n × K × N_DIMS`). |
| Fase 2: acumulación | 71-77 | No | Varios puntos caen en el mismo cluster → escribir `counts[c]`/`new_c` en paralelo daría **race condition**. Y es barata. |
| Fase 2: promedio | 80-90 | No | Solo `K` clusters: trivial, no vale la pena paralelizar. |

**La clave del proyecto**: solo se paraleliza la fase 1 porque es la única (a) costosa y (b) sin dependencias entre iteraciones. El speedup mejora con N grande porque entonces el trabajo paralelo de la fase 1 domina sobre el overhead de crear hilos y sobre la fase 2 serial (Ley de Amdahl).
```
