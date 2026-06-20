# Explicación línea por línea — Criba de Eratóstenes Distribuida con MPI

Esta es una **referencia exhaustiva, línea a línea**, de todo el código del proyecto. A diferencia de `GUIA_ESTUDIO.md` (tutorial conceptual "desde cero"), aquí recorremos **cada línea** de los tres archivos fuente sin saltarnos ninguna, explicando:

- **Qué hace** la línea.
- **Por qué está ahí** (la intención).
- **Sus implicaciones** (rendimiento, memoria, comunicación, riesgos).
- **En qué contexto se ejecuta**: qué procesos la ejecutan (todos / solo rank 0 / solo trabajadores), dentro de qué loop o paso, y si es una llamada de comunicación MPI.

> **Idea central que distingue a este proyecto del de OpenMP**: aquí el paralelismo es **MPI (Message Passing Interface)**. En lugar de hilos que comparten memoria, hay **varios procesos independientes**, cada uno con su **propia memoria separada**, que se coordinan **enviándose mensajes** (`MPI_Send`/`MPI_Recv`). Todos ejecutan el **mismo programa** (modelo SPMD: *Single Program, Multiple Data*) pero se diferencian por su **rank** (identificador). Como no comparten memoria, **no existen las race conditions** del estilo OpenMP; el reto es otro: **repartir el trabajo y comunicar los resultados**.

---

## Tabla de Contenidos

- [Leyenda de etiquetas de contexto](#leyenda-de-etiquetas-de-contexto)
- [Conceptos MPI mínimos](#conceptos-mpi-mínimos)
- [Mapa del proyecto](#mapa-del-proyecto)
- [1. `src/criba.c`](#1-srccribac)
  - [1.1 Cabecera (líneas 1-5)](#11-cabecera-líneas-1-5)
  - [1.2 Includes (líneas 7-10)](#12-includes-líneas-7-10)
  - [1.3 Constante N (líneas 12-13)](#13-constante-n-líneas-12-13)
  - [1.4 PASO 1 — Arranque de MPI (líneas 15-24)](#14-paso-1--arranque-de-mpi-líneas-15-24)
  - [1.5 PASO 2 — Primos base hasta √N (líneas 26-57)](#15-paso-2--primos-base-hasta-n-líneas-26-57)
  - [1.6 PASO 3 — Difundir los primos base (líneas 59-75)](#16-paso-3--difundir-los-primos-base-líneas-59-75)
  - [1.7 PASO 4 — Repartir el rango (líneas 77-87)](#17-paso-4--repartir-el-rango-líneas-77-87)
  - [1.8 PASO 5 — Cribar el segmento local (líneas 89-101)](#18-paso-5--cribar-el-segmento-local-líneas-89-101)
  - [1.9 PASO 6 — Contar primos locales (líneas 103-106)](#19-paso-6--contar-primos-locales-líneas-103-106)
  - [1.10 PASO 7 — Recolectar conteos (líneas 108-118)](#110-paso-7--recolectar-conteos-líneas-108-118)
  - [1.11 PASO 8 — Resumen y resultados.dat (líneas 120-167)](#111-paso-8--resumen-y-resultadosdat-líneas-120-167)
  - [1.12 PASO 9 — Recolectar y escribir primos.dat (líneas 169-213)](#112-paso-9--recolectar-y-escribir-primosdat-líneas-169-213)
  - [1.13 Cierre (líneas 215-219)](#113-cierre-líneas-215-219)
- [2. `Makefile`](#2-makefile)
- [3. `graficar.gp`](#3-graficargp)
- [Resumen de la comunicación MPI](#resumen-de-la-comunicación-mpi)

---

## Leyenda de etiquetas de contexto

A lo largo del documento verás estas etiquetas antes de explicar una línea, para que siempre sepas **quién** la ejecuta y **dónde**:

| Etiqueta | Significado |
|----------|-------------|
| `[TODOS]` | Lo ejecutan **todos** los procesos (cada uno con su copia). |
| `[RANK 0]` | Solo el proceso **maestro** (rank 0). |
| `[RANK != 0]` | Solo los procesos **trabajadores** (rank 1, 2, ...). |
| `[MPI]` | Llamada de comunicación o entorno MPI. |
| `[GLOBAL]` | Fuera de la función (preprocesador). |
| `[LOOP ...]` | Dentro de un bucle (se nombra cuál: criba base, segmento, primos...). |

Cuando una línea cumple varias, se combinan: p. ej. `[RANK 0 · LOOP dest · MPI]` = dentro de la rama del maestro, en el bucle sobre destinos, ejecutando una llamada MPI.

---

## Conceptos MPI mínimos

Para seguir el código conviene tener claros cinco términos:

- **Proceso / rank**: cada copia del programa en ejecución. `rank` es su número (0, 1, 2, ...). `size` es cuántos hay en total. Se lanzan con `mpirun -np 4 ./criba` (4 procesos).
- **SPMD**: todos corren el **mismo código**; las diferencias de comportamiento se logran con `if (rank == 0) ... else ...`.
- **Memoria separada**: una variable `x` en el rank 0 **no es** la misma que `x` en el rank 1. Por eso para compartir datos hay que **enviarlos explícitamente**.
- **`MPI_Send` / `MPI_Recv`**: comunicación **punto a punto** (un emisor, un receptor). Cada par lleva un **tag** (etiqueta entera) que permite distinguir mensajes distintos entre la misma pareja.
- **Descomposición de dominio**: la estrategia de este programa. El rango `[2, N]` se parte en trozos (`[low, high]`), uno por proceso; cada uno criba su trozo en paralelo.

---

## Mapa del proyecto

| Archivo | Qué hace |
|---------|----------|
| `src/criba.c` | El programa MPI: reparte el rango `[2, N]` entre procesos, cada uno criba su trozo, y el rank 0 recolecta conteos y la lista de primos. |
| `Makefile` | Compila con `mpicc` y ofrece atajos (`make run`, `make benchmark`, `make graficas`, `make clean`). |
| `graficar.gp` | Script de Gnuplot: dibuja la distribución de primos por proceso y la curva de speedup. |

Flujo: `make` compila → `make run` (o `mpirun -np 4 ./criba`) genera `resultados/resultados.dat` y `resultados/primos.dat` → `make benchmark` mide 1/2/4 procesos y crea `speedup.dat` → `make graficas` produce los PNG.

---

## 1. `src/criba.c`

### 1.1 Cabecera (líneas 1-5)

```c
/*
 * criba.c — Criba de Eratostenes distribuida con MPI.
 * Compila:  mpicc -O2 -Wall -o criba criba.c -lm
 * Ejecuta:  mpirun -np 4 ./criba
 */
```

`[GLOBAL]`

- **Línea 1** `/*` — abre un comentario de bloque (lo ignora el compilador).
- **Línea 2** `* criba.c — ...` — nombre y propósito: implementa la Criba de Eratóstenes (algoritmo clásico para hallar primos) repartida con MPI.
- **Línea 3** `* Compila: mpicc -O2 -Wall -o criba criba.c -lm` — documenta cómo compilar:
  - `mpicc` es el **wrapper** del compilador para MPI: por dentro llama a `gcc` pero añade automáticamente las rutas y librerías de MPI. No se usa `gcc` directamente.
  - `-O2` optimización nivel 2.
  - `-Wall` activa **todos los warnings** comunes (buena práctica para detectar errores).
  - `-lm` enlaza la librería matemática (se usa `sqrt`), y va al final.
- **Línea 4** `* Ejecuta: mpirun -np 4 ./criba` — cómo lanzarlo: `mpirun` arranca varios procesos; `-np 4` pide **4 procesos**. Cambiar ese número cambia el paralelismo.
- **Línea 5** `*/` — cierra el comentario.

### 1.2 Includes (líneas 7-10)

```c
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
```

`[GLOBAL]`

- **Línea 7** `#include <mpi.h>` — la cabecera de MPI. Aporta `MPI_Init`, `MPI_Send`, `MPI_Recv`, los tipos (`MPI_INT`, `MPI_LONG_LONG`), constantes (`MPI_COMM_WORLD`, `MPI_STATUS_IGNORE`) y funciones (`MPI_Wtime`). Por convención se pone **primero**.
- **Línea 8** `#include <stdio.h>` — E/S estándar: `printf`, `fprintf`, `fopen`, `fclose`, `FILE`.
- **Línea 9** `#include <stdlib.h>` — memoria dinámica y utilidades: `malloc`, `calloc`, `free`.
- **Línea 10** `#include <math.h>` — funciones matemáticas: aquí solo `sqrt` (para calcular √N). Por esto hace falta `-lm` al enlazar.

### 1.3 Constante N (líneas 12-13)

```c
/* Limite superior de la busqueda */
#define N 10000000
```

`[GLOBAL]`

- **Línea 12** `/* Limite superior de la busqueda */` — comentario: `N` es hasta dónde buscamos primos.
- **Línea 13** `#define N 10000000` — define `N` = **10 millones**. Macro de preprocesador (sustitución textual, sin coste en ejecución). **Implicación clave**: con N tan grande, los **índices y sumas** pueden superar el rango de un `int` de 32 bits (~2.100 millones); por eso a lo largo del código se usan `long long` (64 bits) y casts `(long long)` en las multiplicaciones, para evitar **overflow**. Subir N a 10^9 multiplicaría memoria y tiempo, y haría `primos.dat` enorme (ver PASO 9).

### 1.4 PASO 1 — Arranque de MPI (líneas 15-24)

```c
int main(int argc, char *argv[])
{
    int rank, size;

    /* PASO 1: iniciar MPI (rank = id del proceso, size = total de procesos) */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double t_inicio = MPI_Wtime();
```

- **Línea 15** `int main(int argc, char *argv[])` `[TODOS]` — entrada del programa. Recibe `argc`/`argv` (argumentos de línea de comandos) porque **MPI los necesita** (se pasan a `MPI_Init`). Recuerda: este `main` lo ejecutan **todos** los procesos a la vez.
- **Línea 16** `{` `[TODOS]` — abre el cuerpo de `main`.
- **Línea 17** `int rank, size;` `[TODOS]` — declara dos enteros: `rank` (el id de **este** proceso) y `size` (cuántos procesos hay en total). Cada proceso tendrá su propio valor de `rank`, pero el mismo `size`.
- **Línea 19** `/* PASO 1: iniciar MPI ... */` `[TODOS]` — comentario que rotula el primer paso.
- **Línea 20** `MPI_Init(&argc, &argv);` `[TODOS · MPI]` — **inicializa el entorno MPI**. Debe ser la primera llamada MPI. A partir de aquí los procesos pueden comunicarse. Se le pasan las direcciones de `argc`/`argv` para que MPI procese sus propios argumentos internos.
- **Línea 21** `MPI_Comm_rank(MPI_COMM_WORLD, &rank);` `[TODOS · MPI]` — pregunta **"¿cuál es mi número?"** y lo guarda en `rank`. `MPI_COMM_WORLD` es el **comunicador** que agrupa a **todos** los procesos. Tras esta línea, cada proceso conoce su identidad (0, 1, 2, ...).
- **Línea 22** `MPI_Comm_size(MPI_COMM_WORLD, &size);` `[TODOS · MPI]` — pregunta **"¿cuántos somos?"** y lo guarda en `size`. Con `mpirun -np 4`, `size` será 4 en todos los procesos.
- **Línea 24** `double t_inicio = MPI_Wtime();` `[TODOS]` — marca el tiempo de inicio (reloj de pared, en segundos). `MPI_Wtime` es el cronómetro de MPI. Cada proceso mide su propio tiempo; el que importa para el informe es el del rank 0 (línea 120).

### 1.5 PASO 2 — Primos base hasta √N (líneas 26-57)

**Idea**: para cribar el rango `[2, N]` por trozos, basta con conocer los primos hasta **√N**. ¿Por qué? Porque todo número compuesto `c ≤ N` tiene al menos un factor primo `≤ √c ≤ √N`. Así que esos "primos base" son los únicos que hacen falta para marcar compuestos en cualquier segmento. **Solo el rank 0** los calcula (con una criba clásica pequeña) y luego los reparte.

```c
    /* PASO 2: rank 0 calcula los primos base hasta sqrt(N) (los unicos que
     * pueden marcar compuestos en cualquier segmento) */
    int raiz_n = (int)sqrt((double)N);
    int base_count = 0;
    int *base_primes = NULL;
```

- **Líneas 26-27** `/* PASO 2: ... */` `[TODOS]` — comentario que explica la justificación matemática (los primos hasta √N bastan).
- **Línea 28** `int raiz_n = (int)sqrt((double)N);` `[TODOS]` — calcula **√N** y lo trunca a entero. `(double)N` convierte N a flotante para `sqrt`; `(int)` recorta los decimales. Para N=10^7, `raiz_n` ≈ 3162. **Lo calculan todos** (es barato y cada proceso lo necesitará para cribar), aunque la criba base en sí solo la hace el rank 0.
- **Línea 29** `int base_count = 0;` `[TODOS]` — contador de cuántos primos base hay. Arranca en 0. En el rank 0 se llenará calculándolo; en los demás, se recibirá por mensaje (PASO 3).
- **Línea 30** `int *base_primes = NULL;` `[TODOS]` — puntero al arreglo de primos base, inicializado a `NULL`. Cada proceso lo reservará: el rank 0 al calcularlos (línea 49), los trabajadores al recibirlos (línea 70).

```c
    if (rank == 0) {
        char *marcado = calloc(raiz_n + 1, sizeof(char)); /* 0=primo, 1=compuesto */
        if (!marcado) { MPI_Abort(MPI_COMM_WORLD, 1); }

        marcado[0] = 1; /* 0 y 1 no son primos */
        marcado[1] = 1;
```

- **Línea 32** `if (rank == 0) {` `[RANK 0]` — **abre la rama exclusiva del maestro**. Todo lo de dentro (hasta la línea 57) lo ejecuta **solo** el rank 0. Aquí se ve el patrón SPMD: el mismo código, distinto comportamiento según el rank.
- **Línea 33** `char *marcado = calloc(raiz_n + 1, sizeof(char));` `[RANK 0]` — reserva un arreglo de `raiz_n + 1` bytes para marcar números del 0 al `raiz_n`. **`calloc`** (a diferencia de `malloc`) **inicializa todo a cero**. Convención: `0 = primo`, `1 = compuesto`. Se usa `char` (1 byte) en vez de `int` para ahorrar memoria: solo necesitamos un flag.
- **Línea 34** `if (!marcado) { MPI_Abort(MPI_COMM_WORLD, 1); }` `[RANK 0 · MPI]` — si `calloc` falló (sin memoria), aborta. **`MPI_Abort`** mata **todos** los procesos del comunicador (no solo este): en MPI no sirve que un proceso siga si otro murió, porque quedarían esperando mensajes que nunca llegan (deadlock). El `1` es el código de error.
- **Línea 36** `marcado[0] = 1;` `[RANK 0]` — marca el 0 como compuesto (no es primo por definición).
- **Línea 37** `marcado[1] = 1;` `[RANK 0]` — marca el 1 como compuesto (tampoco es primo). El comentario de la línea 36 cubre ambos casos.

```c
        for (int i = 2; (long long)i * i <= raiz_n; i++) {
            if (marcado[i] == 0) {
                for (int j = i * i; j <= raiz_n; j += i)
                    marcado[j] = 1;
            }
        }
```

Este es el **corazón de la criba clásica**: tachar los múltiplos de cada primo.

- **Línea 39** `for (int i = 2; (long long)i * i <= raiz_n; i++) {` `[RANK 0 · LOOP criba-base i]` — recorre los candidatos desde 2. La condición `(long long)i * i <= raiz_n` significa "mientras `i² ≤ raiz_n`": basta cribar hasta la raíz de `raiz_n`, porque los compuestos mayores ya quedarán tachados por un factor menor. El cast `(long long)` evita que `i*i` desborde el `int` (defensivo).
- **Línea 40** `if (marcado[i] == 0) {` `[RANK 0 · LOOP criba-base i]` — si `i` **sigue sin marcar**, es primo → hay que tachar sus múltiplos. Si ya estaba marcado, se salta (sus múltiplos ya los tachó un primo menor).
- **Línea 41** `for (int j = i * i; j <= raiz_n; j += i)` `[RANK 0 · LOOP criba-base i · LOOP j]` — tacha los múltiplos de `i`. **Empieza en `i*i`** (no en `2*i`) porque los múltiplos menores (`2i`, `3i`, ..., `(i-1)i`) ya fueron tachados por primos más pequeños. Avanza de `i` en `i`. Es una **optimización estándar** de la criba.
- **Línea 42** `marcado[j] = 1;` `[RANK 0 · LOOP criba-base i · LOOP j]` — marca `j` como compuesto.
- **Línea 43** `}` `[RANK 0 · LOOP criba-base i]` — cierra el `if`.
- **Línea 44** `}` `[RANK 0]` — cierra el bucle de la criba base. Al salir, `marcado[k]==0` para todo primo `k ≤ raiz_n`.

```c
        for (int i = 2; i <= raiz_n; i++)
            if (marcado[i] == 0) base_count++;

        base_primes = malloc(base_count * sizeof(int));
        if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }

        int idx = 0;
        for (int i = 2; i <= raiz_n; i++)
            if (marcado[i] == 0) base_primes[idx++] = i;

        free(marcado);
    }
```

- **Línea 46** `for (int i = 2; i <= raiz_n; i++)` `[RANK 0 · LOOP conteo i]` — primer recorrido: **contar** cuántos primos base hay (para saber cuánta memoria reservar).
- **Línea 47** `if (marcado[i] == 0) base_count++;` `[RANK 0 · LOOP conteo i]` — cada casilla sin marcar es un primo: incrementa el contador.
- **Línea 49** `base_primes = malloc(base_count * sizeof(int));` `[RANK 0]` — reserva el arreglo exacto para los `base_count` primos. Aquí basta `malloc` (no hace falta `calloc`) porque lo llenaremos enseguida.
- **Línea 50** `if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }` `[RANK 0 · MPI]` — chequeo de memoria, aborta si falla.
- **Línea 52** `int idx = 0;` `[RANK 0]` — índice de escritura en `base_primes`.
- **Línea 53** `for (int i = 2; i <= raiz_n; i++)` `[RANK 0 · LOOP extracción i]` — segundo recorrido: **extraer** los primos al arreglo.
- **Línea 54** `if (marcado[i] == 0) base_primes[idx++] = i;` `[RANK 0 · LOOP extracción i]` — guarda cada primo y avanza `idx`. El `idx++` (post-incremento) usa el valor actual y luego suma 1.
- **Línea 56** `free(marcado);` `[RANK 0]` — libera el arreglo de marcas: ya extrajimos lo que servía (la lista de primos), no se necesita más.
- **Línea 57** `}` `[RANK 0]` — **cierra la rama del rank 0**. A partir de aquí vuelve el código común... pero con una asimetría: solo el rank 0 tiene `base_primes` lleno. El PASO 3 lo soluciona.

### 1.6 PASO 3 — Difundir los primos base (líneas 59-75)

Ahora el rank 0 debe **enviar** los primos base a los demás (ellos los necesitan para cribar su trozo). Como su memoria es separada, hay que mandarlos por mensaje. Se usan **dos mensajes** por destino: primero **cuántos** son (para que el receptor reserve memoria), luego **el arreglo**.

```c
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
```

- **Líneas 59-60** `/* PASO 3: ... tag 0 ... tag 1 ... */` `[TODOS]` — comentario que documenta el protocolo: los **tags** distinguen los dos mensajes (0 = el conteo, 1 = el arreglo).
- **Línea 61** `if (rank == 0) {` `[RANK 0]` — rama del emisor.
- **Línea 62** `for (int dest = 1; dest < size; dest++) {` `[RANK 0 · LOOP dest]` — recorre todos los procesos **excepto el 0** (empieza en `dest=1`). El rank 0 ya tiene los primos, no se manda a sí mismo.
- **Línea 63** `MPI_Send(&base_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);` `[RANK 0 · LOOP dest · MPI]` — **envía el conteo**. Argumentos de `MPI_Send(buffer, cantidad, tipo, destino, tag, comunicador)`: enviar 1 entero (`MPI_INT`) desde `&base_count` al proceso `dest`, con **tag 0**.
- **Línea 64** `MPI_Send(base_primes, base_count, MPI_INT, dest, 1, MPI_COMM_WORLD);` `[RANK 0 · LOOP dest · MPI]` — **envía el arreglo**: `base_count` enteros desde `base_primes` a `dest`, con **tag 1**. El receptor ya supo (por el mensaje anterior) cuántos esperar.
- **Línea 65** `}` `[RANK 0]` — cierra el bucle de destinos.
- **Línea 66** `} else {` `[RANK != 0]` — **rama de los receptores** (todos los trabajadores).
- **Líneas 67-68** `MPI_Recv(&base_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);` `[RANK != 0 · MPI]` — **recibe el conteo**. `MPI_Recv(buffer, cantidad, tipo, origen, tag, comunicador, estado)`: recibe 1 entero del proceso `0` con **tag 0** y lo guarda en `base_count`. **`MPI_STATUS_IGNORE`** indica que no nos interesa el objeto de estado (quién envió, etc.), ya lo sabemos. **Importante**: `MPI_Recv` es **bloqueante** — el proceso espera aquí hasta que el mensaje llegue.
- **Línea 70** `base_primes = malloc(base_count * sizeof(int));` `[RANK != 0]` — ahora que sabe cuántos primos vienen, **reserva** el arreglo del tamaño justo. Por esto se mandó primero el conteo: para poder reservar antes de recibir el arreglo.
- **Línea 71** `if (!base_primes) { MPI_Abort(MPI_COMM_WORLD, 1); }` `[RANK != 0 · MPI]` — chequeo de memoria.
- **Líneas 73-74** `MPI_Recv(base_primes, base_count, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);` `[RANK != 0 · MPI]` — **recibe el arreglo** de primos base: `base_count` enteros del rank 0 con **tag 1**. Tras esto, todos los procesos tienen la misma lista de primos base.
- **Línea 75** `}` `[TODOS]` — cierra el `if/else`. **Nota de diseño**: este reparto manual con un bucle de `MPI_Send` equivale conceptualmente a un `MPI_Bcast` (difusión), pero hecho a mano, punto a punto, con fines didácticos.

### 1.7 PASO 4 — Repartir el rango (líneas 77-87)

**Descomposición de dominio**: cada proceso se queda con un trozo contiguo `[low, high]` del rango `[2, N]`, calculado a partir de su `rank`. Esto se hace **sin comunicación**: cada uno calcula su trozo con una fórmula, porque todos conocen `rank`, `size` y `N`.

```c
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
```

- **Línea 77** `/* PASO 4: ... */` `[TODOS]` — comentario.
- **Línea 78** `long long low = 2 + (long long)rank * (N - 1) / size;` `[TODOS]` — calcula el **inicio** del trozo de este proceso. El rango útil es `[2, N]`, que tiene `N-1` números. Se reparte proporcionalmente: el rank 0 empieza en 2, el siguiente más arriba, etc. El cast `(long long)rank` fuerza la multiplicación `rank * (N-1)` a 64 bits (con N=10^7 y rank alto podría rozar el límite de `int`). El tipo `long long` de `low` permite manejar índices grandes.
- **Línea 79** `long long high;` `[TODOS]` — declara el **fin** del trozo (se asigna abajo según el caso).
- **Línea 80** `if (rank == size - 1)` `[TODOS]` — ¿es este el **último** proceso?
- **Línea 81** `high = N;` `[TODOS]` — sí: el último cubre hasta `N` exacto. Esto **absorbe el resto** de la división entera: si `N-1` no se reparte exacto entre `size`, los números sobrantes caen en el último proceso. Garantiza que **ningún número quede sin cubrir**.
- **Línea 82** `else` `[TODOS]` — para los demás procesos...
- **Línea 83** `high = 2 + (long long)(rank + 1) * (N - 1) / size - 1;` `[TODOS]` — su fin es **justo antes** del inicio del proceso siguiente (el `low` del `rank+1`, menos 1). Así los trozos son contiguos y no se solapan.
- **Línea 84** `long long seg_size = high - low + 1;` `[TODOS]` — tamaño del trozo (cantidad de números). El `+1` es porque el rango es **inclusivo** en ambos extremos.
- **Línea 86** `char *segmento = calloc(seg_size, sizeof(char));` `[TODOS]` — reserva el arreglo local de marcas para **este** trozo, inicializado a cero (`calloc`). De nuevo `char` para ahorrar memoria, con la convención `0=primo, 1=compuesto`. **Cada proceso tiene su propio `segmento`** en su memoria.
- **Línea 87** `if (!segmento) { MPI_Abort(MPI_COMM_WORLD, 1); }` `[TODOS · MPI]` — chequeo de memoria, aborta todo si falla.

### 1.8 PASO 5 — Cribar el segmento local (líneas 89-101)

Ahora **cada proceso, en paralelo y sin comunicarse**, tacha en su `segmento` los múltiplos de cada primo base. Esta es la parte **realmente paralela**: el trabajo pesado se reparte. Como cada uno escribe solo en su propia memoria, **no hay race conditions**.

```c
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
```

- **Línea 89** `/* PASO 5: ... */` `[TODOS]` — comentario.
- **Línea 90** `for (int i = 0; i < base_count; i++) {` `[TODOS · LOOP primos-base i]` — recorre **cada primo base**. Por cada uno se tacharán sus múltiplos dentro del trozo local.
- **Línea 91** `long long p = base_primes[i];` `[TODOS · LOOP primos-base i]` — toma el primo actual. Se copia a `long long` para que las cuentas de múltiplos (que pueden llegar hasta `N`) no desborden.
- **Línea 93** `/* primer multiplo de p >= low */` `[TODOS · LOOP primos-base i]` — comentario: hay que hallar dónde empezar a tachar dentro de **este** trozo.
- **Línea 94** `long long start = low;` `[TODOS · LOOP primos-base i]` — tentativamente, empieza en `low`.
- **Línea 95** `if (low % p != 0)` `[TODOS · LOOP primos-base i]` — si `low` **no** es múltiplo de `p`...
- **Línea 96** `start = low + (p - low % p);` `[TODOS · LOOP primos-base i]` — ...avanza hasta el **primer múltiplo de `p` que sea ≥ low**. La cuenta `p - (low % p)` es cuánto falta desde `low` hasta el siguiente múltiplo. (Si `low` ya era múltiplo, `start` se queda en `low`.) Este cálculo es lo que hace que la criba **segmentada** funcione: cada proceso empieza a tachar justo donde le toca, sin recorrer desde 0.
- **Línea 97** `if (start == p) start += p;` `[TODOS · LOOP primos-base i]` — **guarda importante**: si el primer múltiplo resultó ser `p` mismo (puede pasar en el rank 0, cuyo `low=2` incluye a los primos base), lo saltamos sumándole `p`. **Por qué**: `p` es primo, no debe marcarse como compuesto; solo queremos tachar `2p, 3p, ...`.
- **Línea 99** `for (long long j = start; j <= high; j += p)` `[TODOS · LOOP primos-base i · LOOP múltiplos j]` — recorre los múltiplos de `p` dentro del trozo, de `p` en `p`.
- **Línea 100** `segmento[j - low] = 1;` `[TODOS · LOOP primos-base i · LOOP múltiplos j]` — marca el múltiplo como compuesto. **El índice es `j - low`**, no `j`: el arreglo `segmento` es local y arranca en 0, así que el número `j` ocupa la posición `j - low`. Esta traslación de índices es típica de la criba segmentada.
- **Línea 101** `}` `[TODOS]` — cierra el bucle de primos base. Al salir, `segmento` tiene marcado todos los compuestos del trozo; los ceros restantes son **primos**.

### 1.9 PASO 6 — Contar primos locales (líneas 103-106)

```c
    /* PASO 6: contar los primos del segmento (celdas con valor 0) */
    long long conteo_local = 0;
    for (long long i = 0; i < seg_size; i++)
        if (segmento[i] == 0) conteo_local++;
```

- **Línea 103** `/* PASO 6: ... */` `[TODOS]` — comentario.
- **Línea 104** `long long conteo_local = 0;` `[TODOS]` — cuántos primos halló **este** proceso en su trozo. `long long` por si el conteo es grande.
- **Línea 105** `for (long long i = 0; i < seg_size; i++)` `[TODOS · LOOP conteo-local i]` — recorre todo el segmento local.
- **Línea 106** `if (segmento[i] == 0) conteo_local++;` `[TODOS · LOOP conteo-local i]` — cada casilla en 0 (no marcada) es un primo: lo cuenta.

### 1.10 PASO 7 — Recolectar conteos (líneas 108-118)

Cada proceso conoce su `conteo_local`, pero el rank 0 necesita **todos** para sumar el total. Se recolectan punto a punto con **tag 2**.

```c
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
```

- **Línea 108** `/* PASO 7: ... tag 2 */` `[TODOS]` — comentario.
- **Línea 109** `long long *conteos = NULL;` `[TODOS]` — puntero al arreglo de conteos. Solo el rank 0 lo usará; en los demás se queda en `NULL`.
- **Línea 110** `if (rank == 0) {` `[RANK 0]` — rama del recolector.
- **Línea 111** `conteos = malloc(size * sizeof(long long));` `[RANK 0]` — reserva un hueco por proceso (`size` enteros largos).
- **Línea 112** `conteos[0] = conteo_local;` `[RANK 0]` — el rank 0 coloca **su propio** conteo en la posición 0 (no se manda mensajes a sí mismo).
- **Línea 113** `for (int src = 1; src < size; src++)` `[RANK 0 · LOOP src]` — recorre los demás procesos para recibir su conteo.
- **Líneas 114-115** `MPI_Recv(&conteos[src], 1, MPI_LONG_LONG, src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);` `[RANK 0 · LOOP src · MPI]` — recibe 1 `long long` del proceso `src` con **tag 2** y lo guarda en `conteos[src]`. **`MPI_LONG_LONG`** es el tipo MPI que corresponde al `long long` de C.
- **Línea 116** `} else {` `[RANK != 0]` — rama de los emisores.
- **Línea 117** `MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 2, MPI_COMM_WORLD);` `[RANK != 0 · MPI]` — cada trabajador **envía** su conteo al rank 0 con **tag 2**.
- **Línea 118** `}` `[TODOS]` — cierra el `if/else`.

```c
    double t_fin = MPI_Wtime();
    double tiempo = t_fin - t_inicio;
```

- **Línea 120** `double t_fin = MPI_Wtime();` `[TODOS]` — marca el tiempo **final**. Se toma aquí (tras el cómputo y la recolección de conteos) porque lo que sigue (PASO 8 y 9) es **E/S** (imprimir y escribir archivos), que no es parte del algoritmo paralelo que queremos medir.
- **Línea 121** `double tiempo = t_fin - t_inicio;` `[TODOS]` — duración total del trabajo. El valor relevante es el del rank 0, que es quien lo reporta.

### 1.11 PASO 8 — Resumen y resultados.dat (líneas 120-167)

Solo el rank 0 imprime el resumen en pantalla y escribe `resultados.dat`.

```c
    /* PASO 8: rank 0 muestra el resumen y escribe resultados.dat */
    if (rank == 0) {
        long long total = 0;
        for (int i = 0; i < size; i++) total += conteos[i];

        printf("\n=== Criba de Eratostenes Distribuida (N=%lld) ===\n", (long long)N);
        printf("Procesos MPI        : %d\n", size);
        printf("Primos encontrados  : %lld\n", total);
        printf("Tiempo de ejecucion : %.6f segundos\n\n", tiempo);
```

- **Línea 123** `/* PASO 8: ... */` `[TODOS]` — comentario.
- **Línea 124** `if (rank == 0) {` `[RANK 0]` — abre la rama del maestro para reportar.
- **Línea 125** `long long total = 0;` `[RANK 0]` — acumulador del total de primos.
- **Línea 126** `for (int i = 0; i < size; i++) total += conteos[i];` `[RANK 0 · LOOP suma i]` — suma los conteos de todos los procesos → total de primos hasta N.
- **Línea 128** `printf("\n=== Criba de Eratostenes Distribuida (N=%lld) ===\n", (long long)N);` `[RANK 0]` — imprime el encabezado. `%lld` es el formato para `long long`; el cast `(long long)N` asegura que `N` (que es un `int` por el `#define`) se imprima con el formato correcto, evitando un desajuste tipo/formato.
- **Línea 129** `printf("Procesos MPI : %d\n", size);` `[RANK 0]` — cuántos procesos participaron.
- **Línea 130** `printf("Primos encontrados : %lld\n", total);` `[RANK 0]` — el total calculado (para N=10^7 deben ser 664.579 primos).
- **Línea 131** `printf("Tiempo de ejecucion : %.6f segundos\n\n", tiempo);` `[RANK 0]` — el tiempo medido, con 6 decimales.

```c
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
```

- **Líneas 133-134** `printf("%-8s %-12s %-12s %-10s\n", "Proceso", "Inicio", "Fin", "Primos");` `[RANK 0]` — **cabecera de la tabla**. Los formatos `%-8s`, `%-12s`, etc. alinean cadenas a la izquierda en anchos fijos, creando columnas.
- **Líneas 135-136** `printf(... "-------", ...)` `[RANK 0]` — línea separadora de guiones bajo la cabecera.
- **Línea 137** `for (int i = 0; i < size; i++) {` `[RANK 0 · LOOP tabla i]` — una fila por proceso.
- **Línea 138** `long long lo = 2 + (long long)i * (N - 1) / size;` `[RANK 0 · LOOP tabla i]` — **recalcula** el inicio del trozo del proceso `i`, con la **misma fórmula** de la línea 78. Lo recalcula (en vez de habérselo pedido por mensaje) porque la fórmula es determinista: el rank 0 puede reproducir el reparto de cualquier proceso conociendo solo `i`, `N` y `size`. Ahorra comunicación.
- **Línea 139** `long long hi;` `[RANK 0 · LOOP tabla i]` — fin del trozo de `i`.
- **Líneas 140-143** `if (i == size - 1) hi = N; else hi = ... - 1;` `[RANK 0 · LOOP tabla i]` — mismo cálculo de `high` de las líneas 80-83 (último proceso llega a N; los demás, justo antes del siguiente).
- **Línea 144** `printf("%-8d %-12lld %-12lld %-10lld\n", i, lo, hi, conteos[i]);` `[RANK 0 · LOOP tabla i]` — imprime la fila: proceso, inicio, fin, primos. `%-12lld` es `long long` alineado a la izquierda en ancho 12.
- **Línea 145** `}` `[RANK 0]` — cierra el bucle de la tabla.
- **Línea 147** `printf("\nTIEMPO: %.6f\n", tiempo);` `[RANK 0]` — imprime una línea con el prefijo **`TIEMPO:`**. Esto es **deliberado**: el `Makefile` (target `benchmark`, líneas 50-52) busca esta línea con `grep "^TIEMPO:"` para extraer el tiempo y calcular el speedup. Es un "contrato" entre el programa y el script.

```c
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
```

- **Línea 149** `FILE *fp = fopen("resultados/resultados.dat", "w");` `[RANK 0]` — abre el archivo de salida en modo escritura (`"w"`, trunca si existía). La carpeta `resultados/` debe existir. Este archivo lo lee Gnuplot para la gráfica de distribución.
- **Línea 150** `if (fp) {` `[RANK 0]` — solo escribe si la apertura tuvo éxito (manejo de error silencioso: si falla, simplemente no escribe).
- **Línea 151** `fprintf(fp, "# Criba ... N=%lld ... Procesos=%d\n", (long long)N, size);` `[RANK 0]` — primera línea de comentario (`#`) del archivo, con metadatos. Gnuplot ignora las líneas con `#`.
- **Línea 152** `fprintf(fp, "# Proceso  Inicio  Fin  Primos\n");` `[RANK 0]` — comentario que documenta las columnas.
- **Línea 153** `for (int i = 0; i < size; i++) {` `[RANK 0 · LOOP archivo i]` — una fila de datos por proceso.
- **Líneas 154-159** — **recalculan** `lo` y `hi` exactamente igual que en la tabla de pantalla (líneas 138-143). Es código repetido a propósito para escribir las mismas columnas al archivo.
- **Línea 160** `fprintf(fp, "%d %lld %lld %lld\n", i, lo, hi, conteos[i]);` `[RANK 0 · LOOP archivo i]` — escribe la fila en formato simple separado por espacios (lo que Gnuplot espera): proceso, inicio, fin, primos.
- **Línea 162** `fclose(fp);` `[RANK 0]` — cierra el archivo (vacía el buffer a disco).
- **Línea 163** `printf("Resultados guardados en resultados.dat\n");` `[RANK 0]` — confirma por pantalla.
- **Línea 164** `}` `[RANK 0]` — cierra el `if (fp)`.
- **Línea 166** `free(conteos);` `[RANK 0]` — libera el arreglo de conteos (ya no se usa).
- **Línea 167** `}` `[RANK 0]` — cierra la rama del PASO 8.

### 1.12 PASO 9 — Recolectar y escribir primos.dat (líneas 169-213)

Además del conteo, se quiere la **lista completa de primos** en un archivo, **ordenada**. Como cada proceso tiene solo su trozo, el rank 0 los recolecta **en orden de rank** (tags 3 y 4) y los va escribiendo, de modo que el archivo queda ordenado de menor a mayor.

```c
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
```

- **Líneas 169-172** `/* PASO 9: ... OJO: para N grande este archivo es enorme ... */` `[TODOS]` — comentario que explica el protocolo (tags 3 y 4) y **advierte** del coste: la lista de primos puede tener decenas de millones de líneas para N grande. Es un aviso de rendimiento/espacio en disco.
- **Línea 173** `if (rank == 0) {` `[RANK 0]` — rama del que escribe el archivo.
- **Línea 174** `FILE *fp_p = fopen("resultados/primos.dat", "w");` `[RANK 0]` — abre `primos.dat` para escritura.
- **Línea 175** `if (fp_p) {` `[RANK 0]` — solo procede si abrió bien.
- **Línea 176** `fprintf(fp_p, "# Lista de primos hasta N=%lld\n", (long long)N);` `[RANK 0]` — comentario de cabecera del archivo.
- **Línea 178** `for (long long i = 0; i < seg_size; i++)` `[RANK 0 · LOOP primos-propios i]` — recorre el segmento **del propio rank 0** primero (sus primos son los más pequeños, deben ir al principio para que el archivo quede ordenado).
- **Línea 179** `if (segmento[i] == 0) fprintf(fp_p, "%lld\n", low + i);` `[RANK 0 · LOOP primos-propios i]` — por cada casilla no marcada, escribe el **número real** `low + i` (reconstruido a partir del índice local + el inicio del trozo). Un primo por línea.

```c
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
    }
```

- **Línea 181** `for (int src = 1; src < size; src++) {` `[RANK 0 · LOOP src]` — recibe de los trabajadores **en orden de rank** (1, 2, 3...). Ese orden es lo que garantiza que el archivo final quede **ordenado**, porque los trozos son contiguos y crecientes.
- **Línea 182** `long long n_src;` `[RANK 0 · LOOP src]` — cuántos primos enviará este `src`.
- **Líneas 183-184** `MPI_Recv(&n_src, 1, MPI_LONG_LONG, src, 3, ...);` `[RANK 0 · LOOP src · MPI]` — **recibe el conteo** del proceso `src` con **tag 3** (para saber cuánto reservar).
- **Línea 186** `long long *buf = malloc(n_src * sizeof(long long));` `[RANK 0 · LOOP src]` — reserva un buffer temporal para los primos de este `src`.
- **Línea 187** `if (!buf) { MPI_Abort(MPI_COMM_WORLD, 1); }` `[RANK 0 · LOOP src · MPI]` — chequeo de memoria.
- **Líneas 189-190** `MPI_Recv(buf, n_src, MPI_LONG_LONG, src, 4, ...);` `[RANK 0 · LOOP src · MPI]` — **recibe el arreglo** de primos de `src` con **tag 4**.
- **Línea 192** `for (long long i = 0; i < n_src; i++)` `[RANK 0 · LOOP src · LOOP escritura i]` — recorre los primos recibidos.
- **Línea 193** `fprintf(fp_p, "%lld\n", buf[i]);` `[RANK 0 · LOOP src · LOOP escritura i]` — los escribe al archivo, uno por línea. Como ya vienen como números reales (los empaquetó el emisor en la línea 207), se escriben tal cual.
- **Línea 195** `free(buf);` `[RANK 0 · LOOP src]` — libera el buffer temporal antes de pasar al siguiente `src` (evita acumular memoria).
- **Línea 196** `}` `[RANK 0]` — cierra el bucle de fuentes.
- **Línea 198** `fclose(fp_p);` `[RANK 0]` — cierra `primos.dat`.
- **Línea 199** `printf("Lista de primos guardada en resultados/primos.dat\n");` `[RANK 0]` — confirma.
- **Línea 200** `}` `[RANK 0]` — cierra el `if (fp_p)`.
- **Línea 201** `}` `[RANK 0]` — cierra la rama del rank 0 del PASO 9.

```c
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
```

- **Línea 201** `} else {` `[RANK != 0]` — la **rama de los trabajadores** del PASO 9 (empaquetar y enviar sus primos). (La llave `}` de esta línea cierra el `if (rank==0)` de la línea 173; el `else` abre la rama de los demás.)
- **Línea 202** `long long *buf = malloc(conteo_local * sizeof(long long));` `[RANK != 0]` — reserva un buffer del tamaño exacto de **sus** primos (ya conoce `conteo_local` del PASO 6).
- **Línea 203** `if (!buf) { MPI_Abort(MPI_COMM_WORLD, 1); }` `[RANK != 0 · MPI]` — chequeo de memoria.
- **Línea 205** `long long k = 0;` `[RANK != 0]` — índice de escritura en `buf`.
- **Línea 206** `for (long long i = 0; i < seg_size; i++)` `[RANK != 0 · LOOP empaquetado i]` — recorre su segmento.
- **Línea 207** `if (segmento[i] == 0) buf[k++] = low + i;` `[RANK != 0 · LOOP empaquetado i]` — guarda en `buf` cada primo como **número real** (`low + i`). Esto "comprime" el segmento (que tiene huecos de compuestos) en una lista densa solo de primos, lista para enviar.
- **Línea 209** `MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 3, MPI_COMM_WORLD);` `[RANK != 0 · MPI]` — **envía cuántos** primos tiene al rank 0 con **tag 3** (coincide con el `Recv` de la línea 183).
- **Línea 210** `MPI_Send(buf, conteo_local, MPI_LONG_LONG, 0, 4, MPI_COMM_WORLD);` `[RANK != 0 · MPI]` — **envía el arreglo** de primos con **tag 4** (coincide con el `Recv` de la línea 189).
- **Línea 212** `free(buf);` `[RANK != 0]` — libera el buffer tras enviarlo.
- **Línea 213** `}` `[TODOS]` — cierra el `if/else` del PASO 9.

### 1.13 Cierre (líneas 215-219)

```c
    free(base_primes);
    free(segmento);
    MPI_Finalize();
    return 0;
}
```

- **Línea 215** `free(base_primes);` `[TODOS]` — libera el arreglo de primos base. **Todos** los procesos lo tienen reservado (el rank 0 lo calculó, los demás lo recibieron), así que todos lo liberan.
- **Línea 216** `free(segmento);` `[TODOS]` — libera el segmento local. También lo tienen todos.
- **Línea 217** `MPI_Finalize();` `[TODOS · MPI]` — **cierra el entorno MPI**. Debe ser la última llamada MPI; tras ella ya no se puede comunicar. Es obligatoria para una terminación limpia (libera recursos de MPI).
- **Línea 218** `return 0;` `[TODOS]` — termina con código 0 (éxito). Cada proceso retorna por su cuenta.
- **Línea 219** `}` `[GLOBAL]` — cierra `main`.

---

## 2. `Makefile`

Automatiza compilación y ejecución. En un Makefile las recetas se indentan con **TAB**, `#` inicia comentario, `@` al inicio de un comando lo ejecuta en silencio (sin imprimirlo), y `\` continúa una línea.

```makefile
# =============================================================
# Makefile — Criba de Eratóstenes Distribuida con MPI
#
# SETUP DEL ENTORNO
# -----------------
# OBLIGATORIO (compilar y ejecutar el programa):
#   Ubuntu/Debian : sudo apt install libmpich-dev mpich
#   Fedora/RHEL   : sudo dnf install mpich-devel          <-- ya instalado en este sistema
#   macOS         : brew install open-mpi
#   Windows WSL2  : sudo apt install libopenmpi-dev openmpi-bin
#
# OPCIONAL (gráficas):
#   Ubuntu/Debian : sudo apt install gnuplot
#   Fedora/RHEL   : sudo dnf install gnuplot
#   macOS         : brew install gnuplot
#
# OPCIONAL (informe PDF):
#   Ubuntu/Debian : sudo apt install texlive-full
#   Fedora/RHEL   : sudo dnf install texlive-scheme-full  <-- ya instalado en este sistema
#   macOS         : brew install --cask mactex
# =============================================================
```

- **Líneas 1-21** `[comentarios]` — bloque de documentación. No afecta a la compilación; lista qué paquetes instalar en cada sistema operativo para (a) **obligatorio**: MPI (`mpich`/`openmpi`), necesario para compilar y ejecutar; (b) **opcional**: `gnuplot` (gráficas) y `texlive` (informe PDF). Las anotaciones `<-- ya instalado en este sistema` indican que en Fedora/RHEL ya está listo. Útil para que cualquiera reproduzca el entorno.

```makefile
CC     = mpicc
CFLAGS = -O2 -Wall
TARGET = criba
```

- **Línea 23** `CC = mpicc` — el compilador es **`mpicc`** (el wrapper de MPI), no `gcc`. Esta es la diferencia clave con un proyecto no-MPI.
- **Línea 24** `CFLAGS = -O2 -Wall` — flags: optimización nivel 2 y todos los warnings.
- **Línea 25** `TARGET = criba` — nombre del ejecutable.

```makefile
# ── Compilar ──────────────────────────────────────────────────
all: $(TARGET)

$(TARGET): src/criba.c
	$(CC) $(CFLAGS) -o $(TARGET) src/criba.c -lm
```

- **Línea 27** `# ── Compilar ── ...` — comentario separador.
- **Línea 28** `all: $(TARGET)` — target por defecto (`make` sin argumentos): depende de `criba`. Sin comandos propios.
- **Línea 30** `$(TARGET): src/criba.c` — regla para construir `criba`; depende del fuente (recompila solo si cambió).
- **Línea 31** `	$(CC) $(CFLAGS) -o $(TARGET) src/criba.c -lm` — comando (con TAB). Expande a `mpicc -O2 -Wall -o criba src/criba.c -lm`. El `-lm` va **al final** (tras el fuente) por la resolución de símbolos izq→der del enlazador.

```makefile
# ── Ejecución rápida (4 procesos) ─────────────────────────────
run: $(TARGET)
	mpirun -np 4 ./$(TARGET)
```

- **Línea 33** `# ── Ejecución rápida (4 procesos) ── ...` — comentario.
- **Línea 34** `run: $(TARGET)` — target `run`, depende de `criba` (compila si hace falta).
- **Línea 35** `	mpirun -np 4 ./$(TARGET)` — lanza **4 procesos** con `mpirun -np 4 ./criba`. Es la ejecución estándar.

```makefile
# ── Benchmark: 1, 2 y 4 procesos → tabla de speedup ──────────
# Crea speedup.dat (leído por graficar.gp) y muestra la tabla.
benchmark: $(TARGET)
	@echo ""
	@echo ">>> Corriendo con 1 proceso..."
	@mpirun -np 1 ./$(TARGET) | tee /tmp/_criba_t1.txt
	@echo ""
	@echo ">>> Corriendo con 2 procesos..."
	@mpirun -np 2 ./$(TARGET) | tee /tmp/_criba_t2.txt
	@echo ""
	@echo ">>> Corriendo con 4 procesos..."
	@mpirun -np 4 ./$(TARGET) | tee /tmp/_criba_t4.txt
	@echo ""
```

- **Líneas 37-38** `# ── Benchmark ... # Crea speedup.dat ...` — comentarios que explican el target: corre con 1, 2 y 4 procesos y produce la tabla de speedup y el archivo `speedup.dat`.
- **Línea 39** `benchmark: $(TARGET)` — target `benchmark`, depende de `criba`.
- **Línea 40** `	@echo ""` — imprime una línea en blanco. El **`@`** evita que Make muestre el comando antes de ejecutarlo (salida más limpia).
- **Línea 41** `	@echo ">>> Corriendo con 1 proceso..."` — mensaje de progreso.
- **Línea 42** `	@mpirun -np 1 ./$(TARGET) | tee /tmp/_criba_t1.txt` — corre con **1 proceso** y usa **`tee`**: muestra la salida en pantalla **y** a la vez la guarda en `/tmp/_criba_t1.txt`. Se guarda para luego extraer el tiempo.
- **Líneas 43-45** — repiten el patrón con **2 procesos**, guardando en `_criba_t2.txt`.
- **Líneas 46-48** — repiten con **4 procesos**, guardando en `_criba_t4.txt`.
- **Línea 49** `	@echo ""` — línea en blanco separadora.

```makefile
	@T1=$$(grep "^TIEMPO:" /tmp/_criba_t1.txt | awk '{print $$2}'); \
	 T2=$$(grep "^TIEMPO:" /tmp/_criba_t2.txt | awk '{print $$2}'); \
	 T4=$$(grep "^TIEMPO:" /tmp/_criba_t4.txt | awk '{print $$2}'); \
	 awk -v t1="$$T1" -v t2="$$T2" -v t4="$$T4" 'BEGIN { \
	   print ""; \
	   print "=== TABLA DE SPEEDUP ==="; \
	   print "+----------+-------------+---------+"; \
	   print "| Procesos |    Tiempo   | Speedup |"; \
	   print "+----------+-------------+---------+"; \
	   printf "| %8d | %9.4f s  | %7.3f |\n", 1, t1, 1.0;      \
	   printf "| %8d | %9.4f s  | %7.3f |\n", 2, t2, t1/t2;    \
	   printf "| %8d | %9.4f s  | %7.3f |\n", 4, t4, t1/t4;    \
	   print "+----------+-------------+---------+"; \
	   print ""; \
	   printf "# nprocs  tiempo    speedup\n"         > "resultados/speedup.dat"; \
	   printf "1 %f 1.000\n",       t1               >> "resultados/speedup.dat"; \
	   printf "2 %f %f\n",    t2, t1/t2              >> "resultados/speedup.dat"; \
	   printf "4 %f %f\n",    t4, t1/t4              >> "resultados/speedup.dat"; \
	   print "speedup.dat actualizado."; \
	 }'
```

Este bloque es **un solo comando de shell** repartido en muchas líneas con `\` al final de cada una. **Punto clave de Make**: un `$` literal para el shell se escribe **`$$`** en el Makefile (porque Make interpreta `$` como inicio de sus propias variables). Por eso ves `$$(...)` y `$$2`.

- **Líneas 50-52** `@T1=$$(grep "^TIEMPO:" .../_criba_t1.txt | awk '{print $$2}'); \` — define variables de shell `T1`, `T2`, `T4`. Cada una: busca con `grep` la línea que empieza por `TIEMPO:` (la que imprime el programa en la línea 147 de `criba.c`) en el archivo correspondiente, y con `awk '{print $2}'` extrae el **segundo campo** (el número del tiempo). El `\` encadena las tres asignaciones y el `awk` siguiente como un único comando, de modo que las variables `T1/T2/T4` sigan vivas (cada línea de receta de Make corre en un shell nuevo; por eso todo va junto).
- **Línea 53** `awk -v t1="$$T1" -v t2="$$T2" -v t4="$$T4" 'BEGIN { \` — invoca `awk` pasándole los tres tiempos como variables (`-v`). El bloque `BEGIN { ... }` se ejecuta una vez (sin necesidad de archivo de entrada), solo para **formatear e imprimir**.
- **Líneas 54-58** `print "..."` — imprimen el encabezado de la tabla de speedup (bordes ASCII y títulos de columna).
- **Línea 59** `printf "| %8d | %9.4f s  | %7.3f |\n", 1, t1, 1.0;` — fila de **1 proceso**: speedup es 1.0 por definición (es la base).
- **Línea 60** `printf "... \n", 2, t2, t1/t2;` — fila de **2 procesos**: speedup = `t1/t2` (tiempo base entre tiempo con 2).
- **Línea 61** `printf "... \n", 4, t4, t1/t4;` — fila de **4 procesos**: speedup = `t1/t4`.
- **Línea 62** `print "+----------+-------------+---------+";` — borde inferior de la tabla.
- **Línea 63** `print "";` — línea en blanco.
- **Línea 64** `printf "# nprocs  tiempo    speedup\n" > "resultados/speedup.dat";` — empieza a **escribir el archivo** `speedup.dat`: el `>` **crea/sobrescribe** el archivo con la línea de cabecera. Este archivo lo leerá Gnuplot para la gráfica de speedup.
- **Línea 65** `printf "1 %f 1.000\n", t1 >> "resultados/speedup.dat";` — **añade** (`>>`, append) la fila de 1 proceso: nprocs, tiempo, speedup=1.000.
- **Línea 66** `printf "2 %f %f\n", t2, t1/t2 >> "resultados/speedup.dat";` — añade la fila de 2 procesos con su speedup calculado.
- **Línea 67** `printf "4 %f %f\n", t4, t1/t4 >> "resultados/speedup.dat";` — añade la fila de 4 procesos.
- **Línea 68** `print "speedup.dat actualizado.";` — mensaje de confirmación.
- **Línea 69** `}'` — cierra el bloque `BEGIN` y la invocación de `awk`. Fin del comando del target `benchmark`.

```makefile
# ── Generar gráficas PNG (requiere gnuplot) ───────────────────
graficas:
	gnuplot graficar.gp
	@echo "Generadas: distribucion.png  speedup.png"
```

- **Línea 71** `# ── Generar gráficas PNG ... ──` — comentario.
- **Línea 72** `graficas:` — target sin dependencias.
- **Línea 73** `	gnuplot graficar.gp` — ejecuta Gnuplot con el script, produciendo los PNG. (Requiere que existan `resultados.dat` y `speedup.dat`.)
- **Línea 74** `	@echo "Generadas: distribucion.png  speedup.png"` — confirma en silencio (`@`).

```makefile
# ── Compilar informe IEEE ─────────────────────────────────────
informe:
	cd informe && pdflatex informe_ieee.tex && pdflatex informe_ieee.tex
```

- **Línea 76** `# ── Compilar informe IEEE ──` — comentario.
- **Línea 77** `informe:` — target del informe.
- **Línea 78** `	cd informe && pdflatex informe_ieee.tex && pdflatex informe_ieee.tex` — entra a `informe/` y compila el LaTeX **dos veces** (la primera genera referencias/citas, la segunda las resuelve). `&&` encadena: cada comando solo corre si el anterior tuvo éxito.

```makefile
# ── Compilar referencias ──────────────────────────────────────
referencias:
	cd informe && pdflatex referencias.tex
```

- **Línea 80** `# ── Compilar referencias ──` — comentario.
- **Línea 81** `referencias:` — target.
- **Línea 82** `	cd informe && pdflatex referencias.tex` — compila `referencias.tex` (una pasada).

```makefile
# ── Limpiar archivos generados ────────────────────────────────
clean:
	rm -f $(TARGET)
	rm -f resultados/resultados.dat resultados/speedup.dat
	rm -f resultados/distribucion.png resultados/speedup.png
	rm -f informe/*.aux informe/*.log informe/*.out

.PHONY: all run benchmark graficas informe referencias clean
```

- **Línea 84** `# ── Limpiar archivos generados ──` — comentario.
- **Línea 85** `clean:` — target de limpieza.
- **Línea 86** `	rm -f $(TARGET)` — borra el ejecutable. `-f` (force) no se queja si no existe.
- **Línea 87** `	rm -f resultados/resultados.dat resultados/speedup.dat` — borra los archivos de datos generados.
- **Línea 88** `	rm -f resultados/distribucion.png resultados/speedup.png` — borra las gráficas.
- **Línea 89** `	rm -f informe/*.aux informe/*.log informe/*.out` — borra los auxiliares de LaTeX.
- **Línea 91** `.PHONY: all run benchmark graficas informe referencias clean` — declara estos targets como **"phony"** (acciones, no archivos). Evita que, si existiera un archivo con el nombre de un target (p. ej. `clean`), Make lo considere "ya hecho" y no ejecute la receta. Garantiza que siempre se ejecuten.

---

## 3. `graficar.gp`

Script de Gnuplot que produce dos PNG a partir de `resultados.dat` y `speedup.dat`. `#` inicia comentario; `set` configura opciones.

```gnuplot
# =============================================================
# graficar.gp — Gnuplot: gráficas de la Criba de Eratóstenes
#
# Uso:
#   gnuplot graficar.gp
#
# Requiere:
#   resultados.dat  (generado por: mpirun -np 4 ./criba)
#   speedup.dat     (generado por: make benchmark)
#
# Genera:
#   distribucion.png  — barras: primos encontrados por proceso
#   speedup.png       — líneas: speedup real vs. ideal
# =============================================================
```

- **Líneas 1-14** `[comentarios]` — cabecera que documenta el uso (`gnuplot graficar.gp`), los **prerrequisitos** (los dos `.dat` y cómo se generan) y las **salidas** (los dos PNG y qué muestra cada uno).

### Gráfica 1: Distribución de primos por proceso (líneas 16-48)

```gnuplot
# ─────────────────────────────────────────────────────────────
# GRÁFICA 1: Distribución de primos por proceso
# ─────────────────────────────────────────────────────────────
# Leemos resultados.dat con formato:
#   proceso  inicio  fin  cantidad_primos
# Columna 1 = rank del proceso (eje X)
# Columna 4 = primos encontrados (eje Y)
# ─────────────────────────────────────────────────────────────
```

- **Líneas 16-23** `[comentarios]` — describen la primera gráfica (un histograma/barras: cuántos primos halló cada proceso) y el formato del archivo `resultados.dat`: columna 1 = rank (eje X), columna 4 = primos (eje Y). Sirve para visualizar si la carga quedó **balanceada** entre procesos.

```gnuplot
set terminal pngcairo size 900,600 enhanced font 'Helvetica,13'
set output 'resultados/distribucion.png'

set title 'Distribución de Primos por Proceso MPI (N = 10,000,000)' font ',15'
set xlabel 'Proceso MPI (rank)'
set ylabel 'Cantidad de Primos'
```

- **Línea 25** `set terminal pngcairo size 900,600 enhanced font 'Helvetica,13'` — formato de salida: PNG con motor `pngcairo`, 900×600 px, texto `enhanced`, fuente Helvetica 13.
- **Línea 26** `set output 'resultados/distribucion.png'` — archivo de salida de esta gráfica.
- **Línea 28** `set title 'Distribución de Primos por Proceso MPI (N = 10,000,000)' font ',15'` — título. `font ',15'` deja la fuente por defecto pero a tamaño 15.
- **Línea 29** `set xlabel 'Proceso MPI (rank)'` — etiqueta del eje X.
- **Línea 30** `set ylabel 'Cantidad de Primos'` — etiqueta del eje Y.

```gnuplot
set style data histograms
set style histogram clustered gap 1
set style fill solid 0.80 border -1
set boxwidth 0.6
```

- **Línea 32** `set style data histograms` — dibuja los datos como **histograma** (barras).
- **Línea 33** `set style histogram clustered gap 1` — estilo de histograma "agrupado" con un hueco (`gap 1`) entre grupos de barras.
- **Línea 34** `set style fill solid 0.80 border -1` — rellena las barras con color sólido al 80% de opacidad y un borde (`border -1` = color por defecto).
- **Línea 35** `set boxwidth 0.6` — ancho de cada barra (0.6, deja espacio entre ellas).

```gnuplot
set yrange [0:*]
set xrange [-0.5:*]
set grid y
set key off

set format y '%g'
```

- **Línea 37** `set yrange [0:*]` — eje Y desde 0 hasta automático (`*`).
- **Línea 38** `set xrange [-0.5:*]` — eje X desde -0.5 (un margen para que la primera barra no quede pegada al borde) hasta automático.
- **Línea 39** `set grid y` — cuadrícula solo en el eje Y (líneas horizontales), para leer mejor las alturas.
- **Línea 40** `set key off` — **oculta la leyenda** (con una sola serie de barras no aporta).
- **Línea 42** `set format y '%g'` — formato de las etiquetas del eje Y: `%g` muestra los números de forma compacta (sin ceros innecesarios).

```gnuplot
# Comentario de datos: column(4) = primos, xticlabels(1) = rank
plot 'resultados/resultados.dat' \
     using 4:xticlabels(1) \
     title 'Primos por proceso' \
     linecolor rgb '#4472C4'
```

- **Línea 44** `# Comentario de datos: ...` — comentario que recuerda qué columnas se usan.
- **Línea 45** `plot 'resultados/resultados.dat' \` — empieza el `plot`, leyendo `resultados.dat`. El `\` continúa en la siguiente línea.
- **Línea 46** `using 4:xticlabels(1)` — usa la **columna 4** (cantidad de primos) como altura de la barra, y la **columna 1** (rank) como **etiqueta del eje X** de cada barra (`xticlabels`).
- **Línea 47** `title 'Primos por proceso'` — nombre de la serie (aunque la leyenda está apagada).
- **Línea 48** `linecolor rgb '#4472C4'` — color de las barras (azul).

### Gráfica 2: Speedup real vs ideal (líneas 51-87)

```gnuplot
# ─────────────────────────────────────────────────────────────
# GRÁFICA 2: Speedup real vs. ideal
# ─────────────────────────────────────────────────────────────
# Leemos speedup.dat con formato:
#   nprocs  tiempo  speedup
# Columna 1 = número de procesos (eje X)
# Columna 3 = speedup medido (eje Y)
# La línea ideal es simplemente y = x (speedup perfecto = n)
# ─────────────────────────────────────────────────────────────
```

- **Líneas 51-59** `[comentarios]` — describen la segunda gráfica (speedup medido vs el ideal) y el formato de `speedup.dat`: columna 1 = nº de procesos (X), columna 3 = speedup (Y). El **speedup ideal** es `y = x` (con `p` procesos, ir `p` veces más rápido).

```gnuplot
set terminal pngcairo size 900,600 enhanced font 'Helvetica,13'
set output 'resultados/speedup.png'

set title 'Speedup — Criba de Eratóstenes con MPI' font ',15'
set xlabel 'Número de Procesos'
set ylabel 'Speedup'
```

- **Línea 61** `set terminal pngcairo size 900,600 enhanced font 'Helvetica,13'` — vuelve a fijar el terminal (mismo formato). Repetirlo es inofensivo y deja claro el contexto de esta segunda gráfica.
- **Línea 62** `set output 'resultados/speedup.png'` — **cambia la salida** al segundo PNG. A partir de aquí lo que se grafique va a `speedup.png`.
- **Línea 64** `set title 'Speedup — Criba de Eratóstenes con MPI' font ',15'` — nuevo título (sobrescribe el de la gráfica 1).
- **Línea 65** `set xlabel 'Número de Procesos'` — etiqueta X.
- **Línea 66** `set ylabel 'Speedup'` — etiqueta Y.

```gnuplot
set style data linespoints
unset style histogram

set xrange [0.5:4.5]
set yrange [0:5]
set xtics (1, 2, 4)
set grid
```

- **Línea 68** `set style data linespoints` — cambia el estilo de dibujo a **líneas con puntos** (apropiado para una curva, ya no barras).
- **Línea 69** `unset style histogram` — **desactiva** la configuración de histograma que quedó de la gráfica 1. Limpieza necesaria para que no interfiera.
- **Línea 71** `set xrange [0.5:4.5]` — eje X de 0.5 a 4.5, con margen alrededor de los valores 1, 2, 4.
- **Línea 72** `set yrange [0:5]` — eje Y de 0 a 5 (el speedup ideal con 4 procesos es 4; deja algo de margen).
- **Línea 73** `set xtics (1, 2, 4)` — marcas del eje X **solo** en 1, 2 y 4 (los números de procesos realmente probados).
- **Línea 74** `set grid` — cuadrícula en ambos ejes.

```gnuplot
# Función de speedup ideal: S(p) = p  (escala lineal perfecta)
ideal(x) = x

plot 'resultados/speedup.dat' \
       using 1:3 \
       with linespoints \
       title 'Speedup real' \
       lw 2 pt 7 ps 1.5 lc rgb '#C0392B', \
     ideal(x) \
       with lines \
       title 'Speedup ideal (S=p)' \
       lw 1.5 dt 2 lc rgb '#7F8C8D'
```

- **Línea 76** `# Función de speedup ideal: ...` — comentario.
- **Línea 77** `ideal(x) = x` — define la función `ideal` (la recta `y = x`), que representa el **speedup perfecto** (con `p` procesos, speedup `p`). Sirve de referencia para comparar.
- **Línea 79** `plot 'resultados/speedup.dat' \` — empieza el `plot` con dos curvas separadas por coma.
- **Línea 80** `using 1:3` — primera curva (**speedup real**): columna 1 (nº procesos) en X, columna 3 (speedup medido) en Y.
- **Línea 81** `with linespoints` — dibuja línea **y** puntos.
- **Línea 82** `title 'Speedup real'` — nombre en la leyenda.
- **Línea 83** `lw 2 pt 7 ps 1.5 lc rgb '#C0392B', ` — estilo: grosor 2 (`lw`), punto tipo 7 círculo relleno (`pt`), tamaño 1.5 (`ps`), color rojo (`lc rgb`). La coma final separa de la segunda curva.
- **Línea 84** `ideal(x) \` — segunda curva: la función ideal.
- **Línea 85** `with lines` — solo línea (sin puntos).
- **Línea 86** `title 'Speedup ideal (S=p)'` — nombre en la leyenda.
- **Línea 87** `lw 1.5 dt 2 lc rgb '#7F8C8D'` — estilo: grosor 1.5, trazo discontinuo (`dt 2`, dashed), color gris. Sin coma: es la última curva. Comparar la curva roja (real) con la gris discontinua (ideal) muestra cuán lejos del speedup perfecto queda la implementación.

---

## Resumen de la comunicación MPI

Toda la coordinación del programa se hace con cinco intercambios punto a punto, identificados por **tag**:

| Tag | Mensaje | Dirección | Cuándo (paso) |
|-----|---------|-----------|---------------|
| 0 | `base_count` (cuántos primos base) | rank 0 → cada trabajador | PASO 3 |
| 1 | `base_primes` (el arreglo de primos base) | rank 0 → cada trabajador | PASO 3 |
| 2 | `conteo_local` (primos hallados por cada uno) | cada trabajador → rank 0 | PASO 7 |
| 3 | `conteo_local` (cuántos primos enviará) | cada trabajador → rank 0 | PASO 9 |
| 4 | `buf` (la lista de primos del trabajador) | cada trabajador → rank 0 | PASO 9 |

**Ideas clave del proyecto**:
- **Descomposición de dominio**: el rango `[2, N]` se parte en trozos contiguos, uno por proceso (PASO 4). El cómputo pesado (PASO 5) ocurre **en paralelo y sin comunicación**.
- **Patrón maestro-trabajador**: el rank 0 prepara los datos comunes (primos base), los reparte, y al final **recolecta** conteos y primos.
- **Sin race conditions**: como cada proceso tiene **memoria separada**, nunca dos escriben en la misma celda. El precio que se paga en MPI es el de **comunicar** (los `Send`/`Recv`), no el de sincronizar accesos a memoria compartida.
- **`long long` por todas partes**: con N grande, los índices y conteos superan el rango de `int`; usar 64 bits evita overflow.
- **Determinismo del reparto**: como `low`/`high` se calculan con una fórmula a partir de `rank`, el rank 0 puede reconstruir el trozo de cualquier proceso sin preguntar (lo aprovecha en el PASO 8).
```
