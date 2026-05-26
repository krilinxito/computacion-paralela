# Guía de Estudio: K-Means Paralelo con OpenMP
### Explicado desde cero, para alguien que nunca lo ha visto

**Proyecto:** K-Means Clustering paralelizado con OpenMP — Dataset Iris  
**Nivel:** Sin conocimiento previo de K-Means ni de paralelismo

---

## Tabla de Contenidos

- [Parte 0 — El Dataset Iris: de qué estamos hablando](#parte-0--el-dataset-iris-de-qué-estamos-hablando)
- [Parte 1 — K-Means desde cero](#parte-1--k-means-desde-cero)
- [Parte 2 — Cómputo paralelo desde cero](#parte-2--cómputo-paralelo-desde-cero)
- [Parte 3 — OpenMP desde cero](#parte-3--openmp-desde-cero)
- [Parte 4 — El código C explicado línea a línea](#parte-4--el-código-c-explicado-línea-a-línea)
- [Parte 5 — ¿Por qué solo paralelizamos la asignación?](#parte-5--por-qué-solo-paralelizamos-la-asignación)
- [Parte 6 — ¿Por qué el speedup es menor a 1?](#parte-6--por-qué-el-speedup-es-menor-a-1)
- [Parte 7 — Compilar y ejecutar](#parte-7--compilar-y-ejecutar)
- [Parte 8 — Gráficas con Gnuplot](#parte-8--gráficas-con-gnuplot)
- [Parte 9 — Glosario completo](#parte-9--glosario-completo)

---

## Parte 0 — El Dataset Iris: de qué estamos hablando

### ¿Qué son "datos"?

Antes de hablar de algoritmos, hay que entender con qué trabajamos.

Un **dato** es simplemente información que podemos medir y escribir como número. Por ejemplo, si quisieras describir una flor con números, podrías medir:

- Cuánto mide el sépalo (la parte verde que rodea el pétalo) en largo
- Cuánto mide el sépalo en ancho
- Cuánto mide el pétalo en largo
- Cuánto mide el pétalo en ancho

Eso es exactamente lo que tiene el dataset Iris: **150 flores**, cada una descrita con **4 medidas** en centímetros.

### Cómo se ve el archivo iris.csv

El archivo `iris.csv` es una tabla de texto plano. Las primeras 5 líneas se ven así:

```
5.1,3.5,1.4,0.2,Iris-setosa
4.9,3.0,1.4,0.2,Iris-setosa
4.7,3.2,1.3,0.2,Iris-setosa
4.6,3.1,1.5,0.2,Iris-setosa
5.0,3.6,1.4,0.2,Iris-setosa
```

Cada línea es una flor. Las columnas son:

| Columna 1 | Columna 2 | Columna 3 | Columna 4 | Columna 5 |
|-----------|-----------|-----------|-----------|-----------|
| largo_sépalo | ancho_sépalo | largo_pétalo | ancho_pétalo | especie |
| 5.1 | 3.5 | 1.4 | 0.2 | Iris-setosa |

### Las 3 especies

El dataset tiene exactamente 150 flores divididas en 3 grupos de 50:

- **Filas 1-50:** Iris-setosa (pétalos muy pequeños)
- **Filas 51-100:** Iris-versicolor (pétalos medianos)
- **Filas 101-150:** Iris-virginica (pétalos grandes)

La columna 5 (la especie) se llama **etiqueta**: es la respuesta "correcta" de a qué grupo pertenece cada flor.

### Lo importante: K-Means IGNORA las etiquetas

El algoritmo K-Means solo usa las **4 columnas numéricas**. No le decimos "esta flor es setosa". Le damos solo los números y le preguntamos: *¿puedes encontrar grupos naturales tú solo?*

Eso se llama **aprendizaje no supervisado**: aprender sin que alguien te dé las respuestas correctas.

---

## Parte 1 — K-Means desde cero

### 1.1 El problema: tenemos una mezcla, queremos grupos

Imagina que tienes 150 canicas de colores mezcladas en una caja. Hay rojas, azules y verdes, pero están todas revueltas. Tu trabajo es separarlas en 3 montones, pero... **tienes los ojos cerrados**.

Lo único que puedes hacer es agarrar una canica y sentir si es pesada, grande, rugosa, etc. Con esas características (los "datos"), tienes que decidir a qué montón va.

K-Means hace exactamente eso: agrupa datos **sin ver** a qué clase verdadera pertenecen, solo mirando las características numéricas.

### 1.2 ¿Qué es la distancia entre dos puntos?

Para saber qué tan parecidas son dos flores, necesitamos medir qué tan "lejos" están entre sí en el espacio de sus características.

**Primero en 2D (dos medidas):**

Si tienes la flor A con pétalo largo=1.4, pétalo ancho=0.2, y la flor B con pétalo largo=4.7, pétalo ancho=1.4, la distancia entre ellas en un plano (como la distancia entre dos puntos en un mapa) se calcula con el **Teorema de Pitágoras**:

```
distancia = √( (4.7 - 1.4)² + (1.4 - 0.2)² )
          = √( (3.3)²        + (1.2)²        )
          = √( 10.89         + 1.44          )
          = √( 12.33 )
          = 3.51 cm
```

Visualmente: es la línea recta que une esos dos puntos en un gráfico.

**En 4D (cuatro medidas):**

Con 4 características, el principio es el mismo, solo se suman 4 términos:

```
distancia = √( (x₁-y₁)² + (x₂-y₂)² + (x₃-y₃)² + (x₄-y₄)² )
```

No podemos visualizarlo en 4 dimensiones, pero la matemática funciona igual.

**¿Por qué usamos el cuadrado de la distancia en el código?**

En el código verás que usamos `dist += diff * diff` y nunca calculamos la raíz cuadrada. ¿Por qué?

Porque solo queremos **comparar** distancias, no medirlas exactamente. Si flor A está a distancia 3.51 de un centroide y flor B está a distancia 5.20, la que está más cerca es A. Pero si comparamos los cuadrados: 12.33 vs 27.04, la conclusión es la misma (12.33 < 27.04, entonces A sigue ganando). El orden no cambia, y la raíz cuadrada es una operación costosa. Ahorramos tiempo sin perder precisión.

### 1.3 ¿Qué es un centroide?

Un **centroide** es el "punto central" de un grupo. Se calcula como el promedio de todas las flores del grupo.

**Ejemplo manual con 3 flores en el cluster:**

| Flor | largo_pétalo | ancho_pétalo |
|------|-------------|-------------|
| Flor 1 | 1.4 | 0.2 |
| Flor 2 | 1.4 | 0.2 |
| Flor 3 | 1.3 | 0.2 |

```
centroide_largo  = (1.4 + 1.4 + 1.3) / 3 = 4.1 / 3 = 1.367
centroide_ancho  = (0.2 + 0.2 + 0.2) / 3 = 0.6 / 3 = 0.200
```

El centroide de este grupo es el punto (1.367, 0.200). No es una flor real, es el promedio del grupo.

### 1.4 El algoritmo K-Means paso a paso

Ahora que entendemos distancia y centroide, el algoritmo completo es:

```
PASO 0: Elige K=3 flores al azar como centroides iniciales

REPETIR hasta convergencia:
  PASO 1 (Asignación):
    Para cada flor:
      - Calcula su distancia a los 3 centroides
      - Asígnala al centroide más cercano

  PASO 2 (Recálculo):
    Para cada grupo:
      - Calcula el promedio de todas las flores del grupo
      - Ese promedio es el nuevo centroide

  PASO 3 (Convergencia):
    ¿Se movieron mucho los centroides?
    Si no → terminamos
    Si sí → repetir desde PASO 1
```

**Traza manual con 6 flores inventadas:**

Usamos solo 2 características (largo_pétalo, ancho_pétalo) para poder mostrar los números:

```
Flores:
  F1=(1.0, 0.2)   F2=(1.2, 0.3)   F3=(1.1, 0.2)   ← grupo "pequeño"
  F4=(4.5, 1.5)   F5=(4.8, 1.4)   F6=(4.6, 1.6)   ← grupo "grande"

Centroides iniciales (elegidos al azar):
  C1=(1.0, 0.2)   ← tomamos F1
  C2=(4.5, 1.5)   ← tomamos F4
  C3=(1.2, 0.3)   ← tomamos F2
```

**Iteración 1, PASO 1 — Asignación:**

Para F3=(1.1, 0.2):
```
dist(F3, C1) = (1.1-1.0)² + (0.2-0.2)² = 0.01 + 0 = 0.01
dist(F3, C2) = (1.1-4.5)² + (0.2-1.5)² = 11.56 + 1.69 = 13.25
dist(F3, C3) = (1.1-1.2)² + (0.2-0.3)² = 0.01 + 0.01 = 0.02
→ Más cercano: C1 (distancia 0.01)  → F3 va al cluster 1
```

(Haciendo lo mismo para todas las flores obtenemos: F1,F2,F3 → cluster 1 o 3; F4,F5,F6 → cluster 2)

**Iteración 1, PASO 2 — Recálculo:**

Supón que después de asignar: F1,F2,F3 están en el cluster 1 y F4,F5,F6 en el cluster 2:
```
Nuevo C1 = promedio(F1, F2, F3) = ((1.0+1.2+1.1)/3, (0.2+0.3+0.2)/3)
         = (1.1, 0.233)

Nuevo C2 = promedio(F4, F5, F6) = ((4.5+4.8+4.6)/3, (1.5+1.4+1.6)/3)
         = (4.633, 1.5)
```

**¿Convergió?** El centroide C1 se movió de (1.0, 0.2) a (1.1, 0.233). Desplazamiento = 0.1. Como 0.1 > 0.001 (nuestro umbral), continuamos.

Después de otra iteración, los centroides se estabilizan y el programa para. Eso es convergencia.

### 1.5 ¿Por qué "K" en K-Means?

"K" es el número de grupos que queremos. Nosotros elegimos K=3 porque sabemos (por conocimiento externo del dataset) que hay 3 especies. En la práctica, elegir el K correcto es todo un arte, pero para este proyecto lo sabemos de antemano.

"Means" en inglés significa "promedios". K-Means = K promedios = K centroides que representan K grupos.

---

## Parte 2 — Cómputo paralelo desde cero

### 2.1 ¿Qué hace una CPU?

Una **CPU** (el procesador de tu computadora) ejecuta instrucciones, una por una, extremadamente rápido. "Instrucciones" son cosas simples como: sumar dos números, comparar dos valores, guardar un resultado en memoria.

Las CPUs modernas tienen varios **cores** (núcleos). Cada core puede ejecutar instrucciones de forma independiente. Una CPU con 4 cores puede ejecutar 4 flujos de instrucciones al mismo tiempo.

Piénsalo como 4 calculadoras que trabajan en paralelo en la misma mesa.

### 2.2 ¿Qué es un proceso?

Cuando ejecutas `./kmeans`, el sistema operativo crea un **proceso**: un programa cargado en memoria con su propio espacio de trabajo (variables, instrucciones, pila de llamadas).

Cada proceso tiene:
- Sus propias variables (no puede ver las variables de otros procesos)
- Su propio contador de instrucciones (sabe qué línea está ejecutando)
- Sus propios recursos (archivos abiertos, etc.)

### 2.3 ¿Qué es un hilo (thread)?

Un **hilo** (thread en inglés) es una "sub-tarea" dentro de un proceso. A diferencia de los procesos, los hilos de un mismo proceso **comparten la memoria**.

Imagina un proceso como una oficina y los hilos como los empleados. Todos trabajan en la misma oficina (comparten los archivos, las calculadoras, la fotocopiadora), pero cada uno tiene su propio escritorio (su propia pila de trabajo actual).

En nuestro programa:
- El proceso tiene la tabla `data[150][4]` en memoria
- Cuando creamos 4 hilos, los 4 pueden leer esa tabla sin copiarla
- Eso es eficiente: no gastamos tiempo copiando 150 flores para cada hilo

### 2.4 ¿Qué es paralelismo?

**Paralelismo** es hacer varias cosas al mismo tiempo, en cores distintos de la CPU.

Sin paralelismo (serial):
```
Core 0: Flor 1 → Flor 2 → Flor 3 → ... → Flor 150
Core 1: inactivo
Core 2: inactivo
Core 3: inactivo
```

Con paralelismo (4 hilos):
```
Core 0: Flor 1  → Flor 2  → ... → Flor 37   (flores 1-37)
Core 1: Flor 38 → Flor 39 → ... → Flor 75   (flores 38-75)
Core 2: Flor 76 → Flor 77 → ... → Flor 112  (flores 76-112)
Core 3: Flor 113→ Flor 114→ ... → Flor 150  (flores 113-150)
```

Si cada flor tarda el mismo tiempo, el resultado es ~4 veces más rápido. En teoría.

### 2.5 El problema de compartir datos: race conditions

Imagina dos cajeros en un banco con la misma cuenta:

```
Cuenta: $100

Cajero A: lee saldo → $100, suma $50, escribe → $150
Cajero B: lee saldo → $100, suma $30, escribe → $130

Resultado final: $130  ← ¡MAL! Debería ser $180
```

Esto es una **race condition** (condición de carrera): dos hilos leen el mismo dato, cada uno lo modifica en su propia calculadora, y el último en escribir borra el trabajo del primero.

En programación:

```c
// DOS HILOS ejecutan esto SIMULTÁNEAMENTE:
contador = contador + 1;

// Lo que realmente pasa:
// Hilo 1: lee contador (valor 5)
// Hilo 2: lee contador (valor 5)
// Hilo 1: calcula 5+1=6, escribe 6
// Hilo 2: calcula 5+1=6, escribe 6
// Resultado: contador=6, ¡pero debería ser 7!
```

Las race conditions producen resultados incorrectos e impredecibles. **El programador es responsable de evitarlas**.

### 2.6 El overhead: el costo de contratar trabajadores

Crear hilos no es gratis. El sistema operativo necesita:
1. Reservar memoria para la pila de cada hilo
2. Registrar el hilo en el planificador del SO
3. Distribuir el trabajo entre los hilos
4. Esperar que todos terminen (sincronización)
5. Destruir los hilos y liberar memoria

Todo esto toma tiempo: típicamente entre 10 y 50 microsegundos por hilo. Eso suena poco, pero si el trabajo que vas a paralelizar solo dura 1 microsegundo, estás gastando 10-50x más tiempo en preparativos que en trabajo real.

Esto explica por qué en nuestro proyecto con N=150 flores, usar más hilos hace el programa **más lento**, no más rápido.

---

## Parte 3 — OpenMP desde cero

### 3.1 ¿Qué es OpenMP?

**OpenMP** (Open Multi-Processing) es un estándar de la industria para escribir programas paralelos en C, C++ y Fortran. Consiste en:

1. Un conjunto de **directivas** (instrucciones especiales para el compilador)
2. Una **biblioteca de funciones** en tiempo de ejecución
3. Variables de **entorno** para configurar el comportamiento

Lo especial de OpenMP es que es muy fácil de usar: tomas un bucle serial existente y agregas **una sola línea** para paralelizarlo.

OpenMP no es un lenguaje nuevo. Es una extensión de C/C++/Fortran que el compilador entiende.

### 3.2 ¿Qué es un `#pragma`?

Un `#pragma` es una instrucción al compilador que no afecta la lógica del programa por sí misma, sino que le dice al compilador qué hacer con el código.

```c
#pragma omp parallel for
for (int i = 0; i < 150; i++) {
    /* código */
}
```

Si compilas sin `-fopenmp`, el compilador simplemente **ignora** el `#pragma` y ejecuta el bucle normalmente (como serial). Esto es genial: el mismo código puede correr serial o paralelo según cómo se compile.

Con `-fopenmp`, el compilador **transforma** ese bucle en código que crea hilos, divide las iteraciones y las ejecuta en paralelo.

### 3.3 El modelo de ejecución de OpenMP

OpenMP usa el modelo **fork-join**:

```
Hilo principal (serial)
        │
        │ ←─── #pragma omp parallel
        ▼
   ┌────┴────┐
   │  Fork   │ ← Se crean N hilos
   ▼    ▼    ▼
 Hilo0 Hilo1 Hilo2   ← Todos trabajan en paralelo
   │    │    │
   └────┴────┘
        │ ←─── Fin de la zona paralela (todos esperan aquí)
        │        Esto se llama "barrera implícita"
        ▼
Hilo principal (serial)
```

Fuera de la zona `#pragma`, hay un solo hilo. Dentro, hay N hilos trabajando a la vez.

### 3.4 `omp_set_num_threads(n)` — cuántos hilos usar

```c
omp_set_num_threads(4);  // "de ahora en adelante usa 4 hilos"
```

Esto le dice a OpenMP cuántos hilos crear en la próxima zona paralela. En nuestro código lo llamamos antes de cada experimento (1 hilo, 2 hilos, 4 hilos).

También se puede controlar con la variable de entorno `OMP_NUM_THREADS`:
```bash
OMP_NUM_THREADS=4 ./kmeans
```

### 3.5 `#pragma omp parallel for schedule(static)` — el reparto del trabajo

Esta es la directiva más importante de nuestro código. Vamos parte por parte:

**`parallel`:** Crea un equipo de hilos. A partir de aquí, todos los hilos ejecutan el código.

**`for`:** Le dice a OpenMP que el bucle que sigue se debe **dividir** entre los hilos. Sin esto, cada hilo ejecutaría el bucle completo (¡eso no es lo que queremos!).

**`schedule(static)`:** Define cómo se dividen las iteraciones.

- **static** = división estática en bloques iguales, asignados antes de empezar.

Con 4 hilos y 150 iteraciones:
```
Hilo 0 → iteraciones i = 0, 1, 2, ..., 37   (37-38 iteraciones)
Hilo 1 → iteraciones i = 38, 39, ..., 75
Hilo 2 → iteraciones i = 76, 77, ..., 112
Hilo 3 → iteraciones i = 113, 114, ..., 149
```

Cada hilo trabaja en su rango sin pisarse con los demás.

**¿Por qué `schedule(static)` y no otras opciones?**

| schedule | Comportamiento | Cuándo usarlo |
|----------|---------------|---------------|
| static | Bloques iguales, asignados antes | Cuando todas las iteraciones tardan igual |
| dynamic | Iteraciones asignadas "al vuelo" | Cuando algunas tardan más que otras |
| guided | Bloques grandes al inicio, pequeños al final | Mezcla balanceada |

Para K-Means, cada flor tarda exactamente el mismo tiempo en procesar, así que `static` es perfecto.

### 3.6 Variables privadas vs compartidas

En una zona paralela, las variables pueden ser:

**Compartidas (shared):** Todos los hilos ven el mismo espacio de memoria. Por defecto, las variables declaradas **fuera** del `#pragma` son compartidas.

```c
// FUERA del #pragma → compartidas → todos los hilos las ven
double data[150][4];      // todos los hilos leen esto
int assignments[150];     // cada hilo escribe en su propio índice
double centroids[3][4];   // todos los hilos leen esto
```

**Privadas (private):** Cada hilo tiene su propia copia. Las variables declaradas **dentro** del bucle son privadas automáticamente.

```c
#pragma omp parallel for schedule(static)
for (int i = 0; i < N_POINTS; i++) {
    double min_dist = 1e18;  // ← PRIVADA: cada hilo tiene su propia min_dist
    int    best     = 0;     // ← PRIVADA: cada hilo tiene su propio best
    // ...
}
```

La regla: si dos hilos pueden escribir en el mismo lugar de memoria al mismo tiempo, tenemos un problema. Si cada hilo escribe en lugares diferentes (como `assignments[i]` donde cada hilo tiene un `i` diferente), no hay conflicto.

### 3.7 `omp_get_wtime()` — el cronómetro

```c
double t_start = omp_get_wtime();  // marca tiempo de inicio
// ... código a medir ...
double t_end = omp_get_wtime();    // marca tiempo de fin
double duracion = t_end - t_start; // diferencia en segundos
```

`omp_get_wtime()` retorna el tiempo actual en segundos como número decimal. La diferencia entre dos llamadas da el tiempo transcurrido.

Es como presionar "inicio" y "fin" en un cronómetro. La ventaja sobre `clock()` de C es que `omp_get_wtime()` mide **tiempo de reloj real** (wall time), no tiempo de CPU, por lo que refleja mejor el rendimiento real del programa paralelo.

### 3.8 Por qué `-fopenmp` al compilar

```bash
gcc -O2 -fopenmp -o kmeans kmeans.c -lm
```

Sin `-fopenmp`:
- El compilador ignora todos los `#pragma omp`
- Las funciones `omp_*` no existen → el programa no compila
- El bucle se ejecuta serial

Con `-fopenmp`:
- El compilador entiende los `#pragma omp` y genera código paralelo
- Enlaza automáticamente la biblioteca `libgomp` (la implementación de OpenMP de GCC)
- El bucle se paralela automáticamente

---

## Parte 4 — El código C explicado línea a línea

### 4.1 Las cabeceras (`#include`)

```c
#include <stdio.h>    // entrada/salida: printf, fprintf, fscanf, fopen, fclose
#include <stdlib.h>   // funciones de sistema: exit()
#include <math.h>     // funciones matemáticas: fabs() (valor absoluto)
#include <omp.h>      // funciones de OpenMP: omp_get_wtime(), omp_set_num_threads()
```

Un `#include` es como decirle al compilador: "necesito usar funciones de esta biblioteca". Sin `<stdio.h>`, `printf` no existiría. Sin `<omp.h>`, ninguna función de OpenMP estaría disponible.

### 4.2 Las constantes (`#define`)

```c
#define N_POINTS   150   // cuántas flores tiene el dataset
#define N_DIMS     4     // cuántas medidas tiene cada flor
#define K          3     // cuántos clusters queremos
#define MAX_ITER   100   // máximo de vueltas del algoritmo
#define THRESHOLD  0.001 // "suficientemente quieto" para parar
#define LABEL_LEN  64    // tamaño máximo del nombre de especie
```

`#define` crea un **sustituto de texto**: antes de compilar, el compilador reemplaza cada `N_POINTS` por `150` en todo el código. No es una variable; es como hacer "buscar y reemplazar" en el código fuente.

¿Por qué usar `#define` en lugar de escribir `150` directamente?

1. **Legibilidad**: `N_POINTS` es más claro que `150`
2. **Mantenibilidad**: si el dataset cambia a 200 flores, solo cambias una línea
3. **Evita errores**: si escribes `150` en 20 lugares y luego necesitas cambiarlo, puedes olvidarte de alguno

### 4.3 Los arrays globales

```c
double data[N_POINTS][N_DIMS];          // la tabla de datos (150 filas, 4 columnas)
int    assignments[N_POINTS];           // a qué cluster pertenece cada flor
double centroids[K][N_DIMS];            // posición actual de los 3 centroides
double init_centroids[K][N_DIMS];       // copia guardada para repetir el experimento
```

**¿Qué es un array 2D?**

`data[150][4]` es como una tabla con 150 filas y 4 columnas:

```
data[0][0]=5.1  data[0][1]=3.5  data[0][2]=1.4  data[0][3]=0.2  ← Flor 0
data[1][0]=4.9  data[1][1]=3.0  data[1][2]=1.4  data[1][3]=0.2  ← Flor 1
data[2][0]=4.7  data[2][1]=3.2  data[2][2]=1.3  data[2][3]=0.2  ← Flor 2
...
```

Para acceder a la medida `d` de la flor `i`: `data[i][d]`

**¿Por qué variables globales?**

Las variables globales existen mientras el programa corre y son accesibles desde cualquier función. Las variables locales existen solo dentro de su función y desaparecen cuando la función termina.

Usamos globales aquí porque:
1. Los arrays son grandes (150×4 = 600 doubles = 4800 bytes) — si fueran locales de `main`, habría que pasarlos como parámetros a cada función
2. OpenMP puede acceder a variables globales desde todos los hilos sin problema
3. Simplifica el código: no hay que pasar punteros por todas partes

**`double` vs `int`:**
- `double`: número decimal de doble precisión (64 bits). Ej: 5.1, 3.14159, 0.001
- `int`: número entero (32 bits). Ej: 0, 1, 2, 150

Usamos `double` para medidas (pueden tener decimales) e `int` para índices y contadores (siempre enteros).

### 4.4 La función `load_iris` — leer el CSV

```c
int load_iris(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: no se pudo abrir '%s'\n", filename);
        exit(1);
    }
```

- `FILE *f = fopen(filename, "r")`: abre el archivo en modo lectura. `fopen` retorna un puntero al archivo, o `NULL` si falló.
- `if (!f)`: si `f` es NULL (el archivo no existe o no se puede abrir), mostramos error y salimos.
- `fprintf(stderr, ...)`: escribe en la salida de errores (stderr), no en la salida normal (stdout). Esto es buena práctica para mensajes de error.
- `exit(1)`: termina el programa con código de error 1 (0 = éxito, cualquier otro número = error).

```c
    int n = 0;
    char label[LABEL_LEN];  // buffer temporal para el nombre de especie

    while (n < N_POINTS &&
           fscanf(f, "%lf,%lf,%lf,%lf,%63s",
                  &data[n][0], &data[n][1],
                  &data[n][2], &data[n][3],
                  label) == 5) {
        n++;
    }
```

`fscanf` es como `scanf` pero lee de un archivo en lugar de del teclado.

El formato `"%lf,%lf,%lf,%lf,%63s"` significa:
- `%lf`: lee un número decimal (`double`) — el `l` es de "long", para `double`
- `,`: espera una coma literal en el archivo
- `%63s`: lee una cadena de texto de máximo 63 caracteres (el nombre de la especie)

El `== 5` verifica que se leyeron exactamente **5 campos**. Si el archivo termina o hay una línea mal formada, `fscanf` retorna menos de 5 y el bucle para.

`&data[n][0]` — el `&` significa "dame la **dirección de memoria** de esta variable". `fscanf` necesita saber dónde guardar el número que lee, no el valor actual de la variable.

`label` se usa para leer el nombre de la especie, pero nunca lo usamos después. Solo lo leemos para que `fscanf` avance al siguiente registro.

```c
    fclose(f);

    if (n != N_POINTS) {
        fprintf(stderr, "Advertencia: se esperaban %d filas, se leyeron %d\n",
                N_POINTS, n);
    }
    return n;
}
```

`fclose(f)`: siempre hay que cerrar los archivos que abrimos. Si no, el sistema operativo puede quedarse con el archivo "bloqueado".

La función retorna `n` (cuántas flores se leyeron) para que `main` pueda verificarlo.

### 4.5 La función `init_centroids_from_data` — elegir centroides iniciales

```c
void init_centroids_from_data(int seed_indices[K]) {
    for (int c = 0; c < K; c++) {
        for (int d = 0; d < N_DIMS; d++) {
            centroids[c][d]      = data[seed_indices[c]][d];
            init_centroids[c][d] = data[seed_indices[c]][d];
        }
    }
}
```

Esta función toma un array de índices (`seed_indices = {0, 50, 100}`) y copia las coordenadas de esas flores a los centroides.

El bucle doble:
- El bucle externo `c` recorre los 3 clusters (c = 0, 1, 2)
- El bucle interno `d` recorre las 4 dimensiones (d = 0, 1, 2, 3)

Para c=0, d=0: `centroids[0][0] = data[0][0]` → El centroide 0, dimensión 0 = la flor 0, dimensión 0 = 5.1

También guarda una copia en `init_centroids`. ¿Por qué? Porque el programa ejecuta K-Means 3 veces (con 1, 2 y 4 hilos) y necesita que todas las ejecuciones partan del **mismo punto inicial** para que la comparación de tiempos sea justa. Si no restauráramos los centroides, la segunda ejecución empezaría con los centroides ya convergidos de la primera, terminaría en 1 iteración y parecería más rápida por razones incorrectas.

Usamos índices `{0, 50, 100}` porque coinciden con la primera flor de cada especie, garantizando que los 3 centroides iniciales estén bien separados y el algoritmo converja correctamente.

### 4.6 La función `run_kmeans` — el corazón del programa

Esta es la función más importante. La vamos a analizar en bloques.

**Bloque 1: Preparación**

```c
double run_kmeans(int num_threads) {

    // Restaurar centroides iniciales para comparación justa
    for (int c = 0; c < K; c++)
        for (int d = 0; d < N_DIMS; d++)
            centroids[c][d] = init_centroids[c][d];

    // Configurar cuántos hilos usará OpenMP
    omp_set_num_threads(num_threads);

    // Iniciar el cronómetro
    double t_start = omp_get_wtime();
```

Primero restauramos los centroides (explicado arriba). Luego configuramos los hilos. Luego activamos el cronómetro.

**Bloque 2: El bucle principal**

```c
    for (int iter = 0; iter < MAX_ITER; iter++) {
```

Este bucle se ejecuta hasta 100 veces. Dentro de él, las dos fases del K-Means.

**Bloque 3: Fase de asignación (paralela)**

```c
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N_POINTS; i++) {
            double min_dist = 1e18;
            int    best     = 0;
```

`1e18` es la notación científica para 10^18: un número gigantesco. Lo usamos como "infinito" inicial para que cualquier distancia real sea menor.

`min_dist` y `best` son variables **locales al bucle** → son privadas por hilo → no hay race condition.

```c
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
```

Para cada flor `i`, calculamos la distancia al cuadrado a cada centroide `c`:
- Bucle externo (c): recorre los 3 centroides
- Bucle interno (d): suma los cuadrados de diferencias en cada dimensión
- `diff = data[i][d] - centroids[c][d]`: diferencia en dimensión d
- `dist += diff * diff`: acumula el cuadrado de la diferencia
- Al final del bucle c, si `dist` es menor que la mínima encontrada, actualizamos `best`

```c
            assignments[i] = best;
        }
```

Guardamos el cluster ganador. `assignments[i]` es escrito solo por el hilo que procesa la iteración `i`. No hay dos hilos escribiendo en el mismo índice. No hay race condition.

**Bloque 4: Fase de recálculo (serial)**

```c
        double new_c[K][N_DIMS];
        int    counts[K];

        for (int c = 0; c < K; c++) {
            counts[c] = 0;
            for (int d = 0; d < N_DIMS; d++)
                new_c[c][d] = 0.0;
        }
```

Inicializamos acumuladores a cero. `new_c[c][d]` va a acumular la suma de la dimensión `d` de todas las flores del cluster `c`. `counts[c]` va a contar cuántas flores hay en el cluster `c`.

```c
        for (int i = 0; i < N_POINTS; i++) {
            int c = assignments[i];
            counts[c]++;
            for (int d = 0; d < N_DIMS; d++)
                new_c[c][d] += data[i][d];
        }
```

Para cada flor, miramos a qué cluster fue asignada, incrementamos su contador, y sumamos sus coordenadas a los acumuladores del cluster.

```c
        double max_shift = 0.0;
        for (int c = 0; c < K; c++) {
            if (counts[c] == 0) continue;
            for (int d = 0; d < N_DIMS; d++) {
                double new_val = new_c[c][d] / counts[c];
                double shift   = fabs(new_val - centroids[c][d]);
                if (shift > max_shift) max_shift = shift;
                centroids[c][d] = new_val;
            }
        }
```

Para cada cluster, calculamos el nuevo centroide dividiendo la suma por el conteo (eso es el promedio). Medimos cuánto se movió cada coordenada del centroide con `fabs()` (valor absoluto). Si el mayor movimiento es `max_shift`, lo comparamos con el umbral.

`if (counts[c] == 0) continue;`: si un cluster quedó vacío (ninguna flor le fue asignada), saltamos. De lo contrario dividiríamos por cero, lo cual es un error matemático.

**Bloque 5: Convergencia y fin**

```c
        if (max_shift < THRESHOLD) break;

    } // fin del bucle K-Means

    double t_end = omp_get_wtime();
    return t_end - t_start;
}
```

Si ningún centroide se movió más de 0.001, el algoritmo convergió: salimos del bucle anticipadamente con `break`. Si llegamos a 100 iteraciones sin converger, el bucle termina naturalmente.

Finalmente, calculamos el tiempo transcurrido y lo retornamos.

### 4.7 La función `export_results` — guardar para gnuplot

```c
void export_results(const char *filename) {
    FILE *f = fopen(filename, "w");  // "w" = escritura (crea o sobreescribe)
    if (!f) {
        fprintf(stderr, "Error: no se pudo crear '%s'\n", filename);
        return;
    }
    for (int i = 0; i < N_POINTS; i++) {
        fprintf(f, "%.2f %.2f %d\n",
                data[i][2],            // largo_pétalo (columna 2)
                data[i][3],            // ancho_pétalo (columna 3)
                assignments[i] + 1);   // cluster 1-indexado
    }
    fclose(f);
}
```

Escribimos las dimensiones 2 y 3 (largo y ancho del pétalo) porque estas dos características son las que mejor separan visualmente las 3 especies en un gráfico 2D.

`assignments[i] + 1`: los clusters internamente son 0, 1, 2. Le sumamos 1 para que sean 1, 2, 3 en el archivo. Gnuplot usa esta convención en su sintaxis de filtrado.

`%.2f`: formato de punto flotante con 2 decimales. `%d`: entero sin decimales.

El archivo resultante (`resultados.dat`) tiene este aspecto:
```
1.40 0.20 1
1.40 0.20 1
...
4.70 1.40 2
...
```

### 4.8 La función `export_speedup` — guardar métricas

```c
void export_speedup(const char *filename,
                    int threads[], double times[], int n) {
    FILE *f = fopen(filename, "w");
    if (!f) { fprintf(stderr, "Error: ...\n"); return; }

    double t1 = times[0];  // tiempo base: el de 1 hilo
    for (int i = 0; i < n; i++) {
        double speedup = (times[i] > 0.0) ? (t1 / times[i]) : 1.0;
        fprintf(f, "%d %.4f\n", threads[i], speedup);
    }
    fclose(f);
}
```

**Fórmula del speedup:**

```
Speedup(N hilos) = Tiempo con 1 hilo / Tiempo con N hilos
```

Si el programa tarda 1 segundo con 1 hilo y 0.5 segundos con 2 hilos, el speedup es 2.0 (el doble de rápido). El speedup ideal con N hilos sería N.

`(times[i] > 0.0) ? (t1 / times[i]) : 1.0`: el operador ternario en C. Significa: "si `times[i]` es mayor que 0, calcula `t1/times[i]`; si no, usa 1.0". Es una protección contra división por cero.

### 4.9 La función `main` — el director de orquesta

```c
int main(void) {

    // PASO 1: Cargar dataset
    printf("Cargando iris.csv...\n");
    int n = load_iris("iris.csv");
    printf("  %d puntos cargados, %d dimensiones\n\n", n, N_DIMS);
```

`main` es la función que ejecuta el sistema operativo cuando corres `./kmeans`. Llama a las demás funciones en el orden correcto.

```c
    // PASO 2: Inicializar centroides
    int seeds[K] = {0, 50, 100};
    init_centroids_from_data(seeds);
```

Elegimos la flor 0 (primera setosa), flor 50 (primera versicolor) y flor 100 (primera virginica) como centroides iniciales.

```c
    // PASO 3: Ejecutar K-Means con 1, 2 y 4 hilos
    int    thread_counts[3] = {1, 2, 4};
    double times[3];

    printf("Ejecutando K-Means (K=%d, max_iter=%d)...\n\n", K, MAX_ITER);

    for (int i = 0; i < 3; i++) {
        times[i] = run_kmeans(thread_counts[i]);
    }
```

Ejecutamos K-Means 3 veces con diferente número de hilos. Los tiempos se guardan en `times[]`.

```c
    // PASO 4: Mostrar tabla de resultados
    printf("+---------+-------------+---------+\n");
    printf("| Hilos   | Tiempo (s)  | Speedup |\n");
    printf("+---------+-------------+---------+\n");
    for (int i = 0; i < 3; i++) {
        double speedup = times[0] / times[i];
        printf("| %7d | %11.6f | %7.3f |\n",
               thread_counts[i], times[i], speedup);
    }
    printf("+---------+-------------+---------+\n");
```

Imprime la tabla de resultados con formato de caja ASCII.

`%7d`: entero alineado a la derecha en un campo de 7 caracteres.  
`%11.6f`: decimal con 6 cifras después del punto, campo de 11 caracteres.  
`%7.3f`: decimal con 3 cifras después del punto, campo de 7 caracteres.

```c
    // PASO 5: Exportar datos para gnuplot
    export_results("resultados.dat");
    export_speedup("speedup.dat", thread_counts, times, 3);

    printf("\nEjecutar 'gnuplot graficar.gp' para generar las graficas.\n");
    return 0;
}
```

`return 0`: convención en C para indicar que el programa terminó exitosamente. El sistema operativo recibe este valor.

---

## Parte 5 — ¿Por qué solo paralelizamos la asignación?

### 5.1 El análisis de operaciones

Contemos cuántas operaciones hace cada fase en cada iteración del K-Means:

**Fase de asignación:** Para cada una de las N=150 flores, calculamos la distancia a cada uno de los K=3 centroides, y cada distancia requiere D=4 multiplicaciones y sumas:

```
Operaciones = N × K × D = 150 × 3 × 4 = 1800 multiplicaciones/sumas
```

**Fase de recálculo:** Para cada una de las N=150 flores, sumamos sus D=4 coordenadas al acumulador de su cluster:

```
Operaciones = N × D = 150 × 4 = 600 sumas
```

La asignación hace 3 veces más trabajo que el recálculo. Paralizar la asignación tiene mayor impacto potencial.

### 5.2 La condición de independencia

Para paralelizar un bucle con OpenMP sin errores, las iteraciones deben ser **independientes**: la iteración `i` no puede depender del resultado de la iteración `j`.

**En la fase de asignación:**

```c
// Iteración i:
//   LEE: data[i][...] (solo su propio dato)
//   LEE: centroids[...][...] (todos leen, nadie escribe)
//   ESCRIBE: assignments[i] (su propio índice, nadie más escribe aquí)
```

Cada iteración trabaja en un índice diferente de `assignments`. No hay conflicto. Es **seguro paralelizar**.

**En la fase de recálculo:**

```c
// Iteración i (sin paralelismo):
int c = assignments[i];
counts[c]++;              // ← PROBLEMA
new_c[c][d] += data[i][d]; // ← PROBLEMA
```

Si dos flores i=5 e i=42 pertenecen al mismo cluster c=0, ambas iteraciones intentarían modificar `counts[0]` y `new_c[0][d]` al mismo tiempo. Race condition garantizada.

**¿Se podría paralizar con `reduction`?**

Sí, OpenMP tiene la directiva `reduction` para sumar acumuladores de forma segura:
```c
// Esta versión sería correcta pero más compleja:
#pragma omp parallel for reduction(+:counts) schedule(static)
for (int i = 0; i < N_POINTS; i++) { ... }
```

Sin embargo, el recálculo hace 600 operaciones vs 1800 de la asignación. Añadir complejidad (y overhead) para paralizar el 25% del trabajo cuando el 75% ya está paralelizado no vale la pena. La solución simple (serial) es la correcta aquí.

---

## Parte 6 — ¿Por qué el speedup es menor a 1?

### 6.1 Los resultados reales

Nuestras mediciones (guardadas en `speedup.dat`):

```
1 hilo  → speedup = 1.000  (referencia)
2 hilos → speedup = 0.210  (5 veces MÁS LENTO que 1 hilo)
4 hilos → speedup = 0.150  (6.7 veces MÁS LENTO que 1 hilo)
```

Esto parece un fallo, pero es correcto. ¿Por qué?

### 6.2 Anatomía del tiempo de ejecución

Cuando el programa corre con 4 hilos, el tiempo total se divide así:

| Componente | Tiempo aproximado |
|------------|-----------------|
| Crear 4 hilos | ~20-40 μs |
| Distribuir trabajo (150 flores) | ~1-5 μs |
| Calcular distancias de 150 flores | ~0.5-2 μs |
| Esperar que todos terminen (barrera) | ~5-10 μs |
| Destruir hilos | ~10-20 μs |
| **Total overhead** | ~36-77 μs |
| **Trabajo real** | ~0.5-2 μs |

El overhead de gestión de hilos es **10-100 veces mayor** que el trabajo que hacemos. Es como contratar 4 albañiles para poner un solo ladrillo: el tiempo de llamarlos, explicarles la tarea y despedirlos supera el tiempo de poner el ladrillo.

### 6.3 La Ley de Amdahl — desde cero

La **Ley de Amdahl** (formulada por Gene Amdahl en 1967) dice que el speedup máximo de un programa depende de cuánto de él es paralelizable.

**La intuición:** Si el 90% de tu programa es paralelo y el 10% es serial, con infinitos procesadores el programa nunca puede ser más de 10x más rápido, porque el 10% serial siempre tarda lo mismo.

**La fórmula:**

```
S(N) = 1 / ( (1 - P) + P/N )

Donde:
  S(N) = speedup con N hilos
  P    = fracción del programa que es paralela (entre 0 y 1)
  N    = número de hilos
```

**Con nuestros datos:**

El trabajo paralelo (asignación) dura ~1800 ops. El trabajo serial (resto) incluye recálculo (~600 ops) + overhead de OpenMP. Con N=150 y overhead dominante, la fracción paralela efectiva es muy pequeña.

Si la fracción paralela es P=0.01 (muy pequeña por el overhead):

```
S(2) = 1 / ( (1-0.01) + 0.01/2 ) = 1 / (0.99 + 0.005) = 1 / 0.995 ≈ 1.005
S(4) = 1 / ( (1-0.01) + 0.01/4 ) = 1 / (0.99 + 0.0025) ≈ 1.007
```

Con fracción paralela tan pequeña, agregar más hilos prácticamente no ayuda. De hecho, si el overhead aumenta con N hilos, puede ir hacia atrás.

### 6.4 ¿Cuándo sí valdría la pena?

Con N=1,000,000 de flores, el cálculo de distancias dominaría completamente:
- Cálculo de distancias: ~12,000,000 operaciones (~12 ms)
- Overhead de hilos: ~0.05 ms

```
Fracción paralela: P ≈ 12 / 12.05 ≈ 0.996
S(4) = 1 / (0.004 + 0.996/4) ≈ 3.7  → ¡Casi 4x más rápido!
```

**Conclusión:** OpenMP es poderoso, pero solo cuando el trabajo supera al overhead. Para N=150, el dataset es demasiado pequeño. Para N>100,000, el speedup se acerca al número de hilos.

---

## Parte 7 — Compilar y ejecutar

### 7.1 Verificar que tengas lo necesario

```bash
gcc --version
```

Debe mostrar algo como `gcc (GCC) 13.x.x` o similar.

```bash
make --version
```

Debe mostrar `GNU Make 4.x`.

Si no tienes GCC con OpenMP:

**Fedora/RHEL (como este sistema):**
```bash
sudo dnf install gcc make
```

**Ubuntu/Debian:**
```bash
sudo apt install gcc make
```

**macOS:**
```bash
brew install gcc
# Nota: en macOS el gcc que viene con Xcode no tiene OpenMP.
# Usar el de Homebrew: gcc-13 en lugar de gcc
```

### 7.2 Compilar con Make (recomendado)

```bash
cd /ruta/a/proyKMeans
make
```

Internamente ejecuta:
```bash
gcc -O2 -fopenmp -o kmeans kmeans.c -lm
```

**¿Qué significa cada parte?**

| Flag | Significado |
|------|-------------|
| `gcc` | El compilador de C |
| `-O2` | Optimización nivel 2: el compilador reordena instrucciones para ser más rápido. No cambia el comportamiento del programa, solo lo hace más eficiente |
| `-fopenmp` | Activa soporte OpenMP: el compilador entiende `#pragma omp` y enlaza `libgomp` |
| `-o kmeans` | El archivo ejecutable se llamará `kmeans` |
| `kmeans.c` | El archivo fuente a compilar |
| `-lm` | Enlaza la biblioteca matemática (necesaria para `fabs()`) |

**¿Por qué `-lm` va al final?**

GCC resuelve las bibliotecas de izquierda a derecha. Si escribieras `-lm kmeans.c`, la biblioteca matemática se procesaría antes de que el compilador sepa que necesita `fabs`, y no la enlazaría. Al ponerla al final, el compilador ya sabe qué funciones necesita y las busca en `-lm` correctamente.

### 7.3 Compilar manualmente (sin Make)

Si no tienes `make`:
```bash
gcc -O2 -fopenmp -o kmeans kmeans.c -lm
```

### 7.4 Ejecutar el programa

```bash
./kmeans
```

O con Make:
```bash
make run
```

### 7.5 Interpretar la salida

```
Cargando iris.csv...
  150 puntos cargados, 4 dimensiones

Ejecutando K-Means (K=3, max_iter=100)...

+---------+-------------+---------+
| Hilos   | Tiempo (s)  | Speedup |
+---------+-------------+---------+
|       1 |    0.000124 |   1.000 |
|       2 |    0.000592 |   0.210 |
|       4 |    0.000827 |   0.150 |
+---------+-------------+---------+

Nota: con N=150 el overhead de hilos puede superar la ganancia.
Para N > 100000 el speedup se acerca al numero de hilos.

Resultados exportados a: resultados.dat
Speedup exportado a:     speedup.dat

Ejecutar 'gnuplot graficar.gp' para generar las graficas.
```

- **Tiempo con 1 hilo:** ~0.000124 segundos = ~0.124 milisegundos = ~124 microsegundos
- **Tiempo con 2 hilos:** más lento (overhead de crear y sincronizar 2 hilos)
- **Speedup = 0.210:** con 2 hilos el programa es 0.21x la velocidad de 1 hilo, o sea, ~5x más lento

Los valores exactos varían según la máquina, carga del sistema y velocidad de la CPU.

---

## Parte 8 — Gráficas con Gnuplot

### 8.1 ¿Qué es Gnuplot?

Gnuplot es un programa de línea de comandos para crear gráficas a partir de archivos de datos. Lee un script (`.gp`) con instrucciones sobre qué graficar y cómo formatearlo, y produce imágenes PNG, PDF, SVG, etc.

### 8.2 Instalar Gnuplot

**Fedora/RHEL:**
```bash
sudo dnf install gnuplot
```

**Ubuntu/Debian:**
```bash
sudo apt install gnuplot
```

**macOS:**
```bash
brew install gnuplot
```

### 8.3 Generar las gráficas

Primero ejecuta el programa para crear los archivos de datos:
```bash
./kmeans
```

Luego genera las gráficas:
```bash
gnuplot graficar.gp
```

O en un solo comando con Make:
```bash
make all-plots
```

### 8.4 Las dos gráficas generadas

**`clusters.png` — Scatter plot de flores por cluster**

Muestra las 150 flores como puntos en un gráfico 2D donde:
- Eje X = largo del pétalo (cm)
- Eje Y = ancho del pétalo (cm)
- Color rojo = Cluster 1 (Iris-setosa: pétalos pequeños, abajo-izquierda)
- Color azul = Cluster 2 (Iris-versicolor: pétalos medianos)
- Color verde = Cluster 3 (Iris-virginica: pétalos grandes, arriba-derecha)

Los 3 grupos deben verse claramente separados. Esto demuestra que K-Means encontró los grupos naturales correctamente sin conocer las etiquetas.

**`speedup.png` — Curva de speedup**

Muestra dos líneas:
- Línea gris punteada: speedup **ideal** (lineal: 2 hilos = 2x, 4 hilos = 4x)
- Línea roja: speedup **real** medido

La línea real debería estar muy por debajo de la ideal, demostrando el efecto del overhead de OpenMP con N=150.

### 8.5 Ver las imágenes

```bash
# Linux con entorno gráfico:
eog clusters.png         # Eye of GNOME
xdg-open clusters.png    # abre con el visor predeterminado

# macOS:
open clusters.png

# En el terminal sin interfaz gráfica, copia el archivo a tu máquina local
```

---

## Parte 9 — Glosario completo

### Conceptos de K-Means

| Término | Definición |
|---------|------------|
| **Clustering** | Agrupar datos similares sin conocer las etiquetas correctas de antemano. |
| **Centroide** | El punto "representante" de un cluster, calculado como el promedio de todos los puntos del grupo. |
| **Feature / Característica** | Una columna del dataset: una medida numérica que describe cada elemento. En Iris: largo_sépalo, ancho_sépalo, largo_pétalo, ancho_pétalo. |
| **Distancia Euclidiana** | La "línea recta" entre dos puntos, calculada con el Teorema de Pitágoras extendido a N dimensiones. |
| **Convergencia** | Cuando el algoritmo "se estabiliza": los centroides dejan de moverse más allá de un umbral pequeño. |
| **Aprendizaje no supervisado** | Aprender patrones en datos sin usar etiquetas ni respuestas correctas. K-Means es un ejemplo. |
| **Iteración** | Una vuelta completa del bucle K-Means: asignar todos los puntos + recalcular todos los centroides. |
| **THRESHOLD** | El umbral de convergencia (0.001 en nuestro código). Si ningún centroide se mueve más de este valor, paramos. |
| **K** | El número de clusters que queremos encontrar. En K-Means hay que especificarlo de antemano. |

### Conceptos de sistemas y paralelismo

| Término | Definición |
|---------|------------|
| **CPU / Procesador** | El componente del computador que ejecuta instrucciones del programa. |
| **Core / Núcleo** | Unidad de procesamiento dentro de la CPU. Una CPU de 4 cores puede ejecutar 4 instrucciones simultáneamente. |
| **Proceso** | Un programa en ejecución con su propio espacio de memoria. |
| **Hilo (Thread)** | Una sub-tarea dentro de un proceso. Los hilos de un mismo proceso comparten la memoria. |
| **Paralelismo** | Ejecutar varias tareas simultáneamente en diferentes cores de la CPU. |
| **Race condition** | Error de concurrencia: dos hilos leen y modifican el mismo dato simultáneamente, produciendo resultados incorrectos. |
| **Overhead** | Trabajo extra que no es parte del algoritmo principal. En OpenMP: crear hilos, sincronizarlos, distribuir trabajo. |
| **Speedup** | Cuántas veces más rápido es el programa con N hilos vs. con 1 hilo. Fórmula: S = T₁ / Tₙ. |
| **Ley de Amdahl** | El speedup máximo posible está limitado por la fracción serial del programa. Fórmula: S(N) = 1/((1-P) + P/N). |
| **Embarrassingly parallel** | Un problema donde las tareas son totalmente independientes y no necesitan comunicarse. La fase de asignación es un ejemplo. |
| **Barrera (barrier)** | Punto de sincronización donde todos los hilos esperan hasta que el último termina su trabajo. OpenMP añade una al final de cada `#pragma omp parallel for`. |

### Conceptos de OpenMP

| Término | Definición |
|---------|------------|
| **OpenMP** | Estándar para programación paralela en C/C++/Fortran usando directivas `#pragma`. |
| **`#pragma`** | Instrucción al compilador que no modifica la lógica del programa pero le dice cómo generarlo. |
| **`parallel for`** | Directiva que divide las iteraciones de un bucle entre los hilos disponibles. |
| **`schedule(static)`** | Tipo de reparto: divide las iteraciones en bloques iguales asignados antes de empezar. |
| **`omp_set_num_threads(n)`** | Configura cuántos hilos usará la próxima zona paralela. |
| **`omp_get_wtime()`** | Retorna el tiempo actual en segundos. Útil para medir duración de código. |
| **Modelo fork-join** | El programa corre serial → crea hilos (fork) → trabajan en paralelo → se reúnen (join) → continúa serial. |
| **Variable privada** | Cada hilo tiene su propia copia. Las variables declaradas dentro del bucle son privadas automáticamente. |
| **Variable compartida** | Todos los hilos comparten la misma posición de memoria. Las variables globales son compartidas por defecto. |

### Conceptos de C

| Término | Definición |
|---------|------------|
| **`#include`** | Instrucción para incluir una biblioteca de funciones en el programa. |
| **`#define`** | Define una constante simbólica. El compilador reemplaza el nombre por el valor antes de compilar. |
| **`double`** | Tipo de dato: número decimal de 64 bits (alta precisión). Ejemplos: 5.1, 3.14159. |
| **`int`** | Tipo de dato: número entero de 32 bits. Ejemplos: 0, 1, 150, -5. |
| **Array 2D** | Tabla de datos accesible con dos índices: `data[fila][columna]`. |
| **`fopen` / `fclose`** | Abrir y cerrar archivos. Siempre hay que cerrar lo que se abre. |
| **`fscanf`** | Leer datos formateados de un archivo (análogo a `scanf` pero de archivo). |
| **`fprintf`** | Escribir datos formateados en un archivo (análogo a `printf` pero a archivo). |
| **`fabs(x)`** | Valor absoluto de x en punto flotante. Requiere `#include <math.h>`. |
| **`exit(1)`** | Termina el programa inmediatamente con código de error. |
| **Variable global** | Declarada fuera de toda función. Existe durante toda la ejecución del programa. |
| **Variable local** | Declarada dentro de una función. Existe solo mientras esa función ejecuta. |
| **Puntero (`*`)** | Variable que guarda la dirección de memoria de otra variable, no el valor directamente. |
| **`&variable`** | Operador de dirección: "dame la dirección de memoria de esta variable". |

---

*Proyecto K-Means con OpenMP — Dataset Iris*  
*Guía de estudio para principiantes*
