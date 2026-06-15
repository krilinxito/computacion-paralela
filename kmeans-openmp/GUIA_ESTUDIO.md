# Guía de Estudio: K-Means Paralelo con OpenMP
### Explicado desde cero, para alguien que nunca lo ha visto

**Proyecto:** K-Means Clustering paralelizado con OpenMP — Dataset sintético masivo  
**Nivel:** Sin conocimiento previo de K-Means ni de paralelismo

---

## Tabla de Contenidos

- [Parte 0 — Los datos: de qué estamos hablando](#parte-0--los-datos-de-qué-estamos-hablando)
- [Parte 0.5 — Por qué Iris no basta: el problema de escala](#parte-05--por-qué-iris-no-basta-el-problema-de-escala)
- [Parte 1 — K-Means desde cero](#parte-1--k-means-desde-cero)
- [Parte 2 — Cómputo paralelo desde cero](#parte-2--cómputo-paralelo-desde-cero)
- [Parte 3 — OpenMP desde cero](#parte-3--openmp-desde-cero)
- [Parte 4 — El código C explicado línea a línea](#parte-4--el-código-c-explicado-línea-a-línea)
- [Parte 5 — ¿Por qué solo paralelizamos la asignación?](#parte-5--por-qué-solo-paralelizamos-la-asignación)
- [Parte 6 — El speedup: por qué ahora sí funciona](#parte-6--el-speedup-por-qué-ahora-sí-funciona)
- [Parte 7 — Compilar y ejecutar](#parte-7--compilar-y-ejecutar)
- [Parte 8 — Gráficas con Gnuplot](#parte-8--gráficas-con-gnuplot)
- [Parte 9 — Glosario completo](#parte-9--glosario-completo)

---

## Parte 0 — Los datos: de qué estamos hablando

### ¿Qué son "datos"?

Antes de hablar de algoritmos, hay que entender con qué trabajamos.

Un **dato** es simplemente información que podemos medir y escribir como número. Por ejemplo, si quisieras describir una flor con números, podrías medir el largo y ancho de su sépalo y de su pétalo: 4 medidas. Ese es el famoso **dataset Iris**: 150 flores, cada una con 4 medidas. Es el "hola mundo" del clustering y lo usaremos como **analogía didáctica** a lo largo de esta guía porque es fácil de imaginar.

Pero hay un problema: **150 flores es ridículamente poco** para medir si la paralelización sirve (lo explicamos en la [Parte 0.5](#parte-05--por-qué-iris-no-basta-el-problema-de-escala)). Por eso este proyecto **no carga Iris desde un archivo**: en su lugar **genera datos sintéticos** — puntos artificiales creados por el propio programa — y puede generar tantos como queramos: 10 mil, 100 mil, ¡hasta 1 millón!

### ¿Qué es un punto y qué es una dimensión?

Olvidemos las flores por un momento y pensemos en general.

- Un **punto** es una cosa que queremos agrupar (una flor, un cliente, un píxel...).
- Una **dimensión** es una medida de ese punto (largo del pétalo, edad del cliente, color del píxel...).

Si un punto tiene 2 dimensiones, lo podemos dibujar en una hoja de papel (eje X, eje Y). Si tiene 3, en un cubo. Si tiene 16 dimensiones — como en este proyecto — ya no lo podemos dibujar, pero la matemática funciona exactamente igual: solo son 16 números por punto.

En este proyecto cada punto tiene:
- **D = 16 dimensiones** (16 medidas numéricas por punto)
- Pertenece a uno de **K = 10 grupos** (clusters)

### Cómo se generan los datos sintéticos

En lugar de leer un CSV, el programa **fabrica** los puntos con la función `generate_dataset()`. La idea es sencilla:

1. Queremos K = 10 grupos. Le damos a cada grupo un "centro" diferente: el grupo 0 vive cerca de la coordenada `(0, 0, ..., 0)`, el grupo 1 cerca de `(10, 10, ..., 10)`, el grupo 2 cerca de `(20, 20, ..., 20)`, etc.
2. Para cada punto, elegimos a qué grupo pertenece y lo colocamos **cerca** del centro de ese grupo, pero no exactamente encima: le sumamos un poco de **ruido aleatorio** (una pequeña desviación al azar).

El resultado son 10 "nubes" de puntos bien separadas. Justo lo que K-Means debe redescubrir.

### ¿Qué es el "ruido gaussiano" y Box-Muller?

El ruido que sumamos sigue una **distribución gaussiana** (la famosa "campana de Gauss"): valores cercanos a cero son muy probables, valores grandes son raros. Es el ruido más natural: si midieras la misma flor 1000 veces con una regla, tus errores formarían una campana.

El problema es que la función `rand()` de C solo da números **uniformes** (todos igual de probables entre 0 y 1), no gaussianos. Para convertir uniformes en gaussianos usamos la **transformación de Box-Muller**, una fórmula matemática que toma dos números uniformes `u1, u2` y produce un número gaussiano `g`:

```
g = √(-2 · ln(u1)) · cos(2π · u2)
```

No hace falta entender la fórmula a fondo; basta saber que **convierte aleatoriedad uniforme en aleatoriedad con forma de campana**, que es lo que necesitamos para que cada grupo sea una nube realista.

### Semilla fija: experimentos reproducibles

El programa llama a `srand(42)` antes de generar los datos. La **semilla** (`42`) hace que la secuencia de números "aleatorios" sea siempre la misma. Así, cada vez que ejecutas el programa obtienes **exactamente los mismos datos**, y las comparaciones de tiempo entre 1, 2 y 4 hilos son justas (todos procesan los mismos puntos).

### Lo importante: K-Means IGNORA a qué grupo pertenece cada punto

Aunque nosotros **sabemos** a qué grupo pusimos cada punto al generarlo, K-Means **no lo usa**. Le damos solo los 16 números de cada punto y le preguntamos: *¿puedes encontrar las 10 nubes tú solo?*

Eso se llama **aprendizaje no supervisado**: aprender sin que alguien te dé las respuestas correctas. (Con Iris sería lo mismo: K-Means no ve la columna "especie".)

---

## Parte 0.5 — Por qué Iris no basta: el problema de escala

### El enunciado del proyecto

El proyecto pedía: *"Algoritmo de K-Means Clustering. Optimización de la fase de asignación de centroides para conjuntos de datos masivos."*

La palabra clave es **masivos**. El dataset Iris tiene 150 puntos. Eso **no es masivo** — es minúsculo. Y resulta que con datos minúsculos la paralelización no solo no ayuda: **hace el programa más lento**.

### La intuición: contratar albañiles para un solo ladrillo

Imagina que tienes que poner un ladrillo. Lo haces tú solo en 2 segundos.

Ahora imagina que, para "ir más rápido", llamas a 4 albañiles. Tienes que:
1. Llamarlos por teléfono (tardan en llegar)
2. Explicarles la tarea
3. Repartir... ¿un solo ladrillo entre 4 personas?
4. Esperar a que todos digan "listo"
5. Despedirlos

Todo ese trámite tarda **mucho más** que los 2 segundos de poner el ladrillo tú mismo. Has hecho el trabajo **más lento** por usar más trabajadores.

Eso es exactamente lo que pasa con N = 150: el trabajo real (calcular distancias de 150 puntos) tarda menos de lo que cuesta **crear, coordinar y destruir** los hilos.

### Los números del problema original (Iris, N=150)

Con el proyecto original en Iris, las mediciones daban algo así:

```
1 hilo  → speedup = 1.00  (referencia)
2 hilos → speedup ≈ 0.21  (¡5 veces MÁS LENTO!)
4 hilos → speedup ≈ 0.15  (¡6.7 veces MÁS LENTO!)
```

El `speedup` es "cuántas veces más rápido". Un speedup **menor a 1** significa que empeoramos. Con Iris, paralelizar es contraproducente.

### La solución: datos masivos sintéticos

Para que la fase de asignación realmente domine el tiempo —y la paralelización valga la pena— necesitamos **muchos** puntos. Por eso generamos datos sintéticos con N de 10.000 a 1.000.000.

Con N grande, el cálculo de distancias se vuelve tan pesado que el trámite de los hilos se vuelve insignificante en comparación. Ahí sí, repartir el trabajo entre cores hace el programa **realmente más rápido**:

```
N = 1.000.000:
1 hilo  → speedup = 1.00
2 hilos → speedup ≈ 1.74
4 hilos → speedup ≈ 2.85   ¡Casi 3 veces más rápido!
```

### El concepto clave: el umbral de paralelización

Existe un **tamaño mínimo de problema** por debajo del cual paralelizar no compensa, y por encima del cual sí. A ese punto de equilibrio lo llamamos **umbral de paralelización**.

- **Por debajo del umbral** (datos pequeños): el *overhead* de los hilos domina → speedup < 1 → la paralelización empeora.
- **Por encima del umbral** (datos masivos): el cómputo domina → speedup > 1 → la paralelización ayuda.

Toda esta guía gira en torno a **demostrar empíricamente dónde está ese umbral** y por qué existe. La Ley de Amdahl (Parte 6) lo formaliza con matemáticas.

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

**En D dimensiones (en nuestro código, D=16):**

Con D características, el principio es el mismo, solo se suman D términos:

```
distancia = √( (x₁-y₁)² + (x₂-y₂)² + ... + (x_D-y_D)² )
```

No podemos visualizar 16 dimensiones, pero la matemática funciona igual. En este proyecto cada punto tiene 16 dimensiones, así que sumamos 16 términos. (Con el dataset Iris serían 4.)

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

**¿Convergió?** El centroide C1 se movió de (1.0, 0.2) a (1.1, 0.233). Desplazamiento = 0.1. Como 0.1 > 0.0001 (nuestro umbral, `THRESHOLD`), continuamos.

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
- El proceso tiene la tabla `data` (hasta 1.000.000 × 16 números) en memoria
- Cuando creamos 4 hilos, los 4 pueden leer esa tabla sin copiarla
- Eso es eficiente: no gastamos tiempo copiando millones de puntos para cada hilo

### 2.4 ¿Qué es paralelismo?

**Paralelismo** es hacer varias cosas al mismo tiempo, en cores distintos de la CPU.

Supongamos N = 1.000.000 puntos.

Sin paralelismo (serial):
```
Core 0: Punto 1 → Punto 2 → Punto 3 → ... → Punto 1.000.000
Core 1: inactivo
Core 2: inactivo
Core 3: inactivo
```

Con paralelismo (4 hilos):
```
Core 0: Puntos 1        → ... → 250.000
Core 1: Puntos 250.001  → ... → 500.000
Core 2: Puntos 500.001  → ... → 750.000
Core 3: Puntos 750.001  → ... → 1.000.000
```

Si cada punto tarda el mismo tiempo, el resultado es ~4 veces más rápido. En teoría (en la práctica veremos ~2.8x, por las razones de la Parte 6).

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

Esto explica el **umbral de paralelización** de la Parte 0.5: con datos pequeños (como Iris, N=150) el overhead hace el programa **más lento**. Solo cuando el trabajo es masivo (N de cientos de miles o millones) el cómputo supera al overhead y la paralelización **acelera** de verdad. Por eso este proyecto genera datos sintéticos a gran escala.

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
for (int i = 0; i < n; i++) {
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

Con 4 hilos y 1.000.000 de iteraciones:
```
Hilo 0 → iteraciones i = 0          ... 249.999
Hilo 1 → iteraciones i = 250.000    ... 499.999
Hilo 2 → iteraciones i = 500.000    ... 749.999
Hilo 3 → iteraciones i = 750.000    ... 999.999
```

Cada hilo trabaja en su rango sin pisarse con los demás.

**¿Por qué `schedule(static)` y no otras opciones?**

| schedule | Comportamiento | Cuándo usarlo |
|----------|---------------|---------------|
| static | Bloques iguales, asignados antes | Cuando todas las iteraciones tardan igual |
| dynamic | Iteraciones asignadas "al vuelo" | Cuando algunas tardan más que otras |
| guided | Bloques grandes al inicio, pequeños al final | Mezcla balanceada |

Para K-Means, cada punto tarda exactamente el mismo tiempo en procesar (siempre se compara contra los mismos K centroides), así que `static` es perfecto y además evita el pequeño overhead de repartir trabajo "al vuelo" de `dynamic`.

### 3.6 Variables privadas vs compartidas

En una zona paralela, las variables pueden ser:

**Compartidas (shared):** Todos los hilos ven el mismo espacio de memoria. Por defecto, las variables declaradas **fuera** del `#pragma` son compartidas.

```c
// FUERA del #pragma → compartidas → todos los hilos las ven
double *data;          // todos los hilos leen esto (n × 16 números)
int    *assignments;   // cada hilo escribe en su propio índice
double *centroids;     // todos los hilos leen esto (10 × 16 números)
```

**Privadas (private):** Cada hilo tiene su propia copia. Las variables declaradas **dentro** del bucle son privadas automáticamente.

```c
#pragma omp parallel for schedule(static)
for (int i = 0; i < n; i++) {
    double min_dist = 1e300;  // ← PRIVADA: cada hilo tiene su propia min_dist
    int    best     = 0;      // ← PRIVADA: cada hilo tiene su propio best
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
gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm
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
#include <stdio.h>    // entrada/salida: printf, fprintf, fopen, fclose
#include <stdlib.h>   // memoria y sistema: malloc, free, rand, srand, exit
#include <string.h>   // memcpy, memset (copiar y poner a cero bloques de memoria)
#include <math.h>     // funciones matemáticas: sqrt, log, cos, fabs, M_PI
#include <omp.h>      // OpenMP: omp_get_wtime(), omp_set_num_threads()
```

Un `#include` es como decirle al compilador: "necesito usar funciones de esta biblioteca". Como ahora generamos datos (en vez de leerlos de un archivo) y reservamos memoria dinámica, aparecen `<stdlib.h>` (para `malloc`/`rand`) y `<string.h>` (para `memcpy`/`memset`).

### 4.2 Las constantes (`#define`)

```c
#define N_DIMS    16    // cuántas dimensiones (medidas) tiene cada punto
#define K         10    // cuántos clusters queremos
#define MAX_ITER  100   // máximo de vueltas del algoritmo
#define THRESHOLD 1e-4  // "suficientemente quieto" para parar
```

`#define` crea un **sustituto de texto**: antes de compilar, el compilador reemplaza cada `N_DIMS` por `16` en todo el código. No es una variable; es como hacer "buscar y reemplazar" en el código fuente.

**Diferencia importante con la versión original (Iris):**

Ya **no existe** `N_POINTS`. ¿Por qué? Porque el número de puntos ya no es fijo: el programa prueba varios tamaños (10.000, 100.000, 500.000, 1.000.000) en la misma ejecución. El tamaño actual se pasa como una **variable** `n` a las funciones, no como una constante de compilación.

Subimos `N_DIMS` de 4 a **16** y `K` de 3 a **10**. ¿Por qué? Para que cada distancia sea más costosa de calcular (16 restas y multiplicaciones contra 10 centroides = 160 operaciones por punto), de modo que la fase de asignación domine claramente el tiempo y la paralelización tenga algo sustancial que repartir.

### 4.3 Memoria dinámica: punteros en el heap

En la versión Iris, los datos vivían en arrays globales de tamaño fijo: `double data[150][4]`. Eso funciona para 150 puntos (4.800 bytes), pero **no** para 1.000.000 de puntos en 16 dimensiones:

```
1.000.000 × 16 × 8 bytes = 128.000.000 bytes ≈ 128 MB
```

128 MB no caben en la memoria estática ni en la pila (*stack*) — reventarían el programa. La solución es pedir esa memoria en tiempo de ejecución con `malloc`, que la reserva en una zona llamada **heap** (montón), pensada para bloques grandes:

```c
double *data        = malloc((size_t)n * N_DIMS * sizeof(double));
int    *assignments = malloc((size_t)n * sizeof(int));
double *centroids   = malloc((size_t)K * N_DIMS * sizeof(double));
double *init_c      = malloc((size_t)K * N_DIMS * sizeof(double));
```

**¿Qué es un puntero?** Una variable que guarda la **dirección de memoria** donde empieza un bloque. `double *data` significa "data apunta al primer `double` de un bloque de doubles".

**¿Qué hace `malloc`?** Reserva un bloque de bytes y devuelve su dirección. `malloc(n * N_DIMS * sizeof(double))` pide espacio para `n × 16` números decimales. `sizeof(double)` es 8 (bytes por double). El `(size_t)` evita que la multiplicación se desborde con números grandes.

**Arreglo "aplanado" (1D que simula 2D):**

Como `data` es un puntero plano, no podemos escribir `data[i][d]`. En su lugar calculamos el índice a mano:

```c
data[i * N_DIMS + d]   // punto i, dimensión d
```

La fila `i` empieza en la posición `i * 16`, y dentro de ella la columna `d` está en `+ d`. Es exactamente la misma tabla, solo que escrita como una larga tira de números:

```
[ punto0: d0 d1 ... d15 | punto1: d0 d1 ... d15 | punto2: ... ]
   índices 0..15            índices 16..31
```

**Liberar la memoria:** Todo lo que se pide con `malloc` se debe devolver con `free` cuando ya no se usa, o el programa "gotea" memoria (*memory leak*). Al final de cada tamaño `n`:

```c
free(data); free(assignments); free(centroids); free(init_c);
```

### 4.4 La función `generate_dataset` — fabricar los puntos

En lugar de `load_iris`, ahora **generamos** los datos:

```c
static void generate_dataset(int n, double *data) {
    srand(42);
    for (int i = 0; i < n; i++) {
        int cluster = i % K;
        for (int d = 0; d < N_DIMS; d++) {
            double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
            double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
            double g  = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
            data[i * N_DIMS + d] = cluster * 10.0 + g;
        }
    }
}
```

Línea por línea:

- `srand(42)`: fija la **semilla** del generador aleatorio. Garantiza que cada ejecución produzca los mismos puntos (reproducibilidad).
- `for (int i = 0; i < n; i++)`: recorre los `n` puntos.
- `int cluster = i % K;`: el operador `%` es el **residuo**. `i % 10` da 0,1,2,...,9,0,1,2,... Así asignamos cada punto a uno de los K=10 grupos de forma cíclica (puntos balanceados entre clusters).
- `u1`, `u2`: dos números **uniformes** en el rango (0, 1]. El `+1.0` y `(RAND_MAX+1.0)` evitan obtener exactamente 0 (que rompería el `log`).
- `g = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2)`: la **transformación de Box-Muller**, que convierte `u1, u2` (uniformes) en `g` (gaussiano con media 0 y desviación 1). `M_PI` es π, definido en `<math.h>`.
- `data[i*N_DIMS + d] = cluster * 10.0 + g;`: coloca el punto. El centro del grupo es `cluster * 10` (grupo 0 → 0, grupo 1 → 10, ...) y le sumamos el ruido `g`. Resultado: una nube gaussiana alrededor del centro del grupo.

`static` antes de la función significa "esta función solo se usa dentro de este archivo" — buena práctica para funciones auxiliares.

### 4.5 La función `run_kmeans` — el corazón del programa

Ahora recibe **parámetros** (en vez de usar globales): el tamaño `n`, el número de hilos y punteros a los arreglos.

**Bloque 1: Preparación**

```c
static double run_kmeans(int n, int num_threads,
                         const double *data, int *assignments,
                         double *centroids, const double *init_c) {
    // Restaurar centroides iniciales para comparación justa entre runs
    memcpy(centroids, init_c, (size_t)K * N_DIMS * sizeof(double));

    omp_set_num_threads(num_threads);
    double t_start = omp_get_wtime();
```

- `const double *data`: el `const` promete que la función **no modificará** los datos, solo los leerá. Es una garantía de seguridad y documenta la intención.
- `memcpy(centroids, init_c, ...)`: copia los centroides iniciales de golpe (bloque de bytes). Restauramos el estado inicial para que las corridas con 1, 2 y 4 hilos partan del **mismo punto** y la comparación de tiempos sea justa. Si no lo hiciéramos, la segunda corrida empezaría con centroides ya convergidos y terminaría antes, falseando el speedup.
- `omp_set_num_threads(num_threads)`: configura cuántos hilos usará la próxima zona paralela.
- `omp_get_wtime()`: arranca el cronómetro.

**Bloque 2: Fase de asignación (paralela)**

```c
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
```

- `#pragma omp parallel for schedule(static)`: reparte el bucle de `n` puntos entre los hilos (ver Parte 3.5).
- `1e300`: un número gigantesco usado como "infinito" inicial, para que cualquier distancia real sea menor. (Subimos de `1e18` a `1e300` por seguridad; ambos funcionan.)
- `min_dist` y `best`: locales al bucle → **privadas** por hilo → sin race condition.
- El doble bucle `c`/`d` calcula la distancia al cuadrado del punto `i` a cada centroide `c`. Fíjate en la aritmética de índices aplanados: `data[i*N_DIMS + d]` y `centroids[c*N_DIMS + d]`.
- `assignments[i] = best`: cada hilo escribe solo en su propio índice `i` → sin conflicto.

**Bloque 3: Fase de recálculo (serial)**

```c
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
```

- `new_c[K * N_DIMS]` y `counts[K]`: acumuladores. Como K=10 y N_DIMS=16 son **pequeños y fijos**, sí caben en la pila (no necesitan `malloc`). `new_c` tiene 160 doubles; `counts`, 10 enteros.
- `memset(..., 0, sizeof(...))`: pone todo el bloque a cero de una sola vez (más limpio que un bucle).
- El bucle suma las coordenadas de cada punto al acumulador de su cluster y cuenta cuántos puntos tiene cada cluster.

```c
        double max_shift = 0.0;
        for (int c = 0; c < K; c++) {
            if (counts[c] == 0) continue;
            for (int d = 0; d < N_DIMS; d++) {
                double new_val = new_c[c * N_DIMS + d] / counts[c];
                double shift   = fabs(new_val - centroids[c * N_DIMS + d]);
                if (shift > max_shift) max_shift = shift;
                centroids[c * N_DIMS + d] = new_val;
            }
        }
```

Para cada cluster, el nuevo centroide es el **promedio** (suma / conteo). Medimos cuánto se movió cada coordenada con `fabs()` (valor absoluto) y guardamos el mayor desplazamiento en `max_shift`.

`if (counts[c] == 0) continue;`: si un cluster quedó vacío, saltamos para no dividir por cero.

**Bloque 4: Convergencia y fin**

```c
        if (max_shift < THRESHOLD) break;
    }

    return omp_get_wtime() - t_start;
}
```

Si ningún centroide se movió más de `THRESHOLD` (1e-4), convergimos y salimos con `break`. Finalmente devolvemos el tiempo transcurrido (cronómetro de fin menos cronómetro de inicio).

### 4.6 La función `main` — el director de orquesta

Ahora `main` recorre **varios tamaños** y, para cada uno, prueba **varios números de hilos**.

```c
int main(void) {
    int sizes[]   = {10000, 100000, 500000, 1000000};
    int threads[] = {1, 2, 4};
    int n_sizes   = 4;
    int n_threads = 3;
```

`sizes[]` son los tamaños de dataset a probar (el "sweep" o barrido). `threads[]` son las configuraciones de hilos.

```c
    FILE *fout = fopen("resultados/benchmark.dat", "w");
    if (!fout) {
        fprintf(stderr, "Error: no se pudo crear resultados/benchmark.dat\n");
        return 1;
    }
    fprintf(fout, "# N  hilos  tiempo_s  speedup\n");
```

Abrimos el archivo de salida `benchmark.dat` y escribimos una línea de cabecera (la `#` hace que gnuplot la ignore como comentario).

```c
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
```

Para cada tamaño `n`: reservamos memoria (ver 4.3) y verificamos que `malloc` no devolvió `NULL` (señal de que no hubo memoria).

```c
        generate_dataset(n, data);
        memcpy(init_c, data, (size_t)K * N_DIMS * sizeof(double));
```

Generamos los `n` puntos y tomamos los **primeros K puntos** como centroides iniciales (copiándolos a `init_c`). Es una inicialización simple y reproducible: como los puntos se generan en orden de cluster (0,1,2,...), los primeros K caen cada uno en un grupo distinto.

```c
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
```

Para este `n`, corremos K-Means con 1, 2 y 4 hilos:
- `t1` guarda el tiempo con 1 hilo (la referencia para el speedup).
- `speedup = t1 / t`: cuántas veces más rápido que con 1 hilo.
- Imprimimos en pantalla **y** escribimos la línea en `benchmark.dat`.

`%-10d` significa entero alineado a la **izquierda** en 10 caracteres (el `-` es la alineación izquierda).

```c
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
```

Liberamos la memoria de este tamaño antes de pasar al siguiente (para no acumular cientos de MB). Al terminar todos los tamaños, cerramos el archivo y avisamos.

`return 0`: convención en C para indicar que el programa terminó exitosamente.

**El formato de `benchmark.dat`** queda así (4 columnas: N, hilos, tiempo, speedup):

```
# N  hilos  tiempo_s  speedup
10000 1 0.002402 1.0000
10000 2 0.001337 1.7960
10000 4 0.000921 2.6080
100000 1 0.013071 1.0000
...
```

---

## Parte 5 — ¿Por qué solo paralelizamos la asignación?

### 5.1 El análisis de operaciones

Contemos cuántas operaciones hace cada fase en cada iteración del K-Means, con nuestros parámetros (K=10, D=16) y el caso grande N=1.000.000:

**Fase de asignación:** Para cada uno de los N puntos, calculamos la distancia a cada uno de los K=10 centroides, y cada distancia requiere D=16 multiplicaciones y sumas:

```
Operaciones = N × K × D = 1.000.000 × 10 × 16 = 160.000.000 ops
```

**Fase de recálculo:** Para cada uno de los N puntos, sumamos sus D=16 coordenadas al acumulador de su cluster:

```
Operaciones = N × D = 1.000.000 × 16 = 16.000.000 sumas
```

La asignación hace **10 veces más trabajo** que el recálculo (un factor exacto de K). Es el cuello de botella absoluto: paralelizarla es donde está toda la ganancia. Por eso el enunciado del proyecto pedía específicamente *"optimización de la fase de asignación de centroides"*.

### 5.2 La condición de independencia

Para paralelizar un bucle con OpenMP sin errores, las iteraciones deben ser **independientes**: la iteración `i` no puede depender del resultado de la iteración `j`.

**En la fase de asignación:**

```c
// Iteración i:
//   LEE: data[i*N_DIMS + ...] (solo su propio punto)
//   LEE: centroids[...] (todos leen, nadie escribe)
//   ESCRIBE: assignments[i] (su propio índice, nadie más escribe aquí)
```

Cada iteración trabaja en un índice diferente de `assignments`. No hay conflicto. Es **seguro paralelizar**.

**En la fase de recálculo:**

```c
// Iteración i (sin paralelismo):
int c = assignments[i];
counts[c]++;                          // ← PROBLEMA
new_c[c*N_DIMS + d] += data[i*N_DIMS + d]; // ← PROBLEMA
```

Si dos puntos i=5 e i=42 pertenecen al mismo cluster c=0, ambas iteraciones intentarían modificar `counts[0]` y `new_c[0*N_DIMS + d]` al mismo tiempo. Race condition garantizada.

**¿Se podría paralizar con `reduction`?**

Sí, OpenMP tiene la directiva `reduction` para sumar acumuladores de forma segura:
```c
// Esta versión sería correcta pero más compleja:
#pragma omp parallel for reduction(+:counts) schedule(static)
for (int i = 0; i < n; i++) { ... }
```

Sin embargo, el recálculo hace ~16M operaciones vs ~160M de la asignación (con N=1M): es solo el **~9%** del trabajo, mientras el **~91%** ya está paralelizado. Añadir complejidad (y overhead) para paralizar esa pequeña fracción no vale la pena. La solución simple (serial) es la correcta aquí.

---

## Parte 6 — El speedup: por qué ahora sí funciona

### 6.1 Los resultados reales

Nuestras mediciones (guardadas en `benchmark.dat`, en una máquina de 16 cores) muestran **speedup mayor a 1** para todos los tamaños — la paralelización funciona:

```
                 2 hilos    4 hilos
N = 10.000    →   1.80x  →   2.61x
N = 100.000   →   1.91x  →   2.76x
N = 500.000   →   1.31x  →   2.89x
N = 1.000.000 →   1.74x  →   2.85x
```

Con 4 hilos el programa corre consistentemente **~2.6 a 2.9 veces más rápido**. Compáralo con el desastre del dataset Iris original (Parte 0.5), donde 4 hilos daba speedup 0.15 (¡6.7x más lento!). La diferencia es **el tamaño del problema**.

(Los valores exactos varían según la máquina, su número de cores, su carga y la velocidad de la CPU. Lo importante es la tendencia: speedup > 1, creciente con los hilos.)

### 6.2 Anatomía del tiempo de ejecución (N grande)

Cuando el programa corre con 4 hilos sobre N=1.000.000, el reparto del tiempo es muy distinto al de Iris:

| Componente | Tiempo aproximado |
|------------|-----------------|
| Crear 4 hilos | ~20-40 μs |
| Distribuir trabajo (1M puntos) | ~1-5 μs |
| **Calcular distancias de 1M puntos** | **~130.000 μs (130 ms)** |
| Esperar que todos terminen (barrera) | ~5-10 μs |
| Destruir hilos | ~10-20 μs |
| **Total overhead** | ~36-77 μs |
| **Trabajo real** | ~130.000 μs |

Ahora el **trabajo real es ~2000 veces mayor que el overhead**. El trámite de los hilos (que con Iris arruinaba todo) se vuelve insignificante. Es como contratar 4 albañiles para construir una casa entera: el rato que tardan en organizarse no importa frente a las semanas de trabajo.

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

**Con nuestros datos (N grande):**

La fase paralela (asignación) hace ~160M ops; la serial (recálculo + overhead) hace ~16M ops + un overhead despreciable. La fracción paralela es aproximadamente:

```
P ≈ 160 / (160 + 16) ≈ 0.91
```

Aplicando Amdahl:

```
S(2) = 1 / ( (1-0.91) + 0.91/2 ) = 1 / (0.09 + 0.455) = 1 / 0.545 ≈ 1.83
S(4) = 1 / ( (1-0.91) + 0.91/4 ) = 1 / (0.09 + 0.2275) = 1 / 0.3175 ≈ 3.15
```

¡Estos valores teóricos (1.83 y 3.15) coinciden muy bien con los medidos (~1.8 y ~2.8)! La pequeña diferencia con 4 hilos (3.15 teórico vs ~2.8 real) se debe a factores no contemplados por Amdahl: contención de memoria (varios cores compitiendo por el bus de RAM), efectos de caché y el overhead que sí crece un poco con más hilos.

### 6.4 El umbral de paralelización, cuantificado

Amdahl explica el **umbral** de la Parte 0.5. La fracción paralela efectiva `P` depende del tamaño del problema:

| Dataset | Trabajo cómputo | Overhead relativo | P efectiva | Speedup 4 hilos |
|---------|----------------|-------------------|-----------|-----------------|
| Iris (N=150) | ~μs | domina | ≈ 0.01 | ≈ 0.15 (¡peor!) |
| N = 10.000 | ~ms | pequeño | ≈ 0.85 | ≈ 2.6 |
| N = 1.000.000 | ~130 ms | despreciable | ≈ 0.91 | ≈ 2.85 |

**Conclusión:** OpenMP es poderoso, pero solo cuando el trabajo supera al overhead. Con N=150 (Iris) el dataset es demasiado pequeño y paralelizar empeora. A partir de decenas de miles de puntos, el speedup se hace claramente positivo y se acerca (sin llegar) al número de hilos. Esa es exactamente la lección que el proyecto buscaba demostrar con "conjuntos de datos masivos".

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
gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm
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
gcc -O2 -fopenmp -o kmeans src/kmeans.c -lm
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
K-Means Paralelo con OpenMP
Dataset sintetico gaussiano: D=16 dimensiones, K=10 clusters

N          Hilos    Tiempo (s)     Speedup 
---------- -------- -------------- --------
  [generando N=10000 puntos...]
10000      1        0.002402       1.000   
10000      2        0.001337       1.796   
10000      4        0.000921       2.608   

  [generando N=100000 puntos...]
100000     1        0.013071       1.000   
100000     2        0.006859       1.906   
100000     4        0.004736       2.760   

  [generando N=500000 puntos...]
500000     1        0.066961       1.000   
500000     2        0.051024       1.312   
500000     4        0.023205       2.886   

  [generando N=1000000 puntos...]
1000000    1        0.133134       1.000   
1000000    2        0.076526       1.740   
1000000    4        0.046690       2.851   

Benchmark exportado a: resultados/benchmark.dat
Ejecutar 'gnuplot graficar.gp' para generar las graficas.
```

- **Tiempo con 1 hilo (N=1M):** ~0.133 segundos = 133 milisegundos
- **Tiempo con 4 hilos (N=1M):** ~0.047 segundos → mucho más rápido
- **Speedup = 2.851:** con 4 hilos el programa es ~2.85x la velocidad de 1 hilo

Fíjate cómo el speedup con 4 hilos crece y se mantiene en ~2.6-2.9x para todos los tamaños grandes. Los valores exactos varían según la máquina, su número de cores, carga del sistema y velocidad de la CPU.

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

Ambas se construyen a partir de `resultados/benchmark.dat` (4 columnas: N, hilos, tiempo, speedup).

**`speedup.png` — Speedup vs número de hilos (una curva por N)**

Muestra cómo escala el speedup según el tamaño del dataset:
- Eje X = número de hilos (1, 2, 4)
- Eje Y = speedup (T₁ / Tₙ)
- Línea gris punteada = speedup **ideal** (lineal: 2 hilos = 2x, 4 hilos = 4x)
- Una línea de color por cada N (10k, 100k, 500k, 1M)

Todas las curvas reales deben estar **por encima de 1** y subir con los hilos, acercándose (sin tocar) la línea ideal. Esto demuestra que con datos masivos la paralelización **sí acelera** — lo opuesto al caso Iris.

**`scaling.png` — Tiempo de ejecución vs tamaño N (escala log-log)**

Muestra cómo crece el tiempo al aumentar N y cuánto ayuda cada número de hilos:
- Eje X = N (escala logarítmica: 10⁴, 10⁵, 10⁶)
- Eje Y = tiempo en segundos (escala logarítmica)
- Una línea por número de hilos (1, 2, 4)

La línea de 4 hilos (verde) está siempre **por debajo** de la de 1 hilo (roja): tarda menos para el mismo N. Las tres líneas son aproximadamente rectas y paralelas, lo que en escala log-log confirma el crecimiento **lineal** O(N) del algoritmo.

> Nota: el filtrado por número de hilos en `scaling.png` se hace con `awk` vía pipe (`"< awk '$2==1' ..."`), porque las filas de un mismo hilo no son contiguas en el archivo y el truco habitual de gnuplot (`$2==1 ? $1 : 1/0`) rompería las líneas.

### 8.5 Ver las imágenes

```bash
# Linux con entorno gráfico:
eog resultados/speedup.png         # Eye of GNOME
xdg-open resultados/scaling.png    # abre con el visor predeterminado

# macOS:
open resultados/speedup.png

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
| **THRESHOLD** | El umbral de convergencia (1e-4 en nuestro código). Si ningún centroide se mueve más de este valor, paramos. |
| **K** | El número de clusters que queremos encontrar. En K-Means hay que especificarlo de antemano. (En este proyecto K=10.) |
| **Datos sintéticos** | Puntos generados artificialmente por el programa (no leídos de un archivo real), lo que permite elegir el tamaño N libremente. |
| **Distribución gaussiana** | La "campana de Gauss": ruido aleatorio donde los valores cercanos a la media son los más probables. |
| **Box-Muller** | Fórmula que convierte dos números aleatorios uniformes en uno gaussiano. La usamos para fabricar las nubes de puntos. |
| **Umbral de paralelización** | Tamaño mínimo de problema a partir del cual paralelizar acelera (speedup > 1). Por debajo, el overhead domina y empeora. |

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
| **Heap (montón)** | Zona de memoria para bloques grandes pedidos en tiempo de ejecución. Aquí viven nuestros datos de hasta 128 MB. |
| **`malloc` / `free`** | Reservar y liberar memoria en el heap. Todo lo que se pide con `malloc` se debe devolver con `free`. |
| **Arreglo aplanado** | Un puntero 1D que simula una tabla 2D: el elemento `[fila][col]` se accede como `arr[fila*ANCHO + col]`. |
| **`memcpy` / `memset`** | Copiar un bloque de bytes / poner un bloque de bytes a un valor (típicamente 0). Más rápido que un bucle. |
| **`rand` / `srand`** | Generador de números pseudoaleatorios y su semilla. `srand(42)` fija la semilla para reproducibilidad. |
| **`fopen` / `fclose`** | Abrir y cerrar archivos. Siempre hay que cerrar lo que se abre. |
| **`fprintf`** | Escribir datos formateados en un archivo (análogo a `printf` pero a archivo). |
| **`fabs(x)`** | Valor absoluto de x en punto flotante. Requiere `#include <math.h>`. |
| **`exit(1)`** | Termina el programa inmediatamente con código de error. |
| **Variable local** | Declarada dentro de una función. Existe solo mientras esa función ejecuta. |
| **Puntero (`*`)** | Variable que guarda la dirección de memoria de otra variable (o el inicio de un bloque), no el valor directamente. |
| **`const`** | Promesa de que una función no modificará el dato apuntado (solo lo lee). Documenta la intención y previene errores. |

---

*Proyecto K-Means con OpenMP — Dataset sintético masivo*  
*Guía de estudio para principiantes*
