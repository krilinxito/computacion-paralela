# Guía de Estudio: Criba de Eratóstenes Distribuida con MPI

> Versión **para dummies** — no se asume conocimiento previo de MPI, algoritmos ni parallelismo.
> Se explica todo desde lo más básico. Paso a paso. Con ejemplos numéricos reales.

---

## COMANDOS DE SETUP

Antes de empezar, asegúrate de tener instalado lo necesario.

### Ubuntu / Debian (y WSL2 con Ubuntu)
```bash
# OBLIGATORIO — para compilar y ejecutar el programa
sudo apt install libmpich-dev mpich

# OPCIONAL — para generar las gráficas PNG
sudo apt install gnuplot

# OPCIONAL — para compilar el informe IEEE en PDF
sudo apt install texlive-full
```

### Fedora / RHEL / Rocky Linux (este sistema)
```bash
# mpicc y mpirun ya están instalados en este sistema
# Solo falta gnuplot si quieres las gráficas:
sudo dnf install gnuplot

# pdflatex ya está instalado en este sistema
```

### macOS con Homebrew
```bash
brew install open-mpi          # OBLIGATORIO
brew install gnuplot           # OPCIONAL
brew install --cask mactex     # OPCIONAL (LaTeX)
```

### Windows con WSL2
```bash
# Abre una terminal Ubuntu en WSL2 y ejecuta:
sudo apt install libopenmpi-dev openmpi-bin   # OBLIGATORIO
sudo apt install gnuplot                       # OPCIONAL
sudo apt install texlive-full                  # OPCIONAL
```

---

## PARTE 1: Los números primos (empezamos desde cero)

### 1.1 ¿Qué es un número primo?

Un número primo es un número entero mayor que 1 que **solo puede dividirse exactamente
entre 1 y él mismo**. No tiene otros divisores.

```
¿Es primo el 7?
  7 ÷ 2 = 3.5   → no exacto
  7 ÷ 3 = 2.33  → no exacto
  7 ÷ 4 = 1.75  → no exacto
  7 ÷ 5 = 1.4   → no exacto
  7 ÷ 6 = 1.16  → no exacto
  Solo 7 ÷ 1 = 7 y 7 ÷ 7 = 1 son exactos.
  → 7 ES primo.

¿Es primo el 12?
  12 ÷ 2 = 6    → exacto! tiene un divisor además de 1 y 12
  → 12 NO es primo. Es compuesto (= 2 × 6 = 3 × 4).
```

Ejemplos clasificados:

```
PRIMOS:    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 ...
COMPUESTOS: 4=2×2,  6=2×3,  8=2×4,  9=3×3,  10=2×5,  12=3×4 ...
ESPECIAL:   1 → NO es primo (por definición matemática, el 1 se excluye)
```

Analogía: los primos son como los "átomos" de los números. Todo número compuesto se puede
construir multiplicando primos entre sí:
```
12 = 2 × 2 × 3
100 = 2 × 2 × 5 × 5
35 = 5 × 7
```

### 1.2 El problema: encontrar TODOS los primos hasta 10,000,000

Nuestro programa busca todos los números primos entre 2 y 10,000,000 (diez millones).

**¿Por qué no simplemente probar cada número?**

El método más obvio sería: para cada número X entre 2 y 10,000,000, probar si tiene
algún divisor revisando todos los números desde 2 hasta X-1.

El problema es que eso son demasiadas operaciones:
```
10,000,000 números × ~3,162 pruebas por número ≈ 31,620,000,000 operaciones
```
Eso tardaría minutos en una computadora moderna. La **Criba de Eratóstenes** hace el
mismo trabajo en milisegundos, usando una estrategia completamente diferente.

### 1.3 La Criba de Eratóstenes — ejemplo completo con N=30

Vamos a encontrar todos los primos hasta **30**. El procedimiento es idéntico al que
usa el programa para 10,000,000, solo con números más chicos para poder verlo bien.

---

**INICIO: Escribimos todos los candidatos. Ninguno está descartado todavía.**

```
 2  3  4  5  6  7  8  9 10
11 12 13 14 15 16 17 18 19 20
21 22 23 24 25 26 27 28 29 30
```

---

**PASO 1: El 2 es el primer número de la lista. Como nadie lo descartó, ES PRIMO.**

Ahora tachamos todos sus múltiplos: 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30.
(Los números tachados se muestran con ✗)

```
 2  3  ✗  5  ✗  7  ✗  9  ✗
11  ✗ 13  ✗ 15  ✗ 17  ✗ 19  ✗
21  ✗ 23  ✗ 25  ✗ 27  ✗ 29  ✗
```

*Razonamiento: si 4 = 2×2, entonces 4 NO puede ser primo. Lo mismo para 6, 8, etc.*

---

**PASO 2: El siguiente número SIN tachar es el 3. ES PRIMO.**

Tachamos los múltiplos de 3 que aún no están tachados: 9, 15, 21, 27.
(El 6, 12, 18, 24, 30 ya estaban tachados por el 2, no importa tacharlos de nuevo.)

```
 2  3  ✗  5  ✗  7  ✗  ✗  ✗
11  ✗ 13  ✗  ✗  ✗ 17  ✗ 19  ✗
 ✗  ✗ 23  ✗ 25  ✗  ✗  ✗ 29  ✗
```

---

**PASO 3: El siguiente sin tachar es el 5. ES PRIMO.**

Tachamos los múltiplos de 5 que quedan: solo el 25.
(El 10, 15, 20, 30 ya estaban tachados.)

```
 2  3  ✗  5  ✗  7  ✗  ✗  ✗
11  ✗ 13  ✗  ✗  ✗ 17  ✗ 19  ✗
 ✗  ✗ 23  ✗  ✗  ✗  ✗  ✗ 29  ✗
```

---

**PARAMOS AQUÍ.** El siguiente primo sin tachar es el 7, pero 7 > √30 ≈ 5.47.

¿Por qué paramos en √N? Lo explicamos a continuación.

---

**¿Por qué alcanza con llegar hasta √N?**

Regla: si un número compuesto C tiene algún factor, ese factor siempre es ≤ √C.

Ejemplo: C = 25.
- Sus factores son 5 × 5.
- 5 ≤ √25 = 5. ✓

Ejemplo: C = 21.
- Sus factores son 3 × 7.
- 3 ≤ √21 ≈ 4.58. ✓ (El factor más chico siempre es ≤ √C)

Entonces: si un compuesto C (≤ N) existe, ya fue tachado por algún primo p ≤ √N.
No hay ningún compuesto que "se escape" si trabajamos hasta √N.

Para N = 10,000,000: solo necesitamos los primos hasta √10,000,000 ≈ **3,162**.
Esos son exactamente 446 primos. El programa los llama "primos base".

---

**RESULTADO FINAL:** lo que no está tachado es primo.

```
 2  3  .  5  .  7  .  .  .
11  . 13  .  .  . 17  . 19  .
 .  . 23  .  .  .  .  . 29  .
```

Primos del 2 al 30: **2, 3, 5, 7, 11, 13, 17, 19, 23, 29** (10 primos).

---

## PARTE 2: ¿Qué es ejecutar en paralelo?

### 2.1 ¿Qué es un proceso?

Cuando abres un programa (Chrome, Word, etc.), el sistema operativo crea un **proceso**.
Un proceso es:

- Un programa que está corriendo en ese momento
- Tiene su **propia memoria privada** — sus propios datos
- Ningún otro programa puede leer o escribir en esa memoria

Imagina trabajadores en cubículos separados, cada uno con su propio escritorio cerrado:

```
┌───────────────┐   ┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│  Proceso 0    │   │  Proceso 1    │   │  Proceso 2    │   │  Proceso 3    │
│               │   │               │   │               │   │               │
│ RAM: [datos]  │   │ RAM: [datos]  │   │ RAM: [datos]  │   │ RAM: [datos]  │
│               │   │               │   │               │   │               │
│ (privado,     │   │ (privado,     │   │ (privado,     │   │ (privado,     │
│  solo él lo   │   │  solo él lo   │   │  solo él lo   │   │  solo él lo   │
│  puede ver)   │   │  puede ver)   │   │  puede ver)   │   │  puede ver)   │
└───────────────┘   └───────────────┘   └───────────────┘   └───────────────┘
```

Cuando usamos `mpirun -np 4 ./criba`, se crean 4 procesos separados en tu computador.
No necesitas 4 computadores — los 4 procesos se reparten entre los núcleos de tu CPU.

### 2.2 ¿Qué gana el paralelismo?

Si tienes que cavar una zanja de 4 kilómetros:

```
1 persona sola:
Tiempo → → → → → → → → → → → → → → → →
[============ km 0 a 4 ============]   4 horas

4 personas en paralelo:
Tiempo → → → → →
Persona 0: [= km 0 a 1 =]
Persona 1: [= km 1 a 2 =]   ← trabajan al mismo tiempo
Persona 2: [= km 2 a 3 =]
Persona 3: [= km 3 a 4 =]
                               ~1 hora
```

Con nuestro programa pasa lo mismo. En vez de kilómetros, son rangos de números:

```
1 proceso (secuencial):
[======== del 2 al 10,000,000 (todo) ========]   tiempo: ~0.026 s

4 procesos (paralelo):
[= 2 a 2.5M =]
[= 2.5M a 5M =]   ← los 4 se ejecutan al mismo tiempo
[= 5M a 7.5M =]
[= 7.5M a 10M =]   tiempo: ~0.010 s
```

El resultado es el mismo (664,579 primos), pero más rápido.

### 2.3 El problema: los procesos no comparten memoria

Aquí está el desafío de la paralelización.

El Proceso 0 necesita calcular los "primos base" (los 446 primos hasta 3,162) antes
de poder trabajar. Pero los Procesos 1, 2 y 3 también los necesitan para tamizan
sus propios segmentos.

Como cada proceso tiene **su propia memoria privada**, el Proceso 1 no puede
simplemente "mirar" los datos que calculó el Proceso 0.

Para compartir información, deben **enviarse mensajes explícitamente**. Eso es
exactamente lo que hace **MPI**.

---

## PARTE 3: MPI desde cero

### 3.1 ¿Qué es MPI?

**MPI** (Message Passing Interface) es un estándar que define cómo los procesos se
comunican **enviando y recibiendo mensajes**.

Como el servicio postal entre oficinas separadas:
- Si el jefe (Proceso 0) necesita darle datos al empleado (Proceso 1), le manda un mensaje
- El empleado recibe el mensaje y ahora tiene esos datos en su propia memoria
- Nadie puede "entrar a la oficina del jefe" — solo se comunican por mensajes

MPI puede funcionar:
- En **1 computador** con varios núcleos de CPU (nuestro caso)
- En **varios computadores** conectados por red (clusters de supercomputadores)

### 3.2 Los conceptos fundamentales de MPI

**rank** — el número de identidad de cada proceso.

Con `mpirun -np 4`, los procesos son rank 0, rank 1, rank 2, rank 3.
El rank 0 es el "jefe" o "root" por convención del programa.

**size** — la cantidad total de procesos lanzados.

Con `mpirun -np 4`, `size = 4`. Con `mpirun -np 2`, `size = 2`.

**MPI_COMM_WORLD** — el "comunicador" que incluye a absolutamente todos los procesos.

Cada función MPI necesita saber a qué grupo de procesos le habla. `MPI_COMM_WORLD`
es el grupo que los incluye a todos. En nuestro programa siempre usamos este.

**Clave:** todos los procesos ejecutan **el mismo código** (`criba.c`), pero su
`rank` diferente hace que se comporten diferente. Es como dar el mismo manual a
4 personas donde el manual dice:
```
"Si eres la persona #0: haz A, luego manda los datos a los demás"
"Si eres cualquier otra persona: espera los datos y luego haz B"
```

```
mpirun -np 4 ./criba
     │
     └── El SO lanza 4 copias del programa
         ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌─────────────┐
         │ ./criba     │ │ ./criba     │ │ ./criba     │ │ ./criba     │
         │ rank=0      │ │ rank=1      │ │ rank=2      │ │ rank=3      │
         │ (jefe)      │ │             │ │             │ │             │
         └─────────────┘ └─────────────┘ └─────────────┘ └─────────────┘
           mismo código, pero rank distinto → comportamiento distinto
```

### 3.3 Las funciones de inicialización

```c
MPI_Init(&argc, &argv);                      // Enciende MPI. Va al principio.
MPI_Comm_rank(MPI_COMM_WORLD, &rank);        // "¿Cuál soy yo?" → guarda en rank
MPI_Comm_size(MPI_COMM_WORLD, &size);        // "¿Cuántos somos?" → guarda en size
double t = MPI_Wtime();                      // Devuelve el tiempo actual en segundos
MPI_Finalize();                              // Apaga MPI. Va al final siempre.
```

Regla importante: si un proceso llama `MPI_Init`, **debe** llamar `MPI_Finalize` antes
de terminar. Si no, el programa puede colgarse o dar errores raros.

### 3.4 MPI_Send y MPI_Recv — "mandar y recibir una carta certificada"

Estas son las dos funciones de comunicación **punto a punto** de MPI: un proceso
**envía** (`MPI_Send`) un mensaje a OTRO proceso específico, y ese otro lo
**recibe** (`MPI_Recv`). Siempre van en pareja: por cada `Send` debe haber un
`Recv` que lo espere.

Piensa en **cartas certificadas**: el Proceso 0 escribe una carta y la manda a un
destinatario concreto (Proceso 1, luego Proceso 2, luego Proceso 3). Cada
destinatario tiene que estar esperando para recibir SU carta.

**Los argumentos:**

```c
MPI_Send(datos, count, tipo, dest,   tag, comm);
//       ↑qué    ↑cuántos ↑tipo ↑a quién ↑etiqueta ↑comunicador

MPI_Recv(datos, count, tipo, source, tag, comm, status);
//       ↑dónde guardo   ↑tipo ↑de quién ↑etiqueta      ↑info del mensaje
```

- **dest / source**: el *rank* del proceso con quien me comunico.
- **tag**: una etiqueta numérica para distinguir un mensaje de otro entre el mismo
  par de procesos (p.ej. tag 0 = "el conteo", tag 1 = "el arreglo"). El `Recv`
  debe pedir el MISMO tag que usó el `Send`.
- **status**: información extra del mensaje recibido. Si no la necesitas, pasas
  `MPI_STATUS_IGNORE`.

**Situación:** el Proceso 0 calculó los primos base. Los demás los necesitan pero
no los tienen. Como no hay una función "emitir a todos", rank 0 le manda una copia
a cada proceso **uno por uno**, dentro de un bucle:

```c
if (rank == 0) {
    for (int dest = 1; dest < size; dest++) {       // a cada compañero...
        MPI_Send(&base_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);   // tag 0
        MPI_Send(base_primes, base_count, MPI_INT, dest, 1, MPI_COMM_WORLD); // tag 1
    }
} else {
    MPI_Recv(&base_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    base_primes = malloc(base_count * sizeof(int));   // ahora sé cuánto reservar
    MPI_Recv(base_primes, base_count, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
```

```
ANTES de los Send/Recv:
  Proceso 0: base_primes = [2, 3, 5, 7, 11, 13, ...]   ← tiene los datos
  Proceso 1: base_primes = ??? (memoria sin inicializar)  ← no sabe nada
  Proceso 2: base_primes = ???
  Proceso 3: base_primes = ???

         ↓ rank 0 hace MPI_Send a 1, luego a 2, luego a 3
         ↓ cada uno hace MPI_Recv y reserva memoria

DESPUÉS de los Send/Recv:
  Proceso 0: base_primes = [2, 3, 5, 7, 11, 13, ...]   ← igual que antes
  Proceso 1: base_primes = [2, 3, 5, 7, 11, 13, ...]   ← recibió copia exacta!
  Proceso 2: base_primes = [2, 3, 5, 7, 11, 13, ...]   ← recibió copia exacta!
  Proceso 3: base_primes = [2, 3, 5, 7, 11, 13, ...]   ← recibió copia exacta!
```

**El resultado final es el mismo que daría un "broadcast"**: todos terminan con la
misma copia. La diferencia es que aquí lo construimos nosotros con un bucle de
envíos punto a punto, viendo explícitamente cada par emisor→receptor.

**Regla clave:** cada `MPI_Send` necesita un `MPI_Recv` que lo espere (mismo
`source`, `tag` y `comm`). Si un proceso hace `Recv` y nadie le manda nada, se
queda colgado esperando para siempre — como un empleado sentado esperando una
carta que nunca llega. Por eso rank 0 NO se manda a sí mismo: ya tiene los datos,
y el bucle arranca en `dest = 1`.

**¿Por qué dos mensajes seguidos (conteo y luego arreglo)?**

Los Procesos 1, 2, 3 necesitan reservar memoria para recibir los primos. Para saber
cuánta memoria reservar, primero necesitan saber cuántos primos hay (`base_count`).
Por eso mandamos el conteo primero (tag 0), reservamos memoria con `malloc`, y luego
mandamos el arreglo (tag 1).

### 3.5 Recolectar conteos con Send/Recv — "cada quien reporta al jefe"

Es la operación inversa al reparto: ahora cada proceso le **envía** su resultado a
rank 0, y rank 0 lo **recibe** de todos, uno por uno, guardándolos en un arreglo.

**Situación:** cada proceso terminó de contar cuántos primos hay en su segmento.
Rank 0 necesita sumarlos todos para dar el total.

```c
if (rank == 0) {
    conteos[0] = conteo_local;                  // rank 0 ya tiene el suyo
    for (int src = 1; src < size; src++)        // recibe del 1, del 2, del 3...
        MPI_Recv(&conteos[src], 1, MPI_LONG_LONG, src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
} else {
    MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 2, MPI_COMM_WORLD);  // tag 2
}
```

```
ANTES de recolectar:
  Proceso 0: conteo_local = 183,072   (primos del 2 al 2,500,000)
  Proceso 1: conteo_local = 165,441   (primos del 2,500,001 al 5,000,000)
  Proceso 2: conteo_local = 159,748   (primos del 5,000,001 al 7,500,000)
  Proceso 3: conteo_local = 156,318   (primos del 7,500,001 al 10,000,000)

         ↓ procesos 1,2,3 hacen MPI_Send; rank 0 hace MPI_Recv en bucle

DESPUÉS de recolectar (solo el Proceso 0 tiene esto):
  Proceso 0: conteos[0] = 183,072
             conteos[1] = 165,441
             conteos[2] = 159,748
             conteos[3] = 156,318

  Ahora rank 0 suma: 183072 + 165441 + 159748 + 156318 = 664,579 primos ✓
```

Los Procesos 1, 2 y 3 no reciben el arreglo completo — solo rank 0 lo tiene.

---

## PARTE 4: La estrategia de paralelización

### 4.1 Distribuir el rango entre procesos

Tenemos 10,000,000 números y P procesos. La estrategia: darle a cada proceso
aproximadamente 1/P del rango.

```
Todos los números: [2 ...................................... 10,000,000]

Con 4 procesos:
  Proceso 0: [2 ............... 2,500,000]   ~2.5 millones de números
  Proceso 1: [2,500,001 ....... 5,000,000]   ~2.5 millones de números
  Proceso 2: [5,000,001 ....... 7,500,000]   ~2.5 millones de números
  Proceso 3: [7,500,001 ..... 10,000,000]   ~2.5 millones de números
```

Cada proceso solo necesita memoria para su porción: ~2.5 MB en vez de los 10 MB
que necesitaría un proceso secuencial. Con N=1,000,000,000 esto se vuelve crítico:
un proceso secuencial necesitaría 1 GB, pero 100 procesos solo 10 MB cada uno.

### 4.2 La fórmula de distribución — paso a paso

El código calcula `low` (inicio) y `high` (fin) del segmento de cada proceso:

```c
long long low  = 2 + (long long)rank * (N - 1) / size;
long long high = (rank == size - 1)
                 ? N
                 : 2 + (long long)(rank + 1) * (N - 1) / size - 1;
```

Aplicando con rank=1, size=4, N=10,000,000:

```
low  = 2 + 1 * (10,000,000 - 1) / 4
     = 2 + 1 * 9,999,999 / 4          ← división entera
     = 2 + 2,499,999
     = 2,500,001

high = 2 + (1+1) * 9,999,999 / 4 - 1
     = 2 + 2 * 9,999,999 / 4 - 1      ← 2 * 9999999 = 19999998
     = 2 + 19,999,998 / 4 - 1         ← división entera: 19999998/4 = 4999999
     = 2 + 4,999,999 - 1
     = 5,000,000
```

El Proceso 1 trabaja en el rango [2,500,001 ... 5,000,000].

El `? :` en el código es un "if compacto":
- Si soy el último proceso (`rank == size - 1`), mi `high = N` (llego hasta el final)
- Si no, uso la fórmula (para no dejar un hueco entre mi fin y el inicio del siguiente)

### 4.3 El truco del índice local — el concepto más importante

Aquí viene el concepto que más confunde, así que lo explicamos con mucho cuidado.

**El problema:** el Proceso 2 trabaja con los números del 5,000,001 al 7,500,000.
Necesita un arreglo para marcar cuáles son primos y cuáles compuestos.

Podría crear un arreglo gigante de 10,000,001 celdas (una por cada número del 0 al
10,000,000), pero eso desperdiciaría memoria. En cambio, crea un arreglo de solo
2,500,000 celdas (una por cada número de su segmento).

**La clave:** si mi segmento empieza en `low = 5,000,001`, entonces:
```
El número 5,000,001 → se guarda en segmento[0]      (= 5,000,001 - 5,000,001)
El número 5,000,002 → se guarda en segmento[1]      (= 5,000,002 - 5,000,001)
El número 5,000,003 → se guarda en segmento[2]      (= 5,000,003 - 5,000,001)
...
El número 5,000,010 → se guarda en segmento[9]      (= 5,000,010 - 5,000,001)
...
El número 7,500,000 → se guarda en segmento[2,499,999]
```

La fórmula para convertir **número global → índice local**:
```
índice_local = número_global - low
```

En el código esto aparece como `segmento[j - low] = 1` cuando queremos marcar
el número `j` como compuesto.

**Ejemplo concreto:**
Si queremos marcar el 5,000,009 como compuesto (en el Proceso 2 con low=5,000,001):
```c
segmento[5000009 - 5000001] = 1;   // segmento[8] = 1
```

### 4.4 Encontrar el primer múltiplo de p en el segmento

Este es el truco matemático central del tamizado.

**Situación:** estoy en el Proceso 2 (rango [5,000,001 ... 7,500,000]) y quiero
marcar todos los múltiplos del primo p=3.

Los múltiplos de 3 son: ..., 4,999,998, 5,000,001, 5,000,004, 5,000,007, ...
Necesito el primero que cae dentro de mi rango.

**La fórmula:**
```
primer_múltiplo = ceil(low / p) * p
```

Donde `ceil` = "techo" = redondear hacia arriba.

Ejemplo con low=5,000,001, p=3:
```
1. Dividir:     5,000,001 / 3 = 1,666,667.0   (exacto)
2. Aplicar ceil: ceil(1,666,667.0) = 1,666,667
3. Multiplicar:  1,666,667 × 3 = 5,000,001
```
El primer múltiplo de 3 en nuestro segmento es 5,000,001 mismo.

Otro ejemplo con low=5,000,002, p=3:
```
1. Dividir:     5,000,002 / 3 = 1,666,667.33...
2. Aplicar ceil: ceil(1,666,667.33) = 1,666,668
3. Multiplicar:  1,666,668 × 3 = 5,000,004
```
El primer múltiplo de 3 que cae en [5,000,002 ...] es 5,000,004.

**¿Cómo se calcula ceil en aritmética entera de C?**

C no tiene `ceil` para enteros (solo para flotantes, y son imprecisos para números
grandes). La equivalencia en enteros es:

```c
ceil(a / b) = (a + b - 1) / b    // división entera
```

Por eso en el código la fórmula es:
```c
long long start = ((low + p - 1) / p) * p;
```

Verificación: low=5,000,002, p=3:
```
start = ((5000002 + 3 - 1) / 3) * 3
      = (5000004 / 3) * 3           ← 5000004 / 3 = 1666668 exacto
      = 1666668 * 3
      = 5,000,004 ✓
```

**El caso especial: no marcar al primo mismo como compuesto**

Si el primo p está dentro de nuestro segmento, la fórmula lo devuelve como
"primer múltiplo de p" (porque p = 1×p). Pero p es primo, ¡no debemos marcarlo!

El código lo maneja con:
```c
if (start == p) start += p;   // salta al siguiente múltiplo real (2p)
```

Ejemplo: si p=5 y nuestro segmento empieza en low=2:
```
start = ((2 + 4) / 5) * 5 = (6 / 5) * 5 = 1 * 5 = 5
start == p (ambos son 5), entonces:
start = 5 + 5 = 10   ← empezamos a marcar desde 10, no desde 5
```

---

## PARTE 5: El código explicado línea por línea

Ahora leemos `criba.c` de principio a fin. Para cada bloque explicamos qué hace,
por qué existe, y qué significan las líneas más crípticas.

### Los `#include` y el `#define`

```c
#include <mpi.h>      // Agrega las definiciones de MPI_Init, MPI_Send, etc.
#include <stdio.h>    // Agrega printf, fprintf, fopen, fclose
#include <stdlib.h>   // Agrega malloc, calloc, free
#include <math.h>     // Agrega sqrt()
#include <string.h>   // Agrega funciones de manejo de memoria (usada internamente)

#define N 10000000    // Constante: límite superior de la búsqueda
```

`#include` le dice al compilador "agrega las definiciones de esta librería". Sin
`<mpi.h>`, el compilador no sabría qué es `MPI_Init` y daría error.

`#define N 10000000` es una **macro**: cada vez que el código escribe `N`, el
compilador lo reemplaza por `10000000` antes de compilar. Es más legible y fácil
de cambiar que tener el número "mágico" desperdigado por todo el código.

### PASO 1: Inicializar MPI

```c
int rank, size;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
double t_inicio = MPI_Wtime();
```

**¿Qué hace?** Enciende MPI, pregunta "¿cuál soy yo?" y "¿cuántos somos?",
y arranca el cronómetro.

**Desglose:**

- `int rank, size;` — declara dos variables enteras. Todavía están vacías.

- `MPI_Init(&argc, &argv)` — inicializa el entorno de comunicación MPI. El `&` significa
  "la dirección de" (un puntero). MPI puede necesitar modificar los argumentos del programa,
  por eso recibe punteros y no copias. **Debe ser siempre la primera llamada MPI.**

- `MPI_Comm_rank(MPI_COMM_WORLD, &rank)` — pregunta al sistema: "¿cuál es mi rank?".
  El resultado se guarda en la variable `rank`. Después de esta línea, si `rank = 2`,
  este proceso sabe que es el Proceso 2.

- `MPI_Comm_size(MPI_COMM_WORLD, &size)` — pregunta "¿cuántos procesos hay en total?".
  Si corriste `mpirun -np 4`, `size = 4`.

- `MPI_Wtime()` — devuelve el tiempo actual en segundos con microsegundos de precisión.
  Guardamos este valor para luego calcular cuánto tardó el programa.

Después de estas 4 líneas, cada proceso sabe: *"Yo soy el proceso número `rank`
de un total de `size` procesos."*

### PASO 2: Calcular primos base (solo rank 0)

```c
int raiz_n = (int)sqrt((double)N);   // ≈ 3162
int base_count = 0;
int *base_primes = NULL;

if (rank == 0) {
    char *marcado = calloc(raiz_n + 1, sizeof(char));
    marcado[0] = marcado[1] = 1;   // 0 y 1 no son primos

    for (int i = 2; (long long)i * i <= raiz_n; i++) {
        if (!marcado[i]) {
            for (int j = i * i; j <= raiz_n; j += i)
                marcado[j] = 1;
        }
    }

    for (int i = 2; i <= raiz_n; i++)
        if (!marcado[i]) base_count++;

    base_primes = malloc(base_count * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= raiz_n; i++)
        if (!marcado[i]) base_primes[idx++] = i;

    free(marcado);
}
```

**¿Qué hace?** Solo el Proceso 0 hace una criba de Eratóstenes secuencial sobre
los números del 2 al 3162. El resultado es el arreglo `base_primes` con los 446
primos menores o iguales a 3162.

**¿Por qué solo rank 0?** Porque alguien tiene que hacerlo una vez. Rank 0 lo hace
y luego lo reparte a todos enviándolo con `MPI_Send`. Todos los procesos podrían
calcularlo independientemente (son pocos números), pero lo hace solo rank 0 y lo
reparte con Send/Recv para demostrar el patrón de comunicación punto a punto de MPI.

**Desglose línea a línea:**

- `calloc(raiz_n + 1, sizeof(char))` — reserva un arreglo de 3163 celdas
  (una por cada número del 0 al 3162) y las pone **todas en 0**.
  `calloc` = "allocate and clear". El convenio es: 0 = primo, 1 = compuesto.

  ¿Por qué `char` y no `int`? Un `char` ocupa 1 byte; un `int` ocupa 4 bytes.
  Como solo necesitamos guardar 0 o 1, usar `char` ahorra 4 veces la memoria.

- `marcado[0] = marcado[1] = 1;` — El 0 y el 1 no son primos por definición,
  así que los marcamos como compuestos desde el inicio.

- El primer `for` (criba): para cada i que no esté marcado (es decir, primo),
  marca todos sus múltiplos empezando en `i*i`.
  ¿Por qué empezar en `i*i` y no en `2*i`? Porque los múltiplos menores
  (2i, 3i, ..., (i-1)*i) ya fueron marcados por primos anteriores.

- El segundo `for` cuenta cuántos primos hay → `base_count`.

- `malloc(base_count * sizeof(int))` — reserva memoria para exactamente
  `base_count` enteros. A diferencia de `calloc`, `malloc` no pone ceros
  (el contenido inicial es "basura"), pero lo llenamos en el siguiente for.

- El tercer `for` copia los primos al arreglo `base_primes`.
  `base_primes[idx++] = i;` guarda `i` en la posición `idx` y luego incrementa `idx`.

- `free(marcado)` — devuelve la memoria del arreglo `marcado` al sistema operativo.
  Ya no la necesitamos; los primos están en `base_primes`.

### PASO 3: Enviar los primos base a todos (Send/Recv)

```c
if (rank == 0) {
    for (int dest = 1; dest < size; dest++) {
        MPI_Send(&base_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        MPI_Send(base_primes, base_count, MPI_INT, dest, 1, MPI_COMM_WORLD);
    }
} else {
    MPI_Recv(&base_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    base_primes = malloc(base_count * sizeof(int));
    MPI_Recv(base_primes, base_count, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
```

**¿Qué hace?** Rank 0 envía los primos base a cada proceso, uno por uno, con
comunicación punto a punto (sin usar un broadcast colectivo).

**Desglose:**

1. **`if (rank == 0)`** — solo rank 0 entra al bucle de envíos. Recorre los
   destinos `dest = 1, 2, 3` y a cada uno le manda DOS mensajes:
   - `MPI_Send(&base_count, ..., dest, 0, ...)` (tag 0): cuántos primos hay.
   - `MPI_Send(base_primes, ..., dest, 1, ...)` (tag 1): el arreglo completo.

   Rank 0 NO se manda a sí mismo: ya tiene `base_primes` del paso anterior, por
   eso el bucle arranca en `dest = 1`.

2. **`else` (procesos 1, 2, 3)** — cada uno recibe en el mismo orden:
   - `MPI_Recv(&base_count, ..., 0, 0, ...)` recibe el conteo (tag 0).
   - `malloc(...)` reserva memoria ahora que ya sabe cuántos primos hay.
   - `MPI_Recv(base_primes, ..., 0, 1, ...)` recibe el arreglo (tag 1).

   `MPI_STATUS_IGNORE` indica que no nos interesa la info extra del mensaje.

Después de este bloque, **todos los procesos** tienen exactamente el mismo
arreglo `base_primes` con los 446 primos hasta 3162 — el mismo resultado que
daría un broadcast, pero construido con envíos punto a punto.

### PASO 4: Determinar el segmento de este proceso

```c
long long low  = 2 + (long long)rank * (N - 1) / size;
long long high = (rank == size - 1)
                 ? N
                 : 2 + (long long)(rank + 1) * (N - 1) / size - 1;
long long seg_size = high - low + 1;

char *segmento = calloc(seg_size, sizeof(char));
```

**¿Qué hace?** Cada proceso calcula su rango personal [low, high] y crea el
arreglo local donde marcará los compuestos.

**Desglose:**

- `long long` — tipo de dato de C que ocupa 8 bytes y puede representar números
  hasta ~9,200,000,000,000,000,000. Usamos `long long` porque aunque N=10,000,000
  cabe en un `int`, las multiplicaciones intermedias como `rank * (N-1)` podrían
  desbordar un `int` si N fuera más grande.

- `(long long)rank * (N-1) / size` — el cast `(long long)` convierte `rank` a
  `long long` antes de multiplicar, evitando overflow.

- El operador `? :` es un if-else compacto:
  ```
  condición ? valor_si_verdadero : valor_si_falso
  ```
  El último proceso recibe `high = N` para que llegue exactamente hasta el final.

- `calloc(seg_size, sizeof(char))` — crea el arreglo del segmento, con todo en 0.
  Convenio: 0 = primo (aún no descartado), 1 = compuesto (marcado).

Para N=10,000,000 con 4 procesos, los segmentos son:
```
Proceso 0: low=2,         high=2,500,000   (tamaño: 2,499,999)
Proceso 1: low=2,500,001, high=5,000,000   (tamaño: 2,500,000)
Proceso 2: low=5,000,001, high=7,500,000   (tamaño: 2,500,000)
Proceso 3: low=7,500,001, high=10,000,000  (tamaño: 2,500,000)
```

### PASO 5: Tamizan el segmento local

Esta es la parte central del algoritmo. Cada proceso trabaja independientemente
en su propio segmento — no hay comunicación MPI aquí.

```c
for (int i = 0; i < base_count; i++) {
    long long p = base_primes[i];

    long long start = ((low + p - 1) / p) * p;
    if (start == p) start += p;

    for (long long j = start; j <= high; j += p)
        segmento[j - low] = 1;
}
```

**¿Qué hace?** Para cada primo base p, marca todos los múltiplos de p que están
dentro de [low, high] como compuestos (1).

**Desglose:**

- `long long p = base_primes[i];` — toma el i-ésimo primo base.

- `((low + p - 1) / p) * p` — calcula el primer múltiplo de p que es ≥ low.
  Es `ceil(low/p)*p` expresado en aritmética entera (sin punto flotante).
  Ver explicación detallada en la Parte 4.4.

- `if (start == p) start += p;` — si el primo p mismo está en nuestro segmento
  (es decir, p ≥ low y p ≤ high), saltamos al siguiente múltiplo para no marcarlo.

- `for (long long j = start; j <= high; j += p)` — recorre todos los múltiplos
  de p en [start, high]. Los pasos son de p en p.

- `segmento[j - low] = 1;` — marca el número j como compuesto. El índice local
  es `j - low` (ver explicación en Parte 4.3).

**Ejemplo concreto** con Proceso 1 (low=2,500,001), p=7:
```
start = ((2500001 + 6) / 7) * 7 = (2500007 / 7) * 7 = 357144 * 7 = 2,500,008
j = 2,500,008 → segmento[2500008 - 2500001] = segmento[7] = 1  (marcamos 2,500,008)
j = 2,500,015 → segmento[14] = 1
j = 2,500,022 → segmento[21] = 1
... y así hasta llegar a high=5,000,000
```

### PASO 6: Contar los primos en el segmento

```c
long long conteo_local = 0;
for (long long i = 0; i < seg_size; i++)
    if (!segmento[i]) conteo_local++;
```

**¿Qué hace?** Recorre el arreglo local y cuenta las celdas que siguen en 0
(no fueron marcadas → son primos).

`!segmento[i]` en C significa "si segmento[i] es 0 (falso), esta condición
es verdadera". Es equivalente a escribir `segmento[i] == 0`.

`long long conteo_local` usa 64 bits porque aunque para N=10M los conteos son
del orden de 150,000–180,000, para N mucho más grande podrían superar los 2,147,483,647
que caben en un `int` de 32 bits.

### PASO 7: Recolectar conteos en rank 0

```c
long long *conteos = NULL;
if (rank == 0) {
    conteos = malloc(size * sizeof(long long));
    conteos[0] = conteo_local;          // rank 0 ya tiene el suyo
    for (int src = 1; src < size; src++)
        MPI_Recv(&conteos[src], 1, MPI_LONG_LONG, src, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
} else {
    MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 2, MPI_COMM_WORLD);
}

double t_fin  = MPI_Wtime();
double tiempo = t_fin - t_inicio;
```

**¿Qué hace?** Rank 0 recibe el conteo de primos de cada proceso y los guarda
en el arreglo `conteos`, usando comunicación punto a punto. También se detiene
el cronómetro.

**Desglose:**
```
Procesos 1, 2, 3:  MPI_Send(&conteo_local, 1, MPI_LONG_LONG, 0, 2, comm)
                   → cada uno manda su conteo a rank 0 (tag 2).

Rank 0:  conteos[0] = conteo_local      → su propio conteo, sin enviarlo
         for src = 1..size-1:
             MPI_Recv(&conteos[src], 1, MPI_LONG_LONG, src, 2, comm, IGNORE)
         → recibe uno por uno y los guarda en conteos[src] según el emisor.
```

Resultado: `conteos[0]` tiene el conteo del Proceso 0, `conteos[1]` del Proceso 1, etc.
Es la operación inversa al PASO 3: allí rank 0 enviaba a todos; aquí todos le envían
a rank 0.

**`MPI_LONG_LONG`** — constante MPI que indica "el dato es de tipo long long de C".
MPI necesita saber el tipo para saber cuántos bytes enviar.

### PASO 8: Rank 0 muestra resultados y escribe archivos

```c
if (rank == 0) {
    long long total = 0;
    for (int i = 0; i < size; i++) total += conteos[i];

    printf("\n=== Criba de Eratostenes Distribuida (N=%d) ===\n", N);
    printf("Procesos MPI        : %d\n", size);
    printf("Primos encontrados  : %lld\n", total);
    printf("Tiempo de ejecucion : %.6f segundos\n\n", tiempo);

    // ... tabla por proceso ...

    printf("\nTIEMPO: %.6f\n", tiempo);   // para el script de benchmark

    FILE *fp = fopen("resultados.dat", "w");
    // ... escribe el archivo de datos para Gnuplot ...
    fclose(fp);

    free(conteos);
}
```

**¿Por qué solo rank 0?** Porque solo rank 0 tiene el arreglo `conteos` completo.
Los otros procesos solo conocen su propio `conteo_local`.

**`%lld`** — especificador de formato de `printf` para `long long`. El doble `l`
es necesario para distinguirlo de `int` (%d) y `long` (%ld).

**`%.6f`** — formato para `double` con 6 decimales (microsegundos de precisión).

**`printf("\nTIEMPO: %.6f\n")`** — esta línea en particular la lee el script del
Makefile cuando corre el benchmark, para extraer el tiempo y calcular el speedup.

### Limpieza final

```c
free(base_primes);
free(segmento);
MPI_Finalize();
return 0;
```

**`free()`** — devuelve al sistema la memoria reservada con `malloc`/`calloc`.
Si no hacemos esto (lo que se llama un "memory leak"), la memoria queda ocupada
hasta que el proceso termina. En nuestro caso el proceso termina justo después,
así que no causa problemas prácticos, pero es buena práctica siempre.

**`MPI_Finalize()`** — apaga el entorno MPI limpiamente. Libera recursos internos,
cierra canales de comunicación, etc. Después de esta llamada no se puede usar
ninguna función MPI. **Siempre va al final.**

---

## PARTE 6: Compilar y ejecutar el programa

### ¿Qué es `mpicc`?

`mpicc` no es un compilador diferente. Es un **wrapper** (envoltorio) alrededor
de `gcc` que agrega automáticamente las flags necesarias para compilar con MPI:

```
mpicc criba.c  →  internamente ejecuta algo equivalente a:
gcc criba.c -I/usr/include/mpich -L/usr/lib64/mpich/lib -lmpi -o criba
```

Tú no necesitas recordar esas flags — `mpicc` las agrega por ti.

### ¿Qué hace `mpirun -np 4`?

`mpirun` no requiere 4 computadores. Crea 4 procesos en tu computador actual y
los conecta con MPI. Los 4 procesos se reparten entre los núcleos de tu CPU.

```
mpirun -np 4 ./criba
  │
  └── Lanza 4 procesos:
      ┌──────────────┐ ┌──────────────┐ ┌──────────────┐ ┌──────────────┐
      │ ./criba      │ │ ./criba      │ │ ./criba      │ │ ./criba      │
      │ rank=0       │ │ rank=1       │ │ rank=2       │ │ rank=3       │
      │ (en CPU núcleo 0)│(en núcleo 1)│ (en núcleo 2)│ (en núcleo 3)│
      └──────────────┘ └──────────────┘ └──────────────┘ └──────────────┘
```

### Compilar

```bash
# Opción A: usar el Makefile (recomendado)
make

# Opción B: compilar directamente
mpicc -O2 -Wall -o criba criba.c -lm
```

Las flags del compilador:
- `-O2` — optimizaciones de velocidad (nivel 2, buen balance velocidad/tiempo de compilación)
- `-Wall` — activa todos los warnings ("advertencias") del compilador — ayuda a detectar errores
- `-lm` — enlaza la librería matemática (necesaria para `sqrt()`)

Si `mpicc` no está en PATH, en Fedora prueba:
```bash
/usr/lib64/mpich/bin/mpicc -O2 -Wall -o criba criba.c -lm
```

### Ejecutar

```bash
mpirun -np 1 ./criba   # 1 proceso (secuencial — útil para comparar tiempos)
mpirun -np 2 ./criba   # 2 procesos en paralelo
mpirun -np 4 ./criba   # 4 procesos en paralelo
```

### Verificar que el resultado es correcto

**El programa SIEMPRE debe imprimir:**
```
Primos encontrados  : 664579
```

Ese es el número exacto de primos entre 2 y 10,000,000 (verificado matemáticamente).
Si el número es diferente, hay un bug — ya sea en el algoritmo o en la comunicación MPI.

### Correr el benchmark completo

```bash
make benchmark
```

Corre el programa con 1, 2 y 4 procesos, mide los tiempos e imprime:

```
+----------+-------------+---------+
| Procesos |    Tiempo   | Speedup |
+----------+-------------+---------+
|        1 |    0.0234 s |   1.000 |
|        2 |    0.0117 s |   1.997 |
|        4 |    0.0094 s |   2.474 |
+----------+-------------+---------+
```

También genera el archivo `speedup.dat` para Gnuplot.

### ¿Cómo interpretar el speedup?

```
Speedup = Tiempo con 1 proceso / Tiempo con N procesos
```

- **Speedup = 1.997** con 2 procesos significa: el programa es casi 2 veces más rápido
  que la versión de 1 proceso (muy cerca del ideal).

- **Speedup ideal** con N procesos = N (si 2 procesos hacen exactamente la mitad
  del trabajo cada uno, el tiempo se divide exactamente a la mitad).

- **Speedup real** siempre es menor que el ideal porque hay **overhead de comunicación**:
  el tiempo que tardan los `MPI_Send` y `MPI_Recv`, más el overhead de lanzar procesos,
  sincronización, etc.

```
Eficiencia = (Speedup real / Procesos) × 100%

Con 2 procesos: (1.997 / 2) × 100% = 100% ← casi perfecto
Con 4 procesos: (2.474 / 4) × 100% = 62%  ← cae al subir procesos (más comunicación)
```

---

## PARTE 7: Generar las gráficas con Gnuplot

### Paso 1: Generar los archivos de datos

```bash
# Genera resultados.dat (cuántos primos encontró cada proceso)
mpirun -np 4 ./criba

# Genera speedup.dat (tabla de tiempos para 1, 2 y 4 procesos)
make benchmark
```

### Paso 2: Generar las imágenes PNG

```bash
gnuplot graficar.gp
```

Produce dos imágenes en el directorio actual:
- `distribucion.png` — gráfica de barras: cuántos primos encontró cada proceso
- `speedup.png` — speedup real (línea roja) vs speedup ideal (S=p, línea gris)

Si gnuplot no está instalado:
```bash
sudo dnf install gnuplot    # Fedora/RHEL
sudo apt install gnuplot    # Ubuntu/Debian
```

---

## PARTE 8: Compilar los archivos LaTeX para PDFs

### Informe IEEE

```bash
pdflatex informe_ieee.tex   # Primera pasada
pdflatex informe_ieee.tex   # Segunda pasada (necesaria para referencias cruzadas)
# Resultado: informe_ieee.pdf
```

¿Por qué dos pasadas? LaTeX procesa el documento linealmente. En la primera pasada
genera un archivo auxiliar con los números de página de cada sección. En la segunda
pasada usa ese archivo para que las referencias internas ("ver Sección 3") apunten
a los números correctos.

### Archivo de referencias

```bash
pdflatex referencias.tex
# Resultado: referencias.pdf
```

### Nota sobre las imágenes

Las imágenes `distribucion.png` y `speedup.png` deben existir en el mismo directorio
antes de compilar el informe. Si no existen, LaTeX lanza un warning pero compila igual
(las imágenes aparecen como marcos vacíos).

### Limpiar archivos temporales

```bash
make clean
# Borra el ejecutable criba, los .dat, los .png, y los temporales de LaTeX (.aux, .log, etc.)
```

---

## PARTE 9: Aplicaciones en la vida real

Hasta aquí aprendimos a encontrar primos en paralelo. Pero, ¿esto para qué sirve
**de verdad**? Resulta que tanto los **números primos** como la **técnica de
repartir trabajo entre procesos** están detrás de cosas que usas todos los días.
Vamos por partes.

### 9.1 ¿Para qué sirven los números primos?

Los primos parecen un juego matemático, pero son la base de la seguridad de
internet. Aquí los usos más importantes:

#### a) Criptografía RSA — la seguridad de tu banco y WhatsApp

Cada vez que entras a una página con candado (HTTPS), compras en línea o mandas
un mensaje cifrado, hay números primos trabajando por detrás. El sistema más
famoso se llama **RSA**, y su idea es sorprendentemente simple:

```
Multiplicar es FÁCIL:       61 × 53 = 3233    (lo haces en segundos)
Factorizar es DIFÍCIL:      3233 = ¿? × ¿?    (¿qué dos primos lo forman?)
```

Con números chicos como 3233 es fácil adivinar (61 y 53). Pero RSA usa primos de
**300 dígitos o más**. Multiplicarlos es instantáneo, pero hacer el camino
inverso —descubrir qué dos primos se multiplicaron— le tomaría a la
computadora más potente del mundo **millones de años**.

```
Analogía: es como mezclar dos colores de pintura.
  Mezclar azul + amarillo = verde    → instantáneo
  "Des-mezclar" el verde de vuelta a azul + amarillo → casi imposible
```

Esa asimetría (fácil en un sentido, imposible en el otro) es lo que mantiene
seguros tus datos. Y para crear las claves se necesitan **primos grandes**, que
se buscan con técnicas emparentadas con la criba que programamos.

#### b) Firmas digitales y certificados

La misma matemática de primos permite **firmar** documentos digitalmente: probar
que un mensaje viene de quien dice venir y que nadie lo alteró. Los certificados
que validan que "esta página es realmente la de tu banco" usan firmas basadas en
primos.

#### c) Funciones hash y tablas hash

Una **tabla hash** es una estructura para guardar y buscar datos muy rápido
(la usan bases de datos, diccionarios de programación, cachés). Para repartir los
datos de forma pareja y evitar choques, es común usar un **tamaño primo** de
tabla:

```
Si la tabla tiene tamaño 100 (no primo), muchas claves caen en las mismas casillas.
Si la tabla tiene tamaño 101 (primo), las claves se reparten mucho más parejo.
```

#### d) Generadores de números aleatorios

Los videojuegos, simulaciones y sorteos necesitan números "al azar". Muchos
generadores usan **módulos primos** en sus fórmulas para que la secuencia de
números tarde lo máximo posible en repetirse.

**Resumen 9.1:** los primos no son un juego — son el cimiento de la criptografía
que protege bancos, mensajería y comercio electrónico.

### 9.2 Cazadores de primos gigantes: la criba a escala mundial

Lo más valioso que aprendiste **no es la criba en sí**, sino el **patrón** de
trabajo en paralelo. Lo hiciste en 3 movimientos:

```
1. REPARTIR  → rank 0 envió los primos base con MPI_Send a cada proceso  (send)
2. CALCULAR  → cada uno tamizó su pedazo solo, sin molestar a nadie      (compute)
3. JUNTAR    → rank 0 recibió los conteos con MPI_Recv de cada proceso   (recv)
```

Pues bien: existen proyectos **reales** de investigación matemática que hacen
exactamente esto, pero con miles o millones de computadoras y números
descomunales. Tú lo hiciste hasta 10 millones; ellos llegan hasta números con
**millones de dígitos**. El problema es tan grande que **obligatoriamente** hay
que repartirlo y apoyarse en cribas.

#### a) GIMPS — la cacería de los primos de Mersenne

Un **primo de Mersenne** es un primo de la forma `2^p − 1` (2 elevado a una
potencia, menos 1):

```
2^2  − 1 = 3      ✓ primo
2^3  − 1 = 7      ✓ primo
2^5  − 1 = 31     ✓ primo
2^7  − 1 = 127    ✓ primo
2^11 − 1 = 2047 = 23 × 89   ✗ NO primo (ojo: no todo 2^p−1 es primo)
```

Resulta que **los números primos más grandes que conoce la humanidad son de
Mersenne**. El récord actual tiene **más de 41 millones de dígitos** (si lo
imprimieras, llenarías miles de libros).

¿Cómo se encuentran? Con **GIMPS** (*Great Internet Mersenne Prime Search*), un
proyecto mundial donde **miles de voluntarios** prestan sus computadoras. El
sistema central reparte exponentes candidatos `p` entre todos —igual que tú
repartiste rangos— y cada quien trabaja el suyo.

```
Tu criba:    repartiste rangos de números     [2..2.5M] [2.5M..5M] ...
GIMPS:       reparte exponentes candidatos      p=1000   p=1001   p=1002 ...
Mismo patrón: repartir → calcular → reportar.
```

**¿Dónde entra la criba?** Antes de hacer el test caro, se usan cribas para
**descartar candidatos** rápidamente: los exponentes `p` tienen que ser primos
(¡se generan con una criba como la tuya!) y se eliminan los que tienen factores
pequeños (*trial factoring*). El veredicto final de primalidad se da con un test
especial (Lucas–Lehmer), pero **la criba es el filtro que ahorra el 99% del
trabajo**.

#### b) PrimeGrid — la criba distribuida en estado puro

**PrimeGrid** es otro proyecto mundial, y aquí la conexión es todavía más directa:
tiene subproyectos que se llaman literalmente **"sieving"** (tamizado/criba). Su
trabajo es **exactamente el tuyo**, pero repartido entre computadoras de todo el
planeta: cribar rangos enormes para tachar los candidatos con factores pequeños,
antes de probar cuáles son realmente primos. Así buscan primos gemelos, primos de
Sophie Germain y otras familias especiales.

> Si entendiste la criba de este proyecto, entendiste el corazón de lo que hace
> PrimeGrid. La única diferencia es la escala.

#### c) Probar conjeturas famosas hasta números gigantes

Hay preguntas sobre primos que llevan **siglos** sin demostrarse. Mientras los
matemáticos buscan la prueba, las computadoras las **verifican** hasta números
enormes:

- **Conjetura de Goldbach** ("todo par > 2 es suma de dos primos"): verificada
  por computadora hasta **4,000,000,000,000,000,000** (4 trillones).
- **Conjetura de los primos gemelos** (parejas como 11 y 13, 17 y 19): verificada
  hasta cotas igual de enormes.

Para comprobar esto **hay que generar TODOS los primos** hasta esos límites... y
eso se hace con **cribas segmentadas en paralelo**, exactamente la técnica de
este proyecto. Sin repartir el trabajo, sería imposible.

#### d) Romper la criptografía... también es una criba

¿Recuerdas que RSA (sección 9.1) es seguro porque factorizar es difícil? Pues los
algoritmos **más potentes** para factorizar números enormes se llaman, otra vez,
**cribas**: la *criba del cuerpo de números* (GNFS) y la *criba cuadrática*. Son
cribas mucho más sofisticadas que la de Eratóstenes, pero comparten la idea de
tachar candidatos en paralelo. Con ellas se logró factorizar **RSA-250** (un
número de 250 dígitos) en 2020, usando cientos de computadoras durante meses.

```
Los primos PROTEGEN la información   → construir claves RSA (sección 9.1)
Las cribas distribuidas la ATACAN    → factorizar para romper RSA
Las dos caras de la misma moneda, ambas con cribas y paralelismo.
```

---

**La gran conclusión:** programaste una criba de primos hasta 10 millones, pero
en realidad aprendiste el **mismo patrón** (`MPI_Send` para repartir + cálculo
local + `MPI_Recv` para juntar) que usan los proyectos que cazan los primos más grandes
del mundo, verifican conjeturas centenarias y ponen a prueba la criptografía. Lo
tuyo es la versión de bolsillo —pequeña y verificable ($\pi(10^7)=664{,}579$)— de
algo que se ejecuta a escala planetaria.

---

## Glosario completo

| Término | Significado en lenguaje simple |
|---|---|
| **proceso** | Un programa en ejecución con su propia memoria privada. Con `mpirun -np 4` se crean 4 procesos. |
| **rank** | El número de identidad de un proceso (0, 1, 2, ...). El proceso rank 0 es el "jefe" por convención. |
| **size** | La cantidad total de procesos lanzados. Con `mpirun -np 4`, size=4. |
| **MPI_COMM_WORLD** | El "comunicador" que incluye a absolutamente todos los procesos del programa. |
| **punto a punto** | Comunicación entre DOS procesos concretos: uno envía (`MPI_Send`) y otro recibe (`MPI_Recv`). Es la base de toda la comunicación en MPI. |
| **MPI_Send** | Envía un mensaje a UN proceso destino. Argumentos: datos, count, tipo, dest, tag, comunicador. Bloquea hasta poder enviar. |
| **MPI_Recv** | Recibe un mensaje de UN proceso origen. Debe coincidir con su `Send` en tipo, tag y comunicador. Bloquea hasta que el mensaje llega. |
| **tag** | Etiqueta numérica que distingue mensajes entre el mismo par de procesos. El `Recv` debe pedir el mismo tag que usó el `Send`. |
| **segmento** | La porción del rango [2, N] que le toca procesar a un proceso. Cada uno trabaja en su propio segmento. |
| **primo** | Número entero > 1 que solo se divide exactamente entre 1 y sí mismo: 2, 3, 5, 7, 11, 13, ... |
| **compuesto** | Número que tiene divisores además de 1 y sí mismo: 4=2×2, 6=2×3, 9=3×3, 10=2×5, ... |
| **criba** | Proceso de "filtrar": eliminar los compuestos de una lista para quedarse solo con los primos. |
| **primos base** | Los primos hasta √N (≈ 3162 para N=10M). Son los únicos que se necesitan para tamizan cualquier segmento. |
| **índice local** | La posición de un número dentro del arreglo `segmento`. Fórmula: `índice = número - low`. |
| **calloc** | Función de C que reserva memoria Y la pone toda en 0. `calloc(n, s)` = n celdas de s bytes cada una. |
| **malloc** | Función de C que reserva memoria sin inicializarla (el contenido inicial es "basura"). Más rápido que calloc. |
| **free** | Función de C que devuelve memoria reservada con malloc/calloc al sistema operativo. |
| **char** | Tipo de dato de C que ocupa 1 byte (8 bits). Puede guardar valores 0–255. Aquí lo usamos para guardar 0 (primo) o 1 (compuesto). |
| **long long** | Tipo de dato de C que ocupa 8 bytes. Puede guardar números hasta ~9,200,000,000,000,000,000. |
| **ceil** | Función matemática "techo": redondea hacia arriba. ceil(3.2)=4, ceil(3.0)=3, ceil(-2.7)=-2. |
| **overflow** | Desbordamiento: cuando un número es demasiado grande para el tipo de dato que lo guarda. Ej: un int no puede guardar 3,000,000,000. |
| **speedup** | Cuánto más rápido es el programa paralelo vs el secuencial. Speedup=2.47 = es 2.47 veces más rápido. |
| **eficiencia** | Speedup dividido entre el número de procesos × 100%. Mide qué tan bien se aprovechan los procesos. |
| **overhead de comunicación** | El tiempo "extra" que consumen los mensajes MPI (los `Send`/`Recv`). Hace que el speedup real sea menor que el ideal. |
| **speedup ideal** | El speedup máximo teórico: con N procesos, el ideal es N (trabajo perfectamente dividido, sin overhead). |
| **MPI_Wtime** | Función que devuelve el tiempo actual en segundos con alta precisión. Se usa como cronómetro. |
| **MPI_LONG_LONG** | Constante MPI que representa el tipo `long long` de C. Los mensajes MPI necesitan saber el tipo del dato. |
| **mpicc** | Compilador para programas MPI. Es un wrapper de gcc que agrega automáticamente las flags de MPI. |
| **mpirun** | Lanzador de procesos MPI. `mpirun -np 4 ./criba` crea 4 procesos en el computador actual. |
| **wrapper** | Un programa que "envuelve" a otro agregando funcionalidad extra. `mpicc` envuelve a `gcc`. |
| **memory leak** | Error de programación: reservar memoria con malloc/calloc y no librarla con free. La memoria queda "perdida". |
| **RSA** | Sistema de criptografía de clave pública basado en lo difícil que es factorizar el producto de dos primos grandes. Protege HTTPS, bancos y mensajería. |
| **criptografía** | Técnicas para proteger información (cifrarla, firmarla). La de clave pública moderna depende de números primos grandes. |
| **cluster** | Conjunto de computadoras conectadas por red que trabajan juntas como una sola máquina grande. MPI es el estándar para programarlas. |
| **scatter / gather** | Patrón base del cómputo paralelo: *scatter* = repartir pedazos del trabajo entre procesos; *gather* = recolectar los resultados parciales. |
| **primo de Mersenne** | Primo de la forma `2^p − 1` (ej: 3, 7, 31, 127). Los primos más grandes conocidos son de Mersenne; el récord supera los 41 millones de dígitos. |
| **GIMPS** | *Great Internet Mersenne Prime Search*: proyecto mundial de cómputo distribuido que busca primos de Mersenne repartiendo candidatos entre miles de voluntarios. |
| **PrimeGrid** | Proyecto de cómputo distribuido que ejecuta cribas (*sieving*) a escala mundial para buscar primos gemelos, de Sophie Germain y otras familias. |
| **criba del cuerpo de números (GNFS)** | El algoritmo más potente para factorizar enteros enormes. Es una "criba" sofisticada y masivamente paralela; se usó para romper RSA-250. |
| **conjetura de Goldbach** | Afirmación (no demostrada) de que todo número par > 2 es suma de dos primos. Verificada por computadora hasta 4×10¹⁸ usando cribas paralelas. |
