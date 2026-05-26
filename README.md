# Computación Paralela — Proyectos con MPI y OpenMP

Repositorio de proyectos de la materia **Computación Paralela** — Universidad Mayor de San Andrés.

Cada proyecto implementa un algoritmo clásico usando una tecnología de paralelismo distinta,
incluye informe IEEE, guía de estudio desde cero y resultados experimentales reproducibles.

---

## Proyectos

### 1. [Criba de Eratóstenes con MPI](./criba-mpi/)

Encuentra todos los números primos hasta 10,000,000 distribuyendo el trabajo entre múltiples
procesos con **MPI (Message Passing Interface)**.

| Tecnología | Lenguaje | Paralelismo |
|-----------|----------|-------------|
| MPI (MPICH) | C | Memoria distribuida |

**Compilar y ejecutar:**
```bash
cd criba-mpi
make          # compila
make run      # ejecuta con 4 procesos
make benchmark  # mide speedup con 1, 2 y 4 procesos
```

**Resultados incluidos:** `resultados.dat`, `speedup.dat`, `distribucion.png`, `speedup.png`

---

### 2. [K-Means Clustering con OpenMP](./kmeans-openmp/)

Agrupa las 150 flores del dataset Iris en 3 clusters paralelizando la fase de asignación
con **OpenMP (Open Multi-Processing)**.

| Tecnología | Lenguaje | Dataset | Paralelismo |
|-----------|----------|---------|-------------|
| OpenMP | C | Iris (150 flores, 4 features) | Memoria compartida |

**Compilar y ejecutar:**
```bash
cd kmeans-openmp
make          # compila
make run      # ejecuta (mide 1, 2 y 4 hilos)
make all-plots  # ejecuta y genera gráficas con gnuplot
```

**Resultados incluidos:** `resultados.dat`, `speedup.dat`, `clusters.png`, `speedup.png`

---

## Estructura del repositorio

```
.
├── criba-mpi/
│   ├── src/
│   │   └── criba.c              ← código fuente
│   ├── informe/
│   │   ├── informe_ieee.tex     ← informe LaTeX
│   │   ├── informe_ieee.pdf     ← informe compilado
│   │   ├── referencias.tex
│   │   ├── referencias.pdf
│   │   └── IEEEtran.cls
│   ├── resultados/
│   │   ├── resultados.dat       ← primos por proceso
│   │   ├── speedup.dat          ← mediciones de speedup
│   │   ├── distribucion.png     ← gráfica de distribución
│   │   └── speedup.png          ← gráfica de speedup
│   ├── Makefile
│   ├── graficar.gp              ← script gnuplot
│   └── GUIA_ESTUDIO.md          ← guía para principiantes
└── kmeans-openmp/
    ├── src/
    │   └── kmeans.c             ← código fuente
    ├── datos/
    │   └── iris.csv             ← dataset Iris (UCI)
    ├── informe/
    │   ├── informe_ieee.tex     ← informe LaTeX
    │   ├── informe_ieee.pdf     ← informe compilado
    │   ├── referencias.tex
    │   ├── referencias.pdf
    │   └── IEEEtran.cls
    ├── resultados/
    │   ├── resultados.dat       ← clusters asignados
    │   ├── speedup.dat          ← mediciones de speedup
    │   ├── clusters.png         ← scatter plot de clusters
    │   └── speedup.png          ← gráfica de speedup
    ├── Makefile
    ├── graficar.gp              ← script gnuplot
    └── GUIA_ESTUDIO.md          ← guía para principiantes
```

---

## Dependencias

### Proyecto criba-mpi
```bash
# Fedora/RHEL
sudo dnf install mpich-devel make gnuplot texlive-scheme-full

# Ubuntu/Debian
sudo apt install libmpich-dev mpich make gnuplot texlive-full
```

### Proyecto kmeans-openmp
```bash
# Fedora/RHEL
sudo dnf install gcc make gnuplot texlive-scheme-full

# Ubuntu/Debian
sudo apt install gcc make gnuplot texlive-full
```

> OpenMP viene incluido con GCC — no requiere instalación adicional.
