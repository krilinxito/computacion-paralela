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

# ─────────────────────────────────────────────────────────────
# GRÁFICA 1: Distribución de primos por proceso
# ─────────────────────────────────────────────────────────────
# Leemos resultados.dat con formato:
#   proceso  inicio  fin  cantidad_primos
# Columna 1 = rank del proceso (eje X)
# Columna 4 = primos encontrados (eje Y)
# ─────────────────────────────────────────────────────────────

set terminal pngcairo size 900,600 enhanced font 'Helvetica,13'
set output 'resultados/distribucion.png'

set title 'Distribución de Primos por Proceso MPI (N = 10,000,000)' font ',15'
set xlabel 'Proceso MPI (rank)'
set ylabel 'Cantidad de Primos'

set style data histograms
set style histogram clustered gap 1
set style fill solid 0.80 border -1
set boxwidth 0.6

set yrange [0:*]
set xrange [-0.5:*]
set grid y
set key off

set format y '%g'

# Comentario de datos: column(4) = primos, xticlabels(1) = rank
plot 'resultados/resultados.dat' \
     using 4:xticlabels(1) \
     title 'Primos por proceso' \
     linecolor rgb '#4472C4'


# ─────────────────────────────────────────────────────────────
# GRÁFICA 2: Speedup real vs. ideal
# ─────────────────────────────────────────────────────────────
# Leemos speedup.dat con formato:
#   nprocs  tiempo  speedup
# Columna 1 = número de procesos (eje X)
# Columna 3 = speedup medido (eje Y)
# La línea ideal es simplemente y = x (speedup perfecto = n)
# ─────────────────────────────────────────────────────────────

set terminal pngcairo size 900,600 enhanced font 'Helvetica,13'
set output 'resultados/speedup.png'

set title 'Speedup — Criba de Eratóstenes con MPI' font ',15'
set xlabel 'Número de Procesos'
set ylabel 'Speedup'

set style data linespoints
unset style histogram

set xrange [0.5:4.5]
set yrange [0:5]
set xtics (1, 2, 4)
set grid

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
