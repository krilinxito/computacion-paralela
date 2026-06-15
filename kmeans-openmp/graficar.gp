# graficar.gp — Genera graficas del benchmark K-Means con OpenMP
# Uso: gnuplot graficar.gp
# Requiere: resultados/benchmark.dat (generado por ./kmeans)
#
# benchmark.dat tiene 4 columnas:
#   $1 = N (numero de puntos)
#   $2 = hilos
#   $3 = tiempo en segundos
#   $4 = speedup (T_1 / T_N)

# ===================================================================
# GRAFICA 1: Speedup vs Numero de Hilos (una linea por N)
# Muestra como escala el speedup segun el tamano del dataset.
# Para N pequeno el overhead domina; para N grande el speedup
# se acerca al numero de hilos.
# ===================================================================

set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "resultados/speedup.png"

set title "Speedup de K-Means con OpenMP\n(Dataset sintetico, D=16, K=10)" \
          font "Arial,14"
set xlabel "Numero de Hilos"       font "Arial,12"
set ylabel "Speedup (T_1 / T_N)"   font "Arial,12"

set xrange [0.5:5]
set yrange [0:4.5]
set xtics (1, 2, 4)
set ytics 0.5
set key top left
set grid lc rgb "#DDDDDD"
set border lw 1.5

# Colores: N=10k azul, 100k verde, 500k naranja, 1M rojo
# Filtrado: ($1==N ? $2 : 1/0) retorna col2 si N coincide, NaN si no.
# Gnuplot dibuja NaN como punto ausente, logrando una linea por N.

ideal(x) = x

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

print "Generado: speedup.png"

# ===================================================================
# GRAFICA 2: Tiempo de Ejecucion vs N (escala log-log)
# Muestra escalabilidad: como crece el tiempo al aumentar N,
# y cuanto ayuda cada numero de hilos.
# En escala log-log, relacion lineal = crecimiento O(N).
# ===================================================================

set output "resultados/scaling.png"

set title "Tiempo de Ejecucion vs Tamano del Dataset\n(Dataset sintetico, D=16, K=10)" \
          font "Arial,14"
set xlabel "N (numero de puntos)"  font "Arial,12"
set ylabel "Tiempo (segundos)"     font "Arial,12"

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

# Filtrado por numero de hilos con awk via pipe. Las filas de un mismo
# numero de hilos NO son contiguas en el archivo (estan intercaladas por
# N), lo que rompe las lineas si se filtra con el truco '1/0'. Con awk
# extraemos solo las filas del hilo deseado, quedando contiguas y la
# linea conecta limpio.
plot \
  "< awk '$2==1' resultados/benchmark.dat" using 1:3 \
      title "1 hilo"  with linespoints pt 7 ps 1.4 lw 2 lc rgb "#E74C3C", \
  "< awk '$2==2' resultados/benchmark.dat" using 1:3 \
      title "2 hilos" with linespoints pt 9 ps 1.4 lw 2 lc rgb "#3498DB", \
  "< awk '$2==4' resultados/benchmark.dat" using 1:3 \
      title "4 hilos" with linespoints pt 5 ps 1.4 lw 2 lc rgb "#2ECC71"

print "Generado: scaling.png"
print ""
print "Listo. Abre speedup.png y scaling.png para ver las graficas."
