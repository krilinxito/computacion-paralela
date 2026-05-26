# graficar.gp — Genera graficas del experimento K-Means con Iris
# Uso: gnuplot graficar.gp
# Requiere: resultados.dat y speedup.dat (generados por ./kmeans)

# ===================================================================
# GRAFICA 1: Scatter plot de clusters
# Muestra cada flor del dataset Iris como un punto coloreado segun
# el cluster asignado por K-Means
# Ejes: petal_length (col 1) vs petal_width (col 2)
# ===================================================================

set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "resultados/clusters.png"

set title "K-Means Clustering — Dataset Iris (K=3)" font "Arial,14"
set xlabel "Longitud del Petalo (cm)" font "Arial,12"
set ylabel "Ancho del Petalo (cm)"    font "Arial,12"
set key top left
set grid lc rgb "#DDDDDD"
set border lw 1.5

# Tecnica de filtrado por cluster:
# ($3==N ? $1 : 1/0) retorna el valor de la columna 1 si el cluster
# es N, o 1/0 (=NaN) que gnuplot dibuja como punto ausente.
# Esto permite graficar cada cluster con un color distinto.

plot \
  'resultados/resultados.dat' \
      using ($3==1 ? $1 : 1/0):($3==1 ? $2 : 1/0) \
      title "Cluster 1 (Setosa)"      pt 7 ps 1.3 lc rgb "#E74C3C", \
  'resultados/resultados.dat' \
      using ($3==2 ? $1 : 1/0):($3==2 ? $2 : 1/0) \
      title "Cluster 2 (Versicolor)"  pt 7 ps 1.3 lc rgb "#3498DB", \
  'resultados/resultados.dat' \
      using ($3==3 ? $1 : 1/0):($3==3 ? $2 : 1/0) \
      title "Cluster 3 (Virginica)"   pt 7 ps 1.3 lc rgb "#2ECC71"

print "Generado: clusters.png"

# ===================================================================
# GRAFICA 2: Speedup real vs ideal
# Muestra cuanto mas rapido corre el programa con mas hilos en
# comparacion con la version de 1 hilo (T_1 / T_N)
# La linea ideal seria speedup = num_hilos (perfecto escalado lineal)
# ===================================================================

set output "resultados/speedup.png"

set title "Speedup de K-Means con OpenMP — Dataset Iris" font "Arial,14"
set xlabel "Numero de Hilos"  font "Arial,12"
set ylabel "Speedup (T_1 / T_N)" font "Arial,12"

set xrange [0.5:5]
set yrange [0:3]
set xtics (1, 2, 4)
set ytics 0.5
set key top left
set grid lc rgb "#DDDDDD"

# Funcion de speedup ideal (lineal perfecto)
ideal(x) = x

plot \
  ideal(x) \
      title "Speedup ideal"   lw 2 lc rgb "#AAAAAA" dt 2, \
  'resultados/speedup.dat' \
      using 1:2 \
      title "Speedup real"    with linespoints \
      pt 7 ps 1.5 lw 2 lc rgb "#E74C3C"

print "Generado: speedup.png"
print ""
print "Listo. Abre clusters.png y speedup.png para ver las graficas."
