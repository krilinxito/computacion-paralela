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

CC     = mpicc
CFLAGS = -O2 -Wall
TARGET = criba

# ── Compilar ──────────────────────────────────────────────────
all: $(TARGET)

$(TARGET): src/criba.c
	$(CC) $(CFLAGS) -o $(TARGET) src/criba.c -lm

# ── Ejecución rápida (4 procesos) ─────────────────────────────
run: $(TARGET)
	mpirun -np 4 ./$(TARGET)

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

# ── Generar gráficas PNG (requiere gnuplot) ───────────────────
graficas:
	gnuplot graficar.gp
	@echo "Generadas: distribucion.png  speedup.png"

# ── Compilar informe IEEE ─────────────────────────────────────
informe:
	cd informe && pdflatex informe_ieee.tex && pdflatex informe_ieee.tex

# ── Compilar referencias ──────────────────────────────────────
referencias:
	cd informe && pdflatex referencias.tex

# ── Limpiar archivos generados ────────────────────────────────
clean:
	rm -f $(TARGET)
	rm -f resultados/resultados.dat resultados/speedup.dat
	rm -f resultados/distribucion.png resultados/speedup.png
	rm -f informe/*.aux informe/*.log informe/*.out

.PHONY: all run benchmark graficas informe referencias clean
