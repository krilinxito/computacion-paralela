CC     = gcc
CFLAGS = -O2 -fopenmp
LIBS   = -lm
TARGET = kmeans

# -lm debe ir DESPUES del fuente: GCC resuelve librerias de izq a der
all: $(TARGET)

$(TARGET): src/kmeans.c
	$(CC) $(CFLAGS) -o $(TARGET) src/kmeans.c $(LIBS)

run: $(TARGET)
	./$(TARGET)

# Genera las graficas (requiere gnuplot y que ./kmeans ya se haya ejecutado)
plots: resultados/benchmark.dat
	gnuplot graficar.gp

# Ejecuta el programa y luego genera las graficas en un solo paso
all-plots: run plots

resultados/benchmark.dat: $(TARGET)
	./$(TARGET)

# ── Compilar informe IEEE ────────────────────────────────────
informe:
	cd informe && pdflatex informe_ieee.tex && pdflatex informe_ieee.tex

referencias:
	cd informe && pdflatex referencias.tex

clean:
	rm -f $(TARGET)
	rm -f resultados/benchmark.dat
	rm -f resultados/speedup.png resultados/scaling.png
	rm -f informe/*.aux informe/*.log informe/*.out

.PHONY: all run plots all-plots clean informe referencias
