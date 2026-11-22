# Proyecto final: tobera convergente-divergente 2D

Aerodinámica II (AE604), UACH. Sandra Valeria Tarango Meléndez (365216) y Karyme Hernández Terán (360725), grupo 9A1.

Simulamos el flujo compresible no viscoso en una tobera plana con relación de áreas Ae/At = 1.6, con Gmsh y OpenFOAM-8 (rhoCentralFoam), y lo comparamos con la teoría cuasiunidimensional usando el área real A(x) de nuestra geometría. Los casos son tres: subsónico sin estrangulamiento (S), supersónico adaptado (A) y con choque normal dentro del divergente (C).

## Qué hay en cada carpeta

| Carpeta o archivo | Qué contiene |
| --- | --- |
| `reporte_final.pdf` y `reporte_final/` | El reporte en PDF y su código en LaTeX (con sus figuras y tablas) |
| `memoria_calculo/` | La memoria de cálculo analítico: ecuaciones, ramas, Newton y bisección, con su código en Python |
| `tobera.geo` | La geometría y la malla en Gmsh. El parámetro `r` cambia el refinamiento (r = 1, 2, 4) |
| `calidad_malla/` | La salida de `checkMesh` de las tres mallas |
| `casos/` | Los 29 casos de OpenFOAM de la tabla de corridas |
| `datos/` | Lo que sacamos de cada corrida: perfiles por sección (`axis_*.npz`), historias de gasto de entrada y salida y presión de salida (`.csv`) |
| `scripts/` | La teoría, el posproceso, las figuras y las tablas |
| `figuras/` y `tablas/` | Todas las figuras y tablas del reporte, y la tabla de corridas |
| `tabla_contribuciones.md` | Qué hizo cada una |

El Anexo B del reporte dice, para cada figura y tabla, de qué caso sale, con qué ventana de tiempo, con qué script y de qué archivo de datos. Ahí está la trazabilidad completa, por eso no la repetimos aquí.

## Versiones

- Gmsh 4.15.2
- OpenFOAM-8 (blueCFD-Core en Windows)
- Python 3 con NumPy, SciPy y Matplotlib

## Cómo se arma un caso

1. Generar la malla (en la terminal de blueCFD), con `r` = 1, 2 o 4:
   ```
   gmsh tobera.geo -3 -format msh2 -setnumber r 1 -setnumber withExterior 0 -o tobera_r1.msh
   ```
2. Dentro de la carpeta del caso, convertirla a OpenFOAM:
   ```
   gmshToFoam tobera_r1.msh
   ```
   Después hay que editar `constant/polyMesh/boundary`: el parche `symmetry` debe ser `symmetryPlane` y `frontAndBack` debe ser `empty`. `inlet`, `outlet` y `wall` se quedan como `patch` (la pared es deslizante porque en `0/U` usa `slip`).
3. Revisar la malla con `checkMesh`.
4. Correr en serie: `rhoCentralFoam > log 2>&1`
5. Sacar los datos y hacer figuras y tablas con los scripts de `scripts/` (`resultados.py` hace casi todo; `dtstudy.py` y `orden_courant.py` son del estudio del paso de tiempo; `auditoria.py` compara ventanas consecutivas).

Lo que cambia entre casos está en `0/p` (la presión de salida, o `zeroGradient` en A), en el tamaño de la malla y en `system/controlDict` (el Courant máximo y el tiempo final).

## Lo que usamos en los casos

- Aire como gas ideal: γ = 1.4, R = 287 J/(kg K). Sin viscosidad.
- Entrada con p0 = 300 kPa y T0 = 300 K. Garganta de 0.025 m de alto, salida de 0.04 m, entrada de 0.075 m. Convergente de 0.25 m y divergente de 0.20 m. Profundidad de 0.01 m.
- Simulamos solo la mitad (plano de simetría) y una celda en z.
- Esquema de Kurganov, limitador van Leer, Euler explícito y paso de tiempo adaptable con Courant máximo de 0.1 en las corridas finales.

Contrapresiones: S con 276 kPa, C con 200 kPa. En A la salida es supersónica y no le imponemos presión.

| Caso | Mallas | Tiempo final | Ventana de promedio |
| --- | --- | --- | --- |
| S | r1, r2 | 0.10 s (r1), 0.05 s (r2) | 0.04–0.10 s, 0.03–0.05 s |
| A | r1, r2 | 0.06 s (r1), 0.04 s (r2) | 0.03–0.06 s, 0.02–0.04 s |
| C | r1, r2, r4 | 0.06 s, 0.04 s, 0.02 s | 0.03–0.06 s, 0.02–0.04 s, 0.012–0.02 s |

Las mallas r1, r2 y r4 tienen 2600, 10 400 y 41 600 celdas.

## Cómo decidimos que una corrida ya convergió

Con choques los residuales no llegan a cero, así que usamos cantidades físicas promediadas por ventanas de tiempo. Pedimos tres cosas: que los promedios cambien menos de 1 % entre ventanas consecutivas, que el gasto de entrada y el de salida coincidan dentro de 1 %, y que en C el choque se mueva menos de 0.5 % de Ld.

## Cosas que no salieron a la primera (y que están en las carpetas)

- **El paso de tiempo importó más que la malla.** Las primeras corridas las hicimos con Courant 0.3 a 0.5 y salieron errores grandes: en S, hasta −13 % en el gasto, y en C el choque quedaba 6 a 9 % de Ld más arriba de lo que da la teoría. Con refinar la malla no mejoraba. Al bajar a Courant 0.1 los errores quedaron en décimas de por ciento. Por eso hay carpetas con sufijos `_co01`, `_dt_...` y las originales con Courant alto: las dejamos como evidencia del estudio.
- **Cambiamos la contrapresión de S de 290 a 276 kPa.** A 290 kPa el Mach en la salida es de 0.22 y el flujo es tan lento que el error numérico se amplifica mucho. Las carpetas `p273`, `p280` y la de 290 kPa son de esa decisión.
- **Pruebas de otros esquemas y de la salida** (`minmod`, `tadmor`, `wt005`, `wt020`, `wt045`) para ver si algo distinto al Courant arreglaba S. Ninguno mejoró lo que el Courant 0.1.
- **S desde un campo inicial uniforme** (`caso_dt_S_r1_unif`) para comprobar que el resultado no depende de la condición inicial.

## Limitaciones

- La ejecución en paralelo con MS-MPI de blueCFD nos falló (error 1726), así que todas las corridas fueron en serie.
- En la malla r4 la ventana de promedio es más corta que en r1 y r2 por el tiempo de cómputo. Lo discutimos en el reporte.
- Los resultados son de un modelo no viscoso en 2D y se comparan con teoría ideal, no con una tobera real.
- Por tamaño, cada carpeta de `casos/` guarda `0/`, `constant/`, `system/`, el último instante, `postProcessing/` y los logs. Los instantes intermedios no se incluyen; lo que se extrajo de ellos está en `datos/`. Los logs muy grandes van comprimidos (`.gz`).
- No hicimos el dominio exterior (descarga a la atmósfera), solo el interior.
