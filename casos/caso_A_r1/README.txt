CASO A (supersonico adaptado) - tobera sola, malla r1
Hecho con OpenFOAM-8 (blueCFD) y Gmsh 4.15.2. Solver: rhoCentralFoam.

Modelo: adiabatico e inviscido (mu = 0, pared slip). Entrada: p0 = 300 kPa, T0 = 300 K.
Salida: libre (zeroGradient), sin imponer presion estatica. Aire, gamma = 1.4, R = 287.
Malla: tobera_interna_r1.msh (grupos: inlet, outlet, wall, symmetry, frontAndBack).
Teoria: pe/p01 = 0.141 -> pe = 42.3 kPa, Me = 1.935.

Pasos (en la terminal de blueCFD):
  cd /c/tobera/caso_A_r1
  checkMesh > log.checkMesh
  decomposePar > log.decomposePar
  mpirun -np 4 rhoCentralFoam -parallel > log.A_r1 2>&1
  reconstructPar > log.reconstructPar
Resultados de gasto/presion de salida: postProcessing/
