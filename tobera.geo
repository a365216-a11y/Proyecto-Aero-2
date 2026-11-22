SetFactory("Built-in");
Mesh.MshFileVersion = 2.2;

DefineConstant[ r = {1, Name "r"} ];
DefineConstant[ withExterior = {1, Name "withExterior"} ];

AR = 1.6;
ht = 0.025;
he = AR*ht;
hi = 3.0*ht;
Lc = 0.25;
Ld = 0.20;
Lx = 1.00;
Ly = 0.30;
dz = 0.01;

nxc = 50*r;
nxd = 80*r;
ny  = 20*r;
nxj = 100*r;
nyo = 25*r;
pj  = 1.02^(1.0/r);
pu  = 1.05^(1.0/r);

Nsp = 40;

Point(1) = {-Lc, 0, 0, 1};
Point(2) = { 0,  0, 0, 1};
Point(3) = { Ld, 0, 0, 1};

For i In {0:Nsp}
  x = -Lc + Lc*i/Nsp;
  Point(100+i) = {x, ht + (hi-ht)*(1+Cos(Pi*(x+Lc)/Lc))/2, 0, 1};
EndFor

For i In {1:Nsp}
  x = Ld*i/Nsp;
  Point(200+i) = {x, ht + (he-ht)*(1-Cos(Pi*x/Ld))/2, 0, 1};
EndFor

Spline(1) = {100:100+Nsp};
Spline(2) = {100+Nsp, 201:200+Nsp};

Line(10) = {1, 2};
Line(11) = {2, 3};
Line(13) = {1, 100};
Line(14) = {2, 100+Nsp};
Line(15) = {3, 200+Nsp};

Curve Loop(1) = {10, 14, -1, -13};
Plane Surface(1) = {1};
Curve Loop(2) = {11, 15, -2, -14};
Plane Surface(2) = {2};

If (withExterior == 1)
  Point(4) = {Ld+Lx, 0,  0, 1};
  Point(5) = {Ld,    Ly, 0, 1};
  Point(6) = {Ld+Lx, Ly, 0, 1};
  Point(7) = {Ld+Lx, he, 0, 1};

  Line(12) = {3, 4};
  Line(16) = {200+Nsp, 7};
  Line(17) = {4, 7};
  Line(18) = {200+Nsp, 5};
  Line(19) = {5, 6};
  Line(20) = {7, 6};

  Curve Loop(3) = {12, 17, -16, -15};
  Plane Surface(3) = {3};
  Curve Loop(4) = {16, 20, -19, -18};
  Plane Surface(4) = {4};
EndIf

Transfinite Curve {1, 10}      = nxc+1 Using Progression 1;
Transfinite Curve {2, 11}      = nxd+1 Using Progression 1;
Transfinite Curve {13, 14, 15} = ny+1  Using Progression 1;
Transfinite Surface {1};
Transfinite Surface {2};
Recombine Surface {1, 2};

If (withExterior == 1)
  Transfinite Curve {12, 16, 19} = nxj+1 Using Progression pj;
  Transfinite Curve {17}         = ny+1  Using Progression 1;
  Transfinite Curve {18, 20}     = nyo+1 Using Progression pu;
  Transfinite Surface {3};
  Transfinite Surface {4};
  Recombine Surface {3, 4};
EndIf

out1[] = Extrude {0, 0, dz} { Surface{1}; Layers{1}; Recombine; };
out2[] = Extrude {0, 0, dz} { Surface{2}; Layers{1}; Recombine; };

vol[]    = {out1[1], out2[1]};
frente[] = {1, 2, out1[0], out2[0]};
sym[]    = {out1[2], out2[2]};

If (withExterior == 1)
  out3[] = Extrude {0, 0, dz} { Surface{3}; Layers{1}; Recombine; };
  out4[] = Extrude {0, 0, dz} { Surface{4}; Layers{1}; Recombine; };
  vol[]    += {out3[1], out4[1]};
  frente[] += {3, 4, out3[0], out4[0]};
  sym[]    += {out3[2]};
EndIf

Physical Surface("inlet")        = {out1[5]};
Physical Surface("wall")         = {out1[4], out2[4]};
Physical Surface("symmetry")     = sym[];
Physical Surface("frontAndBack") = frente[];
Physical Volume("fluido")        = vol[];

If (withExterior == 1)
  Physical Surface("plate")   = {out4[5]};
  Physical Surface("ambient") = {out3[3], out4[3], out4[4]};
Else
  Physical Surface("outlet")  = {out2[3]};
EndIf
