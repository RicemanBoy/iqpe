OPENQASM 3.0;
include "stdgates.inc";
bit[16] c;
qubit[33] q;
h q[0];
cx q[0], q[1];
u3(pi, -pi, 0) q[0];
h q[2];
cx q[2], q[3];
h q[5];
cx q[5], q[1];
cx q[5], q[6];
cx q[5], q[2];
h q[7];
cx q[7], q[6];
h q[8];
cx q[8], q[4];
u2(-pi, -pi) q[4];
cx q[8], q[9];
cx q[8], q[5];
u2(-pi, -pi) q[8];
h q[10];
cx q[10], q[9];
cx q[7], q[11];
h q[13];
cx q[13], q[12];
u2(-pi, -pi) q[12];
h q[15];
cx q[15], q[14];
cx q[10], q[14];
cx q[10], q[13];
cx q[7], q[10];
h q[16];
cx q[16], q[17];
u2(-pi, 0) q[16];
h q[18];
cx q[18], q[19];
z q[19];
h q[21];
cx q[21], q[17];
z q[17];
cx q[21], q[22];
cx q[21], q[18];
z q[18];
h q[23];
cx q[23], q[22];
h q[24];
cx q[24], q[20];
u2(-pi, -pi) q[20];
cx q[24], q[25];
cx q[24], q[21];
u2(-pi, -pi) q[24];
h q[26];
cx q[26], q[25];
cx q[23], q[27];
h q[29];
cx q[29], q[28];
u2(-pi, -pi) q[28];
h q[31];
cx q[31], q[30];
cx q[26], q[30];
cx q[26], q[29];
cx q[23], q[26];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[0];
  z q[4];
  z q[8];
  z q[12];
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
h q[1];
h q[2];
h q[3];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
cz q[8], q[24];
cz q[12], q[28];
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
cz q[0], q[16];
h q[16];
h q[17];
cz q[1], q[17];
h q[17];
h q[18];
cz q[2], q[18];
h q[18];
h q[19];
cz q[3], q[19];
h q[19];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
u2(0, 0) q[17];
cz q[1], q[17];
h q[1];
u2(-pi, -pi) q[17];
u2(0, 0) q[18];
cz q[2], q[18];
h q[2];
u2(-pi, -pi) q[18];
u2(0, 0) q[19];
cz q[3], q[19];
h q[3];
u2(-pi, -pi) q[19];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
z q[4];
cz q[8], q[24];
z q[8];
cz q[12], q[28];
z q[12];
u3(pi, -pi, 0) q[16];
cz q[0], q[16];
u2(0, 0) q[0];
u3(pi, -pi, 0) q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[0];
  z q[4];
  z q[8];
  z q[12];
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
h q[1];
h q[2];
h q[3];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
cz q[8], q[24];
cz q[12], q[28];
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
cz q[0], q[16];
h q[16];
h q[17];
cz q[1], q[17];
h q[17];
h q[18];
cz q[2], q[18];
h q[18];
h q[19];
cz q[3], q[19];
h q[19];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
u2(0, 0) q[17];
cz q[1], q[17];
h q[1];
u2(-pi, -pi) q[17];
u2(0, 0) q[18];
cz q[2], q[18];
h q[2];
u2(-pi, -pi) q[18];
u2(0, 0) q[19];
cz q[3], q[19];
h q[3];
u2(-pi, -pi) q[19];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
z q[4];
cz q[8], q[24];
z q[8];
cz q[12], q[28];
z q[12];
u3(pi, -pi, 0) q[16];
cz q[0], q[16];
u2(0, 0) q[0];
u3(pi, -pi, 0) q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[0];
  z q[4];
  z q[8];
  z q[12];
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
h q[1];
h q[2];
h q[3];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
cz q[8], q[24];
cz q[12], q[28];
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
cz q[0], q[16];
h q[16];
h q[17];
cz q[1], q[17];
h q[17];
h q[18];
cz q[2], q[18];
h q[18];
h q[19];
cz q[3], q[19];
h q[19];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
u2(0, 0) q[17];
cz q[1], q[17];
h q[1];
u2(-pi, -pi) q[17];
u2(0, 0) q[18];
cz q[2], q[18];
h q[2];
u2(-pi, -pi) q[18];
u2(0, 0) q[19];
cz q[3], q[19];
h q[3];
u2(-pi, -pi) q[19];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
z q[4];
cz q[8], q[24];
z q[8];
cz q[12], q[28];
z q[12];
u3(pi, -pi, 0) q[16];
cz q[0], q[16];
u2(0, 0) q[0];
u3(pi, -pi, 0) q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[0];
  z q[4];
  z q[8];
  z q[12];
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
h q[0];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[0], q[32];
cx q[1], q[32];
cx q[2], q[32];
cx q[3], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[1], q[32];
  cx q[2], q[32];
  cx q[3], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[1];
    z q[2];
    z q[3];
  }
}
h q[0];
h q[1];
h q[2];
h q[3];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[0], q[32];
cx q[4], q[32];
cx q[8], q[32];
cx q[12], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[0], q[32];
  cx q[4], q[32];
  cx q[8], q[32];
  cx q[12], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[0];
    z q[4];
    z q[8];
    z q[12];
  }
}
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
cz q[8], q[24];
cz q[12], q[28];
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
cz q[0], q[16];
h q[16];
h q[17];
cz q[1], q[17];
h q[17];
h q[18];
cz q[2], q[18];
h q[18];
h q[19];
cz q[3], q[19];
h q[19];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/2, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  z q[16];
  z q[17];
  z q[18];
  z q[19];
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
reset q[32];
u2(pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
reset q[32];
if (c[0]) {
  u2(pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
h q[16];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[17], q[32];
cx q[18], q[32];
cx q[19], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[17], q[32];
  cx q[18], q[32];
  cx q[19], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[17];
    z q[18];
    z q[19];
  }
}
h q[16];
u2(0, 0) q[17];
cz q[1], q[17];
h q[1];
h q[17];
u2(0, 0) q[18];
cz q[2], q[18];
h q[2];
h q[18];
u2(0, 0) q[19];
cz q[3], q[19];
h q[3];
h q[19];
reset q[32];
u2(-pi/4, -pi) q[32];
cx q[16], q[32];
cx q[20], q[32];
cx q[24], q[32];
cx q[28], q[32];
c[0] = measure q[32];
if (c[0]) {
  reset q[32];
  u2(-pi/2, -pi) q[32];
  cx q[16], q[32];
  cx q[20], q[32];
  cx q[24], q[32];
  cx q[28], q[32];
  c[0] = measure q[32];
  if (c[0]) {
    z q[16];
    z q[20];
    z q[24];
    z q[28];
  }
}
cz q[4], q[20];
h q[4];
cz q[8], q[24];
h q[8];
cz q[12], q[28];
h q[12];
u3(pi, -pi, 0) q[16];
cz q[0], q[16];
h q[0];
c[0] = measure q[0];
h q[16];
h q[20];
h q[24];
h q[28];
c[1] = measure q[1];
c[2] = measure q[2];
c[3] = measure q[3];
c[4] = measure q[4];
c[5] = measure q[5];
c[6] = measure q[6];
c[7] = measure q[7];
c[8] = measure q[8];
c[9] = measure q[9];
c[10] = measure q[10];
c[11] = measure q[11];
c[12] = measure q[12];
c[13] = measure q[13];
c[14] = measure q[14];
c[15] = measure q[15];
