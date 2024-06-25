OPENQASM 2.0;
include "qelib1.inc";

qreg q[3];
cx q[1],q[0];
rz(0.25*pi) q[0];
cx q[1],q[0];
