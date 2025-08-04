function q_ab = quat_multiply(a,b)

q0 = a(1);
q1 = a(2);
q2 = a(3);
q3 = a(4);

r0 = b(1);
r1 = b(2);
r2 = b(3);
r3 = b(4);

ab0 = r0*q0 - r1*q1 - r2*q2 - r3*q3;
ab1 = r0*q1 + r1*q0 - r2*q3 + r3*q2;
ab2 = r0*q2 + r1*q3 + r2*q0 - r3*q1;
ab3 = r0*q3 - r1*q2 + r2*q1 + r3*q0;
q_ab = [ab0 ab1 ab2 ab3];