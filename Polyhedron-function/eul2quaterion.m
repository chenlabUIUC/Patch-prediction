function quat = eul2quaterion(psi,theta,phi)

qs = cos(0.5*(phi+psi))*cos(theta/2);
qi = cos(0.5*(phi-psi))*sin(theta/2);
qj = sin(0.5*(phi-psi))*sin(theta/2);
qk = sin(0.5*(phi+psi))*cos(theta/2);
quat = [qs qi qj qk];
end