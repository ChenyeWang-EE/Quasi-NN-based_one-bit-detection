function hHeadBLS = Bussgang_LeastSquare_ISI(xk,z,antennaNum)

x1 = xk(2:end);
x2 = xk(1:end-1);
x = [x1; x2];

[~,xLen] =size(x);
zvec = z(:);
% noiseVec = noise(:);
s = kron(x.',eye(antennaNum));
ss = diag( kron(x.'*conj(x),eye(antennaNum))+eye(antennaNum*xLen) ).^(-1/2) ;
A = sqrt(2/pi).*real(diag(ss));
ATilde = A*s;
hHeadBLStmp =  (ATilde'*ATilde)\ATilde'*zvec;

hHeadBLS = reshape(hHeadBLStmp,[antennaNum,length(hHeadBLStmp)/antennaNum]);
end