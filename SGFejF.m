%Function to find SFG EJ
function [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SGFejF(k0, er, ...
    kx, ky, vTM, vTE, iTM, iTE, zeta0, zetaS)
    kRho = sqrt(kx.^2 + ky.^2);
    Gxx = -(((vTM.*(kx.^2)) + (vTE.*(ky.^2)))./(kRho.^2));
    Gyx = ((vTE - vTM).*kx.*ky)./(kRho.^2);
    Gzx = (zeta0.*iTM.*kx)./k0;   
    Gyy = -(((vTE.*(kx.^2)) + (vTM.*(ky.^2)))./(kRho.^2));
    Gxy = ((vTE - vTM).*kx.*ky)./(kRho.^2);
    Gzy = (zeta0.*iTM.*ky)./k0;
end