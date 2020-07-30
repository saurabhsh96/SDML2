%GEM
function [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(k0, er, kx, ky, vTM, vTE, iTM, iTE, zeta0, zetaS)
    zetaX = zeta0./sqrt(er);
    kRho = sqrt(kx.^2 + ky.^2);
    Gxx = (vTM - vTE).*kx.*ky./(kRho.^2); 
    Gyx = (vTE.*kx.^2 + vTM.*ky.^2)./(kRho.^2);
    Gzx = -zetaX.*ky.*iTM./k0;
    
    Gxy = -(vTE.*ky.^2 + vTM.*kx.^2)./(kRho.^2);
    Gyy = (vTM - vTE).*kx.*ky./(kRho.^2); 
    Gzy = -zetaX.*kx.*iTM./k0;
end