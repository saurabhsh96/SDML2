%TRX Superstrate
function [vTM, vTE, iTM, iTE] = trxline_SuperTest(k0, er, h, hs, zeta0, zetaS, kRho, z)

    %In Slab
    ks = sqrt(er)*k0;
    
    %kZ
    kz0 = -1j*sqrt(-((k0^2)-(kRho.^2)));
    kzs = -1j*sqrt(-((ks^2)-(kRho.^2)));

    %Air TE TM impedance
    Z0TE = (zeta0.*k0)./kz0;
    Z0TM = (zeta0.*kz0)./k0;

    %Slab TE TM impedance
    ZsTE = (zetaS.*ks)./kzs;
    ZsTM = (zetaS.*kzs)./ks;
    
    %ZL from voltage
    ZlTM = ZsTM.*(Z0TM + 1j.*ZsTM.*tan(kzs*hs))./...
        (ZsTM + 1j.*Z0TM.*tan(kzs*hs));
    ZlTE = ZsTE.*(Z0TE + 1j.*ZsTE.*tan(kzs*hs))./...
        (ZsTE + 1j.*Z0TE.*tan(kzs*hs));
    
    %Gamma
    gamma1TE = (ZlTE - Z0TE)./(ZlTE + Z0TE);
    gamma1TM = (ZlTM - Z0TM)./(ZlTM + Z0TM);
    
    gamma2TE = (Z0TE - ZsTE)./(Z0TE + ZsTE);
    gamma2TM = (Z0TM - ZsTM)./(Z0TM + ZsTM);
    
    %Constants
    %Substrate and Air
    V0_pTM = 1./(1+gamma1TM.*exp(-2*1j*kz0*(h)));
    V0_pTE = 1./(1+gamma1TE.*exp(-2*1j*kz0*(h)));
    
    V0TM = V0_pTM.*(exp(-1j.*kz0.*h).*(1+gamma1TM));
    V0TE = V0_pTE.*(exp(-1j.*kz0.*h).*(1+gamma1TE));
    
    VsPTM = V0TM./(1+gamma2TM.*exp(-2*1j*kzs*(hs)));
    VsPTE = V0TE./(1+gamma2TE.*exp(-2*1j*kzs*(hs)));
    
    VsTM = VsPTM.*(exp(-1j.*kzs.*hs).*(1+gamma2TM));
    VsTE = VsPTE.*(exp(-1j.*kzs.*hs).*(1+gamma2TE));
    
    vTE = VsTE.*exp(-1j.*kz0.*(z-(h+hs)));
    vTM = VsTM.*exp(-1j.*kz0.*(z-(h+hs)));
    iTE = vTE./Z0TE;
    iTM = vTM./Z0TM;
end