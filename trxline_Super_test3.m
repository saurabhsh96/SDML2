function [vTM, vTE, iTM, iTE] = trxline_Super_test3(k0,er,h,kro,z)

%Characteristic impedance
zeta = 120*pi;
zetas = zeta/sqrt(er);

% Propagation along z
ks = sqrt(er)*k0;
kz = -1j*sqrt(-((k0^2)-(kro.^2)));
kzs = -1j*sqrt(-((ks^2)-(kro.^2)));

%Impedances:
Z0TE = (zeta.*k0)./kz;
Z0TM = (zeta.*kz)./k0;

ZsTE = (zetas.*ks)./kzs;
ZsTM = (zetas.*kzs)./ks;

TauTE = ((ZsTE - Z0TE))./(ZsTE + Z0TE); 
TauTM = ((ZsTM - Z0TM))./(ZsTM + Z0TM);

V01TE = 1./(1+TauTE.*exp(-2i.*kz.*h));
V01TM = 1./(1+TauTM.*exp(-2i.*kz.*h));

v0TE = V01TE.*exp(-1i.*kz.*h).*(1+TauTE); % (TauTE.*exp(1i.*kz.*h)).*exp(1i.*kzs.*h).*exp(-1i.*kzs.*z)

vTE = v0TE.*exp(-1i.*kzs.*(z-h));
iTE = vTE./ZsTE;

v0TM = V01TM.*exp(-1i.*kz.*h).*(1+TauTM);

vTM = v0TM.*exp(-1i.*kzs.*(z-h));
iTM = vTM./ZsTM;
end