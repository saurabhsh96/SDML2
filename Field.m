%Field calculations
function [Eth, Eph, Emag, Emax] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
    Gzx, Gxy, Gyy, Gzy, J, z, z_dash)
    %Finding Electric field; Assuming obs pi in far field
    EFx = 1j.*kzs.*(Gxx.*squeeze(J(1,:,:)) ...
        + Gxy.*squeeze(J(2,:,:))).*(exp(1j.*kzs.*(z - z_dash))).*(exp(-1j*k0*r)./(2*pi*r));
    
    EFy = 1j.*kzs.*(Gyx.*squeeze(J(1,:,:)) ...
        + Gyy.*squeeze(J(2,:,:))).*(exp(1j.*kzs.*(z - z_dash))).*(exp(-1j*k0*r)./(2*pi*r));
    
    EFz = 1j.*kzs.*(Gzx.*squeeze(J(1,:,:)) ...
        + Gzy.*squeeze(J(2,:,:))).*(exp(1j.*kzs.*(z - z_dash))).*(exp(-1j*k0*r)./(2*pi*r));
    
    %Getting Spherical Coordinates
    Er = (EFx.*sin(thi).*cos(phi)) + (EFy.*sin(thi).*sin(phi)) + (EFz.*cos(thi));
    Eth = (EFx.*cos(thi).*cos(phi)) + (EFy.*cos(thi).*sin(phi)) - (EFz.*sin(thi));
    Eph = (-EFx.*sin(phi)) + (EFy.*cos(phi));
    
    Emag = sqrt(abs(EFx).^2 + abs(EFy).^2 + abs(EFz).^2); %Magnitude of E
    Emax = max(max(Emag));
end