function [Eth_fs,Eph_fs] = farfied_fs(k01,r,th,ph,SGF,Jx1)

    kz = k01.*cos(th);

    Exx = 1j.*kz.*(squeeze(SGF(1,1,:,:)).*(Jx1)).*(exp(-1j*k01*r)./(2*pi*r));
    
    Eyy = 1j.*kz.*(squeeze(SGF(2,1,:,:)).*(Jx1)).*(exp(-1j*k01*r)./(2*pi*r));
    
    Ezz = 1j.*kz.*(squeeze(SGF(3,1,:,:)).*(Jx1)).*(exp(-1j*k01*r)./(2*pi*r));
    
    %Spherical coordinates
    Er_fs = (Exx.*sin(th).*cos(ph)) + (Eyy.*sin(th).*sin(ph)) + (Ezz.*cos(th));
    Eth_fs = (Exx.*cos(th).*cos(ph)) + (Eyy.*cos(th).*sin(ph)) - (Ezz.*sin(th));
    Eph_fs = (-Exx.*sin(ph)) + (Eyy.*cos(ph));
end
