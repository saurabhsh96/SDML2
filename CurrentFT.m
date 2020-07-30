% Routine to calculate the Fourier Transform of a current
function [currFT] = CurrentFT(k0, kl, kw, L, W, J)
    %kl represents the orientation of Length of dipole
    %kw represents the orientation of Width of dipole
    %Assuming PWS as current Dist.
    %Assuming I0 = 1
    %J = [1; 0; 0];
    Tx = sinc(kw.*W./2./pi); %Why divide by pi?
    Lx = (2.*k0.*(cos(kl.*L/2) - cos(k0.*L/2))./((k0.^2 - kl.^2).*sin(k0.*L/2)));
    
    currFT(1,:,:) = Lx.*Tx.*J(1);
    currFT(2,:,:) = Lx.*Tx.*J(2);
    currFT(3,:,:) = Lx.*Tx.*J(3);
    %someVal = (2*k0.*(cos(kl*L/2) - cos(k0*L/2)).*sinc(kw*W/2))./((k0^2 - kl.^2).*sin(k0*L/2));
    %CurrentFT(1,:,:) = someVal;
    %CurrentFT(2,:,:) = 0;
    %CurrentFT(3,:,:) = 0;
    %CurrentFT = (2*k0.*(cos(kl*L/2) - cos(k0*L/2)).*sinc(kw*W/2))./((k0^2 - kl.^2).*sin(k0*L/2));
end