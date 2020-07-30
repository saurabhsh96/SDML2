%Function to calculate the directivity of the antenna using the E-Field
function [Dir, Prad] = DirectivityF(Emag, er, r, thi, dth, dph)
    %Finding zeta
    %eps_0 = 8.854187817e-12;
    %mu_0 = 1.2566370614e-6;
    zeta0 = 120*pi;
    zetaS = zeta0./sqrt(er);
        
    %Calculate field intensity
    Intensity = r^2.*(abs(Emag).^2)./(2*zetaS);
    
    %Calculate prad
    Prad = (sum(Intensity.*sin(thi), 'all'))*dth*dph;
    
    %Calculate directivity
    Dir = (4*pi).*Intensity/Prad; 
end