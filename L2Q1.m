%Q1 Lecture 2: Far field of Ground Slab thing
clear;
close all;

%% Inputs

%Dimensions
L = 1e-3;
W = 1e-3;
h = 2e-3;

%EM
er = 10;
c = 3e8;

%Impedance (in Ohm)
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = (sqrt(mu_0/(eps_0*1))); 
zetaS = (sqrt(mu_0/(eps_0*er)));

%Current
J = [1,0,0];

%Meshgrid
drad = pi/180;
th = linspace(eps, pi/2-drad, 90);
ph = linspace(eps, 2*pi-drad, 360);
[thi, phi] = meshgrid(th, ph);
dth = thi(1, 2) - thi(1, 1);
dph = phi(2, 1) - phi(1, 1);

%% Q1-1

%EM
freq = 15e9;
lam = c/freq;
k0 = 2*pi/lam;
ks = k0*sqrt(er);

kxs = k0.*sin(thi).*cos(phi);
kys = k0.*sin(thi).*sin(phi);
kzs = k0.*cos(thi);
kRho = sqrt(kxs.^2 + kys.^2); 

%Observation point
r = 100*lam;
z = r.*cos(thi);
z_dash = h;
%z = r;
%z_dash = 0;

%Tx-line equivalent circuits
[vTM, vTE, iTM, iTE] = trxline_GroundSlab(k0, er, h, zeta0, zetaS, kRho, z);

%Green's function calculations
[Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SGFejF(k0, er, ...
    kxs, kys, vTM, vTE, iTM, iTE, zeta0, zetaS);

%JFT of the current
keq = (k0 + ks)./2;
JFT = CurrentFT(keq, kxs, kys, L, W, J);

%Calculating Electric field
[Eth, Eph, Emag, Emax] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
    Gzx, Gxy, Gyy, Gzy, JFT, z, z_dash);

%Plotting Electric Field
figure(1);
th_vec = [-thi(1, size(thi, 2):-1:1) thi(1, :)].*180/pi;
EmagNorm = mag2db(Emag./Emax);
EVec = [EmagNorm(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1) EmagNorm(1, :)];
plot(th_vec, EVec, 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
EVec = [EmagNorm(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1) EmagNorm(91, :)];
plot(th_vec, EVec, 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
title('Normalized Far-Field vs. \theta, for differenct \phi');
xlabel('\theta (deg)');
ylabel('Normalized Far-Field E(\theta, \phi) (dB)');
ylim([-30, 0]);
legend show;
grid on;

%% Q1-2

%Freq variation
fr = 1e9:0.5e9:25e9;
%Spectral
PradS = zeros(size(fr));
DirS = zeros([size(thi) size(fr, 2)]);
%FS
PradF = zeros(size(fr));
DirF = zeros([size(thi) size(fr, 2)]);

for ind = 1:size(PradS, 2)
    freq = fr(ind);
    lam = c/freq;
    k0 = 2*pi/lam;
    ks = k0*sqrt(er);

    kxs = k0.*sin(thi).*cos(phi);
    kys = k0.*sin(thi).*sin(phi);
    kzs = k0.*cos(thi);
    kRho = sqrt(kxs.^2 + kys.^2); 

    %Observation point
    r = 100*lam;
    z = r.*cos(thi);
    z_dash = h;
    %z = r;
    %z_dash = 0;

    %Tx-line equivalent circuits
    [vTM, vTE, iTM, iTE] = trxline_GroundSlab(k0, er, h, zeta0, zetaS, kRho, z);

    %Green's function calculations
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SGFejF(k0, er, ...
        kxs, kys, vTM, vTE, iTM, iTE, zeta0, zetaS);

    %JFT of the current
    keq = (k0 + ks)./2;
    JFT = CurrentFT(keq, kxs, kys, L, W, J);

    %Calculating Electric field
    [Eth, Eph, Emag, Emax] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
        Gzx, Gxy, Gyy, Gzy, JFT, z, z_dash);

    %Spectral version of directivity function
    [DirS(:,:,ind), PradS(1, ind)] = DirectivityF(Emag, 1, r, thi, dth, dph);
    
    %Free space -> Assuming no stratified medium at all, not even the
    %ground plane, As taking till theta = 90 deg, prad is doubled in below
    %calculations
    [DirF(:,:,ind), PradF(1, ind)] = Directivity(freq, L, W, 1, r, thi, phi);
    
    %Assuming ground plane is there, then no need of multiplying Prad by 2
    %in below calculations
    %[DirF(:,:,ind), PradF(1, ind)] = DirectivityH(freq, L, W, 1, r, thi, phi, h);
end

figure();
Prad = PradS./(2*PradF);
plot(fr./(10^9), Prad, 'LineWidth', 1.5);
title('P_{rad} vs. Frequency (Normalized) w.r.t. P_{rad} in FS');
xlabel('Frequency (in GHz)');
ylabel('Normalized P_{rad}');
grid on;

%% Q1-3
%Taking freq as 12.5 GHz as PradNorm is max there and then chaning the
%height of the substrate
h = eps:0.1e-3:15e-3;

freq = 12.5e9;
lam = c/freq;
k0 = 2*pi/lam;
ks = k0*sqrt(er);

kxs = k0.*sin(thi).*cos(phi);
kys = k0.*sin(thi).*sin(phi);
kzs = k0.*cos(thi);
kRho = sqrt(kxs.^2 + kys.^2); 

%Observation point
r = 100*lam;
z = r.*cos(thi);

%Spectral
PradS = zeros(size(h));
DirS = zeros([size(thi) size(h, 2)]);
%FS
PradF = zeros(size(h));
DirF = zeros([size(thi) size(h, 2)]);

for ind = 1:size(PradS, 2)
    %freq = fr(ind);
    z_dash = h(ind);

    %Tx-line equivalent circuits
    [vTM, vTE, iTM, iTE] = trxline_GroundSlab(k0, er, h(ind), zeta0, zetaS, kRho, z);

    %Green's function calculations
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SGFejF(k0, er, ...
        kxs, kys, vTM, vTE, iTM, iTE, zeta0, zetaS);

    %JFT of the current
    keq = (k0 + ks)./2;
    JFT = CurrentFT(keq, kxs, kys, L, W, J);

    %Calculating Electric field
    [Eth, Eph, Emag, Emax] = Field(k0, kzs, r, thi, phi, Gxx, Gyx,...
        Gzx, Gxy, Gyy, Gzy, JFT, z, z_dash);

    %Spectral version of directivity function
    [DirS(:,:,ind), PradS(1, ind)] = DirectivityF(Emag, 1, r, thi, dth, dph);
    
    %Free space -> Assuming no stratified medium at all, not even the
    %ground plane, As taking till theta = 90 deg, prad is doubled in below
    %calculations
    [DirF(:,:,ind), PradF(1, ind)] = Directivity(freq, L, W, 1, r, thi, phi);
    
    %Assuming ground plane is there, then no need of multiplying Prad by 2
    %in below calculations
    %[DirF(:,:,ind), PradF(1, ind)] = DirectivityH(freq, L, W, 1, r, thi, phi, h);
end

figure();
Prad = PradS./(2*PradF);
plot(h./(10^-3), Prad, 'LineWidth', 1.5);
title('P_{rad} vs. Height (Normalized) w.r.t. P_{rad} in FS');
xlabel('Height (in mm)');
ylabel('Normalized P_{rad}');
grid on;
