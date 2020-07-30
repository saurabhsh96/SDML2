%Question 2, assignment 2
clear;
close all;

%% Input Definition

%Dimensions
h = 15e-3;
hs = 2.1e-3;

%EM
er = 12;
freq = 9e9:1e9:11e9;
c = 3e8;
lam = c./freq;
k0 = 2*pi./lam;
L = lam./2;
W = lam./20;

%Impedance (in Ohm)
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = (sqrt(mu_0/(eps_0*1))); 
zetaS = (sqrt(mu_0/(eps_0*er)));

%Current
M = [1,0,0];

%Meshgrid
% drad = pi/180;
% th = linspace(eps, pi/2-drad, 90);
% ph = linspace(eps, 2*pi-drad, 360);
% [thi, phi] = meshgrid(th, ph);
% dth = thi(1, 2) - thi(1, 1);
% dph = phi(2, 1) - phi(1, 1);


ph = (0:5:360)*pi/180;
th = linspace(eps,89,181)*pi/180;
dth = th(2)-th(1);
dph = ph(2)-ph(1);
[thi,phi] = meshgrid(th,ph);

%Observation point
r = 1;
z = r.*cos(thi);
z_dash = 0;

%% Q1-1 and Q1-2

%Allocating memory
Eth = zeros([size(k0, 2) size(thi)]);
Eph =  zeros([size(k0, 2) size(thi)]);
Emag =  zeros([size(k0, 2) size(thi)]);
Emax = zeros(size(k0));

%Free space
EthFS = zeros([size(k0, 2) size(thi)]);
EphFS =  zeros([size(k0, 2) size(thi)]);
EmagFS =  zeros([size(k0, 2) size(thi)]);
EmaxFS = zeros(size(k0));

%Plotting vec
th_vec = [-thi(1, size(thi, 2):-1:1), thi(1, :)].*180/pi;

for ind = 1:size(k0, 2)
    %Propagation constants
    ks = k0(ind).*sqrt(er);
    kxs = k0(ind).*sin(thi).*cos(phi);
    kys = k0(ind).*sin(thi).*sin(phi);
    kzs = k0(ind).*cos(thi);
    kRho = sqrt(kxs.^2 + kys.^2); 

    %Tx Line Equivalence
    [vTM, vTE, iTM, iTE] = trxline_SuperStrate(k0(ind), er, h, hs,...
        zeta0, zetaS, kRho, z);
    
    %JFT of the current
    %keq = (k0(ind) + ks)./2;
    MFT = CurrentFT(k0(ind), kxs, kys, L(ind), W(ind), M);

    %Calling SGF
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(k0(ind), 1, kxs, kys, vTM, ...
        vTE, iTM, iTE, zeta0, zetaS);

    %Calling Field function
    [Eth(ind,:,:), Eph(ind,:,:), Emag(ind,:,:), Emax(ind)] = Field(k0(ind), ...
        kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);

    %Free space Field
    [EFx, EFy, EFz] = FF(freq(ind), 1, L(ind), W(ind), r, thi, phi);
    
    %ErFS = (EFx.*sin(thi).*cos(phi)) + (EFy.*sin(thi).*sin(phi)) + (EFz.*cos(thi));
    EthFS(ind, :, :) = (EFx.*cos(thi).*cos(phi)) + ...
        (EFy.*cos(thi).*sin(phi)) - (EFz.*sin(thi));
    EphFS(ind, :, :) = (-EFx.*sin(phi)) + (EFy.*cos(phi));
    EmagFS(ind, :, :) = sqrt((abs(EFx)).^2 + (abs(EFy)).^2 + (abs(EFz)).^2); %Magnitude of E
    EmaxFS(ind) = max(max(EmagFS(ind,:,:)));
    EmagNormFS = squeeze((EmagFS(ind,:,:)./EmaxFS(ind)));
    EVecFS = [EmagNormFS(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1), EmagNormFS(1, :)];
    EVec1FS = [EmagNormFS(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1), EmagNormFS(91, :)];
        
    figure();
    EmagNorm = squeeze((Emag(ind,:,:)./Emax(ind)));
    EVec = [EmagNorm(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1) EmagNorm(1, :)];
    plot(th_vec, mag2db(EVec), 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
    EVec1 = [EmagNorm(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1) EmagNorm(91, :)];
    plot(th_vec, mag2db(EVec1), 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
    title(['Normalized Far-Field vs. \theta, Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    xlabel('\theta (deg)');
    ylabel('Normalized Far-Field E(\theta, \phi) (dB)');
    ylim([-35, 0]);
    legend show;
    hold off;
    grid on;
    
    figure(); %Comparison 0 deg
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], [EmagNormFS(1,size(EmagNormFS,2):-1:1), ...
        EmagNormFS(1,:)], 'LineWidth', 1.5, 'DisplayName', 'Free Space');
    hold on;
    %polarplot(EVec, 'LineWidth', 1.5, 'DisplayName', 'Spectral Field'); hold on;
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec...
       , 'LineWidth', 1.5, 'DisplayName', 'Spectral');
    title(['Comparison FF FS and Spectral \phi= 0 deg Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    legend show;
    hold off;
    
    figure(); %Comparison 90 deg
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], [EmagNormFS(91,size(EmagNormFS,2):-1:1), ...
        EmagNormFS(91,:)], 'LineWidth', 1.5, 'DisplayName', 'Free Space');
    hold on;
    %polarplot(EVec, 'LineWidth', 1.5, 'DisplayName', 'Spectral Field'); hold on;
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec1...
       , 'LineWidth', 1.5, 'DisplayName', 'Spectral');
    title(['Comparison FF FS and Spectral \phi= 90 deg Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    legend show;
    hold off;
    

    figure(); %FF
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec, 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec1, 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
    title(['Normalized Far-Field, Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    legend show;
    hold off;
end

%% Q2-3
%Change in directivity vs. Er

freq = 10e9;
lam = c./freq;
k0 = 2.*pi./lam;
L = lam./2;
W = lam./20;

kxs = k0.*sin(thi).*cos(phi);
kys = k0.*sin(thi).*sin(phi);
kzs = k0.*cos(thi);
kRho = sqrt(kxs.^2 + kys.^2); 

er = 1.2:0.2:25;
hs = lam./(4.*sqrt(er));

Dir = zeros([size(er, 2) size(thi)]);
DirReq = zeros(size(er));

for ind = 1:size(hs, 2)

    zetaS = (sqrt(mu_0/(eps_0*er(ind))));
    ks = k0.*sqrt(er(ind));
    
    %JFT of the current
    %keq = (k0 + ks)./2;
    MFT = CurrentFT(k0, kxs, kys, L, W, M);

    %Tx Line Equivalence
    [vTM, vTE, iTM, iTE] = trxline_SuperTest(k0, er(ind), h, hs(ind),...
        zeta0, zetaS, kRho, z);
    
    %Calling SGF
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(k0, 1, kxs, kys, vTM, ...
        vTE, iTM, iTE, zeta0, zetaS);

    %Calling Field function
    %[Eth_fs,Eph_fs] = farfied_fs(k01,r,th,ph,SGF,MFT);

    [Eth1, Eph1, Emag1, Emax1] = Field(k0, ...
        kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);

    %Spectral version of directivity function
    [Dir(ind,:,:), Prad] = DirectivityF(Emag1, 1, r, thi, dth, dph);
   
    DirReq(ind) = max(Dir(ind, :, 1));
end
figure();
plot(er, DirReq, 'LineWidth', 1.5);
title('Directivity vs. \epsilon_r (Linear Scale)');
xlabel('\epsilon_r');
ylabel('Directivity');
grid on;

figure();
plot(er, pow2db(DirReq), 'LineWidth', 1.5);
title('Directivity vs. \epsilon_r (Log Scale)');
xlabel('\epsilon_r');
ylabel('Directivity (dB)');
grid on;