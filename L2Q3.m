%L2Q3
clear;
close all;

%% Defining Inputs

%Dim
h = 15e-3;

%EM
freq = 10e9;
c = 3e8;
lam = c/freq;
L = lam./2;
W = lam./20;
k0 = 2*pi./lam;

%Impedance (in Ohm)
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = (sqrt(mu_0/(eps_0*1))); 

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

%% Q3.1

er = 12;

%Propagation const
ks = k0.*sqrt(er);
kxs = ks.*sin(thi).*cos(phi);
kys = ks.*sin(thi).*sin(phi);
kzs = ks.*cos(thi);
kRho = sqrt(kxs.^2 + kys.^2); 
zetaS = (sqrt(mu_0/(eps_0*er)));

%Tx Line Equivalence
[vTM, vTE, iTM, iTE] = trxline_SuperStrate3(k0, ...
    er, h, zeta0, zetaS, kRho, z);

%JFT of the current
keq = (k0 + ks)./2;
MFT = CurrentFT(k0, kxs, kys, L, W, M); %k0 or ks

%Calling SGF
[Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(ks, er, kxs, kys, vTM, ...
    vTE, iTM, iTE, zeta0, zetaS);

%Calling Field function
[Eth, Eph, Emag, Emax] = Field(ks, ...
    kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);

%Plotting
% th_vec = [-thi(1, size(thi, 2):-1:1), thi(1, :)].*180/pi;
% 
% figure();
% EmagNorm = squeeze((Emag(:,:)./Emax));
% EVec = [(EmagNorm(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1))...
%     (EmagNorm(1, :))];
% plot(th_vec, mag2db(EVec), 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
% EVec1 = [(EmagNorm(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1))...
%     (EmagNorm(91, :))];
% plot(th_vec, mag2db(EVec1), 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
% title(['Normalized Far-Field vs. \theta, Freq = ', num2str(freq./10^9), ' GHz']);
% xlabel('\theta (deg)');
% ylabel('Normalized Far-Field E(\theta, \phi) (dB)');
% ylim([-50, 0]);
% legend show;
% hold off;
% 
% figure();
% polarplot(th_vec./180*pi, EVec, 'LineWidth', 1.5);

%% Q3.2
%Chaning er and finding directivity accordingly
er = 1.2:0.2:25;
%hs = lam./(4.*sqrt(er));

Dir = zeros([size(er, 2) size(thi)]);
DirReq = zeros(size(er));

for ind = 1:size(er, 2)

    zetaS = (sqrt(mu_0/(eps_0*er(ind))));
    ks = k0.*sqrt(er(ind));
    kxs = ks.*sin(thi).*cos(phi);
    kys = ks.*sin(thi).*sin(phi);
    kzs = ks.*cos(thi);
    kRho = sqrt(kxs.^2 + kys.^2); 

    [vTM, vTE, iTM, iTE] = trxline_Super_test3(k0,er(ind),h,kRho,z);
    
    %trxline_SuperStrate3(k0,er(ind), h, zeta0, zetaS, kRho, z);

    %JFT of the current
    %keq = (k0 + ks)./2;
    MFT = CurrentFT(k0, kxs, kys, L, W, M); %k0 or ks

    %Calling SGF
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(ks, er(ind), kxs, kys, vTM, ...
        vTE, iTM, iTE, zeta0, zetaS);

    %Calling Field function
    [Eth, Eph, Emag, Emax] = Field(ks, ...
        kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);

    %Spectral version of directivity function
    [Dir(ind,:,:), Prad] = DirectivityF(Emag, er(ind), r, thi, dth, dph);
   
    DirReq(ind) = max(Dir(ind, :, 1));
end
figure();
plot(er, DirReq, 'LineWidth', 1.5);
title('Directivity vs. \epsilon_r (Linear Scale)');
xlabel('\epsilon_r');
ylabel('Directivity');
grid on;

figure();
%plot(er, pow2db(abs(DirReq)), 'LineWidth', 1.5);
plot(er, pow2db(DirReq), 'LineWidth', 1.5);
title('Directivity vs. \epsilon_r (Log Scale)');
xlabel('\epsilon_r');
ylabel('Directivity (dB)');
grid on;

%% Q3-3
DirReqLast = load('N:\MASTERS\Quarter 4\Spectral Domain Techniques\Matlab\Lecture 2\Pictures\matlab.mat');

figure();
plot(er, DirReq, 'LineWidth', 1.5, 'DisplayName', 'Dir-Q3'); hold on;
plot(er, DirReqLast.DirReq, 'LineWidth', 1.5, 'DisplayName', 'Dir-Q2');
title('Comparison of Directivity variation of Q2 and Q3 (Linear Scale)');
xlabel('\epsilon_r');
ylabel('Directivity');
grid on;
legend show;

figure();
%plot(er, pow2db(abs(DirReq)), 'LineWidth', 1.5);
plot(er, pow2db(DirReq), 'LineWidth', 1.5, 'DisplayName', 'Dir-Q3'); hold on;
plot(er, pow2db(DirReqLast.DirReq), 'LineWidth', 1.5, 'DisplayName', 'Dir-Q2');
title('Comparison of Directivity variation of Q2 and Q3 (Log Scale)');
xlabel('\epsilon_r');
ylabel('Directivity (dB)');
grid on;
legend show;