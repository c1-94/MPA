close all;
clc

N = 5e6;    % number density
v0 = 500e3; % bulk speed
mi = 1.672649E-27;  % ion mass
me = 9.11e-31;  % electron mass
m = mi; % christmas
q = 1.602176565E-19;    % elementary charge
kb = 1.38064852e-23;    % boltzman constant
T = 100000; % 10⁵K
v_step = 1
v = 1:v_step:2.5e6; % velocity range
vt = (2*kb*T/m)^(1/2);  % thermal velocity

r = 1e-3;
a = 2*r;
l = 1e-3;
pixel_size = l;
A = pi*r^2;  % aperture area
Phi0 = 360;  % azymuthal aperture angle
Theta0 = 5;    % elevation aperture angle
tau = 5e-3; % acquisition time

y_maxwell = v.^3.*exp(-(v-v0).^2/vt^2); % velocity part of the spherical 3D distribution
y_1D_maxwell = v.*exp(-(v-v0).^2/vt^2);  % velocity part of the 1D distribution
kappa = 6;
y_kappa = v.^3.*(1+(2/(2*kappa-3)).*((v-v0).^2/vt^2)).^(-kappa-1);
K_kappa = (1/vt^3)*(2/(pi*(2*kappa-3)))^(3/2)*gamma(kappa+1)/gamma(kappa-1/2);
S_maxwell = 0;
S_1D_maxwell = 0;
S_kappa = 0;
for i=1:length(v)-1
    S_maxwell = S_maxwell + (v(i+1)-v(i))*(y_maxwell(i+1)+y_maxwell(i))/2;  % integration over velocities
    S_1D_maxwell = S_1D_maxwell + (v(i+1)-v(i))*(y_1D_maxwell(i+1)+y_1D_maxwell(i))/2;  % integration over velocities
    S_kappa = S_kappa + (v(i+1)-v(i))*(y_kappa(i+1)+y_kappa(i))/2;  % integration over velocities
end

C_maxwell = N*tau*A*(m/(2*pi*kb*T))^(3/2)*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180))*S_maxwell     % total number of counts over all velocities
C_1D_maxwell = N*tau*A*(m/(2*pi*kb*T))^(1/2)*S_1D_maxwell     % total number of counts over all velocities 1D
Cv_maxwell = tau*A*(N/(pi^(3/2)*vt^3))*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180)).*y_maxwell;  % number of counts as a function of velocity
Cv_1D_maxwell = N*tau*A*v*(m/(2*pi*kb*T))^(1/2).*exp(-(v-v0).^2/vt^2);   % number of counts as a function of velocity 1D
C_kappa = N*tau*A*K_kappa*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180))*S_kappa
Cv_kappa = N*tau*A*K_kappa*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180)).*y_kappa;

figure()
semilogy(v, Cv_maxwell);
hold on;
semilogy(v, Cv_1D_maxwell);
semilogy(v, Cv_kappa);
legend("C(v) maxwellian", "C(v) 1D maxwellian", "C(v) kappa");
xlim([0.3 0.7].*1e6);
xlabel("velocity (m/s)");
ylabel("number of particles");
title("Number of particles as a function of velocity");

B0 = 0.1;
distance = 2*v./(q*B0/m);
step = distance(2) - distance(1);
S_maxwell = 0;
S_kappa = 0;
Count_pixel = zeros(2,length(distance));    % maxwell | kappa
floor(pixel_size/step)
for i=1:floor(pixel_size/step):length(distance)-floor(pixel_size/step)
    for j=i:i+floor(pixel_size/step)
        S_maxwell = S_maxwell + Cv_maxwell(j);
        S_kappa = S_kappa + Cv_kappa(j);
    end
    Count_pixel(1,i) = S_maxwell*v_step;
    Count_pixel(2,i) = S_kappa*v_step;
    S_maxwell = 0;
    S_kappa = 0;
end

Count_pixel(Count_pixel == 0) = 1e-10;

figure()
subplot(211)
semilogy(distance, Count_pixel(1,:));
xlabel("hit distance (m)");
ylabel("number of counts per pixel");
legend("Maxwellian");
ylim([1e-4 max(Count_pixel(1,:))*1.5]);
xlim([0 0.3]);
title("Number of counts per pixel Maxwellian distribution");
subplot(212)
semilogy(distance, Count_pixel(2,:));
xlabel("hit distance (m)");
ylabel("number of counts per pixel");
legend("Kappa");
ylim([1e-4 max(Count_pixel(2,:))*1.5]);
xlim([0 0.3]);
title("Number of counts per pixel Kappa distribution");

Counts_maxwell = 0;
Counts_kappa = 0;
v_counts_maxwell = 0;
v_counts_kappa = 0;
j = 1;
k = 1;
for i=1:length(Count_pixel(1,:))
    if (Count_pixel(1,i) >= 1)
        Counts_maxwell(j) = Count_pixel(1,i);
        v_counts_maxwell(j) = v(i);
        j = j+1;
    end
    if (Count_pixel(2,i) >= 1)
        Counts_kappa(k) = Count_pixel(2,i);
        v_counts_kappa(k) = v(i);
        k = k+1;
    end
end
distance_counts = 2*v_counts_maxwell./(q*B0/m);
figure()
hAx = subplot(211);
e = errorbar(distance_counts, Counts_maxwell, sqrt(Counts_maxwell)./2);
hAx.YScale = 'log';
hAx.YScale = 'log';
xlim([0 0.2]);
%xlim([0 0.15]);
ylim([1 max(Counts_maxwell)*2]);
e.CapSize = 10;
%ylim([1e-4 max(Count_pixel(1,:))*1.5]);
title("Number of particles as a function of hit distance");
xlabel("hit distance (m)");
ylabel("number of counts");
subplot(212);
semilogy(distance_counts, 1-(Counts_maxwell - sqrt(Counts_maxwell))./Counts_maxwell);
hold on;
ylim([0 1]);
xlim([0 0.2]);
line([min(xlim);max(xlim)],[0.1;0.1], 'color', 'g', 'LineWidth', 1);
line([min(xlim);max(xlim)],[0.2;0.2], 'color', 'r', 'LineWidth', 1);
title("relative error in percent");
xlabel("hit distance (m)");
ylabel("e_r (%)");

yneg = sqrt(Counts_maxwell);
ypos = sqrt(Counts_maxwell);
xneg = zeros(1,j-1);
xpos = zeros(1,j-1);
for i=1:j-1
    xneg(i) = (r+l)*(q*B0)/(4*m);
    xpos(i) = (r+l)*(q*B0)/(4*m);
end
figure()
hold on;
hAx = subplot(211);
% semilogy(v_counts_kappa, Counts_kappa, 'g', 'LineWidth', 2)
%e = errorbar(v_counts_maxwell, Counts_maxwell, yneg, ypos, xneg, xpos, 'LineStyle','none', 'Color', 'k','linewidth', 1.5);
e = errorbar(v_counts_maxwell, Counts_maxwell, yneg, ypos, 'LineStyle','none', 'Color', 'k','linewidth', 1.5);
hAx.YScale = 'log';
set(hAx, 'FontSize', 14)
hold on;
semilogy(v_counts_maxwell, Counts_maxwell, 'b', 'LineWidth', 2)
%xlim([0 0.15]);
%ylim([1 max(Counts)*2]);
e.CapSize = 10;
%ylim([1e-4 max(Count_pixel(1,:))*1.5]);
title("Number of counts as a function of velocity");
xlabel("velocity (m/s)");
ylabel("number of counts");
ylim([0.5*min(Counts_maxwell) 1.5*max(Counts_maxwell)]);
xlim([0.98*min(v_counts_maxwell) 1.02*max(v_counts_maxwell)]);
fig = subplot(212);
set(fig, 'FontSize', 14);
hold on;
semilogy(v_counts_maxwell, (1-(Counts_maxwell - sqrt(Counts_maxwell))./Counts_maxwell).*100, 'b', 'LineWidth', 2);
ylim([0 100]);
xlim([0.98*min(v_counts_maxwell) 1.02*max(v_counts_maxwell)]);
%xlim([0.06 0.15]);
line([min(xlim);max(xlim)],[10;10], 'color', 'g', 'LineWidth', 1);
line([min(xlim);max(xlim)],[20;20], 'color', 'r', 'LineWidth', 1);
title("Relative Poisson error in percent");
xlabel("velocity (m/s)");
ylabel("e_r (%)");







% differential number of counts and statistical errors with wide pixels not
% really relevent
% Cv_maxwell_norm = Cv_maxwell./max(Cv_maxwell);
% Cv_kappa_norm = Cv_kappa./max(Cv_kappa);
% dif_distrib = abs(Cv_maxwell_norm - Cv_kappa_norm);
% B0 = 0.1;
% distance = 2*v./(q*B0/m);
% step = distance(2) - distance(1);
% pixel_size = 3e-3;
% S_dif = 0;
% Count_pixel_dif = 0;
% distance_counts_dif = 0;
% k=1;
% for i=1:floor(pixel_size/step):length(distance)-floor(pixel_size/step)
%     for j=i:i+floor(pixel_size/step)
%         S_dif = S_dif + dif_distrib(j);
%     end
%     Count_pixel_dif(k) = S_dif;
%     S_dif = 0;
%     distance_counts_dif(k) = distance(i);
%     k = k + 1;
% end
% figure()
% hAx = subplot(211);
% hAx.YScale = 'log';
% e = errorbar(distance_counts_dif, Count_pixel_dif, sqrt(Count_pixel_dif)./2);
% e.CapSize = 10;
% xlim([0.06 0.15]);
% %ylim([1e-4 max(Count_pixel(1,:))*1.5]);
% title("number of particles as a function of velocity");
% xlabel("hit distance (m)");
% ylabel("number of counts");
% subplot(212);
% semilogy(distance_counts_dif, 1-(Count_pixel_dif - sqrt(Count_pixel_dif))./Count_pixel_dif);
% hold on;
% ylim([0 1]);
% xlim([0.06 0.15]);
% line([min(xlim);max(xlim)],[0.1;0.1], 'color', 'g', 'LineWidth', 1);
% line([min(xlim);max(xlim)],[0.2;0.2], 'color', 'r', 'LineWidth', 1);
% title("relative error in percent");
% xlabel("hit distance (m)");
% ylabel("e_r (%)");

% S = zeros(length(v)-1,1);
% Sj = 0;
% for j=1:length(v)-1
%     for i=1:j
%         Sj = Sj + (v(i+1)-v(i))*(y_maxwell(i+1)+y_maxwell(i))/2;
%     end
%     S(j) = Sj;
%     Sj = 0;
% end
% figure()
% plot(v(1:end-1), S);
% hold on;
% plot(v, Cv_maxwell);
% legend("integral count number", "count number per velocity");
% xlabel("velocity (m/s)");
% ylabel("number of counts");

