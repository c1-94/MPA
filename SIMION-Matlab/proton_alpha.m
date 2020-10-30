clc
%close all;

file = csvread("proton_alpha.csv");

y_center = 15.0;    % careful, can change with simion fly2 definition file

n_a = (3.5/100)*5e6;    % alpha density
n = 5e6-n_a;    % proton density
v0 = 500e3; % bulk speed
mi = 1.672649e-27;  % ion mass
me = 9.11e-31;  % electron mass
ma = 6.64e-27;  % alpha particle mass
m = mi; % christmas
q = 1.602176565E-19;    % elementary charge
kb = 1.38064852e-23;    % boltzman constant
T = 1e5; % 10⁵K
T_a = 1.5e5;
vt = (2*kb*T/m)^(1/2);  % thermal velocity
vt_a = (2*kb*T_a/(4*m))^(1/2);
radius = 1e-3;  % aperture radius
pixel_size = 1e-3;
instrument_length = 0.4;
B0 = 0.1;
delta_v = (q*B0*pixel_size)/(2*m);
delta_v_a = (2*q*B0*pixel_size)/(2*ma);
A = pi*radius^2;
nb_pix = floor(instrument_length/pixel_size);
tau = 5e-3;
Phi0 = 360;  % azymuthal aperture angle
Theta0 = 5;    % elevation aperture angle
G = tau*A*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180));

final_pos = (file(:,2)-y_center).*1e-3;
final_pos = sort(final_pos);
Count_pixel = zeros(1,nb_pix);
v_counts = zeros(1,nb_pix);
j = 1;
k = 1;
count_alpha = [];
v_counts_alpha = [];
count_proton = [];
v_counts_proton = [];
for i=1:nb_pix
    upper = i*pixel_size;
    lower = (i-1)*pixel_size;
    Count_pixel(i) = sum(final_pos < upper) - sum(final_pos < lower);
    v_counts(i) = ((i*pixel_size*q*B0)/(2*m) + ((i-1)*pixel_size*q*B0)/(2*m))/2;
    if v_counts(i) > 1.5*v0
        count_alpha(j) = Count_pixel(i);
        v_counts_alpha(j) = v_counts(i);
        j = j + 1;
    else
        count_proton(k) = Count_pixel(i);
        v_counts_proton(k) = v_counts(i);
        k = k + 1;
    end
    
end


figure()
semilogy(1000*2*m*v_counts_proton./(q*B0), count_proton, 'b-*', 'LineWidth', 2);
hold on;
set(gca,'FontSize',14);
semilogy(1000*2*m*v_counts_alpha./(q*B0), count_alpha, 'g-*', 'LineWidth', 2);
%title([{"Number of counts per pixel"}, {"as a function of hit distance"}]);
xlabel("hit distance (mm)");
ylabel("number of counts");
legend("protons", "alphas");

v_counts_alpha = v_counts_alpha./2;

f_est_proton = count_proton./(G*delta_v*v_counts_proton.^3);
csvwrite('fitting_data_maxwell_SIMION_proton.txt', [v_counts_proton; f_est_proton]);
f_est_alpha = count_alpha./(G*delta_v_a*v_counts_alpha.^3);
csvwrite('fitting_data_maxwell_SIMION_alpha.txt', [v_counts_alpha; f_est_alpha]);

figure()
semilogy(v_counts_proton, f_est_proton, 'r-*', 'LineWidth', 2);
hold on;
set(gca,'FontSize',14);
semilogy(v_counts_alpha, f_est_alpha, 'g-o', 'LineWidth', 2);
semilogy(v_counts_proton, (n/(pi^(3/2)*vt^3)).*exp(-(v_counts_proton-v0).^2/vt^2), 'b--', 'LineWidth', 2);
semilogy(v_counts_alpha, (n_a/(pi^(3/2)*vt_a^3)).*exp(-(v_counts_alpha-v0).^2/vt_a^2), 'k--', 'LineWidth', 2);
ylim([1e-11 1e-7])
legend("f_{estimated} proton", "f_{estimated} alpha", "f_{input} proton", "f_{input} alpha");
xlabel("velocity (m/s)");
ylabel("f (m^{-3} (m/s)^{-3}");
%title("    Estimated and theoretical distributions");



