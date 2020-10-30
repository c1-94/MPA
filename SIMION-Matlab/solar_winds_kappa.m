clc 
close all;


fast = csvread("flight_kappa_700e3_3e6_2e5.csv");
bulk = csvread("flight_kappa_500e3_5e6_1e5.csv");
slow = csvread("flight_kappa_400e3_5_5e6_8e4.csv");

y_center = 15.0;
radius = 1e-3;
n = 5e6;
v0 = 500e3; % bulk speed
mi = 1.672649e-27;  % ion mass
me = 9.11e-31;  % electron mass
ma = 6.64e-27;
m = mi; % christmas
q = 1.602176565E-19;    % elementary charge
kb = 1.38064852e-23;    % boltzman constant
T = 1e5; % 10⁵K
kappa = 3;
vt = (2*kb*T/m)^(1/2);  % thermal velocity
pixel_size = 1e-3;
instrument_length = 0.4;
B0 = 0.1;
delta_v = (q*B0*pixel_size)/(2*m);
A = pi*radius^2;
nb_pix = floor(instrument_length/pixel_size);
tau = 5e-3;
Phi0 = 360;  % azymuthal aperture angle
Theta0 = 5;    % elevation aperture angle
G = tau*A*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180));

final_posbulk = (bulk(:,2)-y_center).*1e-3;
final_posbulk = sort(final_posbulk);
Count_pixel_bulk = zeros(1,nb_pix);
v_counts_bulk = zeros(1,nb_pix);
for i=1:nb_pix
    upper = i*pixel_size;
    lower = (i-1)*pixel_size;
    Count_pixel_bulk(i) = sum(final_posbulk < upper) - sum(final_posbulk < lower);
    v_counts_bulk(i) = ((i*pixel_size*q*B0)/(2*m) + ((i-1)*pixel_size*q*B0)/(2*m))/2;
end

final_pos_fast = (fast(:,2)-y_center).*1e-3;
final_pos_fast = sort(final_pos_fast);
Count_pixel_fast = zeros(1,nb_pix);
v_counts_fast = zeros(1,nb_pix);
for i=1:nb_pix
    upper = i*pixel_size;
    lower = (i-1)*pixel_size;
    Count_pixel_fast(i) = sum(final_pos_fast < upper) - sum(final_pos_fast < lower);
    v_counts_fast(i) = ((i*pixel_size*q*B0)/(2*m) + ((i-1)*pixel_size*q*B0)/(2*m))/2;
end

final_pos_slow = (slow(:,2)-y_center).*1e-3;
final_pos_slow = sort(final_pos_slow);
Count_pixel_slow = zeros(1,nb_pix);
v_counts_slow = zeros(1,nb_pix);
for i=1:nb_pix
    upper = i*pixel_size;
    lower = (i-1)*pixel_size;
    Count_pixel_slow(i) = sum(final_pos_slow < upper) - sum(final_pos_slow < lower);
    v_counts_slow(i) = ((i*pixel_size*q*B0)/(2*m) + ((i-1)*pixel_size*q*B0)/(2*m))/2;
end

v_counts_bulk = v_counts_bulk(Count_pixel_bulk ~= 0);
Count_pixel_bulk = Count_pixel_bulk(Count_pixel_bulk ~= 0);
v_counts_fast = v_counts_fast(Count_pixel_fast ~= 0);
Count_pixel_fast = Count_pixel_fast(Count_pixel_fast ~= 0);
v_counts_slow = v_counts_slow(Count_pixel_slow ~= 0);
Count_pixel_slow = Count_pixel_slow(Count_pixel_slow ~= 0);


figure()
semilogy(1000*2*m*v_counts_bulk./(q*B0), Count_pixel_bulk, 'b-*', 'LineWidth', 2);
hold on;
set(gca,'FontSize',14);
semilogy(1000*2*m*v_counts_fast./(q*B0), Count_pixel_fast, 'g-*', 'LineWidth', 2);
semilogy(1000*2*m*v_counts_slow./(q*B0), Count_pixel_slow, 'r-*', 'LineWidth', 2);
%line([min(xlim);max(xlim)],[100;100], 'color', 'k', 'LineWidth', 3);

%title([{"Number of counts per pixel"},{"as a function of hit distance"}]);
xlabel("hit distance (mm)");
ylabel("number of counts");
legend("intermediate", "fast", "slow");


f_est_kappa_bulk = Count_pixel_bulk./(G*delta_v*v_counts_bulk.^3);
f_est_kappa_fast = Count_pixel_fast./(G*delta_v*v_counts_fast.^3);
f_est_kappa_slow = Count_pixel_slow./(G*delta_v*v_counts_slow.^3);
csvwrite('fitting_data_kappa_SIMION.txt', [v_counts_bulk; f_est_kappa_bulk]);

vt_fast = ((2*kb*2e5)/m)^(1/2);
vt_slow = ((2*kb*0.8e5)/m)^(1/2);
vt_bulk = ((2*kb*1e5)/m)^(1/2);

K_kappa = (2/(pi*(2*kappa-3)))^(3/2)*gamma(kappa+1)/gamma(kappa-1/2);

figure()
semilogy(v_counts_bulk, f_est_kappa_bulk, 'r-*', 'LineWidth', 2);
hold on;
set(gca,'FontSize',14);
semilogy(v_counts_bulk, (5e6/(vt_bulk^3))*K_kappa.*(1+(2/(2*kappa-3)).*((v_counts_bulk-500e3).^2/vt_bulk^2)).^(-kappa-1), 'b--', 'LineWidth', 2);
semilogy(v_counts_fast, f_est_kappa_fast, 'r-*', 'LineWidth', 2);
semilogy(v_counts_fast, (3e6/(vt_fast^3))*K_kappa.*(1+(2/(2*kappa-3)).*((v_counts_fast-700e3).^2/vt_fast^2)).^(-kappa-1), 'b--', 'LineWidth', 2);
semilogy(v_counts_slow, f_est_kappa_slow, 'r-*', 'LineWidth', 2);
semilogy(v_counts_slow, (5.5e6/(vt_slow^3))*K_kappa.*(1+(2/(2*kappa-3)).*((v_counts_slow-400e3).^2/vt_slow^2)).^(-kappa-1), 'b--', 'LineWidth', 2);

ylim([1e-11 1e-7])
legend("f_{estimated}", "f_{input}");
xlabel("velocity (m/s)");
ylabel("f (m^{-3} (m/s)^{-3})");
%title("Estimated and theoretical distributions");

csvwrite("fitting_data_kappa_SIMION_bulk.txt", [v_counts_bulk; f_est_kappa_bulk]);
csvwrite("fitting_data_kappa_SIMION_fast.txt", [v_counts_fast; f_est_kappa_fast]);
csvwrite("fitting_data_kappa_SIMION_slow.txt", [v_counts_slow; f_est_kappa_slow]);


