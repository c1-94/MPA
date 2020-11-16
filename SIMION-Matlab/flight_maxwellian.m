clc 
close all;

file = csvread("file_name_from_simion_output.csv");

y_center = 20.0;    % careful, can change with simion fly2 definition file

n = 5e6;    % density
v0 = 500e3; % bulk speed
mi = 1.672649e-27;  % ion mass
me = 9.11e-31;  % electron mass
ma = 6.64e-27;  % alpha particle mass
m = mi; % christmas
q = 1.602176565E-19;    % elementary charge
kb = 1.38064852e-23;    % boltzman constant
T = 1e5; % 10⁵K
vt = (2*kb*T/m)^(1/2);  % thermal velocity
radius = 1e-3;  % aperture radius
pixel_size = 1e-3;
instrument_length = 0.4;
B0 = 0.1;   % magnetic field strength
delta_v = (q*B0*pixel_size)/(2*m);
A = pi*radius^2;
nb_pix = floor(instrument_length/pixel_size);
tau = 5e-3;
Phi0 = 360;  % azymuthal aperture angle
Theta0 = 5;    % elevation aperture angle
G = tau*A*(Phi0/4*pi/180)*(1-cos(Theta0*pi/180));

final_pos = (file(:,2) - y_center).*1e-3;
final_pos = sort(final_pos);
Count_pixel = zeros(1,nb_pix);
v_counts = zeros(1,nb_pix);
for i=1:nb_pix
    upper = i*pixel_size;
    lower = (i-1)*pixel_size;
    Count_pixel(i) = sum(final_pos < upper) - sum(final_pos < lower);
    v_counts(i) = ((i*pixel_size*q*B0)/(2*m) + ((i-1)*pixel_size*q*B0)/(2*m))/2;
end

figure()
semilogy(1000*2*m*v_counts./(q*B0), Count_pixel, 'b-*', 'LineWidth', 2);
hold on;
line([min(xlim);max(xlim)],[100;100], 'color', 'g', 'LineWidth', 1);
title("Number of counts per pixel as a function of hit distance");
xlabel("hit distance (mm)");
ylabel("number of counts");


f_est_maxwell = Count_pixel./(G*delta_v*v_counts.^3);
csvwrite('fitting_data_maxwell_SIMION.txt', [v_counts; f_est_maxwell]);

figure()
semilogy(v_counts, f_est_maxwell, 'r-*', 'LineWidth', 2);
hold on;
semilogy(v_counts, (n/(pi^(3/2)*vt^3)).*exp(-(v_counts-v0).^2/vt^2), 'b--', 'LineWidth', 2);
ylim([1e-11 1e-7])
legend("f_{estimated}", "f");
xlabel("velocity (m/s)");
ylabel("f (m^{-3} (m/s)^{-3})");
title("Estimated and theoretical distributions");


