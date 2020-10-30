clc;
close all;

N = 5e6;    % number density
v0 = 500e3; % bulk speed
mi = 1.672649e-27;  % ion mass
me = 9.11e-31;  % electron mass
m = mi; % christmas
q = 1.602176565E-19;    % elementary charge
kb = 1.38064852e-23;    % boltzman constant
T = 1e5; % 10⁵K
v_step = 1e3;
v = 1:v_step:2.5e6; % velocity range
vt = (2*kb*T/m)^(1/2);  % thermal velocity
B0 = 0.1;

instrument_length = 0.3;    %length of the instrument in meter
Phi0 = 360;  % azymuthal aperture angle
Theta0 = 5;    % elevation aperture angle
tau = 5e-3; % acquisition time
kappa = 3;


y_maxwell = v.^3.*exp(-(v-v0).^2/vt^2); % velocity part of the spherical 3D distribution
y_kappa = v.^3.*(1+(2/(2*kappa-3)).*((v-v0).^2/vt^2)).^(-kappa-1);
K_kappa = (2/(pi*(2*kappa-3)))^(3/2)*gamma(kappa+1)/gamma(kappa-1/2);

list_a = (0.1:0.2:3).*1e-3; % radius
list_l = [20 15 10 7.5 5 2 1.5 1 0.75 0.5 0.2 0.1 0.075 0.05 0.025].*1e-3;

for p_a=1:length(list_a)
    p_a
    for p_l=1:length(list_l)
        a = list_a(p_a);
        D = 2*a;
        A = pi*a^2;
        pixel_size = list_l(p_l);
        nb_pix = floor(instrument_length/pixel_size);
        G = tau*A*(Phi0*pi/180)*(1-cos(Theta0/2*pi/180));
        Cv_maxwell = N*G*(1/(pi^(3/2)*vt^3)).*y_maxwell;  % number of counts as a function of velocity
        Cv_kappa = N*G*(1/vt^3)*K_kappa.*y_kappa;
        
        distance_maxwell = [];
        distance_kappa = [];
        k_maxwell = 1;
        k_kappa = 1;
        for i=1:length(v)
            j = 0;
            count_lim_maxwell = Cv_maxwell(i)*v_step;
            while (j < floor(count_lim_maxwell))
                distance_maxwell(1,k_maxwell) = 2*v(i)./(q*B0/m) + D*rand() - D/2;
                k_maxwell = k_maxwell + 1;
                j = j + 1;
            end
            j = 0;
            count_lim_kappa = Cv_kappa(i)*v_step;
            while (j < floor(count_lim_kappa))
                distance_kappa(1,k_kappa) = 2*v(i)./(q*B0/m) + D*rand() - D/2;
                k_kappa = k_kappa + 1;
                j = j + 1;
            end
        end
        distance_maxwell = sort(distance_maxwell);
        distance_kappa = sort(distance_kappa);
        Count_pixel = zeros(2,nb_pix);    % maxwell | kappa
        v_counts = zeros(1,nb_pix);
        for i=1:nb_pix
            upper = (i)*pixel_size+D/2;
            lower = (i-1)*pixel_size+D/2;
            Count_pixel(1,i) = sum(distance_maxwell < upper) - sum(distance_maxwell < lower);
            Count_pixel(2,i) = sum(distance_kappa < upper) - sum(distance_kappa < lower);
            v_counts(i) = (((i-0.5)*pixel_size+D/2)*q*B0)/(2*m);
        end
        
        Counts_maxwell = Count_pixel(1,:);
        Counts_maxwell = Counts_maxwell(Counts_maxwell ~= 0);
        v_counts_maxwell = v_counts(Count_pixel(1,:) ~= 0);
        Counts_kappa = Count_pixel(2,:);
        Counts_kappa = Counts_kappa(Counts_kappa ~= 0);
        v_counts_kappa = v_counts(Count_pixel(2,:) ~= 0);
        
        delta_v = (q*B0*pixel_size)/(2*m);
        f_est_maxwell = Counts_maxwell./(G*delta_v*v_counts_maxwell.^3);
        csvwrite("fitting_data_maxwell_error_v_a_" + list_a(p_a) + "_l_" + list_l(p_l) + ".txt", [v_counts_maxwell; f_est_maxwell]);
        f_est_kappa = Counts_kappa./(G*delta_v*v_counts_kappa.^3);
        csvwrite("fitting_data_kappa_error_v_a_" + list_a(p_a) + "_l_" + list_l(p_l) + ".txt", [v_counts_kappa; f_est_kappa]);
    end
end

