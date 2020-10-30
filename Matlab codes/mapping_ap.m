clc
close all;

B_field = (0.01:0.03:0.53)
pixel_length = [20 15 10 7.5 5 2 1.5 1 0.75 0.5 0.2 0.1 0.075 0.05 0.025].*1e-3;

kappa = 3;
N = 5e6;
v0 = 500e3;
T = 1e5; % 10⁵K
mi = 1.672649e-27;  % ion mass
me = 9.11e-31;  % electron mass
m = mi; % christmas
q = 1.602176565E-19;    % elementary charge
kb = 1.38064852e-23;    % boltzman constant
vt = (2*kb*T/m)^(1/2);  % thermal velocity


density_est_maxwell = csvread("density_maxwell.csv");
bulk_speed_est_maxwell = csvread("bulk_speed_maxwell.csv");
vt_est_maxwell = csvread("thermal_speed_maxwell.csv");
kappa_est = csvread("kappa_param.csv");
density_est_kappa = csvread("density_kappa.csv");
bulk_speed_est_kappa = csvread("bulk_speed_kappa.csv");
vt_est_kappa = csvread("thermal_speed_kappa.csv");


err_kappa = 100*abs(kappa - kappa_est)/kappa;
err_density_maxwell = 100*abs(N-density_est_maxwell)/N;
err_density_kappa = 100*abs(N-density_est_kappa)/N;
err_vt_maxwell = 100*abs(vt-vt_est_maxwell)/vt;
err_vt_kappa = 100*abs(vt-vt_est_kappa)/vt;
err_bulk_maxwell = 100*abs(v0-bulk_speed_est_maxwell)/v0;
err_bulk_kappa = 100*abs(v0-bulk_speed_est_kappa)/v0;


err_kappa(err_kappa > 100) = 100;
err_density_maxwell(err_density_maxwell > 100) = 100;
err_density_kappa(err_density_kappa > 100) = 100;
err_vt_maxwell(err_vt_maxwell > 100) = 100;
err_vt_kappa(err_vt_kappa > 100) = 100;
err_bulk_maxwell(err_bulk_maxwell > 100) = 100;
err_bulk_kappa(err_bulk_kappa > 100) = 100;

img_bulk_maxwell = interp2(err_bulk_maxwell, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_bulk_maxwell');
set(gca,'YScale','log','XScale','log')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title({"Relative error in % on the determination","of bulk speed for Maxwellian"});
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;

img_bulk_kappa = interp2(err_bulk_kappa, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_bulk_kappa');
set(gca,'YScale','log','XScale','log')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title({"Relative error in % on the determination", "of bulk speed for Kappa"});
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;


img_kappa = interp2(err_kappa, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_kappa');
set(gca,'YScale','log','XScale','log')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title("Relative error in % on the determination of \kappa");
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;

img_density_maxwell = interp2(err_density_maxwell, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_density_maxwell');
set(gca,'YScale','log','XScale','log')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title("Relative error in % on density for Maxwellian");
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;

img_density_kappa = interp2(err_density_kappa, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_density_kappa');
set(gca,'YScale','log','XScale','log')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title("Relative error in % on density for Kappa");
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;

img_vt_maxwell = interp2(err_vt_maxwell, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_vt_maxwell');
set(gca,'YScale','log','XScale','log')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title("Relative error in % on thermal speed for Maxwellian");
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;

img_vt_kappa = interp2(err_vt_kappa, 3);
figure()
imagesc(B_field', pixel_length'.*1e3, img_vt_kappa');
hold on
set(gca,'YScale','log', 'XScale','log', 'YMinorTick', 'on')
colorbar
axis([min(B_field') max(B_field') min(pixel_length'.*1e3) max(pixel_length'.*1e3)]);
% tt = title("Relative error in % on thermal speed for Kappa");
% tt.FontSize = 14;
xl = xlabel("Magnetic field strength (T)");
yl = ylabel("pixel size (mm)");
xl.FontSize = 14;
yl.FontSize = 14;

