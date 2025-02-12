clc;
clear;
% Parameters for Ag (gold)
epsilon_inf = 1.0;  % High-frequency dielectric constant
sigma_epsilon0 = 1355.01;  % Conductivity over epsilon_0
planck_constant = 6.62607e-34;
light_speed=3e8;
e = 1.60217e-19;

A_p = [-8.577e4, -2.875, -997.6, -1.630];  % Ap [eV]
B_p = [-1.156e4, 0.0, -3090, -4.409];      % Bp [eV^2]
C_p = [5.557e7, 2.079e3, 6.921e5, 26.15];  % Cp [eV^2]
     

% Drude model parameters
omega_p = 5.8;  % Plasma frequency in eV for gold
gamma = 0.047;   % Damping factor (collision frequency) in eV

% Energy (eV) and Angular Frequency (omega)
wavelength= linspace(450,650,1000)*1e-9;
energy_joule = planck_constant*light_speed ./wavelength;  % Energy range from 0.5 to 10 eV
energy = energy_joule/e;
omega = energy;  % Assume h-bar omega ≈ energy in eV

% Initialize real and imaginary parts of the dielectric function
epsilon_real = epsilon_inf * ones(size(omega));
epsilon_imag = -sigma_epsilon0 ./ omega;  % Conductivity term

% Calculate the dielectric function using the Lorentz model
for p = 1:4
    denom = (omega.^2 + B_p(p)).^2 + (A_p(p) .* omega).^2;
    epsilon_real = epsilon_real + C_p(p) .* (omega.^2 + B_p(p)) ./ denom;
    epsilon_imag = epsilon_imag - C_p(p) .* (A_p(p) .* omega) ./ denom;
end

epsilon_lorentz = epsilon_real + epsilon_imag*1i ;

n_Au_lorentz = sqrt((abs(epsilon_lorentz) + epsilon_real) / 2) + ...
       1i * sqrt((abs(epsilon_lorentz) - epsilon_real) / 2);


% Drude model calculation
epsilon_real_drude = epsilon_inf - (omega_p^2 ./ (omega.^2 + gamma^2));
epsilon_imag_drude = (omega_p^2 .* gamma) ./ (omega .* (omega.^2 + gamma^2));


epsilon_drude = epsilon_imag_drude*1i + epsilon_real_drude ;
n_Au_drude = sqrt((abs(epsilon_drude) + epsilon_real_drude) / 2) + ...
       1i * sqrt((abs(epsilon_drude) - epsilon_real_drude) / 2);


% Radius in nm
radius = 30e-9; % nm

% Refractive index
refractive_index_medium = 1.5;

% Create arrays
radius_array = radius * ones(1, 1000);
refractive_index_medium_array = refractive_index_medium * ones(1, 1000);

function [Qext,Qsca,Qabs]=MieScattering2(lambda0,r0,n_m,n_Au)
%% parameters
m=n_Au/n_m;
k=2*pi*n_m/lambda0;
x=k*r0;
z=m*x;
N=round(2+x+4*x^(1/3));
%% computation
j=(1:N);
sqr=sqrt(pi*x/2);
sqrm=sqrt(pi*z/2);
phi=sqr.*besselj(j+0.5,x);
xi=sqr.*(besselj(j+0.5,x)+i*bessely(j+0.5,x));
phim=sqrm.*besselj(j+0.5,z);
phi1=[sin(x), phi(1:N-1)];
phi1m=[sin(z), phim(1:N-1)];
y=sqr*bessely(j+0.5,x);
y1=[-cos(x), y(1:N-1)];


phip=(phi1-j/x.*phi);
phimp=(phi1m-j/z.*phim);
xip=(phi1+i*y1)-j/x.*(phi+i*y);
aj=(m*phim.*phip-phi.*phimp)./(m*phim.*xip-xi.*phimp);
bj=(phim.*phip-m*phi.*phimp)./(phim.*xip-m*xi.*phimp);

Csca=sum( (2*j+1).*(abs(aj).*abs(aj)+abs(bj).*abs(bj)) );
Cext=sum( (2*j+1).*real(aj+bj) );


Qext=Cext*2*pi/(k*k);
Qsca=Csca*2*pi/(k*k);
Qabs=Qext-Qsca;
end

% lorentz loop

for i = 1:length(wavelength)
    n_m = 1.5;  
    % Call the MieScattering2 function
    [Qext, Qsca, Qabs] = MieScattering2(wavelength(i),radius, n_m, n_Au_lorentz(i));

    % Store the results
    Qext_array_lorentz(i) = Qext;
    Qsca_array_lorentz(i) = Qsca;
    Qabs_array_lorentz(i) = Qabs;
end

% drude loop

for i = 1:length(wavelength)
    n_m = 1.5;

    % Call the MieScattering2 function
    [Qext, Qsca, Qabs] = MieScattering2(wavelength(i), radius, n_m, n_Au_drude(i));

    % Store the results
    Qext_array_drude(i) = Qext;
    Qsca_array_drude(i) = Qsca;
    Qabs_array_drude(i) = Qabs;
end

% Plot the results for Lorentz and Drude models
figure;

% Lorentz Model Plot
subplot(2, 1, 1);
plot(wavelength * 1e9, Qext_array_lorentz, 'r-', 'LineWidth', 1.5); hold on;
plot(wavelength * 1e9, Qsca_array_lorentz, 'b-', 'LineWidth', 1.5);
plot(wavelength * 1e9, Qabs_array_lorentz, 'g-', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Cross Section');
legend('Extinction', 'Scattering', 'Absorption');
title('Lorentz Model');
grid on;

% Drude Model Plot
subplot(2, 1, 2);
plot(wavelength * 1e9, Qext_array_drude, 'r', 'LineWidth', 1.5); hold on;
plot(wavelength * 1e9, Qsca_array_drude, 'b', 'LineWidth', 1.5);
plot(wavelength * 1e9, Qabs_array_drude, 'g', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Cross Section');
legend('Extinction', 'Scattering', 'Absorption');
title('Drude Model');
grid on;

% Add a common title
sgtitle('Mie Scattering Cross Sections for Gold Nanoparticle');


    


