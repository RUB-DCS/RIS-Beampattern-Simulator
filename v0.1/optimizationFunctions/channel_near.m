function [channel] = channel_near(dist,freq,c)


wavelength = c/freq;
%dist_size=size(dist);
%channel=zeros(dist_size(1),dist_size(2));


    %Ric = exp(-absorption_koeff*(dist(i,j)))/(1-exp(-absorption_koeff*(dist(i,j))));
    h_los= ( c./(4*pi*freq*dist) ) .* exp(1j*2*pi*(dist./wavelength));
    %h_absor= ( (1- exp(-absorption_koeff*dist(i,j)) )^0.5 ) * ( c/(4*pi*freq*dist(i,j)) ) * exp(1j*2*pi*rand(1,1));
    %channel(i,j)= sqrt(Ric/(Ric+1))*h_los + sqrt(1/(Ric+1))*h_absor;
    channel = h_los;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function absorption_koeff = initialize(freq,c)
% Polynomial approximation for absorption coefficient
% Physical constants
mu_perc = 90.0; %avg. annual humidity for Berlin
d = 10;                           % Tx-Rx distance in m [100-500]
mu = percentage2mu(mu_perc);        % VMR of water vapor
%f = linspace(100e9, 500e9, 4e2+1);    % Frequency grid 100-500GHz
f= linspace(200e9, 500e9, 3e2+1);
% Polynomial parameters

A = 5.159e-5 * (1-mu) * (-6.65e-5 * (1-mu) + 0.0159);
B = (-2.09e-4 * (1-mu) + 0.05) ^ 2;
C = 0.1925 * mu * (0.135 * mu + 0.0318);
D = (0.4241 * mu + 0.0998) ^ 2;
E = 0.2251 * mu * (0.1314 * mu + 0.0297);
F = (0.4127 * mu + 0.0932) ^ 2;
G = 2.053 * mu * (0.1717 * mu + 0.0306);
H = (0.5394 * mu + 0.0961) ^ 2;
I = 0.177 * mu * (0.0832 * mu + 0.0213);
J = (0.2615 * mu + 0.0668) ^ 2;
K = 2.146 * mu * (0.1206 * mu + 0.0277);
L = (0.3789 * mu + 0.0871) ^ 2;

p1 = 3.96e0;
p2 = 6.11e0;
p3 = 10.84e0;
p4 = 12.68e0;
p5 = 14.65e0;
p6 = 14.94e0;

a = 0.915e-112;
b = 9.42;


% Absorption polynomials

y1 = A ./ (B + (f / 100 / c - p1) .^ 2);
y2 = C ./ (D + (f / 100 / c - p2) .^ 2);
y3 = E ./ (F + (f / 100 / c - p3) .^ 2);
y4 = G ./ (H + (f / 100 / c - p4) .^ 2);
y5 = I ./ (J + (f / 100 / c - p5) .^ 2);
y6 = K ./ (L + (f / 100 / c - p6) .^ 2);

g = mu / 0.0157 * (2e-4 + a * f .^ b);


% Polynomial-based absorption + path-loss model

kappa = y1 + y2 + y3 + y4 + y5 + y6 + g;
PL = exp(kappa * d) .* (4 * pi * d * f) .^ 2 / c ^ 2;
Attenuation = exp(kappa * d);

s=find(f==freq);
absorption_koeff = kappa(s);
end

function mu = percentage2mu(mu_perc)
%            v = [0.0031 0.0094 0.0157 0.0220 0.0282]
%water vapor % = [10%     30%     50%    70%    90%]

%line-fitting
s = 3.137500000000000e-04;
t = -3.750000000000021e-05;

mu = s*mu_perc+t;
end
