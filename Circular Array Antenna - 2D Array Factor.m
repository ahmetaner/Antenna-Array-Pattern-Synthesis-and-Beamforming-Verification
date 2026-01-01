clc;
clear;
mu_0 = 4*pi*1e-7;
eps_0 = 8.854187817e-12; % Birimi: F/m (Farad/metre)
w = 2*pi*2.4e9;
k = w * sqrt(mu_0*eps_0); %wave number
degtorad = pi/180; % dereceden radyana geçiş
N = 12; % antenna number
I_n = 1; % antenna phase
a = 2*pi*(1)/k; % radius of antennas
theta_focus = 90*degtorad; % elevation
phi_focus = 100*degtorad; % desired angle, azimuth
for n = 1:N
% her bir antenin x ekseni ile yaptığı açı
phi_antenna(n) = 2*pi*n/N; % indices start from 1 in the matlab
%her bir antenin çalışma fazı(istenilen açıda maksimum kazancı verecek şekilde)
phase_excitation(n) = -k*a*sin(theta_focus)*cos(phi_focus-phi_antenna(n));
end
theta = pi/2; % şu anki analizimiz iki boyutlu olduğu için değiştirmeye gerek yok
for phi = 1:360 %azimuth angle
phi_rad = phi*degtorad;
AF(phi) = 0; % phi açısı için başlangıç AF
for n = 1:N
% her bir döngüde tek bir antenden gelen katkı toplanır
AF(phi) = AF(phi) + I_n*exp(i*[k*a*sin(theta)*cos(phi_rad-phi_antenna(n)) + phase_excitation(n)]);
end
end
degree = 1:1:360;
% çok küçük array faktörü, negatif büyüklükte çok büyük desibeller
% oluşturuyor bu da grafiğin görünüşünü bozuyor. Bu durumdan kaçınmak için
% bu kod bloğu eklenmiştir.
for phi = 1:360
if 20*log10(abs(AF(phi))) < 0
AF(phi) = 1;
end
end
radian = degree * degtorad;
polarplot(radian, 20*log10(abs(AF)));
title("Hüzme Grafiği");
grid on;