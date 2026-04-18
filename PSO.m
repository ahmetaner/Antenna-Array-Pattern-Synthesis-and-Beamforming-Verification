% 8-Elemanlı 5-Bit UCA - Şablon Eşleştirme (Mask Matching) PSO
clear; clc; close all;

%% 1. Kullanıcı Girdisi: 10 Derecelik Adımlarla Hedef Desen (Maske)
% -180'den +180'e 10'ar derecelik adımlarla toplam 37 nokta
phi_hedef = -180:10:180; 

% Hedef güç değerleri (dB cinsinden). 
% Örnek bir desen oluşturuyorum: 
% 0 derecede 0 dB (Ana hüzme)
% +/- 90 derecede -40 dB (Null noktaları)
% Diğer yönlerde -15 dB (SLL)
AF_istenen_dB = -20 * ones(size(phi_hedef)); % Temel SLL seviyesi -15 dB

% Ana hüzme bölgesi (0 derece ve çevresi)
AF_istenen_dB(phi_hedef == 0)   = 0;
AF_istenen_dB(phi_hedef == 10)  = 0;
AF_istenen_dB(phi_hedef == -10) = 0;
AF_istenen_dB(phi_hedef == 20)  = -5;
AF_istenen_dB(phi_hedef == -20) = -5;
AF_istenen_dB(phi_hedef == 30)  = -5;
AF_istenen_dB(phi_hedef == -30) = -5;
AF_istenen_dB(phi_hedef == 40)  = -10;
AF_istenen_dB(phi_hedef == -40) = -10;


% Null (Sıfır) noktaları (+/- 90 derece)
%AF_istenen_dB(phi_hedef == 90)  = -40;
AF_istenen_dB(phi_hedef == -120) = -40;

% Not: Kendi hedef vektörünüzü doğrudan buraya manuel olarak girebilirsiniz.
% Örn: AF_istenen_dB = [-15, -15, -20, ... , 0, ... -40]; (37 elemanlı olmalı)

%% 2. Dizi ve Kuantalama Parametreleri
N = 8;                      
r_lambda = 0.5;             
bit_cuzunurluk = 4;         
q_adim = 360 / (2^bit_cuzunurluk); 

%% 3. PSO Ayarları
num_particles = 150;        
max_iter = 400;             
num_vars = N;               
alt_sinir = 0;              
ust_sinir = 360;            

w = 1; c1 = 2; c2 = 2; 

%% 4. PSO Başlatma ve Cophasal Aşı
pos = alt_sinir + rand(num_particles, num_vars) * (ust_sinir - alt_sinir);
vel = zeros(num_particles, num_vars);

% Ana hüzmenin 0 derecede olduğunu varsayarak analitik aşı (ilk 5 parçacık)
phi_n_rad = deg2rad((0:N-1) * 360 / N);
alfa_teorik = mod(-360 * r_lambda * cos(deg2rad(0) - phi_n_rad), 360);
for i = 1:5
    pos(i, :) = alfa_teorik; 
end

pBest_pos = pos;
pBest_cost = inf(num_particles, 1);
gBest_cost = inf;
gBest_pos = zeros(1, num_vars);

disp('Şablon Eşleştirmeli PSO Başlıyor...');

%% 5. Ana PSO Döngüsü
for iter = 1:max_iter
    for i = 1:num_particles
        cost = uca_mask_fitness(pos(i, :), phi_hedef, AF_istenen_dB, N, r_lambda, q_adim);
        
        if cost < pBest_cost(i)
            pBest_cost(i) = cost;
            pBest_pos(i, :) = pos(i, :);
        end
        if cost < gBest_cost
            gBest_cost = cost;
            gBest_pos = pos(i, :);
        end
    end
    
    for i = 1:num_particles
        r1 = rand(1, num_vars); r2 = rand(1, num_vars);
        vel(i, :) = w * vel(i, :) + c1 * r1 .* (pBest_pos(i, :) - pos(i, :)) + c2 * r2 .* (gBest_pos - pos(i, :));
        pos(i, :) = pos(i, :) + vel(i, :);
        pos(i, :) = max(pos(i, :), alt_sinir);
        pos(i, :) = min(pos(i, :), ust_sinir);
    end
    
    if mod(iter, 50) == 0
        fprintf('İterasyon %d/%d | En Düşük Hata (Cost): %.2f\n', iter, max_iter, gBest_cost);
    end
    
end

%% 6. Sonuçların Kuantalanması ve Çizimi
en_iyi_fazlar_q = round(gBest_pos / q_adim) * q_adim;
en_iyi_fazlar_q(en_iyi_fazlar_q == 360) = 0; 

fprintf('\nBulunan Kuantalanmış Faz Değerleri (Derece):\n');
disp(en_iyi_fazlar_q);

% Sürekli çizim için 1 derecelik hassasiyetle hesapla
phi_cizim = -180:1:180;
AF_cizim = zeros(size(phi_cizim));
for i = 1:N
    AF_cizim = AF_cizim + exp(1j * (2 * pi * r_lambda * cos(deg2rad(phi_cizim) - phi_n_rad(i)) + deg2rad(en_iyi_fazlar_q(i))));
end
AF_cizim_dB = 20*log10(abs(AF_cizim));
AF_cizim_dB = AF_cizim_dB - max(AF_cizim_dB); % Normalize

figure('Name', '8-Elemanlı UCA Mask Matching Sonucu');
plot(phi_cizim, AF_cizim_dB, 'b-', 'LineWidth', 2); hold on;
% İstenen hedef noktaları kırmızı yuvarlaklarla göster
scatter(phi_hedef, AF_istenen_dB, 50, 'r', 'filled'); 
ylim([-40 5]); xlim([-180 180]); grid on;
legend('Sentezlenen Desen', 'İstenen Hedef Noktalar (10° adımlı)', 'Location', 'best');
xlabel('Açı (Derece)'); ylabel('Normalize Dizi Faktörü (dB)');
title('Şablon Eşleştirme (Mask Matching) Performansı');


%% --- GÜNCELLENMİŞ FİTNESS FONKSİYONU ---
function cost = uca_mask_fitness(fazlar, phi_hedef, AF_istenen_dB, N, r_lambda, q_adim)
    % Fazları kuantize et
    fazlar_q = round(fazlar / q_adim) * q_adim;
    
    phi_hedef_rad = deg2rad(phi_hedef);
    phi_n_rad = deg2rad((0:N-1) * 360 / N);
    
    AF = zeros(size(phi_hedef));
    for i = 1:N
        AF = AF + exp(1j * (2 * pi * r_lambda * cos(phi_hedef_rad - phi_n_rad(i)) + deg2rad(fazlar_q(i))));
    end
    AF_dB = 20*log10(abs(AF) + eps);
    AF_dB = AF_dB - max(AF_dB); % Normalize et
    
    hata_vektoru = zeros(size(phi_hedef));
    
    for k = 1:length(phi_hedef)
        aci = phi_hedef(k);
        % Eğer açı Ana Hüzme bölgesindeyse (örneğin -20 ile 20 derece arası)
        % Şekli TAM olarak taklit etmeye çalışsın (Eşitlik Kısıtı)
        if aci >= -20 && aci <= 20 
            hata_vektoru(k) = abs(AF_istenen_dB(k) - AF_dB(k));
            
        % Eğer SLL veya Null bölgesindeyse
        % Sadece hedefi AŞARSA ceza ver (Eşitsizlik Kısıtı)
        else
            if AF_dB(k) > AF_istenen_dB(k)
                hata_vektoru(k) = AF_dB(k) - AF_istenen_dB(k);
            end
            % Hedefin altındaysa hata_vektoru(k) zaten 0 kalır.
        end
    end
    
    % Hataların karesini alarak toplam maliyeti bul
    cost = sum(hata_vektoru.^2); 
end