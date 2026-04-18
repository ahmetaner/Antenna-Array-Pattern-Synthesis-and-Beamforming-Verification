% PSO Hiperparametreleri İçin Lineer Arama (Grid Search)
clear; clc; close all;

disp('--- PSO İçin Lineer Arama (Grid Search) Başlıyor ---');
disp('Bu işlem bilgisayarınızın hızına bağlı olarak 15-20 dakika sürebilir...');

%% 1. Taranacak Katsayı Aralıkları
% Her katsayı için denenecek değerleri vektör olarak tanımlıyoruz.
% Örnek olarak her birinden 10'ar değer belirledik (10x10x10x10 = 10.000 deneme)

w_degerleri      = linspace(0.4, 1.2, 10);   % Atalet ağırlığı w (0.4 ile 1.2 arası 10 adım)
wdamp_degerleri  = linspace(0.90, 0.99, 10); % Atalet sönümleme w_damp (0.90 ile 0.99 arası 10 adım)
c1_degerleri     = linspace(1.0, 3.0, 10);   % Bilişsel öğrenme c1 (1.0 ile 3.0 arası 10 adım)
c2_degerleri     = linspace(1.0, 3.0, 10);   % Sosyal öğrenme c2 (1.0 ile 3.0 arası 10 adım)

toplam_deneme = length(w_degerleri) * length(wdamp_degerleri) * length(c1_degerleri) * length(c2_degerleri);
fprintf('Toplam Denenecek Kombinasyon Sayısı: %d\n\n', toplam_deneme);

%% 2. En İyi Sonuçları Tutacak Değişkenler
en_iyi_hata = inf;
en_iyi_katsayilar = zeros(1, 4); % [w, w_damp, c1, c2]

%% 3. İç İçe Döngülerle (Lineer Arama) Tüm Kombinasyonları Deneme
sayac = 0;
tic; % Süre ölçümünü başlat

for i = 1:length(w_degerleri)
    for j = 1:length(wdamp_degerleri)
        for k = 1:length(c1_degerleri)
            for m = 1:length(c2_degerleri)
                
                sayac = sayac + 1;
                
                % O anki denenecek katsayılar
                w_test      = w_degerleri(i);
                wdamp_test  = wdamp_degerleri(j);
                c1_test     = c1_degerleri(k);
                c2_test     = c2_degerleri(m);
                
                % UCA PSO'yu bu katsayılarla 1 kez çalıştır
                guncel_hata = ic_pso_calistir(w_test, wdamp_test, c1_test, c2_test);
                
                % Eğer bulunan hata şu ana kadarki en iyisiyse, kaydet
                if guncel_hata < en_iyi_hata
                    en_iyi_hata = guncel_hata;
                    en_iyi_katsayilar = [w_test, wdamp_test, c1_test, c2_test];
                end
                
                % Ekrana ilerleme durumunu yazdır (Her 500 adımda bir)
                if mod(sayac, 500) == 0
                    gecen_sure = toc;
                    kalan_sure_tahmini = (gecen_sure / sayac) * (toplam_deneme - sayac);
                    fprintf('İlerleme: %d / %d (%%%0.2f) | Şu anki En İyi Hata: %.2f | Tahmini Kalan Süre: %.1f sn\n', ...
                        sayac, toplam_deneme, (sayac/toplam_deneme)*100, en_iyi_hata, kalan_sure_tahmini);
                end
                
            end
        end
    end
end

gecen_sure_toplam = toc;

%% 4. Sonuçların Gösterimi
disp('======================================================');
fprintf('Lineer Arama (Grid Search) Tamamlandı! Toplam Süre: %.1f saniye\n', gecen_sure_toplam);
disp('BULUNAN EN İYİ (MÜKEMMEL) KATSAYILAR:');
fprintf('w      = %.3f\n', en_iyi_katsayilar(1));
fprintf('w_damp = %.3f\n', en_iyi_katsayilar(2));
fprintf('c1     = %.3f\n', en_iyi_katsayilar(3));
fprintf('c2     = %.3f\n', en_iyi_katsayilar(4));
fprintf('Bu katsayılarla ulaşılan en düşük hata: %.2f\n', en_iyi_hata);
disp('======================================================');


%% --- İÇ FONKSİYON: ASIL UCA PSO ALGORİTMASI ---
% Hızlı test için parametreleri optimize edilmiş asıl PSO kodu
function final_cost = ic_pso_calistir(w, w_damp, c1, c2)
    
    % Hedef Deseni (Orijinal Kod)
    phi_hedef = -180:10:180; 
    AF_istenen_dB = -20 * ones(size(phi_hedef));
    AF_istenen_dB(phi_hedef == 0)   = 0;  AF_istenen_dB(phi_hedef == 10)  = 0;
    AF_istenen_dB(phi_hedef == -10) = 0;  AF_istenen_dB(phi_hedef == 20)  = -5;
    AF_istenen_dB(phi_hedef == -20) = -5; AF_istenen_dB(phi_hedef == 30)  = -5;
    AF_istenen_dB(phi_hedef == -30) = -5; AF_istenen_dB(phi_hedef == 40)  = -10;
    AF_istenen_dB(phi_hedef == -40) = -10; AF_istenen_dB(phi_hedef == -120) = -40;

    N = 8; r_lambda = 0.5; bit_cuzunurluk = 4; q_adim = 360 / (2^bit_cuzunurluk); 
    
    % 10.000 kez çalışacağı için iç simülasyonu makul seviyede tutuyoruz
    num_particles = 40;   
    max_iter = 100;       
    
    num_vars = N; alt_sinir = 0; ust_sinir = 360; 
    
    pos = alt_sinir + rand(num_particles, num_vars) * (ust_sinir - alt_sinir);
    vel = zeros(num_particles, num_vars);

    phi_n_rad = deg2rad((0:N-1) * 360 / N);
    alfa_teorik = mod(-360 * r_lambda * cos(deg2rad(0) - phi_n_rad), 360);
    for i = 1:5
        pos(i, :) = alfa_teorik; 
    end

    pBest_pos = pos;
    pBest_cost = inf(num_particles, 1);
    gBest_cost = inf;
    gBest_pos = zeros(1, num_vars);

    phi_hedef_rad = deg2rad(phi_hedef);

    % PSO Döngüsü
    for iter = 1:max_iter
        for i = 1:num_particles
            % Fitness hesaplama (Orijinal Hata Fonksiyonu)
            fazlar_q = round(pos(i, :) / q_adim) * q_adim;
            AF = zeros(size(phi_hedef));
            for k = 1:N
                AF = AF + exp(1j * (2 * pi * r_lambda * cos(phi_hedef_rad - phi_n_rad(k)) + deg2rad(fazlar_q(k))));
            end
            AF_dB = 20*log10(abs(AF) + eps);
            AF_dB = AF_dB - max(AF_dB);
            
            hata_vektoru = zeros(size(phi_hedef));
            for k = 1:length(phi_hedef)
                aci = phi_hedef(k);
                if aci >= -20 && aci <= 20 
                    hata_vektoru(k) = abs(AF_istenen_dB(k) - AF_dB(k));
                else
                    if AF_dB(k) > AF_istenen_dB(k)
                        hata_vektoru(k) = AF_dB(k) - AF_istenen_dB(k);
                    end
                end
            end
            cost = sum(hata_vektoru.^2); 
            
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
        
        w = w * w_damp;  % Atalet Sönümleme Katsayısını Uygula
    end
    
    final_cost = gBest_cost;
end