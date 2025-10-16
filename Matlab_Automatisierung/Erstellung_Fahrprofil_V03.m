close all
clear

%% =========================================================================
%  Dateien laden
% =========================================================================
DATA_RAW = load("FSG24_Wet_Endu_Lap17.mat");   % Messdaten
Faktor_Abwaermen_Motoren = 0.5;

% AMK Loss-Map (liefert Speed [rpm], Shaft_Torque [Nm], Total_Loss [W])
load("4E4_W_T80_MT_80.mat");

% Maxi Motor Loss-Map aus CSV
filename = 'LossMap_Maxi_Motor_Final.csv';
Maxi_Motor_Loss_Map = readmatrix(filename);

% Sanity-Check: benötigte Felder in DATA_RAW
mustHave = {'Speed_combined_XY', ...
            's_AMK_fr_ActualVelocity','s_AMK_rr_ActualVelocity', ...
            's_AMK_fr_ActualTorque','s_AMK_rr_ActualTorque'};
for k = 1:numel(mustHave)
    assert(isfield(DATA_RAW,mustHave{k}), ...
        "Feld '%s' fehlt in East24_Dry_Endu_Lap20.mat", mustHave{k});
end

%% =========================================================================
%  Geschwindigkeit robust extrahieren (m/s & km/h) + gemeinsame Zeitbasis
% =========================================================================
sp = DATA_RAW.Speed_combined_XY;

if isfield(sp,'Value') && isfield(sp,'Time') && ~isempty(sp.Value)
    speed_mps  = sp.Value(:);
    speed_time = sp.Time(:);
elseif isfield(sp,'X') && isfield(sp,'Y') ...
       && isfield(sp.X,'Value') && isfield(sp.Y,'Value')
    vx         = sp.X.Value(:);
    vy         = sp.Y.Value(:);
    speed_mps  = hypot(vx,vy);
    if isfield(sp.X,'Time');      speed_time = sp.X.Time(:);
    elseif isfield(sp.Y,'Time');  speed_time = sp.Y.Time(:);
    else; error('Speed_combined_XY: keine Zeitachse gefunden.');
    end
else
    error('Speed_combined_XY hat eine unerwartete Struktur.');
end

% säubern & auf 0 beginnend legen
good          = isfinite(speed_mps) & isfinite(speed_time);
speed_mps     = speed_mps(good);
speed_time    = speed_time(good);
speed_time    = speed_time - speed_time(1);

% Heuristik für Einheit (m/s vs. km/h) – bei sehr großen Werten annehmen, es sind km/h:
if median(speed_mps,'omitnan') > 120/3.6    % > ~33 m/s
    speed_kmh = speed_mps;
    speed_mps = speed_kmh/3.6;
else
    speed_kmh = speed_mps*3.6;
end

% gemessene Abtastrate aus der Speed‑Zeitachse
if numel(speed_time) >= 2 && speed_time(end) > speed_time(1)
    Hz = (numel(speed_time)-1)/(speed_time(end)-speed_time(1));
else
    Hz = 100;  % Fallback
end

%% =========================================================================
%  Motor‑Zeitbasis (für Loss‑Map‑Interpolation) und Signale
% =========================================================================
t_mot  = DATA_RAW.s_AMK_fr_ActualVelocity.Time(:);
t_mot  = t_mot - t_mot(1);  % auf 0
velocity_fr = DATA_RAW.s_AMK_fr_ActualVelocity.Value(:).';   % rpm
velocity_rr = DATA_RAW.s_AMK_rr_ActualVelocity.Value(:).';   % rpm
moment_fr   = DATA_RAW.s_AMK_fr_ActualTorque.Value(:).';     % Nm
moment_rr   = DATA_RAW.s_AMK_rr_ActualTorque.Value(:).';     % Nm

% Längen angleichen (vorsichtshalber auf gemeinsame Minimallänge kürzen)
minN = min([numel(t_mot), numel(velocity_fr), numel(moment_fr), ...
            numel(velocity_rr), numel(moment_rr)]);
t_mot        = t_mot(1:minN);
velocity_fr  = velocity_fr(1:minN);
moment_fr    = moment_fr(1:minN);
velocity_rr  = velocity_rr(1:minN);
moment_rr    = moment_rr(1:minN);

%% =========================================================================
%  Interpolation Loss‑Map — Vorderachse (AMK)
% =========================================================================
X = Speed;             % rpm
Y = Shaft_Torque;      % Nm
Z = Total_Loss;        % W

[Xi, Yi] = meshgrid(linspace(min(X(:)), max(X(:)), 100), ...
                    linspace(min(Y(:)), max(Y(:)), 100));
Zi = griddata(X, Y, Z, Xi, Yi, 'cubic');

x_val = velocity_fr;
y_val = abs(moment_fr);
Interpolated_Loss_FR = griddata(Xi, Yi, Zi, x_val, y_val, 'cubic');

%% =========================================================================
%  Interpolation Loss‑Map — Vorderachse (Maxi Motor 94.1)
% =========================================================================
X_M = Maxi_Motor_Loss_Map(1, 2:end);   % rpm
Y_M = Maxi_Motor_Loss_Map(2:end, 1);   % Nm
Z_M = Maxi_Motor_Loss_Map(2:end, 2:end);

[Xm, Ym] = meshgrid(linspace(min(X_M(:)), max(X_M(:)), 50), ...
                    linspace(min(Y_M(:)), max(Y_M(:)), 50));
Zm = griddata(X_M, Y_M, Z_M, Xm, Ym, 'cubic');

Interpolated_Loss_FR_M = griddata(Xm, Ym, Zm, x_val, y_val, 'cubic');

%% =========================================================================
%  Interpolation Loss‑Map — Hinterachse (AMK & Maxi)
% =========================================================================
x_val = velocity_rr;
y_val = abs(moment_rr);

Interpolated_Loss_RR   = griddata(Xi, Yi, Zi, x_val, y_val, 'cubic');

[Xm2, Ym2] = meshgrid(linspace(min(X_M(:)), max(X_M(:)), 100), ...
                      linspace(min(Y_M(:)), max(Y_M(:)), 100));
Zm2 = griddata(X_M, Y_M, Z_M, Xm2, Ym2, 'cubic');
Interpolated_Loss_RR_M = griddata(Xm2, Ym2, Zm2, x_val, y_val, 'cubic');

% NaNs -> 0
Interpolated_Loss_FR( isnan(Interpolated_Loss_FR) ) = 0;
Interpolated_Loss_FR_M(isnan(Interpolated_Loss_FR_M))= 0;
Interpolated_Loss_RR( isnan(Interpolated_Loss_RR) ) = 0;
Interpolated_Loss_RR_M(isnan(Interpolated_Loss_RR_M))= 0;

%% =========================================================================
%  Plots Verluste (Motor‑Zeitbasis)
% =========================================================================
figure(2); tiledlayout(1,2)
nexttile
plot(t_mot, Interpolated_Loss_FR/2); hold on
ybar = mean(Interpolated_Loss_FR)*Faktor_Abwaermen_Motoren;
plot([t_mot(1), t_mot(end)], [ybar, ybar], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung | AMK VL'); grid on
legend('Data', sprintf('Average (%.1f W)', ybar));

nexttile
plot(t_mot, Interpolated_Loss_FR_M); hold on
ybar = mean(Interpolated_Loss_FR_M);
plot([t_mot(1), t_mot(end)], [ybar, ybar], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung | Motor94.1 VL'); grid on
legend('Data', sprintf('Average (%.1f W)', ybar));

figure(3); tiledlayout(1,2)
nexttile
plot(t_mot, Interpolated_Loss_RR/2); hold on
ybar = mean(Interpolated_Loss_RR)*Faktor_Abwaermen_Motoren;
plot([t_mot(1), t_mot(end)], [ybar, ybar], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung | AMK HL'); grid on
legend('Data', sprintf('Average (%.1f W)', ybar));

nexttile
plot(t_mot, Interpolated_Loss_RR_M); hold on
ybar = mean(Interpolated_Loss_RR_M);
plot([t_mot(1), t_mot(end)], [ybar, ybar], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung | Motor94.1 HL'); grid on
legend('Data', sprintf('Average (%.1f W)', ybar));

%% =========================================================================
%  Abwärmeleistung Traktionsinverter (Motor‑Zeitbasis)
% =========================================================================
DATA_Q_TI(1,:) = t_mot.';                                  % Zeit (s)
DATA_Q_TI(2,:) = velocity_fr(1:minN);                      % rpm VA
DATA_Q_TI(3,:) = moment_fr(1:minN);                        % Nm  VA
DATA_Q_TI(4,:) = velocity_rr(1:minN);                      % rpm HA
DATA_Q_TI(5,:) = moment_rr(1:minN);                        % Nm  HA

% FSB-Anteile
Anteil_FSB_VL = sum(DATA_Q_TI(2,:) > 11000)/numel(DATA_Q_TI(2,:));
Anteil_FSB_HL = sum(DATA_Q_TI(4,:) > 11000)/numel(DATA_Q_TI(4,:));

% Durchschnittsmomente (robust)
AVG_Moment_FSB_VL        = sum(abs(DATA_Q_TI(3, DATA_Q_TI(2, :) > 11000))) / max(Anteil_FSB_VL     * numel(DATA_Q_TI(2,:)), eps);
AVG_Moment_nicht_FSB_VL  = sum(abs(DATA_Q_TI(3, DATA_Q_TI(2, :) < 11000))) / max((1-Anteil_FSB_VL) * numel(DATA_Q_TI(2,:)), eps);
AVG_Moment_FSB_HL        = sum(abs(DATA_Q_TI(5, DATA_Q_TI(4, :) > 11000))) / max(Anteil_FSB_HL     * numel(DATA_Q_TI(4,:)), eps);
AVG_Moment_nicht_FSB_HL  = sum(abs(DATA_Q_TI(5, DATA_Q_TI(4, :) < 11000))) / max((1-Anteil_FSB_HL) * numel(DATA_Q_TI(4,:)), eps);

% Ströme (AMK: 1/0.38 ARMS/Nm)
Motorkonstante_AMK = 1/0.38; % ARMS/Nm
DATA_Q_TI(6,:) = abs(DATA_Q_TI(3,:))*Motorkonstante_AMK; % Strom VA
DATA_Q_TI(7,:) = abs(DATA_Q_TI(5,:))*Motorkonstante_AMK; % Strom HA

% Inverter-Verluste
% Infineon FS200R12PT4 (6 Schalter)
DATA_Q_TI(8,:)  = 6*(0.002589.*DATA_Q_TI(6,:).*DATA_Q_TI(6,:) + 0.557347.*DATA_Q_TI(6,:) + 0.945667); 
DATA_Q_TI(9,:)  = 6*(0.002589.*DATA_Q_TI(7,:).*DATA_Q_TI(7,:) + 0.557347.*DATA_Q_TI(7,:) + 0.945667);
% Powermodul Dresden 10-EY126PB011ME-PJ19F18T
DATA_Q_TI(10,:) = (0.0481.*DATA_Q_TI(6,:).*DATA_Q_TI(6,:) - 0.2247.*DATA_Q_TI(6,:) + 8.9979);
DATA_Q_TI(11,:) = (0.0481.*DATA_Q_TI(7,:).*DATA_Q_TI(7,:) - 0.2247.*DATA_Q_TI(7,:) + 8.9979);

% Plots TI
figure(4); tiledlayout(1,2)
nexttile
scatter(DATA_Q_TI(2,:),DATA_Q_TI(3,:), 2); hold on
xline(11000, 'r--'); yline(0, 'k--')
xlabel('Drehzahl in 1/min'); ylabel('Drehmoment in Nm');
xlim([0, 18000]); ylim([-15 30]); grid on; set(gca,'FontSize',14)
title('VA Drehzahl–Moment');

nexttile
scatter(DATA_Q_TI(4,:),DATA_Q_TI(5,:), 2); hold on
xline(11000, 'r--'); yline(0, 'k--')
xlabel('Drehzahl in 1/min'); ylabel('Drehmoment in Nm');
xlim([0, 18000]); ylim([-15 30]); grid on; set(gca,'FontSize',14)
title('HA Drehzahl–Moment');

figure(5); tiledlayout(1,2)
nexttile
plot(DATA_Q_TI(1,:)-DATA_Q_TI(1,1), DATA_Q_TI(8,:)); hold on
plot([t_mot(1), t_mot(end)], [mean(DATA_Q_TI(8,:)), mean(DATA_Q_TI(8,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('Infineon FS200R12PT4 VL'); grid on; set(gca,'FontSize',14)

nexttile
plot(DATA_Q_TI(1,:)-DATA_Q_TI(1,1), DATA_Q_TI(9,:)); hold on
plot([t_mot(1), t_mot(end)], [mean(DATA_Q_TI(9,:)), mean(DATA_Q_TI(9,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('Infineon FS200R12PT4 HL'); grid on; set(gca,'FontSize',14)

figure(6); tiledlayout(1,2)
nexttile
plot(DATA_Q_TI(1,:)-DATA_Q_TI(1,1), DATA_Q_TI(10,:)); hold on
plot([t_mot(1), t_mot(end)], [mean(DATA_Q_TI(10,:)), mean(DATA_Q_TI(10,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('10-EY126PB011ME-PJ19F18T VL'); grid on; set(gca,'FontSize',14)

nexttile
plot(DATA_Q_TI(1,:)-DATA_Q_TI(1,1), DATA_Q_TI(11,:)); hold on
plot([t_mot(1), t_mot(end)], [mean(DATA_Q_TI(11,:)), mean(DATA_Q_TI(11,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('10-EY126PB011ME-PJ19F18T HL'); grid on; set(gca,'FontSize',14)

%% =========================================================================
%  Fahrprofil (auf die Speed‑Zeitachse interpoliert)
% =========================================================================
% Interpolation der Motor‑verlust-Leistungen auf speed_time
interpLoss = @(t,y) interp1(t, y(:), speed_time, 'linear', 'extrap');

Q_FR_AMK   = interpLoss(t_mot, Interpolated_Loss_FR);
Q_RR_AMK   = interpLoss(t_mot, Interpolated_Loss_RR);
Q_FR_MAXI  = interpLoss(t_mot, Interpolated_Loss_FR_M);
Q_RR_MAXI  = interpLoss(t_mot, Interpolated_Loss_RR_M);
Q_TI_VA_Inf= interpLoss(t_mot, DATA_Q_TI(8,:));
Q_TI_HA_Inf= interpLoss(t_mot, DATA_Q_TI(9,:));
Q_TI_VA_SiC= interpLoss(t_mot, DATA_Q_TI(10,:));
Q_TI_HA_SiC= interpLoss(t_mot, DATA_Q_TI(11,:));

Fahrprofil = zeros(10, numel(speed_time));
Fahrprofil(1,:)  = speed_time.';   % Zeit
Fahrprofil(2,:)  = speed_mps.';    % v in m/s
Fahrprofil(3,:)  = Q_FR_AMK.';     % Verluste VL (AMK)
Fahrprofil(4,:)  = Q_RR_AMK.';     % Verluste HL (AMK)
Fahrprofil(5,:)  = Q_TI_VA_Inf.';  % Inverter VA (Infineon)
Fahrprofil(6,:)  = Q_TI_HA_Inf.';  % Inverter HA (Infineon)
Fahrprofil(7,:)  = Q_TI_VA_SiC.';  % Inverter VA (SiC)
Fahrprofil(8,:)  = Q_TI_HA_SiC.';  % Inverter HA (SiC)
Fahrprofil(9,:)  = Q_FR_MAXI.';    % Verluste VL Motor94.1
Fahrprofil(10,:) = Q_RR_MAXI.';    % Verluste HL Motor94.1

% Laps duplizieren (optional)
No_Laps   = 10;
Cols_BSL  = size(Fahrprofil,2);
Laptime   = Fahrprofil(1,Cols_BSL) - Fahrprofil(1,1);

if No_Laps > 1
    for k = 2:No_Laps
        c0 = size(Fahrprofil,2);
        Fahrprofil(2:10, c0+(1:Cols_BSL)) = Fahrprofil(2:10, 1:Cols_BSL);
    end
    % neue Zeit auf speed_time basierend
    Fahrprofil(1,:) = 0:(1/Hz):(No_Laps*Laptime + (1/Hz)*(No_Laps-1));
end
Endtime = No_Laps * Laptime %#ok<NOPTS>

%% =========================================================================
%  Downsampling auf Zielrate (z. B. 20 Hz)
% =========================================================================
target_datarate   = 20;
decimation_factor = max(1, round(Hz / target_datarate));
Fahrprofil_reduziert = Fahrprofil(1:10, 1:decimation_factor:size(Fahrprofil,2));

%% =========================================================================
%  Pumpe (Beispielkennlinie)
% =========================================================================
Pumpenkennlinie_Davis_Craig(1,:) = 0:5:30;                         % l/min
Pumpenkennlinie_Davis_Craig(1,:) = Pumpenkennlinie_Davis_Craig(1,:)/60;    % l/s
Pumpenkennlinie_Davis_Craig(1,:) = Pumpenkennlinie_Davis_Craig(1,:)/1000;  % m^3/s
Pumpenkennlinie_Davis_Craig(2,:) = [0.84321 0.81771 0.76721 0.69171 0.59121 0.46571 0.31521]; % bar

%% =========================================================================
%  Geschwindigkeitsplots (Figure 7/8)
% =========================================================================
% Histogramm
speeds_kmh = speed_kmh(:).';
classEdges = 0:5:100;
[counts, edges] = histcounts(speeds_kmh, classEdges);
relativeFrequencies = (counts / max(1,sum(counts))) * 100;

figure;
bar(edges(1:end-1), relativeFrequencies, 'FaceColor', "#0072BD");
xlabel('Geschwindigkeit in km/h', 'FontSize', 14);
ylabel('Relative Häufigkeit in %', 'FontSize', 14);
title('Geschwindigkeitsverteilung');
set(gca, 'FontSize', 14); grid on
xticks(0:10:100);

% Zeitverlauf (echte Zeit auf x)
figure
plot(speed_time, speeds_kmh, 'LineWidth', 1.1)
xlabel("Zeit in s")
ylabel("Fahrgeschwindigkeit in km/h")
set(gca, 'FontSize', 12);
grid on
axis tight
