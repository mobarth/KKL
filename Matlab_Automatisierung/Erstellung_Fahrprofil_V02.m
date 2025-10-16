close all
clear

%% ------------------------------------------------------------------------
%  Daten laden
% -------------------------------------------------------------------------
DATA_RAW = load("East24_Dry_Endu_Lap20.mat");
load("FS_TMM_Workspace.mat")                  % falls benötigt (Speed, Shaft_Torque, Total_Loss)
Faktor_Abwaermen_Motoren = 0.5;

% Ursprüngliche Loss-Map (AMK) – enthält Speed, Shaft_Torque, Total_Loss
load("4E4_W_T80_MT_80.mat");

% Fahrdatensatz-Name (für eval-Nutzung beibehalten)
name        = 'AutoX_FSG_CM';      % rein interner Name
repetitions = 1;
File        = '.mat';
load ([name, File]);               % lädt Variablen wie s_AMK_fr_ActualVelocity etc. aus dem Workspace-File

% --- Sanity-Check: erwartete Felder
mustHave = {'Speed_combined_XY', ...
            's_AMK_fr_ActualVelocity','s_AMK_rr_ActualVelocity', ...
            's_AMK_fr_ActualTorque','s_AMK_rr_ActualTorque'};
for k = 1:numel(mustHave)
    assert(isfield(DATA_RAW,mustHave{k}), ...
        "Feld '%s' fehlt in East24_Dry_Endu_Lap20.mat", mustHave{k});
end

%% ------------------------------------------------------------------------
%  Vorbereitung der Zeitbasis
% -------------------------------------------------------------------------
% Hz aus dem Front‑Drehzahl‑Signal bestimmen
ersterWert  = s_AMK_fr_ActualVelocity.Time(1);
letzterWert = s_AMK_fr_ActualVelocity.Time(end);
anzahlWerte = length(s_AMK_fr_ActualVelocity.Time);
t_ges       = letzterWert - ersterWert;
Hz          = (anzahlWerte-1)/t_ges;

% einheitlichen Zeitvektor erzeugen
eval([['time_', name], ' = linspace(0, t_ges*repetitions, (1+t_ges*repetitions * Hz));']);


%% ------------------------------------------------------------------------
%  Fahrdaten: Drehzahl / Moment VA/HA für Wiederholungen
% -------------------------------------------------------------------------
eval([['velocity_fr_', name], ' = repmat(s_AMK_fr_ActualVelocity.Value, 1, repetitions);']);
eval([['velocity_rr_', name], ' = repmat(s_AMK_rr_ActualVelocity.Value, 1, repetitions);']);

% !!! Anpassung auf ActualTorque
eval([['moment_fr_', name],   ' = repmat(s_AMK_fr_ActualTorque.Value,   1, repetitions);']);
eval([['moment_rr_', name],   ' = repmat(s_AMK_rr_ActualTorque.Value,   1, repetitions);']);

% CSV-Datei Maxi Motor Loss Map importieren
filename = 'LossMap_Maxi_Motor_Final.csv';
Maxi_Motor_Loss_Map = readmatrix(filename);

%% ------------------------------------------------------------------------
%  Interpolation Loss‑Map — Vorderachse (AMK)
% -------------------------------------------------------------------------
Motor    = 'fr';
time     = eval(['time_', name]);
velocity = eval(['velocity_' ,Motor,'_', name]);      % rpm
moment   = eval(['moment_'   ,Motor,'_', name]);      % Nm

X = Speed;                   % rpm
Y = Shaft_Torque;            % Nm
Z = Total_Loss;              % W

[X_interp, Y_interp] = meshgrid(linspace(min(X(:)), max(X(:)), 100), ...
                                linspace(min(Y(:)), max(Y(:)), 100));
Z_interp = griddata(X, Y, Z, X_interp, Y_interp, 'cubic');

x_value = velocity;
y_value = abs(moment);
Interpolated_Loss_FR = griddata(X_interp, Y_interp, Z_interp, x_value, y_value, 'cubic');

%% ------------------------------------------------------------------------
%  Interpolation Loss‑Map — Vorderachse (Motor94.1 Maxi Map)
% -------------------------------------------------------------------------
X_M = Maxi_Motor_Loss_Map(1, 2:end);   % rpm
Y_M = Maxi_Motor_Loss_Map(2:end, 1);   % Nm
Z_M = Maxi_Motor_Loss_Map(2:end, 2:end);

[X_M_interp, Y_M_interp] = meshgrid(linspace(min(X_M(:)), max(X_M(:)), 50), ...
                                    linspace(min(Y_M(:)), max(Y_M(:)), 50));
Z_M_interp = griddata(X_M, Y_M, Z_M, X_M_interp, Y_M_interp, 'cubic');

% Visualisierung (optional)
figure
surf(X_M_interp, Y_M_interp, Z_M_interp)
ylabel("Motordrehmoment in Nm")
xlabel("Motordrehzahl in 1/min")
zlabel("Verlustabwärme in W")
set(gca, 'FontSize', 20);
xticks([0, 5000, 10000, 15000, 20000, 25000]);
yticks([0, 5, 10, 15, 20, 25, 30, 35]);
zticks([0, 1000, 2000, 3000, 4000, 5000]);
colormap("turbo")

Interpolated_Loss_FR_M = griddata(X_M_interp, Y_M_interp, Z_M_interp, x_value, y_value, 'cubic');

%% ------------------------------------------------------------------------
%  Verlustleistung plotten — Front
% -------------------------------------------------------------------------
Interpolated_Loss_FR(isnan(Interpolated_Loss_FR)) = 0;
Interpolated_Loss_FR_M(isnan(Interpolated_Loss_FR_M)) = 0;

figure(2); tiledlayout(1,2)
nexttile
plot(time, Interpolated_Loss_FR/2); hold on
plot([min(time), max(time)], [mean(Interpolated_Loss_FR)*Faktor_Abwaermen_Motoren, ...
     mean(Interpolated_Loss_FR)*Faktor_Abwaermen_Motoren], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung während FSG AutoX Runde | AMK VL'); grid on
legend('Data', ['Average losses (' num2str(mean(Interpolated_Loss_FR)*Faktor_Abwaermen_Motoren) ')']);

nexttile
plot(time, Interpolated_Loss_FR_M); hold on
plot([min(time), max(time)], [mean(Interpolated_Loss_FR_M), mean(Interpolated_Loss_FR_M)], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung während FSG AutoX Runde | Motor94.1 VL'); grid on
legend('Data', ['Average losses (' num2str(mean(Interpolated_Loss_FR_M)) ')']);

%% ------------------------------------------------------------------------
%  Interpolation Loss‑Map — Hinterachse (AMK & Maxi)
% -------------------------------------------------------------------------
Motor    = 'rr';
velocity = eval(['velocity_' ,Motor,'_', name]);
moment   = eval(['moment_'   ,Motor,'_', name]);

[X_interp, Y_interp] = meshgrid(linspace(min(X(:)), max(X(:)), 100), ...
                                linspace(min(Y(:)), max(Y(:)), 100));
Z_interp = griddata(X, Y, Z, X_interp, Y_interp, 'cubic');

x_value = velocity; y_value = abs(moment);
Interpolated_Loss_RR = griddata(X_interp, Y_interp, Z_interp, x_value, y_value, 'cubic');

[X_M_interp, Y_M_interp] = meshgrid(linspace(min(X_M(:)), max(X_M(:)), 100), ...
                                    linspace(min(Y_M(:)), max(Y_M(:)), 100));
Z_M_interp = griddata(X_M, Y_M, Z_M, X_M_interp, Y_M_interp, 'cubic');

Interpolated_Loss_RR_M = griddata(X_M_interp, Y_M_interp, Z_M_interp, x_value, y_value, 'cubic');

%% ------------------------------------------------------------------------
%  Verlustleistung plotten — Heck
% -------------------------------------------------------------------------
Interpolated_Loss_RR(isnan(Interpolated_Loss_RR))   = 0;
Interpolated_Loss_RR_M(isnan(Interpolated_Loss_RR_M)) = 0;

figure(3); tiledlayout(1,2)
nexttile
plot(time, Interpolated_Loss_RR/2); hold on
plot([min(time), max(time)], [mean(Interpolated_Loss_RR)*Faktor_Abwaermen_Motoren, ...
     mean(Interpolated_Loss_RR)*Faktor_Abwaermen_Motoren], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung während FSG AutoX Runde | AMK HL'); grid on
legend('Daten', ['Average losses (' num2str(mean(Interpolated_Loss_RR)*Faktor_Abwaermen_Motoren) ')']);

nexttile
plot(time, Interpolated_Loss_RR_M); hold on
plot([min(time), max(time)], [mean(Interpolated_Loss_RR_M), mean(Interpolated_Loss_RR_M)], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 2000]);
title('Abwärmeleistung während FSG AutoX Runde | Motor94.1 HL'); grid on
legend('Data', ['Average losses (' num2str(mean(Interpolated_Loss_RR_M)) ')']);

%% ------------------------------------------------------------------------
%  Berechnung Abwärmeleistung Traktionsinverter
% -------------------------------------------------------------------------
% Aufbau Datensatz
DATA_Q_TI(1,:) = DATA_RAW.Speed_combined_XY.Time;                   % Zeit in s
DATA_Q_TI(2,:) = DATA_RAW.s_AMK_fr_ActualVelocity.Value;            % Drehzahl VA 1/min
DATA_Q_TI(3,:) = DATA_RAW.s_AMK_fr_ActualTorque.Value;              % Moment VA Nm
DATA_Q_TI(4,:) = DATA_RAW.s_AMK_rr_ActualVelocity.Value;            % Drehzahl HA 1/min
DATA_Q_TI(5,:) = DATA_RAW.s_AMK_rr_ActualTorque.Value;              % Moment HA Nm

% FSB‑Anteil über eine Runde berechnen
Anteil_FSB_VL = sum(DATA_Q_TI(2,:) > 11000)/numel(DATA_Q_TI(2,:));
Anteil_FSB_HL = sum(DATA_Q_TI(4,:) > 11000)/numel(DATA_Q_TI(4,:));

% Durchschnittsmomente im/außer FSB
AVG_Moment_FSB_VL      = sum(abs(DATA_Q_TI(3, DATA_Q_TI(2, :) > 11000))) / (Anteil_FSB_VL     * numel(DATA_Q_TI(2,:)));
AVG_Moment_nicht_FSB_VL= sum(abs(DATA_Q_TI(3, DATA_Q_TI(2, :) < 11000))) / ((1-Anteil_FSB_VL)* numel(DATA_Q_TI(2,:)));

AVG_Moment_FSB_HL      = sum(abs(DATA_Q_TI(5, DATA_Q_TI(4, :) > 11000))) / (Anteil_FSB_HL     * numel(DATA_Q_TI(4,:)));
AVG_Moment_nicht_FSB_HL= sum(abs(DATA_Q_TI(5, DATA_Q_TI(4, :) < 11000))) / ((1-Anteil_FSB_HL)* numel(DATA_Q_TI(4,:)));

% Ströme aus Momenten (AMK: 1/0.38 ARMS/Nm)
Motorkonstante_AMK = 1/0.38; % ARMS/Nm
DATA_Q_TI(6,:) = abs(DATA_Q_TI(3,:))*Motorkonstante_AMK; % Strom VA
DATA_Q_TI(7,:) = abs(DATA_Q_TI(5,:))*Motorkonstante_AMK; % Strom HA

% Inverter-Verluste (zwei Varianten)
% Infineon FS200R12PT4 (6 Schalter)
DATA_Q_TI(8,:)  = 6*(0.002589.*DATA_Q_TI(6,:).*DATA_Q_TI(6,:) + 0.557347.*DATA_Q_TI(6,:) + 0.945667); 
DATA_Q_TI(9,:)  = 6*(0.002589.*DATA_Q_TI(7,:).*DATA_Q_TI(7,:) + 0.557347.*DATA_Q_TI(7,:) + 0.945667);

% Powermodul Dresden 10-EY126PB011ME-PJ19F18T
DATA_Q_TI(10,:) = (0.0481.*DATA_Q_TI(6,:).*DATA_Q_TI(6,:) - 0.2247.*DATA_Q_TI(6,:) + 8.9979);
DATA_Q_TI(11,:) = (0.0481.*DATA_Q_TI(7,:).*DATA_Q_TI(7,:) - 0.2247.*DATA_Q_TI(7,:) + 8.9979);

%% ------------------------------------------------------------------------
%  Plots TI
% -------------------------------------------------------------------------
figure(4); tiledlayout(1,2)
nexttile
scatter(DATA_Q_TI(2,:),DATA_Q_TI(3,:), 2); hold on
xline(11000, 'r--'); yline(0, 'black--')
xlabel('Drehzahl in 1/min'); ylabel('Drehmoment in Nm');
xlim([0, 18000]); ylim([-15 30]); grid on; set(gca,'FontSize',14)

nexttile
scatter(DATA_Q_TI(4,:),DATA_Q_TI(5,:), 2); hold on
xline(11000, 'r--'); yline(0, 'black--')
xlabel('Drehzahl in 1/min'); ylabel('Drehmoment in Nm');
xlim([0, 18000]); ylim([-15 30]); grid on; set(gca,'FontSize',14)

% Infineon
figure(5); tiledlayout(1,2)
nexttile
plot(DATA_Q_TI(1,:)-DATA_Q_TI(1,1), DATA_Q_TI(8,:)); hold on
plot([min(time), max(time)], [mean(DATA_Q_TI(8,:)), mean(DATA_Q_TI(8,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('Abwärmeleistung Runde | Infineon FS200R12PT4 VL'); grid on; set(gca,'FontSize',14)

nexttile
plot(time, DATA_Q_TI(9,:)); hold on
plot([min(time), max(time)], [mean(DATA_Q_TI(9,:)), mean(DATA_Q_TI(9,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('Abwärmeleistung Runde | Infineon FS200R12PT4 HL'); grid on; set(gca,'FontSize',14)

% SiC‑Modul
figure(6); tiledlayout(1,2)
nexttile
plot(DATA_Q_TI(1,:)-DATA_Q_TI(1,1), DATA_Q_TI(10,:)); hold on
plot([min(time), max(time)], [mean(DATA_Q_TI(10,:)), mean(DATA_Q_TI(10,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('Abwärmeleistung | 10-EY126PB011ME-PJ19F18T VL'); grid on; set(gca,'FontSize',14)

nexttile
plot(time, DATA_Q_TI(11,:)); hold on
plot([min(time), max(time)], [mean(DATA_Q_TI(11,:)), mean(DATA_Q_TI(11,:))], 'r--', 'LineWidth', 2);
xlabel('Zeit in s'); ylabel('Abwärmeleistung in W'); ylim([0, 350]);
title('Abwärmeleistung Runde | 10-EY126PB011ME-PJ19F18T HL'); grid on; set(gca,'FontSize',14)

%% ------------------------------------------------------------------------
%  Assembly Fahrprofil
% -------------------------------------------------------------------------
Fahrprofil(1,:)  = DATA_RAW.Speed_combined_XY.Time;        % Zeit
Fahrprofil(2,:)  = DATA_RAW.Speed_combined_XY.Value;       % Geschwindigkeit (m/s)
Fahrprofil(3,:)  = Interpolated_Loss_FR;                   % Verluste VL (AMK)
Fahrprofil(4,:)  = Interpolated_Loss_RR;                   % Verluste HL (AMK)
Fahrprofil(5,:)  = DATA_Q_TI(8,:);                         % Inverter VA (Infineon)
Fahrprofil(6,:)  = DATA_Q_TI(9,:);                         % Inverter HA (Infineon)
Fahrprofil(7,:)  = DATA_Q_TI(10,:);                        % Inverter VA (SiC)
Fahrprofil(8,:)  = DATA_Q_TI(11,:);                        % Inverter HA (SiC)
Fahrprofil(9,:)  = Interpolated_Loss_FR_M;                 % Verluste VL Motor94.1
Fahrprofil(10,:) = Interpolated_Loss_RR_M;                 % Verluste HL Motor94.1

No_Laps          = 10;
Rows_BSL_Dataset = size(Fahrprofil,1);
Cols_BSL_Dataset = size(Fahrprofil,2);
Laptime          = Fahrprofil(1,Cols_BSL_Dataset) - Fahrprofil(1,1);

if No_Laps > 1
    for index = 2:No_Laps
        Rows_Dataset = size(Fahrprofil,1); %#ok<NASGU>
        Cols_Dataset = size(Fahrprofil,2);

        Fahrprofil(2,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(2,1:Cols_BSL_Dataset);
        Fahrprofil(3,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(3,1:Cols_BSL_Dataset);
        Fahrprofil(4,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(4,1:Cols_BSL_Dataset);
        Fahrprofil(5,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(5,1:Cols_BSL_Dataset);
        Fahrprofil(6,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(6,1:Cols_BSL_Dataset);
        Fahrprofil(7,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(7,1:Cols_BSL_Dataset);
        Fahrprofil(8,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(8,1:Cols_BSL_Dataset);
        Fahrprofil(9,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset))  = Fahrprofil(9,1:Cols_BSL_Dataset);
        Fahrprofil(10,(Cols_Dataset+1):(Cols_Dataset+Cols_BSL_Dataset)) = Fahrprofil(10,1:Cols_BSL_Dataset);
    end
end

% neue Zeit für Fahrprofil setzen
Fahrprofil(1,:) = 0:(1/Hz):(No_Laps*Laptime+(1/Hz)*(No_Laps-1));
Endtime         = No_Laps * Laptime %#ok<NOPTS>

%% ------------------------------------------------------------------------
%  Downsampling auf 20 Hz
% -------------------------------------------------------------------------
original_datarate  = 200;   % Original (200 Hz)
target_datarate    = 20;    % Ziel (20 Hz)
decimation_factor  = original_datarate / target_datarate;

num_cols = size(Fahrprofil,2);
Fahrprofil_reduziert(1:10,:) = zeros(10,round(num_cols/decimation_factor));
Fahrprofil_reduziert(1:10,:) = Fahrprofil(1:10,1:decimation_factor:end);

%% ------------------------------------------------------------------------
%  Pumpe (Beispielkennlinie)
% -------------------------------------------------------------------------
Pumpenkennlinie_Davis_Craig(1,:) = 0:5:30;                       % Volumenstrom l/min
Pumpenkennlinie_Davis_Craig(1,:) = Pumpenkennlinie_Davis_Craig(1,:)/60;  % l/s
Pumpenkennlinie_Davis_Craig(1,:) = Pumpenkennlinie_Davis_Craig(1,:)/1000;% m^3/s
Pumpenkennlinie_Davis_Craig(2,:) = [0.84321 0.81771 0.76721 0.69171 0.59121 0.46571 0.31521]; % bar

%% ------------------------------------------------------------------------
%  Geschwindigkeitsverteilung (Histogramm)
% -------------------------------------------------------------------------
speeds = DATA_RAW.Speed_combined_XY.Value .* 3.6;  % km/h

classEdges = 0:5:100;
[counts, edges] = histcounts(speeds, classEdges);
relativeFrequencies = (counts / sum(counts)) * 100;

figure;
bar(edges(1:end-1), relativeFrequencies, 'FaceColor', "#0072BD");
xlabel('Geschwindigkeit in km/h', 'FontSize', 14);
ylabel('Relative Häufigkeit in %', 'FontSize', 14);
set(gca, 'FontSize', 14); grid on
xticks(0:10:100);

figure
plot(speeds)
xticks([0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000]);
xticklabels({'0','10','20','30','40','50','60','70','80'});
yticks(0:10:100);
xlabel("Zeit in s")
ylabel("Fahrgeschwindigkeit in km/h")
set(gca, 'FontSize', 12);
grid on
axis([0 16000 0 100]);
