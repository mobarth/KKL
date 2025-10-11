%% === Pfad zur MAT-Datei anpassen ===
matPath = 'D:\GT\GT\KKL\Matlab_files\Prüfstandprogramm MLAPP StandMaxi_V02\temp_2025_10_10_13_45_28\canlog_2025_10_10_13_45_28.mat';

%% === Laden und prüfen ===
S = load(matPath);
assert(isfield(S,'TT') && istimetable(S.TT), 'Erwarte timetable TT in der MAT-Datei.');
TT = S.TT;

% Variablen prüfen
varsNeeded = {'MassFlow','T_tc_15'};
for v = varsNeeded
    assert(ismember(v{1}, TT.Properties.VariableNames), ...
        'Variable "%s" nicht in TT gefunden.', v{1});
end

%% === Zeitachse in Sekunden ab t0 ===
tsec = timeToSeconds(TT);

% Gültige Daten filtern (optional, falls NaNs)
valid = ~isnan(TT.MassFlow) & ~isnan(TT.T_tc_15) & ~isnan(tsec);
tsec = tsec(valid);
mass = TT.MassFlow(valid);
tc15 = TT.T_tc_15(valid);

%% === Plot: 2 y-Achsen, unterschiedliche Farben ===
figure('Name','MassFlow & T\_tc\_15 (gemeinsame Zeit, zwei Y-Achsen)','Color','w');
hold on; grid on;

% Farben (RGB) – bei Bedarf anpassen
colMass = [0 0.45 0.74];     % Blau
colTc15 = [0.85 0.33 0.10];  % Orange

yyaxis left
p1 = plot(tsec, mass, 'LineWidth', 1.4, 'Color', colMass);
ylabel(labelWithUnit(TT,'MassFlow','MassFlow'));

yyaxis right
p2 = plot(tsec, tc15, 'LineWidth', 1.4, 'LineStyle','--', 'Color', colTc15);
ylabel(labelWithUnit(TT,'T_tc_15','T\_tc\_15'));

xlabel('Zeit [s]');
title('MassFlow und T\_tc\_15 über Zeit');
xlim([0, max(tsec)]);
legend([p1 p2], {'MassFlow','T\_tc\_15'}, 'Location','best');
hold off;

%% === Hilfsfunktionen ===
function ylab = labelWithUnit(TT,varName,pretty)
    % Baut Achsenlabel inkl. Einheit aus timetable-Properties (falls vorhanden)
    ylab = pretty;
    try
        units = TT.Properties.VariableUnits;
        if ~isempty(units)
            idx = find(strcmp(TT.Properties.VariableNames, varName), 1);
            if ~isempty(idx) && ~isempty(units{idx})
                ylab = sprintf('%s [%s]', pretty, units{idx});
            end
        end
    catch
        % ignorieren, wenn Properties fehlen
    end
end

function tsec = timeToSeconds(TT)
    % Liefert relativen Zeitvektor in Sekunden (t(1) -> 0 s)
    try
        t = TT.Properties.RowTimes;
    catch
        if ismember('Time', TT.Properties.DimensionNames)
            t = TT.Time;
        else
            error('Konnte die Zeitachse der timetable nicht ermitteln.');
        end
    end

    if isdatetime(t) || isduration(t)
        tsec = seconds(t - t(1));
    elseif isnumeric(t)
        tsec = t - t(1);
    else
        error('Nicht unterstützter Zeittyp: %s', class(t));
    end
    tsec = double(tsec(:));
end
