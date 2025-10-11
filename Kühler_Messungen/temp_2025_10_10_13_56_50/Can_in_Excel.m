%% ============================================================
%  MAT -> Excel mit Einheiten (für 13_08_30)
%  - erste Spalte: Time [s] ab 0
%  - optionales Resampling (gleichmäßiges Δt)
%  - Spaltennamen inkl. angenommener Einheiten
%% ============================================================

% >>> Pfad/Datei
matPath = 'D:\GT\GT\KKL\Matlab_files\Prüfstandprogramm MLAPP StandMaxi_V02\temp_2025_10_10_13_08_30\canlog_2025_10_10_13_08_30.mat';

% Optionen
doSync      = true;     % true = retime auf gleichmäßiges Raster
fixedStep_s = [];       % [] = Median-Takt; z.B. 0.1 für 100 ms

assert(exist(matPath,'file')==2, 'Datei nicht gefunden: %s', matPath);
S = load(matPath);
assert(isfield(S,'TT') && istimetable(S.TT), 'Erwarte timetable TT.');
TT = S.TT;

%% ggf. synchronisieren
if doSync
    t = TT.Properties.RowTimes;
    assert(numel(t)>=2, 'Zu wenige Zeitstempel.');
    if isempty(fixedStep_s), dt = median(diff(t)); else, dt = seconds(fixedStep_s); end
    newTime = (t(1):dt:t(end))';

    % numerische vs. andere Spalten OHNE vartype trennen
    V = TT.Properties.VariableNames;
    isNum = false(1,numel(V));
    for i = 1:numel(V), isNum(i) = isnumeric(TT.(V{i})); end

    TTn = retime(TT(:, isNum), newTime, 'linear');       % numerisch: linear
    if any(~isNum)
        TTo = retime(TT(:,~isNum), newTime, 'previous'); % nicht-numerisch: hold
        TT  = [TTn TTo];
    else
        TT  = TTn;
    end
end

%% Zeit in Sekunden relativ
tsec = seconds(TT.Properties.RowTimes - TT.Properties.RowTimes(1));
tsec = double(tsec(:));

%% Tabelle ohne RowTimes + Excel-freundlich machen
Tvars = timetable2table(TT, 'ConvertRowTimes', false);
for i = 1:width(Tvars)
    if isduration(Tvars.(i)), Tvars.(i) = string(Tvars.(i)); end
    if isdatetime(Tvars.(i)), Tvars.(i) = string(Tvars.(i)); end
end

% Header mit Einheiten
varNames = Tvars.Properties.VariableNames;
header = cell(1, numel(varNames)+1); header{1} = 'Time [s]';
for k = 1:numel(varNames)
    u = guessUnit(varNames{k});
    header{k+1} = ifelse(~isempty(u), sprintf('%s [%s]', varNames{k}, u), varNames{k});
end

% Datenzellen
dataCell = [num2cell(tsec), table2cell(Tvars)];

% Schreiben
[outDir, base, ~] = fileparts(matPath);
outXlsx = fullfile(outDir, [base '_export_units.xlsx']);
writecell(header,   outXlsx, 'Sheet','Data',  'Range','A1');
writecell(dataCell, outXlsx, 'Sheet','Data',  'Range','A2');

% Units-Sheet
unitsList = [{'Variable','Unit'}; [varNames(:), cellfun(@guessUnit, varNames(:), 'uni',0)]];
writecell(unitsList, outXlsx, 'Sheet','Units', 'Range','A1');

fprintf('[OK] Excel geschrieben: %s\n', outXlsx);

%% --- Hilfsfunktionen ---
function out = ifelse(cond, a, b), if cond, out=a; else, out=b; end, end
function u = guessUnit(name)
    n = lower(string(name)); u = '';
    if startsWith(n,"t_") || contains(n,"temp"), u = '°C';
    elseif contains(n,"diffpressure") || contains(n,"pressure") || endsWith(n,"_dp"), u = 'Pa';
    elseif contains(n,"massflow") || contains(n,"mdot") || contains(n,"m_air"), u = 'kg/s';
    elseif contains(n,"volumeflow") || contains(n,"vdot") || contains(n,"v_flow"), u = 'm^3/h';
    elseif contains(n,"rpm") || contains(n,"fanrpm") || contains(n,"speed"), u = 'rpm';
    elseif contains(n,"power") || endsWith(n,"_p"), u = 'W';
    elseif contains(n,"voltage") || startsWith(n,"u_") || contains(n,"_u"), u = 'V';
    elseif contains(n,"current") || startsWith(n,"i_") || contains(n,"_i"), u = 'A';
    elseif contains(n,"duty") || contains(n,"percent") || endsWith(n,"_pct"), u = '%';
    elseif contains(n,"humidity") || contains(n,"rh"), u = '%RH';
    elseif contains(n,"flow") && isempty(u), u = 'SI';
    end
end
