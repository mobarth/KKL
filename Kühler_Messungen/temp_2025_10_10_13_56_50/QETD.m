%% ============================================================
%  Q–ETD Kennfeld (mit robustem Auto-Mapping der Variablen)
%  - lädt timetable TT aus MAT (v7.3 ok)
%  - mappt Kanalnamen robust (resolveVar)
%  - Q = m_dot_c * cp_c * (T_c,in - T_c,out)
%  - ETD = T_c,in - T_air,in
%  - Binning nach Luftmassenstrom, lineare Fits pro Bin
%  - Export CSV
%  Autor: (du)
%% ============================================================

%% ======= Pfad zur MAT-Datei =================================
matPath = 'D:\GT\GT\KKL\Matlab_files\Prüfstandprogramm MLAPP StandMaxi_V02\temp_2025_10_10_13_56_50\canlog_2025_10_10_13_56_50.mat';

%% ======= Laden ==============================================
S = load(matPath);
assert(isfield(S,'TT') && istimetable(S.TT), 'Erwarte timetable TT in der MAT-Datei.');
TT = S.TT;

%% ======= Physik-Parameter & Einheiten =======================
rho_coolant = 1020;     % Dichte [kg/m^3] (Wasser-Glykol 50/50 ~1020)
cp_coolant  = 3800;     % Wärmekapazität [J/(kg*K)] (~3.8 kJ/kgK)

air_mdot_unit = 'kg/h'; % Einheit Luftmassenstrom in TT: 'kg/h' oder 'kg/s'
Vdot_unit     = 'L/min';% Einheit Kühlmittel-Volumenstrom (falls genutzt)

%% ======= Robustes Auto-Mapping ==============================
% -> Passe Kandidaten bei Bedarf an deine Kanalnamen an
var_Tc_in   = resolveVar(TT, {'T_tc_15','T_c_in','Tc_in','CoolantIn','T_cool_in'}, 'Kühlmittel-Eintritt');
var_Tc_out  = resolveVar(TT, {'T_tc_9','T_tc_09','T_c_out','Tc_out','CoolantOut','T_cool_out'}, 'Kühlmittel-Austritt');
var_Tair_in = resolveVar(TT, {'T_int_0','T_air_in','T_in_air','AmbientIn','T_amb_in'}, 'Luft-Eintritt');
var_mdot_air= resolveVar(TT, {'MassFlow','m_air','mdot_air','AirMassFlow'}, 'Luftmassenstrom');

% Kühlmittelmassenstrom: entweder direkt...
var_mdot_coolant = tryResolve(TT, {'CoolantMassFlow','m_coolant','mdot_c','m_dot_coolant'});
% ...oder aus Volumenstrom ableiten (nur, wenn oben leer blieb)
var_Vdot_coolant = "";
if var_mdot_coolant == ""
    var_Vdot_coolant = resolveVar(TT, {'CoolantVolumeFlow','V_coolant','Vdot_c','VolFlow_cool'}, 'Kühlmittel-Volumenstrom');
end

%% ======= Daten holen & Einheiten normieren ==================
Tc_in   = TT.(var_Tc_in);
Tc_out  = TT.(var_Tc_out);
Tair_in = TT.(var_Tair_in);

mdot_air = TT.(var_mdot_air);
switch lower(air_mdot_unit)
    case 'kg/h', mdot_air = mdot_air/3600; % -> kg/s
    case 'kg/s' % ok
    otherwise, error('Unbekannte Luftmassenstrom-Einheit: %s', air_mdot_unit);
end

if var_mdot_coolant ~= ""
    mdot_c = TT.(var_mdot_coolant); % kg/s
else
    Vdot = TT.(var_Vdot_coolant);
    switch lower(Vdot_unit)
        case 'l/min',  Vdot_m3s = (Vdot./1000)/60;
        case 'm^3/h',  Vdot_m3s = Vdot/3600;
        case 'm^3/s',  Vdot_m3s = Vdot;
        otherwise, error('Unbekannte Volumenstrom-Einheit: %s', Vdot_unit);
    end
    mdot_c = rho_coolant .* Vdot_m3s; % kg/s
end

%% ======= Q & ETD ============================================
dT_c = Tc_in - Tc_out;           % K
ETD  = Tc_in - Tair_in;          % K
Q    = mdot_c .* cp_coolant .* dT_c;  % W

valid = isfinite(Q) & isfinite(ETD) & isfinite(mdot_air);
Q      = Q(valid);
ETD    = ETD(valid);
md_air = mdot_air(valid);

%% ======= Binning & Fits ====================================
nbins = 5;
edges = quantile(md_air, linspace(0,1,nbins+1));
binIdx = discretize(md_air, edges);

figure('Name','Q–ETD Kennfeld','Color','w'); hold on; grid on;
cmap = lines(nbins);
leg = strings(1,nbins);
fits = struct('a',[],'b',[],'md_air_mean',[]);

for b = 1:nbins
    sel = (binIdx == b);
    if ~any(sel), continue; end
    x = ETD(sel);          % K
    y = Q(sel)/1000;       % kW
    scatter(x, y, 12, cmap(b,:), 'filled', 'MarkerFaceAlpha', 0.25);
    p = polyfit(x, y, 1);  % y = a*x + b
    xf = linspace(min(x), max(x), 100);
    plot(xf, polyval(p, xf), 'Color', cmap(b,:), 'LineWidth', 1.8);

    md_mean = mean(md_air(sel)); % kg/s
    leg(b) = sprintf('m_{air}≈%.3f kg/s | y=%.3f·ETD%+ .3f', md_mean, p(1), p(2));
    fits(b).a = p(1); fits(b).b = p(2); fits(b).md_air_mean = md_mean;
end
xlabel('ETD = T_{c,in} - T_{air,in}  [K]');
ylabel('Q  [kW]');
title('Q–ETD Kennfeld (Kühlmittel-Seite)');
legend(leg(leg~=""), 'Location','northwest');
hold off;

%% ======= Export =============================================
outDir  = fileparts(matPath);
writetable(table(ETD(:), Q(:), md_air(:), ...
    'VariableNames', {'ETD_K','Q_W','m_air_kgps'}), ...
    fullfile(outDir, 'QETD_points.csv'));

a = [fits.a]'; b = [fits.b]'; mdm = [fits.md_air_mean]';
writetable(table(mdm, a, b, edges(1:end-1)', edges(2:end)', ...
    'VariableNames', {'m_air_mean_kgps','slope_kWperK','offset_kW','bin_lo','bin_hi'}), ...
    fullfile(outDir, 'QETD_fits.csv'));

fprintf('[OK] Export:\n  %s\n  %s\n', ...
    fullfile(outDir,'QETD_points.csv'), fullfile(outDir,'QETD_fits.csv'));

%% ======= Hilfsfunktionen ===================================
function name = resolveVar(TT, candidates, hint)
% Wählt die erste vorhandene Variable aus den Kandidaten.
% Gibt bei Nichtfund einen hilfreichen Fehler mit ähnlichen Namen aus.
    V = string(TT.Properties.VariableNames);
    name = "";
    for c = string(candidates)
        if any(V == c), name = c; return; end
    end
    % Ähnliche Namen per contains (grobe Heuristik):
    key = extractBefore(string(candidates(1))+"_","_"); % Stamm bis zum 1. '_'
    if key == "", key = candidates(1); end
    sugg = V(contains(V, key));
    if isempty(sugg)
        % Fallback: zeige alle Variablen, die mit erstem Kandidaten anfangen
        base = string(candidates(1));
        sugg = V(startsWith(V, extractBefore(base+"_","_")));
    end
    error('Variable für %s nicht gefunden. Versucht: %s\nÄhnliche Variablen in TT: %s', ...
        hint, strjoin(string(candidates), ', '), strjoin(sugg, ', '));
end

function name = tryResolve(TT, candidates)
% Wie resolveVar, aber ohne Fehler. Gibt "" zurück, wenn nichts passt.
    V = string(TT.Properties.VariableNames);
    name = "";
    for c = string(candidates)
        if any(V == c), name = c; return; end
    end
end
