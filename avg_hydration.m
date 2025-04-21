files = dir(' .txt'); % import txt file

E_values = [];
hydration_values = [];

n_units = ;  % define chain length

for i = 1:length(files)
    filename = files(i).name;
    fprintf('\n🔍 Processing: %s\n', filename);

    % === Robustly parse E value from filename ===
    % Expect format like: E0_5 → 0.5, E1_75 → 1.75
    match = regexp(filename, 'E(\d+)_?(\d*)', 'tokens');
    if isempty(match)
        warning('⚠️ Could not parse E from: %s', filename);
        continue;
    end
    digits = match{1};
    whole = str2double(digits{1});
    fractional = str2double(digits{2});
    if isnan(fractional)
        fractional = 0;
    end
    % Convert: E1_25 → 1.25
    factor = 10^(-length(digits{2}));
    E_val = whole + fractional * factor;
    fprintf('✅ Parsed E = %.2f\n', E_val);

    % === Read hydration value from inside the file ===
    hydration_val = NaN;
    try
        fileText = fileread(filename);
        match = regexp(fileText, 'max_hydration:\s*([\d\.Ee+-]+)', 'tokens');
        if ~isempty(match)
            hydration_val = str2double(match{1}{1});
            fprintf('✅ Hydration: %.4f\n', hydration_val);
        else
            fprintf('⚠️ max_hydration not found in %s\n', filename);
            continue;
        end
    catch
        fprintf('❌ Error reading: %s\n', filename);
        continue;
    end

    % === Store values ===
    E_values(end+1) = E_val;
    hydration_values(end+1) = hydration_val;
end

% === Final check ===
if isempty(E_values)
    error('❌ No valid hydration data found.');
end

% === Sort and average ===
[E_values_sorted, idx] = sort(E_values);
hydration_sorted = hydration_values(idx);
avg_hydration_per_monomer = hydration_sorted / n_units;

% === Table ===
disp(table(E_values_sorted', avg_hydration_per_monomer', ...
    'VariableNames', {'E', 'AvgHydration'}))

% === Plot ===
figure('Units', 'inches', 'Position', [0, 0, 13, 8], 'Color', 'w');
scatter(E_values_sorted, avg_hydration_per_monomer, 120, 'b', 'filled');
hold on;

% Polynomial fit (3rd degree)
p = polyfit(E_values_sorted, avg_hydration_per_monomer, 3);
x_fit = linspace(min(E_values_sorted), max(E_values_sorted), 200);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'k--', 'LineWidth', 2);

xlabel('Substrate Modulus E (MPa)', 'FontSize', 18);
ylabel('Avg. Hydration per Monomer', 'FontSize', 18);
title('Hydration Scaling with PDMS Modulus (n = 20)', 'FontSize', 20);
grid on;
set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
