% hydration_by_class.m
clear; clc;

% Folder containing the CSV files
folder = './';  % or specify full path

% Define elastic moduli and corresponding suffixes
E_values = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0];
atom_classes = {'N', 'O', 'C', 'H'};  % NH3+, OH, CH, H
colors = [0.4 0 0.6; 1.0 0.5 0; 1.0 0.8 0; 0.1 0.6 1];  % purple, orange, yellow, blue

% Initialize matrix to store total hydration contributions
hydration_matrix = zeros(length(E_values), length(atom_classes));

for i = 1:length(E_values)
    E = E_values(i);
    suffix = strrep(sprintf('n20_E%.1f', E), '.', '_');
    for j = 1:length(atom_classes)
        atom = atom_classes{j};
        fname = fullfile(folder, sprintf('%s_hydration_%s.csv', suffix, atom));
        if isfile(fname)
            data = readmatrix(fname);
            hydration_matrix(i, j) = sum(data(~isnan(data)));
        else
            warning('Missing file: %s', fname);
        end
    end
end

% Create bar chart
figure('Color', 'w', 'Position', [200, 200, 900, 500]);
bar(E_values, hydration_matrix, 'stacked');
colormap(colors);
legend({'NH3^+', 'Free OH', 'Bonded OH', 'Nonpolar'}, 'Location', 'northwest');
xlabel('Substrate Modulus E (MPa)', 'FontSize', 14);
ylabel('Total Hydration Energy', 'FontSize', 14);
title('Hydration Energy Contribution by Site Class', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);

% Optional: save figure
saveas(gcf, 'hydration_by_class_plot.png');