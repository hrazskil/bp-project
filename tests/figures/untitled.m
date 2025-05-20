clc
clear


fig = openfig('fig_norObjectF_rad_obs_points.fig');  % or any of your .fig files
ax = gca;
lines = findall(ax, 'Type', 'line');

% Pick the line you want
X = lines(1).XData(:);
Y = lines(1).YData(:);


% Your Y vector
minY = 0
maxY = 4.5

% Generate 5 linearly spaced tick values
yticks_lin = linspace(minY, maxY, 7);

% Format each number to 4 decimal places
ytick_strs = arrayfun(@(v) sprintf('%.4f', v), yticks_lin, 'UniformOutput', false);

% Join into LaTeX-compatible string
ytick_latex = sprintf('ytick = {%s},', strjoin(ytick_strs, ','));

% Display the result
disp(ytick_latex);


% Your Y vector
minX = min(0)
maxX = max(6000)

% Generate 5 linearly spaced tick values
xticks_lin = linspace(minX, maxX, 7);

% Format each number to 4 decimal places
xtick_strs = arrayfun(@(v) sprintf('%.4f', v), xticks_lin, 'UniformOutput', false);

% Join into LaTeX-compatible string
xtick_latex = sprintf('ytick = {%s},', strjoin(xtick_strs, ','));

% Display the result
disp(xtick_latex);

T = table(X, Y);
% Save it without header
writetable(T, 'C:\Users\kilia\Plocha\gitHub\texStudio\Figures\fig_norObjectF_rad_obs_points.tsv', ...
           'Delimiter', '\t', 'FileType', 'text', 'WriteVariableNames', false);