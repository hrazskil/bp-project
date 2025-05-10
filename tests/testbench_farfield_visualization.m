% Load far-field data
dataH = load('tests/test_structure_2/dataHorizontal.tsv');  % [angle, intensity]
dataV = load('tests/test_structure_2/dataVertical.tsv');    % [angle, intensity]

anglesH = dataH(:,1);   % radians
I_H = dataH(:,2);       % power/intensity

anglesV = dataV(:,1);
I_V = dataV(:,2);


figure;

subplot(1,2,1);
polarplot(anglesH, sqrt(I_H), 'r');  % sqrt for E-field magnitude if needed
title('Horizontal Far-Field Pattern');

subplot(1,2,2);
polarplot(anglesV, sqrt(I_V), 'b');
title('Vertical Far-Field Pattern');


