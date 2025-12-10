clc;
clear;
close all;

% ------------------ SLL=-20 ------------------
F20 = readmatrix("results_SLL20_objectivePareto.xlsx");
I20 = readmatrix("results_SLL20_currentsPareto.xlsx");

[err20, id20] = min(F20(:,1));
bestI20 = I20(id20,:);
D20     = -F20(id20,2);

fprintf("\n=== SLL=-20 dB ===\n");
fprintf(" Error : %8.3e\n", err20);
fprintf(" D[dBi]: %8.3f\n", 10*log10(D20));
fprintf(" Curr  : [%6.3f %6.3f %6.3f %6.3f]\n", bestI20);


% ------------------ SLL=-30 ------------------
F30 = readmatrix("results_SLL30_objectivePareto.xlsx");
I30 = readmatrix("results_SLL30_currentsPareto.xlsx");

[err30, id30] = min(F30(:,1));
bestI30 = I30(id30,:);
D30     = -F30(id30,2);

fprintf("\n=== SLL=-30 dB ===\n");
fprintf(" Error : %8.3e\n", err30);
fprintf(" D[dBi]: %8.3f\n", 10*log10(D30));
fprintf(" Curr  : [%6.3f %6.3f %6.3f %6.3f]\n", bestI30);


% ------------------ SLL=-40 ------------------
F40 = readmatrix("results_SLL40_objectivePareto.xlsx");
I40 = readmatrix("results_SLL40_currentsPareto.xlsx");

[err40, id40] = min(F40(:,1));
bestI40 = I40(id40,:);
D40     = -F40(id40,2);

fprintf("\n=== SLL=-40 dB ===\n");
fprintf(" Error : %8.3e\n", err40);
fprintf(" D[dBi]: %8.3f\n", 10*log10(D40));
fprintf(" Curr  : [%6.3f %6.3f %6.3f %6.3f]\n", bestI40);
