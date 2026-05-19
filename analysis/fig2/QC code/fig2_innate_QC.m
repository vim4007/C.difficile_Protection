%% fig2_innate_QC.m
% Quality control checks for innate immune flow cytometry data
% Run AFTER fig2_innate_analysis.m — requires variables from that script
% Run with MATLAB working directory set to data/raw/flow_cytometry/

%% Overall gating summary
totalCells = height(allMiceWithGates);
didNotPassAnyGate = sum(allMiceWithGates{:, 16:22}, 2) == 0;
nNotGated = sum(didNotPassAnyGate);
nGated = sum(~didNotPassAnyGate);
nWithinRange = height(allCellsWithinRange);

fprintf('\n=== Innate Panel Gating Quality Control ===\n');
fprintf('Total cells loaded:              %d\n', totalCells);
fprintf('Cells passing at least one gate: %d (%.1f%%)\n', nGated, 100*nGated/totalCells);
fprintf('Cells passing NO gate:           %d (%.1f%%)\n', nNotGated, 100*nNotGated/totalCells);
fprintf('Cells within fluorescence range: %d (%.1f%%)\n', nWithinRange, 100*nWithinRange/totalCells);

%% Per mouse breakdown
fprintf('\nPer mouse gating summary:\n');
mice = [1 6 7 11 12 13];
groups = {'UI', 'UI', 'UI', 'Avirulent', 'Avirulent', 'Avirulent'};

for i = 1:length(mice)
    mouseCells = allMiceWithGates(allMiceWithGates.mouse == mice(i), :);
    mouseTotal = height(mouseCells);
    mouseGated = sum(sum(mouseCells{:, 16:22}, 2) > 0);
    fprintf('  Mouse %d (%s): %d total, %d gated (%.1f%%)\n', ...
        mice(i), groups{i}, mouseTotal, mouseGated, 100*mouseGated/mouseTotal);
end

%% Cells per gate per mouse
cellTypeVars = allCellsWithinRange.Properties.VariableNames(16:22);
cellTypeNames = allCellsWithinRange.Properties.VariableNames(16:22);

fprintf('\nCells per gate per mouse:\n');
fprintf('%-10s', 'Mouse');
for j = 1:length(cellTypeNames)
    fprintf('%-18s', cellTypeNames{j});
end
fprintf('\n');

for i = 1:length(mice)
    fprintf('%-10d', mice(i));
    mouseCells = allCellsWithinRange(allCellsWithinRange.mouse == mice(i), :);
    for j = 1:length(cellTypeNames)
        n = sum(mouseCells{:, cellTypeVars{j}});
        fprintf('%-18d', n);
    end
    fprintf('\n');
end

%% Flag low cell count populations
fprintf('\n=== Low Cell Count Warnings (< 30 cells) ===\n');
warned = false;
for i = 1:length(mice)
    mouseCells = allCellsWithinRange(allCellsWithinRange.mouse == mice(i), :);
    for j = 1:length(cellTypeNames)
        n = sum(mouseCells{:, cellTypeVars{j}});
        if n < 30
            fprintf('  WARNING: Mouse %d, %s — only %d cells\n', mice(i), cellTypeNames{j}, n);
            warned = true;
        end
    end
end
