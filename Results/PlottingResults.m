%% code to simulate
clear all
close all

%% load results
S1 = readtable('Scenario1_New.xlsx');
S2 = readtable('Scenario2_New.xlsx');
S3 = readtable('Scenario3_New.xlsx');

%% plot results
figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.MMR, 20, 'Normalization', 'probability', 'FaceColor', 'r')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)
xlim([0.6 1.4])

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S1.MMR, S1.P, 40, S1.C, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('X_{Plagioclase}^{Reacted assemblage}', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%colormap viridis

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S1.MMR, -S1.H/1000, 40, S1.C, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('\Delta H (kJ/kg)', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.MMR, 20, 'Normalization', 'probability', 'FaceColor', 'r')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)
xlim([0.6 1.4])

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S2.MMR, S2.P, 40, S2.C, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('X_{Plagioclase}^{Reacted assemblage}', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%colormap viridis

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S2.MMR, -S2.H/1000, 40, S2.C, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('\Delta H (kJ/kg)', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.MMR, 20, 'Normalization', 'probability', 'FaceColor', 'r')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)
xlim([0.6 1.4])

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S3.MMR, S3.P, 40, S3.O, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('X_{Plagioclase}^{Reacted assemblage}', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%colormap viridis

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S3.MMR, -S3.H/1000, 40, S3.O, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('\Delta H (kJ/kg)', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%% histograms
figure('rend','painters','pos',[20 10 1150 1220])
subaxis(3,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.O,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{React}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.P,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{React}', 'FontSize', 16)

subaxis(3,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.C,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{React}', 'FontSize', 16)

subaxis(3,3,6,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.T,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('T_{initial} (^{o}C)', 'FontSize', 16)

subaxis(3,3,4,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.phi,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('\phi', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,5,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.M,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('M', 'FontSize', 16)

subaxis(3,3,7,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.Oc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{Cryst}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,8,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.Pc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{Cryst}', 'FontSize', 16)

subaxis(3,3,9,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1.Cc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{Cryst}', 'FontSize', 16)

figure('rend','painters','pos',[20 10 1150 1220])
subaxis(3,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.O,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{React}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.P,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{React}', 'FontSize', 16)

subaxis(3,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.C,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{React}', 'FontSize', 16)

subaxis(3,3,6,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.T,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('T_{initial} (^{o}C)', 'FontSize', 16)

subaxis(3,3,4,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.phi,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('\phi', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,5,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.M,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('M', 'FontSize', 16)

subaxis(3,3,7,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.Oc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{Cryst}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,8,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.Pc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{Cryst}', 'FontSize', 16)

subaxis(3,3,9,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S2.Cc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{Cryst}', 'FontSize', 16)

figure('rend','painters','pos',[20 10 1150 1220])
subaxis(3,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.O,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{React}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.P,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{React}', 'FontSize', 16)

subaxis(3,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.C,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{React}', 'FontSize', 16)

subaxis(3,3,6,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.T,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('T_{initial} (^{o}C)', 'FontSize', 16)

subaxis(3,3,4,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.phi,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('\phi', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,5,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.M,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('M', 'FontSize', 16)

subaxis(3,3,7,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.Oc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{Cryst}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,8,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.Pc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{Cryst}', 'FontSize', 16)

subaxis(3,3,9,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.Cc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{Cryst}', 'FontSize', 16)

%% correlation plots
figure('rend','painters','pos',[20 10 450 420])
% sample correlation matrix
r = corrcoef(S1{:,2:end});

% labels
labels = ["MMR", "X_{Olivine}", "X_{Plagioclase}", "X_{Clinopyroxene}", "phi", "M", "T", "X_{OlCryst}", "X_{PlagCryst}", "X_{ClinoCryst}",  "delta H"];
% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(x(:), y(:), 300.*abs(r(:)), r(:), 'filled', 'MarkerFaceAlpha', 1)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text(1:n, (n:-1:1) + 0.5, labels, 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1, n, 1), labels, ...
    'HorizontalAlignment', 'right', 'Rotation', 270)
h = gca;
colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
%colormap('viridis')

figure('rend','painters','pos',[20 10 450 420])
% sample correlation matrix
r = corrcoef(S2{:,2:end});

% labels
labels = ["MMR", "X_{Olivine}", "X_{Plagioclase}", "X_{Clinopyroxene}", "phi", "M", "T", "X_{OlCryst}", "X_{PlagCryst}", "X_{ClinoCryst}",  "delta H"];
% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(x(:), y(:), 300.*abs(r(:)), r(:), 'filled', 'MarkerFaceAlpha', 1)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text(1:n, (n:-1:1) + 0.5, labels, 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1, n, 1), labels, ...
    'HorizontalAlignment', 'right', 'Rotation', 270)
h = gca;
colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
%colormap('viridis')

figure('rend','painters','pos',[20 10 450 420])
% sample correlation matrix
r = corrcoef(S3{:,2:end});

% labels
labels = ["MMR", "X_{Olivine}", "X_{Plagioclase}", "X_{Clinopyroxene}", "phi", "M", "T", "X_{OlCryst}", "X_{PlagCryst}", "X_{ClinoCryst}",  "delta H"];
% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(x(:), y(:), 300.*abs(r(:)), r(:), 'filled', 'MarkerFaceAlpha', 1)
% enclose markers in a grid
xl = [1:n+1;repmat(n+1, 1, n+1)];
xl = [xl(:, 1), xl(:, 1:end-1)];
yl = repmat(n+1:-1:1, 2, 1);
line(xl, yl, 'color', 'k') % horizontal lines
line(yl, xl, 'color', 'k') % vertical lines
% show labels
text(1:n, (n:-1:1) + 0.5, labels, 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1, n, 1), labels, ...
    'HorizontalAlignment', 'right', 'Rotation', 270)
h = gca;
colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.9;
axis(h, 'equal')
%colormap('viridis')

%% clinopyroxene plots
Cpx = readtable('CpxRes_FeOt.xlsx');
CpxFC = readtable('CpxFC.xlsx');
AB = readtable('AB_mineral_Major.xlsx');

figure('rend','painters','pos',[20 10 450 420])
hold on
xlim([0.6 0.95])
ylim([0 2.5])
h=hexscatter((AB.MgO(AB.Mineral == "Cpx")./40.3044)./...
    (AB.MgO(AB.Mineral == "Cpx")./40.3044+AB.FeO(AB.Mineral == "Cpx")./71.844), ...
    AB.TiO2(AB.Mineral == "Cpx"), 'xlim', [0.60 0.95],...
    'ylim', [0.0 2.5], 'showZeros', false, 'res', 25);
colormap(flipud(gray))
plot((CpxFC.MgO_Cpx./40.3044)./(CpxFC.MgO_Cpx./40.3044+CpxFC.FeOt_Cpx./71.844), CpxFC.TiO2_Cpx, 'o', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', 'None', 'MarkerSize', 5)
plot((Cpx.MgO_Cpx./40.3044)./...
    (Cpx.MgO_Cpx./40.3044+Cpx.FeOt_Cpx./71.844),...
    Cpx.TiO2_Cpx, 'o', 'MarkerFaceColor', [0.56 0 0], 'MarkerEdgeColor', 'None', 'MarkerSize', 4)
box on

xlabel('Mg#', 'FontSize', 16)
ylabel('TiO_{2} (wt%)', 'FontSize', 16)

%% Mineral changes plots
figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
h = gobjects(3, 1);
h(1) = plot(S2.O, S2.Oc, 'ok', 'MarkerFaceColor', 'yellow', ...
    'MarkerSize', 5, 'DisplayName', 'Scenario 2');
h(2) = plot(S3.O, S3.Oc, 'ok', 'MarkerFaceColor', 'blue', ...
    'MarkerSize', 5, 'DisplayName', 'Scenario 3');
h(3) = plot(S1.O, S1.Oc, 'ok', 'MarkerFaceColor', 'red', ...
    'MarkerSize', 5, 'DisplayName', 'Scenario 1');
plot([0 1], [0 1], '--k', 'LineWidth', 3)
box on
xlabel('X_{Olivine}^{React}', 'FontSize', 16)
ylabel('X_{Olivine}^{Cryst}', 'FontSize', 16)
legend(h([3 1 2]),'Location', 'southeast')

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
plot(S2.P, S2.Pc, 'ok', 'MarkerFaceColor', 'yellow', ...
    'MarkerSize', 5)
plot(S3.P, S3.Pc, 'ok', 'MarkerFaceColor', 'blue', ...
    'MarkerSize', 5)
plot(S1.P, S1.Pc, 'ok', 'MarkerFaceColor', 'red', ...
    'MarkerSize', 5)
plot([0 1], [0 1], '--k', 'LineWidth', 3)
box on
xlabel('X_{Plagioclase}^{React}', 'FontSize', 16)
ylabel('X_{Plagioclase}^{Cryst}', 'FontSize', 16)

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
plot(S2.C, S2.Cc, 'ok', 'MarkerFaceColor', 'yellow', ...
    'MarkerSize', 5)
plot(S3.C, S3.Cc, 'ok', 'MarkerFaceColor', 'blue', ...
    'MarkerSize', 5)
plot(S1.C, S1.Cc, 'ok', 'MarkerFaceColor', 'red', ...
    'MarkerSize', 5)
plot([0 1], [0 1], '--k', 'LineWidth', 3)
xlabel('X_{Clinopyroxene}^{React}', 'FontSize', 16)
ylabel('X_{Clinopyroxene}^{Cryst}', 'FontSize', 16)
box on

%% load results
S1_Dick = readtable('Scenario1_Dick2000Start.xlsx');
S3_Dick = readtable('Scenario3_Dick2000Start.xlsx');

S1_Dick = rmmissing(S1_Dick);
S3_Dick = rmmissing(S3_Dick);

S3_Dick(S3_Dick.T > 1205, :) = [];

figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.MMR, 20, 'Normalization', 'probability', 'FaceColor', 'r')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)
xlim([0.6 1.4])

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S1_Dick.MMR, S1_Dick.P, 40, S1_Dick.C, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('X_{Plagioclase}^{Reacted assemblage}', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%colormap viridis

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S1_Dick.MMR, -S1_Dick.H/1000, 40, S1_Dick.C, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('\Delta H (kJ/kg)', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

figure('rend','painters','pos',[20 10 1250 320])
subaxis(1,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3_Dick.MMR, 20, 'Normalization', 'probability', 'FaceColor', 'r')
xlabel('Melt Mass Ratio', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)
xlim([0.6 1.4])

subaxis(1,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S3_Dick.MMR, S3_Dick.P, 40, S3_Dick.O, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('X_{Plagioclase}^{Reacted assemblage}', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%colormap viridis

subaxis(1,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
scatter(S3_Dick.MMR, -S3_Dick.H/1000, 40, S3_Dick.O, 'o', 'filled', 'MarkerEdgeColor', 'k')
xlabel('Melt Mass Ratio', 'FontSize', 16)
xlim([0.6 1.4])
ylabel('\Delta H (kJ/kg)', 'FontSize', 16)
a = colorbar;
a.Label.String = 'X_{Clinopyroxene}^{Reacted assemblage}';
a.Label.FontSize = 16;
a.Location = 'northoutside';
a.Position = a.Position + [0, 0.18, 0,0];

%% histograms
figure('rend','painters','pos',[20 10 1150 1220])
subaxis(3,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.O,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{React}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.P,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{React}', 'FontSize', 16)

subaxis(3,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.C,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{React}', 'FontSize', 16)

subaxis(3,3,6,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.T,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('T_{initial} (^{o}C)', 'FontSize', 16)

subaxis(3,3,4,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.phi,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('\phi', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,5,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.M,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('M', 'FontSize', 16)

subaxis(3,3,7,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.Oc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{Cryst}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,8,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.Pc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{Cryst}', 'FontSize', 16)

subaxis(3,3,9,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S1_Dick.Cc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{Cryst}', 'FontSize', 16)

figure('rend','painters','pos',[20 10 1150 1220])
subaxis(3,3,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.O,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{React}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.P,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{React}', 'FontSize', 16)

subaxis(3,3,3,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.C,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{React}', 'FontSize', 16)

subaxis(3,3,6,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.T,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('T_{initial} (^{o}C)', 'FontSize', 16)

subaxis(3,3,4,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.phi,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('\phi', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,5,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.M,15, 'FaceColor', 'blue', 'FaceAlpha', 0.6)
xlabel('M', 'FontSize', 16)

subaxis(3,3,7,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.Oc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Olivine}^{Cryst}', 'FontSize', 16)
ylabel('Count', 'FontSize', 16)

subaxis(3,3,8,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.Pc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Plagioclase}^{Cryst}', 'FontSize', 16)

subaxis(3,3,9,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12)
hold on
box on
histogram(S3.Cc,15, 'FaceColor', 'red', 'FaceAlpha', 0.6)
xlabel('X_{Clinopyroxene}^{Cryst}', 'FontSize', 16)

%% Melt flux calculations
n = 27.39; %Pa S - calculated using the model of Giordano et al. 2008 using the melt composition from the MELTS FC model at 1180oC
p_mush = 3288.9*0.05 + 2661.1*0.55 + 3223.1*0.40; % Ol, Plg, Cpx density and modal abundance in mush
p_melt = 2665.4; % kg/m3 from MELTS model at 1180 oC
g = 9.81; % Gravitational acceleration m/s2
r_mush = linspace(100,5000,1000); % radius of mush system in m
r_channel = linspace(0.01, 0.5, 1000); % radius of melt channel in m
phi = [0.05 0.1 0.2 0.4]; % mush porosity
k = (0.001^2).*0.002.*(phi.^3);

for i = 1:4
    Q_mush(i,:) = ((k(i).*pi().*r_mush.^2)./n).*((p_mush - p_melt).*g);
end

Q_channel = (pi()/8).*((p_mush - p_melt).*g.*(r_channel.^4))./n;

figure('rend','painters','pos',[20 10 800 350])
subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12, 'YScale', 'log')
hold on
box on
plot(r_mush, Q_mush(1,:), '-b')
plot(r_mush, Q_mush(2,:), '-r')
plot(r_mush, Q_mush(3,:), '-c')
plot(r_mush, Q_mush(4,:), '-g')
ylim([10^(-7) 10^2])
xlabel('Mush radius (m)', 'FontSize', 16)
ylabel('Q (m^{3}/s)', 'FontSize', 16)

subaxis(1,2,2,'SpacingVert',0.08,'SpacingHoriz',0.06)
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12, 'YScale', 'log', 'YAxisLocation', 'right')
hold on
box on
plot(r_channel, Q_channel, '-k')
ylim([10^(-7) 10^2])
xlabel('Channel radius (m)', 'FontSize', 16)
ylabel('Q (m^{3}/s)', 'FontSize', 16)