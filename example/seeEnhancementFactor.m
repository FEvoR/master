wpp = importdata('eDist_wk_initial.csv',',',1);
wdata = wpp.data;
%# Crystal, C-Axis (x), C-Axis (y), C-Axis (z), Size (m), Disl. dens. (1/m^2), Last recr. time (s), Size at last recr. (m)
wN = wdata(:,2:4)';

ipp = importdata('eDist_iso_initial.csv',',',1);
idata = ipp.data;
%# Crystal, C-Axis (x), C-Axis (y), C-Axis (z), Size (m), Disl. dens. (1/m^2), Last recr. time (s), Size at last recr. (m)
iN = idata(:,2:4)';

wN = wN';
iN = iN';

% make sure all upper hemisphere
wmsk = wN(:,3) < 0;
imsk = iN(:,3) < 0;

wN(cmsk,:) = -wN(cmsk,:);
iN(fmsk,:) = -iN(fmsk,:);

figure;
Sims.Analyze.equalAreaPoint(wN)
figure;
Sims.Analyze.equalAreaPoint(iN)
