%% Compairison between FEvoR and Thor time step
    % make sure to run FEvoR/bin/compairisonStep before this script. '!' call to
    % it doesn't work.
    
close all, clear all, clc

addpath ~/Documents/Programs/Thor/trunk/

% perform a time step
Test.timing;

iN   = [sin(icdist.theta).*cos(icdist.phi) sin(icdist.theta).*sin(icdist.phi) cos(icdist.theta)];
fN   = [sin(fcdist.theta).*cos(fcdist.phi) sin(fcdist.theta).*sin(fcdist.phi) cos(fcdist.theta)];

% make data files from time step
 idata = [(1:SET.numbcrys)',iN, icdist.size, icdist.dislDens, iSET.to, iSET.Do];
 fdata = [(1:SET.numbcrys)',fN, fcdist.size, fcdist.dislDens, fSET.to, fSET.Do];
    
% write data file
ifilename = 'compairisonDist.csv';
ffilename = 'compairisonDist_stepped_matlab.csv';
cfilename = 'compairisonDist_stepped_cpp.csv';

fid = fopen(ifilename, 'w');
fprintf(fid, '# Crystal, C-Axis (x), C-Axis (y), C-Axis (z), Size (m), Disl. dens. (1/m^2), Last recr. time (s), Size at last recr. (m)\n');
fclose(fid);

dlmwrite(ifilename, idata, '-append', 'precision', '%.6f', 'delimiter', ',');

fid = fopen(ffilename, 'w');
fprintf(fid, '# Crystal, C-Axis (x), C-Axis (y), C-Axis (z), Size (m), Disl. dens. (1/m^2), Last recr. time (s), Size at last recr. (m)\n');
fclose(fid);

dlmwrite(ffilename, fdata, '-append', 'precision', '%.6f', 'delimiter', ',');

% cpp step
% !../bin/compairisonStep > compairisonDist_stepped_cpp.csv
    % stupid this doesn't work.

cpp = importdata(cfilename,',',1);
cdata = cpp.data;

cN = cdata(:,2:4);

%% flip N(:,3) < 0
imsk = iN(:,3) < 0;
fmsk = fN(:,3) < 0;
cmsk = cN(:,3) < 0;

iN(imsk,:) = -iN(imsk,:);
fN(fmsk,:) = -fN(fmsk,:);
cN(cmsk,:) = -cN(cmsk,:);

idata(:,2:4) = iN;
fdata(:,2:4) = fN;
cdata(:,2:4) = cN;


%% plot data

display(bedot);
display(NPOLY);
display(NMIGRE);

MSize = 90;
LWidth = 3;

figure;
scatter3(idata(1:100:end,2),idata(1:100:end,3),idata(1:100:end,4), MSize, 'xb','LineWidth',LWidth)
hold on

scatter3(fdata(1:100:end,2),fdata(1:100:end,3),fdata(1:100:end,4), 'or','LineWidth',LWidth)
hold on

scatter3(cdata(1:100:end,2),cdata(1:100:end,3),cdata(1:100:end,4), 'og','LineWidth',LWidth)
hold off

axis([-1,1,-1,1,-1,1])
view([0,90])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEvoR: Fabric Evolution with Recrystallization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2009-2014  Joseph H Kennedy
%
% This file is part of FEvoR.
%
% FEvoR is free software: you can redistribute it and/or modify it under the 
% terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or (at your option) any later 
% version.
%
% FEvoR is distributed in the hope that it will be useful, but WITHOUT ANY 
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details.
%
% You should have received a copy of the GNU General Public License along 
% with FEvoR.  If not, see <http://www.gnu.org/licenses/>.
%
% Additional permission under GNU GPL version 3 section 7
%
% If you modify FEvoR, or any covered work, to interface with
% other modules (such as MATLAB code and MEX-files) available in a
% MATLAB(R) or comparable environment containing parts covered
% under other licensing terms, the licensors of FEvoR grant
% you additional permission to convey the resulting work.

