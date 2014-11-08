%% script to test Thor results line by line for compairison to FEvoR

temperature = -10; % Celsius
stress = [ 10000,     0, 10000;...
               0,     0,     0;...
           10000,     0,-10000];

dt = 1000; % years
addpath ~/Documents/Programs/Thor/trunk/


theta = pi/3;
phi   = pi/3;

n = 3;
R = 0.008314472; % units: kJ K^{-1} mol^{-1}
beta = 630.0;    % from Thors 2001 paper (pg 510, above eqn 16)

if (temperature > -10)
    Q = 115;
else
    Q = 60;
end
A = 3.5e-25*beta*exp(-(Q/R)*(1.0/(273.15+temperature)-1.0/263.15)); % units: s^{-1} Pa^{-n}

% sines and cosines so calculation only has to be preformed once
st = sin(theta); ct = cos(theta);
sp = sin(phi);   cp = cos(phi); sq3 = sqrt(3);

% Basal plane vectors
B1 = [ct.*cp/3, ct.*sp/3, -st/3]; % -
B2 = [(-ct.*cp - sq3.*sp)/6, (-ct.*sp+sq3.*cp)/6, st/6]; % -
B3 = [(-ct.*cp + sq3.*sp)/6, (-ct.*sp-sq3.*cp)/6, st/6]; % -

% C-axis orientation
cAxis = [st.*cp, st.*sp, ct]; % -

shmidt1 = B1'*cAxis;
shmidt2 = B2'*cAxis;
shmidt3 = B3'*cAxis;

rss1 = sum(sum(shmidt1.*stress));
rss2 = sum(sum(shmidt2.*stress));
rss3 = sum(sum(shmidt3.*stress));
Mrss = sqrt( sum( ( B1*rss1 + B2*rss2+B3*rss3 ).^2 ) );  % Pa

Mbase1 = kron(shmidt1,shmidt1);
Mbase2 = kron(shmidt2,shmidt2);
Mbase3 = kron(shmidt3,shmidt3);
M = A*(Mbase1*rss1^(n-1)+Mbase2*rss2^(n-1)+Mbase3*rss3^(n-1));

G1 = A*rss1.*abs(rss1).^(n-1); % s^{-1}
G2 = A*rss2.*abs(rss2).^(n-1); % s^{-1}
G3 = A*rss3.*abs(rss3).^(n-1); % s^{-1}

% calculate the velocity gradient (size 3x3xN)
Gvel = shmidt1*G1+shmidt2*G2+shmidt3*G3; % s^{-1}
GvelT = Gvel';
           
% calculate the strain rate (size 3x3xN)
Gecdot = Gvel/2+Gvel'/2; % s^{-1}

Gmecdot = sqrt(1/2*sum(sum(Gecdot.^2)));

display('************************************')
display('   Stepping a crystal')
display('************************************')

display(shmidt1)
display(shmidt2)
display(shmidt3)
display(rss1)
display(rss2)
display(rss3)
display(Mrss)
display(M)
display(Gvel)
display(GvelT)
display(Gecdot)
display(Gmecdot)

%%

cpp = importdata('comparisonDist_initial.csv',',',1);
cdata = cpp.data;
%# Crystal, C-Axis (x), C-Axis (y), C-Axis (z), Size (m), Disl. dens. (1/m^2), Last recr. time (s), Size at last recr. (m)
cN = cdata(:,2:4)';

% get new angles
HXY = sqrt(cN(1,:).^2+cN(2,:).^2);
c.theta = atan2(HXY,cN(3,:))';
c.phi   = atan2(cN(2,:),cN(1,:))';

c.size = cdata(:,5);
c.dislDens=cdata(:,6);

cdist = c;

SET = load('/home/joseph/Documents/Programs/Thor/trunk/+Test/timing.mat','A');
SET.nelem = 1;
SET.numbcrys = 8000;
SET.stress = stress; 
SET.xcec = [1, 0];
SET.T = temperature; 
SET.glenexp = n;
SET.width = [20 20 20];
SET.tstep = dt*365*24*60*60;
SET.CONN = 'cube8000.mat';
SET.to = zeros(SET.numbcrys,SET.nelem);
SET.ti = zeros(SET.nelem, 1);
SET.poly = true;
SET.migre = true;
SET.run = 1; 
SET.Do = cdist.size;


eigenMask = ones(8000,1);
ii = 1; 

NPOLY = zeros(SET.nelem,size(eigenMask,2));
NMIGRE = zeros(SET.nelem,size(eigenMask,2));

    % calculate velocity gradients and crystal strain rates
    cdist = Thor.Utilities.vec( cdist, SET, ii);        

    % calculate bulk strain rate
    [bedot] = Thor.Utilities.bedot( cdist );

    % grow the crystals
    [cdist, K] = Thor.Utilities.grow(cdist, SET, ii);

    % calculate new dislocation density
    [cdist, ~] = Thor.Utilities.disl(cdist, SET, ii, K);

    % check for migration recrystallization
    if (SET.migre)
        [cdist, SET, NMIGRE(ii,:)] = Thor.Utilities.migre(cdist, SET, ii, eigenMask);
    end

    % check for polygonization
    if (SET.poly)
        [cdist, SET, NPOLY(ii,:)] = Thor.Utilities.poly(cdist, SET, ii, eigenMask);
    end

%     % check crystal orientation bounds
%     cdist = Thor.Utilities.bound(cdist);

    % rotate the crystals from last time steps calculations
    cdist = Thor.Utilities.rotate(cdist, SET, ii );

display('************************************')
display('   Stepping a disribution')
display('************************************')
    
display(NMIGRE)
display(NPOLY)
display(bedot)

%% plot

% sines and cosines so calculation only has to be preformed once
st = sin(cdist.theta); ct = cos(cdist.theta);
sp = sin(cdist.phi);   cp = cos(cdist.phi); sq3 = sqrt(3);

% C-axis orientation
mN = [st.*cp, st.*sp, ct]; % -

fpp = importdata('comparisonDist_final.csv',',',1);
fdata = fpp.data;
%# Crystal, C-Axis (x), C-Axis (y), C-Axis (z), Size (m), Disl. dens. (1/m^2), Last recr. time (s), Size at last recr. (m)
fN = fdata(:,2:4)';

% get new angles
HXY = sqrt(fN(1,:).^2+fN(2,:).^2);
f.theta = atan2(HXY,fN(3,:))';
f.phi   = atan2(fN(2,:),fN(1,:))';

f.size = fdata(:,5);
f.dislDens=fdata(:,6);

cN = cN';
fN = fN';

% make sure all upper hemisphere
cmsk = cN(:,3) < 0;
fmsk = fN(:,3) < 0;
mmsk = mN(:,3) < 0;

cN(cmsk,:) = -cN(cmsk,:);
fN(fmsk,:) = -fN(fmsk,:);
mN(mmsk,:) = -mN(mmsk,:);

MSize = 200;
LWidth = 3;


figure;
scatter3(cN(1:100:end,1),cN(1:100:end,2),cN(1:100:end,3), MSize*2, 'xb','LineWidth',LWidth)
hold on

scatter3(fN(1:100:end,1),fN(1:100:end,2),fN(1:100:end,3),MSize, 'or','LineWidth',LWidth)
hold on

scatter3(mN(1:100:end,1),mN(1:100:end,2),mN(1:100:end,3), MSize,'+g','LineWidth',LWidth)
hold off

axis([-1,1,-1,1,-1,1])
view([0,90])