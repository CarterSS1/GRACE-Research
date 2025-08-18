clear;
clc;
close all;

%% Read In Data
LM = readtable('Mead_Storage.csv');
t_Mead = decyear(  table2array(LM(23654:end,1))  );
V_Mead = table2array(LM(23654:end,2))*1.2335E-6;     % convert to km3 from acre-feet 

%% Get Numbers from GRACE data
GSFC_Raw = readtable('GSFC_mscn1918.csv');
GSFCx_Raw = GSFC_Raw(13:237,1);
GSFCx = table2array(GSFCx_Raw);
GSFCy_Raw = GSFC_Raw(13:237,2);
GSFCy = table2array(GSFCy_Raw);

%% Create Plot for Volume Data of the Water
%%
% Subtract the starting amount of water from volume data to get water
% equivalent height
%%
% Convert GRACE data from centimeters to Kilometers ^ 3 using the surface area of
% the mascon
figure();
V_Mead = V_Mead - 24.2633;
plot(t_Mead,V_Mead)
hold on;
GSFCy = GSFCy ./ 100000;
GSFCy = GSFCy .* 1232;
GSFCy = GSFCy .* 9;
plot (GSFCx,GSFCy)

title('Water Equivalent Height (km^3) vs years')
legend('Water Equivalent Height from Volume data','Water Equivalent Height from GRACE Data','Location','northwest')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);

%% Plotting Lines of Best Fit
% Create the lines of best fit along with error bars that represent 95%
% accuracy.
figure();
hold on;
[mb1,S1] = polyfit(t_Mead, V_Mead,1);
y1 = mb1(1) .* t_Mead + mb1(2);
plot(t_Mead,y1);

[mb2,S2] = polyfit(GSFCx, GSFCy,1);
y2 = mb2(1) .* GSFCx + mb2(2);
plot(GSFCx,y2);

[extrap_fit1,delta1] = polyval(mb1,t_Mead,S1);
errorMinus1 = extrap_fit1- 2 * delta1;
errorplus1 = extrap_fit1+ 2 * delta1;
plot(t_Mead,errorMinus1,'--');
plot(t_Mead,errorplus1,'--');

[extrap_fit2,delta2] = polyval(mb2,GSFCx,S1);
errorMinus2 = extrap_fit2 - 2 * delta2;
errorplus2 = extrap_fit2 + 2 * delta2;
plot(GSFCx,errorMinus2,'--');
plot(GSFCx,errorplus2,'--');


title('Water Equivalent Height (km^3) vs years')
legend('Line of Best Fit from Volume data',' Line of Best Fit from GRACE Data','Volume Data - 2 * \sigma_x','Volume Data + 2 * \sigma_x','GRACE Data - 2 * \sigma_x','GRACE data + 2 * \sigma_x','Location','northwest')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);

%% Trying to Graph the difference from data point to data point by interpolating but I could not figure it out

% yGSFCInterpolated = interp1(GSFCx, GSFCy, t_Mead);
% yGSFCInterpolated(1:50) = 0;
% badrows=isnan(yGSFCInterpolated);
% 
% N = length(yGSFCInterpolated);
% badrows = double(badrows);
% 
% for i = 1:N
%     if badrows(i) == 1
%         yGSFCInterpolated(i) = [0];
%     end
% end
% yGSFCInterpolated(isnan(yGSFCInterpolated)) = 0;
% 
% GSFCx = GSFCx - 2002;
% t_Mead = t_Mead - 2002;
% 
% n1 = length(V_Mead);
% yGSFCInterpolated = interp1(GSFCx, GSFCy, t_Mead);
% V_MeadDiff = zeros(n1,1);
% for i = 1:(n1-1)
%     V_MeadDiff(i) = V_Mead(i+1) - V_Mead(i);
% end
% 
% figure();
% hold on;
% plot(t_Mead,yGSFCInterpolated);
% plot(t_Mead,V_MeadDiff.*10);
% 
% V_Mead = V_Mead - 24.2633;
% GSFCx = GSFCx - 2002;
% t_Mead = t_Mead - 2002;
% GSFCy = GSFCy ./ 100000;
% GSFCy = GSFCy .* 12321;
% n1 = length(V_Mead);
% yGSFCInterpolated = interp1(GSFCx, GSFCy, t_Mead);
% V_MeadDiff = zeros(n1,1);
% for i = 1:(n1-1)
%     V_MeadDiff(i) = V_Mead(i+1) - V_Mead(i);
% end
% 
% figure();
% hold on;
% plot(t_Mead,yGSFCInterpolated);
% plot(t_Mead,V_MeadDiff.*10);