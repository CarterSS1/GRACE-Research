clear;
clc;
close all;

%% Read In Data
LM = readtable('Mead_Storage.csv');
t_Mead = decyear(  table2array(LM(23654:end,1))  );
V_Mead = table2array(LM(23654:end,2))*1.2335E-6;     % convert to km3 from acre-feet 

%% Get Numbers from GRACE data
GSFC_Raw1895 = readtable('GSFC_mscn1895_1_of_9.csv');
GSFC_Raw1911 = readtable('GSFC_mscn1911_2_of_9.csv');
GSFC_Raw1912 = readtable('GSFC_mscn1912_3_of_9.csv');
GSFC_Raw1854 = readtable('GSFC_mscn1854_4_of_9.csv');
GSFC_Raw1918 = readtable('GSFC_mscn1918_5_of_9.csv');
GSFC_Raw1919 = readtable('GSFC_mscn1919_6_of_9.csv');
GSFC_Raw1859 = readtable('GSFC_mscn1859_7_of_9.csv');
GSFC_Raw1924 = readtable('GSFC_mscn1924_8_of_9.csv');
GSFC_Raw1925 = readtable('GSFC_mscn1925_9_of_9.csv');

GSFCx_Raw = GSFC_Raw1895(13:237,1);
GSFCx = table2array(GSFCx_Raw);

GSFCy_Raw1895 = GSFC_Raw1895(13:237,2);
GSFCy_Raw1911 = GSFC_Raw1911(13:237,2);
GSFCy_Raw1912 = GSFC_Raw1912(13:237,2);
GSFCy_Raw1854 = GSFC_Raw1854(13:237,2);
GSFCy_Raw1918 = GSFC_Raw1918(13:237,2);
GSFCy_Raw1919 = GSFC_Raw1919(13:237,2);
GSFCy_Raw1859 = GSFC_Raw1859(13:237,2);
GSFCy_Raw1924 = GSFC_Raw1924(13:237,2);
GSFCy_Raw1925 = GSFC_Raw1925(13:237,2);

GSFCy1895 = table2array(GSFCy_Raw1895);
GSFCy1911 = table2array(GSFCy_Raw1911);
GSFCy1912 = table2array(GSFCy_Raw1912);
GSFCy1854 = table2array(GSFCy_Raw1854);
GSFCy1918 = table2array(GSFCy_Raw1918);
GSFCy1919 = table2array(GSFCy_Raw1919);
GSFCy1859 = table2array(GSFCy_Raw1859);
GSFCy1924 = table2array(GSFCy_Raw1924);
GSFCy1925 = table2array(GSFCy_Raw1925);

GSFCarea_Raw1895 = GSFC_Raw1895(4,2);
GSFCarea_Raw1911 = GSFC_Raw1911(4,2);
GSFCarea_Raw1912 = GSFC_Raw1912(4,2);
GSFCarea_Raw1854 = GSFC_Raw1854(4,2);
GSFCarea_Raw1918 = GSFC_Raw1918(4,2);
GSFCarea_Raw1919 = GSFC_Raw1919(4,2);
GSFCarea_Raw1859 = GSFC_Raw1859(4,2);
GSFCarea_Raw1924 = GSFC_Raw1924(4,2);
GSFCarea_Raw1925 = GSFC_Raw1925(4,2);

GSFCarea1895 = table2array(GSFCarea_Raw1895);
GSFCarea1911 = table2array(GSFCarea_Raw1911);
GSFCarea1912 = table2array(GSFCarea_Raw1912);
GSFCarea1854 = table2array(GSFCarea_Raw1854);
GSFCarea1918 = table2array(GSFCarea_Raw1918);
GSFCarea1919 = table2array(GSFCarea_Raw1919);
GSFCarea1859 = table2array(GSFCarea_Raw1859);
GSFCarea1924 = table2array(GSFCarea_Raw1924);
GSFCarea1925 = table2array(GSFCarea_Raw1925);

% Net water lost
GSFCy = (GSFCy1895 + GSFCy1911 + GSFCy1912 ...
    + GSFCy1854 + GSFCy1918 + GSFCy1919 + GSFCy1859 + GSFCy1924 + GSFCy1925);
% Mean water lost per mascon
GSFCy = GSFCy / 9;
% Total area
GSFCarea =  GSFCarea1895 + GSFCarea1911 + GSFCarea1912 + GSFCarea1854...
    + GSFCarea1918 + GSFCarea1919 + GSFCarea1859 + GSFCarea1924 + GSFCarea1925;
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
GSFCy = GSFCy .* GSFCarea;

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

