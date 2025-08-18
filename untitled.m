clear;
clc;
close all;

%% Read In Data
LM = readtable('Mead_Storage.csv');
t_Mead = decyear(  table2array(LM(23654:end,1))  );
V_Mead = table2array(LM(23654:end,2))*1.2335E-6;     % convert to km3 from acre-feet 

%% Get Numbers from GRACE data
% Read in data from the 9 GRACE mascons
GSFC_Raw1895 = readtable('GSFC_mscn1895_1_of_9.csv');
GSFC_Raw1911 = readtable('GSFC_mscn1911_2_of_9.csv');
GSFC_Raw1912 = readtable('GSFC_mscn1912_3_of_9.csv');
GSFC_Raw1854 = readtable('GSFC_mscn1854_4_of_9.csv');
GSFC_Raw1918 = readtable('GSFC_mscn1918_5_of_9.csv');
GSFC_Raw1919 = readtable('GSFC_mscn1919_6_of_9.csv');
GSFC_Raw1859 = readtable('GSFC_mscn1859_7_of_9.csv');
GSFC_Raw1924 = readtable('GSFC_mscn1924_8_of_9.csv');
GSFC_Raw1925 = readtable('GSFC_mscn1925_9_of_9.csv');

%Create x axis (Same for all 9 mascons)
GSFCx_Raw = GSFC_Raw1895(13:237,1);
GSFCx = table2array(GSFCx_Raw);

%Measure water loss and gain from each mascon)
GSFCy_Raw1895 = GSFC_Raw1895(13:237,2);
GSFCy_Raw1911 = GSFC_Raw1911(13:237,2);
GSFCy_Raw1912 = GSFC_Raw1912(13:237,2);
GSFCy_Raw1854 = GSFC_Raw1854(13:237,2);
GSFCy_Raw1918 = GSFC_Raw1918(13:237,2);
GSFCy_Raw1919 = GSFC_Raw1919(13:237,2);
GSFCy_Raw1859 = GSFC_Raw1859(13:237,2);
GSFCy_Raw1924 = GSFC_Raw1924(13:237,2);
GSFCy_Raw1925 = GSFC_Raw1925(13:237,2);

% Convert to array
GSFCy1895 = table2array(GSFCy_Raw1895);
GSFCy1911 = table2array(GSFCy_Raw1911);
GSFCy1912 = table2array(GSFCy_Raw1912);
GSFCy1854 = table2array(GSFCy_Raw1854);
GSFCy1918 = table2array(GSFCy_Raw1918);
GSFCy1919 = table2array(GSFCy_Raw1919);
GSFCy1859 = table2array(GSFCy_Raw1859);
GSFCy1924 = table2array(GSFCy_Raw1924);
GSFCy1925 = table2array(GSFCy_Raw1925);

% Pull out area data from each mascon
GSFCarea_Raw1895 = GSFC_Raw1895(4,2);
GSFCarea_Raw1911 = GSFC_Raw1911(4,2);
GSFCarea_Raw1912 = GSFC_Raw1912(4,2);
GSFCarea_Raw1854 = GSFC_Raw1854(4,2);
GSFCarea_Raw1918 = GSFC_Raw1918(4,2);
GSFCarea_Raw1919 = GSFC_Raw1919(4,2);
GSFCarea_Raw1859 = GSFC_Raw1859(4,2);
GSFCarea_Raw1924 = GSFC_Raw1924(4,2);
GSFCarea_Raw1925 = GSFC_Raw1925(4,2);

% Convert to array
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
% Subtract the starting amount of water from volume data to get water
% equivalent height
% Convert GRACE data from centimeters to Kilometers ^ 3 using the surface area of
% the mascon
figure();
% V_Mead = V_Mead - 24.2633;
V_Mead = V_Mead - mean(V_Mead);

plot(t_Mead,V_Mead)
hold on;
GSFCy = GSFCy ./ 100000;
GSFCy = GSFCy .* GSFCarea;
GSFCy = GSFCy - mean(GSFCy);

plot (GSFCx,GSFCy)

title('Water Equivalent Height (km^3) vs years')
legend('Water Equivalent Height from Volume data','Water Equivalent Height from GRACE Data','Location','northwest')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);

% %% Plotting Lines of Best Fit
% % Create the lines of best fit along with error bars that represent 95%
% % accuracy.
% figure();
% hold on;
% [mb1,S1] = polyfit(t_Mead, V_Mead,1);
% y1 = mb1(1) .* t_Mead + mb1(2);
% plot(t_Mead,y1);
% 
% [mb2,S2] = polyfit(GSFCx, GSFCy,1);
% y2 = mb2(1) .* GSFCx + mb2(2);
% plot(GSFCx,y2);
% 
% [extrap_fit1,delta1] = polyval(mb1,t_Mead,S1);
% errorMinus1 = extrap_fit1- 2 * delta1;
% errorplus1 = extrap_fit1+ 2 * delta1;
% plot(t_Mead,errorMinus1,'--');
% plot(t_Mead,errorplus1,'--');
% 
% [extrap_fit2,delta2] = polyval(mb2,GSFCx,S1);
% errorMinus2 = extrap_fit2 - 2 * delta2;
% errorplus2 = extrap_fit2 + 2 * delta2;
% plot(GSFCx,errorMinus2,'--');
% plot(GSFCx,errorplus2,'--');
% 
% 
% title('Water Equivalent Height (km^3) vs years')
% legend('Line of Best Fit from Volume data',' Line of Best Fit from GRACE Data','Volume Data - 2 * \sigma_x','Volume Data + 2 * \sigma_x','GRACE Data - 2 * \sigma_x','GRACE data + 2 * \sigma_x','Location','northwest')
% xlabel('time in years since 2002')
% ylabel('Water Equivalent Height (km^3)')
% xlim([2002 2024]);
% 
% 

%% Read GLDAS
%% Lat1= 36.03;  Lon1= -114.44;

% Re = 6378136.3E-3;
% lat1 = [36.5-0.5     36.5+0.5       36.5+0.5      36.5-0.5];
% lon1 = [-114.5-0.5   -114.5-0.5     -114.5+0.5    -114.5+0.5];
% % figure;plot(lon1,lat1,'*-')
% a1 = 4*pi*Re^2*areaint(lat1,lon1)
% 
% wgs84 = wgs84Ellipsoid("km");
% a2 = areaint(lat1,lon1,wgs84)
% 
% a3 = (Re*1*pi/180)*(cosd(36.5)*Re*1*pi/180)

%% 
files = dir('GLDAS-Noah-TWS/*.nc');
name1=fullfile(files(1).folder,files(1).name);
% ncdisp(name1);
lat_GLDAS = ncread(name1,'lat');
lon_GLDAS = ncread(name1,'lon');
[lonM_GLDAS, latM_GLDAS] = meshgrid(lon_GLDAS, lat_GLDAS);

tws1_GLAD = ncread(name1,'TWS_monthly');

%% 9 GSFC 1-deg mascons covering Lake Mead and its surroundings
% --- 1: center
% Mascon #: 1918
Lat1= 36.03;
Lon1= -114.44;




% -- 2: west
% Mascon #: 1854
Lat2= 36.15;
Lon2= -115.45;

% -- 3: northwest
% Mascon #: 1895
Lat3= 37.02;
Lon3= -115.92;

% -- 4: north
% Mascon #: 1911
Lat4= 36.93;
Lon4= -114.79;

% --- 5: northeast
% Mascon #: 1912
Lat5= 36.97;
Lon5= -113.37;

% --- 6: east
% Mascon #: 1919
Lat6= 36.11;
Lon6= -113.12;

% --- 7: southeast
% Mascon #: 1925
Lat7= 35.07;
Lon7= -113.14;

% --- 8: south
% Mascon #: 1924
Lat8= 34.89;
Lon8= -114.20;

% --- 9: southwest
% Mascon #: 1859
Lat9= 34.82;
Lon9= -115.49;

lat_Int = [Lat1; Lat2; Lat3; Lat4; Lat5; Lat6; Lat7; Lat8; Lat9];
lon_Int = [Lon1; Lon2; Lon3; Lon4; Lon5; Lon6; Lon7; Lon8; Lon9];
mascon_ID = [1918; 1854; 1895; 1911; 1912; 1919; 1925; 1924; 1859]; 

for i = 1:length(mascon_ID)
    [minA, lon_idx(i,1)] = min(abs(lon_GLDAS-lon_Int(i)));
    [minB, lat_idx(i,1)] = min(abs(lat_GLDAS-lat_Int(i)));
    % tws1_GLAD(lon_idx(i,1), lat_idx(i,1))

    % lonM_GLDAS0 = lonM_GLDAS;
    % latM_GLDAS0 = latM_GLDAS;
    % nanGLDAS = isnan(tws1_GLAD);
    % tws1_GLAD(nanGLDAS==1)=[];
    % lonM_GLDAS0(nanGLDAS==1)=[];
    % latM_GLDAS0(nanGLDAS==1)=[];
    % TWS_atGR= griddata(lonM_GLDAS0', latM_GLDAS0',tws1_GLAD',lon_Int(i,1),lat_Int(i,1),'linear');


end


%% 



for k = 1:length(files)

    lonM_GLDAS0 = lonM_GLDAS;
    latM_GLDAS0 = latM_GLDAS;

    namek=fullfile(files(k).folder,files(k).name);
    t0_days = ncread(namek,'time');
    t_GLDAS(k,1) = decyear( datetime(2002,01,01) + t0_days );

    twsi_GLAD = ncread(namek,'TWS_monthly');


    for ll=1:9
        TWS_atGR(1,ll) = twsi_GLAD(lon_idx(ll,1), lat_idx(ll,1));
    end

    % nanGLDAS = isnan(twsi_GLAD);
    % twsi_GLAD(nanGLDAS==1)=[];
    % lonM_GLDAS0(nanGLDAS==1)=[];
    % latM_GLDAS0(nanGLDAS==1)=[];
    % TWS_atGR= griddata(lonM_GLDAS0, latM_GLDAS0,twsi_GLAD,lon_Int',lat_Int','nearest');

    % F = scatteredInterpolant(lonM_GLDAS0', latM_GLDAS0',twsi_GLAD'); 
    % TWS_atGR = F(lon_Int,lat_Int)';

    GLDAS_9mas(k,1:9) = TWS_atGR * 0.1;  % convert to cm


end


avg_GLDAS = mean(GLDAS_9mas')';   % EWH in cm
avg_GLDAS_vol = (avg_GLDAS/1E5) * GSFCarea;  % volume in km^3

avg_GLDAS_vol = avg_GLDAS_vol - mean(avg_GLDAS_vol);

figure(); 

plot(t_GLDAS, avg_GLDAS_vol,'*-'); grid on; hold on;
plot (GSFCx,GSFCy)
title('Measured Water (km^3) vs years')
legend('Water Moisture',' Water equivalent height GRACE Data','Location','northwest')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);
save GLDAS_9mascon t_GLDAS  GLDAS_9mas;


% plottedGLDAS = mean(GLDAS_9mas')';
% GLDAS_soil = [t_GLDAS  GLDAS_9mas];

plottedGLDAS = avg_GLDAS_vol;
GLDAS_soil = [t_GLDAS  avg_GLDAS_vol];

yGLDASInterpolated = interp1(t_GLDAS, plottedGLDAS, GSFCx);

FinalWater = (GSFCy - yGLDASInterpolated);

% GLDAS_soil = [0 mascon_ID'; 
%               0 lat_Int'; 
%               0 lon_Int'; 
%               GLDAS_soil];
% FinalWater = (GSFCy - yGLDASInterpolated)-5.9921;

figure();
hold on
plot(t_Mead,V_Mead,'k');
plot(GSFCx,GSFCy,'b');
plot(GSFCx,FinalWater,'r');


title('Water Equivalent Height (km^3) vs years')
% legend('Line of Best Fit from Volume data',' Line of Best Fit from GRACE Data','Volume Data - 2 * \sigma_x','Volume Data + 2 * \sigma_x','GRACE Data - 2 * \sigma_x','GRACE data + 2 * \sigma_x','Location','northwest')
legend('in situ','GRACE','GRACE - GLDAS')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);
% write data 
% save GLDAS_soil  GLDAS_soil;

% writematrix(GLDAS_soil, 'GLDAS_soil.xls')


%% removing teh groundwater effect 


% t_Mead V_Mead
% GSFCx GSFCy
%  GSFCx FinalWater : GRAACE - soil moisture

load GW_WGHM.mat;
%    t_GW   GW_Vol
GW_Vol = GW_Vol - mean(GW_Vol);

% t_GR = GSFCx

x = GSFCx;
dt = datetime(floor(x), 1, 1) + years(x-floor(x));
[Y_GR, M_GR, D_GR] = ymd(dt);
t_YMD_GR = [Y_GR, M_GR, D_GR];


% 
t_YMD_GR(160, 2) = 3;
t_YMD_GR(166, 2) = 10;
%


% GSFCx GSFCy  FinalWater
Lia = ismember(t_GW(:,1:2), t_YMD_GR(:,1:2)  ,'rows');
% sum(Lia)
idx_common = find(Lia == 1);

t_GW_2 = t_GW(idx_common,:);
GW_Vol_2 = GW_Vol(idx_common);

t_GW_2_decY = decyear(t_GW_2);


%
Lia_3 = ismember(t_YMD_GR(:,1:2), t_GW(:,1:2)  ,'rows');
% sum(Lia_3)
idx_common_3 = find(Lia_3 == 1);

GSFCx_2 = GSFCx(idx_common,1);
GSFCy_2 = GSFCy(idx_common,1);
FinalWater_2 = FinalWater(idx_common,1);

% remove GW effect 
FinalWater_3 = FinalWater_2 - GW_Vol_2;

% 
figure();
hold on
plot(t_Mead,V_Mead,'k','linewidth',3);
plot(t_GLDAS, avg_GLDAS_vol,'g','linewidth',2)
plot(t_GW_2_decY, GW_Vol_2,'m','linewidth',2);
plot(GSFCx,GSFCy,'b');
% plot(t_GW_2_decY,FinalWater_2,'c','linewidth',2);
% plot(t_GW_2_decY,FinalWater_3,'r','linewidth',2);

title('Water Equivalent Height (km^3) vs years')
% legend('Line of Best Fit from Volume data',' Line of Best Fit from GRACE Data','Volume Data - 2 * \sigma_x','Volume Data + 2 * \sigma_x','GRACE Data - 2 * \sigma_x','GRACE data + 2 * \sigma_x','Location','northwest')
legend('in situ','soil moisture','groundwater','GRACE')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);
grid on;



%
figure();
hold on
plot(t_Mead,V_Mead,'k','linewidth',3);
% plot(t_GLDAS, avg_GLDAS_vol,'g')
% plot(t_GW_2_decY, GW_Vol_2,'m');
plot(GSFCx,GSFCy,'b','linewidth',2);
plot(t_GW_2_decY,FinalWater_2,'c','linewidth',2);
plot(t_GW_2_decY,FinalWater_3,'r','linewidth',2);


title('Water Equivalent Height (km^3) vs years')
% legend('Line of Best Fit from Volume data',' Line of Best Fit from GRACE Data','Volume Data - 2 * \sigma_x','Volume Data + 2 * \sigma_x','GRACE Data - 2 * \sigma_x','GRACE data + 2 * \sigma_x','Location','northwest')
legend('in situ','GRACE','GRACE - soil','GRACE - soil - GW')
xlabel('time in years since 2002')
ylabel('Water Equivalent Height (km^3)')
xlim([2002 2024]);
grid on;