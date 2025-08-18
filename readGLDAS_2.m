%% Read GLDAS 

clc
clear
close all


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


figure; plot(t_GLDAS, mean(GLDAS_9mas')','*-'); grid on;

save GLDAS_9mascon t_GLDAS  GLDAS_9mas;

GLDAS_soil = [t_GLDAS  GLDAS_9mas];

GLDAS_soil = [0 mascon_ID'; 
              0 lat_Int'; 
              0 lon_Int'; 
              GLDAS_soil];

% write data 
save GLDAS_soil  GLDAS_soil;

writematrix(GLDAS_soil, 'GLDAS_soil.xls')

