%% read GW data from WGHM model

clc
clear
close all

%% 

dir_data = 'data/watergap22e_gswp3-w5e5_groundwstor_histsoc_monthly_1901_2019.nc';

% get the info
ncdisp(dir_data);

%%  extract data 
% latitude 
lat_WGHM = ncread(dir_data,'lat');
lat_WGHM = double(lat_WGHM);

% long 
lon_WGHM = ncread(dir_data,'lon');
lon_WGHM = double(lon_WGHM);

% time 
t_WGHM = ncread(dir_data,'time');
t_WGHM = double(t_WGHM);


date_WGHM = datetime('1901-01-01') + calmonths(t_WGHM);

YMD_WGHM = datevec(date_WGHM);
YMD_WGHM(:,4:6) = [];

% GW 
gws_WGHM = ncread(dir_data,'groundwstor');  % unit: mm
gws_WGHM = double(gws_WGHM);

gws_WGHM = pagetranspose(gws_WGHM);


%% plot global grid of a particular month
year_i = 2009;
month_i = 11;

idx = find(YMD_WGHM(:,1)==year_i & YMD_WGHM(:,2)== month_i);

gws_month = gws_WGHM(:,:, idx);

% figure; imagesc(gws_1901_jan)
% colorbar; colormap(jet); 
% clim([-300   +300])
% title('Jan. 2019')

% meshgrid 
[lonMESH, latMESH] = meshgrid(lon_WGHM, lat_WGHM);

% 
load coastlines; 

% 
figure;
worldmap World;
geoshow(latMESH,lonMESH,gws_month,'DisplayType','texturemap')
hold on; 
geoshow(coastlat,coastlon,'DisplayType','line','color','k','linewidth',1.5)
colorbar; colormap(jet); 
clim([-300   +300])
title([year_i  month_i])



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

lat_L = lat_Int - 0.5; 
lat_U = lat_Int + 0.5; 

lon_L = lon_Int - 0.6;
lon_R = lon_Int + 0.6;

lat_vec = [lat_L  ; lat_U];
lon_vec = [lon_L  ; lon_R];

% extract boundaries
min_lat = min(lat_vec);
max_lat = max(lat_vec);

min_lon = min(lon_vec);
max_lon = max(lon_vec);

%% extract gw storage from WGHM model in the study area

% latitude bound
idx_lat = find(lat_WGHM>=min_lat & lat_WGHM<=max_lat);
% longitude bound
idx_lon = find(lon_WGHM>=min_lon & lon_WGHM<=max_lon);
% time (GRACE period)
t_begin = [2002  04];  % 2002.296
idx_t_begin = find(YMD_WGHM(:,1)==2002 & YMD_WGHM(:,2)== 4);
t_end = [2019  12];   % 
idx_t_end = find(YMD_WGHM(:,1)==2019 & YMD_WGHM(:,2)== 12);

% time vector 
tvec_3 = ( datetime(YMD_WGHM(idx_t_begin,:)) :calmonths(1): datetime(YMD_WGHM(idx_t_end,:)) )';
tvec_4 = datevec(tvec_3);
tvec_4(:,4:6) = [];
decyr_tvec_4 = decyear(tvec_4);


gws_WGHM_LakeMead = gws_WGHM(idx_lat, idx_lon, idx_t_begin:idx_t_end);

for i = 1:length(decyr_tvec_4)

  ts_gws_LM(i,1) = mean2(gws_WGHM_LakeMead(:,:,i));

end


ts_gws_LM_res = ts_gws_LM - mean(ts_gws_LM);

figure; plot(decyr_tvec_4, ts_gws_LM_res/10)
ylabel('GWS [cm]')
grid on;



GSFCarea =     111477.72; % km^2
ts_gws_LM_res_vol = (ts_gws_LM_res/10)*1E-5*GSFCarea;


figure; plot(decyr_tvec_4, ts_gws_LM_res_vol)
ylabel('GWS [km^3]')
grid on;

