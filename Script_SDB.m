%%% IMPORTANT !
%%% LOAD ALL FILES IN THE "CURRENT FOLDER" WHERE MATLAB OPERATES
%%% MIND NOT TO OVERWRITE FILES IN THE RESULTS FOLDER, CHANGE THE WORKSPACE
%%% NAME OR THE FILENAMES BEFORE NEW RUNS

%%% CLEAR WORKSPACE AND COMMAND WINDOW
close all
clear
clc

%%% READ THE PREPROCESSED SATELLITE IMAGE RASTER FILE WITH THE RESPECTIVE "R" REFERENCE SYSTEM FILE
%%% -> RENAME WITH THE PROPER FILENAME.tif

[A, R] = readgeoraster([pwd '/Data/S2A_MSIL2A_20200319T101021_N0214_R022_T32TQR_20200319T130518_resampled_Subset_BandMaths_Deglint.tif']);

%%% UNCOMMENT TO CREATE AND SHOW THE MERGED RGB IMAGE 
% Image = cat(3, A(:,:,3), A(:,:,2), A(:,:,1));
% RGBImage = double(Image);
% imshow(RGBImage)
% hold off
% close all
%%% RGB BANDS TO MATRICES CORRECTLY ORIENTED
Red= flip(A(:,:,3));
Green= flip(A(:,:,2));
Blue= flip(A(:,:,1));


%%% READ THE GROUND TRUTH VECTOR FILE 
%%% -> RENAME WITH THE PROPER DEPTH GROUND TRUTH POINTS FILENAME.shp 

V = readgeotable([pwd '/Data/P_ICESAT_2020_IN.shp']);

%%% EPSG CODE RETRIEVAL 
CRSstr= R.ProjectedCRS.Name;
zone=extractBetween(CRSstr,19,20);
emi=extractAfter(CRSstr,20);
if emi=="S"
	rootepsg="327" ; 
else
	
	rootepsg="326" ;
end
epsg=strcat(rootepsg,zone);

Lat = V.Shape.Latitude(:,1);
Lon = V.Shape.Longitude(:,1);
p2 = projcrs(str2num(epsg));
[xx,yy] = projfwd(p2,Lat,Lon);

%%% RETRIEVE GROUND TRUTH COORDINATES FROM "T" AND EXTRACT RGB VALUES IN
%%% THAT POINTS

ind=[(xx)-(R.XWorldLimits(1)), (yy)-R.YWorldLimits(1)];
Ind=round(ind./10);    % 10m IS THE SPATIAL RESOLUTION OF THE RASTER

X = Ind(:,1);
Y = Ind(:,2);

%%% UNCOMMENT TO SHOW MAP WITH WITH S-2 AND GROUND TRUTH POINTS AS CONTROL
% pcolor(Red)
% hold on
% shading flat
% hold on
% plot(X, Y,  'r*')
% hold off
% close all
%%% CLOSE PLOT 

%%% CREATE AND FILL MATRICES OF THE PROPER SIZE WITH RGB VALUES
aa=size(X);
Mat_R=zeros(aa(1),1);
Mat_G=zeros(aa(1),1);
Mat_B=zeros(aa(1),1);
s = aa(1);

for c = 1:s
   b=  X(c);
   a=  Y(c);
   Mat_R(c)=Red(a,b);
end

for c = 1:s
   b=  X(c);
   a=  Y(c);
   Mat_G(c)=Green(a,b);
end

for c = 1:s
   b=  X(c);
   a=  Y(c);
   Mat_B(c)=Blue(a,b);
end

Depth = table2array(V(:,2));


%%% CALCULATE RATIOS VALUES (ACCORDING TO STUMPF 2003) IN ARRAYS

Ratio_BR = log(1000*Mat_B)./log(1000*Mat_R);
Ratio_BG = log(1000*Mat_B)./log(1000*Mat_G);

%%% CALCULATE RATIOS VALUES IN  MAPS
Map_ratio_red = log(double(Blue)*1000)./ log(double(Red)*1000);
Map_ratio_green = log(double(Blue)*1000)./ log(double(Green)*1000);

%%% FIND AND REMOVE NAN
N = isnan(Ratio_BR);
Depth = Depth(~N);
Lat = Lat(~N);
Lon= Lon(~N);
Ratio_BR = Ratio_BR(~N);
Ratio_BG= Ratio_BG(~N);

N = isnan(Ratio_BG);
Depth = Depth(~N);
Lat = Lat(~N);
Lon= Lon(~N);
Ratio_BR = Ratio_BR(~N);
Ratio_BG= Ratio_BG(~N);

Sampled_Points = [Lat,Lon,Depth,Ratio_BR, Ratio_BG];

%%% TRAINING AND TESTING POINTS SETS CREATION
[m,n] = size(Sampled_Points) ;
P = 0.4; %%% CHANGE THE SPLIT VALUES ACCORDING TO PREFERENCE
idx = randperm(m)  ;
Training = Sampled_Points(idx(1:round(P*m)),:) ; 
Testing = Sampled_Points(idx(round(P*m)+1:end),:) ;

Depth_Cal = Training(:,3);
Depth_Val = Testing(:,3);

Ratio_BR_Cal = Training(:,4);
Ratio_BR_Val = Testing(:,4);

Ratio_BG_Cal = Training(:,5);
Ratio_BG_Val = Testing(:,5);

%%% CALIBRATION - RED BAND RATIO 

x = Depth_Cal;
y = Ratio_BR_Cal;

%%% LINEAR REGRESSION
p1 = polyfit(x,y,1);
f1=polyval(p1,x);

mdl = fitlm(x,y);
R_squared = mdl.Rsquared.Adjusted;
N = length(x);

%%% LINEAR REGRESSION
m1 = 1/p1(1);
m0 = p1(2)/p1(1);

%%% OUTLIERS REMOVAL AND REGRESSION PLOT
f = fittype('f1*x');
fit1 = fit(x,y,'poly1');
hold on
fdata = feval(fit1,x); 
I = abs(f1 - y) > 3*std(y); 
outliers = excludedata(x,y,'indices',I);
Y= y(~outliers);
X= x(~outliers);
M=max(X);
fit2 = fit(X,Y,'poly1');
hold on
plot(x,y,'w*')
plot(fit2,'r-', x,y,'k.',outliers,'m*')
hold on
plot(x, f1+3*std(y),   'g:')
hold on
plot(x,  f1-3*std(y),  'g:')

hold on
Min_y = 1-min(y)- 3*std(y);
Max_y = 1+max(y)+ 3*std(y);
Min_x = min(x)-0.2;
Max_x = max(x)+0.2;

grid on
axis([Min_x Max_x  Min_y Max_y])
title('Calibration'),
hold on
ylabel('BLUE/RED')
hold on
xlabel('Depth (m)')

p1_clean = polyfit(X,Y,1);
f1_clean=polyval(p1_clean,X);
mdl_clean = fitlm(X,Y);
R_squared_clean = mdl_clean.Rsquared.Adjusted;
N = length(X);

m1_clean = 1/p1_clean(1);
m0_clean = p1_clean(2)/p1_clean(1);

theString1 = sprintf('m1 = %.3f, m0 = %.3f ', m1_clean, m0_clean);
theString2= sprintf('R2 = %0.2f\n',R_squared);
theString3= sprintf(' N = %0.2f\n',N);
theString4= sprintf('R2 clean = %0.2f\n',R_squared_clean);

legend(theString1, theString3, theString4, Location='northwest'); 
hold on
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
hold on
saveas(gcf,[pwd '/Data/Results/CAL_BR.png']) % -> SAVE THE PLOT

CAL_STAS_BR = [theString3, theString1, theString2, theString4]
hold off
close all
writematrix(CAL_STAS_BR, [pwd '/Data/Results/CAL_STATS_BR.csv'] )% SAVE THE STATS


%%% VALIDATION - RED BAND RATIO 

SDB_BR = (Ratio_BR_Val)*m1_clean - m0_clean; %%% ALGORITHM APPLICATION FOR VALIDATION (ACCORDING TO STUMPF 2003)

L = [1 2 3 4 5 6 7 8 9 10]; 
x = Depth_Val;
y = SDB_BR;

p1 = polyfit(x,y,1);
f1=polyval(p1,x);

mdl = fitlm(x,y);

R_squared = mdl.Rsquared.Adjusted;
N = length(x);
BIAS = x - y;
BIAS_STD = std(BIAS);
BIASP = BIAS > 10;
BIAS = BIAS(~BIASP);
BIASN = BIAS < -10;
BIAS = BIAS(~BIASN);

BIAS_BR= BIAS;

%%% VALIDATION STATS
RMSE = sqrt(mean((BIAS).^2));
MAE = mean(abs(BIAS));
BIAS_AV = mean(BIAS);
BIAS_STD = std(BIAS);

res1=polyval(p1, x)-y;
rmse = sprintf('RMSE =   %0.2f\n', RMSE);
r2= sprintf('R2 = %0.2f\n',R_squared);
n= sprintf(' N = %0.2f\n',N);
biasav = sprintf('BIAS AV =  %0.2f\n',BIAS_AV);
mae = sprintf('MAE = %.3f',MAE);
biastd = sprintf('BIAS STD = %.3f',BIAS_STD );

%%% VALIDATION PLOTS
plot(x,y,'k.')
hold on
plot(x,f1,'w-')
hold on
plot(x,f1,'w-')
hold on
plot(x,f1,'w-')
hold on
plot(x,f1,'w-')
hold on
plot(x, f1, 'k-')
hold on
plot(L,'r--')
hold on
grid on
hold on
xlabel('Depth (m)')
ylabel('SDB (m)')
hold on 
legend(n, rmse, mae, biastd, biasav, Location='northwest');
hold on
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
hold on
saveas(gcf,[pwd '/Data/Results/VAL_BR.png'])
hold off
VAL_STATS_BR = [n rmse mae biastd biasav r2]
close all
histogram(BIAS)
hold on
xline(0, 'r--', ... 
      'LineWidth', 3, ...
    'Interpreter', 'latex', ...
    'LabelOrientation', 'horizontal')
hold on
saveas(gcf, [pwd '/Data/Results/BIAS_BR.png'])
hold off
close all

hold on
writematrix( VAL_STATS_BR, [pwd '/Data/Results/VAL_STATS_BR.csv'] )

Map_SDB_BR = -flip((Map_ratio_red)*m1_clean - m0_clean); %ALGORITHM APPLICATION (ACCORDING TO STUMPF 2003), THE MINUS TO SHOW DEPTH NEGATIVE VALUES
imshow(Map_SDB_BR, [-10 0]) %-> CHANGE [MIN MAX] VALUES ACCORDING TO AOI FEATURES
hold on
colormap turbo
hold on
colorbar eastoutside

%%% SAVE MAP IN PNG FORMAT
saveas(gcf,[pwd '/Data/Results/MAP_SDB_BR.png'])

%%% SAVE MAP IN TIFF FORMAT WITH THE PROPER EPSG CODE FOR THE GEOGRAPHICAL
% ZONE OF THE AOI
Epsg= strcat('EPSG:',strcat(rootepsg,zone));
geotiffwrite( [pwd '/Data/Results/Map_SDB_BR'], Map_SDB_BR, R, CoordRefSysCode= Epsg);
hold off
close all



%%% CALIBRATION - GREEN BAND RATIO

x = Depth_Cal;
y = Ratio_BG_Cal;

%%% LINEAR REGRESSION
p1 = polyfit(x,y,1);
f1=polyval(p1,x);

mdl = fitlm(x,y);
R_squared = mdl.Rsquared.Adjusted;
N = length(x);

m1 = 1/p1(1);
m0 = p1(2)/p1(1);

%%% OUTLIERS REMOVAL AND REGRESSION PLOT
f = fittype('f1*x');
fit1 = fit(x,y,'poly1');
hold on
fdata = feval(fit1,x); 
I = abs(f1 - y) > 3*std(y); 
outliers = excludedata(x,y,'indices',I);
Y= y(~outliers);
X= x(~outliers);
M=max(X);
fit2 = fit(X,Y,'poly1');
hold on
plot(x,y,'w*')
plot(fit2,'r-', x,y,'k.',outliers,'m*')
hold on
plot(x, f1+3*std(y),   'g:')
hold on
plot(x,  f1-3*std(y),  'g:')

hold on
Min_y = 1-min(y)- 3*std(y);
Max_y = 1+max(y)+ 3*std(y);
Min_x = min(x)-0.2;
Max_x = max(x)+0.2;

grid on
axis([Min_x Max_x  Min_y Max_y])
title('Calibration'),
hold on
ylabel('BLUE/GREEN')
hold on
xlabel('Depth (m)')

p1_clean = polyfit(X,Y,1);
f1_clean=polyval(p1_clean,X);
mdl_clean = fitlm(X,Y);
R_squared_clean = mdl_clean.Rsquared.Adjusted;
N = length(X);

m1_clean = 1/p1_clean(1);
m0_clean = p1_clean(2)/p1_clean(1);

theString1 = sprintf('m1 = %.3f, m0 = %.3f ', m1_clean, m0_clean);
theString2= sprintf('R2 = %0.2f\n',R_squared);
theString3= sprintf(' N = %0.2f\n',N);
theString4= sprintf('R2 clean = %0.2f\n',R_squared_clean);

legend(theString1, theString3, theString4, Location='northwest'); 
hold on
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
hold on
hold on
saveas(gcf,[pwd '/Data/Results/CAL_BG.png']) 

CAL_STAS_BG = [theString3, theString1, theString2, theString4]
hold off
close all

writematrix(CAL_STAS_BG, [pwd '/Data/Results/CAL_STATS_BG.csv'] ) 

%%% VALIDATION - GREEN BAND RATIO 

SDB_BG = (Ratio_BG_Val)*m1_clean - m0_clean; % %%% ALGORITHM APPLICATION FOR VALIDATION (ACCORDING TO STUMPF 2003)

L = [1 2 3 4 5 6 7 8 9 10];
x = Depth_Val;
y = SDB_BG;

p1 = polyfit(x,y,1);
f1=polyval(p1,x);

mdl = fitlm(x,y);

R_squared = mdl.Rsquared.Adjusted;
N = length(x);
BIAS = x - y;
BIAS_STD = std(BIAS);
BIASP = BIAS > 10;
BIAS = BIAS(~BIASP);
BIASN = BIAS < -10;
BIAS = BIAS(~BIASN);

BIAS_BG= BIAS;

%%% VALIDATION STATS
RMSE = sqrt(mean((BIAS).^2));
MAE = mean(abs(BIAS));
BIAS_AV = mean(BIAS);
BIAS_STD = std(BIAS);

res1=polyval(p1, x)-y;
rmse = sprintf('RMSE =   %0.2f\n', RMSE);
r2= sprintf('R2 = %0.2f\n',R_squared);
n= sprintf(' N = %0.2f\n',N);
biasav = sprintf('BIAS AV =  %0.2f\n',BIAS_AV);
mae = sprintf('MAE = %.3f',MAE);
biastd = sprintf('BIAS STD = %.3f',BIAS_STD );

%%% VALIDATION PLOTS
plot(x,y,'k.')
hold on
plot(x,f1,'w-')
hold on
plot(x,f1,'w-')
hold on
plot(x,f1,'w-')
hold on
plot(x,f1,'w-')
hold on
plot(x, f1, 'k-')
hold on
plot(L,'r--')
hold on
grid on
hold on
xlabel('Depth (m)')
ylabel('SDB (m)')
hold on 
legend(n, rmse, mae, biastd, biasav, Location='northwest');
hold on
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
hold on
saveas(gcf,[pwd '/Data/Results/VAL_BG.png'])
hold off
close all

histogram(BIAS)
hold on
xline(0, 'r--', ... 
      'LineWidth', 3, ...
    'Interpreter', 'latex', ...
    'LabelOrientation', 'horizontal')
hold on
saveas(gcf, [pwd '/Data/Results/BIAS_BG.png'])
hold off
close all

VAL_STATS_BG = [n rmse mae biastd biasav r2]
hold on
writematrix( VAL_STATS_BG, [pwd '/Data/Results/VAL_STATS_BG.csv'] )

Map_SDB_BG = -flip((Map_ratio_green)*m1_clean - m0_clean); % ALGORITHM APPLICATION (ACCORDING TO STUMPF 2003), THE MINUS TO SHOW DEPTH NEGATIVE VALUES
imshow(Map_SDB_BG, [-10 0]) %-> CHANGE [MIN MAX] VALUES ACCORDING TO AOI FEATURES
hold on
colormap turbo
hold on
colorbar eastoutside
%%% SAVE MAP IN PNG FORMAT
saveas(gcf,[pwd '/Data/Results/MAP_SDB_BG.png']) 

%%% SAVE MAP IN TIFF FORMAT WITH THE PROPER EPSG CODE FOR THE GEOGRAPHICAL
% ZONE OF THE AOI

geotiffwrite( [pwd '/Data/Results/Map_SDB_BG'], Map_SDB_BG, R, CoordRefSysCode= Epsg);
close all

%%% SAVE WORKSPACE, CHANGE FILENAME TO AVOID OVERWRITING 
filename = "SDB";
save(filename)