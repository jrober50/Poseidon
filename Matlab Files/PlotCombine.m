function [] = PlotCombine()

% Specific Formating Variables
width = 12;     % Width in inches
height = 8;    % Height in inches
fsz = 15;      % Fontsize


% Define Speed of Light
c = 29979245800;
csqr = c*c;







% Create Filenames
BaseSources_File = "Data/Base_Sources.out";
Results_File = "Data/Results.out";

% Open Base Source File
BaseSource_ID = fopen(BaseSources_File);
if BaseSource_ID == -1
    error('Cannot open file: %s',BaseSources_File)
end
fgets(BaseSource_ID);
BaseSources = fscanf(BaseSource_ID,'%f %f %f',[3 inf]);


% Open Results File
Results_ID = fopen(Results_File);
if Results_ID == -1
    error('Cannot open file: %s',Results_File)
end
fgets(Results_ID);
Results = fscanf(Results_ID,'%f %f %f %f %f %f',[6 inf]);

titlename = "Poseidon Results"; 
fig = figure;
sgt = sgtitle(titlename);
sgt.FontSize = 20;







% Plot rho
subplot(2,2,1);
plot(BaseSources(1,:)/100000,BaseSources(2,:));
hold on;

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Density');
xlabel('Radial Distance (km)');
ylabel('Density (g cm^{-3})');
ylim([0,1e14]);
%xlim([1e1,1e4]);
xlim([0,10e4]);

% Plot Radial Velocity
subplot(2,2,2);
plot(BaseSources(1,:)/100000,BaseSources(3,:)/c);
hold on;

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Velocity');
xlabel('Radial Distance (km)');
ylabel('Velocity (c)');
ylim([-0.15 0]);
xlim([0,10e4]);


% Plot Lapse and Conformal Factor
subplot(2,2,3)
plot(Results(1,:)/100000,Results(2,:));
hold on;
plot(Results(1,:)/100000,Results(3,:));
hold on;
%plot(Results(1,:)/100000,Results(3,:)./Results(2,:));

legend('\psi','\alpha\psi');

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('\psi & \alpha\psi Solver Values');
xlabel('Radial Distance (km)');
ylim([0.75 1.2]);
xlim([0,10e4]);



% Plot Shift1
subplot(2,2,4)
plot(Results(1,:)/100000,Results(4,:));

hold on;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('\beta^1 Solver Value');
xlabel('Radial Distance (km)');
ylabel('Shift Value (cm/s)');
ylim([-3e-10,1e-10]);
xlim([0,10e4]);




end

