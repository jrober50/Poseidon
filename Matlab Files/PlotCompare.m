function [] = PlotCompare()

% Specific Formating Variables
width = 12;     % Width in inches
height = 8;    % Height in inches
fsz = 15;      % Fontsize


% Define Speed of Light
c = 29979245800;
csqr = c*c;







% Create Filenames
Results_File = "Data/Results.out";
Solution_File = "Data/Solution.out";

% Open Results File
Results_ID = fopen(Results_File);
if Results_ID == -1
    error('Cannot open file: %s',Results_File)
end
fgets(Results_ID);
Results = fscanf(Results_ID,'%f %f %f %f %f %f',[6 inf]);


% Open Solution File
Solution_ID = fopen(Solution_File);
if Solution_ID == -1
    error('Cannot open file: %s',Solution_File)
end
fgets(Solution_ID);
Solution = fscanf(Solution_ID,'%f %f %f',[3 inf]);




titlename = "Poseidon Results Compartive"; 
fig = figure;
sgt = sgtitle(titlename);
sgt.FontSize = 20;





% Plot Lapse and Conformal Factor
subplot(2,1,1)
plot(Results(1,:)/100000,Results(2,:));
hold on;
plot(Results(1,:)/100000,Results(3,:));
hold on;
plot(Solution(1,:)/100000,1 - Solution(2,:)/(2*csqr),'--');
hold on;
plot(Solution(1,:)/100000,1 + Solution(2,:)/(2*csqr),'--');

legend('\psi Solver','\alpha\psi Solver','\psi Analytic','\alpha\psi Analytic');

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('\psi & \alpha\psi Solver Values');
xlabel('Radial Distance (km)');
xlim([0,10e4]);



% Plot Shift1
subplot(2,1,2)
plot(Results(1,:)/100000,Results(4,:));
hold on;
plot(Solution(1,:)/100000,Solution(3,:),'--');

hold on;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

legend('\beta^1 Solver','\beta^1 Analytic');

title('\beta^1 Solver Value');
xlabel('Radial Distance (km)');
ylabel('Shift Value (cm/s)');
xlim([0,10e4]);




end

