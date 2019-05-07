function [] = PlotSources()
width = 10;     % Width in inches
height = 6;    % Height in inches
fsz = 15;      % Fontsize


filename = "Data/Sources.out";
fid = fopen(filename);
if fid == -1
    error('Cannot open file: %s',filename)
end
fgets(fid);
data = fscanf(fid,'%f %f %f %f %f',[5 inf]);


% Plot E
figure( 'Name','Source: E, Matter-Energy Density','NumberTitle','off');
plot(data(1,:)/100000,data(2,:));
hold on;

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('E');
xlabel('Radial Distance (km)');
ylabel('Matter-Energy Density g/(cm s^2)');

hold off;


% Plot Si
figure('Name','Source: Si, Radial Momentum Density','NumberTitle','off');
plot(data(1,:)/100000,data(3,:));

hold on;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Radial Momentum Density');
xlabel('Radial Distance (km)');
ylabel('');

hold off;



% Plot S
figure('Name','Source: S, Trace of S^ij','NumberTitle','off');
plot(data(1,:)/100000,data(4,:));

hold on;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Trace of Momentum Density');
xlabel('Radial Distance (km)');
ylabel('');

hold off;

fclose(fid);
end

