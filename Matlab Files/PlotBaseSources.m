function [] = PlotBaseSources()
width = 10;     % Width in inches
height = 6;    % Height in inches
fsz = 15;      % Fontsize


filename = "Data/Base_Sources.out";
fid = fopen(filename);
if fid == -1
    error('Cannot open file: %s',filename)
end
fgets(fid);
data = fscanf(fid,'%f %f %f',[3 inf]);
%data = data';

c = 29979245800;
csqr = c*c;


% Plot rho
figure( 'Name','Base Source: rho, Density','NumberTitle','off');
plot(data(1,:)/100000,data(2,:));
hold on;

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size

set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('E');
xlabel('Radial Distance (km)');
ylabel('Density g cm^-3)');

hold off;


% Plot Velocity
figure('Name','Base Source: Velocity','NumberTitle','off');
plot(data(1,:)/100000,data(3,:)/c);

hold on;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Velocity');
xlabel('Radial Distance (km)');
ylabel('Velocity (c)');

hold off;



fclose(fid);
end

