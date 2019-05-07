function [] = PlotSolution()

width = 10;     % Width in inches
height = 6;    % Height in inches
fsz = 15;      % Fontsize

filename = "Data/Solution.out";
fid = fopen(filename);
if fid == -1
    error('Cannot open file: %s',filename)
end
fgets(fid);
data = fscanf(fid,'%f %f %f',[3 inf]);



c = 29979245800;
csqr = c*c;


% Plot Potential
figure( 'Name','Analytic Conformal Factor and Lapse Function','NumberTitle','off');
plot(data(1,:)/100000,1 - data(2,:)/(2*csqr));
hold on;
plot(data(1,:)/100000,1 + data(2,:)/(2*csqr));

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

legend('Psi','AlphaPsi');

title('Analytic Conformal Factor and Lapse Function');
xlabel('Radial Distance (km)');
ylabel('Potential (cm^2/s^2)');

hold off;


% Plot Shift
figure('Name','Analytic Beta^1 Value','NumberTitle','off');
plot(data(1,:)/100000,data(3,:));

hold on;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Beta^1 Analytic Value');
xlabel('Radial Distance (km)');
ylabel('Shift Value (cm/s)');

hold off;

fclose(fid);
end

