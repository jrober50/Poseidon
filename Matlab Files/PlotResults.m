function [] = PlotResults()

B1on = 1;
B2on = 0;

width = 10;     % Width in inches
height = 6;    % Height in inches
fsz = 15;      % Fontsize



filename = "Data/Results.out";
fid = fopen(filename);
if fid == -1
    error('Cannot open file: %s',filename)
end
fgets(fid);
data = fscanf(fid,'%f %f %f %f %f %f',[6 inf]);
%data = data';


% Plot Potential
titlename = 'Conformal Factor and Lapse Function';
figure( 'Name',titlename,'NumberTitle','off');
plot(data(1,:)/100000,data(2,:));
hold on;
plot(data(1,:)/100000,data(3,:));

legend('Psi','AlphaPsi');

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
set(gca,'XScale','log')
set(gca,'FontSize',fsz,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

title('Conformal Factor and Lapse Function');
xlabel('Radial Distance (km)');
ylabel('Potential (cm^2/s^2)');

hold off;


% Plot Shift1
if (B1on == 1 )
    titlename = 'Beta^1 Value';
    figure('Name',titlename,'NumberTitle','off');
    plot(data(1,:)/100000,data(4,:));

    hold on;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
    set(gca,'XScale','log')
    set(gca,'FontSize',fsz,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

    title('Beta^1 Solver Value');
    xlabel('Radial Distance (km)');
    ylabel('Shift Value (cm/s)');

    hold off;
end 


% Plot Shift2
if ( B2on == 1) 
    titlename = 'Beta^2 Value';
    figure('Name',titlename,'NumberTitle','off');
    plot(data(1,:)/100000,data(5,:));

    hold on;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*80, height*100]); %<- Set size
    set(gca,'XScale','log')
    set(gca,'FontSize',fsz,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',fsz,'fontWeight','bold')

    title('Beta^2 Solver Value');
    xlabel('Radial Distance (km)');
    ylabel('Shift Value (cm/s)');

    hold off;
end 



fclose(fid);
end

