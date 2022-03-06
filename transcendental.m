syms x

sigma = 7/4;
n = 5;

fplot((-x^2-sigma*n^2*x+n^2), "LineWidth", 3, "Color", 'b');
hold on
%fplot((x^2-sigma*n^2*x-n^2), "LineWidth", 1.5, "Color", "r")
ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
ax.GridLineStyle = ':';
ax.MinorGridLineStyle = 'none';

h1 = line([0, 0], [-50, 50], "LineWidth", 1.25, "Color", "k");
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h2 = line([n, n], [-50, 50], "LineWidth", 1.25, "Color", "k");
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h3 = line([-10, 10], [0, 0], "LineWidth", 1.25, "Color", "k");
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%grid on
fplot((x^2-sigma*n^2*x-n^2)*exp(-4*x), [0,10], "LineWidth", 3, "Color", 'g')
fplot((-x^2-sigma*n^2*x+n^2)+(x^2-sigma*n^2*x-n^2)*exp(-4*x), [0,n], "LineWidth", 3, "Color", 'r', "LineStyle", "-")
xlabel("\nu")
legend(["$f_1(\nu)$", "$f_2(\nu)$", "$f(\nu)$"], 'Interpreter', 'latex')
exportgraphics(ax,'s04-n5.pdf','ContentType','vector')

%%%% 
figure;
n = 10000;

fplot((-x^2-sigma*n^2*x+n^2), "LineWidth", 3, "Color", "b")
hold on
ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
ax.GridLineStyle = ':';
ax.MinorGridLineStyle = 'none';

h1 = line([0, 0], [-4*10^11, 10^11], "LineWidth", 1.25, "Color", "k");
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h2 = line([n, n], [-4*10^11, 10^11], "LineWidth", 1.25, "Color", "k");
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h3 = line([-10, n], [0, 0], "LineWidth", 1.25, "Color", "k");
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


fplot((x^2-sigma*n^2*x-n^2)*exp(-4*x), [0,n], "LineWidth", 3, "Color", "g")
fplot((-x^2-sigma*n^2*x+n^2)+(x^2-sigma*n^2*x-n^2)*exp(-4*x), [0,n], "LineWidth", 3, "Color", "r", "LineStyle", "-")
xlabel("\nu")
legend(["$f_1(\nu)$", "$f_2(\nu)$", "$f(\nu)$"], 'Interpreter', 'latex')
exportgraphics(ax,'s04-n10k.pdf','ContentType','vector')