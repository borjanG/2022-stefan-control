syms x

sigma = 5;
n = 1000;
fplot((x^2/n^2+1)*tan(2*x), [0, 7], "LineWidth", 2, "Color", 'b')
hold on
%fplot((x^2-sigma*n^2*x-n^2), "LineWidth", 1.5, "Color", "r")
ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
%set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
ax.GridLineStyle = ':';
ax.MinorGridLineStyle = 'none';
fplot(x*sigma, [0, 7], "LineWidth", 1.5, "Color", "k")
xlabel("\nu")
%legend(["$f_1(\nu)$", "$f_2(\nu)$", "$f(\nu)$"], 'Interpreter', 'latex')
exportgraphics(ax,'tanfix.pdf','ContentType','vector')