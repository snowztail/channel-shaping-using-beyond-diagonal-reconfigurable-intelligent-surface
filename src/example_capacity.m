myfun = @(x,y) x.^3 .* y; 
[X,Y] = meshgrid(linspace(-3, 3, 30));
Z = myfun(X,Y);
Imx = find(Z == max(Z(:)));
Imn = find(Z == min(Z(:)));
figure
s1 = surf(X, Y, Z, 'FaceAlpha',0.8);
hold on
p1 = plot3(X(Imx), Y(Imx), Z(Imx), '^r', 'MarkerFaceColor','r', 'MarkerSize',10);
p2 = plot3(X(Imn), Y(Imn), Z(Imn), 'vg', 'MarkerFaceColor','g', 'MarkerSize',10);
hold off
grid on
xlabel('x')
ylabel('y')
zlabel('z')
legend([s1 p1(1) p2(1)], 'myfun(x,y)', 'Maximum', 'Minimum')
