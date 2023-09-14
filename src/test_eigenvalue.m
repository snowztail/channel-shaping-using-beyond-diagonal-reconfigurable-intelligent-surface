clear; close; setup;

[N_t, N_s, N_r] = deal(2, 2, 2);

H_f = sqrt(0.5) * (randn(N_s, N_t) + 1i * randn(N_s, N_t));
H_b = sqrt(0.5) * (randn(N_r, N_s) + 1i * randn(N_r, N_s));

[U_f, S_f, V_f] = svd(H_f);
[U_b, S_b, V_b] = svd(H_b);

M = 1e2;
[sv_diagonal, sv_unitary] = deal(zeros(M, M, 2));
[X, Y] = meshgrid(2 * pi * (1:M) / M);
for i = 1 : M
    for j = 1 : M
        sv_diagonal(i, j, :) = svd(H_b * ris_diagonal(2 * pi * [i, j] / M) * H_f);
        sv_unitary(i, j, :) = svd(H_b * ris_unitary(2 * pi * [i, j] / M) * H_f);
    end
end

tiledlayout(2, 1);

nexttile;
Z = sv_diagonal(:, :, 1);
[I_max, I_min] = deal(find(Z == max(Z(:))), find(Z == min(Z(:))));
S(1) = mesh(X, Y, Z);
hold on;
datatip(S(1), 'DataIndex', I_max(1));
datatip(S(1), 'DataIndex', I_min(1));
S(1).DataTipTemplate.DataTipRows = S(1).DataTipTemplate.DataTipRows(end);
S(1).DataTipTemplate.DataTipRows.Format = '%.3f';
% plot3(X(I_max), Y(I_max), Z(I_max), '^r');
% plot3(X(I_min), Y(I_min), Z(I_min), 'vg');
Z = sv_diagonal(:, :, 2);
[I_max, I_min] = deal(find(Z == max(Z(:))), find(Z == min(Z(:))));
S(2) = mesh(X, Y, Z);
datatip(S(2), 'DataIndex', I_max(1));
datatip(S(2), 'DataIndex', I_min(1));
S(2).DataTipTemplate.DataTipRows = S(2).DataTipTemplate.DataTipRows(end);
S(2).DataTipTemplate.DataTipRows.Format = '%.3f';
% plot3(X(I_max), Y(I_max), Z(I_max), '^r');
% plot3(X(I_min), Y(I_min), Z(I_min), 'vg');
hold off;
legend(S, '$\sigma_1$', '$\sigma_2$');
xlabel('$\theta_1$');
ylabel('$\theta_2$');
zlabel({'Channel'; 'singular value'});
xticks([-2*pi, -1.5*pi, -pi, -0.5*pi, 0, 0.5*pi, pi, 1.5*pi, 2*pi]);
yticks([-2*pi, -1.5*pi, -pi, -0.5*pi, 0, 0.5*pi, pi, 1.5*pi, 2*pi]);
xticklabels({'-2$\pi$', '-1.5$\pi$', '-$\pi$', '-0.5$\pi$', '0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'});
yticklabels({'-2$\pi$', '-1.5$\pi$', '-$\pi$', '-0.5$\pi$', '0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'});
title('2x2 MIMO + Diagonal 2-Element RIS');


nexttile;
Z = sv_unitary(:, :, 1);
[I_max, I_min] = deal(find(Z == max(Z(:))), find(Z == min(Z(:))));
S(1) = mesh(X, Y, Z);
hold on;
datatip(S(1), 'DataIndex', I_max(1));
datatip(S(1), 'DataIndex', I_min(1));
S(1).DataTipTemplate.DataTipRows = S(1).DataTipTemplate.DataTipRows(end);
S(1).DataTipTemplate.DataTipRows.Format = '%.3f';
% plot3(X(I_max), Y(I_max), Z(I_max), '^r');
% plot3(X(I_min), Y(I_min), Z(I_min), 'vg');
Z = sv_unitary(:, :, 2);
[I_max, I_min] = deal(find(Z == max(Z(:))), find(Z == min(Z(:))));
S(2) = mesh(X, Y, Z);
datatip(S(2), 'DataIndex', I_max(1));
datatip(S(2), 'DataIndex', I_min(1));
S(2).DataTipTemplate.DataTipRows = S(2).DataTipTemplate.DataTipRows(end);
S(2).DataTipTemplate.DataTipRows.Format = '%.3f';
% plot3(X(I_max), Y(I_max), Z(I_max), '^r');
% plot3(X(I_min), Y(I_min), Z(I_min), 'vg');
hold off;
legend(S, '$\sigma_1$', '$\sigma_2$');
xlabel('$\alpha$');
ylabel('$\phi$');
zlabel({'Channel'; 'singular value'});
xticks([-2*pi, -1.5*pi, -pi, -0.5*pi, 0, 0.5*pi, pi, 1.5*pi, 2*pi]);
yticks([-2*pi, -1.5*pi, -pi, -0.5*pi, 0, 0.5*pi, pi, 1.5*pi, 2*pi]);
xticklabels({'-2$\pi$', '-1.5$\pi$', '-$\pi$', '-0.5$\pi$', '0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'});
yticklabels({'-2$\pi$', '-1.5$\pi$', '-$\pi$', '-0.5$\pi$', '0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'});
title('2x2 MIMO + Unitary 2-Element RIS');




function [Theta] = ris_diagonal(arg)
    Theta = diag([exp(1i * arg(1)), exp(1i * arg(2))]);
end

function [Theta] = ris_unitary(arg)
    Theta = [exp(1i * arg(1)) * cos(arg(2)), 1i * sin(arg(2)); 1i * sin(arg(2)), exp(-1i * arg(1)) * cos(arg(2))];
end
