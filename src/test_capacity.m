sigma_F = diag([2, 1]);
sigma_B = diag([4, 3]);

X0 = [1 0; 0 1];
X1 = [0 1; 1 0];
X2 = randn(2) + 1i*randn(2); X2 = X2 * (X2' * X2) ^ (-0.5);

H0 = sigma_B * X0 * sigma_F;
H1 = sigma_B * X1 * sigma_F;
H2 = sigma_B * X2 * sigma_F;

[U0,S0,V0] = svd(H0, 'vector');
[U1,S1,V1] = svd(H1, 'vector');
[U2,S2,V2] = svd(H2, 'vector');

snr = -20 : 5 : 20;
[R0, R1, R2] = deal(zeros(1, length(snr)));

for i = 1 : length(snr)
    P0 = waterfill(1, 1 ./ (S0.^2' * db2pow(snr(i))));
    P1 = waterfill(1, 1 ./ (S1.^2' * db2pow(snr(i))));
    P2 = waterfill(1, 1 ./ (S2.^2' * db2pow(snr(i))));

    R0(i) = sum(log2(1 + P0 .* S0.^2'));
    R1(i) = sum(log2(1 + P1 .* S1.^2'));
    R2(i) = sum(log2(1 + P2 .* S2.^2'));
end

R0
R1
