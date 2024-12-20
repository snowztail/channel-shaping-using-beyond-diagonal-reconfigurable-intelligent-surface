function [Q, S] = takagi_factorization(A)
    if norm(A - A.', 'fro') > 1e-10
        error('Matrix A must be symmetric.');
    end

    [U, S, V] = svd(A);
    t = diag(U' * conj(V));
    phi = angle(t) / 2;
    Q = U * diag(exp(1i * phi));
end
