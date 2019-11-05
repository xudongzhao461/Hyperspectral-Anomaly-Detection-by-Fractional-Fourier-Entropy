function D = RX(X)

[N M] = size(X);

X_mean = mean(X.').';

X = X - repmat(X_mean, [1 M]);

Sigma = (X * X')/M;

Sigma_inv = inv(Sigma);
for m = 1:M
 D(m) = X(:, m)' * Sigma_inv * X(:, m);
end
