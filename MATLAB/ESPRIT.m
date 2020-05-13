% Function Description: implements ESPRIT algorithm for direction of
% arrival estimation.


function thetas = ESPRIT(X, K)
% wavelength:
lambda = 2;
% sensor separation:
dist = 1;

% get dimensions of X:
dim = size(X);
N = dim(1);
T = dim(2);

% estimate autocorrelation matrix:
Rx = (1/T)*(X*X');

% perform eigendecomposition of autocorrelation matrix:
[V, D] = eig(Rx);
% sort eigenvectors in decreasing order of their eigenvalues:
[d,indices] = sort(diag(D), 'descend');
% D_sorted = D(indices, indices);
V_sorted = V(:,indices);
% extract first K eigenvectors (columns):
V1 = V_sorted(:, 1:K);

% W1 and W2 matrices:
W1 = V1(1:N-1, :);
W2 = V1(2:end, :);

% solve optomization by LS:
psi = inv(W1'*W1) * W1' * W2;

% perform eigendecomposition of psi:
[V_psi, D_psi] = eig(psi);

% extract theta_k's:
angles = angle(diag(D_psi));
sines = angles / (-2*pi*(dist/lambda));
thetas = asin(sines);
% convert to degrees:
thetas = (180/pi) * thetas;

end
