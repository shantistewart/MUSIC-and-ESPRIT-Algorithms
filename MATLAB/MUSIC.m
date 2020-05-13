% Function Description: implements MUSIC algorithm for direction of
% arrival estimation.


function Smusic = MUSIC(X, K)
% wavelength:
lambda = 2;
% sensor separation:
dist = 1;
% angles to scan:
theta_scan_deg = -90:0.05:90; 
theta_scan = theta_scan_deg*pi/180;
dim_theta = size(theta_scan);
num_angles = dim_theta(2);

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
% extract last N-K eigenvectors (columns):
V2 = V_sorted(:, K+1:end);
disp(size(V2));

% calculate spectrum:
Smusic = zeros(num_angles, 1);
for k = 1 : num_angles
    % construct a(theta_k):
    a = zeros(N, 1);
    for n = 1 : N
        a(n, 1) = exp(1i*(-2*pi)*(dist/lambda)*(n-1)*sin(theta_scan(k)));
    end
    
    Smusic(k, 1) = 1 / norm(V2'*a);
end

end
