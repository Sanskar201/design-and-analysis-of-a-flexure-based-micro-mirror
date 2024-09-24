% Student Name: Sanskar Anil Nalkande
% The code for calculating the Euler Stiffness Matrix is taken from the
% course modules on Bruinlearn

% Values of Material Constants
E = 169*10^9;
G = 94*10^9;
rho = 2330;

% Constructing the Constraint Matrix
Constraint = [ 0.0005 -0.00025 0.000175 -1  0 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
               0.0005  0.00025 0.000025 -1  0 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
              -0.0005 -0.00025 0.000025  1  0 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
              -0.0005  0.00025 0.000175  1  0 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
             -0.00025 -0.0005  0.000175  0  1 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
              0.00025 -0.0005  0.000025  0  1 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
             -0.00025  0.0005  0.000025  0 -1 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
              0.00025  0.0005  0.000175  0 -1 0 0 0 1 0.001 0.00005 0.00005 E G 1 0;
              1 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0];

% Using EulerStiffnessMatrix.m to calculate KK
[KK] = EulerStiffnessMatrix(Constraint);

% Mass of the stage
M_Stage = rho * 0.001 * 0.001 * 0.0002;

Z = zeros(3);
I = eye(3);

%D elta Matrix
Delta = [Z I; I Z];

% New axes at the center of mass of stage
n1 = [1 0 0].';
n2 = [0 1 0].';
n3 = [0 0 1].';

% Coordinates of the center of mass of the stage
L = [0 0 0.001];

z = [0 0 0].';
% N Matrix
N = [n1 n2 n3 z z z;
    cross(L,n1).' cross(L,n2).' cross(L,n3).' n1 n2 n3];

% Moments of inertia about center of mass
Ix = M_Stage*(0.001 * 0.001 + 0.0002 * 0.0002)/12;
Iy = M_Stage*(0.001 * 0.001 + 0.0002 * 0.0002)/12;
Iz = M_Stage*(0.001 * 0.001 + 0.001 * 0.001)/12;

% Inertia Matrix
In = diag([Ix Iy Iz M_Stage M_Stage M_Stage]);

% Mass Matrix
[MTW] = N * Delta * In / N;

% Calculating the eigenvalues and eigenvectors of inv(MTW)*K
[V,D] = eig(MTW\KK);

% Taking the square root and sorting them in ascending order
[E,X] = sort(diag(D^0.5));
V = V(:, X);

E