%% H2 optimal state feedback controller
clear all;
close all;

%System
m = [100 200 300 400 500];
A = [0 1; 0 0];
B1 = [0; 0];
for i = 1:length(m)
    B2(:,:,i) = [0; 1/m(i)];
end
C1 = [0 0];
C2 = [1 0];
D11 = 0; D12 = 0; D21 = 0; D22 = 0;

yalmip('clear');
Const = [];
eta = .0001;
n = size(A,1);
Z = sdpvar(1,n);
X = sdpvar(n,n);
W = sdpvar(1);
gam2 = sdpvar(1);
Const = [Const; X>=eta*eye(size(X))];
for i = 1:length(m)
    mat1 = [A B2(:,:,i)]*[X; Z]+[X Z']*[A'; B2(:,:,i)']+B1*B1';
    mat2 = [X (C1*X+D12*Z)'; C1*X+D12*Z W];
    
    Const = [Const; mat1 <= -eta*eye(size(mat1))];
    Const = [Const; mat2 >= eta*eye(size(mat2))];
    Const = [Const; trace(W)<= gam2];
end

opt=sdpsettings('solver','sedumi','verbose',0);
optimize(Const, gam2, opt);
ZZ = value(Z);
XX = value(X);
K = ZZ*inv(XX);
gamma = sqrt(value(gam2));
disp(K);
disp(gamma);