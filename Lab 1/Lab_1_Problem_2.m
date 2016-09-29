%% Jeffrey Huang - Linear Systems and Signals Lab 1 - Problem 2
%%

clc; clear all; close all;

%% Part A
%   y'(t) + 2 * y(t) = 3 * x'(t) + 2 * x(t) --> (3s+2)/(s+2)
%           --> a = 2, b0 = 3, b1 = 2, y(0-)= -1
%           --> scale factor = (b1 - a * b0) = 1 - 2 * 3 = -4
%   v'(t) =-2 * v(t) + x(t) --> a = 2
%   y (t) =-4v(t) + 3 * x(t) --> (b1 - a * b0) = -4, b0 = 3
%           --> b1 - 2 * 3 = -4 --> b1 = 2
% Equations (2) implies equation (1)
a  = 2;
b0 = 3;
b1 = 2;
y0 = -1;
v0 = y0/(b1-a*b0);
%% Part B
syms t s Y
% --- determine H(s) and h(t) ---
H = (b0*s+b1)/(s+a);
H = partfrac(H);        % H(s) = 3 - 4/(s + 2)
h = ilaplace(H);        % h(t) = 3*dirac(t) - 4*exp(-2*t)
%% Part C
% The initial value of the state, v(0-) = v(0+) is related to the intial
% value , y(0-) through being the same value intially.
%% Part D
% --- determine complete y(t) using Laplace ---
x = exp(-t);        % input signal
X = laplace(x);     % X(s)
Y1 = solve(s*Y-y0 + a*Y == b0*s*X + b1*X, Y); % y(t) - Laplace method
Y1 = partfrac(Y1,s); % Y(s) = 3/(s + 2) - 1/(s + 1)
y1 = ilaplace(Y1);   % y(t) = 3*exp(-2*t) - exp(-t)
%% Part E
% --- determine complete y(t) using dsolve ---
syms y(t)
x0 = subs(x,t,0);
dy = diff(y,t);
dx = diff(x,t);
y2 = dsolve(dy + a*y == b0*dx + b1*x, y(0) == y0+b0*x0); % y(t) = 3*exp(-2*t) - exp(-t)
yzs1 = ilaplace(H*X); % zero-state component
%% Part F
% -- determine complete V(s), v(t) ---
syms V
V1 = solve(s*V-v0 + a*V == X, V);   % V(s) = (1/(s + 1) + 1/3)/(s + 1)
V1 = partfrac(V1)     % V(s) = 1/(s + 1) - 3/(4*(s + 2))
v1 = ilaplace(V1)     % v(t) = exp(-t)/3 + t*exp(-t)
simplify(y2 + 4*v1 - 3*x)    % verify: y = -4*v(t) - 3*x(t), ans = 0
%% Part G
% --- plotting y(t) and v(t) over the interval 0 <= t <= 5
t1 = linspace(0,5,201);
y1p = subs(y1,t,t1);
v1p = subs(v1,t,t1);
figure;
plot(t1, y1p, 'b-', t1, v1p, 'r--');
axis([0 5 -.25 2.25]);
title('exact outputs, {\ity}(0-) = -1');
legend('{\ity}({\itt})', ' {\itv}({\itt})');
%% Part H
% --- state-space implementation ---
t2 = linspace(0,5,201)';
x1 = exp(-t2);
num = [b0,b1];
den = [1,a];
H = tf(num,den);
[A,B,C,D] = tf2ss(num,den);
S = ss(A,B,C,D);
s = tf('s');
Hv = 1/(s+a);
Sv = ss(Hv);
y3 = lsim(S,x1,t2,v0);
yzs = lsim(S,x1,t2);
v2 = lsim(Sv,x1,t2,v0);
norm(-4*v2+3*x1 - y3)
figure;
plot(t2, y3, 'b-', t2, v2, 'r--');
axis([0 5 -.25 2.25]);
title('lsim outputs, {\ity}(0-) = -1');
legend('{\ity}({\itt})', ' {\itv}({\itt})');
%% Part I
% --- discrete-time implementation ---
syms y
T = 0.01;
A1 = -exp(-a*T);
B1 = b1*(1-exp(-a*T))/a-b0;
tn = 0:T:5;
N = length(tn);
x = exp(-tn);
w = y0;
v = 0;
for n=0:N-1,             
    y(n+1) = -A1*w + b0*x(n+1) + B1*v;
    w = y(n+1);
    v = x(n+1);
end

w = 0; v = 0;
for n=0:N-1,             
    yzs(n+1) = -A1*w + b0*x(n+1) + B1*v;     % (n+1) is MATLAB index
    w = yzs(n+1);
    v = x(n+1);
end

clear v;
A1 = -exp(-a*T);
B1 = (1 - exp(-a*T))/a;
w = v0; u = 0;
for n=0:N-1,             
    v(n+1) = -A1*w + B1*u;
    w = v(n+1);
    u = x(n+1);
end
figure;
plot(tn,y,'b-', tn, v, 'g-', t2, y3, 'r:', t2, v2, 'r:');
title('discrete-time outputs,  {\ity}(0-) = -1, {\itT} = 0.1');
legend('{\ity}({\itt_n})', ' {\itv}({\itt_n})', 'exact');
%% Part J
% --- Everything above from part C to part i where y0 = 0 ---
clear all;
syms t s Y
a  = 2;
b0 = 3;
b1 = 2;
y0 = 0;
v0 = y0/(b1-a*b0);
H = (b0*s+b1)/(s+a);
H = partfrac(H);
h = ilaplace(H);
%% Part J-C
% The initial value of the state, v(0-) = v(0+) is related to the intial
% value , y(0-) through being the same value intially.
%% Part J-D
% --- determine complete y(t) using Laplace ---
x = exp(-t);        % input signal
X = laplace(x);     % X(s)
Y1 = solve(s*Y-y0 + a*Y == b0*s*X + b1*X, Y); % y(t) - Laplace method
Y1 = partfrac(Y1,s); % Y(s) = 3/(s + 2) - 1/(s + 1)
y1 = ilaplace(Y1);   % y(t) = 3*exp(-2*t) - exp(-t)
%% Part J-E
% --- determine complete y(t) using dsolve ---
syms y(t)
x0 = subs(x,t,0);
dy = diff(y,t);
dx = diff(x,t);
y2 = dsolve(dy + a*y == b0*dx + b1*x, y(0) == y0+b0*x0); % y(t) = 3*exp(-2*t) - exp(-t)
yzs1 = ilaplace(H*X); % zero-state component
%% Part J-F
% -- determine complete V(s), v(t) ---
syms V
V1 = solve(s*V-v0 + a*V == X, V);   % V(s) = (1/(s + 1) + 1/3)/(s + 1)
V1 = partfrac(V1)     % V(s) = 1/(s + 1) - 3/(4*(s + 2))
v1 = ilaplace(V1)     % v(t) = exp(-t)/3 + t*exp(-t)
simplify(y2 + 4*v1 - 3*x)    % verify: y = -4*v(t) - 3*x(t), ans = 0
%% Part J-G
% --- plotting y(t) and v(t) over the interval 0 <= t <= 5
t1 = linspace(0,5,201);
y1p = subs(y1,t,t1);
v1p = subs(v1,t,t1);
figure;
plot(t1, y1p, 'b-', t1, v1p, 'r--');
axis([0 5 -.25 3.25]);
title('exact outputs, {\ity}(0-) = 0');
legend('{\ity}({\itt})', ' {\itv}({\itt})');
%% Part J-H
% --- state-space implementation ---
t2 = linspace(0,5,201)';
x1 = exp(-t2);
num = [b0,b1];
den = [1,a];
H = tf(num,den);
[A,B,C,D] = tf2ss(num,den);
S = ss(A,B,C,D);
s = tf('s');
Hv = 1/(s+a);
Sv = ss(Hv);
y3 = lsim(S,x1,t2,v0);
yzs = lsim(S,x1,t2);
v2 = lsim(Sv,x1,t2,v0);
norm(-4*v2+3*x1 - y3)
figure;
plot(t2, y3, 'b-', t2, v2, 'r--');
axis([0 5 -.25 3.25]);
title('lsim outputs, {\ity}(0-) = 0');
legend('{\ity}({\itt})', ' {\itv}({\itt})');
%% Part J-I
% --- discrete-time implementation ---
syms y
T = 0.01;
A1 = -exp(-a*T);
B1 = b1*(1-exp(-a*T))/a-b0;
tn = 0:T:5;
N = length(tn);
x = exp(-tn);
w = y0;
v = 0;
for n=0:N-1,             
    y(n+1) = -A1*w + b0*x(n+1) + B1*v;
    w = y(n+1);
    v = x(n+1);
end

w = 0; v = 0;
for n=0:N-1,             
    yzs(n+1) = -A1*w + b0*x(n+1) + B1*v;     % (n+1) is MATLAB index
    w = yzs(n+1);
    v = x(n+1);
end

clear v;
A1 = -exp(-a*T);
B1 = (1 - exp(-a*T))/a;
w = v0; u = 0;
for n=0:N-1,             
    v(n+1) = -A1*w + B1*u;
    w = v(n+1);
    u = x(n+1);
end
figure;
plot(tn,y,'b-', tn, v, 'g-', t2, y3, 'r:', t2, v2, 'r:');
title('discrete-time outputs,  {\ity}(0-) = 0, {\itT} = 0.1');
legend('{\ity}({\itt_n})', ' {\itv}({\itt_n})', 'exact');
axis([0 5 -.25 3.25]);