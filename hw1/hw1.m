% Load data
data = load("specheat_data.mat")
C = data.Cj
D = data.Dj
T = data.Tj

% X axis range for plotting
t_range = 300:0.1:600

%%%%%%%%%%%%%%%%%%%%% Q1 %%%%%%%%%%%%%%%%%%%%%

% Corresponding data points
x = T([6:6:24]).'
y = C([6:6:24]).'

% Form the Vandermonde matrix
n = length(x);
Vand = ones(n,n);
for i = 1:n
    Vand(:,i) = (x.^(n-i));
end

% Solve
coeffs = Vand \ y

% Plot
figure(1); hold on
plot(t_range, polyval(coeffs, t_range), 'LineWidth',3)
plot(x,y,'kx','MarkerSize',16,'LineWidth',4)

%%%%%%%%%%%%%%%%%%%%% Q2 %%%%%%%%%%%%%%%%%%%%%

% Corresponding data points
% y2's are the coefficients
x2 = T([3:3:24]).'
y2 = C([3:3:24]).'

% Using Barycentric form 
pval = barycentric_interpolate(x2,y2,t_range)

% Plot
figure(2); hold on
plot(t_range,pval,'LineWidth',3)
plot(x2,y2,'kx','MarkerSize',16,'LineWidth',4)

%%%%%%%%%%%%%%%%%%%%% Q3 %%%%%%%%%%%%%%%%%%%%%

% Corresponding data points
x3 = T([8:8:24]).'
y3 = C([8:8:24]).'
d3 = D([8:8:24]).'

n = length(x3);

% Concatenate
y_with_der = [y3.' d3.'].'

% Form the Vandermonde-like matrix for Hermite interpolation
for i = 1:n
    for j = 1:(2*n)
        temp(j) = (x3(i).^(j-1))
        temp2(j) = ((x3(i).^(j-2)))
    end
    Herm(i,:) = temp 
    Herm(i+n,:) = temp2 .* [0:5]
end

% Solve
coeffs3 = Herm \ y_with_der
% Reverse the coeffs
coeffs3 = coeffs3(end:-1:1)

% Plot
figure(3); hold on
% polyval gives the desired polynomial
plot(t_range, polyval(coeffs3, t_range),'LineWidth',3)
plot(x3,y3,'kx','MarkerSize',16,'LineWidth',4)

figure(4); hold on
% polyder gives the derivative
plot(t_range, polyval(polyder(coeffs3), t_range),'LineWidth',3)
plot(x3,d3,'kx','MarkerSize',16,'LineWidth',4)

%%%%%%%%%%%%%%%%%%%%% Q4 %%%%%%%%%%%%%%%%%%%%%

% Corresponding data points
x4 = T([1:5:end]).'
y4 = C([1:5:end]).'
d4 = D([1:5:end]).'

% Repeat every element
x4_rep = repelem(x4,2)
y4_rep = repelem(y4,2)
d4_rep = repelem(d4,2)

% Form the divided difference table
[coeffs4, table] = divdif(x4_rep, y4_rep, d4_rep)

% Plot
figure(5); hold on
plot(t_range, evalnewt(t_range, x4_rep, coeffs4), 'LineWidth',3)
plot(x4,y4,'kx','MarkerSize',16,'LineWidth',4)
 
% New X axis range for the dimension match 
t_range2 = 300:0.1:599.9

% Plot
figure(6); hold on
p = polyfit(t_range,evalnewt(t_range, x4_rep, coeffs4),9)
plot(t_range2, polyval(polyder(p), t_range2),'LineWidth',3) 
plot(x4,d4,'kx','MarkerSize',16,'LineWidth',4)

%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%

function pval = barycentric_interpolate(xj,yj,x)

n = length(xj);

% (1) construction
 
for j = 1:n
    weight = 1;
    for k = 1:n
        if (k ~= j)
            weight = weight*(xj(j) - xj(k));
        end
    end
    
    w(j) = 1/weight;
    yw(j) = yj(j)*w(j);
end

% (2) evaluation

p = length(x);

num = zeros(1,p);
denom = zeros(1,p);
for j = 1:n
    temp = w(j)./(x-xj(j));
    num = num + yj(j)*temp;
    denom = denom + temp;
end

pval = num./denom;

return
end

function pval = evalnewt(x,xj,coef)

k = length(x);
np1 = length(xj);

pval = coef(np1)*ones(1,k);
for j = np1-1:-1:1
    pval = pval.*(x - xj(j)) + coef(j);
end
return;
end

function [coef,table] = divdif(xj,yj,dj)
%
% function [coef,table] = divdif(xi,yi)
%
% Construct a divided difference table based on data points (xi, yi).
% Upon return, the Newton interpolation coefficients are in coef

np1 = length(xj); n = np1-1;

table = zeros(np1,np1);
table(1:np1,1) = yj;

for k = 2:np1
    for j = k:np1
        % forms f[x_{j-k+1}, ... , x_j]
        if (k==2) & rem(j,2) == 0
            table(j,k) = dj(j)
        else
            table(j,k) = (table(j,k-1) - table(j-1,k-1)) / (xj(j) - xj(j-k+1));
        end
    end
end

coef = diag(table);

return;
end

% Final version