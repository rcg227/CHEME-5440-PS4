% CHEME 5440, HW4 Q1
clear all
close all
% first set kappas=0.1, then to 10
K=0.1; 
K1=K; 
K2=K; 
K3=K; 
K4=K;

p=25; % number of points along x axis
%input=linspace(0.01,0.999,p); % where input=(1/KD)
input=linspace(0.01,100,p); % where input=(1/KD)
%input=linspace(0.01,5,p); % where input=(1/KD)


% initialize theta, x_roots, y_roots vectors
theta=zeros(1,p); 
x=zeros(1,p); 
y=zeros(1,p); 

a=zeros(1,p); 
b=zeros(1,p); 
c=zeros(1,p); 
d=zeros(1,p); 
e=zeros(1,p); 
f=zeros(1,p); 

for i=1:p
    %compute theta(i)
    theta(i)=(1/(1+(1/input(i)))); 
   
    % compute x_roots(i,:) [note roots will output a 1x2 vector]
    % Quadratic function: 0 = d(x^2) + ex + f
    d(i) = 5*theta(i) - 1; 
    e(i) = K1 + 1 + 5*theta(i)*K2 - 5*theta(i);
    f(i) = -5*theta(i)*K2; 
    
    %vector of only the positive roots of x 
    x(i)= -e(i) + sqrt(((e(i)^2)-(4*d(i)*f(i))))/(2*d(i)); 
   
    a(i) = -1 + 10*x(i);
    b(i) = -10*x(i) + 10*x(i)*K4 + K3 + 1;
    c(i) = -10*x(i)*K4;
    
    %vector of only the positive roots of y 
    y(i)= -b(i) + sqrt(((b(i)^2)-(4*a(i)*c(i))))/(2*a(i)); 
end
figure(1)
hold on 
subplot(2,2,1)
plot(input,theta,'b')
xlabel("Non-dimensional Input (1/KappaD)")
ylabel("ThetaB")
title("Kappa=0.1")
hold on 
subplot(2,2,2)
plot(input,x,'b')
xlabel("Non-dimensional Input (1/KappaD)")
ylabel("x*")
title("Kappa=0.1")
hold on 
subplot(2,2,3)
plot(input,y,'b')
xlabel("Non-dimensional Input (1/KappaD)")
ylabel("Output, y*")
title("Kappa=0.1")

hold off

%% Now let kappas=10
% CHEME 5440, HW4 Q1
clear all
close all
% first set kappas=0.1, then to 10
K=10; 
K1=K; 
K2=K; 
K3=K; 
K4=K;

p=25; % number of points along x axis
%input=linspace(0.001,0.999,p); % where input=(1/KD)
%input=linspace(0.01,10,p); % where input=(1/KD)
input=linspace(0.01,100,p); % where input=(1/KD)

% initialize theta, x_roots, y_roots vectors
theta=zeros(1,p); 
x=zeros(1,p); 
y=zeros(1,p); 

a=zeros(1,p); 
b=zeros(1,p); 
c=zeros(1,p); 
d=zeros(1,p); 
e=zeros(1,p); 
f=zeros(1,p); 

for i=1:p
    %compute theta(i)
    theta(i)=(1/(1+(1/input(i))));  
   
    % compute x_roots(i,:) [note roots will output a 1x2 vector]
    % Quadratic function: 0 = d(x^2) + ex + f
    d(i) = 5*theta(i) - 1; 
    e(i) = K1 + 1 + 5*theta(i)*K2 - 5*theta(i);
    f(i) = -5*theta(i)*K2; 
    
    %vector of only the positive roots of x 
    x(i)= -e(i) + sqrt(((e(i)^2)-(4*d(i)*f(i))))/(2*d(i)); 
    
    %compute y_roots(:,p)
    a(i) = -1 + 10*x(i);
    b(i) = -10*x(i) + 10*x(i)*K4 + K3 + 1;
    c(i) = -10*x(i)*K4;
    
    %vector of only the positive roots of y 
    y(i)= -b(i) + sqrt(((b(i)^2)-(4*a(i)*c(i))))/(2*a(i)); 
end
figure(2)
hold on 
subplot(2,2,1)
plot(input,theta,'b')
xlabel("Non-dimensional Input (1/KappaD)")
ylabel("ThetaB")
title("Kappa=10")
hold on 
subplot(2,2,2)
plot(input,x,'b')
xlabel("Non-dimensional Input (1/KappaD)")
ylabel("x*")
title("Kappa=10")
hold on 
subplot(2,2,3)
plot(input,y,'b')
xlabel("Non-dimensional Input (1/KappaD)")
ylabel("Output, y*")
title("Kappa=10")

hold off