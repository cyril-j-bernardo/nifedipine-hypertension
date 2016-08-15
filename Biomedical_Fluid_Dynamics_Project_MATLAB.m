%Pressure Waveform for a Healthy Patient

%Since there were a range of values for R,C and Z, this script
%attempts to find the proper set of values that will result in a maximum
%pressure of 120mmHg and a minimum pressure of 80mmHg

clear

Tc = 60/75;                                     % Total time for one cardiac cycle
Ts = (2*Tc) ./ 5;                                  % Total Time for Systolic Phase

i = 1;
t = 0:0.01:Tc;                                  % Discretization of times in one cardiac cycle
q(size(t,2)) = 0;                               % Initializing of the vector 

q0 = 70 ./(Ts/2);                                 % Peak Flow

for x = 0:0.01:Ts                               % Setting the dependent variable
    q(i) = q0 .* (sin(pi*x ./ Ts)) .^ 2;                 % Dependent Variable q
    i = i + 1;
end

Q = [q,q(2:end)];                               % New Q Matrix of two cycles
Q = Q(1:size(Q,2));
T = 0:0.01:2*Tc;                                % Matrix of 2 cardiac cycles
flag = 0;

for R = 0.89:0.01:1.99                          % Various R Values
    for Z = 0.03:0.01:0.06                      % Various Z Values
        for C = 0.47:0.01:1.86                  % Various C Values

            A = -1/(R*C);
            B = 1;
            C_c = (R+Z)/(R*C);
            D = Z;

            sys = ss(A,B,C_c,D);                % State-Space Model

            Y = lsim(sys,Q,T,((80-D)/C_c));     % Time Response given initial condition
            peak = max(Y);                      % Peak Value
            i = i + 1;      
            if((peak > 119) && (peak < 121))    
                peak2 = max(Y(60:end));
                if((peak2 > 119) && (peak2 < 121))
                    flag = 1;
                    break;
                end
            end
        end
       if (flag == 1)
           break;
       end 
    end
    if (flag == 1)
        plot(T,Y)
        title('Pressure Waveform for Healthy Patient')
        xlabel('Time (s)')
        ylabel('Pressure (mmHg)')
        break;
    end
end


fprintf('\nR value calculated: %-3.2f' ,R)
fprintf('\nC value calculated: %-3.2f' ,C)
fprintf('\nZ value calculated: %-3.2f \n' ,Z)

%%

%Pressure Waveform for Patient with Hypertension
   
clear

R = 1.311;
C = 1.234;
Z = 0.069;

Tc = 60/73;
Ts = (2*Tc)/5;

i = 1;
t = 0:0.01:Tc;
q(size(t,2)) = 0;

q0 = 93/(Ts/2);

for x = 0:0.01:Ts
    q(i) = q0*(sin(pi*x/Ts))^2;
    i = i + 1;
end

Q = [q,q(2:end)];
Q = Q(1:size(Q,2));
T = 0:0.01:2*Tc;


A = -1/(R*C);
B = 1;
C_c = (R+Z)/(R*C);
D = Z;

sys = ss(A,B,C_c,D);

P = lsim(sys,Q,T,((134-D)/C_c));

figure

plot(T,P)
title('Pressure Waveform for Patient with Hypertension')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
axis([0,1.64,130,210])

%%

%Pressure Waveform for Patient Using Nifedipine to Treat Hypertension

clear

R = 0.97;
C = 1.776;
Z = 0.069;

Tc = 60/73;
Ts = (2*Tc)/5;

i = 1;
t = 0:0.01:Tc;
q(size(t,2)) = 0;

q0 = 93/(Ts/2);

for x = 0:0.01:Ts
    q(i) = q0*(sin(pi*x/Ts))^2;
    i = i + 1;
end

Q = [q,q(2:end)];
Q = Q(1:size(Q,2));
T = 0:0.01:2*Tc;


A = -1/(R*C);
B = 1;
C_c = (R+Z)/(R*C);
D = Z;

sys = ss(A,B,C_c,D);

figure

P = lsim(sys,Q,T,((104-D)/C_c));
plot(T,P)
title('Pressure Waveform for Patient using Nifedipine to treat Hypertension')
xlabel('Time (s)')
ylabel('Pressure (mmHg)')
axis([0,1.64,100,170])