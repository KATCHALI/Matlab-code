# Matlab-code
Mathematical numerical simulation of the BSFL based fertilizer Model

% Function describing variable dynamics of larval growth and nutrients
function Katchali_MSC()
clear; clc;
%--------------------------------------------------------------------------
% Constants
%--------------------------------------------------------------------------
V_i = 280; %inner production volume [m^3]
k_air = 1.005; % Specific heat capacity of air in [J/g·°C]
k_w = 4.2;     %water specific heat at standard condition in [J/g·°C]
k_wv = 1.996;  % water vapor specific heat
rho_air = 1.293;  % density of air in standard condition
k_D_1 = 6.35; % duration of highest feed assimilation (2nd to 4th instars)
k_D_2 = 2.5; % duration when the assimilation decreases (pre-pupae)
k_D_3 = 2.5; % duration to when the assimilation ended (pupae) in hour
B_max = 0.69; %Maximum individual larvae mass in [mg]
W = 60; % substrate moisture in %
A_m = 0.55*0.38; % Tray surface in [m^2]
k_hl= 2260; %(kJ/kg, The latent heat capacity of water at [100^oC])
W_c = 0.75; %l/kg water content in the substrate
k_Qass = 1.4e3; % specific heat due to assimilation in J/g
k_Qmai = 28e3; % specific heat due to maintenance in J/g
k_Qma = 3.1e3; % sepcific heat due to maturity in J/g
A_ai = 1.98;
n_lar = 200;
N_B = 102.5;
N_F = 15;
P_B = 92.5;
P_F = 15;
K_B = 132.5;
K_F = 15;
S = 10;  
k = 2;  %slop of the temperature rate between optimal and bounds 
T_min = 25;
T_max = 42; 
%--------------------------------------------------------------------------
% Load and process real data from a CSV file
%--------------------------------------------------------------------------
data = readtable('mean_by_timewell.csv'); % Provided data
tableHeaderNames = data.Properties.VariableNames; % Table headers
t_data = data.Time; % Time column
%Normalization
observed_data = [
    (data.Bdry - min(data.Bdry))./ (max(data.Bdry) - min(data.Bdry)),
    (data.Tai - min(data.Tai)) ./ (max(data.Tai) - min(data.Tai)),
    (data.Tsub - min(data.Tsub)) ./ (max(data.Tsub) - min(data.Tsub)),
    (data.Hai - min(data.Hai)) ./ (max(data.Hai) - min(data.Hai)),
    (data.Cai - min(data.Cai)) ./ (max(data.Cai) - min(data.Cai)),
    (data.N - min(data.N)) ./ (max(data.N) - min(data.N)),
    (data.P - min(data.P))./ (max(data.P) - min(data.P)),
    (data.K - min(data.K)). / (max(data.K) - min(data.K))
];

%--------------------------------------------------------------------------
% Parameter estimation
%--------------------------------------------------------------------------
initial_parameter_guesses =[
    2.000;  % k_r_maxgm
    2.128;  % k_r_maxA
    2.500;  % sigma
    0.01877  % k_A_half
    0.00049; % F_dhalf
    1.873;  % r_C_1
    1.873;  % r_C_2
    1.5;    % k_dev
    1.00567; % h_m_a
    1.61e-4; % k_ing
    1.789; % k_mle
    2e3;    % k_sub
    5e6;    % k_le
    0.5762; % k_alpha_ex
    1.0879;  % k_r_maxW
    0.02135; % k_alpha_ass
    5.6779e-4; % k_main
    9.6779e-4; % k_ma
    32;     % T_opt
    5.879;      % k_r_maxT
    3.97;    %x
    44.6;
    69;
    78;
    ];
Y = objectiveFunc(optimal_params, t_data);

% Optimization options for lsqcurvefit
optimization_options = optimoptions(@lsqcurvefit, ...
    'StepTolerance', 1e-10, 'Display', 'iter-detailed', 'FunctionTolerance', 1e-12);
% Specify lower bound (lb) and upper bound (ub) for the model parameters
lb = [0.1,1,0.0065,0.1,0,1,1,2,1,0,1,1.5e3,4e6,0,1,0,1,0,25,0.5,0, 0, 0, 0]; %lower parameters bounds
ub = [2.0,1.5,5.76,1,1,3.873,3.873,4,2,1,2.5,3e3,6e6,1.5,2.5,1,1,1,35,1.5, 1, 100, 1600, 950];%upper parameters bounds

% Use lsqcurvefit to find optimal parameters
[optimal_params,resnorm,residuals,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@objectiveFunc,...
    initial_parameter_guesses,t_data,observed_data, lb, ub, optimization_options);
% Display optimal parameters
fprintf(1,'\tOptimal Parameters:\n')
for i = 1:length(optimal_params)
    fprintf(1, '\t\tparams(%d) = %10.6f\n', i, optimal_params(i))
end
%--------------------------------------------------------------------------
% Solve the system using RK4 (ode45)
%--------------------------------------------------------------------------
t = linspace(min(t_data), max(t_data), 1000);
[~,n] = size(observed_data); iter = 4; m = length(t);
BigMatrix = zeros(m*iter,n); variedParamLegend = {iter};
variedParam = W ; %<-- CHANGE TARGET PARAMETER
variedParamName = 'W '; %<-- CHANGE TARGET PARAMETER
for i = 1:iter
    variedParamLegend{i} = strcat(variedParamName,'=',num2str(variedParam));
    Y = objectiveFunc(optimal_params, t); % Optimal solutions
    BigMatrix(m*(i-1)+1:i*m,:) = Y;
    variedParam = variedParam+5 ;
    W  = variedParam; %<-- CHANGE TARGET PARAMETER
    
end
fprintf('R^2: %.6f\n',R_squared)
fprintf('rmse_N: %.6f\n',rmse_N)
%--------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------- 
color = {'b','k','g','m'};
for i = [1, 6:8]  % Plot for the first variable (index 1) and last three variables (indices 6, 7, 8)
    figure(i)
    % plot(t_data, observed_data(:,i), 'ro', 'DisplayName',[tableHeaderNames{i+1},'-Data'])
    for j = 1:iter
        hold on
        plot(t, BigMatrix(m*(j-1)+1:j*m,i), color{j}, 'LineWidth', 1, 'DisplayName', [tableHeaderNames{i+1}, '-Simulated']);
        grid on
    end
    xlabel('Time')
    ylabel('Values')
    legend([tableHeaderNames{i+1},'-Data'],variedParamLegend{1:iter},'Location','N')
    switch i
        case 1
            ylabel('Larval weight [mg]')
        case 6
            ylabel('Nitrogen')
        case 7
            ylabel('Phosphorus')
        case 8
            ylabel('Potassium')
    end
end

% Magnify
figHandler = figure(8);
magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'magnifierShape', 'ellipse',...
        'initialPositionSecondaryAxes', [326.933 259.189 164.941 102.65],...
        'initialPositionMagnifier',     [174 49 14 37],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 1,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.91 0.91 0.91]... 
            );

%--------------------------------------------------------------------------
% Define the system of ODEs
%--------------------------------------------------------------------------
    function sol = objectiveFunc(params,tspan)
        % Parameters to be estimated
        k_r_maxgm = params(1);  %maximal larvae growth rate for feed availability
        k_r_maxA = params(2);   %maximal larvae growth rate with airflow rate (O2)
        sigma = params(3);  %the mass transfer coefficient for cond-evap
        k_A_half = params(4); %airflow rate leading to half maximal growth rate
        F_dhalf = params(5);  %half feed saturation rate
        r_C_1 = params(6);%specific rate of CO2 production due to assimilation 
        r_C_2 = params(7);%maintenance metabolism
        k_dev = params(8);    %specific development rate
        h_m_a = params(9);%the convective heat transfer coefficient between the air and the growingmedium,
        k_ing = params(10);   %[mg s-1] feed ingestion rate
        k_mle = params(11); %coefficient of leakage
        k_sub = params(12); %substrate heat coeffient
        k_le = params(13);  %Co2 leakage coefficient
        k_alpha_ex= params(14);  %[mg s-1] feed excretion rate
        k_r_maxW = params(15); %maximum growth with respect to substrate moisture 
        k_alpha_ass = params(16); %specific assimilationn rate
        k_main = params(17);     %specific maintenance rate
        k_ma = params(18);       %specific maturity rate
        T_opt = params(19);     % Optimal temperature for growing
        k_r_maxT = params(20);  %maximal larvae growth rate with respect to temperature rate change 
        x= params(21);
        theta= params(22);
        omega=params(23);
        delta=params(24);
        % Nitrogen variation coefficient for C:N ratio
        %--------------------------------------------------------------------------
        % Solve the system
        %--------------------------------------------------------------------------
        % Initial Conditions
        Bdry0 = observed_data(1,1);
        Tai0 = observed_data(1,2);
        Tsub0 = observed_data(1,3);
        Hai0 = observed_data(1,4);
        Cai0 = observed_data(1,5);
        N0 = observed_data(1,6);
        P0 = observed_data(1,7);
        K0 = observed_data(1,8);
        IC = [Bdry0,Tai0,Tsub0,Hai0,Cai0,N0,P0,K0];
        % IC = ones(1,8);
        % Solve the system using ode45()
        options = odeset('RelTol',1e-3,'AbsTol',1e-4);
        
        [T,Y] = ode45(@f,tspan,IC,options);
        
        function dYdt = f(t,Y)
            % Calculation of the rates
            %F_d = S./Y(1); %feed availability
            tp=504;%growing period from second instar to last instar in [h]
            F_d =((S .* 1000) ./ n_lar)./ tp;%feed density
            r_Tsub = k_r_maxT./ (1 + exp((T_opt - Y(3)) ./ (k .* (T_opt - T_min))) + exp((Y(3) - T_opt) ./ (k .* (T_max - T_opt)))); %temperature rate change
            r_Fd = k_r_maxgm.*(F_d./(F_d + F_dhalf)); %feed rate change
            r_A_air = k_r_maxA.*(A_ai./(A_ai + k_A_half)); % airflow rate
            W_sub_percent = W./(W + Y(6)); %substrate moisture in [%]
            r_W = k_r_maxW./(1 + exp(-17.684.*W_sub_percent + 7.0622)); % moisture correction
            phi_Q_le = -k_air.*rho_air.*k_mle.*Y(2); %heat (Temperature) looses due to leakage and air fans
            k_ha = k_air.*rho_air + k_wv.*Y(4); % Total heat capacity in the production air
            k_hs = k_sub.*Y(6) + k_w.*W_c; %total heat capacity in thye substrate
            phi_H_le = k_le.*Y(4); % humidity looses due to leakage
            phi_Q_m_a = h_m_a.*A_m.*(Y(2) - Y(3)); % convective heat flux air./growing medium
            phi_C_le = -k_le.*Y(5); %Carbone dioxy d looses due to leakage
            e_s = 6.125.*exp(17.27.*Y(2)./(Y(2) + 243.03)); % magnus equation for evaporation
            e_a = (Y(4)./42).*e_s;
            phi_W_L_S = sigma.*A_m.*(e_a - e_s) ; %Evaporation condensation from surface of growing meduim
            phi_Q_L_S = k_w.*((100 - Y(3)) + k_hl).*phi_W_L_S;
            D_sum = (r_Tsub./k_r_maxT).*(r_Fd./k_r_maxgm).*(r_A_air./k_r_maxA).*(r_W./k_r_maxW).*k_dev;% development time of BFSL
            if (k_D_1 < D_sum && D_sum < k_D_3)%maturity activation function
                r_B_ma = 1;
            else
                r_B_ma = 0;
            end
            r_main = (r_Tsub./k_r_maxT).*(r_Fd./k_r_maxgm).*(r_A_air./k_r_maxA); %maintenance rate change with respect to environmental variables
            r_ma = r_B_ma.*(r_Tsub./k_r_maxT).*(r_Fd./k_r_maxgm).*(r_A_air./k_r_maxA);%maturity rate change with respect to environmental variables
            if (D_sum < k_D_1)
                r_B_ass = (1 - Y(1)./B_max); %assimilation function with respect to actuel larave size
            elseif (k_D_1 < D_sum && D_sum < k_D_2)
                r_B_ass = (1 - Y(1)./B_max).*(D_sum - k_D_2)./(k_D_1 - k_D_2);
            else
                r_B_ass = 0;
            end
            r_ass = r_B_ass.*(r_Tsub./k_r_maxT).*(r_Fd./k_r_maxgm).*(r_A_air./k_r_maxA).*(r_W./k_r_maxW);
            phi_Q_bio = n_lar.*(k_Qass.*k_alpha_ass.*r_ass.*k_ing + k_Qmai.*r_main.*k_main + k_Qma.*r_ma.*k_ma).*Y(1);%[j s-1]heat produced in the substrate due to the larval+microbiome activity
            phi_C_bio = ((r_C_1.*k_alpha_ass.*r_ass.*k_ing) + (r_C_2.*r_main.*k_main));%Co2 produced due to metabolism activity of larave
            
            %System of first order ODEs
            dYdt = zeros(8,1);
            dYdt(1) = (1 - k_alpha_ex - k_alpha_ass).*r_ass.*k_ing.*Y(1) - (r_main.*k_main + r_ma.*k_ma ).*Y(1);
            dYdt(2) = (phi_Q_le + phi_Q_m_a)./k_ha;
            dYdt(3) = (phi_Q_bio - phi_Q_m_a - phi_Q_L_S)./k_hs;
            dYdt(4) = (- phi_H_le + phi_W_L_S)./V_i;
            dYdt(5) = (phi_C_le + phi_C_bio)./V_i;
            dYdt(6) = theta.*(-(N_B./N0).*(r_ass.*k_ing.*Y(1)) + (N_F./N0).*(k_alpha_ex.*k_ing.*Y(1))+x.*phi_C_bio);
            dYdt(7) = omega.*(-(P_B./P0).*(r_ass.*k_ing.*Y(1)) + (P_F./P0).*(k_alpha_ex.*k_ing.*Y(1)));
            dYdt(8) = delta.*(-(K_B./K0).*(r_ass.*k_ing.*Y(1)) + (K_F./K0).*(k_alpha_ex.*k_ing.*Y(1)));
       end
     sol = Y;
    end
end
