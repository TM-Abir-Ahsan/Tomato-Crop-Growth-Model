clear all; clc;

Data=xlsread('TOMGRO_INPUT.xlsx','Sheet1');

Thin=Data(:,3);
PPFD_in=Data(:,2);

% Initial Variables %Jones (1999)
N = 6.0; 
LAI = 0.006; % [-]
W = 0.28;
Wm = 0.0;
Wf = 0.0;
CO2 = 450; % ppm

% Growth data per day
N_hist = [];
LAI_hist = [];
W_hist = [];
Wm_hist = [];
Wf_hist = [];

% Growth rate per day
delN = [];
delLAI = [];
delW = [];
delWm = [];
delWf = [];

Tday = [];


for i = 1:1:250
    % Reset variables
    dNdt_ = 0;
    Td = 0; % Change!
    Tdaytime = 0;
    ndaytime = 0;
    PPFDd = 0;

    % for daily temperature
    for h = 1:1:24
        Td = Td + Thin(i*h);
        PPFDd = PPFDd + PPFD_in(i*h);
        Tdaytime = Tdaytime + (PPFD_in(i*h) ~=0)*Thin(i*h);
        ndaytime = ndaytime + (PPFD_in(i*h) ~=0);
        fN_ = fN(Thin(i*h));
        dNdt_ = dNdt_ + dNdt(fN_);
    end
    Td = Td / 24; % average daily temperature
    PPFDd = PPFDd / 24; % average daily temperature
    Tdaytime = Tdaytime / ndaytime;

%    fN_ = fN(Td);
%    dNdt_ = dNdt(fN_);

    % d(LAI)/dt
    lambda_ = lambdas(Td);
    dLAIdt_ = dLAIdt(LAI, 3.10, N, lambda_, dNdt_);

    % dWfdt
    fR_ = fR(N);
    LFmax_ = LFmax(CO2);
    PGRED_ = PGRED(Td);
    Pg_ = Pg(LFmax_, PGRED_, PPFDd, LAI);
    Rm_ = Rm(Td, W, Wm);
    GRnet_ = GRnet(Pg_, Rm_, fR_);
    fF_ = fF(Td);
    g_ = g(Tdaytime);
    dWfdt_ = dWfdt(GRnet_, fF_, N, g_);

    % dWdt
    dWdt_ = dWdt(LAI, dWfdt_, GRnet_, 3.10, dNdt_);

    % dWmdt
    Df_ = Df(Td);
    dWmdt_ = dWmdt(Df_, Wf, Wm, N);

    % Reload
    N = N + dNdt_;
    LAI = LAI + dLAIdt_;
    Wf = Wf + dWfdt_;
    W = W + dWdt_;
    Wm = Wm + dWmdt_;

    % Save
    N_hist = [N_hist, N];
    LAI_hist = [LAI_hist, LAI];
    W_hist = [W_hist, W];
    Wf_hist = [Wf_hist, Wf];
    Wm_hist = [Wm_hist, Wm];
    delN = [delN, dNdt_];
    delLAI = [delLAI, dLAIdt_];
    delW = [delW, dWdt_];
    delWf = [delWf, dWfdt_];
    delWm = [delWm, dWmdt_];
    Tday = [Tday, Tdaytime];
end

figure
plot(1:1:250,LAI_hist)
hold on
figure
plot(1:1:250,Wm_hist,1:1:250,Wf_hist,1:1:250,W_hist)
legend('Mature fruit dry weight','Fruit dry weight', 'Aboveground plant weight')
xlabel('Days after transplanation')
ylabel('Dry Weight [gm^-2]')

function result = fN(T)
    % Heuvelink(1994) & Jones(1991)
    if (T > 12 && T <= 28)
        result = 1.0 + 0.0281 * (T - 28);
    elseif (T > 28 && T < 50)
        result = 1.0 - 0.0455 * (T - 28);
    else
        result = 0;
    end
end

function result = dNdt(fN_)
    Nm = 0.5/24; % P Jones(1991) converted from per day to per hour and later integrated on a daily basis
    result = Nm * fN_;
end

function result = lambdas(Td)
    result = 1.0;
end

function result = dLAIdt(LAI, dens, N, lambda_, dNdt_)
    sigma = 0.038; % P Maximum leaf area expansion per node, coefficient in expolinear equation; Jones(1999)
    beta = 0.169; % P Coefficient in expolinear equation; Jones(1999)
    Nb = 16.0; % P Coefficient in expolinear equation, projection of linear segment of LAI vs N to horizontal axis; Jones(1999)
    LAImax = 4.0; % P Jones(1999)

    if (LAI > LAImax)
        result = 0;
    else
        a = exp(beta * (N - Nb));
        result = dens * sigma * lambda_ * a * dNdt_ / (1 + a);
    end
end

function result = dWdt(LAI, dWfdt_, GRnet_, dens, dNdt_)
    LAImax = 4.0; % P Jones(1999)
    if (LAI >= LAImax)
        p1 = 2.0; % P Jones(1999)
    else
        p1 = 0;
    end
    Vmax = 8.0; % P Jones(1999)

    a = dWfdt_ + (Vmax - p1) * dens * dNdt_;
    b = GRnet_ - p1 * dens * dNdt_;
    result = min(a, b);
end

function result = Df(T)
    % The rate of development or aging of fruit at temperature T; Jones(1991)
    if (T > 9 && T <= 28)
        result = 0.0017 * T - 0.015;
    elseif (T > 28 && T <= 35)
        result = 0.032;
    else
        result = 0;
    end
end

function result = dWmdt(Df_, Wf, Wm, N)
    NFF = 22.0; % P Jones(1999)
    kF = 5.0; % P Jones(1999)

    if (N <= NFF + kF)
        result = 0;
    else
        result = Df_ * (Wf - Wm);
    end
end

function result = fR(N)
    % root phenology-dependent fraction; Jones(1991)
    if (N >= 30)
        result = 0.07;
    else
        result = -0.0046 * N + 0.2034;
    end
end

function result = LFmax(CO2)
    % maximum leaf photosynthetic rate; Jones(1991)
    % CO2[ppm]?
    tau = 0.0693; % P carbon dioxide use efficiency; Jones(1991)
    result = tau * CO2;
end

function result = PGRED(T)
    % function to modify Pg under suboptimal daytime temperatures; Jones(1991)
    if (T > 0 && T <= 12)
        result = 1.0 / 12.0 * T;
    elseif (T > 12 && T < 35)
        result = 1.0;
    else
        result = 0;
    end
end

function result = Pg(LFmax_, PGRED_, PPFD, LAI)
    D = 2.593; % P coefficient to convert Pg from CO2 to CH2O; Jones(1991)
    K = 0.58; % P light extinction coefficient; Jones(1991)
    m = 0.1; % P leaf light transmission coefficient; Jones(1991)
    Qe = 0.0645; % P leaf quantum efficiency; Jones(1991)

    a = D * LFmax_ * PGRED_ / K;
    b = log(((1-m) * LFmax_ + Qe * K * PPFD) / ((1-m) * LFmax_ + Qe * K * PPFD * exp(-1 * K * LAI)));
    result = a * b;
end

function result = Rm(T, W, Wm)
    % Jones(1999)
    Q10 = 1.4; % P Jones(1991)
    rm = 0.016; % P Jones(1999)
    result = Q10 ^ ((T-20)/10) * rm * (W - Wm);
end

function result = GRnet(Pg_, Rm_, fR_)
    E = 0.717; % P convert efficiency; Dimokas(2009)
    result = max(0, E * (Pg_ - Rm_) * (1 - fR_));
end

function result = fF(Td)
    % Jones(1991)
    if (Td > 8 && Td <= 28)
        result = 0.0017 * Td - 0.0147;
    elseif (Td > 28)
        result = 0.032;
    else
        result = 0;
    end
end

function result = g(T_daytime)
    T_CRIT = 24.4; % P mean daytime temperature above which fruits abortion start; Jones(1999)
    if (T_daytime <= T_CRIT)
        result = 0;
    else
        result = max(0, 1.0 - 0.154 * (T_daytime - T_CRIT));
    end
end

function result = dWfdt(GRnet_, fF_, N, g_)
    NFF = 22; % 22.0 P Nodes per plant when the first fruit appears; Jones(1999)
    alpha_F = 0.80; % P Maximum partitioning of new growth to fruit; Jones(1999)
    v = 0.135; % P Transition coefficient between vegetative and full fruit growth; Jones(1999)
    fF_ = 0.5; % P ORIGINAL
    if (N <= NFF)
        result = 0;
    else
        result = GRnet_ * alpha_F * fF_ * (1 - exp(v * (NFF - N))) * g_;
    end
end

