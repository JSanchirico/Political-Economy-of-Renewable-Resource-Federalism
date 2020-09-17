%  Main program file for Sanchirico et al. Political-Economy of
%  Renewable Resource Federalism,  Ecological Applications,  2020. Please
%  cite the main article if you use the code/algorithms for the open-loop
%  Nash equilibrium.

%  Please note: Code is developed to generate a solution in the paper for each management 
%  at a base case set of parameters. This code does not
%  recreate each figure in the paper or do the sensitivity analysis. 
%  Information on how to recreate the runs is available in the paper and SI.
%  To run this code, you need to download software from Tomlab tomopt.com/tomlab/
%  (propt and main solver package). While this is a licensed product, you are able to download a free trial.


% Contact Jim at jsanchirico@ucdavis.edu for any questions. 

clear all
close all

% Ecological parameters
r = [.85 .85]; % Growth rates
K = [1 1]; % Carrying capacity
d = .25*min(r); % Dispersal rate
S_0=[.20*K(1) .20*K(2)]; % Initial conditions

% Economic parameters
P =[10 10];  % Price
delta_=0.04;  % Discount rate

% The cii terms are calibrated so that with linear costs the open-access
% biomass is alph_% of K and the quadratic costs are 25% of the linear costs.
alph_=.20; fract=.25;
c = [alph_*P(1)*K(1) fract*alph_*P(1)*K(1);
     alph_*P(2)*K(2) fract*alph_*P(2)*K(2)];

% Tomlab parameters 
% Initial tomlab
T = 50; %  final T
Nset=80; % # of collocation points

%  The period over which to calculate the NPV
TCounter=30; % Note this when the solution is near the turnpike (steady-state values)
SSTime=30; % Time to pull out steady-state values

% creating the parameter set for sensitivity analysis
SenRuns=22; % grid used for sensitivity (not included in current code)
sens=1; % Indicator for case

% Calls a function file that develops the models across the cases
DISPERSAL=1; RR=1;DIS=1; % these were used to index all of the senstivity runs. Turned off for now.
x0=[]; RES=[];SOL=[]; % These are not used in the call
R=r; C=c; % R and C are useful when running sensitivity analysis

% Please note at the current set-up there is no heterogeneity in the
% parameters and so the 1st best=uniform effort=uniform shadow price. You
% need to introduce heterogeneity above for different results. Note that
% the solution can become unstable for very large differences (what large
% is depends on the level of the base paramters). 

% CASE=1; 1st best 
% CASE=2; Nash
% CASE=3: Uniform Effort
% CASE=5; Uniform shadow price

for CASE=[1 2 3 5] % Case 4 not utilized in paper
[solution,T1,S1,S2,h1,h2,rENT1,rENT2,nPV1,nPV2,x1s,x2s,Result,nPV1n,nPV2n] = ...
                Federalism_Final(T,Nset,x0,K,R,DISPERSAL,d,C,P,CASE,delta_,RR,RES,[],S_0,TCounter,SSTime);

            % Storing results of the different runs
             t1{DISPERSAL,CASE,RR,DIS}=T1;
             s1{DISPERSAL,CASE,RR,DIS}=S1;
             s2{DISPERSAL,CASE,RR,DIS}=S2;
             H1{DISPERSAL,CASE,RR,DIS}=h1;
                    H2{DISPERSAL,CASE,RR,DIS}=h2;
                    NPV1{DISPERSAL,CASE,RR,DIS} = nPV1;
                    NPV2{DISPERSAL,CASE,RR,DIS} = nPV2;
                    NPV1n{DISPERSAL,CASE,RR,DIS} = nPV1n;
                    NPV2n{DISPERSAL,CASE,RR,DIS} = nPV2n;
                    rent1{DISPERSAL,CASE,RR,DIS} = rENT1;
                    rent2{DISPERSAL,CASE,RR,DIS} = rENT2;
                    X1s{DISPERSAL,CASE,RR,DIS} =x1s;
                    X2s{DISPERSAL,CASE,RR,DIS} =x2s;
                    RESult{DISPERSAL,CASE,RR,DIS} = Result;
            
end            


%% Switching from one management regime to another
% This code is the anslysis of when you switch from one type of management effort to another. 
% Specifically, this is between uniform effort and game. To switch to the
% shadow price and game, you would change the case from 3 to 5

SCENARIO=1; % Way to store results when running the different cases
RR=1; % Used in the sensitivity analysis as a running parameter that called in different set of parameters
    
    
for CASE=[2 3]  %1=optimal, 2=game, 3=effort, 5=M rent

% This picks the correct counterfactual
    if CASE==2
        cASE=3;
    else
        cASE=2;
    end

    % These initial conditions are the optimal "steady-state" solutions
    % from the above analysis
    S_0=[X1s{DISPERSAL,cASE,RR,DIS} X2s{DISPERSAL,cASE,RR,DIS}];

    [solution,T1,S1,S2,h1,h2,rENT1,rENT2,nPV1,nPV2,x1s,x2s,Result] = ...
        Federalism_Final(T,Nset,x0,K,R,DISPERSAL,d,C,P,CASE,delta_,RR,RES,SOL,S_0,TCounter,SSTime);

    nt1{SCENARIO,DISPERSAL,CASE,RR,DIS}=T1;
    ns1{SCENARIO,DISPERSAL,CASE,RR,DIS}=S1;
    ns2{SCENARIO,DISPERSAL,CASE,RR,DIS}=S2;
    nH1{SCENARIO,DISPERSAL,CASE,RR,DIS}=h1;
    nH2{SCENARIO,DISPERSAL,CASE,RR,DIS}=h2;

    nrent1{SCENARIO,DISPERSAL,CASE,RR,DIS} = rENT1;
    nrent2{SCENARIO,DISPERSAL,CASE,RR,DIS} = rENT2;


end



