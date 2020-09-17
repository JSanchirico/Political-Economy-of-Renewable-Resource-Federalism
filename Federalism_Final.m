% Function file for the different runs. Easier to use when running all of
% the analysis in the paper.

function [solution,Ti,S1s,S2s,h1s,h2s,RENT1,RENT2,NPV1,NPV2,X1s,X2s,Result,NPV1n,NPV2n] = ...
    Federalism_Final(T,Nset,x0,K,R,DISPERSAL,d,c,P,CASE,delta_,RR,RES,SOL,S_0,TCounter,SSTime)

% Tomlab variables
toms t
Phas=tomPhase('Phas',t,0,T,Nset,[],'Gauss'); % Guassian collocation points
setPhase(Phas)
tomStates S1 S2
tomControls h1 h2

% Initial guess
x0 = { icollocate({S1 == .5*K(1)
    S2 == .5*K(2)})
    collocate(h1==.25*R(1)*K(1))
    collocate(h2==.25*R(2)*K(2))};

% Non-negativity constraints
cbox = {0<= icollocate(S1)
    0<= icollocate(S2)
    0<= icollocate(h1)
    0<= icollocate(h2)
    };

cterm = final({}); % If we want to impose a T condition, we can add it here
cbnd = initial({S1==S_0(1), S2==S_0(2)}); % initial conditions

% Dispersal
d11 = -d; d12=d;
d22 = -d; d21=d;
D = [d11 d12 ; d21 d22]; % Dispersal matrix
D1 = D(1,:)*[S1; S2];
D2 = D(2,:)*[S1; S2];

% Cost and rent functions
C1= c(1,1)*h1  + c(1,2)*h1^2;
C2= c(2,1)*h2  + c(2,2)*h2^2;
Rent1 = P(1)*h1*S1 - C1;
Rent2 = P(2)*h2*S2 - C2;

if CASE~=2  % CASE 2 is the nash open-loop solution, this case utilizes code below
    
    if CASE==1  % first-best
        cbox2 = { };
        % State equations
        ceq = collocate({...
            dot(S1)==R(1)*S1*(1-S1/K(1)) + D1 - h1*S1;
            dot(S2)==R(2)*S2*(1-S2/K(2)) + D2 - h2*S2;
            });
        
        
    elseif CASE==3 % effort restriction
        cbox2= {collocate(h1)==collocate(h2)};
        
        % State equations
        ceq = collocate({...
            dot(S1)==R(1)*S1*(1-S1/K(1)) + D1 - h1*S1;
            dot(S2)==R(2)*S2*(1-S2/K(2)) + D2 - h2*S2;
            
            });
        
        
    elseif CASE==4 % Not in use in final version of the paper
        
    else  % Marginal rent restriction
        
        MRent1= (P(1)*S1 - c(1,1) - 2*c(1,2)*h1)./S1;
        MRent2= (P(2)*S2 - c(2,1) - 2*c(2,2)*h2)./S2;
        
        cbox2=collocate({(MRent1-MRent2==0)});
        
        % State equations
        ceq = collocate({...
            dot(S1)==R(1)*S1*(1-S1/K(1)) + D1 - h1*S1;
            dot(S2)==R(2)*S2*(1-S2/K(2)) + D2 - h2*S2;
            
            });
    end
    
    
    objective = -integrate(exp(-delta_*t).*(Rent1+Rent2));
    options=struct;
    options.solver = 'snopt'; % choose solver ('knitro' is other good one)
    [solution, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},x0,options);
    
    % If a run is having a hard time solving, you can go thru these runs
    % that change the initial guess and/or solver used
    counter=0;
    while result.ExitFlag~=0 && counter<3 % limits the loop size
        
        if counter==1
            % Generate new initial guess without constraint
            options.solver = 'snopt';
            [solution1, ~]=ezsolve(objective,{cbnd,ceq,cbox,cterm},x0,options);
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution1,options);
        elseif counter==2
            options.solver = 'knitro';
            [solution1, ~]=ezsolve(objective,{cbnd,ceq,cbox,cterm},x0,options);
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution1,options);
        else
            options.solver = 'knitro';
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution,options);
        end
        
        if result.ExitFlag~=0
            options.solver = 'npsol';
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution1,options);
        end
        
        if result.ExitFlag~=0
            options.solver = 'snopt';
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution1,options);
        end
        disp('Non-game model')
        counter=counter+1
        solution=solution1;
    end
    
    
    % If we have a solution, then we generate some values
    if result.ExitFlag==0
        Result=result.f_k;
        Ti=subs(icollocate(t),  solution);
        S1s=subs(icollocate(S1), solution);
        S2s=subs(icollocate(S2), solution);
        h1s=subs(icollocate(h1), solution);
        h2s=subs(icollocate(h2), solution);
        
        Ind=max(find(round(Ti)==TCounter));
        Time=Ti(1:Ind);
        
        
        RENT1 =  (P(1).*h1s(1:Ind).*S1s(1:Ind) -...
            ( c(1,1).*h1s(1:Ind)  + c(1,2).*h1s(1:Ind).^2));
        RENT2 =( P(2)*h2s(1:Ind).*S2s(1:Ind) -...
            ( c(2,1)*h2s(1:Ind)  + c(2,2)*h2s(1:Ind).^2));
        
        
        % This code fits a function to rents in each patch and uses the
        % function to calculate the net present value of the infinite
        % horizon problem. Otherwise the objective function from the solver
        % is only over the finite time period.
        f = fittype('pchipinterp');
        g1=fit(Time,exp(-delta_.*Time).*RENT1,f);
        int1=integrate(g1,0:1:Time(end),0);
        g2=fit(Time,exp(-delta_.*Time).*RENT2,f);
        int2=integrate(g2,0:1:Time(end),0);
        
        
        % This includes the scrap value of operating at the steady-state
        % value forever (discounted back to the present)
        NPV1n = int1(end)+exp(-delta_.*Time(end)).*RENT1(end)./delta_;
        NPV2n = int2(end)+exp(-delta_.*Time(end)).*RENT2(end)./delta_;
        
        % Another way to calculate the NPV
        NPV1 = trapz(exp(-delta_.*Time).*RENT1)+ exp(-delta_.*Time(end)).*RENT1(end)./delta_;
        NPV2 = trapz(exp(-delta_.*Time).*RENT2)+ exp(-delta_.*Time(end)).*RENT2(end)./delta_;
        
        X1s =S1s(max(find(round(Ti)==TCounter)));
        X2s =S2s(max(find(round(Ti)==TCounter)));
        
        
        
    else
        Result=-999999; % This code is called if the solver did not find a solution. It prohibits the code when running in many loops to put in values corresponding to run that did not converge.
    end
end

% Here is where the game theory code us us called
if CASE==2
    
    [solution,Ti,S1s,S2s,h1s,h2s,RENT1,RENT2,NPV1,NPV2,X1s,X2s,Result,NPV1n,NPV2n] =...
        federalism_game_inside(T,Nset,[],K,R,DISPERSAL,d,c,P,CASE,delta_,RR,RES,SOL,S_0,D,D1,D2,TCounter,SSTime);
end

end

function [solution,Ti,S1s,S2s,h1s,h2s,RENT1,RENT2,NPV1,NPV2,X1s,X2s,Result,NPV1n,NPV2n] =...
    federalism_game_inside(T,Nset,x0,K,R,DISPERSAL,d,c,P,CASE,delta_,RR,RES,SOL,S_0,D,D1,D2,TCounter,SSTime)

% All of the same calls to Tomlab variables
toms t
Phas=tomPhase('Phas',t,0,T,Nset,[],'Gauss'); % Guassian collocation points
setPhase(Phas)
tomStates S1 S2   
tomControls h1 h2

% Solving this part of the problem included some tricks, like using the
% solution from the 1st best as a starting guess for the game solution. I
% have left in some the code that was utilized to help the solver along.

if isempty(RES)
    x0 = { icollocate({S1 == .5*K(1)
        S2 == .5*K(2)})
        collocate(h1==.25*R(1)*K(1))
        collocate(h2==.25*R(2)*K(2))};
    % Setting the harvest fixed in patch 2 initially
    h2_init = collocate(h2==.25*R(2)*K(2));
    
else
    
    % takes initial guess from outside or inside the program depending on
    % the run
    Ti=subs(icollocate(t),  SOL);
    s1=subs(icollocate(S1), SOL);
    s2=subs(icollocate(S2), SOL);
    H1=subs(icollocate(h1), SOL);
    H2=subs(icollocate(h2), SOL);
    x0 = {
        icollocate({S1 == s1
        S2 == s2})
        collocate(h1==H1)
        collocate(h2==H2)};
    
    % Setting the harvest fixed in patch 2 initially
    h2_init = subs(collocate(h2),SOL);
    
end

cbox = {0<= icollocate(S1)
    0<= icollocate(S2)
    0<= icollocate(h1)
    0<= icollocate(h2)
    };

cterm = final({}); % If we want to impose a T condition, we can add it here
cbnd = initial({S1==S_0(1), S2==S_0(2)}); % initial conditions

% Cost and rent functions
C1= c(1,1)*h1  + c(1,2)*h1^2;
C2= c(2,1)*h2  + c(2,2)*h2^2;
Rent1 = P(1)*h1*S1 - C1;
Rent2 = P(2)*h2*S2 - C2;

iter=25;
options=struct;
ERR=ones(iter,1);
F=zeros(iter,2);


% Here is where we are running the game algorithm
for k=1:iter
    
    options.solver = 'snopt'; % choose solver ('knitro' is other good one)
    
    % Solving H1 given H2
    if k==1
        cbox2 = {collocate(h2)==(h2_init)};
    else
        cbox2 = {collocate(h2)==collocate(h2_init)};
    end
    
    ceq = collocate({...
        dot(S1)==R(1)*S1*(1-S1/K(1)) + D1 - h1*S1;
        dot(S2)==R(2)*S2*(1-S2/K(2)) + D2 - h2*S2;
        });
    
    objective = -integrate(exp(-delta_*t).*(Rent1));  
    [solution, gresult1]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},x0,options);
    GAMESOlution{k,1}=solution;
    
    % Below this helps when not finding a solution in the initial run
    counter=0;
    while gresult1.ExitFlag~=0 && counter<10 % limits the lopp size
        options.solver = 'knitro';
        if counter==1
            [solution, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},x0,options);
        else
            [solution, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution,options);
        end
        
        if gresult1.ExitFlag~=0
            options.solver = 'npsol';
            [solution, gresult1]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution,options);
        end
        
        if gresult1.ExitFlag~=0
            options.solver = 'snopt';
            [solution, gresult1]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution,options);
        end
        
        disp('game model')
        counter=counter+1
        options.solver
    end
    
    % Find the level of patch 1 harvest from the above solution to impose
    % on patch 2
    h1_init = subs(h1,solution);
   
    
    %% Patch 2's optimization taking the solution from patch 1's optimization
    %solving H2 given H1
    cbox2 = {collocate(h1)==collocate(h1_init)};
    
    ceq = collocate({...
        dot(S1)==R(1)*S1*(1-S1/K(1)) + D1 - h1*S1;
        dot(S2)==R(2)*S2*(1-S2/K(2)) + D2 - h2*S2;
        });
    
    options.solver = 'snopt';
    objective = -integrate(exp(-delta_*t).*(Rent2)); %-exp(-delta_.*T).*(final(Rent2)./delta_);
    [solution, gresult2]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},x0,options);
    GAMESOlution{k,2}=solution;
    
    counter=0;
    while gresult2.ExitFlag~=0 && counter<10 % limits the lopp size
        options.solver = 'knitro';
        if counter==1
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},x0,options);
        else
            [solution1, result]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution,options);
        end
        
        if gresult2.ExitFlag~=0
            options.solver = 'npsol';
            [solution, gresult2]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution1,options);
        end
        
        if gresult2.ExitFlag~=0
            options.solver = 'snopt';
            [solution, gresult2]=ezsolve(objective,{cbnd,ceq,cbox,cterm,cbox2},solution1,options);
        end
        
        disp('game model 2')
        CASE
        counter=counter+1
        options.solver
    end
    
    
    
    s1_init = subs(S1,solution);
    s2_init = subs(S2,solution);
    h1_init = subs(h1,solution);
    h2_init = subs(h2,solution);
    
    
    x0 = {
        icollocate({S1 == s1_init
        S2 == s2_init
        })
        collocate(h1==h1_init)
        collocate(h2==h2_init)
        };
    
    % here is where we are doing the convergence of the two optimizations
    if gresult1.ExitFlag==0 && gresult2.ExitFlag==0 % makes sure both are converging to count as solution
        F(k,1)=gresult1.f_k; % Value of the objective function in patch 1's optimization holding patch 2 harvest fixed
        F(k,2)=gresult2.f_k; % Value of the objective function in patch 2's optimization holding patch 1 harvest fixed
        
    else
        F(k,1)=-99999; % setting these as such should preculde exiting with runs that did not converge
        F(k,2)=-99999;
    end
    
    Result=F(k,1)+F(k,2);
    ERR(k)=abs(F(k,1))+abs(F(k,2));
    
    if k>2
        Error1 = (abs(F(k,1))-abs(F(k-1,1)))^2; % Calculating the difference in objective function across runs for each patch's optimization
        Error2 = (abs(F(k,2))-abs(F(k-1,2)))^2;
        Err2=Error1+Error2; 
    end
    
    if k>3 && Err2<1e-10 % Convergence criteria on sum of the squared differences 
        
        % with a solution we do the same calculations as above
        Ti=subs(icollocate(t),  solution);
        S1s=subs(icollocate(S1), solution);
        S2s=subs(icollocate(S2), solution);
        h1s=subs(icollocate(h1), solution);
        h2s=subs(icollocate(h2), solution);
        
        
        Ind=max(find(round(Ti)==TCounter));
        Time=Ti(1:Ind);
        
        
        RENT1 =  (P(1).*h1s(1:Ind).*S1s(1:Ind) -...
            ( c(1,1).*h1s(1:Ind)  + c(1,2).*h1s(1:Ind).^2));
        RENT2 =( P(2)*h2s(1:Ind).*S2s(1:Ind) -...
            ( c(2,1)*h2s(1:Ind)  + c(2,2)*h2s(1:Ind).^2));
        
       
        
        f = fittype('pchipinterp');
        g1=fit(Time,exp(-delta_.*Time).*RENT1,f);
        int1=integrate(g1,0:1:Time(end),0);
        g2=fit(Time,exp(-delta_.*Time).*RENT2,f);
        int2=integrate(g2,0:1:Time(end),0);
        NPV1n = int1(end)+exp(-delta_.*Time(end)).*RENT1(end)./delta_;
        NPV2n = int2(end)+exp(-delta_.*Time(end)).*RENT2(end)./delta_;
        
        NPV1 = trapz(exp(-delta_.*Time).*RENT1)+ exp(-delta_.*Time(end)).*RENT1(end)./delta_;
        NPV2 = trapz(exp(-delta_.*Time).*RENT2)+ exp(-delta_.*Time(end)).*RENT2(end)./delta_;
        
        X1s =S1s(max(find(round(Ti)==TCounter)));
        X2s =S2s(max(find(round(Ti)==TCounter)));
        
 
        
        ITER=k;
        ERRj=Err2;
        
        break;
    end
    
    if k==iter && Err2>1e-9 % Report nothing when solver does not work, which will create an error later in the code
        
        Ti=[];
        S1s=[];
        S2s=[];
        h1s=[];
        h2s=[];
        X1s =[];
        X2s =[];
        RENT1 =  [];
        RENT2 = [];
        
        NPV1 = [];
        NPV2 = [];
        NPV1n = [];
        NPV2n = [];
        rent1 = [];
        rent2 = [];
    end
    
end
end


