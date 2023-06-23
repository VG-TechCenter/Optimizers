



%%% Designed and Developed by Pavel Trojovský and Mohammad Dehghani %%%


function[Best_score,Best_pos,POA_curve]=POA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)

lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

%% INITIALIZATION
for i=1:dimension
    X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
end

for i =1:SearchAgents
    L=X(i,:);
    fit(i)=fitness(L);
end
%%

for t=1:Max_iterations
    %% update the best condidate solution
    [best , location]=min(fit);
    if t==1
        Xbest=X(location,:);                                           % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        Xbest=X(location,:);
    end
    
    %% UPDATE location of food
    
    X_FOOD=[];
    k=randperm(SearchAgents,1);
    X_FOOD=X(k,:);
    F_FOOD=fit(k);
    
    %%
    for i=1:SearchAgents
        
        %% PHASE 1: Moving towards prey (exploration phase)
        I=round(1+rand(1,1));
        if fit(i)> F_FOOD
            X_new=X(i,:)+ rand(1,1).*(X_FOOD-I.* X(i,:)); %Eq(4)
        else
            X_new=X(i,:)+ rand(1,1).*(X(i,:)-1.*X_FOOD); %Eq(4)
        end
        X_new= max(X_new,lowerbound);X_new = min(X_new,upperbound);
        
        % Updating X_i using (5)
        f_new = fitness(X_new);
        if f_new <= fit (i)
            X(i,:) = X_new;
            fit (i)=f_new;
        end
        %% END PHASE 1: Moving towards prey (exploration phase)
        
        %% PHASE 2: Winging on the water surface (exploitation phase)
        X_new=X(i,:)+0.2*(1-t/Max_iterations).*(2*rand(1,dimension)-1).*X(i,:);% Eq(6)
        X_new= max(X_new,lowerbound);X_new = min(X_new,upperbound);
        
        % Updating X_i using (7)
        f_new = fitness(X_new);
        if f_new <= fit (i)
            X(i,:) = X_new;
            fit (i)=f_new;
        end
        %% END PHASE 2: Winging on the water surface (exploitation phase)
    end

    best_so_far(t)=fbest;
    average(t) = mean (fit);
    
end
Best_score=fbest;
Best_pos=Xbest;
POA_curve=best_so_far;
end
%%

