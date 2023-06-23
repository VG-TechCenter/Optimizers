%_______________________________________________________________________________________%
%  Dwarf Mongoose Optimization Algorithm source codes (version 1.0)                     %
%                                                                                       %
%  Developed in MATLAB R2015a (7.13)                                                    %
%  Author and programmer: Jeffrey O. Agushaka and Absalom E. Ezugwu and Laith Abualigah %
%         e-Mail:  EzugwuA@ukzn.ac.za                                                   %
%                                                                                       %
%   Main paper:                                                                         %
%__Dwarf Mongoose Optimization Algorithm: A new nature-inspired metaheuristic optimizer %
%__Main paper: please, cite it as follws:_______________________________________________%
%_______________________________________________________________________________________%
%_______________________________________________________________________________________%
function [BEF,BEP,BestCost]=DMOA(nPop,MaxIt,VarMin,VarMax,nVar,F_obj)



%nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

%VarMin=-10;         % Decision Variables Lower Bound
%VarMax= 10;         % Decision Variables Upper Bound

%% ABC Settings

% MaxIt=1000;              % Maximum Number of Iterations

% nPop=100;               % Population Size (Family Size)

nBabysitter= 3;         % Number of babysitters

nAlphaGroup=nPop-nBabysitter;         % Number of Alpha group

nScout=nAlphaGroup;         % Number of Scouts

L=round(0.6*nVar*nBabysitter); % Babysitter Exchange Parameter 

peep=2;             % Alpha female痴 vocalization 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty Mongoose Structure
empty_mongoose.Position=[];
empty_mongoose.Cost=[];

% Initialize Population Array
pop=repmat(empty_mongoose,nAlphaGroup,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;
tau=inf;
Iter=1;
sm=inf(nAlphaGroup,1);

% Create Initial Population
for i=1:nAlphaGroup
        for j=1:VarSize
        BB=unifrnd(VarMin(j),VarMax(j),VarSize);
        AA(j,:)=BB(1,:);
       end
    pop(i).Position=AA;
    pop(i).Cost=F_obj(pop(i).Position);
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nAlphaGroup,1);
CF=(1-Iter/MaxIt)^(2*Iter/MaxIt);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% DMOA Main Loop

for it=1:MaxIt
    
    % Alpha group
     F=zeros(nAlphaGroup,1);
     MeanCost = mean([pop.Cost]);
    for i=1:nAlphaGroup
        
        % Calculate Fitness Values and Selection of Alpha
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
        P=F/sum(F);
      % Foraging led by Alpha female
    for m=1:nAlphaGroup
        
        % Select Alpha female
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to Alpha
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);
        
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newpop.Cost=F_obj(newpop.Position);
        
        % Comparision
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
        else
            C(i)=C(i)+1;
        end
        
    end   
    
    % Scout group
    for i=1:nScout
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nAlphaGroup];
        k=K(randi([1 numel(K)]));
        
        % Define Vocalization Coeff.
        phi=(peep/2)*unifrnd(-1,+1,VarSize);
        
        % New Mongoose Position
        newpop.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newpop.Cost=F_obj(newpop.Position);
        
        % Sleeping mould
        sm(i)=(newpop.Cost-pop(i).Cost)/max(newpop.Cost,pop(i).Cost);
        
        % Comparision
        if newpop.Cost<=pop(i).Cost
            pop(i)=newpop;
        else
            C(i)=C(i)+1;
        end
        
    end    
    % Babysitters
    for i=1:nBabysitter
        if C(i)>=L
            for j=1:VarSize
            BB=unifrnd(VarMin(j),VarMax(j),VarSize);
            AA(j,:)=BB(1,:);
            end
%             pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            pop(i).Position=AA;
            pop(i).Cost=F_obj(pop(i).Position);
            C(i)=0;
        end
    end    
     % Update Best Solution Ever Found
    for i=1:nAlphaGroup
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end    
        
   % Next Mongoose Position
   newtau=mean(sm);
   for i=1:nScout
        M=(pop(i).Position.*sm(i))/pop(i).Position;
        if newtau>tau
           newpop.Position=pop(i).Position-CF*phi*rand.*(pop(i).Position-M);
        else
           newpop.Position=pop(i).Position+CF*phi*rand.*(pop(i).Position-M);
        end
        tau=newtau;
   end
       
   % Update Best Solution Ever Found
    for i=1:nAlphaGroup
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BEF=BestSol.Cost;
    BEP=BestSol.Position;
    % Display Iteration Information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    

end
end

%%

function i=RouletteWheelSelection(P)

    r=rand;
    
    C=cumsum(P);
    
    i=find(r<=C,1,'first');

end
%%
function z=Sphere(x)

    z=sum(x.^2);

end
%%

