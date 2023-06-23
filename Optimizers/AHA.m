function [BestF,BestX,HisBestFit,VisitTable]=AHA(nPop,MaxIt,Low,Up,Dim,BenFunctions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FunIndex: The index of function.                    %
% MaxIt: The maximum number of iterations.            %
% nPop: The size of hummingbird population.           %
% PopPos: The position of population.                 %
% PopFit: The fitness of population.                  %
% Dim: The dimensionality of prloblem.                %
% BestX: The best solution found so far.              %
% BestF: The best fitness corresponding to BestX.     %
% HisBestFit: History best fitness over iterations.   %
% Low: The low boundary of search space               %
% Up: The up boundary of search space.                %
% VisitTable: The visit table.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [Low,Up,Dim]=FunRange(FunIndex);
    PopPos=zeros(nPop,Dim);
    PopFit=zeros(1,nPop);
    for i=1:nPop
        PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
        PopFit(i)=BenFunctions(PopPos(i,:));
    end

    BestF=inf;
    BestX=[];

    for i=1:nPop
        if PopFit(i)<=BestF
            BestF=PopFit(i);
            BestX=PopPos(i,:);
        end
    end

    % Initialize visit table
    HisBestFit=zeros(MaxIt,1);
    VisitTable=zeros(nPop) ;
    VisitTable(logical(eye(nPop)))=NaN;    
    
    for It=1:MaxIt
        DirectVector=zeros(nPop,Dim);% Direction vector/matrix

        for i=1:nPop
            r=rand;
            if r<1/3     % Diagonal flight
                RandDim=randperm(Dim);
                if Dim>=3
                    RandNum=ceil(rand*(Dim-2)+1);
                else
                    RandNum=ceil(rand*(Dim-1)+1);
                end
                DirectVector(i,RandDim(1:RandNum))=1;
            else
                if r>2/3  % Omnidirectional flight
                    DirectVector(i,:)=1;
                else  % Axial flight
                    RandNum=ceil(rand*Dim);
                    DirectVector(i,RandNum)=1;
                end
            end

            if rand<0.5   % Guided foraging
                [MaxUnvisitedTime,TargetFoodIndex]=max(VisitTable(i,:));
                MUT_Index=find(VisitTable(i,:)==MaxUnvisitedTime);
                if length(MUT_Index)>1
                    [~,Ind]= min(PopFit(MUT_Index));
                    TargetFoodIndex=MUT_Index(Ind);
                end

                newPopPos=PopPos(TargetFoodIndex,:)+randn*DirectVector(i,:).*...
                    (PopPos(i,:)-PopPos(TargetFoodIndex,:));
                newPopPos=SpaceBound(newPopPos,Up,Low);
                newPopFit=BenFunctions(newPopPos);
                
                if newPopFit<PopFit(i)
                    PopFit(i)=newPopFit;
                    PopPos(i,:)=newPopPos;
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(i,TargetFoodIndex)=0;
                    VisitTable(:,i)=max(VisitTable,[],2)+1;
                    VisitTable(i,i)=NaN;
                else
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(i,TargetFoodIndex)=0;
                end
            else    % Territorial foraging
                newPopPos= PopPos(i,:)+randn*DirectVector(i,:).*PopPos(i,:);
                newPopPos=SpaceBound(newPopPos,Up,Low);
                newPopFit=BenFunctions(newPopPos);
                if newPopFit<PopFit(i)
                    PopFit(i)=newPopFit;
                    PopPos(i,:)=newPopPos;
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(:,i)=max(VisitTable,[],2)+1;
                    VisitTable(i,i)=NaN;
                else
                    VisitTable(i,:)=VisitTable(i,:)+1;
                end
            end
        end

        if mod(It,2*nPop)==0 % Migration foraging
            [~, MigrationIndex]=max(PopFit);
            PopPos(MigrationIndex,:) =rand(1,Dim).*(Up-Low)+Low;
            PopFit(MigrationIndex)=BenFunctions(PopPos(MigrationIndex,:));
            VisitTable(MigrationIndex,:)=VisitTable(MigrationIndex,:)+1;
            VisitTable(:,MigrationIndex)=max(VisitTable,[],2)+1;
            VisitTable(MigrationIndex,MigrationIndex)=NaN;            
        end

        for i=1:nPop
            if PopFit(i)<BestF
                BestF=PopFit(i);
                BestX=PopPos(i,:);
            end
        end

        HisBestFit(It)=BestF;
    end
end
%%
function  X=SpaceBound(X,Up,Low)


    Dim=length(X);
    S=(X>Up)+(X<Low);    
    X=(rand(1,Dim).*(Up-Low)+Low).*S+X.*(~S);
end
