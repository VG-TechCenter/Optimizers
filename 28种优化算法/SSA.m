

function [fMin , bestX,Convergence_curve ] = SSA(pop, M,c,d,dim,fobj  )
        
 P_percent = 0.2;    % The population size of producers accounts for "P_percent" percent of the total population size       
   %生产者占所有种群的0.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNum = round( pop *  P_percent );    % The population size of the producers    
%生产者数量取整

lb= c.*ones( 1,dim );    % Lower limit/bounds/     a vector    约束上限
ub= d.*ones( 1,dim );    % Upper limit/bounds/     a vector  约束下限
%Initialization
for i = 1 : pop    
    x( i, : ) = lb + (ub - lb) .* rand( 1, dim );   %随机初始化n个种群
    fit( i ) = fobj( x( i, : ) ) ;               %计算所有群体的适应情况，如果求最小值的,越小代表适应度越好              
end
% 以下找到最小值对应的麻雀群
pFit = fit;                      
pX = x;                            % The individual's best position corresponding to the pFit
[ fMin, bestI ] = min( fit );      % fMin denotes the global optimum fitness value
bestX = x( bestI, : );             % bestX denotes the global optimum position corresponding to fMin
 
 % Start updating the solutions.
for t = 1 : M    
      
  [ ans, sortIndex ] = sort( pFit );% Sort.
     
  [fmax,B]=max( pFit );
   worse= x(B,:);         %找到最差的个体      

   r2=rand(1);      %产生随机数    感觉没有啥科学依据，就是随机数
   %大概意思就是在0.8概率内原来种群乘一个小于1的数，种群整体数值缩小了
   %大概意思就是在0.2概率内原来种群乘一个小于1的数，种群整体数值+1
   %变化后种群的数值还是要限制在约束里面
   %对前pNum适应度最好的进行变化 ，即生产者进行变化，可见生产者是挑最好的
if(r2<0.8)
    for i = 1 : pNum                                                        % Equation (3)
         r1=rand(1);
        x( sortIndex( i ), : ) = pX( sortIndex( i ), : )*exp(-(i)/(r1*M)); %将种群按适应度排序后更新 
        x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );  %将种群限制在约束范围内
        fit( sortIndex( i ) ) = fobj( x( sortIndex( i ), : ) );   
    end
  else
  for i = 1 : pNum            
  x( sortIndex( i ), : ) = pX( sortIndex( i ), : )+randn(1)*ones(1,dim);
  x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );
  fit( sortIndex( i ) ) = fobj( x( sortIndex( i ), : ) );       
  end      
end
 %把经过变化后最好的种群记录下来  
 [ fMMin, bestII ] = min( fit );      
  bestXX = x( bestII, : );            
  
 %下面是乞讨者 
   for i = ( pNum + 1 ) : pop                     % Equation (4)
         A=floor(rand(1,dim)*2)*2-1;           %产生1和-1的随机数
         
          if( i>(pop/2))
           %如果i>种群的一半，代表遍历到适应度靠后的一段，代表这些序列的种群可能在挨饿
           x( sortIndex(i ), : )=randn(1)*exp((worse-pX( sortIndex( i ), : ))/(i)^2);  
           %适应度不好的，即靠后的麻雀乘了一个大于1的数，向外拓展
          else
        x( sortIndex( i ), : )=bestXX+(abs(( pX( sortIndex( i ), : )-bestXX)))*(A'*(A*A')^(-1))*ones(1,dim);  
           %这是适应度出去介于生产者之后，又在种群的前半段的,去竞争生产者的食物，在前面最好种群的基础上
           %再进行变化一次，在原来的基础上减一些值或者加一些值
         end  
        x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );  %更新后种群的限制在变量范围
        fit( sortIndex( i ) ) = fobj( x( sortIndex( i ), : ) );                    %更新过后重新计算适应度
   end
   %在全部种群中找可以意识到危险的麻雀
  c=randperm(numel(sortIndex));
   b=sortIndex(c(1:20));
  for j =  1  : length(b)      % Equation (5)
    if( pFit( sortIndex( b(j) ) )>(fMin) )
         %如果适应度比最开始最小适应度差的话，就在原来的最好种群上增长一部分值
        x( sortIndex( b(j) ), : )=bestX+(randn(1,dim)).*(abs(( pX( sortIndex( b(j) ), : ) -bestX)));
        
        else
        %如果适应度达到开始最小的适应度值，就在原来的最好种群上随机增长或减小一部分
        x( sortIndex( b(j) ), : ) =pX( sortIndex( b(j) ), : )+(2*rand(1)-1)*(abs(pX( sortIndex( b(j) ), : )-worse))/ ( pFit( sortIndex( b(j) ) )-fmax+1e-50);

          end
        x( sortIndex(b(j) ), : ) = Bounds( x( sortIndex(b(j) ), : ), lb, ub );
       
       fit( sortIndex( b(j) ) ) = fobj( x( sortIndex( b(j) ), : ) );
 end
    for i = 1 : pop 
        %如果哪个种群适应度好了，就把变化的替换掉原来的种群
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin )   %最优值以及最优值位置看是否变化
           fMin= pFit( i );
            bestX = pX( i, : );
         
            
        end
    end
  
    Convergence_curve(t)=fMin;
  
end


% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
  % Apply the lower bound vector
  temp = s;
  I = temp < Lb;
  temp(I) = Lb(I);
  
  % Apply the upper bound vector 
  J = temp > Ub;
  temp(J) = Ub(J);
  % Update this new move 
  s = temp;

%---------------------------------------------------------------------------------------------------------------------------
