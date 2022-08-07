clc;
clear;
close all;
%%Initialize for DPSO
pcnt = 8 ; %粒子数目
maxLoop = 50 ; %最大迭代次数
gb = 1:maxLoop
w = 0.9 ;
c1 = 0.8 ;
c2 = 0.7 ;
n = 15 ;%车位数目
m = 10 ;%车辆数目
load("cost.mat")
%T = round(100* rand((m+n),2)) ;
x = T(:,1) ;
y = T(:,2) ;
plot(x(1:n),y(1:n),'b*',x(n+1:n+m),y(n+1:n+m),'go')
save('cost.mat',"T","x","y")
dis = ones(n + m, n + m); 
for i = 1:n + m
    for j = 1:n + m
        dis(i, j) = ((T(i, 1) - T(j, 1))^2 + (T(i, 2) - T(j, 2))^2)^0.5;
        %前n行n列为车位所在位置
    end
end
%%Initialize particles
p.X = zeros(1,m);
p.V = zeros(1,m);
p.pBestX = p.X ;
p.pBest = 0;
p = repmat(p,[1 pcnt]);
gb = p.pBestX ;
%X：粒子当前位置
%XL：粒子自身最佳位置
for k = 1:pcnt
    p(k).X = randperm(n,m);%从1-10中随机抽取三个数
    p(k).V = randperm(m);
    p(k).pBestX = p(k).X;
    p(k).pBest = fitness(p(k).X ,  m , n, dis );
end

gBest = fitness(p(1).X , m , n , dis );
for i = 1:pcnt
    if fitness(p(i).X ,m , n ,dis ) <= gBest 
        gBest = fitness(p(i).X , m , n , dis  );
        gBestX = p(i).X;
    end
end

%%Update particles
% Xi = w.*Vi-1 + c1.* r1.*(pBest-Xi-1) + c2.*r2.*(gBest - Xi-1)
for loop = 1 : maxLoop
    for i = 1 : pcnt
        r = rand ;
        p(i).V = w .*p(i).V + c1 .* r .* ( p(i).pBestX - p(i).X ) +c2 .* r .* ( gBestX - p(i).X );
%         p(i).V = p(i).pBest - p(i).X ;
%         p(i).V = gBestX - p(i).X ;
        if p(i).V > (n/2)
            p(i).V = (n/2) ;
        elseif p(i).V < -(n/2)
            p(i).V = -(n/2) ;
        end
        p(i).V = floor( p(i).V ) ;

        %粒子位置变化
        for j = 1 : m 
            p(i).X(j) = p(i).X(j) + p(i).V(j) ;
            p(i).X(j) = mod( p(i).X(j) , n );%取余
            if ismember(0 , p(i).X) %判断是否有0
                ind = find( 0 == p(i).X) ;
                p(i).X(ind) = n ;
            end
            k = j ;
            while k > 0
                %修改第j维
                repeat = length( p(i).X(1:j) )- length( unique( p(i).X(1:j) )) ;
                %判断是否被分配到同一车位
                while logical( repeat )
                    %根据粒子速度是否为正进行判断
                    if p(i).V(j) > 0
                        p(i).X(j) = p(i).X(j)  - repeat;
                    elseif p(i).V(j) < 0
                        p(i).X(j) = p(i).X(j)  + repeat;
                    end
                    repeat = repeat - 1 ;
                end                
                k = k - 1 ;
            end    
            if ismember(0 , p(i).X) %判断是否有0
                ind = find( 0 == p(i).X) ;
                p(i).X(ind) = n ;
            end
            if p(i).X(j) >= n
                p(i).X(j) = n ;            
            elseif p(i).X(j) <= 1
                p(i).X(j) = 1 ;
            end
        end

        %更新pBest与gBest
        fitnessValue = fitness( p(i).X , m , n , dis ) ;
        if fitnessValue < p(i).pBest
            p(i).pBest = fitnessValue;
            p(i).pBestX = p(i).X ;
        end
        if fitnessValue < gBest
            gBest = fitnessValue;
            gBestX = p(i).X ;
        end
        gb(loop) = gBest;
    end
    
end
x = 1:maxLoop
plot(x,gb)


function y = fitness(val,m,n,dis)
    gBest = 0;
    for i = 1 : m
        
        x = dis(val(i),i+n);
        gBest = x + gBest ;
    end
    y = gBest ;
end
