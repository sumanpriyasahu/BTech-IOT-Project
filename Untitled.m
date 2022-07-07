close all;
clear;
clc;

xm=100;
ym=100;

sink.x=0.5*xm;
sink.y=0.5*ym;

n=100;

p=0.1;
Eo=0.5;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10*0.00000000001;
Emp=0.0013*0.00000000001;
EDA=5*0.000000001;

m=0.1;
a=1;
rmax=4000;

do=sqrt(Efs/Emp);
%figure(1);

for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    S(i).type='N';
    temp_rnd0=i;
    
    if(temp_rnd0>=m*n+1)
        S(i).E=Eo;
        S(i).ENERGY=0;
       % plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    
     if(temp_rnd0<m*n+1)
        S(i).E=Eo;
        S(i).ENERGY=1;
       % plot(S(i).xd,S(i).yd,'+');
        hold on;
     end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%plot(S(n+1).xd,S(n+1).yd,'x');
%figure(1);

countCHs=0;
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax
    r
if(mod(r,round(1/p))==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
end
hold off;

dead=0;
dead_a=0;
dead_n=0;

packets_TO_BS=0;
packets_TO_CH=0;

PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;
%figure(1);

for i=1:1:n
    if(S(i).E<=0)
        %plot(S(i).xd,S(i).yd,'red.');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;
    end
    
    if(S(i).E>0)
        S(i).type='N';
        if(S(i).ENERGY==0)
           % plot(S(i).xd,S(i).yd,'o');
        end
        if(S(i).ENERGY==1)
           % plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
    
    tototalenergy(i)=S(i).E;
end
tototal_energy=sum(tototalenergy);
STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;

if r+1==1095
    S(i).E=0;
end
if(dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
    if(S(i).E>0)
        temp_rand=rand;
        if((S(i).G)<=0)
        if(temp_rand<=(p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            %plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt((S(i).xd-(S(n+1).xd))^2+(S(i).yd-(S(n+1).yd))^2);
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            distance;
            
            if(distance>do)
                S(i).E-((ETX+EDA)*(4000)+Emp*4000*(distance*distance*distance*distance));
            end
            if(distance<=do)
                S(i).E-((ETX+EDA)*(4000)+Efs*4000*(distance*distance));
            end
        end
    end
end
end

STATISTICS(r+1).CLUSTERHERADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

for i=1:1:n
    if(S(i).type=='N' && S(i).E>0)
        if(cluster-1>=1)
            min_dis=sqrt((S(i).xd-S(n+1).xd)^2+(S(i).yd-S(n+1).yd)^2);
            min_dis_cluster=1;
            for c=1:1:cluster-1
                temp=min(min_dis,sqrt((S(i).xd-C(c).xd)^2+(S(i).yd-C(c).yd)^2));
                
                if(temp<min_dis)
                    min_dis=temp;
                    min_dis_cluster=c;
                end
            end
            min_dis;
            if(min_dis>do)
                S(i).E=S(i).E-(ETX*(4000)+Emp*4000*(min_dis*min_dis*min_dis*min_dis));
            end
            if(min_dis<=do)
                S(i).E=S(i).E-(ETX*4000)+Efs*4000*(min_dis*min_dis);
            end
            if(min_dis>0)
                S(C(min_dis_cluster).id).E=S(C(min_dis_cluster).id).E-((ERX+EDA)*4000);
                PACKETS_TO_CH(r+1)=n-dead-cluster+1;
            end
            S(i).min_dis=min_dis;
            S(i).min_dis_cluster=min_dis_cluster;
        end
    end
end
hold on;
countCHs;
rcountCHs=countCHs+countCHs;

end

AlliveLeach=100-DEAD;
plot(1:4001,AlliveLeach,'-b')


                
            
            
            
    
    
    
    