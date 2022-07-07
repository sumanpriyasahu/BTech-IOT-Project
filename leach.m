clear;

xm=100;
ym=100;
sink.x=0.5*xm;
sink.y=0.5*ym;
n=200;
allive=n;

p=0.1;
sv=0;
b=0;
ph=0.5;
Eo=0.5;

ETX=50*0.000000001;
ERx=50*0.000000001;
Efs=10*0.000000000001;
Emp=10*0.000000000001;
EDA=5*0.000000001;

h=100;
s=2;

tempi=50;
tempf=200;

m=0.1;
a=1;
u=a/2;
rmax=50;
allive=n;
E_adv=Eo*(1+a);
packets_TO_BS=0;
do=sqrt(Efs/Emp);
rcountCHs=0;




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
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    
     if(temp_rnd0<m*n+1)
         S(i).E=Eo;
         S(i).ENERGY=1;
         plot(S(i).xd,S(i).yd,'+');
         hold on;
     end  
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
countCHs=0;
cluster=1;
flag_first_dead=0;

dead=0;
first_dead=0;
allive=n;

packets_TO_BS=0;
packets_TO_CH=0;

for r=0:1:rmax
    r
    cv=tempi+(tempf-tempi)*rand(1,1);
    pnrm=(p/(1+a*m+b*u));
    padv=(p*(1+a)/(1+a*m+b*u));
    
    if(mod(r,round(1/pnrm))==0)
        for i=1:1:n
            if(S(i).ENERGY==0)
                S(i).G=0;
                S(i).cl=0;
            end
        end
    end
    
    if(mod(r,round(1/padv))==0)
        for i=1:1:n
            if(S(i).ENERGY==1)
                S(i).G=0;
                S(i).cl=0;
            end
        end
    end
    hold off;
    
    dead=0;
    dead_a=0;
    dead_n=0;
    
    packets_TO_CH=0;
    
    figure(1);
    
    for i=1:1:n
        if(S(i).E<=0)
            plot(S(i).xd,S(i).yd,'red');
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
                    plot(S(i).xd,S(i).yd,'o');
                end
                
                if(S(i).ENERGY==1)
                    plot(S(i).xd,S(i).yd,'+');
                end
                hold on;
            end
            totalenergy(i)=S(i).E;
    end
    figure(1)
    plot(S(n+1).xd,S(n+1).yd,'x');
    
    total_energy=sum(totalenergy);
    
    STATISTICS(r+1).DEAD=dead;
    DEAD(r+1)=dead;
    DEAD_N(r+1)=dead_n;
    DEAD_A(r+1)=dead_a;
    
    alive(r+1)=allive-dead;
    
    if(dead==1)
        if(flag_first_dead==0)
            first_dead=r;
            flag_first_dead=1;
        end
    end
    
    countChs=0;
    cluster=1;
    for i=1:1:n
        if(S(i).E>0)
            temp_rand=rand;
            if((S(i).G)<=0)
                if((S(i).ENERGY==0 && (temp_rand<=(pnrm/(1-pnrm*mod(r,round(1/pnrm)))))*(S(i).E/total_energy)))
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    S(i).type='C';
                    S(i).G=100;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    plot(S(i).xd,S(i).yd,'k*');
                    
                    
                    distance=sqrt((S(i).xd-(S(n+1).xd))^2+(S(i).yd-(S(n+1).yd))^2);
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    
                    
                    distance;
                    if(cv>=h)
                        test=cv-sv;
                        if(test>=s)
                            distance;
                            if(distance>do)
                                S(i).E=S(i).E-((ETX+EDA)*(4000)+Emp*4000*(distance*distance*distance*distance));
                            end
                            
                            if(distance<=do)
                                S(i).E=S(i).E-((ETX+EDA)*(4000)+Efs*4000*(distance*distance));
                            end
                        end
                    end
                end
                
                
                if((S(i).ENERGY==1 && (temp_rand<=(padv/(1-padv*mod(r,round(1/padv)))))*(S(i).E/total_energy)))
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    S(i).type='C';
                    S(i).G=100;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    plot(S(i).xd,S(i).yd,'k*');
                    
                    
                    distance=sqrt((S(i).xd-(S(n+1).xd))^2+(S(i).yd-(S(n+1).yd))^2);
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    
                    
                    distance;
                    if(cv>=h)
                        test=cv-sv;
                        if(test>=s)
                            if(distance>do)
                                S(i).E=S(i).E-((ETX+EDA)*(4000)+Emp*4000*(distance*distance*distance*distance));
                            end
                            
                            if(distance<=do)
                                S(i).E=S(i).E-((ETX+EDA)*(4000)+Efs*4000*(distance*distance));
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    STATISTICS(r+1).CLUSTERHEADS=cluster-1;
    CLUSTERHS(r+1)=cluster-1;
    
    
    for i=1:1:n
        if(S(i).type=='N'&&S(i).E>0)
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
                
                if(cv>=h)
                    test=cv-sv;
                    if(test>=s)
                        if(min_dis>do)
                            S(i).E=S(i).E-(ETX*(4000)+Emp*4000*(min_dis*min_dis*min_dis*min_dis));
                        end
                        if(min_dis<=do)
                            S(i).E=S(i).E-(ETX*(4000)+Efs*4000*(min_dis*min_dis));
                        end
                        
                        
                        
                        if(min_dis>0)
                          S(C(min_dis_cluster).id).E=S(C(min_dis_cluster).id).E-((ERx+EDA)*4000);
                          PACKETS_TO_CH(r+1)=n-dead-cluster+1;
                        end
                    end
                end
                S(i).min_dis=min_dis;
                S(i).min_dis_cluster=min_dis_cluster;
                
            end
        end
    end
    hold on;
   
                            
   countCHs;
   rcountCHs=rcountCHs+countCHs;
   sv=cv;
   inisv=sv;
   
   if(cv>=180&&r<=2500)
       sv=inisv;
   end
   
   
   [vx,vy]=voronoi(X,Y);
   %plot(X,Y,'r*',vx,vy,'b-');
   %hold on;
   voronoi(X,Y);
   axis([0 xm 0 ym]);
   
   P.PACKETS_TO_BS(r+1)=packets_TO_BS;
   pause(0.1)
end
r=0:50;
AlliveLeach=100-DEAD;


figure(2)
plot(r,DEAD,'--b');
legend('LEACH');
xlabel('Number of Rounds');
ylabel('Dead Nodes');

figure(3)
plot(r,alive,'--b');
legend('LEACH');
xlabel('Number of rounds');
ylabel('Alive Nodes');

figure(4)
plot(r,P.PACKETS_TO_BS,'--b');
legend('LEACH');
xlabel('Number of Rounds');
ylabel('Throughput');

figure(5)
plot(r,AlliveLeach,'--b');
legend('LEACH');
xlabel('Number of Rounds');
ylabel('Total number of Nodes');
