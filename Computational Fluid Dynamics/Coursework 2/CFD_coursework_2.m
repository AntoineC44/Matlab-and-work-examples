%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              CFD Coursework 2                           %
%                                                                         %
%                       Antoine Collier - CID 01145965                    %                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Question 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear
close all

%Custom color matrix for prettiers plots, taken from the linspecer
%library which is not included with this work to hand back a single '.m' file
Colors_4 = [0.3467 0.5360 0.6907;0.9153 0.2816 0.2878;0.4416 0.7490 0.4322;1.0000 0.5984 0.2000];


%%%%%NOTE
%We could have used matlab's angle and abs functions to plot the dispersion and diffusion errors from
%the complex values of G and G~ but I chose to use the analytic expressions
%since using those functions gave me discontinuities in the dispersion error plotting.
%%%%%%%%%
figure
hold on
phi=(0:0.01:pi); %phi phase vector
sigma_vector =[0.1 0.5 1 1.5];

%We now plot the amplification factor for each value of sigma
for i = 1:size(sigma_vector,2)
    Amplification_factor = (cos(phi/2).^2+(sigma_vector(i)*sin(phi/2)).^2).^2;
    
    plot(phi,Amplification_factor,'color',Colors_4(i,:),'LineWidth',2);
    legendgraph{i} = ['\sigma = ' num2str(sigma_vector(i))]; %We store the legend in a vector  
end 

%Tweaking of the figure
xlabel('\Phi','FontSize', 16);
ylabel('\midG\mid','FontSize', 16);
legend(legendgraph,'Location','northwest');
%set(gca,'xlim',[0 pi()],'xtick',[0:pi()/4:pi()],'xticklabel',{'0' 'p/4' 'p/2' '3p/4' 'p'})
title('\midG\mid vs \Phi','FontSize', 16);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Question 3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%Custom colors from 'linspecer.m' file
Colors_4 = [0.3467 0.5360 0.6907;0.9153 0.2816 0.2878;0.4416 0.7490 0.4322;1.0000 0.5984 0.2000];


%Loading the data from the .dat file created by 'Shock_tube_analytic.m'

analytic_data=importdata('shock_analytic.dat',' ');
xx = analytic_data.data(:,1);
rho_analytic = analytic_data.data(:,2);
p_analytic = analytic_data.data(:,3);
velocity_analytic = analytic_data.data(:,4);
mach_analytic = analytic_data.data(:,5);
entropy_analytic = analytic_data.data(:,6);

%Creating the plots
figure(1); hold on; plot (xx,velocity_analytic,'color',Colors_4(1,:),'LineStyle','-','LineWidth', 1.5); title('Velocity vs x'); xlabel('x'); ylabel ('u'); grid on; 
figure(2); hold on; plot (xx,p_analytic,'color',Colors_4(1,:),'LineStyle','-','LineWidth', 1.5); title('Pressure vs x'); xlabel('x'); ylabel ('p'); grid on; 
figure(3); hold on; plot (xx,rho_analytic,'color',Colors_4(1,:),'LineStyle','-','LineWidth', 1.5); title('Density vs x'); xlabel('x'); ylabel ('\rho'); grid on;
figure(4); hold on; plot (xx,mach_analytic,'color',Colors_4(1,:),'LineStyle','-','LineWidth', 1.5); title('Mach vs x'); xlabel('x'); ylabel ('M'); grid on;
figure(5); hold on; plot (xx,entropy_analytic,'color',Colors_4(1,:),'LineStyle','-','LineWidth', 1.5); title('Entropy vs x'); xlabel('x'); ylabel ('s'); grid on;


%Setting flow parameters
p2=10;
rho2=8;
p1=p2/10;
rho1=rho2/8;
gamma=1.4;

%Time parameter
T=0.5;

%The 3 values for N told by the coursework sheet
for N=[101 201 301];
    
    delta_x=4/(N+1);
    
    %Using the CFL condition with the quicket component of the flow
    u_max = max(velocity_analytic); 
    c_max = max(sqrt(1.4*p_analytic./rho_analytic));
    delta_t = delta_x/(u_max+c_max);
    
    lambda_m=delta_t/delta_x;

    x=transpose(-2:delta_x:2); 

    %Creation of the initial u vector
    U0=zeros(3,N+2);

        %Setting initial value for the first component of U at time 0
        U0(1,(1:1:(N+1)/2+1))=rho2*ones(1,(N+1)/2+1);
        U0(1,((N+1)/2)+2:1:N+2)=rho1*ones(1,(N+1)/2);

        %No need to set the second component since it is zero on all the domain

        %Setting initial value for the third component of U at time 0
        U0(3,(1:1:(N+1)/2+1))=p2*ones(1,(N+1)/2+1)/(gamma-1);
        U0(3,((N+1)/2)+2:1:N+2)=p1*ones(1,(N+1)/2)/(gamma-1);


    %Initialization
    Un=U0;
    Un1_2=zeros(3,N+1);
    Un1=zeros(3,N);
    
    %Iteration
    for n=1:T/delta_t
        
    %We first apply the (n) to (n+1/2) part of the two-step Lax-Friedrich scheme
        for i=1:N+1

            %We calculate the composants of F to have a more readable code
            F_Un=[Un(2,:);...
                (3-gamma)/2*Un(2,:).^2./Un(1,:)+(gamma-1)*Un(3,:);...
                gamma*Un(2,:).*Un(3,:)./Un(1,:)-(gamma-1)/2*Un(2,:).^3./(Un(1,:).^2)];

            Un1_2(:,i)=0.5*(Un(:,i+1)+Un(:,i)) - lambda_m/2*(F_Un(:,i+1)-F_Un(:,i));

        end

    %We then apply the (n+1/2) to (n+1) part of the two-step Lax-Friedrich scheme
        for i=1:N

            F_Un1_2=[Un1_2(2,:);...
                (3-gamma)/2*Un1_2(2,:).^2./Un1_2(1,:)+(gamma-1)*Un1_2(3,:);...
                gamma*Un1_2(2,:).*Un1_2(3,:)./Un1_2(1,:)-(gamma-1)/2*Un1_2(2,:).^3./(Un1_2(1,:).^2)];

           Un1(:,i)=0.5*(Un1_2(:,i+1)+Un1_2(:,i)) - lambda_m/2*(F_Un1_2(:,i+1)-F_Un1_2(:,i));

        end

    %We update the time (n) vecotr U with the values of U at time (n+1) while adding the boundary conditions U0=U1 and UN+1=UN
        Un=[Un1(:,1) Un1 Un1(:,N)]; 
    end
    
    %We now create the relevant values from the final Un vector
    u_num = Un(2,:)./Un(1,:);
    p_num = 0.4*(Un(3,:)-Un(2,:).^2./Un(1,:)/2);
    rho_num = Un(1,:);
    mach_num = u_num./sqrt(1.4*p_num./rho_num);
    entropy_num = log(p_num./(rho_num.^1.4));
    
    %We use a different color for each N
    if N==101 
        color = Colors_4(2,:);
    elseif N==201
        color = Colors_4(3,:);
    else
        color = Colors_4(4,:);
    end;
    
    figure(1), plot(x,u_num,'--','Color',color,'Linewidth',1.5)
    figure(2), plot(x,p_num,'--','Color',color,'Linewidth',1.5)
    figure(3), plot(x,rho_num,'--','Color',color,'Linewidth',1.5)
    figure(4), plot(x,mach_num,'--','Color',color,'Linewidth',1.5)
    figure(5), plot(x,entropy_num,'--','Color',color,'Linewidth',1.5)

end

%Legend
figure(1), legend('Analytic solution','N = 101','N = 201','N = 301')
figure(2), legend('Analytic solution','N = 101','N = 201','N = 301')
figure(3), legend('Analytic solution','N = 101','N = 201','N = 301')
figure(4), legend('Analytic solution','N = 101','N = 201','N = 301')
figure(5), legend('Analytic solution','N = 101','N = 201','N = 301','Location','NorthWest')



