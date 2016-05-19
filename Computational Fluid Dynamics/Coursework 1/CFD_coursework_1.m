%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              CFD Coursework                             %
%                                                                         %
%                       Antoine Collier - CID 01145965                    %                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Question 4   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%Custom color matrix for prettiers plots, taken from the linspecer
%library which is not included with this work to hand back a single '.m' file
Colors_4 = [0.3467 0.5360 0.6907;0.9153 0.2816 0.2878;0.4416 0.7490 0.4322;1.0000 0.5984 0.2000];


%%%%%NOTE
%We could have used matlab's angle and abs functions to plot the dispersion and diffusion errors from
%the complex values of G and G~ but I chose to use the analytic expressions
%since using those functions gave me discontinuities in the dispersion error plotting.
%%%%%%%%%

phi=(0:0.01:pi); %phi phase vector

        %%%%%%%%%%%%%%%%%%%%FOR BETA VARYING%%%%%%

figure(1);

%Beta varying, sigma fixed
sigma=0.2;
beta_vector=[0.2 0.3 0.4 0.5];

top = subplot(211);
hold on 
bottom = subplot(212);
hold on %each subplot need its own hold on

%We now plot the diffusion and dispersion errors for each value of beta
for i = 1:size(beta_vector,2)
    
    Diffusion_error=exp(beta_vector(i)*phi.^2)./sqrt((1+2*beta_vector(i)*(1-cos(phi))).^2+(2*sigma*sin(phi)).^2);
    plot(top,phi,Diffusion_error,'color',Colors_4(i,:));
    legendInfoTop{i} = ['\beta = ' num2str(beta_vector(i))]; %We store the legend in a vector
    
    Dispersion_error=atan(sigma*sin(phi)./(1+2*beta_vector(i)*(1-cos(phi))))./(sigma.*phi);
    plot(bottom,phi,Dispersion_error,'color',Colors_4(i,:));
    legendInfoBottom{i} = ['\beta = ' num2str(beta_vector(i))]; %We store the legend in a vector
end 

%Tweaking of the figure
xlabel(top,'\phi','FontSize', 16);
ylabel(top,'\epsilon_D','FontSize', 16);
legend(top,legendInfoTop);
title(top,'Diffusion error','FontSize', 16);

xlabel(bottom,'\phi','FontSize', 16);
ylabel(bottom,'\epsilon_\phi','FontSize', 16);
legend(bottom,legendInfoBottom);
title(bottom,'Dispersion error','FontSize', 16);




          %%%%%%%%%%%%%%%%%%%%FOR SIGMA VARYING%%%%%%
          
figure(2);
%Sigma varying, beta fixed
beta=0.2;
sigma_vector=[0.5 1 2 3];

top = subplot(211);
hold on 
bottom = subplot(212);
hold on %each subplot need its own hold on


%We now plot the diffusion and dispersion errors for each value of sigma
for i = 1:size(sigma_vector,2) 
    
    Diffusion_error=exp(beta*phi.^2)./sqrt((1+2*beta*(1-cos(phi))).^2+(2*sigma_vector(i)*sin(phi)).^2); 
    plot(top,phi,Diffusion_error,'color',Colors_4(i,:));
    legendInfoTop{i} = ['\sigma = ' num2str(sigma_vector(i))]; %We store the legend in a vector
    
    Dispersion_error=atan(sigma_vector(i)*sin(phi)./(1+2*beta*(1-cos(phi))))./(sigma_vector(i).*phi);
    plot(bottom,phi,Dispersion_error,'color',Colors_4(i,:));
    legendInfoBottom{i} = ['\sigma = ' num2str(sigma_vector(i))]; %We store the legend in a vector
end 

%Tweaking of the figure
xlabel(top,'\phi','FontSize', 16);
ylabel(top,'\epsilon_D','FontSize', 16);
legend(top,legendInfoTop);
title(top,'Diffusion error','FontSize', 16);

xlabel(bottom,'\phi','FontSize', 16)
ylabel(bottom,'\epsilon_\phi','FontSize', 16);
legend(bottom,legendInfoBottom);
title(bottom,'Dispersion error','FontSize', 16);





%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Question 5   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Custom color matrix for prettiers plots, taken from the linspecer
%library which is not included with this work to hand back a single '.m' file
Colors_3 = [0.3467 0.5360 0.6907;0.9153 0.2816 0.2878;0.4416 0.7490 0.4322];

%Constants given by the Coursework sheet
a=1;
alpha=0.005;
delta_x=0.01;
delta_t=0.001;
T=0.5;
%Calulation of some values
sigma=a*delta_t/delta_x;
beta=alpha*delta_t/delta_x^2;
N=1/delta_x;


x_vector=transpose(0:0.01:1); %creation of the x vector 

%creation of the initial u vector
syms m;
U_0=10*(symsum((-1)^m*4*sin((2*m+1).*pi*x_vector)/((2*m+1)^2*pi^2), m, 0, 4)); %U_0 from the initial conditions given
U_0=double(U_0); %double function to convert U_0 from symbolic matrix to a numeric one


     %%%%%%%%%%%%%%%%%%%%%CASE 1 : CENTRED ADVECTION TERM %%%%%%%%%%

            %>>>>>>Creation of the matrix C

%Creation of the diagonal vectors taking place in the matrix
lower_diagonal=(beta+sigma/2)*ones(N,1);
diagonal=(1-2*beta)*ones(N+1,1);
upper_diagonal=(beta-sigma/2)*ones(N,1);

%Creation of the last row vector
Last_row=zeros(1,N+1);
Last_row(N)=1;


%Finally, creating C
C=diag(lower_diagonal,-1)+diag(diagonal)+diag(upper_diagonal,1); %Adding of the diagonals to the right place using the diag function
C(1,:)=zeros(1,N+1); %The first row is set to zero
C(N+1,:)=Last_row; %The last row is replaced by the propper one


            %>>>>>>Now we increment the solution using the matrix C until we reach the time T
            %(each multiplication account for a delta_t step)


U_n=U_0; %initialization

for i=0:T/delta_t
    U_n1=C*U_n;
    U_n=U_n1; 
end

Colors=linspecer(2);


%creation of u plot (CASE 1)
figure
hold on
plot(x_vector,U_0,':','color',Colors_3(1,:),'LineWidth',1.4);
plot(x_vector,U_n,'color',Colors_3(2,:),'LineWidth',1);



        %%%%%%%%%%%%%%%%%%%CASE 2 : UPWIND ADVECTION TERM %%%%%%%%%%

                %>>>>>>Creation of the matrix C

%Creation of the diagonal vectors taking place in the matrix
lower_diagonal=(beta+sigma)*ones(N,1);
diagonal=(1-2*beta-sigma)*ones(N+1,1);
upper_diagonal=(beta)*ones(N,1);

%Creation of the last row vector
Last_row=zeros(1,N+1);
Last_row(N)=1;

%Finally, creating C
C=diag(lower_diagonal,-1)+diag(diagonal)+diag(upper_diagonal,1); %Adding of the diagonals to the right place using the diag function
C(1,:)=zeros(1,N+1); %The first row is set to zero
C(N+1,:)=Last_row; %The last row is replaced by the propper one




                %>>>>>>Now we increment the solution using he matrix C until we reach the time T
                %(each multiplication account for a delta_t step)


U_n=U_0; %initialization


for i=0:T/delta_t
    U_n1=C*U_n;
    U_n=U_n1; 
end


%creation of u plot (CASE 2)

hold on

plot(x_vector,U_n,'color',Colors_3(3,:),'LineWidth',1);
xlabel('x','FontSize', 13);
lgd=legend('u(x,0)','centred numerical scheme','upwind numerical scheme','Location','Northwest');
