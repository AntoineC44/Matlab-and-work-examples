%AFD Coursework
clear all
close all

%Moment on the studied section
%8 cases, matrix with the sign coefficient
Weight_load = 32000;
Tail_load_1 = 15000;
Tail_load_2 = 5000;
Vertical_stabilizer_load = 5000;
cases=[1 1 -1 1;1 1 -1 -1;1 -1 1 1;1 -1 1 -1;-1 1 -1 1;-1 1 -1 -1;-1 -1 1 1;-1 -1 1 -1];

M=zeros(3,8);
for k=1:8
    M(1,k) = cases(k,1)*(cases(k,2)*Tail_load_1 + cases(k,3)*Tail_load_2) - 0.9*cases(k,4)*Vertical_stabilizer_load;
    M(2,k) = -Weight_load*2-cases(k,1)*4*(Tail_load_1 + Tail_load_2);
    M(3,k) = cases(k,4)*3.5*Vertical_stabilizer_load;
end



%Entry parameters
N_vector=(6:1:80); %Nb stringers
t_vector=(0.0001:0.0001:0.004); %skin thickness
alpha_offset=0; %offset

%%Moments
Mx=M(1,2);
My=M(2,2);
Mz=M(3,2);

%%Section parameters
D = 1.2; %Section diameter
E_skin = 72.4*10^9; %Modulus of elasticity skin (Aluminum 2024-T6)

%%Stringers parameters
dh_ratio = 0.4; %Stringer ratio d/h
E_stringer = 71.7*10^9; %Modulus of elasticity stringer (Aluminum 7075-T6)

%Loop speed
Total_area=zeros(size(N_vector,2),size(t_vector,2));
Failure=zeros(size(N_vector,2),size(t_vector,2));
Failure_criterion=zeros(size(N_vector,2),size(t_vector,2));

for i=1:size(N_vector,2)
    for j=1:size(t_vector,2)

            N=N_vector(i);
            t=t_vector(j);
            
            %Relevent dimensions
            b = pi*D/N; %Beam distance
            p=b-30*t; %effective beam distance
            As = b*t; %Stringer area
            
            %Location of stringers
            Y = zeros(1,N);
            Z = zeros(1,N);
            for k=1:size(Y,2)
                Y(k)=cos(alpha_offset+(N/4+1-k)*pi/(N/2))*D/2; % Ycoordinate stringers
                Z(k)=sin(alpha_offset+(N/4+1-k)*pi/(N/2))*D/2; % Zcoordinate stringers
            end
            
            %For loop speed
            A = zeros(1,N);
            sigmaxx = zeros(1,N);
                     
            %Creation of equivalent area of the skin/stringer distribution
            A(1)=As ...
                    + t*p*(2+ (My*Z(2)-Mz*Y(2))/(My*Z(1)-Mz*Y(1)))/6 ...
                    + t*p*(2+ (My*Z(1)-Mz*Y(1))/(My*Z(N)-Mz*Y(N)))/6;

            for k=2:size(A,2)-1
                A(k)=As ...
                    + t*p*(2+ (My*Z(k+1)-Mz*Y(k+1))/(My*Z(k)-Mz*Y(k)))/6 ...
                    + t*p*(2+ (My*Z(k)-Mz*Y(k))/(My*Z(k-1)-Mz*Y(k-1)))/6;
            end

             A(N)=As ...
                    + t*p*(2+ (My*Z(1)-Mz*Y(1))/(My*Z(N)-Mz*Y(N)))/6 ...
                    + t*p*(2+ (My*Z(N)-Mz*Y(N))/(My*Z(N-1)-Mz*Z(N-1)))/6;

             %Creation of the second moment of inertia
             Iyy=0;
             Izz=0;

             for k=1:N
                Iyy=Iyy+A(k)*Y(k)^2;
                Izz=Izz+A(k)*Z(k)^2;
             end    

             for k=1:N
                sigmaxx(k)=My*Z(k)/Iyy - Mz*Y(k)/Izz;
             end


             %stringers-skin panel compression
             K=5.2; %Skin stringer panel correction factor
             sigma_cr=K*E_stringer*(t/b)^2; %Initial buckling stress
             Rc=max(abs(sigmaxx))/sigma_cr;


             %Uniform skin shear
             qtwist=Mx/2/pi/(D/2)^2;
             tau=qtwist/t;
             T=(t+As/p); %Effective thickness
             tau_cr=4.86*E_skin*(T/p)^2;
             Rs=abs(tau)/tau_cr;

             Total_area(i,j)=(N*As+pi*(((D+2*t)/2)^2-(D/2)^2))*10^6;
             Failure_criterion(i,j)=Rs^2+Rc;
             %Failure test
             if Rs^2+Rc>0.99;
                 Failure(i,j)=1;
             end
                 
    end
end   

%%%%%%%%%Failure plot
figure1 = figure('PaperSize',[20.98404194812 29.67743169791],...
    'Colormap',[0.0666666701436043 0.858823537826538 0.0666666701436043;0.0814814865589142 0.845191419124603 0.0656084716320038;0.0962963029742241 0.831559300422668 0.0645502656698227;0.111111111938953 0.817927181720734 0.0634920671582222;0.125925928354263 0.804295063018799 0.0624338649213314;0.140740737318993 0.790662944316864 0.0613756664097309;0.155555561184883 0.777030825614929 0.0603174641728401;0.170370370149612 0.763398706912994 0.0592592619359493;0.185185194015503 0.74976658821106 0.0582010596990585;0.200000002980232 0.736134469509125 0.057142861187458;0.214814811944962 0.72250235080719 0.0560846589505672;0.229629635810852 0.708870232105255 0.0550264567136765;0.244444444775581 0.69523811340332 0.053968258202076;0.259259253740311 0.681605994701385 0.0529100559651852;0.274074077606201 0.667973875999451 0.0518518537282944;0.288888901472092 0.654341757297516 0.0507936552166939;0.30370369553566 0.640709638595581 0.0497354529798031;0.31851851940155 0.627077519893646 0.0486772507429123;0.333333343267441 0.613445401191711 0.0476190485060215;0.348148137331009 0.599813282489777 0.046560849994421;0.362962961196899 0.586181163787842 0.0455026477575302;0.37777778506279 0.572549045085907 0.0444444455206394;0.39259260892868 0.558916926383972 0.0433862470090389;0.407407402992249 0.545284807682037 0.0423280447721481;0.422222226858139 0.531652688980103 0.0412698425352573;0.43703705072403 0.518020570278168 0.0402116440236568;0.451851844787598 0.504388451576233 0.0391534417867661;0.466666668653488 0.490756303071976 0.0380952395498753;0.481481492519379 0.477124184370041 0.0370370373129845;0.496296286582947 0.463492065668106 0.035978838801384;0.51111114025116 0.449859946966171 0.0349206365644932;0.525925934314728 0.436227828264236 0.0338624343276024;0.540740728378296 0.422595709562302 0.0328042358160019;0.555555582046509 0.408963590860367 0.0317460335791111;0.570370376110077 0.395331472158432 0.0306878332048655;0.585185170173645 0.381699353456497 0.0296296309679747;0.600000023841858 0.368067234754562 0.028571430593729;0.614814817905426 0.354435116052628 0.0275132283568382;0.629629611968994 0.340802997350693 0.0264550279825926;0.644444465637207 0.327170878648758 0.0253968276083469;0.659259259700775 0.313538759946823 0.0243386253714561;0.674074053764343 0.299906641244888 0.0232804249972105;0.688888907432556 0.286274522542953 0.0222222227603197;0.703703701496124 0.272642403841019 0.0211640223860741;0.718518495559692 0.259010285139084 0.0201058220118284;0.733333349227905 0.245378151535988 0.0190476197749376;0.748148143291473 0.231746032834053 0.017989419400692;0.762962937355042 0.218113914132118 0.0169312171638012;0.777777791023254 0.204481795430183 0.0158730167895555;0.792592585086823 0.190849676728249 0.0148148154839873;0.807407379150391 0.177217558026314 0.0137566141784191;0.822222232818604 0.163585439324379 0.0126984138041735;0.837037026882172 0.149953320622444 0.0116402124986053;0.851851880550385 0.136321201920509 0.010582011193037;0.866666674613953 0.122689075767994 0.00952380988746881;0.881481468677521 0.109056957066059 0.0084656085819006;0.896296322345734 0.0954248383641243 0.00740740774199367;0.911111116409302 0.0817927196621895 0.00634920690208673;0.92592591047287 0.0681606009602547 0.00529100559651852;0.940740764141083 0.0545284785330296 0.0042328042909503;0.955555558204651 0.0408963598310947 0.00317460345104337;0.970370352268219 0.0272642392665148 0.00211640214547515;0.985185205936432 0.0136321196332574 0.00105820107273757;1 0 0]);

% Create axes
axes1 = axes('Parent',figure1,'FontSize',14,'CLim',[0 1]);
ylim(axes1,[0 75]);
grid(axes1,'on');
hold(axes1,'all');

% Create mesh
mesh(Failure,'Parent',axes1,'MarkerSize',3,'Marker','o','LineStyle','none');
xlabel('10*t (mm)','FontSize',14);
ylabel('N','FontSize',14);





%%%%%%%%%Total area
figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2,'FontSize',14);
ylim(axes1,[0 75]);
grid(axes1,'on');
hold(axes1,'all');

% Create mesh
mesh(Total_area,'Parent',axes1,'MarkerSize',3,'Marker','o','LineStyle','none');
xlabel('10*t (mm)','FontSize',14);
ylabel('N','FontSize',14);



%%%%%%%%%Area stringers plot
N=36; %Chosen number of stringers
b=pi*D/N; %beam distance
figure3 = figure;

% Create axes
axes1 = axes('Parent',figure3,'FontSize',14);
ylim(axes1,[0 3]);
box(axes1,'on');
hold(axes1,'all');

plot(t_vector,Failure_criterion(N,:))
hold on
plot(t_vector,0.99*ones(size(t_vector,2)),'r')
legend({'Failure criterion','Limit'})

% Create xlabel
xlabel('t (mm)','FontSize',14);
ylabel('Failure criterion','FontSize',14);
legend(axes1,'show');

t=0.0015;
ts=t*1.25;
As=t*b;
h=As/(ts+2*0.4*ts);
d=h*0.4;

%%
    


%%Sketch
figure1 = figure;
hold on

axes1 = axes('Visible','Off','Parent',figure1);
xlim(axes1,[-0.6 0.6]);
ylim(axes1,[-0.6 0.6]);
box(axes1,'on');
hold(axes1,'all');

%Generating the skin
nbpoints=1000;
for k=1:nbpoints
    Y1(k)=cos(alpha_offset+k*2*pi/nbpoints)*D/2; % Ycoordinate
    Z1(k)=sin(alpha_offset+k*2*pi/nbpoints)*D/2; % Zcoordinate
end

%Generating the stringers
stringer0 = [D/2,ts/2 ; D/2-h+ts,ts/2 ; D/2-h+ts,ts/2+d ; D/2-h,ts/2+d ; D/2-h,-ts/2 ; D/2-ts,-ts/2 ; D/2-ts,-ts/2-d ; D/2,-ts/2-d]';
stringer = zeros(2,8);

for n = 1:N
    theta = 2*pi*n/N;
     for i=1:8
        stringer(:,i) = [cos(theta) , -sin(theta) ; sin(theta) , cos(theta)]*stringer0(:,i);
    end;
    for i=1:7
        plot([stringer(1,i),stringer(1,i+1)],[stringer(2,i),stringer(2,i+1)],'k')
    end;
        plot([stringer(1,8),stringer(1,1)],[stringer(2,8),stringer(2,1)],'k')
end;

plot(Y1,Z1) %inner skin layer

plot(Y1*(D/2+t)/(D/2),Z1*(D/2+t)/(D/2)) %outer skin layer
daspect([1 1 1])


hold on
axis equal
axis off






