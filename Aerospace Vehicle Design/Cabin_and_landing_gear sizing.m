%%Aerospace Vehicle Design
clear all
close all

%*********** Parameters for the problem

%Weight parameters
AverageLuggageSizePerPassenger=17;
Male_FemaleRatio=0.7;
AverageFemaleWeight=75; %with hand baggage and infants below 2 years of age
AverageMaleWeight=94; %with hand baggage and infants below 2 years of age

%specific fuel consumption (mg/Ns)
Sfc_cruise=14.1;
Sfc_loiter=11.3;

%Earth gravity field accelaration
g=9.81;

%Max lift drag ratio parameters
Kld=15.5; %given, for a civil jet
%[Aspect_ratio,W0]=meshgrid(5:0.1:10,100000:100:1000000);
Aspect_ratio=7.5; %arbitrary between B747 and B777
Wetted_aspect_ratio=6; %from B747 in graph in lecture

%Lift/Drag optimal ratio factor
lift_drag_loiter_factor=1;
lift_drag_cruise_factor=0.866;

%Composite structure correction factor
composite_factor=0.95;

%***********


%Payload weight
NbPassengers=590;
W_Payload = NbPassengers * (AverageLuggageSizePerPassenger + Male_FemaleRatio * AverageMaleWeight + (1 - Male_FemaleRatio) * AverageFemaleWeight);

%%Crew Weight
NbCrew=12;
W_Crew = NbCrew * (AverageLuggageSizePerPassenger + Male_FemaleRatio * AverageMaleWeight + (1 - Male_FemaleRatio) * AverageFemaleWeight);

%Mission Profile
W_1_0=0.970;
W_2_1=1.0065-0.0325*0.82; %Climbing at 10 km
%W_3_2 to determine
W_4_3=0.995; %"Landing"  at 0 km height
%W_5_4=0.985; %Climbing at 4km
%W_6_5 to determine
%%Maybe an antoher step from 4km to loiter here
%W_7_6 to determine
W_8_7=0.995; %Depending on the height chosen
W_9_8=0.995; %Same

Lift_drag_ratio_max=Kld*sqrt(Aspect_ratio./Wetted_aspect_ratio); %Lift drag ratio (needed for Breguet Range equations)

%%%%W_3_2 :

T_10000=15.05-0.00649*10000+273.15 %Assuming air is a perfect gas
M_3_2=0.82;
V_3_2=M_3_2*sqrt(1.4*287*T_10000);
W_3_2=exp(-11000000*g*Sfc_cruise/1000000./(V_3_2*lift_drag_cruise_factor*Lift_drag_ratio_max)); %Breguet range equation

%%%%W_6_5 :
altitude=5000;
Mach_diversion=0.5;
W_5_4=1.0065-0.0325*Mach_diversion;
T_6741=15.05-0.00649*altitude+273.15;
V_6_5=Mach_diversion*sqrt(1.4*287*T_6741);

W_6_5=exp(-370000*g*Sfc_cruise/1000000./(V_6_5*lift_drag_cruise_factor*Lift_drag_ratio_max)); %Breguet range equation
W_6_5
6741*cos(8*3.14/360)/sin(8*3.14/360)

%%%%W_7_6 :

W_7_6=exp(-45*60*Sfc_loiter/100000./(lift_drag_loiter_factor*Lift_drag_ratio_max)); %Breguet endurance equation

Fuel_fraction=1.02*(1-W_1_0.*W_2_1.*W_3_2.*W_4_3.*W_5_4.*W_6_5.*W_7_6.*W_8_7.*W_9_8); %Total fuel fraction with correction

A=0.97;
c=-0.06;
W0=linspace(0,1000000,100);
Empty_ratio=composite_factor*A*W0.^c; %factor to account for a composite aircraft

plot(W0,W0)

hold on

W02=(W_Crew+W_Payload)./(1-Fuel_fraction-Empty_ratio);
plot(W0,W02)



% for i=1:size(W02,1)
%     for j=1:size(W02,2)
%         if (abs(W02(i,j)-W0(i,j))>1000)
%             W02(i,j)=0;
%         end
%     end
% end
% mesh(W0,Aspect_ratio,W02)
    



% oswald_efficiency_factor=1/(1.05+0.007*3.14*Aspect_ratio); %http://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
% K=1/(3.14*Aspect_ratio*oswald_efficiency_factor); %Oswald Span Efficiency Method
% S=510; %B747 %better to find one ourselves ?
% C_fe=0.0026; %given
% C_D0=C_fe*Wetted_aspect_ratio;
% T=273.15-17.5;
% pho=1.292*273.15./T;
% W=W0*W_1_0*W_2_1*W_3_2*W_4_3*W_5_4;
% 
% Vmd=(K/C_D0)^(1/4)*sqrt(2*W./(pho*S))/sqrt(1.4*287*T)

%%
%undercarriage constraint diagram
close all
clear

cg_position_to_the_fuselage_bottom=2.3;
cg_height_to_the_ground=[0:0.1:4]+cg_position_to_the_fuselage_bottom;
dihedral_angle=pi*5/180;
engine_distance_to_fuselage=12;
engine_height=3.1;
fuselage_width=5.922;
track_distance=9;
wing_fixation_point_in_reference_to_cg=-1;
distance_undercarriage_to_engine=engine_distance_to_fuselage+fuselage_width/2-track_distance/2;

distance_nose_cg=36;
distance_nose_mg=40.6;
distance_nose_ng=4.5;

most_aft_cg=38;
most_forward_cg=34;

minimum_vertical_distance_engine_to_ground=distance_undercarriage_to_engine*tan(pi*5/180)+6*2.54/cos(5*pi/180)/100;
minimum_engine_clearance_angle=atan(minimum_vertical_distance_engine_to_ground/distance_undercarriage_to_engine)*ones(size(cg_height_to_the_ground));

engine_height_to_the_ground=cg_height_to_the_ground+wing_fixation_point_in_reference_to_cg+engine_distance_to_fuselage*sin(dihedral_angle)-engine_height;
engine_clearance_angle=atan(engine_height_to_the_ground/distance_undercarriage_to_engine);

%
angle_nose_main_gear=atan(track_distance/2/(distance_nose_mg-distance_nose_ng));
cg_distance_to_nose_main_gear_axis=sin(angle_nose_main_gear)*(distance_nose_cg-distance_nose_ng);

forward_turn_angle=atan(cg_height_to_the_ground./cg_distance_to_nose_main_gear_axis);
max_overturn_angle=63*ones(size(cg_height_to_the_ground));
%
static_tail_angle=20*ones(size(cg_height_to_the_ground));


vertical_cg_angle_main_aft=atan((distance_nose_mg-most_aft_cg)./cg_height_to_the_ground);
vertical_cg_angle_main_forward=atan((distance_nose_mg-most_forward_cg)./cg_height_to_the_ground);

%%%
Colors=linspecer(7);

figure
hold on
p1=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,max_overturn_angle,'--','Color',Colors(1,:),'LineWidth',1);
p2=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,forward_turn_angle*180/pi,'Color',Colors(1,:),'LineWidth',2);
p3=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,minimum_engine_clearance_angle*180/pi,'--','color',Colors(3,:),'LineWidth',1);
p4=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,engine_clearance_angle*180/pi,'color',Colors(3,:),'LineWidth',2);
p5=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,vertical_cg_angle_main_aft*180/pi,'color',Colors(5,:),'LineWidth',2);
p6=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,vertical_cg_angle_main_forward*180/pi,'.','color',Colors(5,:),'LineWidth',2);
p7=plot(cg_height_to_the_ground-cg_position_to_the_fuselage_bottom,static_tail_angle,'--','color',Colors(5,:),'LineWidth',1);

xlabel('Height fuselage above ground (m)');
ylabel('Angle (deg)');
legend([p1 p2 p3 p4 p5 p6 p7], {'Max overturn angle', 'Overturn angle','Minimum engine clearance angle','Engine clearance angle','Most aft CG vertical angle','Most forward CG vertical angle','Static tail angle'});

%%
syms beta;
syms gamma;

xcg=36;
xng=4.5;
xmg=40.6;
xrg=41.6;
eq1=xcg-0.1*xng-beta*xmg-gamma*xrg==0;
eq2=0.1+beta+gamma==1;

[A,B] = equationsToMatrix([eq1 eq2],[beta gamma])






