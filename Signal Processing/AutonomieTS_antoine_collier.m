%
% STI tc0 - Autonomie de Traitement du Signal - annee 2012-2013 - S6
% Estimation du milieu oceanique
%
%
% Elements relatifs a matlab fournis pour ce travail :
%   - le present fichier 'AutonomieTS.m', squelette du travail a realiser
%   - le fichier simulink 'Ocean.mdl', simulateur utilise dans la partie 3)
%   - l'image 'schemaSourceRecepteurSmall.jpg' utilisee par 'Ocean.mdl'
%   - le fichier de donnees 'data.mat' pour la partie 4)
%

% Nom : Collier
% Prenom : Antoine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Utiliser "Run and advance" plut�t que run pour voir les r�ponses
%%%%%%%%%%%%%% question par question (pr�sence d'un close all; au d�but de
%%%%%%%%%%%%%% chaque question)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------%
% 2.1) Construction du signal source deterministe
%
%% 

% Question 7 : Tracez et interpretez la fonction d'autocorrelation d'un signal MLF
% pour differentes valeurs de T puis de B.

close all;
clear;

%%Autocorr�lation pour diff�rentes valeurs de T

figure('name','Autocorr�lation d''un signal MLF en fonction de T');
hold on;
B=10;
k=1;
tau=[-1:1/1000:1];

color=['r','g','b','m','k'];

for T=[0.01,0.1,1,100,1000];
  R_MLF=T.*tripuls(tau/T,2).*sinc(B*tripuls(tau/T,2).*tau);
  subplot(5,1,k)
  plot(tau,R_MLF,color(k));
  k=k+1;
  legend(['T=' num2str(T)]);
end
hold on


%%On constate que plus T diminue, plus on se rapproche d'un dirac. C'est
%%normal, plus on diminue le support temporel, plus les fr�quences entre 0 et B
%%sont balay�es en un temps bref, et donc plus on se rapproche du dirac th�orique (infinit� des fr�quences balay�es en un temps infiniment court).
  

figure('name','Autocorr�lation d''un signal MLF en fonction de B');
hold on;
T=10;
k=1;
for B=[0.1,1,10,100,1000];
  R_MLF=T.*tripuls(tau/T,2).*sinc(B*tripuls(tau/T,2).*tau);
  subplot(5,1,k)
  plot(tau,R_MLF,color(k));
  k=k+1;
  legend(['B=' num2str(B)]);
end

%%Plus B augmente plus on se rapproche d'un dirac. Normal, th�oriquement un
%%dirac balaye toutes les fr�quences en un temps infiniment court. Plus B
%%est grand, plus on va dans ce sens
 

%%


% Question 8 : Construisez le vecteur des echantillons {xk}_{k=0:N-1} associe au 
% signal MLF x dont vous fixerez les parametres T et B en accord avec une frequence
% d'echantillonnage Ts de 10^{-4} sec.


Ts = 10^(-4); % p�riode d'echantillonnage (sec) % 
B=4000; %On prend B inf�rieur � B_max (voir apr�s)
T=1;  
kTs=0:Ts:T-Ts;
x_kTs=sin(pi*B*kTs.*kTs/T).*rectangularPulse((kTs-T/2)/T); 
figure('name','Allure du signal x �chantillonn� sur son support [0:T]');
stem(kTs, x_kTs,'.')

%Valeur max pour B : Th de shanon 2*nu_max <nu_s => soit avec le concept de
%fr�quence instant�e, qui est max en t=T : B_max=5000

%% 

% Question 9 : Calculez N echantillons de la transformee de Fourier X du signal
% continu x. Calculez sa fonction d'autocorrelation Rx et evaluez la densite spectrale
% d'energie Sx du signal de deux fa�ons differentes.



figure('name','Transform�e de Fourrier de x');
N=length(kTs);
X=Ts*fftshift(fft(x_kTs));
nu=-1/(2*Ts):1/(T-Ts):1/(2*Ts); 

plot(nu,abs(X))

R_x=xcorr(x_kTs,'unbiased'); 

tau=(-N+1:N-1)*Ts;

figure('name','Autocorr�lation du signal');
plot(tau,R_x)

%densit� spectrale d'�nergie m�thode 1 : on passe par le module du spectre
%de x au carr� :

figure('name','Densit� spectrale de x : 2 m�thodes');
subplot(211);
plot(nu,abs(X).^2)
title('Calcul avec S_x=|X|�');

%densit� spectrale d'�nergie m�thode 2 : on passe par la transform�e de
%Fourier de l'autocorr�lation

nu_2=-1/2/Ts:1/2/(T-Ts):1/2/Ts; %%Le vecteur doit avoir une taille de 2N-1
S_x=Ts*fftshift(fft(R_x));
subplot(212)
plot(nu_2,abs(S_x))
title('Calcul avec S_x=F[R_x]');



%% 

%-----------------------------------------------------------------------------------%
% 2.2) Construction du signal source aleatoire
%

% Question 12 : construisez le vecteur des echantillons {bk}_{k=0:N-1} associe a
% une realisation b d'un bruit blanc, gaussien, echantillonne a la
% periode Ts=10^{-4} sec et de duree T.


Ts = 10^(-4); % periode d'echantillonnage (sec) %
b_k=randn(1,N);
figure('name','Signal discret b* extrait d''une r�alisation b d''un signal de type bruit blanc');
stem(kTs,b_k);

%% 


% Question 13 : Calculez la transformee de Fourier, la fonction d'autocorrelation
% et la densite spectrale de cette realisation b.



figure('name','Transform�e de Fourrier de b �chantillonn�');
X=Ts*fftshift(fft(b_k)); %On multiplie par Ts car c'est Xn qui est renvoy� par fft(xkTs)
nu=-1/(2*Ts):1/(T-Ts):1/(2*Ts); 

plot(nu,abs(X))

R_b=xcorr(b_k); 

tau=(-N+1:N-1)*Ts;
figure('name','Autocorr�lation du signal b �chantillonn� ');
plot(tau,R_b)


figure('name','Densit� spectrale du signal b �chantillonn�');
plot(nu,abs(X).^2)
title('Calcul avec S_x=|X|�');




%% 

% Question 14 : Realisez un filtre et appliquez le au signal b afin de limiter
% son support spectral a l'intervalle [0;B] Hz. Representez le signal x obtenu en
% sortie du filtre, sa transformee de Fourier, sa fonction d'autocorrelation et sa
% densite spectrale d'energie.

%Le signal b_k enregistr� dans Matlab �tant le signal b �chantillonn�, on va donc lui appliquer un
%filtrage num�rique. Le filtre choisi est un passe bas de type butterworth
%(donc � RII) qui va �tre facile � construire avec un gabarit.


%On choisit les valeurs d'att�nuations max en bande passante et min en
%bande r�jection
R_p=0.5;
R_s=20;

nu_pass=B-B/20; % Fr�quence de fin de la bande passante
nu_stop=B; % Fr�quence de d�but de la bande de r�jection. Ainsi, le support spectral du signal sera bien [-B,B]

%On retranscrit ces caract�ristiques en temps continu
nu_contpass=1/pi/Ts*tan(pi*Ts*nu_pass);
nu_contstop=1/pi/Ts*tan(pi*Ts*nu_stop);

%Ordre et fr�quence de coupure du filtre de butterworth correspondant en
%temps continu

[n, wn]=buttord(2*pi*nu_contpass,2*pi*nu_contstop,R_p,R_s,'s');


%On calcule le num�rateur et le d�nominateur correspondant en temps continu

[B_c,A_c]=butter(n,wn,'s');

%On transpose avec la transformation bilin�aire (ou de Tustin)
[B_b,A_b]=bilinear(B_c,A_c,1/Ts);

%Calcul de la r�ponse fr�quentielle
H_b=freqz(B_b, A_b, nu, 1/Ts);

%trac� du module de la r�ponse fr�quentielle
figure('name','Module de la r�ponse fr�quentielle du filtre num�rique');
plot(nu,abs(H_b))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Comparaison avant/apr�s filtrage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Temporel

b_k_filtre=filtfilt(B_b,A_b,b_k);
figure('name','Comparaison temporelle signal avant puis apr�s filtrage')
subplot(241)
stem(kTs,b_k)
title('Signal b* g�n�r� par Matlab')

subplot(245)
stem(kTs,b_k_filtre)
title('Signal b* filtr� num�riquement')

%FFT

X=Ts*fftshift(fft(b_k)); 
nu=-1/(2*Ts):1/(T-Ts):1/(2*Ts);

subplot(242)
plot(nu,abs(X))
title('FFT du signal b �chantillonn�')

X_filtre=Ts*fftshift(fft(b_k_filtre));

subplot(246)
plot(nu,abs(X_filtre))
title('FFT du signal b �chantillonn� puis filtr�')

%Autocorr�lation
R_b=xcorr(b_k); 

subplot(243)
plot(tau,R_b)
title('Autocorr�lation du signal b �chantillonn�')

R_b_filtre=xcorr(b_k_filtre);
subplot(247)
plot(tau,R_b_filtre)
title('Autocorr�lation du signal b �chantillonn� puis filtr�')


%Densit� spectrale d'�nergie

subplot(244)
plot(nu,abs(X).^2)
title('DSE signal b echantillonn�');

subplot(248)
plot(nu,abs(X_filtre).^2)
title('DSE signal b echantillonn� puis filtr�');

%% 
%Question BONUS

figure('name','Comparaison de la DSE et de l''autocorr�lation des deux signaux sources �tudi�s')
subplot(221)
plot(tau,R_x)
title('Autocorr�lation du signal x');

subplot(222)
X=Ts*fftshift(fft(x_kTs));
plot(nu,abs(X).^2)
title('Calcul avec S_x=|X|�');


subplot(223)
R_b_filtre=xcorr(b_k_filtre);
plot(tau,R_b_filtre)
title('Autocorr�lation du signal b �chantillonn� puis filtr�')

subplot(224)
plot(nu,abs(X_filtre).^2)
title('DSE signal b echantillonn� puis filtr�')

%%Il n'y a pas photo, en terme de constance de la DSE sur l'intervale voulu
%%comme de d'absence de corr�lation en dehors de 0, le signal MLF est
%%meilleur. 
%% 


%-----------------------------------------------------------------------------------%
% 3) Mise en place de la methode
%

% simulation de la propagation %
% nom du signal source : xk (obligatoire)
% nom du signal re�u : yk
Ts = 10^(-4);   % periode d'echantillonnage (sec) %



%Etant donn� la r�ponse � la question bonus pr�c�dente, on choisi le signal
%MLF en source, qui est le plus proche du dirac.
xk=x_kTs;
sim('Ocean');   % lance le modele simulink 'Ocean.mdl' % 
y_kTs = yk(:)';    % mise sous forme de vecteur ligne %

%% 

% Question 15 : Observez le signal re�u apres propagation ainsi que son spectre.

N=length(x_kTs);
M=length(y_kTs);

kTs_2=0:Ts:(T-Ts)*(M-1)/(N-1);

figure('nam','Etude de la r�ponse impulsionnelle estim�e de l''oc�an')
subplot(311)
plot(kTs,x_kTs)
title('Signal �mis par la source')
subplot(312)
plot(kTs_2,y_kTs)
title('Signal y_kTs re�u par le r�cepteur')
subplot(313)
Y=Ts*fftshift(fft(y_kTs));
nu=-1/2/Ts:1/(T-Ts)*(N-1)/(M-1):1/2/Ts;
plot(nu,abs(Y))
title('Spectre du signal re�u par le r�cepteur')
%% 

% Question 16 : Implementez la methode d'estimation de la reponse impulsionnelle
% par intercorrelation. Representez la reponse impulsionnelle estimee 'hhat' ainsi
% que son spectre.

%On construit l'intercorr�lation entre y_kTs et x_kTs (cf question 3)


tau=(-M+1:M-1)*Ts;
hhat=xcorr(y_kTs,x_kTs);

figure('name','R�ponse impulsionnelle estim�e du milieu')
subplot(211)


plot(tau,hhat)
title('R�ponse impulsionnelle estim�e du milieu')

Hhat=Ts*fftshift(fft(hhat));
nu=-1/2/Ts:1/2/(T-Ts)*(N-1)/(M-1):1/2/Ts;

subplot(212)
plot(nu,abs(Hhat))
title('Spectre de la r�ponse impulsionnelle estim�e du milieu')


%% 

% Question 17 : Estimez la reponse impulsionnelle du milieu pour quelques valeurs
% des parametres T et B du signal. Observez, comparez et interpretez les resultats.







figure('name','R�ponse impulsionnelle estim�e du milieu pour diff�rentes valeur de B et T')
k=1;
B=1;
for T=[0.01 1 10];


kTs=0:Ts:T-Ts; 
x_kTs=sin(pi*B*kTs.*kTs/T).*rectangularPulse((kTs-T/2)/T);

xk=x_kTs;
sim('Ocean');   % lance le modele simulink 'Ocean.mdl' % 
y_kTs = yk(:)';

M=length(y_kTs);
tau=(-M+1:M-1)*Ts;
hhat=xcorr(y_kTs,x_kTs);

subplot(3,3,k)
%legend(['T=' num2str(T)],['B=' num2str(B)]);

plot(tau,hhat,'b')
legend(['T=' num2str(T)]);
if(k==1)
    title('B=1')
end
    
k=k+3;
end

k=2;
B=200;
for T=[0.01 1 10];


kTs=0:Ts:T-Ts; 
x_kTs=sin(pi*B*kTs.*kTs/T).*rectangularPulse((kTs-T/2)/T);

xk=x_kTs;
sim('Ocean');   % lance le modele simulink 'Ocean.mdl' % 
y_kTs = yk(:)';

M=length(y_kTs);
tau=(-M+1:M-1)*Ts;
hhat=xcorr(y_kTs,x_kTs);

subplot(3,3,k)
%legend(['T=' num2str(T)],['B=' num2str(B)]);

plot(tau,hhat,'r')
legend(['T=' num2str(T)]);
if(k==2)
    title('B=200')
end
    
k=k+3;
end

k=3;
B=4000;
for T=[0.01 1 10];


kTs=0:Ts:T-Ts; 
x_kTs=sin(pi*B*kTs.*kTs/T).*rectangularPulse((kTs-T/2)/T);

xk=x_kTs;
sim('Ocean');   % lance le modele simulink 'Ocean.mdl' % 
y_kTs = yk(:)';

M=length(y_kTs);
tau=(-M+1:M-1)*Ts;
hhat=xcorr(y_kTs,x_kTs);

subplot(3,3,k)
%legend(['T=' num2str(T)],['B=' num2str(B)]);

plot(tau,hhat,'k')
legend(['T=' num2str(T)]);
if(k==3)
    title('B=4000')
end
    
k=k+3;
end

%% 

% Question 19 D�terminez la valeur de Z

%Param�tres de la zone �tudi�e
c=1510;
R=1000;
Z_s=30;
Z_r=50;

%Temps des contributions des trajets direct, par r�flexion sur le fond
%puis sur la surface(rep�r�s gr�ce � leur signe) obtenus sur la r�ponse
%impulsionnelle
tau_1=0.01; %Trajet direct
tau_2=0.0119; %R�flexion sur la surface
tau_3=0.0256; %r�flexion sur le fond



syms Z;

%Distance points de r�flexion sur la surface et le fond
x=R*Z_s/(Z_r+Z_s); %Sur la surface
y=R*(Z-Z_s)/(2*Z-Z_s-Z_r); %Sur le fond

%Distance parcourue lors des 3 trajets
d1=sqrt(R^2+(Z_r-Z_s)^2); %Direct
d2=sqrt(Z_s^2+x^2)+sqrt(Z_r^2+(R-x)^2); %R�flexion sur la surface
d3=sqrt((Z-Z_s)^2+y^2)+sqrt((Z-Z_r)^2+(R-y)^2); %R�flexion sur le fond


%Temps des trajets
t1=d1/c;
t2=d2/c;
t3=d3/c;

sol=solve(t3==(tau_3-tau_1)/(tau_2-tau_1)*(t2-t1)+t1,'Z');
Z=eval(sol(1))

%La profondeur estim�e est donc Z=152 m
%% 

%-----------------------------------------------------------------------------------%
% 4) Determination des parametres du milieu
%

% Question 20 : on a mesure le signal y_experience2 qui est donne dans le
% fichier data.mat en reponse a un signal source MLF de duree T = 0.2 sec
% et de bande B=2000 Hz. Determinez la vitesse de propagation c2 de l'onde
% acoustique dans le milieu.

load data;

figure('name','R�ponse impulsionnelle estim�e du milieu dans l''exp�rience 2'); 

T=0.2;               
t=0:Ts:T-Ts;           
N=length(t); 
M=length(y_experience2);
B=2000;              

x=sin(pi*B*t.*t/T).*rectangularPulse((t-T/2)/T);





R_xy=xcorr(y_experience2,x);      
tau=(-M+1:M-1)*Ts;  


plot(tau,R_xy, 'b');


y=R*(Z-Z_s)/(2*Z-Z_s-Z_r);
d3=sqrt((Z-Z_s)^2+y^2)+sqrt((Z-Z_r)^2+(R-y)^2);
tau_1=0.01;
tau_2=0.012;
tau_3=0.0263;

c2=(d3-d1)/(tau_3-tau_1)



