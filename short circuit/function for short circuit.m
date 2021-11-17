function fval = libfunc(~,T)
xsei=T(1);
z=T(2);
xne=T(3);
alp=T(4);
soc=T(5);
TT=T(6);

T0=425;
%% calculating for the sei layer 
%constants

asei=1.67e15;
Esei=2.24e-19;
kb=1.381e-23;
msei=8.1e-3;
Hsei=257000;

%equations

    fval(1,1)=-xsei*asei*exp(-Esei/(kb*TT));
    Qsei=-msei*Hsei*(-xsei*asei*exp(-Esei/(kb*TT)));


%% calculations for the anode
%constants

ane=2.5e13;
Ene=2.24e-19;
z0=0.033;
mne=8.1e-3;
Hne=1714000;

%equations

    fval(2,1)=xne*ane*exp(-Ene/(kb*TT))*exp(-z/z0);
    fval(3,1)=-xne*ane*exp(-Ene/(kb*TT))*exp(-z/z0);
    Qne=-mne*Hne*(-xne*ane*exp(-Ene/(kb*TT))*exp(-z/z0));



%% calculations for the cathode
% constants
mpe=18.3E-3;
ape=6.67E11;
Hpe=314000;
Epe=2.03E-19;
% equations

    fval(4,1)=alp*(1-alp)*ape*exp(-Epe/(kb*TT));
    Qpe=mpe*Hpe*(alp*(1-alp)*ape*exp(-Epe/(kb*TT)));

%% calculations for the short circuit
% constants
asoc=3.37E12;
Esoc=1.58E-19;
V=4.2;
C=2.4;
% equations
if(TT<330.15)
    isc=0;
else
    isc=1;
end    
fval(5,1)= -isc*soc*asoc*exp(-Esoc/(kb*TT));
Hsoc=V*C*0.28*3600;
Qsoc=-Hsoc*(-isc*soc*asoc*exp(-Esoc/(kb*TT)));
%% calculations for the convection
%constants
hconv=7.17;
asurf=3.5E-3;

%equation
Qconv=hconv*asurf*(TT-T0);
%% calculations for the radiation
%constants 
eff=0.8;
sig=5.67E-8;
%equations
tpo=TT^4;
tpo1=T0^4;
Qrad=eff*sig*asurf*(tpo-tpo1);
%% calculations for the total heat generation
% constants
vcell=1.663e-5;
rojr=2580;
cpjr=830;
if(TT<443.15)
    Qheater=29;
else
    Qheater=0;
end
%equations
prod=vcell*rojr*cpjr;
fval(6,1)=(Qsei+Qne+Qpe+Qsoc+Qheater-Qconv-Qrad)/prod;
end