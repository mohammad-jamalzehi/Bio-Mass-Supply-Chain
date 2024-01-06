sets
t set of all periods from 1 to 365 where each time period coreesponds to one day/1*365/
b set of all biomass available/wc,wp,bq/

parameters
C(b) cost of biomass type on dry basis ($dry tonne)/wc 20,wp 90,bq 136/
Cg Cost of heat produced from natural gas ($GJ)/37/
HD(t) Heat demand from the campus on day t/1*365 12624/
EE(b) gasifictaion existing system/wc 0.735,wp 0.763,bq 0.761/
EN(b) gasifictaion new system/wc 0.768,wp 0.811,bq 0.807/
N Maximum number of truckloads received per day/6/
S Maximum storage capacity (in m3)/412/
Sp Maximum storage capacity of the silo (in m cube) for storing wood pellets/192/
VC(b) Vehicle capacity for biomass in terms of volume/wc 110,wp 110,bq 110/
HCE Maximum heat capacity of the existing system in GJ/518.4/
HCN Maximum heat capacity of the new system in GJ/1555.2/
GHGwood Emission ratio for wood in the existing system (tonnes of CO2-eq. per GJ of heat)/0.96/
GHGng Emission ratio for natural gas (tonnes of CO2-eq. per GJ of heat)/0.5/
p carbon price per tonne Co2-eq in carbon pricing policies/50/
InitialAllowance Initial emission allowance tonnes of CO2-eq allocated in the carbon cap-and-trade policy/2000/
EIApositive Positive deviation of total emissions from the initial allowance./1/
EIAnegative Negative deviation of total emissions from the initial allowance./1/
ComplianceTarget Emissions compliance target (tonnes of CO2-eq.) assigned in the carbon offset policy/2000/
ECTPositive deviation of total emissions from the compliance target. This is the number of carbon offsets to be purchased/1/
ECTnegative Negative deviation of total emissions from the compliance target/1/
table HV(b,t) Heating value of biomass b on day t(GJ dry tonne)
         1*365
wc       12.5

wp       18

bq       19

table MC(b,t) Moisture content of biomass b on day t  (% weight)
         1*365
wc       0.3

wp       0.1

bq       0.12

table D(b,t) Density of biomass b at period t (delivered tonne m^3)
         1*365
wc       0.15

wp       0.75

bq       1.2
;

variables
z1 Total feedstock cost corresponding to a feasible solution x of the model
z2 Total emissions corresponding to a feasible solution x of the model
z3 Bi-objective optimization model with carbon tax policy

positive variables
xE(b,t) Quantity of biomass type b(m3) burned on day t T in the existing system from what is received on dayt
xN(b,t) Quantity of biomass type b (m3) burned on day t T in the new system from what is received on dayt
xES(b,t) Quantity of biomass type b (m3) burned on day t T in the existing system from what is stored on dayt 1
xNS(b,t) Quantity of biomass type b (m3) burned on day t T in the new system from what is stored on dayt 1
sb(b,t) Quantity of biomass of type b stored at the end of day t T (m3)
g(t) Quantity of energy required to be produced from natural gas on day t T (GJ)

integer variables
nt(b,t) Number of truckloads of biomass b received on dayt T;

equations
objectivefunction1
objectivefunction2
co1(t)
co2(b,t)
co3(b,t)
co4(b,t)
co5(t)
co6(t)
co7(b,t)
co8(t)
co9(t);


objectivefunction1       ..z1 =e= sum((t,b),(nt(b,t)*VC(b)*D(b,t)*(1-MC(b,t))*C(b)))+sum(t,g(t)*Cg);
objectivefunction2       ..z2 =e= sum((t,b),D(b,t)*(1-MC(b,t))*HV(b,t)*GHGwood *((xE(b,t)*EE(b))+(xN(b,t)*EN(b))))+sum((t,b),D(b,t-1)*(1-MC(b,t-1))*HV(b,t-1)*GHGwood*((xES(b,t)*EE(b))+(xNS(b,t)*EN(b)))+(g(t)*GHGng));
co1(t)                   ..sum(b,xE(b,t)*(1-MC(b,t))*HV(b,t)*EE(b))+sum(b,xN(b,t)*D(b,t)*(1-MC(b,t))*HV(b,t)*EN(b)) + sum(b,xES(b,t)*D(b,t-1)*(1-MC(b,t-1))*HV(b,t-1)*EE(b)) + sum(b,xNS(b,t)*D(b,t-1)*(1-MC(b,t-1))*HV(b,t-1)*EN(b))+ g(t)=e=HD(t);
co2(b,t)                 ..xNS(b,t) + xES(b,t)=e=sb(b,t);
co3(b,t)                 ..nt(b,t)*VC(b) - xE(b,t) - xN(b,t)=e=sb(b,t);
co4(b,t)                 ..sb(b,t-1)+nt(b,t)*VC(b)-xE(b,t)-xN(b,t)=e=sb(b,t);
co5(t)                   ..sum(b,nt(b,t))=l=N ;
co6(t)                   ..sum(b,sb(b,t))=l=S;
co7(b,t)                 ..sb(b,t)=l=Sp;
co8(t)                   ..sum(b,xN(b,t)*D(b,t)*(1-MC(b,t))*HV(b,t)*EN(b))+sum(b,xNS(b,t)*D(b,t-1)*(1-MC(b,t-1))*HV(b,t-1)*EN(b))=l=HCN;
co9(t)                   ..sum(b,xE(b,t)*D(b,t)*(1-MC(b,t))*HV(b,t)*EE(b))+sum(b,xES(b,t)*D(b,t-1)*(1-MC(b,t-1))*HV(b,t-1)*EE(b))=l=HCE;

Model FirstObjective /objectivefunction1,co1,co2,co3,co4,co5,co6,co7,co8,co9/;
Option Mip=Cplex,optca=0,Optcr=0;
solve FirstObjective using MIP min z1;


Model SecondObjective /objectivefunction2,co1,co2,co3,co4,co5,co6,co7,co8,co9/;
Option Mip=Cplex,optca=0,Optcr=0;
solve SecondObjective using MIP min z2;


display z1.l,z2.l,xE.l,xN.l,xES.l,xNS.l,sb.l,g.l,nt.l;
