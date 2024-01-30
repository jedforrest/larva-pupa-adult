
restart:
with(LinearAlgebra):
with(Statistics):
with(Student[NumericalAnalysis]):
with(VectorCalculus):
# In this implementation, the parameters in the LPA system were renamed as follows:
# a=1-\mu_L, b=b, c=1-\mu_A, d=c_{EL}, f=c_{EA}, g=c_{PA}
# 
# Below are the 'answer' parameter values, which we are trying to reobtain.
aans:= 1-0.2055;
bans:= 6.598;
cans:= 1-7.629*10^(-3);
dans:= 1.209*10^(-2);
fans:= 1.155*10^(-2);
gans:= 4.7*10^(-3);
# In the full code, parameter intervals are randomly generated. However this code matches the example in the paper, hence, we have fixed the values from one particular iteration of the experiment.
# 
# Since our application will need the endpoints, this is how we define the intervals.
amin := 0.7535053327;amax := 0.8354946674;bmin := 6.266830333;bmax := 6.929169668;cmin := 0.9414827824;cmax := 1.043259218;dmin := 0.01021583228;dmax := 0.01396416777;fmin := 0.009702832180;fmax := 0.01339716787;gmin := 0.003195332091;gmax := 0.006204668001;

# We generate simulated population (observed) data from fixed initial populations and using the model without noise, as well as the answer parameters. In this experiment we rounded up.
L(0):= 107 ;
P(0):= 73 ;
A(0):= 214 ;
for t from 0 to 5 do
L(t+1):=ceil(bans*A(t)*exp(-dans*L(t)-fans*A(t)));
P(t+1):=ceil(L(t)*aans);
A(t+1):=ceil(P(t)*exp(-gans*A(t))+A(t)*cans);
end do
;
# Next we use the population data and parameter intervals to find the expansion intervals. This is the algorithm Expansion Ranges
for t from 0 to 6 do
RLmin(t):=-dmax*L(t)-fmax*A(t): RLmax(t):=-dmin*L(t)-fmin*A(t):
RAmin(t):=-gmax*A(t): RAmax(t):=-gmin*A(t):
end do
;
# In the LPA model we use midpoints as the expansion point so we compute these next.
for t from 0 to 6 do
expansionL(t):=(RLmax(t)+RLmin(t))/2:
expansionA(t):=(RAmax(t)+RAmin(t))/2:
end do
;
# Next we need our system and its prolongations. For this implementation, we computed 6 time prolongations for each population L(t), P(t), A(t), although this turned out to be redundant as the system has solutions by the third prolongation.
# 
# Note also the initial populations are just numerical and have no data on the parameters, hence, are of no interest here.
Lequations:=[];
for i from 1 to 6 do
TaylorPolynomial(exp(x),x=expansionL(i),order=2):
EquationL:=L(i)=b*A(i-1)*subs(x=-d*L(i-1)-f*A(i-1),%):
Lequations:=[op(Lequations),EquationL]
end do
;
Aequations:=[];
for i from 1 to 6 do
TaylorPolynomial(exp(x),x=expansionA(i),order=2):
EquationA:=A(i)=P(i-1)*subs(x=-g*A(i-1),%)+A(i-1)*c:
Aequations:=[op(Aequations),EquationA]
end do
;
Pequations:=[];
for i from 1 to 6 do
EquationP:=P(i)=a*L(i-1):
Pequations:=[op(Pequations),EquationP]
end do
;
# For illustration purpose, we show the first three L(t) and A(t) equations as used in the paper (but not rounded) out of these lists:
Lequations[1];Lequations[2];Lequations[3];Aequations[1];Aequations[2];Aequations[3];
# As Maple is not the main implementation project, we do not have automated implementation of the Square System algorithm. We will demonstrate the Jacobian establishing the need to replace the P(2) equation as in the paper.
ReducedRowEchelonForm(Jacobian([(rhs(Aequations[1])-lhs(Aequations[1])),(rhs(Pequations[1])-lhs(Pequations[1])),(rhs(Lequations[1])-lhs(Lequations[1])),(rhs(Aequations[2])-lhs(Aequations[2])),(rhs(Pequations[2])-lhs(Pequations[2])),(rhs(Lequations[2])-lhs(Lequations[2]))],[a,b,c,d,f,g]))
;
# Finally we demonstrate the answer as in the paper:
solve({Aequations[1],Aequations[2],Pequations[1],Lequations[1],Lequations[2],Lequations[3]},[a,b,c,d,f,g])
;

