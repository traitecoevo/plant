#params <- function(){

#canopy shape parameters from Yokozawa et al 1995
p.eta=12;

#ratio leaf area to sapwood area
p.theta=4669;
p.theta=10000;


#height  - leaf mass scaling
p.a1=5.44;
p.B1=0.306;

#root leaf scaling
p.a3=0.07;

#scaling of leaf turnover(/yr) to LMA
p.a4=0.4565855*0.5;  #TROPICAL RATe
p.B4=1.71;

p.b = 0.17;

#nitrogen concentrations & photosynthesis
p.n_area=1.87E-3;    #leaf kg/m2
p.c_p1=150.36;
p.c_p2=0.19;

#respiration rates
p.c_Rl = 2.1E4; #mol / kg / yr
p.c_Rs = 4012;  #mol / m3 / yr
p.c_Rr = 217;   #mol / kg / yr
p.c_Rb = 2*p.c_Rs;

#carbon conversion parameter
p.Y = 0.7;
p.c_bio= 12E-3/0.49;

#turnover
p.k_b = 0.2;
p.k_r = 1.0;
p.k_s = 0.2;

#REPRODUCTION
p.c_r1=1;
p.c_r2=50;
p.c_acc=4;
