#params <- function(){

#shading_spline shape parameters from Yokozawa et al 1995
p.eta=12;

#ratio leaf area to sapwood area
p.theta=1/4669;
p.theta=1/10000;


#height  - leaf mass scaling
p.a_l1=5.44;
p.a_l2=0.306;

#root leaf scaling
p.a_r1=0.07;

#scaling of leaf turnover(/yr) to LMA
p.k_l=0.4565855*0.5;  #TROPICAL RATe

#scaling of stem turnover(/yr) to wood desnity
p.k_s = 0.2;

p.a_b1= 0.17;

#nitrogen concentrations & photosynthesis
p.n_area=1.87E-3;    #leaf kg/m2
p.a_p1=151.177775377968;
p.a_p2=0.204716166503633;

#respiration rates
p.r_l = 2.1e4 * 1.87e-3 / 0.1978791; #mol / kg / yr
p.r_s = 4012/608;  #mol / m3 / yr
p.r_r = 217;   #mol / kg / yr
p.r_b = 2*p.r_s;

#carbon conversion parameter
p.a_y = 0.7;
p.a_bio= 12E-3/0.49;

#turnover
p.k_b = 0.2;
p.k_r = 1.0;

#REPRODUCTION
p.a_f1=1;
p.a_f2=50;
p.a_f3=3*3.8e-5;

# light interception
p.k_I=0.5;
