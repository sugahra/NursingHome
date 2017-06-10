const decl aday="0608";
//const decl iburn=10000;
//const decl irepeat=100000;
const decl iburn=0;
const decl irepeat=100000;
const decl iK_d=7; // number of explanatory variables for demand, including F and f
const decl iK_s=6; // number of explanatory variables for supply
const decl ivars_file=17;
const decl iK_gamma_T=2; // number of T in gamma
const decl iK_gamma_f=2; // number of f in gamma
const decl iK_gamma=iK_gamma_T+iK_gamma_f+2;
const decl iM=38; //# Markets
const decl iH=1265; // Total # of homes 10*20+5*20+20*10, 

const decl in_param=iK_d+iK_s+iK_gamma+3;
//beta_d, beta_s, gamma, alpha, sigma2_omega, sigma2_xi

//Hyperparameters

const decl c_dmu_beta_d0=0;
const decl c_dmu_beta_s0=0;
const decl c_dmu_gamma0=0;
const decl c_dmu_logalpha0=0;
const decl c_dmu_alpha0=1;

const decl c_dsigma2_beta_d0=1000;
const decl c_dsigma2_beta_s0=1000;
const decl c_dsigma2_gamma0=10;
const decl c_dsigma2_logalpha0=10;
const decl c_dsigma2_alpha0=10;
const decl c_dalpha_tau_omega_10=3;
const decl c_dalpha_tau_omega_20=10;
const decl c_dalpha_tau_xi_10=3;
const decl c_dalpha_tau_xi_20=10;


//Apaptive MH
const decl c_dweight=0.05; 
//beta in Rosenthal's papers. Common for all parameters

//Polya trees
// const decl iJ=10; */
// const decl iN=1024;  */
// const decl iL=2046; */
 const decl iJ=5; 
 const decl iN=32; 
 const decl iL=62; 
//iJ~=round(log_2 (iN)). log_x(Y)=log_e(Y)/log_e(x)
//iM=2^J, # distinct 
//iL=# distinct values of \bm{epsilon}= 2+4+8+...+iM
const decl c_dc=10;

//prediction
const decl iR=100;
