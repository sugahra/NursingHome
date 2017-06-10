#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<maximize>
#import<solvenle>
#include"myconst.h"

static decl s_vbeta_d,s_dalpha, s_vbeta_s,s_vgamma,
  s_dtau_omega,s_dtau_xi,
  s_vq,s_vp,s_mX,s_mW,s_vf,s_vT,s_mTf,s_vGamma_T,s_vs,
  s_ck,
  s_vmu_beta_d0,s_vmu_beta_s0,s_vmu_gamma0,
  s_minv_Sigma_beta_d0,s_minv_Sigma_beta_s0,s_minv_Sigma_gamma0,
  s_vh_m,s_vj,s_mijj2,
  s_mX_logit,s_varname,s_vconst,
  s_vomega,s_vxi,s_vlog_one_minus_sum_s,s_cloop;

//Functions Section
fn_varname(){
  s_varname=new array[in_param];
  s_varname[0]="$\\beta_{d}$: Constant";
  s_varname[1]="$\\beta_{d}$: \# Residents per worker";
  s_varname[2]="\\beta_{d}: Years from opening";
  s_varname[3]="\\beta_{d}: Occupancy rate";
  s_varname[4]="\\beta_{d}:Chain";
  s_varname[5]="\\alpha_{F}";
  s_varname[6]="\\alpha_{f}";
  s_varname[7]="\\gamma_{0}";    
  s_varname[8]="\\gamma_{T1}";	     
  s_varname[9]="\\gamma_{T2}";
  s_varname[10]="\\gamma_{f1}";    
  s_varname[11]="\\gamma_{f2}";	     
  s_varname[12]="\\gamma_{fT}";	     
  s_varname[13]="\\alpha";     
  s_varname[14]="\\beta_{s}:Constant"; 
  s_varname[15]="\\beta_{s}:Worker_Ratio"; 
  s_varname[16]="\\beta_{s}:Years from opening";    
  s_varname[17]="\\beta_{s}:Chain";
  s_varname[18]="\\beta_{s}:Local_Average_Rent";  
  s_varname[19]="\\beta_{s}:Local_Average_Wage"; 
  s_varname[20]="\\tau_{\\omega}"; 
  s_varname[21]="\\tau_{\\xi}"; 
  return 1;
}


fn_Generate_Q(const Func, const vq)
{
  decl vfunc,dsumexpq,vs;
  vs=exp(vq);
  dsumexpq=sumc(vs);
  vfunc=1+(vs-1).* (-vq +log(1-dsumexpq)+s_vconst);
  Func[0]=vfunc;
  return 1;
}

fn_Gamma_T(const vgamma)
{
  decl vTgamma,vexpTgamma,vGamma_T,ch;
  vTgamma=s_mTf*vgamma;
  vexpTgamma=exp(vTgamma);
  vGamma_T=vexpTgamma./(1+vexpTgamma);
  for(ch=0;ch<iH;ch++)
    {
      if(s_mTf[ch][1]==0)
	{
	  vGamma_T[ch]=0;
	}
      else if( isdotinf(vexpTgamma[ch])==1  && vTgamma[ch]>0)
	{
	  vGamma_T[ch]=1;
	}
      else if(isdotinf(vexpTgamma[ch])==1  && vTgamma[ch]<0)
	{

	  vGamma_T[ch]=0;
	}
    }
      return vGamma_T;
}


fn_randtrann_R(B)
{
  // Generate a right-truncated normal r.v.
  decl mu;
  mu=ranu(sizer(B),sizec(B)).*(probn(B));
  return quann(mu);
}

fn_randtrann(B,iP)
{
  // Truncated Normal: 
    //The first argment(matrix) is the truncation point, 
    //The second argment(int) is the truncation oriantation:
    // iP=1: Right truncation (generate an r.v. smaller than B)
    // iP=3: Left truncation (generate an r.v. larger than B)
    // iP=2 must be the both truncation(not yet coded) 
  decl mx;
  if(iP==1)
    {return fn_randtrann_R(B);}
  if(iP==3)
    {
      mx=fn_randtrann_R(-B);
      return B+fabs(mx-(-B));
    }
}

fn_polya_floor(const dtau, const vW)
{
  decl mW,mW_floor,mn_index, mn,mf_num,mf_den,ci,cj,cl,
    isum, df;
 
  mW=probn(sqrt(dtau)*vW)*s_vj;
  mW_floor=floor(mW);
  mn_index=zeros(iH,iL);
  mn=zeros(iH,iJ);
  mf_num=mf_den=zeros(iH-1,iJ);


  //ci=0
  //cj=0;
  isum=0;
  for(cj=0;cj<iJ;cj++)
    {	  
      cl=isum+mW_floor[0][cj];
      if(cl==isum+2^(cj+1)){cl=cl-1;}
      mn_index[0][cl]=1;
      isum+=2^(cj+1);
    }

  //end of ci=0

  for(ci=1;ci<iH;ci++)
    {
      isum=0;
      cl=isum+mW_floor[ci][0];
      if(cl==2){cl=cl-1;}
      mn_index[ci][cl]=1;
      mn[ci][0]=sumc(mn_index[0:ci-1][cl]);
      mf_num[ci-1][0]=s_mijj2[0][0]+mn[ci][0];
      mf_den[ci-1][0]=s_mijj2[0][1]+ci;
      isum+=2;

      for(cj=1;cj<iJ;cj++)
	{
	  cl=isum+mW_floor[ci][cj];
	  if(cl==isum+2^(cj+1)){cl=cl-1;}
	  mn_index[ci][cl]=1;
	  mn[ci][cj]=sumc(mn_index[0:ci-1][cl]);
	  mf_num[ci-1][cj]=s_mijj2[cj][0]+mn[ci][cj];
	  mf_den[ci-1][cj]=s_mijj2[cj][1]+mn[ci][cj-1];
	  isum+=2^(cj+1);
	}
    }
  return sumc(sumr(log(mf_num)-log(mf_den)));
}

fn_omega(const dalpha, const vbeta_s, const vGamma_T)
{
  return log(s_vp+ s_vf.*vGamma_T+(ones(iH,1)./(dalpha*(s_vs-1)))) 
    -s_mW*vbeta_s;
}

fn_xi(const dalpha, const vbeta_d )
{
  return s_vq -s_mX*vbeta_d + dalpha*s_vp-s_vlog_one_minus_sum_s;
}

fn_sumZinv(const vGamma_T, const dalpha)
{
  return sumc(log(s_vp+s_vf.*vGamma_T+(ones(iH,1)./ (dalpha*(s_vs-1)))));
}

fn_logdensity_normal(const dsigma2, const vx)
{
  return -0.5*iH*log(dsigma2)-0.5*sumsqrc(vx)/dsigma2;
}

fn_beta_d(const Beta, const Func, const Score, const Hess)
   {
     decl vbeta,dlikelihood,vxi,dpolya_xi,dg_xi;
     vbeta=s_vbeta_d;
     vbeta[s_ck]=Beta;
     vxi=fn_xi(s_dalpha, vbeta);
     dpolya_xi=fn_polya_floor(s_dtau_xi,vxi);
     dg_xi=fn_logdensity_normal((1/s_dtau_xi), vxi);
     dlikelihood=dpolya_xi+dg_xi;
     Func[0]=-0.5 *(vbeta-s_vmu_beta_d0)'*s_minv_Sigma_beta_d0
       *(vbeta-s_vmu_beta_d0)
       +dlikelihood;
     return 1;
   }

fn_beta_s(const Beta, const Func, const Score, const Hess)
   {
     decl vbeta,dlikelihood,vomega,dpolya_omega,dg_omega;
     vbeta=s_vbeta_s;
     vbeta[s_ck]=Beta;
     vomega=fn_omega(s_dalpha, vbeta, s_vGamma_T);
     dpolya_omega=fn_polya_floor(s_dtau_omega, vomega);
     dg_omega=fn_logdensity_normal((1/s_dtau_omega), vomega);
     dlikelihood=dpolya_omega+dg_omega;
     Func[0]=-0.5 *(vbeta-s_vmu_beta_s0)'*s_minv_Sigma_beta_s0
       *(vbeta-s_vmu_beta_s0)
       +dlikelihood;
     return 1;
   }
fn_gamma(const Gamma, const Func, const Score, const Hess)
   {
     decl vgamma,dlikelihood,vGamma_T,dsumZinv,vomega,dpolya_omega,dg_omega;
     vgamma=s_vgamma;
     vgamma[s_ck]=Gamma;
     vGamma_T=fn_Gamma_T(vgamma);
     dsumZinv=fn_sumZinv(vGamma_T,s_dalpha);
     vomega=fn_omega(s_dalpha, s_vbeta_s, vGamma_T);
     // if(s_cloop==23)
     //   {
     // 	 println(vGamma_T~vomega);
     //   }

     dpolya_omega=fn_polya_floor(s_dtau_omega, vomega);
     dg_omega=fn_logdensity_normal((1/s_dtau_omega), vomega);
     dlikelihood=-dsumZinv+dpolya_omega+dg_omega;
     Func[0]=-0.5 *(vgamma-s_vmu_gamma0)'*s_minv_Sigma_gamma0
       *(vgamma-s_vmu_gamma0)
       +dlikelihood;
     return 1;
   }

fn_alpha(const LogAlpha, const Func, const Score, const Hess)
   {
     decl dlogalpha,dalpha,dlikelihood,dsumZinv,
       vomega,dpolya_omega,dg_omega,vxi,dpolya_xi,dg_xi;
     dlogalpha=LogAlpha;
     dalpha=exp(dlogalpha);
     vomega=fn_omega(dalpha, s_vbeta_s, s_vGamma_T);
     vxi=fn_xi(dalpha, s_vbeta_d);
     dpolya_omega=fn_polya_floor(s_dtau_omega, vomega);
     dpolya_xi=fn_polya_floor(s_dtau_xi,vxi);
     dg_omega=fn_logdensity_normal((1/s_dtau_omega), vomega);
     dg_xi=fn_logdensity_normal((1/s_dtau_xi), vxi);
     dsumZinv=fn_sumZinv(s_vGamma_T,dalpha);
     dlikelihood=-dsumZinv+dpolya_omega+dpolya_xi+dg_omega+dg_xi;
     Func[0]=-0.5 * ((dalpha-c_dmu_alpha0)^(2)) /c_dsigma2_alpha0
       +dlikelihood;
     return 1;
   }

fn_tau_omega(const Tau, const Func, const Score, const Hess)
   {
     decl dtau_omega,dlikelihood,dpolya_omega,dg_omega;
     dtau_omega=Tau;
     dpolya_omega=fn_polya_floor(dtau_omega, s_vomega);
     dg_omega=fn_logdensity_normal((1/dtau_omega), s_vomega);
     dlikelihood=dpolya_omega+dg_omega;
     Func[0]=(c_dalpha_tau_omega_10-1)*log(dtau_omega)
       -c_dalpha_tau_omega_20*dtau_omega
       +dlikelihood;
     return 1;
   }

fn_tau_xi(const Tau, const Func, const Score, const Hess)
   {
     decl dtau_xi,dlikelihood,dpolya_xi,dg_xi;
     dtau_xi=Tau;
     dpolya_xi=fn_polya_floor(dtau_xi,s_vxi);
     dg_xi=fn_logdensity_normal((1/dtau_xi), s_vxi);
     dlikelihood=dpolya_xi+dg_xi;
     Func[0]=(c_dalpha_tau_xi_10-1)*log(dtau_xi)
       -c_dalpha_tau_xi_20*dtau_xi
       +dlikelihood;
     return 1;
   }



main()
{
  decl time,mSigma,dmh,cloop,file,ck,
    dold,dnew,dtarget_old,dtarget_new,
    mbeta_s_sample, mbeta_d_sample, mgamma_sample,
    valpha_sample, vtau_omega_sample, vtau_xi_sample,
    msample,
    vmh_beta_d, vmh_beta_s, vmh_gamma,
    imh_alpha,imh_tau_xi,imh_tau_omega,
    mdata,mdata_market,cj,
    vy_s,dopt,maxi,
    dlog_L_alpha,vL_G,vlogL_G_frac,dL_gamma,vL_gamma,
    vf_nozero,vindex_nozero_f,ch,vindex_L_G,vL_G_one,
    du, vindex_gamma,vT_nozero,dproposal_old,dproposal_new,
    vtheta,vTgamma,dx,solve,hess,vGamma_old,vgamma_new,
    vmean,mcredible_01,mcredible_05,mcredible_10,astar,vse,vmh,
    vone_minus_s,vlog_one_minus_s,vp_f_nozero,vone_minus_s_f_nozero,
    ih_past,cm,iHm,vchain,mTf_nozero;

  //Preset
  time=timer();
  file=fopen("../data/market_index0918.dat","r");
  fscan(file,"%#m",iM,2,&mdata_market);
  fclose(file);

  file=fopen("../data/nursing0918.dat");
  fscan(file,"%#m",iH,ivars_file,&mdata);
  fclose(file);

  s_vp=mdata[][1];
  s_vq=mdata[][2];
  s_vs=exp(s_vq);
  vchain=mdata[][7]+mdata[][8]+mdata[][9]+mdata[][10]+mdata[][11]+mdata[][12];
  s_mX=ones(iH,1)~mdata[][3:5]~vchain;
  s_mW=ones(iH,1)~mdata[][3:4]~vchain~mdata[][13:14];
  s_vf=mdata[][15]/12;
  s_vT=mdata[][16]*12;
  s_mX=s_mX~((s_vf.*s_vT)/1000)~s_vf;
  s_vh_m=mdata_market[][1];
  s_mX[][2]=(s_mX[][2]-meanc(s_mX[][2]))/sqrt(varc(s_mX[][2]));
  s_mW[][2]=(s_mW[][2]-meanc(s_mW[][2]))/sqrt(varc(s_mW[][2]));
  s_mW[][iK_s-2]=(s_mW[][iK_s-2]-meanc(s_mW[][iK_s-2]))
    /sqrt(varc(s_mW[][iK_s-2]));
  s_mW[][iK_s-1]=(s_mW[][iK_s-1]-meanc(s_mW[][iK_s-1]))
    /sqrt(varc(s_mW[][iK_s-1]));
  s_vlog_one_minus_sum_s=zeros(iH,1);
  ih_past=0;
  for(cm=0;cm<iM;cm++)
    {
      iHm=s_vh_m[cm];      
      s_vlog_one_minus_sum_s[ih_past:ih_past+iHm-1]
	=log(1-sumc(s_vs[ih_past:ih_past+iHm-1]))*ones(iHm,1);
      ih_past+=iHm;
    }

  vf_nozero=s_vf;
  vT_nozero=s_vT;
  vindex_nozero_f=ones(iH,1);

  for(ch=0;ch<iH;ch++)
    {
      if(s_vT[ch]==0 || s_vf[ch]==0)
	{
	  s_vT[ch]=s_vf[ch]=0;
	  vT_nozero[ch]=1;
	  vf_nozero[ch]=1;
	  vindex_nozero_f[ch]=0;
	}
    }

  s_mTf=ones(iH,1);  
  mTf_nozero=ones(iH,1);
  for(ck=1;ck<iK_gamma_T+1;ck++)
    {
      s_mTf=s_mTf~(s_vT.^(ck));
      mTf_nozero=mTf_nozero~(vT_nozero.^(ck));
    }
  for(ck=1;ck<iK_gamma_f+1;ck++)
    {
      s_mTf=s_mTf~(s_vf.^(ck));
      mTf_nozero=mTf_nozero~(vf_nozero.^(ck));
    }
  s_mTf=s_mTf~(s_vf.*s_vT);
  mTf_nozero=mTf_nozero~(vf_nozero.*vT_nozero);

  s_mTf[][1]=s_mTf[][1]/10;
  s_mTf[][2]=s_mTf[][2]/1000;
  s_mTf[][3]=s_mTf[][3]/10;
  s_mTf[][4]=s_mTf[][4]/100;
  s_mTf[][5]=s_mTf[][5]/1000;

  mTf_nozero[][1]=mTf_nozero[][1]/10;
  mTf_nozero[][2]=mTf_nozero[][2]/1000;
  mTf_nozero[][3]=mTf_nozero[][3]/10;
  mTf_nozero[][4]=mTf_nozero[][4]/100;
  mTf_nozero[][5]=mTf_nozero[][5]/1000;

  s_vj=zeros(iJ,1);
  for(cj=0;cj<iJ;cj++)
    {
      s_vj[cj]=2^(cj+1);
    }
  s_vj=s_vj';
  s_mijj2=zeros(iJ,2);
  for(cj=0;cj<iJ;cj++)
    {
      s_mijj2[cj][0]=c_dc*(cj+1)*(cj+1);
    }
  s_mijj2[][1]=2*s_mijj2[][0];

  //Set sizes of matrices
  mbeta_d_sample=zeros(iK_d,irepeat);
  mbeta_s_sample=zeros(iK_s,irepeat);
  mgamma_sample=zeros(iK_gamma,irepeat);
  valpha_sample=zeros(irepeat,1);
  vtau_omega_sample=zeros(irepeat,1);
  vtau_xi_sample=zeros(irepeat,1);

  //Hypermarameters
  s_vmu_beta_d0=c_dmu_beta_d0*ones(iK_d,1);
  s_vmu_beta_s0=c_dmu_beta_s0*ones(iK_s,1);
  s_vmu_gamma0=c_dmu_gamma0*ones(iK_gamma,1);
  s_minv_Sigma_beta_d0=(1/c_dsigma2_beta_d0)*unit(iK_d);
  s_minv_Sigma_beta_s0=(1/c_dsigma2_beta_s0)*unit(iK_s);
  s_minv_Sigma_gamma0=(1/c_dsigma2_gamma0)*unit(iK_gamma);

  //Initial values
  // vtheta=<-4.30008,  -0.74125,   0.03405,   0.67931,   0.35624,   0.05144,  0, 0.42536,  
  // -2.73817,   0.55788,  -3.83525,  -0.80675,   4.12104,   0.13906,   2.73134,  
  // -0.32992,  -0.09607,   0.16706,   0.21921,  -0.03727,   3.03464,   0.90808>';
  vtheta=<-4.26406,     -0.75568, -0.01696,  0.73783,  0.38768,       0.01027,   0.00262, -0.96932, 
    -2.43534, -1.02059,    -1.21184, -4.38132,   4.83881,   0.13849,  3.05343, 
    -0.46827,  -0.10023,  0.24151,  0.18374,  0.01514,  2.52629,  0.98116>';
  s_vbeta_d=vtheta[0:iK_d-1];
  s_vgamma=vtheta[iK_d:iK_d+iK_gamma-1];
  s_dalpha=vtheta[iK_d+iK_gamma];
  s_vbeta_s=vtheta[iK_d+iK_gamma+1:iK_d+iK_gamma+1+iK_s-1];
  s_dtau_omega=vtheta[iK_d+iK_gamma+iK_s+1];
  s_dtau_xi=vtheta[iK_d+iK_gamma+iK_s+2];

  s_vGamma_T=fn_Gamma_T(s_vgamma);
  vone_minus_s=(1-s_vs);
  vlog_one_minus_s=log(vone_minus_s);
  vone_minus_s_f_nozero=vone_minus_s.*vf_nozero;
  vp_f_nozero=s_vp./vf_nozero;

  dlog_L_alpha=max(-log(s_vf.*s_vGamma_T+s_vp)-vlog_one_minus_s);
  vL_G=ones(iH,1)./(s_dalpha*vone_minus_s_f_nozero) - vp_f_nozero;
  vindex_L_G=(vL_G.>0);
  vindex_gamma=(minc((vindex_nozero_f~vindex_L_G)'))';
  vlogL_G_frac=log((vL_G./(1-vL_G)).^vindex_gamma);

  //Set to null

  decl    c_vsigma_beta_d_rw,  c_vsigma_beta_s_rw, c_vsigma_gamma_rw,
  c_dsigma_logalpha_rw,c_dsigma_logtau_omega_rw,c_dsigma_logtau_xi_rw;
  // decl c_vsigma_beta_d_rw=0.1*ones(iK_d,1);
  // decl c_vsigma_beta_s_rw=0.1*ones(iK_s,1);

  c_vsigma_beta_d_rw=0.8*ones(iK_d,1);
  c_vsigma_beta_d_rw[0]=0.045;
  c_vsigma_beta_d_rw[1]=0.015;
  c_vsigma_beta_d_rw[2]=0.05;
  c_vsigma_beta_d_rw[3]=0.055;
  c_vsigma_beta_d_rw[4]=0.082;
  c_vsigma_beta_d_rw[5]=0.031;
  c_vsigma_beta_d_rw[6]=0.003;
  c_vsigma_gamma_rw=<0.4;0.3;0.5;0.7;1;0.9>;
  c_vsigma_gamma_rw[0]=4;
  c_vsigma_gamma_rw[1]=2.5;
  c_vsigma_gamma_rw[2]=3.2;
  c_vsigma_gamma_rw[3]=4;
  c_vsigma_gamma_rw[4]=3;
  c_vsigma_gamma_rw[5]=3.9;
  c_dsigma_logalpha_rw=0.002;
  c_vsigma_beta_s_rw=0.001*ones(iK_s,1);
  c_vsigma_beta_s_rw[0]=0.003;
  c_vsigma_beta_s_rw[1]=0.0015;
  c_vsigma_beta_s_rw[2]=0.007;
  c_vsigma_beta_s_rw[3]=0.015;
  c_vsigma_beta_s_rw[4]=0.0065;
  c_vsigma_beta_s_rw[5]=0.0065;
  c_dsigma_logtau_omega_rw=0.059;
  c_dsigma_logtau_xi_rw=0.089;

  vmh_beta_d=zeros(iK_d,1);
  vmh_beta_s=zeros(iK_s,1);
  vmh_gamma=zeros(iK_gamma,1);
  imh_alpha=imh_tau_xi=imh_tau_omega=0;

  //Iteration start
  for(cloop=-iburn;cloop<irepeat;cloop++)
    {
      println(cloop);
      //println(s_vbeta_d'~s_vgamma'~s_dalpha~s_vbeta_s'~s_dtau_omega~s_dtau_xi);
      //For demand side
      //Sample beta_d
      s_cloop=cloop;
     for(s_ck=0;s_ck<iK_d;s_ck++)
       {
	 dold=s_vbeta_d[s_ck];
	 fn_beta_d(dold,&dtarget_old,0,0);
	 dnew=dold+c_vsigma_beta_d_rw[s_ck]*rann(1,1);
	 fn_beta_d(dnew,&dtarget_new,0,0);
	 dmh=exp(dtarget_new-dtarget_old);
	 if(ranu(1,1)<=dmh)
	   {     
	     s_vbeta_d[s_ck]=dnew;
	     if(cloop>=0){vmh_beta_d[s_ck]++;}
	   }
       }

     //Sample alpha
     dold=log(s_dalpha);
     fn_alpha(dold,&dtarget_old,0,0);
     dnew=dold+c_dsigma_logalpha_rw
       *fn_randtrann((dlog_L_alpha-dold)/c_dsigma_logalpha_rw,3);
     fn_alpha(dnew,&dtarget_new,0,0);
     dproposal_old=-dold-log(probn((dnew-dlog_L_alpha)/c_dsigma_logalpha_rw));
     dproposal_new=-dnew-log(probn((dold-dlog_L_alpha)/c_dsigma_logalpha_rw));
     dmh=exp(dtarget_new+dproposal_old-dtarget_old-dproposal_new);
     if(ranu(1,1)<=dmh)
       {     
	 s_dalpha=exp(dnew);
	 if(cloop>=0){imh_alpha++;}
	 vL_G=ones(iH,1)./(s_dalpha*vone_minus_s_f_nozero) - vp_f_nozero;
	 vindex_L_G=(vL_G.>0);
	 vindex_gamma=(minc((vindex_nozero_f~vindex_L_G)'))';
	 vlogL_G_frac=log((vL_G./(1-vL_G)).^vindex_gamma);
       }

     //Sample gamma
     for(s_ck=0;s_ck<iK_gamma;s_ck++)
       {
	 dold=s_vgamma[s_ck];
	 fn_gamma(dold,&dtarget_old,0,0);
	 vL_gamma=(vlogL_G_frac-s_mTf*s_vgamma+dold*s_mTf[][s_ck])
	   ./(mTf_nozero[][s_ck]).*(vindex_gamma)
	   -100000000*(1-vindex_gamma);
	 dL_gamma=max(vL_gamma);
	 du=fn_randtrann((dL_gamma-dold)/c_vsigma_gamma_rw[s_ck],3);
	 dnew=dold+c_vsigma_gamma_rw[s_ck]*du;
	 //println(dL_gamma~du~dnew);
	 fn_gamma(dnew,&dtarget_new,0,0);
	 dproposal_old=-log(probn((dnew-dL_gamma)/c_vsigma_gamma_rw[s_ck]));
	 dproposal_new=-log(probn((dold-dL_gamma)/c_vsigma_gamma_rw[s_ck]));
	 dmh=exp(dtarget_new+dproposal_old-dtarget_old-dproposal_new);
	 if(ranu(1,1)<=dmh)
	   {     
	     s_vgamma[s_ck]=dnew;
	     if(cloop>=0){vmh_gamma[s_ck]++;}
	     s_vGamma_T=fn_Gamma_T(s_vgamma);
	     dlog_L_alpha=max(-log(s_vf.*s_vGamma_T+s_vp)-vlog_one_minus_s);
	   }
       }
     //Supply side

     //Sample beta_s
     for(s_ck=0;s_ck<iK_s;s_ck++)
       {
	 dold=s_vbeta_s[s_ck];
	 fn_beta_s(dold,&dtarget_old,0,0);
	 dnew=dold+c_vsigma_beta_s_rw[s_ck]*rann(1,1);
	 fn_beta_s(dnew,&dtarget_new,0,0);
	 dmh=exp(dtarget_new-dtarget_old);
	 if(ranu(1,1)<=dmh)
	   {     
	     s_vbeta_s[s_ck]=dnew;
	     if(cloop>=0){vmh_beta_s[s_ck]++;}
	   }
       }

     //Variance components

     //Sample tau_omega

     s_vomega=fn_omega(s_dalpha, s_vbeta_s, s_vGamma_T);
     dold=log(s_dtau_omega);
     fn_tau_omega(exp(dold),&dtarget_old,0,0);
     dnew=dold+c_dsigma_logtau_omega_rw*rann(1,1);
     fn_tau_omega(exp(dnew),&dtarget_new,0,0);
     dproposal_old=-dold;
     dproposal_new=-dnew;
     dmh=exp(dtarget_new+dproposal_old-dtarget_old-dproposal_new);
     if(ranu(1,1)<=dmh)
       {     
	 s_dtau_omega=exp(dnew);
	 if(cloop>=0){imh_tau_omega++;}
       }

     //Sample tau_xi

     s_vxi=fn_xi(s_dalpha, s_vbeta_d);
     dold=log(s_dtau_xi);
     fn_tau_xi(exp(dold),&dtarget_old,0,0);
     dnew=dold+c_dsigma_logtau_xi_rw*rann(1,1);
     fn_tau_xi(exp(dnew),&dtarget_new,0,0);
     dproposal_old=-dold;
     dproposal_new=-dnew;
     dmh=exp(dtarget_new+dproposal_old-dtarget_old-dproposal_new);
     if(ranu(1,1)<=dmh)
       {     
	 s_dtau_xi=exp(dnew);
	 if(cloop>=0){imh_tau_xi++;}
       }
      

     if(cloop>=0)
       {
	 mbeta_d_sample[][cloop]=s_vbeta_d;
	 mbeta_s_sample[][cloop]=s_vbeta_s;
	 mgamma_sample[][cloop]=s_vgamma;
	 valpha_sample[cloop]=s_dalpha;
	 vtau_omega_sample[cloop]=s_dtau_omega;
	 vtau_xi_sample[cloop]=s_dtau_xi;	  
       }
    }       

  //Save Output
  msample=mbeta_d_sample|mgamma_sample|valpha_sample'
    |mbeta_s_sample|vtau_omega_sample'| vtau_xi_sample';
  msample=msample';
  vmh=vmh_beta_d|vmh_gamma|imh_alpha|vmh_beta_s|imh_tau_omega|imh_tau_xi;
  vmh=vmh/irepeat;
  println("time lapsed=",timespan(time));

  format(2000);
  file=fopen("msample.dat","w");
  fprint(file,"%15.5f",msample);
  fclose(file);

  fn_varname();
  vmean=meanc(msample)';
  vse=sqrt(varc(msample)');
  mcredible_01=quantilec(msample,<0.005,0.995>)';
  mcredible_05=quantilec(msample,<0.025,0.975>)';
  mcredible_10=quantilec(msample,<0.05,0.95>)';
  astar=new array[in_param];

  for(ck=0;ck<in_param;ck++)
    {
      if(mcredible_01[ck][0]>0 || mcredible_01[ck][1]<0)
	    {
	      astar[ck]="***";
	    }
      else if(mcredible_05[ck][0]>0 || mcredible_05[ck][1]<0)
	{
	  astar[ck]="**";
	}
      else if(mcredible_10[ck][0]>0 || mcredible_10[ck][1]<0)
	{
	  astar[ck]="*";
	}
      else
	{
	  astar[ck]=" ";
	}
    }

  for(ck=0;ck<in_param;ck++)
    {
      println(s_varname[ck], "%10.3f",vmean[ck] ,astar[ck], 
	      "   (",
	      "%5.3f",vse[ck],")   ", " [ ",mcredible_05[ck][0], "  ,  ", mcredible_05[ck][1]," ] "
	      ,"%10.3f", vmh[ck]);  
    }

}
