#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<maximize>
#import<solvenle>
#include"myconst.h"
//15 is 20, 20 is 30. Sorry!

static decl s_vbeta_d,s_dalpha, s_vbeta_s,s_vgamma,
  s_dtau_omega,s_dtau_xi,s_dsigma2_omega,s_dsigma2_xi,
  s_vj,s_mijj2,
  s_vconst,s_dconst,s_vone_minus_sum_s,
  s_vn_xi,s_vn_omega,s_iHm,s_vXB,s_vWB,s_dWB,s_dXB,s_vq_others,
  s_dp,s_dq,s_ds,s_dsigma_omega,s_dsigma_xi,
  s_vq,s_vp,s_mX,s_mW,
  s_vGamma_T,s_vs_all,s_chm,s_dlog_one_minus_sum_s,
  s_vf,s_vT;

//Functions Section

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

fn_xi_prediction_p(const dp)
{
  return s_dq -s_dXB + s_dalpha*dp-s_dlog_one_minus_sum_s;
}

fn_xi_prediction_q(const vq)
{
  return vq-s_vXB+s_dalpha*s_vp-log(1-sumc(exp(vq)));
}

fn_omega_prediction(const dp, const ds)
{
  return log(dp+1/(s_dalpha*(ds-1)))-s_dWB;
}

fn_polya_n(const dtau, const vW)
{
  decl mW, mW_floor, mn_index,cj,cl,ci,
    isum, df, vn,vs,in_all;
  mW=probn(sqrt(dtau)*vW)*s_vj;
  mW_floor=floor(mW);
  mn_index=zeros(iH,iL);

  isum=0;
  for(cj=0;cj<iJ;cj++)
    {	  
      cl=isum+mW_floor[0][cj];
      if(cl==isum+2^(cj+1))
	{cl=cl-1;}
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
      isum+=2;
      for(cj=1;cj<iJ;cj++)
	{
	  cl=isum+mW_floor[ci][cj];
	  if(cl==isum+2^(cj+1)){cl=cl-1;}
	  mn_index[ci][cl]=1;
	  isum+=2^(cj+1);
	}
    }
  vn=(sumc(mn_index))';
  return vn; 
}

fn_predictive_dist(const dW,const dtau, const vn)
{
  decl vprob,idenom,vW,vW_floor,isum,cl,cj;
  vprob=zeros(iJ,1);
  idenom=iH;
  vW=probn(sqrt(dtau)*dW)*s_vj;
  vW_floor=floor(vW);
  isum=0;
  for(cj=0;cj<iJ;cj++)
    {
      cl=isum+vW_floor[cj];
      if(cl==isum+2^(cj+1)){cl=cl-1;}
      vprob[cj]=(s_mijj2[cj][0]+vn[cl]) /(s_mijj2[cj][1]+idenom);
      idenom=vn[cl];
      isum+=2^(cj+1);
    }
  return vprob;
}

fn_detD(vq)
{
  decl ds,done_minus_sums,vs,mD,chm;
  done_minus_sums=1-sumc(exp(vq));
  mD=zeros(s_iHm,s_iHm);
  for(chm=0;chm<s_iHm;chm++)
    {
      ds=exp(vq[chm]);
      mD[][chm]=(ds/done_minus_sums)*ones(s_iHm,1);
      mD[chm][chm]=1 + ds/done_minus_sums + ds/((1-ds)^(2));
    }
  return determinant(mD);
}

fn_logZ(const dp, const ds)
{
  return -log(dp+(1/ (s_dalpha*(ds-1))));
}


fn_q(const dq, const Func, const Score, const Hess)
{
  decl domega,vxi,ds,vdist_xi,vq,chm;
  ds=exp(dq);
  vq=s_vq;
  vq[s_chm]=dq;
  domega=fn_omega_prediction(s_dp, ds);
  vxi=fn_xi_prediction_q(vq);
  vdist_xi=zeros(s_iHm,1);
  for(chm=0;chm<s_iHm;chm++)
    {
      vdist_xi=sumc(log(fn_predictive_dist(vxi[chm],s_dtau_xi, s_vn_xi)));
    }
  Func[0]=log(fn_detD(vq))+fn_logZ(s_dp,ds)
    -0.5*domega*domega/s_dsigma2_omega
    +sumc(log(fn_predictive_dist(domega,s_dtau_omega, s_vn_omega)))
    -0.5*(vxi'*vxi)/s_dsigma2_xi;
    +sumc(vdist_xi);
    return 1;
}

fn_p(const dp, const Func, const Score, const Hess)
{
  decl domega,dxi;
  dxi=fn_xi_prediction_p(dp);
  domega=fn_omega_prediction(dp, s_ds);
  Func[0]=fn_logZ(dp,s_ds)
    -0.5*domega*domega/s_dsigma2_omega
    -0.5*dxi*dxi/s_dsigma2_xi
    +sumc(log(fn_predictive_dist(dxi,s_dtau_xi, s_vn_xi)))
    +sumc(log(fn_predictive_dist(domega,s_dtau_omega, s_vn_omega)));
  return 1;
}

main()
{
  decl time,file,cj,ck,cm,cloop,
    vchain,
    msample,mdata,mdata_market,
    ih_past,iHm,
    vmarket,ch,
    vindex_nozero_f,
    vxi_all,vomega_all,
    mq_sample,mp_sample,vWB,vXB,vq_others,
    mW,chm,
    vmh_q,vmh_p,dold,dnew,dmh,
    vq_all,vp_all,mX_all,mW_all,vf_all,vT_all,mTf_all,
    vGamma_T,vs_all,vlog_one_minus_sum_s_all,vh_m,
    vTgamma,vexpTgamma,
    vsigma_p_rw,vsigma_q_rw,dtarget_old,dtarget_new,mmh,vq_new,
    dU_q,dU_q_one,dU_q_two,dL_p,dproposal_old,dproposal_new,
    vtotalpay,mp_new,vp_new_mean,vtotalpay_new_mean,mprediction,
    iT,vid,mprediction_noT,mprediction_withT,ch_noT,ch_withT,
    vprediction_5,mprediction_5,mprediction_5_withT,itau_5,
    vtotalpay_5,vtotalpay_new_mean_5,
    vprediction_10,mprediction_10,mprediction_10_withT,
    itau_10,vtotalpay_10,vtotalpay_new_mean_10,
    vprediction_15,mprediction_15,mprediction_15_withT,itau_15,
    vtotalpay_15,vtotalpay_new_mean_15,
    vprediction_20,mprediction_20,mprediction_20_withT,itau_20,
    vtotalpay_20,vtotalpay_new_mean_20,
    cp,mmean_p,vstep,cc,mprediction_all,vT_withT;

  //Preset
  time=timer();
  file=fopen("../data/market_index0918.dat","r");
  fscan(file,"%#m",iM,2,&mdata_market);
  fclose(file);

  file=fopen("../data/nursing0918.dat");
  fscan(file,"%#m",iH,ivars_file,&mdata);
  fclose(file);

  vp_all=mdata[][1];
  vq_all=mdata[][2];
  vs_all=exp(vq_all);
  vchain=mdata[][7]+mdata[][8]+mdata[][9]+mdata[][10]+mdata[][11]+mdata[][12];
  mX_all=ones(iH,1)~mdata[][3:5]~vchain;
  mW_all=ones(iH,1)~mdata[][3:4]~vchain~mdata[][13:14];
  vf_all=mdata[][15]/12;
  vT_all=mdata[][16]*12;
  mX_all=mX_all~((vf_all.*vT_all)/1000);
  vh_m=mdata_market[][1];
  mX_all[][2]=(mX_all[][2]-meanc(mX_all[][2]))
    /sqrt(varc(mX_all[][2]));
  mW_all[][2]=(mW_all[][2]-meanc(mW_all[][2]))
    /sqrt(varc(mW_all[][2]));
  mW_all[][iK_s-2]=(mW_all[][iK_s-2]-meanc(mW_all[][iK_s-2]))
    /sqrt(varc(mW_all[][iK_s-2]));
  mW_all[][iK_s-1]=(mW_all[][iK_s-1]-meanc(mW_all[][iK_s-1]))
    /sqrt(varc(mW_all[][iK_s-1]));
  vlog_one_minus_sum_s_all=zeros(iH,1);
  ih_past=0;
  for(cm=0;cm<iM;cm++)
    {
      iHm=vh_m[cm];      
      vlog_one_minus_sum_s_all[ih_past:ih_past+iHm-1]
	=log(1-sumc(vs_all[ih_past:ih_past+iHm-1]))*ones(iHm,1);
      ih_past+=iHm;
    }

  for(ch=0;ch<iH;ch++)
    {
      if(vT_all[ch]==0 || vf_all[ch]==0)
	{
	  vT_all[ch]=vf_all[ch]=0;
	}
    }

  mTf_all=ones(iH,1);  
  for(ck=1;ck<iK_gamma_T+1;ck++)
    {
      mTf_all=mTf_all~(vT_all.^(ck));
    }
  for(ck=1;ck<iK_gamma_f+1;ck++)
    {
      mTf_all=mTf_all~(vf_all.^(ck));
    }
  mTf_all=mTf_all~(vf_all.*vT_all);
  mTf_all[][1]=mTf_all[][1]/10;
  mTf_all[][2]=mTf_all[][2]/1000;
  mTf_all[][3]=mTf_all[][3]/10;
  mTf_all[][4]=mTf_all[][4]/100;
  mTf_all[][5]=mTf_all[][5]/1000;

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

  //Only Shizuoka
  cm=18;
  s_iHm=vh_m[cm];
  ih_past=sumc(vh_m[0:cm-1]);

  s_mX=mX_all[ih_past:ih_past+s_iHm-1][];
  s_mX[][iK_d-1]=zeros(s_iHm,1);
  s_mW=mW_all[ih_past:ih_past+s_iHm-1][];
  s_vp=vp_all[ih_past:ih_past+s_iHm-1];
  s_vq=vq_all[ih_past:ih_past+s_iHm-1];
  s_vf=vf_all[ih_past:ih_past+s_iHm-1];
  s_vT=vT_all[ih_past:ih_past+s_iHm-1];
  vindex_nozero_f=ones(s_iHm,1);

  for(ch=0;ch<s_iHm;ch++)
    {
      if(s_vf[ch]==0)
	{
	  vindex_nozero_f[ch]=0;
	}
    }
  vstep=range(1, 5000, 5);
  cc =  idiv(5000, 5);

  //Prediction
  mmean_p=zeros(s_iHm,iprediction_end);
  for(cp=0;cp<iprediction_end;cp++)
    {
      file=fopen(sprint(cp,"p_shizuoka.dat"),"r");
      fscan(file,"%#m",irepeat_prediction,s_iHm,&mp_new);
      fclose(file);
      for(ch=0;ch<s_iHm;ch++)
	{
	  DrawXMatrix(ch,thinc(mp_new[][ch]' , cc),ch,vstep,"",0,1);
	}
      SaveDrawWindow(sprint(cp,"p.eps"));
      CloseDrawWindow();
      mmean_p[][cp]=meanc(mp_new)';
    }

  vp_new_mean=meanr(mmean_p);

  itau_5=5*12;
  itau_10=10*12;
  itau_15=20*12;
  itau_20=30*12;
  vtotalpay_5=s_vp*itau_5+s_vf.*s_vT;
  vtotalpay_new_mean_5=vp_new_mean*itau_5;
  vtotalpay_10=s_vp*itau_10+s_vf.*s_vT;
  vtotalpay_new_mean_10=vp_new_mean*itau_10;
  vtotalpay_15=s_vp*itau_15+s_vf.*s_vT;
  vtotalpay_new_mean_15=vp_new_mean*itau_15;
  vtotalpay_20=s_vp*itau_20+s_vf.*s_vT;
  vtotalpay_new_mean_20=vp_new_mean*itau_20;

  iT=sumc(vindex_nozero_f);
  vid=cumulate(ones(s_iHm,1));
  mprediction=vid~(s_vp+s_vf)~vp_new_mean;
  mprediction_5=vid~s_vT~vtotalpay_5~vtotalpay_new_mean_5;
  mprediction_10=vid~s_vT~vtotalpay_10~vtotalpay_new_mean_10;
  mprediction_15=vid~s_vT~vtotalpay_15~vtotalpay_new_mean_15;
  mprediction_20=vid~s_vT~vtotalpay_20~vtotalpay_new_mean_20;
  mprediction_noT=zeros(s_iHm-iT,3);
  mprediction_withT=zeros(iT,3);
  mprediction_5_withT=zeros(iT,4);
  mprediction_10_withT=zeros(iT,4);
  mprediction_15_withT=zeros(iT,4);
  mprediction_20_withT=zeros(iT,4);
  ch_noT=ch_withT=0;
  vT_withT=zeros(iT,1);
  for(chm=0;chm<s_iHm;chm++)
    {
      if(vindex_nozero_f[chm]==1)
	{
	  mprediction_withT[ch_withT][]=mprediction[chm][];
	  mprediction_5_withT[ch_withT][]=mprediction_5[chm][];
	  mprediction_10_withT[ch_withT][]=mprediction_10[chm][];
	  mprediction_15_withT[ch_withT][]=mprediction_15[chm][];
	  mprediction_20_withT[ch_withT][]=mprediction_20[chm][];
	  vT_withT[ch_withT]=s_vT[chm];
	  ch_withT++;
	}
      else
	{
	  mprediction_noT[ch_noT][]=mprediction[chm][];
	  ch_noT++;
	}
    }

  mprediction_withT=sortbyc(mprediction_withT,1);
  mprediction_noT=sortbyc(mprediction_noT,1);
  mprediction_5_withT=sortbyc(mprediction_5_withT,2);
  mprediction_10_withT=sortbyc(mprediction_10_withT,2);
  mprediction_15_withT=sortbyc(mprediction_15_withT,2);
  mprediction_20_withT=sortbyc(mprediction_20_withT,2);

  DrawMatrix(0,mprediction_5_withT[][2]',"Before: 5 years" ,1,1);
  DrawMatrix(0,mprediction_5_withT[][3]',"After: 5 years" ,1,1,2,3);
  DrawMatrix(1,mprediction_10_withT[][2]',"Before: 10 years" ,1,1);
  DrawMatrix(1,mprediction_10_withT[][3]',"After: 10 years" ,1,1,2,3);
  DrawMatrix(2,mprediction_15_withT[][2]',"Before: 20 years" ,1,1);
  DrawMatrix(2,mprediction_15_withT[][3]',"After: 20 years" ,1,1,2,3);
  DrawMatrix(3,mprediction_20_withT[][2]',"Before: 30 years" ,1,1);
  DrawMatrix(3,mprediction_20_withT[][3]',"After: 30 years" ,1,1,2,3);
  SaveDrawWindow("Shizuoka_prediction_T_star_withT.eps");
  CloseDrawWindow();

  DrawMatrix(0,mprediction_withT[][1]',"Before: $p_{h}+f_{h}$" ,1,1);
  DrawMatrix(0,mprediction_withT[][2]',"After: $p_{h}^{new}$" ,1,1,2,3);
  SaveDrawWindow("Shizuoka_prediction_withT.eps");
  CloseDrawWindow();

  DrawMatrix(0,mprediction_noT[][1]',"Before: $p_{h}$" ,1,1);
  DrawMatrix(0,mprediction_noT[][2]',"After: $p_{h}^{new}$" ,1,1,2,3);
  SaveDrawWindow("Shizuoka_prediction_noT.eps");
  CloseDrawWindow();


  format(200000);
  file=fopen("shizuoka_predictors.dat","w");
  fprint(file,"%15.5f",mprediction);
  fclose(file);

  format(200000);
  file=fopen("shizuoka_predictors_Tstar.dat","w");
  fprint(file,"%15.5f",mprediction_5~mprediction_10);
  fclose(file);

  mprediction_all=zeros(iT,11);
  
  mprediction_all[][0]=mprediction_withT[][1];
  mprediction_all[][1]=mprediction_withT[][2];
  mprediction_all[][2]=mprediction_5_withT[][2];
  mprediction_all[][3]=mprediction_5_withT[][3];
  mprediction_all[][4]=mprediction_10_withT[][2];
  mprediction_all[][5]=mprediction_10_withT[][3];
  mprediction_all[][6]=mprediction_15_withT[][2];
  mprediction_all[][7]=mprediction_15_withT[][3];
  mprediction_all[][8]=mprediction_20_withT[][2];
  mprediction_all[][9]=mprediction_20_withT[][3];
  mprediction_all[][10]=vT_withT;

  format(200000);
  file=fopen("0610shizuoka.dat","w");
  fprint(file,"%15.5f",mprediction_all);
  fclose(file);
  
}
