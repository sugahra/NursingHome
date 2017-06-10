#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include <oxdraw.h>
#import<maximize>
#import<solvenle>
#include"myconst.h"


static decl s_vbeta_d,s_dalpha, s_vbeta_s,s_vgamma,
  s_dtau_omega,s_dtau_xi,s_dsigma2_omega,s_dsigma2_xi,
  s_vj,s_mijj2,
  s_vconst,s_dconst,s_vone_minus_sum_s,
  s_vn_xi,s_vn_omega,s_iHm,s_vXB,s_vWB,s_dWB,s_dXB,s_vq_others,
  s_dp,s_dq,s_ds,s_dsigma_omega,s_dsigma_xi,
  s_vq,s_vp,s_mX,s_mW,
  s_vGamma_T,s_vs_all,s_chm,s_dlog_one_minus_sum_s;

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
    cp,ip,vmean_p;

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

  file=fopen("msample.dat");
  fscan(file,"%#m",irepeat,in_param,&msample);
  fclose(file);
  msample=msample';

  //Only Shizuoka
  cm=18;
  s_iHm=vh_m[cm];
  ih_past=sumc(vh_m[0:cm-1]);
  s_mX=mX_all[ih_past:ih_past+s_iHm-1][];
  s_mX[][iK_d-1]=zeros(s_iHm,1);
  s_mW=mW_all[ih_past:ih_past+s_iHm-1][];

  //Prediction
  vmean_p=zeros(iprediction_end,1);
  for(cp=0;cp<iprediction_end;cp++)
    {
      ip=cp*iprediction_interval;
      s_vbeta_d=msample[0:iK_d-1][ip];
      s_vgamma=msample[iK_d:iK_d+iK_gamma-1][ip];
      s_dalpha=msample[iK_d+iK_gamma][ip];
      s_vbeta_s=msample[iK_d+iK_gamma+1:iK_d+iK_gamma+1+iK_s-1][ip];
      s_dtau_omega=msample[iK_d+iK_gamma+iK_s+1][ip];
      s_dtau_xi=msample[iK_d+iK_gamma+iK_s+2][ip];
      s_dsigma2_omega=1/s_dtau_omega;
      s_dsigma2_xi=1/s_dtau_xi;
      s_dsigma_omega=sqrt(s_dsigma2_omega);
      s_dsigma2_xi=sqrt(s_dsigma2_xi);

      vTgamma=mTf_all*s_vgamma;
      vexpTgamma=exp(vTgamma);
      vGamma_T=vexpTgamma./(1+vexpTgamma);
      for(ch=0;ch<iH;ch++)
	{
	  if(mTf_all[ch][1]==0)
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

      vxi_all=vq_all -mX_all*s_vbeta_d + s_dalpha*vp_all-vlog_one_minus_sum_s_all;
      s_vn_xi=fn_polya_n(s_dtau_xi,vxi_all);
      vomega_all=log(vp_all+ vf_all.*vGamma_T+(ones(iH,1)./(s_dalpha*(vs_all-1)))) 
	-mW_all*s_vbeta_s;
      s_vn_omega=fn_polya_n(s_dtau_omega,vomega_all);

      s_vp=vp_all[ih_past:ih_past+s_iHm-1]
	+vf_all[ih_past:ih_past+s_iHm-1].*vGamma_T[ih_past:ih_past+s_iHm-1];
      s_vq=vq_all[ih_past:ih_past+s_iHm-1];

      mq_sample=mp_sample=zeros(s_iHm,irepeat_prediction);
      vsigma_p_rw=0.5*ones(s_iHm,1);
      vsigma_q_rw=2*ones(s_iHm,1);
      vmh_p=vmh_q=zeros(s_iHm,1);
      s_vXB=s_mX*s_vbeta_d;
      s_vWB=s_mW*s_vbeta_s;

      //Iteration start
      for(cloop=-iburn_prediction;cloop<irepeat_prediction;cloop++)
	{
	  //println(cloop);      
	  //Sample q	  
	  for(chm=0;chm<s_iHm;chm++)
	    {
	      //println(chm);      
	      dold=s_vq[chm];
	      s_dp=s_vp[chm];
	      s_dWB=s_vWB[chm];
	      s_chm=chm;	  
	      dU_q_one=log(1-sumc(exp(s_vq))+exp(dold));	  
	      dU_q_two=log(1-1/(s_dalpha*s_dp));
	      dU_q=min(dU_q_one,dU_q_two);
	      fn_q(dold,&dtarget_old,0,0);
	      dnew=dold+vsigma_q_rw[chm]*fn_randtrann((dU_q-dold)/vsigma_q_rw[chm],1);
	      fn_q(dnew,&dtarget_new,0,0);
	      dproposal_old=-log(probn((dU_q-dnew)/vsigma_q_rw[chm]));
	      dproposal_new=-log(probn((dU_q-dold)/vsigma_q_rw[chm]));
	      dmh=exp(dtarget_new+dproposal_old-dtarget_old-dproposal_new);
	      if(ranu(1,1)<=dmh)
		{     
		  s_vq[chm]=dnew;
		  if(cloop>=0){vmh_q[chm]++;}
		}
	    }

	  //Sample p
	  s_dlog_one_minus_sum_s=log(1-sumc(exp(s_vq)));
	  for(chm=0;chm<s_iHm;chm++)
	    {
	      //println(chm);      
	      dold=log(s_vp[chm]);
	      s_dWB=s_vWB[chm];
	      s_dXB=s_vXB[chm];
	      s_dq=s_vq[chm];
	      s_ds=exp(s_dq);
	      dL_p=-log(s_dalpha*(1-s_ds));
	      fn_p(exp(dold),&dtarget_old,0,0);
	      dnew=dold+vsigma_p_rw[chm]*fn_randtrann((dL_p-dold)/vsigma_p_rw[chm],3);
	      fn_p(exp(dnew),&dtarget_new,0,0);
	      dproposal_old=-dold-log(probn((dnew-dL_p)/vsigma_p_rw[chm]));
	      dproposal_new=-dnew-log(probn((dold-dL_p)/vsigma_p_rw[chm]));
	      dmh=exp(dtarget_new+dproposal_old-dtarget_old-dproposal_new);

	      if(ranu(1,1)<=dmh)
		{     
		  s_vp[chm]=exp(dnew);
		  if(cloop>=0){vmh_p[chm]++;}
		}
	    }

	  if(cloop>=0)
	    {
	      mp_sample[][cloop]=s_vp;
	      mq_sample[][cloop]=s_vq;
	    }
	}
      format(200000);
      file=fopen(sprint(cp,"p_shizuoka.dat"),"w");
      fprint(file,"%15.5f",mp_sample');
      fclose(file);

      format(200000);
      file=fopen(sprint(cp,"q_shizuoka.dat"),"w");
      fprint(file,"%15.5f",mq_sample');
      fclose(file);
    }

    

}
