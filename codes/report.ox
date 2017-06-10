#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<ReportMCMC.ox>
#include"myconst.h"

static decl s_varname,s_vtrue,file;
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

main()
{
  decl aday,dbm,
    moutput,ci,ipath,mpath,out,
    mcredible_01,mcredible_05,mcredible_10,astar,ck,vmean,msample;
  file=fopen("msample.dat","r");
  fscan(file,"%#m",irepeat,in_param,&moutput);
  fclose(file);
  fn_varname();
  out=new ReportMCMC(moutput);
  out.SetVarNames(s_varname);
  out.Report();
  delete out;
}