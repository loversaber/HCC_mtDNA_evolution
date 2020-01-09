#python3 cprN1AndTumor_v2.py caddxlsfils.txt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal,mannwhitneyu,ranksums,ttest_ind

testMethods=[kruskal,mannwhitneyu,ranksums,ttest_ind]
def test_Tumor_2_type_mut(testm,df,mut1,mut2):
 statistic1,pvalue1=testm(df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']==mut1],\
			  df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']==mut2])
 return pvalue1

def tests(testm,df,syb,NT):
 statistic1,pvalue1=testm(df['heteroplasmy'][df['Normal&Tumor']==NT][df['Code&NonCode']=='NC'],\
			  df['heteroplasmy'][df['Normal&Tumor']==NT][df['Code&NonCode']=='C'])
 if syb=="2":
  statistic2,pvalue2=testm(df['heteroplasmy'][df['Normal&Tumor']==NT][df['NSS&NC']=='NS'],\
			   df['heteroplasmy'][df['Normal&Tumor']==NT][df['NSS&NC']=='S'],\
			   df['heteroplasmy'][df['Normal&Tumor']==NT][df['NSS&NC']=='NC'])
  return pvalue1,pvalue2
 else:
  return pvalue1

def posnum(df):
 n2c_maf=df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['Code&NonCode']=="C"]
 n2nc_maf=df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['Code&NonCode']=="NC"]
 t2c_maf=df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['Code&NonCode']=="C"]
 t2nc_maf=df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['Code&NonCode']=="NC"]
 n2c_mean=n2c_maf.mean();n2nc_mean=n2nc_maf.mean();t2c_mean=t2c_maf.mean();t2nc_mean=t2nc_maf.mean()
 
 n3s_maf=df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['NSS&NC']=="S"]
 n3ns_maf=df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['NSS&NC']=="NS"]
 n3nc_maf=df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['NSS&NC']=="NC"]
 t3s_maf=df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']=="S"]
 t3ns_maf=df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']=="NS"]
 t3nc_maf=df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']=="NC"]
 n3s_mean=n3s_maf.mean();n3ns_mean=n3ns_maf.mean();n3nc_mean=n3nc_maf.mean()
 t3s_mean=t3s_maf.mean();t3ns_mean=t3ns_maf.mean();t3nc_mean=t3nc_maf.mean() 

 normal2_C=n2c_maf.count();normal2_NC=n2nc_maf.count()
 tumor2_C=t2c_maf.count()
 t2c=tumor2_C/23 
 tumor2_NC=t2nc_maf.count()
 t2nc=tumor2_NC/23

 normal3_S=n3s_maf.count();normal3_NS=n3ns_maf.count();normal3_NC=n3nc_maf.count()
 tumor3_S=t3s_maf.count()
 t3s=tumor3_S/23
 tumor3_NS=t3ns_maf.count()
 t3ns=tumor3_NS/23
 tumor3_NC=t3nc_maf.count()
 t3nc=tumor3_NC/23

 n2c_median=n2c_maf.median();n2nc_median=n2nc_maf.median();t2c_median=t2c_maf.median();t2nc_median=t2nc_maf.median()
 n3s_median=n3s_maf.median();n3ns_median=n3ns_maf.median();n3nc_median=n3nc_maf.median()
 t3s_median=t3s_maf.median();t3ns_median=t3ns_maf.median();t3nc_median=t3nc_maf.median()

 return [normal2_C,tumor2_C,normal2_NC,tumor2_NC,normal3_S,tumor3_S,normal3_NS,tumor3_NS,normal3_NC,tumor3_NC,t2c,t2nc,t3s,t3ns,t3nc],\
        [n2c_mean,t2c_mean,n2nc_mean,t2nc_mean,n3s_mean,t3s_mean,n3ns_mean,t3ns_mean,n3nc_mean,t3nc_mean],\
        [n2c_median,t2c_median,n2nc_median,t2nc_median,n3s_median,t3s_median,n3ns_median,t3ns_median,n3nc_median,t3nc_median]

def testsNT(testm,df,syb,prosyb):
 if syb=="2":#CNC
  statistic1,pvalue1=testm(df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['Code&NonCode']==prosyb],\
                          df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['Code&NonCode']==prosyb])
  return pvalue1
 elif syb=="3":#NSSNC
  statistic2,pvalue2=testm(df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['NSS&NC']==prosyb],\
                           df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']==prosyb])
  return pvalue2

df1,df2=pd.DataFrame(),pd.DataFrame()
with open(sys.argv[1],"r")as fils:
 for fil in fils:
  fil=fil.strip();name=fil.split("_")[0]
  fildf=pd.read_table(fil,sep="\t",header=0)
  fildf.loc[:,"filname"]=name
  df1=df1.append(fildf, ignore_index=True)#23 tumor  ignore_index=True should be or index will be 01 012 0123
 def indexNormalTumor(idx):
  if idx=="N1":
   return "N1"
  else:
   return "Tumor" 
 def codeNoncode(Code_NonCode):
  if Code_NonCode =="Code":
   return "C"
  else:
   return "NC"
 df1=df1.copy()
 df1.loc[:,'Normal&Tumor']=df1["filname"].map(indexNormalTumor)
 df1['Code&NonCode']=df1['Code&NonCode'].map(codeNoncode)

 lposnum,lmean,lmedian=posnum(df1)#get N1 pos num tumor avg num and mean of maf
 print(lposnum)
 print(lmean)
 lmean2f=["%.4f"%i for i in lmean]
 print(lmean2f)
 print(lmedian)
 lmedian2f=["%.4f"%i for i in lmedian]
 print(lmedian2f)
 print(df1.head(15))#work well
 print(df1["NSS&NC"].unique())
 #########cpr NS NC S in Tumor
 dpvalue_Tumor={}
 #test_Tumor_2_type_mut(testm,df,mut1,mut2) 
 def test_Tumor_3(df,mut12name):
  mut12=mut12name.split("_")
  mut1=mut12[0]
  mut2=mut12[1]
  dpvalue_Tumor[mut12name]={}
  for test in testMethods:
   dpvalue_Tumor[mut12name][test.__name__]=test_Tumor_2_type_mut(test,df,mut1,mut2)
 ########use Mann test to check tumor ns is different with S NC in tumor 
 test_Tumor_3(df1,"S_NS")
 test_Tumor_3(df1,"S_NC")
 test_Tumor_3(df1,"NS_NC")
 ########check delete how much outliner Ns in tumor will be not significant in tumor NS S NC
 def mann_check(df,mut1,mut2="NS"): 
  c1=list(df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']==mut1])
  c1.sort()
  c2=list(df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']==mut2])
  c2.sort()
  print(len(c1),len(c2))
  pvalue_99=0.01
  while pvalue_99<0.05:
   statistic_99,pvalue_99=mannwhitneyu(c1,c2)
   print("&"*40,pvalue_99)
   c2.remove(max(c2))
  print("After Delete:",len(c2),pvalue_99)
 mann_check(df1,"NC","NS")  
 mann_check(df1,"S","NS") 
 #########check delete how much outliner Ns in tumor will be not significant in tumor NS and N1 NS
 def mann_check_NS_in_N1_Tumor(df): 
  c1=list(df['heteroplasmy'][df['Normal&Tumor']=="N1"][df['NSS&NC']=="NS"])
  c1.sort()
  #c1=c1[:-2]
  #c1.remove(max(c1))
  #c1.remove(max(c1))
  c2=list(df['heteroplasmy'][df['Normal&Tumor']=="Tumor"][df['NSS&NC']=="NS"])
  c2.sort()
  print(len(c1),len(c2))
  pvalue_99=0.01
  while pvalue_99<0.05:
   statistic_99,pvalue_99=mannwhitneyu(c1,c2)
   print("&"*40,pvalue_99)
   c2.remove(max(c2))
  print("After Delete:",len(c2),pvalue_99)
 mann_check_NS_in_N1_Tumor(df1) 

 dpvalue={}
 def testdf(df,dfname,NT): 
  dpvalue[dfname]={}
  for test in testMethods:
   if test.__name__=="kruskal":
    p1,p2=tests(test,df,"2",NT)#NSSNC
    dpvalue[dfname]['kruskal1']=p1;dpvalue[dfname]['kruskal2']=p2
   else:
    p1=tests(test,df,"1",NT)#only use ttest eg. to test NC C
    dpvalue[dfname][test.__name__]=p1
 testdf(df1,"N1_MAF",'N1')#N1 pvalues
 testdf(df1,"Tumor_MAF",'Tumor')#tumors pvalue
 
 dpvalue2={}
 def testdfNT(df,prosname,syb):
  dpvalue2[prosname]={}
  if syb=="2":
   for pro in["C","NC"]: #["Code","NonCode"]:
    dpvalue2[prosname][pro]={}
    for testm in testMethods:
     pvalue=testsNT(testm,df,"2",pro)
     dpvalue2[prosname][pro][testm.__name__]=pvalue
  elif syb=="3":
   for pro in ["S","NS","NC"]:
    dpvalue2[prosname][pro]={}
    for testm in testMethods:
     pvalue=testsNT(testm,df,"3",pro)
     dpvalue2[prosname][pro][testm.__name__]=pvalue
 testdfNT(df1,"CNC","2");testdfNT(df1,"NSSNC","3") 
 #print(dpvalue2)#store the cpr of N1 and Tumor each mut MAF work well
 #print(dpvalue)
 def dict2df(dpvalue,s):
  dfpvalue=pd.DataFrame(dpvalue).T
  dfpvalue.to_csv(s+".xls")
  return(dfpvalue)

 df_dpvalue_Tumor=dict2df(dpvalue_Tumor,"Tumor_mut3_0528")
 print("*"*30)
 print("^"*30,"MANN Whitney U test")
 print(df_dpvalue_Tumor)
 print("*"*30)
 dfpvalue1=dict2df(dpvalue,"pvalue_N1_Tumor_0528")
 print(dpvalue)
 print(dfpvalue1)
 print("##########################")
 dfpvalue2CNC=dict2df(dpvalue2['CNC'],"pvalue_N1_Tumor_cpr_CNC_0528")
 dfpvalue2NSSNC=dict2df(dpvalue2['NSSNC'],"pvalue_N1_Tumor_cpr_NSSNC_0528")
 print(dfpvalue2CNC)
 print("##########################")
 print(dfpvalue2NSSNC)
 
 def plotDf(df,dfpvalue,picname,ec,lw,sz,alp):
  plt.style.use("ggplot")
  fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))#;fig2=plt.figure(figsize=(20,10))
  #ax1=axs[0,0];ax2=axs[0,1]
  newPal   = dict(N1= "#4878CF",Tumor= "#6ACC65")
  newPal2=dict(N1="#000080",Tumor="#008000")
  sns.stripplot(x="Code&NonCode",y='heteroplasmy',hue="Normal&Tumor",hue_order=['N1','Tumor'],palette=newPal2,data=df,dodge=True,order=['C','NC'],edgecolor=ec,linewidth=lw,size=sz,alpha=alp,ax=ax1)
  #ax1=ax1.legend_.remove()
  sns.violinplot(x="Code&NonCode",y='heteroplasmy',hue="Normal&Tumor",hue_order=['N1','Tumor'],data=df,palette=newPal,inner=None,order=['C','NC'],ax=ax1)#,**{showmedians:True})
  handles, labels = ax1.get_legend_handles_labels()
  ax1.legend(handles[0:2], labels[0:2],loc=1,framealpha=0,bbox_to_anchor=(1, 0.96))
  #l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2
  sns.stripplot(x="NSS&NC",y='heteroplasmy',hue="Normal&Tumor",hue_order=['N1','Tumor'],palette=newPal2,data=df,dodge=True,order=['S','NS','NC'],edgecolor=ec,linewidth=lw,size=sz,alpha=alp,ax=ax2)
  sns.violinplot(x="NSS&NC",y='heteroplasmy',hue="Normal&Tumor",hue_order=['N1','Tumor'],data=df,palette=newPal,inner=None,order=['S','NS','NC'],ax=ax2)
  handles2, labels2 = ax1.get_legend_handles_labels()
  ax2.legend(handles2[0:2], labels2[0:2],loc=1,framealpha=0,bbox_to_anchor=(1, 0.96))
#[normal2_C,tumor2_C,normal2_NC,tumor2_NC,normal3_S,tumor3_S,normal3_NS,tumor3_NS,normal3_NC,tumor3_NC,10:t2c,t2nc,t3s,t3ns,t3nc]

  title1="Pvalue:N1-%.2f Tumor-%.2f C-%.2f NC-%.2f"%(dfpvalue.loc["N1_MAF",'kruskal1'],dfpvalue.loc["Tumor_MAF",'kruskal1'],\
				  dfpvalue2CNC.loc["C","kruskal"],dfpvalue2CNC.loc["NC","kruskal"])
  title2="Pvalue:N1-%.2f Tumor-%.2f S-%.2f NS-%.2f NC-%.2f"%(dfpvalue.loc["N1_MAF",'kruskal2'],dfpvalue.loc["Tumor_MAF",'kruskal2'],\
					dfpvalue2NSSNC.loc["S","kruskal"],\
					dfpvalue2NSSNC.loc["NS","kruskal"],\
					dfpvalue2NSSNC.loc["NC","kruskal"])
  
  def axlabel(ax,xlabel,title,size,xcoor,x_l):
   #ax.legend(loc=1,framealpha=0,bbox_to_anchor=(x_l, 0.96))#up right 1;7 is right center #make the legend in the position#framealpha=0
   ax.set_ylim([-0.1, 0.46])
   #ax.text(xcoor,0.4,title,horizontalalignment='left', color='black',fontsize=16,weight='semibold')
   #ax.set_title("",fontsize=15,fontweight='bold',loc='left')
   ax.set_xlabel(xlabel,fontsize=size,fontweight='semibold')
   ax.set_ylabel("MAF",fontsize=size,fontweight='semibold') 
   ax.tick_params(axis="both",labelsize=size)
   plt.setp(ax.get_legend().get_texts(), fontsize=size)  # for legend text
   plt.setp(ax.get_legend().get_title(), fontsize=size)
  axlabel(ax1,"",title1,18,-0.4,1);axlabel(ax2,"",title2,18,-0.4,1)
  #########draw significant syb
  #y1=0.4;y2=0.42
  #text1=0.42;text2=0.45
  y1=0.39;y2=0.41;y3=0.435
  text1=0.41;text2=0.43;text3=0.455
  ax2.annotate("",xy=(0.8,y1),xycoords="data",xytext=(1.2,y1),\
textcoords="data",arrowprops=dict(arrowstyle="-", ec='k',connectionstyle="bar,fraction=0.2",lw=3),fontsize=16,weight='semibold')

  ax2.annotate("",xy=(0.2,y3),xycoords="data",xytext=(1.2,y3),\
textcoords="data",arrowprops=dict(arrowstyle="-", ec='k',connectionstyle="bar,fraction=0.1",lw=3),fontsize=16,weight='semibold')
  ax2.annotate("",xy=(1.2,y2),xycoords="data",xytext=(2.2,y2),\
textcoords="data",arrowprops=dict(arrowstyle="-", ec='k',connectionstyle="bar,fraction=0.1",lw=3),fontsize=16,weight='semibold')

  ax2.text(1,text1,"*",horizontalalignment='center',verticalalignment="center",color='black',fontsize=22,weight='semibold')
  ax2.text(0.7,text3,"*",horizontalalignment='center',verticalalignment="center",color='black',fontsize=22,weight='semibold')
  ax2.text(1.7,text2,"*",horizontalalignment='center',verticalalignment="center",color='black',fontsize=22,weight='semibold')
  #plt.title("The Difference of MAF between Different Types of heteroplasmic mutations")
  fig.tight_layout()
  fig.savefig(picname+".pdf",dpi=200,bbox_inches='tight')
 
 ec,lw,sz,alp='w',0.75,6,0.75
 plotDf(df1,dfpvalue1,"N1_Tumor_0916",ec,lw,sz,alp)
