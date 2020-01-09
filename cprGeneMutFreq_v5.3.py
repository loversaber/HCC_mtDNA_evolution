import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import numpy as np
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
plt.style.use('ggplot')

dfmitocart=pd.read_excel("Human.MitoCarta2.0.xls",sheet_name="A Human MitoCarta2.0",header=0)
dfsnv269=pd.read_excel("snv269_1.xlsx",header=[0,1],sheet_name="Sheet2")

mutfrq=dfsnv269.xs('MutFreq.',level=1,axis=1)

muttype=pd.DataFrame(dfsnv269['MutType']['Missense/Silent/Nonsense/Splice site '])
dfGeneMutType=dfsnv269.ix[:,dfsnv269.columns.get_level_values(1).isin({"MutFreq.", "Missense/Silent/Nonsense/Splice site "})]#just the 23 tumors not N1 N1 mutfrq is 0

gene12=np.intersect1d(dfGeneMutType.index,dfmitocart['Symbol'])
def indexType(ix):
 if ix in gene12:
  return "M"
 else:
  return "NM"
mutfrq=mutfrq.copy()
mutfrq['Gene']=mutfrq.index
mutfrq.loc[:,'GeneSyb']=mutfrq['Gene'].apply(lambda x: indexType(x))

mutfrq12=mutfrq[mutfrq.index.map(lambda x:x in gene12)]
mutfrq257=mutfrq[mutfrq.index.map(lambda x :x not in gene12)]

fig=plt.figure(figsize=(20,10))

se12,se257=pd.Series(),pd.Series()
#print(mutfrq)
def sizeLabel(ax,title,size,syb,pvalue_syb):
 ax.set_ylim([-0.3, 1.6])
 #if syb=="red":
  #ax.set_title(title,fontsize=size,loc='left',color="red",weight='bold')
 #else:
  #ax.set_title(title,fontsize=size,loc='left',weight='bold')
 ax.set_xlabel("",fontsize=size)#,weight='bold')
 ax.set_ylabel("MutFreq",fontsize=size)
 #ax.text(0.5,1.2,pvalue_syb,horizontalalignment='center', color='black',fontsize=18,weight='semibold')
 ax.tick_params(axis="both",labelsize=size)
 ax.text(0.05, 0.95, title, transform=ax.transAxes,\
      fontsize=16, fontweight='bold', va='top')
dpvalue={}
for i,col in enumerate(mutfrq):
 if col !="Gene" and col != "GeneSyb":
  k,u=kruskal(mutfrq12[col],mutfrq257[col],nan_policy='omit')
  dpvalue[col]={}
  dpvalue[col]['Hpvalue']=u
  se12=se12.append(mutfrq12[col]);se257=se257.append(mutfrq257[col])
  ax1=fig.add_subplot(5,5,i+1)
  sns.violinplot(x="GeneSyb",y=col,data=mutfrq,palette="muted",inner="box",ax=ax1)
  #ax1.set_ylim([-0.3, 1.3])
  title1=col#+" Pvalue:"+"%.2f"%u
  sizeLabel(ax1,title1,18,"bl","Pvalue:"+"%.2f"%u)
  #print(col,u)
ax0=fig.add_subplot(5,5,25)
m1=mutfrq.set_index(['Gene','GeneSyb'])
m2 = m1.stack().reset_index()
sns.violinplot(x="GeneSyb",y=0,data=m2,palette="muted",inner="box",ax=ax0)
#ax0.set_ylim([-0.3, 1.3])
kall,uall=kruskal(se12,se257,nan_policy='omit')

dpvalue["ALL"]={};dpvalue["ALL"]['Hpvalue']=uall
dfpvalue=pd.DataFrame(dpvalue).T
dfpvalue['AdjustedP_BH']=multipletests(dfpvalue['Hpvalue'],method='fdr_bh')[1]
dfpvalue['AdjustedP_BY']=multipletests(dfpvalue['Hpvalue'],method='fdr_by')[1]
dfpvalue['AdjustedP_Bonf']=multipletests(dfpvalue['Hpvalue'],method='bonferroni')[1]
dfpvalue=dfpvalue.sort_values(by='AdjustedP_BH',ascending=True)
#dfpvalue.to_csv("HtestPvalueAdjusted_1221.xls",sep="\t",index=True)

kall1,uall1=kruskal(m2[0][m2['GeneSyb']=="M"],m2[0][m2['GeneSyb']=="NM"],nan_policy='omit')
title0="ALL"# Pvalue:"+"%.2f"%uall#+"%.2f"%uall1
sizeLabel(ax0,title0,18,"bl","Pvalue:"+"%.2f"%uall)
#print("ALL:%g:%g"%(uall,uall1))

fig.tight_layout()
fig.savefig("MitoNoMito_0214.pdf",dpi=200)

#print(dfpvalue)
#print(dfpvalue.index)
#print(dfpvalue.columns)

#se12_2,se257_2=pd.Series(),pd.Series()
fig2=plt.figure(figsize=(20,10))
for i,col in enumerate(mutfrq):
 if col !="Gene" and col != "GeneSyb":
  k2,u2=kruskal(mutfrq12[col],mutfrq257[col],nan_policy='omit')
  #se12_2=se12_2.append(mutfrq12[col]);se257_2=se257_2.append(mutfrq257[col])
  ax2=fig2.add_subplot(5,5,i+1)
  sns.violinplot(x="GeneSyb",y=col,data=mutfrq,palette="muted",inner="box",ax=ax2)
  #ax2.set_ylim([-0.3, 1.3])
  title2=col#+" Pvalue:"+"%.2f"%u2+"/"+"%.2f"%dfpvalue.loc[col,"AdjustedP_BH"]
  sizeLabel(ax2,title2,18,"black","Pvalue:"+"%.2f"%u2+"/"+"%.2f"%dfpvalue.loc[col,"AdjustedP_BH"])
  print(col,u2,dfpvalue.loc[col,"AdjustedP_BH"])
ax02=fig2.add_subplot(5,5,25)

mutfrq1=mutfrq.drop('23TMD',axis=1)
m12=mutfrq1.set_index(['Gene','GeneSyb'])
m22 = m12.stack().reset_index()
sns.violinplot(x="GeneSyb",y=0,data=m22,palette="muted",inner="box",ax=ax02)
#ax02.set_ylim([-0.3, 1.3])
#kall,uall=kruskal(se12,se257,nan_policy='omit')


kall2,uall2=kruskal(m22[0][m22['GeneSyb']=="M"],m22[0][m22['GeneSyb']=="NM"],nan_policy='omit')
####U-test
uall_1,pall_1=mannwhitneyu(m22[0][m22['GeneSyb']=="M"],m22[0][m22['GeneSyb']=="NM"])
print("##########U-test ALL:",pall_1)

title02="ALL"# Pvalue:"+"%.2f"%uall2+"/%.2f"%dfpvalue.loc["ALL","AdjustedP_BH"]
sizeLabel(ax02,title02,18,"black","Pvalue:"+"%.2f"%uall2+"/%.2f"%dfpvalue.loc["ALL","AdjustedP_BH"])
ax02.annotate("",xy=(0,1.15),xycoords="data",xytext=(1,1.15),\
textcoords="data",arrowprops=dict(arrowstyle="-", ec='k',connectionstyle="bar,fraction=0.1",lw=3),fontsize=16,weight='semibold') 

ax02.text(0.5,1.4,"*",horizontalalignment='center',verticalalignment="center",color='black',fontsize=22,weight='semibold')
print("ALL:%g:%g"%(uall2,dfpvalue.loc["ALL","AdjustedP_BH"]))

fig2.tight_layout()
fig2.savefig("MitoNoMitoAdjustPvalue_0412.pdf",dpi=200)

fig3,ax_all=plt.subplots(figsize=(6,5))
newPal   = dict(NM= "#4878CF",M= "#6ACC65")
newPal2=dict(NM="#000080",M="#008000")
ec,lw,sz,alp='w',0.75,2,0.75
sns.swarmplot(x="GeneSyb",y=0,data=m22,palette=newPal2,edgecolor=ec,linewidth=lw,size=sz,alpha=alp,ax=ax_all)
#sns.stripplot(x="GeneSyb",y=0,data=m22,palette=newPal2,edgecolor=ec,linewidth=lw,size=sz,alpha=alp,ax=ax_all)
sns.violinplot(x="GeneSyb",y=0,data=m22,palette=newPal,inner=None,ax=ax_all)
sizeLabel(ax_all,"",18,"black","Pvalue:"+"%.2f"%uall2+"/%.2f"%dfpvalue.loc["ALL","AdjustedP_BH"])
ax_all.annotate("",xy=(0,1.15),xycoords="data",xytext=(1,1.15),\
textcoords="data",arrowprops=dict(arrowstyle="-", ec='k',connectionstyle="bar,fraction=0.1",lw=3),fontsize=16,weight='semibold')

ax_all.text(0.5,1.4,"*",horizontalalignment='center',verticalalignment="center",color='black',fontsize=22,weight='semibold')
ax_all.set_ylabel("Mutation Frequency",fontsize=18,fontweight="semibold")
fig3.tight_layout()
fig3.savefig("ALL_0916.pdf",dpi=200,bbox_inches='tight')
