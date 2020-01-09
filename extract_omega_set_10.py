import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from Bio.Phylo.PAML import codeml
from multiprocessing import Pool
import statsmodels.stats.api as sms
from scipy.stats import norm

#print(Tumor_result_fils)
def get_KaKs(fil):
 results=codeml.read(fil)
 omega_value=results.get("NSsites").get(0).get("parameters").get("omega")
 return omega_value

def get_list(fils):
 p=Pool(50)
 multiRecs=[p.apply_async(get_KaKs,args=(fil,)) for fil in fils]
 kaks_list=[rec.get() for rec in multiRecs if rec.get()<10]
 kaks_list=list(set(kaks_list))
 print(len(kaks_list))
 return kaks_list

def cal_95CI(l):
 min_l,max_l=min(l),max(l)
 CI=sms.DescrStatsW(l).tconfint_mean()
 CI_1=CI[0]
 CI_2=CI[1]
 return CI_1,CI_2,min_l,max_l

def hist_list(list_l,optpdf):
 plt.style.use("ggplot")
 fig=plt.figure()
 #n_hist,bins_h,patch=plt.hist(list_l,bins=20,density=False,rwidth=100,alpha=1,edgecolor='white')
 sns.distplot(list_l, bins=20, fit=norm, norm_hist=False,color='b', kde=False,fit_kws={'label': 'Fitted Norm', "color": "r"}, hist_kws={"edgecolor": "white", "alpha": 1})
 plt.xlabel('Permutation Ka/Ks',fontsize=15,fontweight='semibold')
 #plt.ylabel('Frequency',fontsize=15,fontweight='semibold')
 plt.legend(loc='best')
 plt.tight_layout()
 fig.savefig(optpdf,dpi=200)

def pro_fils(fils,optfil):
 c=get_list(fils)
 hist_list(c,optfil)
 CI_1,CI_2,min_l,max_l=cal_95CI(c)
 print(CI_1,CI_2,min_l,max_l)

if __name__=="__main__":
 Tumor_result_fils=glob.glob("Tumor_codeml_*.results")
 N1_result_fils=glob.glob("N1_codeml_*.results")
 pro_fils(Tumor_result_fils,"T_hist.pdf")
 pro_fils(N1_result_fils,"N_hist.pdf")
