import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#python3 generate13cds.py
import pandas as pd,seaborn as sns
from Bio import SeqIO
from Bio.Phylo.PAML import codeml
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys,random,os,shutil
import statsmodels.stats.api as sms
from scipy.stats import norm

coor=pd.read_csv("./coor_cds_no_stop",sep="\t",header=None,names=['gene','start','end','strand'])#coor_cds
coorCDS=coor.loc[:12,:]
def get_gene_by_pos(df,pos):
 if df[(pos>=df["start"])&(pos<=df["end"])].empty==False:
  gene=df[(pos>=df["start"])&(pos<=df["end"])]["gene"].values[0]
 else:
  gene="***"
 return gene

def return_cds_mut(d):#delete Het pos not in CDS
 d1={}
 for pos in d.keys():
  if coor[(pos>=coor["start"])&(pos<=coor["end"])].empty==False:
   d1[pos]=d[pos]
 return d1

mito=SeqIO.read("/home/liuqi/work/tumor/bwamem/nDNAGeneKaKs/N1_ref_kaks/pro_script/mito_38.fa","fasta")
heter=pd.read_csv("/home/liuqi/work/tumor/bwamem/nDNAGeneKaKs/N1_ref_kaks/pro_script/het_matrix_mut.xls",header=0,index_col=0,sep="\t")
c=heter.loc["N1",:][heter.loc["N1",:].notnull()]
c.index=c.index.map(int)
n1Pos0=c.to_dict()
n1Pos=return_cds_mut(n1Pos0)

heter1=heter.copy()
heter1=heter1.drop("N1",axis=0)
heter1=heter1.loc[:,heter1.notnull().any(axis=0)]
tumorPos0={331:"A/C",2465:"T/A",2778:"T/G",3424:"G/A",3955:"G/A",4937:"T/C",5530:"C/T",6480:"G/C",6565:"A/C",6581:"A/G",6587:"C/A",6593:"A/C",8190:"C/T",8523:"C/A",10844:"A/C",14804:"G/A",16117:"T/C"}#T/C/T 7985:"C/A" is deleted
tumorPos=return_cds_mut(tumorPos0)

print("n1Pos_length:%d tumorPos_length:%d"%(len(n1Pos.keys()),len(tumorPos.keys())))
sub_matrix=pd.read_csv("/home/liuqi/work/tumor/bwamem/nDNAGeneKaKs/N1_ref_kaks/pro_script/sub_matrix.xls",header=0,index_col=0,sep="\t")
d_sub_N10=sub_matrix.loc["N1"].dropna().to_dict()
d_sub_N1=return_cds_mut(d_sub_N10)
print("^^^^^^^^^^^^^^^^^^N1 Substutions Number:%d"%len(d_sub_N1.keys()))#27

def getRef_N1(d):#d_sub_N1 make N1 major Ref
 mutseq=MutableSeq(str(mito.seq))
 n_sub_N1=0
 for i in d.keys():
  i0=int(i)
  ref_alt=re.split(r"\d+",d[i])
  ref=ref_alt[0] 
  alt=ref_alt[1]
  if mutseq[i0-1]==ref:
   mutseq[i0-1]=alt
   n_sub_N1+=1
 if str(mutseq)==str(mito.seq):
  print("!!!!!!!!!Nothing Changed IN N1 Ref")
 print("^^^^^^^^^^^^^^^^^^^^^N1 %d Substutions"%n_sub_N1)
 return n_sub_N1,mutseq

def getMito_N1(d):#change the seq by heter pos
 #print(d)
 mutseq=MutableSeq(str(mito.seq))
 n=0
 for i in d.keys():
  i0=int(i)
  gene_N1=get_gene_by_pos(coorCDS,i0)
  ref_AssU_Ualt_maf=d[i].split(":")
  #print(ref_AssU_Ualt_maf)
  ref_AssU_Ualt=ref_AssU_Ualt_maf[0].split("/")
  maf=ref_AssU_Ualt_maf[1]
  ref=ref_AssU_Ualt[0];AssU=ref_AssU_Ualt[1];Ualt=ref_AssU_Ualt[2]
  #print(i0,ref,alt)
  if mutseq[i0-1]==ref:#16117 not change
   mutseq[i0-1]=Ualt
   n+=1
 if str(mutseq)==str(mito.seq):
  print("!!!!!!!!!You Have Changed Nothing")
 print("##################N1 %d Subs"%n)
 return n,mutseq

def getMito_Tumor(d):#change the seq by heter pos
 #print(d)
 mutseq=MutableSeq(str(mito.seq))
 n=0
 for i in d.keys():
  i0=int(i)
  gene_Tumor=get_gene_by_pos(coorCDS,i0)
  ref_alt=d[i].split("/")
  ref=ref_alt[0];alt=ref_alt[1]
  #print(i0,ref,alt)
  if mutseq[i0-1]==ref:#16117 not change
   mutseq[i0-1]=alt
   n+=1
 print("##################Tumor %d Subs"%n)
 return n,mutseq

#N1,Tumors seq cut the 13 CDS
def getSeq(seq,gene,st,ed,syb):#get the gene seq by start end position from nucleotide_trans_coor
 seq=seq[st-1:ed]
 if syb==1:
  #print(gene)
  seq1=str(seq)
 elif syb==-1:
  #print(gene,"-1")
  seq=SeqRecord(seq)
  seq1=str(seq.reverse_complement().seq)
 return seq1

def getSeq13CDS(seq):
 s13cds="";s13pro=""
 for i in coorCDS.index:
  gene=coorCDS.loc[i,'gene'];st=coorCDS.loc[i,'start'];ed=coorCDS.loc[i,'end']
  syb=coorCDS.loc[i,'strand']
  seq_i=getSeq(seq,gene,st,ed,syb)
  pro_i=Seq(seq_i).translate(table="Vertebrate Mitochondrial")#translate the seq
  s13cds+=seq_i;s13pro+=pro_i
 return s13cds,s13pro

def generate_seq(d_origin,getMito,n_mut,id_syb):#generate mut 13CDS seq by sampling dict n_mut items
 mut_list=random.sample(list(d_origin.keys()),k=n_mut)
 mut_dict={i:d_origin[i] for i in mut_list}
 n_mut_true,mut_seq=getMito(mut_dict)#replace nucl
 mut_seq_cds,mut_seq_pro=getSeq13CDS(mut_seq)
 des=str(n_mut_true)+"_h"
 rec=SeqRecord(Seq(str(mut_seq_cds)),id=id_syb,description=des)
 return rec 

def change_2_I_phylip(fil):#change muscle phylip to paml format in codeml seqfile
 name=fil.split(".")[0]#MT_N1
 opt_change_fil=name+"_change.afa"
 f=open(fil,"r")
 line1=f.readline()
 line1=line1.split("\n")[0]+" I"+"\n"
 line2=f.readline()
 line2=line2[:10]+" "+line2[10:]
 line3=f.readline()
 line3=line3[:10]+" "+line3[10:]
 
 to_file = open(opt_change_fil,mode="w")
 to_file.write(line1)
 to_file.write(line2)
 to_file.write(line3)
 shutil.copyfileobj(f, to_file)
 return opt_change_fil

def change_2_tree(fil):
 name=fil.split(".")[0]
 opt_change_fil=name+"_change.tree"
 f=open(fil,"r")
 line1=f.readline()
 line2=f.readline()
 line2=line2.split("_")[1]
 line3=f.readline()
 line4=f.readline()
 line4=line4.split("_")[1]
 line5=f.readline()
 to_file = open(opt_change_fil,mode="w")
 to_file.write(line1)
 to_file.write(line2)
 to_file.write(line3)
 to_file.write(line4)
 to_file.write(line5)
 to_file.write(";")
 return opt_change_fil

def muscle_pro(records,records_file,optfil):#muscle and change phylip file to PAML file
 SeqIO.write(records, records_file, "fasta")
 muscle_exe=r"/usr/bin/mafft"
 muscle_cline = MafftCommandline(muscle_exe, input=records_file,auto=True,phylipout=True,treeout=True,quiet=True) #out=optfil,tree1=tree_fil,phyi=True)
 stdout, stderr=muscle_cline()
 with open(optfil, "w") as handle:
  handle.write(stdout)
 tree_fil1=records_file+".tree"
 for_codeml_phylip_file=change_2_I_phylip(optfil)#G12_100.afa
 tree_fil=change_2_tree(tree_fil1)
 return for_codeml_phylip_file,tree_fil

def muscle_pro_notwork_in_BIG(records,records_file,optfil,tree_fil):#muscle and change phylip file to PAML file
 SeqIO.write(records, records_file, "fasta")
 muscle_exe=r"/share/nas1/pub/pubtools/muscle/muscle3.8.31_i86linux64"
 muscle_cline = MuscleCommandline(muscle_exe, input=records_file, out=optfil,tree1=tree_fil,phyi=True)
 stdout, stderr=muscle_cline()
 for_codeml_phylip_file=change_2_I_phylip(optfil)
 return for_codeml_phylip_file,tree_fil

def codeml_pro1(n_pro,align_fil,tree_fil,opt_fil):
 print("######CODEML PRO %d"%n_pro)
 name=opt_fil.split("_")[0]
 codeml_path=str(n_pro)+"_"+name+"_codeml_pro"
 os.mkdir(codeml_path)
 cml = codeml.Codeml(alignment = align_fil, tree = tree_fil,out_file = opt_fil, working_dir = codeml_path)
 #cml.read_ctl_file("codeml1.ctl")
 cml.set_options(noisy = 0,verbose = 0,runmode = 0,seqtype = 1,CodonFreq = 2,clock = 0,aaDist = 0,model = 0,NSsites =[0],icode = 1,Mgene = 0,fix_kappa = 0,kappa = 2,fix_omega = 0,omega = .4,fix_alpha = 1,alpha = 0.,Malpha = 0,ncatG = 3,fix_rho = 1,rho = 0.,getSE = 0,RateAncestor = 0)
 #results = cml.run()
 results = cml.run(command="/usr/bin/codeml",verbose = True)
 omega_value=results.get("NSsites").get(0).get("parameters").get("omega")
 print("######CODEML %d PRO OVER!!!!"%n_pro)
 return omega_value

def codeml_pro(n_mut_N1,n_mut_Tumor,n_pro):
 #1:get rec
 #rCRS_cds,rCRS_pro=getSeq13CDS(MutableSeq(str(mito.seq)))#get rCRS
 #rec_mito_cds=SeqRecord(Seq(str(rCRS_cds)),id="MT",description="rCRS")
 
 #rec_N1_Ref=generate_seq(d_sub_N1,getRef_N1,)
 
 rec_N1=generate_seq(n1Pos,getMito_N1,n_mut_N1,"N1")
 rec_Tumor=generate_seq(tumorPos,getMito_Tumor,n_mut_Tumor,"Tumor")
 #2:muscle and replace phylip
 print("####Muscle Pro N1")#get muscle change afa file
 N1_for_codeml_phylip_file,N1_tree=muscle_pro([rec_N1_Ref,rec_N1],"MT_N1_"+str(n_pro)+".recs","MT_N1_"+str(n_pro)+".afa")
 print("####Muscle Pro Tumor")
 Tumor_for_codeml_phylip_file,Tumor_tree=muscle_pro([rec_N1_Ref,rec_Tumor],"MT_Tumor_"+str(n_pro)+".recs","MT_Tumor_"+str(n_pro)+".afa")
 #3:codeml Pro
 N1_omega_value=codeml_pro1(n_pro,N1_for_codeml_phylip_file,N1_tree,"N1_codeml_"+str(n_pro)+".results")
 Tumor_omega_value=codeml_pro1(n_pro,Tumor_for_codeml_phylip_file,Tumor_tree,"Tumor_codeml_"+str(n_pro)+".results")
 return N1_omega_value,Tumor_omega_value

def cal_95CI(l):
 CI=sms.DescrStatsW(l).tconfint_mean()
 CI_1=CI[0]
 CI_2=CI[1]
 return CI_1,CI_2

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

def check_999(l,n_pro,kaks_value):
 if int(kaks_value)==999:
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 999 PRO %d"%n_pro)
 else:
  l.append(kaks_value)

if __name__=="__main__":
 rec_N1_Ref=generate_seq(d_sub_N1,getRef_N1,27,"N1REF")
 print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  0 Pro")
 N1_omega_value,Tumor_omega_value=codeml_pro(12, 12,0)
 N1_list=[]
 Tumor_list=[]
 n=1
 while n<=1000:
  print("****************************************************** %d Pro"%n)
  N1_omega_value_tmp,Tumor_omega_value_tmp=codeml_pro(8, 8,n)
  check_999(N1_list,n,N1_omega_value_tmp)
  check_999(Tumor_list,n,Tumor_omega_value_tmp)
  n+=1

 #N1_list=remove_999(N1_list)
 #Tumor_list=remove_999(Tumor_list)

 N1_CI_1,N1_CI_2=cal_95CI(N1_list)
 Tumor_CI_1,Tumor_CI_2=cal_95CI(Tumor_list)
 print(N1_omega_value,N1_CI_1,N1_CI_2,min(N1_list),max(N1_list))
 print(Tumor_omega_value,Tumor_CI_1,Tumor_CI_2,min(Tumor_list),max(Tumor_list))

 hist_list(N1_list,"N1_CI_hist.pdf")
 hist_list(Tumor_list,"Tumor_CI_hist.pdf")
