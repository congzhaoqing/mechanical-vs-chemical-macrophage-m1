exec(open('/tmp/sc_all.py').read().split('# GSE275018')[0])
import glob,os,pandas as pd,numpy as np,scanpy as sc
ads=[]
for f in sorted(glob.glob("sc/GSE186185/*Expression.csv.gz")):
    tag=os.path.basename(f).split("_")[1]
    d=pd.read_csv(f)
    sym=d["GeneSymbol"].astype(str).values
    m=d.drop(columns=["GeneSymbol","GeneEnsembl"]).values.astype("float32")
    ad=sc.AnnData(m.T); ad.var_names=sym; ad.var_names_make_unique()
    ad.obs_names=[f"{tag}_{c}" for c in d.columns[2:]]; ad.obs["sample"]=tag
    ads.append(ad); print(" ",tag,ad.shape,flush=True)
pipeline(sc.concat(ads,join="inner"),"GSE186185")
