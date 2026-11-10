exec(open('/tmp/sc_all.py').read().split('# GSE275018')[0])
import scipy.io as sio, glob, os, pandas as pd, numpy as np, scanpy as sc
# GSE274828 — mechanical ventilation, 4 samples (YNV/YV/ANV/AV)
D="sc/GSE274828"
pref=sorted({os.path.basename(f).split("_matrix")[0] for f in glob.glob(D+"/*_matrix.mtx.gz")})
ads=[]
for p in pref:
    tag=p.split("_",1)[1]
    X=sio.mmread(f"{D}/{p}_matrix.mtx.gz").T.tocsr()
    feat=pd.read_csv(f"{D}/{p}_features.tsv.gz",sep="\t",header=None)
    bc=pd.read_csv(f"{D}/{p}_barcodes.tsv.gz",sep="\t",header=None)[0].values
    ad=sc.AnnData(X); ad.var_names=feat[1].astype(str).values if feat.shape[1]>1 else feat[0].astype(str).values
    ad.var_names_make_unique(); ad.obs_names=[f"{tag}_{b}" for b in bc]; ad.obs["sample"]=tag
    ads.append(ad); print(" ",tag,ad.shape,flush=True)
pipeline(sc.concat(ads,join="inner"),"GSE274828")
# GSE186185 — orthodontic force, 2 CSV matrices
D="sc/GSE186185"; ads=[]
for f in sorted(glob.glob(D+"/*Expression.csv.gz")):
    tag=os.path.basename(f).split("_")[1]
    d=pd.read_csv(f,index_col=0)
    ad=sc.AnnData(d.T.values.astype("float32")); ad.var_names=d.index.astype(str); ad.var_names_make_unique()
    ad.obs_names=[f"{tag}_{c}" for c in d.columns]; ad.obs["sample"]=tag
    ads.append(ad); print(" ",tag,ad.shape,flush=True)
pipeline(sc.concat(ads,join="inner"),"GSE186185")
