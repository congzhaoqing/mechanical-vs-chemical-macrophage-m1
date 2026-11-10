
import os, glob, gzip, numpy as np, pandas as pd, scanpy as sc, scipy.io as sio
import scipy.sparse as sp
sc.settings.verbosity=1
if not hasattr(sp.csr_array,"getnnz"): sp.csr_array.getnnz = sp.csr_array._getnnz
if not hasattr(sp.csc_array,"getnnz"): sp.csc_array.getnnz = sp.csc_array._getnnz
D="data/GSE275018"
pref = sorted({os.path.basename(f).split("-matrix")[0] for f in glob.glob(D+"/*-matrix.mtx.gz")})
arm = {"A":"PBS","B":"MNM","C":"RMF"}
ads=[]
for p in pref:
    tag = p.split("_",1)[1]          # A_1
    a = arm[tag.split("_")[0]]
    X = sio.mmread(f"{D}/{p}-matrix.mtx.gz").T.tocsr()
    feat = pd.read_csv(f"{D}/{p}-features.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(f"{D}/{p}-barcodes.tsv.gz", sep="\t", header=None)[0].values
    ad = sc.AnnData(X)
    ad.var_names = feat[1].astype(str).values if feat.shape[1]>1 else feat[0].astype(str).values
    ad.var_names_make_unique()
    ad.obs_names = [f"{tag}_{b}" for b in bc]
    ad.obs["sample"]=tag; ad.obs["arm"]=a; ad.obs["gsm"]=p.split("_")[0]
    ads.append(ad)
    print(tag, a, ad.shape, flush=True)
ad = sc.concat(ads, join="inner")
ad.layers["counts"]=ad.X.copy()
print("concat", ad.shape, flush=True)

ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
n0=ad.n_obs
ad = ad[(ad.obs.n_genes_by_counts>=200)&(ad.obs.pct_counts_mt<20)].copy()
sc.pp.filter_genes(ad, min_cells=3)
print(f"QC min_genes>=200 pct_mt<20: {n0} -> {ad.n_obs}", flush=True)

ad.raw=None
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
ad.layers["logn"]=ad.X.copy()
sc.pp.highly_variable_genes(ad, n_top_genes=2000)
sc.pp.pca(ad, n_comps=50, mask_var="highly_variable")
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
sc.tl.leiden(ad, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
print("clusters", ad.obs.leiden.nunique(), flush=True)

M = {
 "Myeloid":["Lyz2","Csf1r","Itgam","Adgre1","Cd68","Fcgr1","Ctss","Tyrobp","Fcer1g"],
 "Neutrophil":["S100a8","S100a9","Retnlg","Ly6g"],
 "T_NK":["Cd3e","Cd3d","Cd3g","Trbc2","Nkg7","Gzmb","Klrb1c"],
 "B":["Cd79a","Cd79b","Ms4a1","Ighm"],
 "Tumor_Epi":["Epcam","Krt8","Krt18","Mlana","Pmel","Dct","Tyr"],
 "Fibroblast":["Col1a1","Col1a2","Dcn","Pdgfrb"],
 "Endothelial":["Pecam1","Cdh5","Flt1"],
}
for k,g in M.items():
    gg=[x for x in g if x in ad.var_names]
    sc.tl.score_genes(ad, gg, score_name="sc_"+k)
    print("markers",k,gg, flush=True)
sm = ad.obs.groupby("leiden", observed=True)[["sc_"+k for k in M]].mean()
sm.columns=[c[3:] for c in sm.columns]
lab = sm.idxmax(axis=1)
ad.obs["celltype"]=ad.obs.leiden.map(lab).astype(str)
sm["label"]=lab; sm["n_cells"]=ad.obs.leiden.value_counts()
sm.to_csv("out/GSE275018_cluster_markers.csv")
print(sm.round(3).to_string(), flush=True)

for subset,name in [(["Myeloid","Neutrophil"],"myeloid"), (["Myeloid"],"mnp")]:
    sel = ad.obs.celltype.isin(subset)
    sub = ad[sel]
    cnt = sub.obs.groupby(["sample","arm"], observed=True).size()
    print(name, "cells", int(sel.sum()), "\n", cnt.to_string(), flush=True)
    pb={}
    for s in sorted(sub.obs["sample"].unique()):
        m = sub[sub.obs["sample"]==s].layers["counts"]
        pb[s]=np.asarray(m.sum(axis=0)).ravel()
    pbdf=pd.DataFrame(pb, index=sub.var_names)
    pbdf.to_csv(f"out/GSE275018_pseudobulk_{name}.csv")
    cnt.to_frame("n_cells").to_csv(f"out/GSE275018_ncells_{name}.csv")
print("SCDONE", flush=True)
