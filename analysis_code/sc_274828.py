import os, glob, numpy as np, pandas as pd, scanpy as sc, scipy.io as sio
import scipy.sparse as sp
sc.settings.verbosity=1
if not hasattr(sp.csr_array,"getnnz"): sp.csr_array.getnnz = sp.csr_array._getnnz
if not hasattr(sp.csc_array,"getnnz"): sp.csc_array.getnnz = sp.csc_array._getnnz
D="data/GSE274828"
pref = sorted({os.path.basename(f).split("_matrix")[0] for f in glob.glob(D+"/*_matrix.mtx.gz")})
ads=[]
for p in pref:
    gsm,tag = p.split("_",1)                       # GSM8458655, YNV
    age = "young" if tag.startswith("Y") else "aged"
    vent = 1 if tag.endswith("V") and not tag.endswith("NV") else 0
    X = sio.mmread(f"{D}/{p}_matrix.mtx.gz").T.tocsr()
    feat = pd.read_csv(f"{D}/{p}_features.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(f"{D}/{p}_barcodes.tsv.gz", sep="\t", header=None)[0].values
    ad = sc.AnnData(X)
    ad.var_names = feat[1].astype(str).values if feat.shape[1]>1 else feat[0].astype(str).values
    ad.var_names_make_unique()
    ad.obs_names=[f"{tag}_{b}" for b in bc]
    ad.obs["sample"]=tag; ad.obs["age"]=age; ad.obs["vent"]=vent; ad.obs["gsm"]=gsm
    ads.append(ad); print(tag, age, "vent="+str(vent), ad.shape, flush=True)
ad = sc.concat(ads, join="inner"); ad.layers["counts"]=ad.X.copy()
print("concat", ad.shape, flush=True)
ad.var["mt"]=ad.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
n0=ad.n_obs
ad = ad[(ad.obs.n_genes_by_counts>=200)&(ad.obs.pct_counts_mt<20)].copy()
sc.pp.filter_genes(ad, min_cells=3)
print(f"QC min_genes>=200 pct_mt<20: {n0} -> {ad.n_obs}", flush=True)
ad.raw=None
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
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
 "Tumor_Epi":["Epcam","Krt8","Krt18","Sftpc","Scgb1a1","Ager"],
 "Fibroblast":["Col1a1","Col1a2","Dcn","Pdgfrb"],
 "Endothelial":["Pecam1","Cdh5","Flt1"],
}
for k,g in M.items():
    gg=[x for x in g if x in ad.var_names]; sc.tl.score_genes(ad, gg, score_name="sc_"+k)
sm = ad.obs.groupby("leiden", observed=True)[["sc_"+k for k in M]].mean()
sm.columns=[c[3:] for c in sm.columns]; lab=sm.idxmax(axis=1)
ad.obs["celltype"]=ad.obs.leiden.map(lab).astype(str)
sm["label"]=lab; sm["n_cells"]=ad.obs.leiden.value_counts()
sm.to_csv("out/GSE274828_cluster_markers.csv"); print(sm.round(3).to_string(), flush=True)
sel = ad.obs.celltype.isin(["Myeloid","Neutrophil"]); sub = ad[sel]
cnt = sub.obs.groupby(["sample","age","vent"], observed=True).size()
print("myeloid cells", int(sel.sum()), "\n", cnt.to_string(), flush=True)
pb={s: np.asarray(sub[sub.obs["sample"]==s].layers["counts"].sum(axis=0)).ravel() for s in sorted(sub.obs["sample"].unique())}
pd.DataFrame(pb, index=sub.var_names).to_csv("out/GSE274828_pseudobulk_myeloid.csv")
cnt.to_frame("n_cells").to_csv("out/GSE274828_ncells_myeloid.csv")
print("SCDONE", flush=True)
