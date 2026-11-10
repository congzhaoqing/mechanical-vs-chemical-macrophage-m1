import os,glob,gzip,numpy as np,pandas as pd,scanpy as sc,scipy.io as sio,scipy.sparse as sp
sc.settings.verbosity=0
for a in (sp.csr_array,sp.csc_array):
    if not hasattr(a,"getnnz"): a.getnnz=a._getnnz
MK={"Myeloid":["Lyz2","Csf1r","Itgam","Adgre1","Cd68","Fcgr1","Ctss","Tyrobp","Fcer1g"],
 "Neutrophil":["S100a8","S100a9","Retnlg","Ly6g"],
 "T_NK":["Cd3e","Cd3d","Cd3g","Trbc2","Nkg7","Gzmb","Klrb1c"],
 "B":["Cd79a","Cd79b","Ms4a1","Ighm"],
 "Tumor_Epi":["Epcam","Krt8","Krt18","Mlana","Pmel","Dct","Tyr"],
 "Fibroblast":["Col1a1","Col1a2","Dcn","Pdgfrb"],
 "Endothelial":["Pecam1","Cdh5","Flt1"],
 "Epithelial_lung":["Sftpc","Sftpb","Scgb1a1","Ager"]}
def pipeline(ad,tag):
    ad.layers["counts"]=ad.X.copy()
    ad.var["mt"]=ad.var_names.str.startswith(("mt-","MT-"))
    sc.pp.calculate_qc_metrics(ad,qc_vars=["mt"],inplace=True,percent_top=None,log1p=False)
    n0=ad.n_obs
    ad=ad[(ad.obs.n_genes_by_counts>=200)&(ad.obs.pct_counts_mt<20)].copy()
    sc.pp.filter_genes(ad,min_cells=3)
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad,n_top_genes=2000)
    sc.pp.pca(ad,n_comps=min(50,ad.n_obs-1,ad.n_vars-1),mask_var="highly_variable")
    sc.pp.neighbors(ad,n_neighbors=15,n_pcs=30)
    sc.tl.leiden(ad,resolution=1.0,flavor="igraph",n_iterations=2,directed=False)
    for k,g in MK.items():
        gg=[x for x in g if x in ad.var_names]
        if gg: sc.tl.score_genes(ad,gg,score_name="sc_"+k)
    cols=[c for c in ad.obs.columns if c.startswith("sc_")]
    smm=ad.obs.groupby("leiden",observed=True)[cols].mean(); smm.columns=[c[3:] for c in smm.columns]
    ad.obs["celltype"]=ad.obs.leiden.map(smm.idxmax(axis=1)).astype(str)
    sel=ad.obs.celltype.isin(["Myeloid","Neutrophil"])
    print(f"{tag}: {n0} -> {ad.n_obs} cells QC; myeloid {int(sel.sum())}",flush=True)
    sub=ad[sel]
    pb={s:np.asarray(sub[sub.obs['sample']==s].layers["counts"].sum(axis=0)).ravel()
        for s in sorted(sub.obs['sample'].unique())}
    pd.DataFrame(pb,index=sub.var_names).to_csv(f"out/{tag}_pseudobulk_myeloid.csv")
    sub.obs.groupby(["sample"],observed=True).size().to_frame("n_cells").to_csv(f"out/{tag}_ncells.csv")

# GSE275018 — 9 samples, PBS/MNM/RMF
D="sc/GSE275018"
pref=sorted({os.path.basename(f).split("-matrix")[0] for f in glob.glob(D+"/*-matrix.mtx.gz")})
ads=[]
for p in pref:
    tag=p.split("_",1)[1]
    X=sio.mmread(f"{D}/{p}-matrix.mtx.gz").T.tocsr()
    feat=pd.read_csv(f"{D}/{p}-features.tsv.gz",sep="\t",header=None)
    bc=pd.read_csv(f"{D}/{p}-barcodes.tsv.gz",sep="\t",header=None)[0].values
    ad=sc.AnnData(X); ad.var_names=feat[1].astype(str).values if feat.shape[1]>1 else feat[0].astype(str).values
    ad.var_names_make_unique(); ad.obs_names=[f"{tag}_{b}" for b in bc]; ad.obs["sample"]=tag
    ads.append(ad)
pipeline(sc.concat(ads,join="inner"),"GSE275018")
