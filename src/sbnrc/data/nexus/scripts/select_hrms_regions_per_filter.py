#!/usr/bin/env python3
"""Select sparse HRMS regions independently for each filter."""
import argparse, json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt
from scipy.ndimage import binary_dilation, gaussian_filter

def product(c,f,s):
    m=sorted(c.glob(f"JWST_NRC_*_{f}_60mas_{s}.fits"))
    if len(m)!=1: raise RuntimeError(f"Expected one {f} {s}, found {len(m)}")
    return m[0]

def filters(c):
    out=[p.name.split("_60mas_sci.fits")[0].split("_")[-1] for p in sorted(c.glob("JWST_NRC_*_60mas_sci.fits"))]
    if not out: raise RuntimeError(f"No 60mas SCI products in {c}")
    return out

def nonoverlap(items,n,d):
    out=[]
    for q in sorted(items,key=lambda x:x["score"]):
        y,x=q["center"]
        if all(np.hypot(y-r["center"][0],x-r["center"][1])>=d for r in out): out.append(q)
        if len(out)==n: break
    return out

def preview(a,items,box,path):
    rows=int(np.ceil(len(items)/5)); fig,axes=plt.subplots(rows,5,figsize=(20,3.5*rows)); axes=np.atleast_1d(axes).ravel()
    finite=a[np.isfinite(a)]; sample=finite[::max(1,finite.size//500000)]; lo,hi=np.nanpercentile(sample,[2,98])
    for i,(ax,q) in enumerate(zip(axes,items)):
        y,x=q["corner"]; ax.imshow(a[y:y+box,x:x+box],origin="lower",cmap="gray_r",vmin=lo,vmax=hi)
        ax.set_title(f"{i}: score={q['score']:.3f}\nsrc={q['occupancy']:.3f}"); ax.set_xticks([]); ax.set_yticks([])
    for ax in axes[len(items):]: ax.axis("off")
    fig.tight_layout(); fig.savefig(path,dpi=150); plt.close(fig)

def select(c,f,box,stride,target,minimum,outdir):
    sci=np.asarray(fits.getdata(product(c,f,"sci"),memmap=True)); rms=np.asarray(fits.getdata(product(c,f,"rms"),memmap=True)); wht=np.asarray(fits.getdata(product(c,f,"wht"),memmap=True))
    valid=np.isfinite(sci)&np.isfinite(rms)&(rms>0)&np.isfinite(wht)&(wht>0)
    sample=sci[::8,::8][valid[::8,::8]]; _,med,std=sigma_clipped_stats(sample,sigma=3,maxiters=8)
    typical=np.nanmedian(rms[valid]); floor=min(typical,std) if std>0 else typical
    sn=np.where(valid,(sci-med)/np.maximum(rms,floor),0.0); smooth=gaussian_filter(sn,sigma=1.2)
    source=binary_dilation((sn>2.5)|(smooth>1.8),iterations=5); half=box//2; candidates=[]
    for cy in range(half,sci.shape[0]-half,stride):
        for cx in range(half,sci.shape[1]-half,stride):
            ys=slice(cy-half,cy-half+box); xs=slice(cx-half,cx-half+box); good=valid[ys,xs]; frac=float(good.mean())
            if frac<0.995: continue
            occ=float(source[ys,xs].mean()); vals=sn[ys,xs][good]; vmed=np.nanmedian(vals)
            scatter=float(1.4826*np.nanmedian(np.abs(vals-vmed))); tail=float(np.nanpercentile(vals,99)); score=12*occ+0.35*scatter+0.025*max(tail,0)
            candidates.append({"center":[cy,cx],"corner":[cy-half,cx-half],"score":score,"occupancy":occ,"scatter":scatter,"tail99":tail,"valid_fraction":frac})
    ranked=nonoverlap(candidates,20,box*1.15); chosen=ranked[:min(target,len(ranked))]
    if len(chosen)<minimum: raise RuntimeError(f"{f}: only {len(chosen)} valid boxes; minimum {minimum}")
    preview(sci,ranked,box,outdir/f"candidate_regions_{f}.png"); preview(sci,chosen,box,outdir/f"selected_regions_{f}.png")
    return {"box_count":len(chosen),"corners":[q["corner"] for q in chosen],"selected_metrics":chosen,"candidate_metrics":ranked}

def main():
    p=argparse.ArgumentParser(); p.add_argument("--collection",type=Path,default=Path("collection/60mas")); p.add_argument("--box-size",type=int,default=301); p.add_argument("--stride",type=int,default=75); p.add_argument("--target",type=int,default=10); p.add_argument("--minimum",type=int,default=5); p.add_argument("--output-dir",type=Path,default=Path(".")); a=p.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    payload={"method":"per-filter sparse-region fallback","box_size_60mas":a.box_size,"target_count":a.target,"minimum_count":a.minimum,"filters":{}}
    for f in filters(a.collection): print(f"Selecting {f}",flush=True); payload["filters"][f]=select(a.collection,f,a.box_size,a.stride,a.target,a.minimum,a.output_dir)
    (a.output_dir/"selected_regions_per_filter.json").write_text(json.dumps(payload,indent=2)); print(json.dumps({f:q["box_count"] for f,q in payload["filters"].items()},sort_keys=True))
if __name__=="__main__": main()
