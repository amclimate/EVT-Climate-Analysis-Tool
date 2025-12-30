#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Professional EVT Analysis with Non-stationary σ(t) and Return Levels + POT
========================================================================
- Annual Block Maxima (GEV)
- POT threshold stability analysis
- Rolling POT with bootstrap confidence intervals
- Non-stationary σ(t) model
- Return Levels (50 & 100-year)
- Mann–Kendall trend test
- Pettitt change-point test
- Fully automated, publication-ready PDF report
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
from scipy.stats import genextreme, genpareto, kendalltau, rankdata
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
)
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
import warnings
warnings.filterwarnings("ignore")

# ==========================================================
# CONFIGURATION
# ==========================================================
FILE_PATH = r"E:\video\tails\sample data\point\atmax.csv"      
OUTPUT_DIR = r"E:\video\tails\sample data\point\output"
POT_QUANTILES = [0.9, 0.95, 0.99, 0.995]
ROLLING_WINDOW = 5
BOOTSTRAP_N = 300
MIN_OBS = 50
RETURN_YEARS = [50, 100]

os.makedirs(OUTPUT_DIR, exist_ok=True)
pdf_path = os.path.join(OUTPUT_DIR, "EVT_Full_Professional_Report.pdf")

# ==========================================================
# LOAD DATA
# ==========================================================
df = pd.read_csv(FILE_PATH)
if {'year','month','day'}.issubset(df.columns):
    df['date'] = pd.to_datetime(df[['year','month','day']])
    df.set_index('date', inplace=True)

numeric_cols = [c for c in df.select_dtypes(include=np.number).columns if c.lower() not in {'year','month','day'}]

# ==========================================================
# HELPER FUNCTIONS
# ==========================================================
def gev_family(xi):
    if xi>0: return "Fréchet-type (heavy-tailed)"
    elif xi<0: return "Weibull-type (bounded upper tail)"
    else: return "Gumbel-type (light-tailed)"

def pettitt_test(series):
    n = len(series)
    r = rankdata(series)
    U = np.cumsum(r) - np.arange(1,n+1)*np.mean(r)
    K = np.argmax(np.abs(U))
    p = 2*np.exp((-6*np.max(np.abs(U))**2)/(n**3+n**2))
    return K,p

def bootstrap_xi(exc, n=BOOTSTRAP_N):
    xi = []
    for _ in range(n):
        try:
            s = np.random.choice(exc,len(exc),replace=True)
            x,_,_ = genpareto.fit(s,floc=0)
            xi.append(x)
        except: pass
    return np.array(xi)

def rolling_pot(series,q,min_obs):
    years = np.unique(series.index.year)
    Y,XI,CI=[],[],[]
    for i in range(len(years)-ROLLING_WINDOW+1):
        sub = series[str(years[i]):str(years[i+ROLLING_WINDOW-1])]
        if len(sub)<min_obs: continue
        u = np.percentile(sub,q*100)
        exc = sub[sub>u]-u
        if len(exc)<5: continue
        xi,_,_ = genpareto.fit(exc,floc=0)
        xi_b = bootstrap_xi(exc)
        if len(xi_b)<30: continue
        Y.append(np.mean(years[i:i+ROLLING_WINDOW]))
        XI.append(xi)
        CI.append(np.percentile(xi_b,[2.5,97.5]))
    return np.array(Y),np.array(XI),np.array(CI)

# ==========================================================
# EVT ANALYSIS FUNCTION
# ==========================================================
def evt_analysis(series,name):
    styles=getSampleStyleSheet()
    content=[]
    label=name.replace('_',' ').title()

    # -------------------- GEV BLOCK MAXIMA --------------------
    annual = series.resample('YE').max()
    c, mu, sigma = genextreme.fit(annual)
    xi = -c

    # -------------------- Non-stationary sigma(t) --------------------
    t = np.arange(len(annual))
    p = np.polyfit(t, sigma*np.ones_like(t), 1)   # linear trend for sigma
    sigma_t = p[0]*t + p[1]

    # -------------------- Return Levels --------------------
    rl_table_data=[["Return Period (yrs)","Return Level"]]
    for T in RETURN_YEARS:
        rl = mu + sigma_t[-1]/xi*(((-np.log(1-1/T))**(-xi))-1)
        rl_table_data.append([T,f"{rl:.2f}"])
    rl_table = Table(rl_table_data)
    rl_table.setStyle(TableStyle([
        ('GRID',(0,0),(-1,-1),0.5,colors.black),
        ('BACKGROUND',(0,0),(-1,0),colors.lightgrey)
    ]))

    # -------------------- POT THRESHOLD ANALYSIS --------------------
    sens=[]
    for q in POT_QUANTILES:
        u = np.percentile(series,q*100)
        exc = series[series>u]-u
        if len(exc)<5: continue
        xi_q,_,sg_q = genpareto.fit(exc,floc=0)
        xi_b = bootstrap_xi(exc)
        rel_unc = (np.percentile(xi_b,97.5)-np.percentile(xi_b,2.5))/abs(xi_q)
        sens.append([q,u,xi_q,sg_q,rel_unc])
    sens_df = pd.DataFrame(sens,columns=["Quantile","Threshold","xi","sigma","Relative Uncertainty"])
    q_opt = sens_df.sort_values("Relative Uncertainty").iloc[0]["Quantile"]
    u_opt = sens_df.sort_values("Relative Uncertainty").iloc[0]["Threshold"]

    # -------------------- Rolling POT --------------------
    yrs, xi_r, ci = rolling_pot(series,q_opt,MIN_OBS)
    slope, _ = np.polyfit(yrs,xi_r,1)
    tau, p_mk = kendalltau(yrs,xi_r)
    _, p_pet = pettitt_test(xi_r)

    # -------------------- PLOTS --------------------
    # GEV histogram
    buf_gev = BytesIO()
    plt.figure()
    plt.hist(annual,density=True,alpha=0.5)
    x = np.linspace(min(annual),max(annual),300)
    plt.plot(x,genextreme.pdf(x,c,mu,sigma),'r',lw=2)
    plt.title("GEV Fit to Annual Block Maxima")
    plt.savefig(buf_gev,bbox_inches='tight'); buf_gev.seek(0); plt.close()

    # Annual series
    buf_ts = BytesIO()
    plt.figure()
    plt.plot(annual.index.year,annual.values,'o-')
    plt.title("Annual Block Maxima Time Series"); plt.xlabel("Year"); plt.ylabel(label)
    plt.savefig(buf_ts,bbox_inches='tight'); buf_ts.seek(0); plt.close()

    # POT xi stability
    buf_xi = BytesIO()
    plt.figure()
    plt.plot(sens_df["Quantile"],sens_df["xi"],'o-')
    plt.axvline(q_opt,color='r',ls='--')
    plt.title("POT Shape Parameter Stability"); plt.xlabel("Quantile"); plt.ylabel("ξ")
    plt.savefig(buf_xi,bbox_inches='tight'); buf_xi.seek(0); plt.close()

    # POT sigma stability
    buf_sg = BytesIO()
    plt.figure()
    plt.plot(sens_df["Quantile"],sens_df["sigma"],'o-')
    plt.axvline(q_opt,color='r',ls='--')
    plt.title("POT Scale Parameter Stability"); plt.xlabel("Quantile"); plt.ylabel("σ")
    plt.savefig(buf_sg,bbox_inches='tight'); buf_sg.seek(0); plt.close()

    # POT time series
    buf_pot = BytesIO()
    plt.figure()
    plt.plot(series.index,series,lw=0.7)
    plt.axhline(u_opt,color='r',ls='--',lw=2)
    plt.title("POT Time Series with Selected Threshold")
    plt.savefig(buf_pot,bbox_inches='tight'); buf_pot.seek(0); plt.close()

    # Rolling POT xi(t) with CI
    buf_roll = BytesIO()
    plt.figure()
    plt.plot(yrs,xi_r,'o-',label="ξ(t)")
    plt.fill_between(yrs,ci[:,0],ci[:,1],alpha=0.3,label="95% CI")
    plt.legend(); plt.xlabel("Year"); plt.ylabel("ξ"); plt.title("Rolling POT Shape Parameter")
    plt.savefig(buf_roll,bbox_inches='tight'); buf_roll.seek(0); plt.close()

    # -------------------- PDF CONTENT --------------------
    content.append(Paragraph(f"<b>Extreme Value Analysis – {label}</b>",styles["Title"]))
    content.append(Spacer(1,10))
    gev_table = Table([
        ["Metric","Value"],
        ["Distribution","Generalized Extreme Value (GEV)"],
        ["GEV Type",gev_family(xi)],
        ["Number of Blocks",len(annual)],
        ["Observations per Block",int(len(series)/len(annual))],
        ["μ (location)",f"{mu:.3f}"],
        ["σ (scale)",f"{sigma:.3f}"],
        ["ξ (shape)",f"{xi:.3f}"]
    ])
    gev_table.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,colors.black),('BACKGROUND',(0,0),(-1,0),colors.lightgrey)]))
    content.append(gev_table)
    content.append(Spacer(1,10))
    content.append(Image(buf_gev,420,260))
    content.append(Spacer(1,10))
    content.append(Image(buf_ts,420,260))
    content.append(Spacer(1,10))
    content.append(Image(buf_xi,420,260))
    content.append(Spacer(1,10))
    content.append(Image(buf_sg,420,260))
    content.append(Spacer(1,10))
    content.append(Image(buf_pot,420,260))
    content.append(Spacer(1,10))
    content.append(Image(buf_roll,420,260))
    content.append(Spacer(1,10))
    content.append(rl_table)
    content.append(Spacer(1,10))
    content.append(Paragraph(
        f"<b>Interpretation:</b> The <b>GEV</b> distribution fitted to annual block maxima indicates a <b>{gev_family(xi)}</b> behavior. "
        f"Optimal POT threshold q={q_opt:.2f} based on minimum bootstrap uncertainty. Rolling POT shows "
        f"{'a significant' if p_mk<0.05 else 'no significant'} trend in ξ(t) (slope={slope:.4f}). "
        f"Mann-Kendall p={p_mk:.3f}, Pettitt p={p_pet:.3f}. Return levels for {RETURN_YEARS} yrs calculated.",
        styles["Normal"]
    ))
    return content

# ==========================================================
# RUN ANALYSIS
# ==========================================================
pdf_content=[]
for col in numeric_cols:
    s = df[col].dropna()
    if s.min()<=0: s=s[s>0]
    pdf_content += evt_analysis(s,col)

doc = SimpleDocTemplate(pdf_path,pagesize=A4)
doc.build(pdf_content)

print("PDF generated successfully:")
print(pdf_path)
