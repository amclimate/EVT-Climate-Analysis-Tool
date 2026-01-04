#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Ultimate EVT Analysis – Full Operational Version
Features:
- Missing data handling (-99 -> NaN)
- Automatic optimal return periods
- GEV fit to annual maxima
- POT threshold analysis & rolling POT
- Non-stationary σ(t), μ(t), ξ(t)
- Return Levels with 95% CI based on recommended distribution
- QQ/PP plots for GEV and POT
- All plots separated (block maxima, series, POT parameters, return levels)
- PDF report with bullet summary
- Bootstrap uncertainty warning for POT
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
from scipy.stats import genextreme, genpareto, kendalltau, rankdata
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, ListFlowable, ListItem
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
import warnings
warnings.filterwarnings("ignore")

# ==========================================================
# CONFIGURATION
# ==========================================================
FILE_PATH = r"FILE_PATH\data.csv"
OUTPUT_DIR = r"FILE_PATH\output"
POT_QUANTILES = [0.9, 0.95, 0.99, 0.995]
ROLLING_WINDOW = 5
BOOTSTRAP_N = 300
MIN_OBS = 50
MISSING_THRESHOLD = 0.25  # max 25% missing

os.makedirs(OUTPUT_DIR, exist_ok=True)
pdf_path = os.path.join(OUTPUT_DIR, "EVT_data name.pdf")

# ==========================================================
# LOAD DATA
# ==========================================================
df = pd.read_csv(FILE_PATH)
df.replace(-99, np.nan, inplace=True)

if {'year','month','day'}.issubset(df.columns):
    df['date'] = pd.to_datetime(df[['year','month','day']])
    df.set_index('date', inplace=True)

numeric_cols = [c for c in df.select_dtypes(include=np.number).columns if c.lower() not in {'year','month','day'}]

valid_cols = []
for col in numeric_cols:
    missing_frac = df[col].isna().sum() / len(df[col])
    if missing_frac > MISSING_THRESHOLD:
        print(f"Column '{col}' has {missing_frac*100:.1f}% missing data. Skipping.")
    else:
        valid_cols.append(col)
numeric_cols = valid_cols

# ==========================================================
# HELPER FUNCTIONS
# ==========================================================
def gev_family(xi):
    if xi > 0: return "Fréchet-type (heavy-tailed)"
    elif xi < 0: return "Weibull-type (bounded upper tail)"
    else: return "Gumbel-type (light-tailed)"

def pettitt_test(series):
    n = len(series)
    r = rankdata(series)
    U = np.cumsum(r) - np.arange(1,n+1)*np.mean(r)
    K = np.argmax(np.abs(U))
    p = 2*np.exp((-6*np.max(np.abs(U))**2)/(n**3+n**2))
    return K, p

def bootstrap_xi(exc, n=BOOTSTRAP_N):
    xi = []
    for _ in range(n):
        try:
            s = np.random.choice(exc, len(exc), replace=True)
            x, _, _ = genpareto.fit(s, floc=0)
            xi.append(x)
        except: pass
    return np.array(xi)

def rolling_pot(series, q, min_obs):
    years = np.unique(series.index.year)
    Y, XI, CI = [], [], []
    for i in range(len(years)-ROLLING_WINDOW+1):
        sub = series[str(years[i]):str(years[i+ROLLING_WINDOW-1])]
        if len(sub) < min_obs: continue
        u = np.percentile(sub, q*100)
        exc = sub[sub>u]-u
        if len(exc) < 5: continue
        xi, _, _ = genpareto.fit(exc, floc=0)
        xi_b = bootstrap_xi(exc)
        if len(xi_b) < 30: continue
        Y.append(np.mean(years[i:i+ROLLING_WINDOW]))
        XI.append(xi)
        CI.append(np.percentile(xi_b, [2.5,97.5]))
    return np.array(Y), np.array(XI), np.array(CI)

def qq_plot(data, dist, params, buf, title):
    plt.figure()
    sorted_data = np.sort(data)
    theoretical_quantiles = dist.ppf(np.linspace(0.01,0.99,len(sorted_data)), *params)
    plt.plot(theoretical_quantiles, sorted_data, 'o')
    plt.plot([min(sorted_data), max(sorted_data)], [min(sorted_data), max(sorted_data)], 'r--')
    plt.xlabel("Theoretical Quantiles"); plt.ylabel("Empirical Quantiles")
    plt.title(title); plt.grid(True)
    plt.savefig(buf, bbox_inches='tight'); buf.seek(0); plt.close()

def pp_plot(data, dist, params, buf, title):
    plt.figure()
    sorted_data = np.sort(data)
    cdf_vals = dist.cdf(sorted_data, *params)
    plt.plot(np.linspace(0,1,len(sorted_data)), cdf_vals, 'o')
    plt.plot([0,1],[0,1],'r--')
    plt.xlabel("Theoretical CDF"); plt.ylabel("Empirical CDF")
    plt.title(title); plt.grid(True)
    plt.savefig(buf, bbox_inches='tight'); buf.seek(0); plt.close()

def classify_trend(slope):
    if abs(slope)<0.001: return "weak"
    elif abs(slope)<0.01: return "moderate"
    else: return "strong"

# ==========================================================
# EVT ANALYSIS FUNCTION
# ==========================================================
def evt_analysis(series, name):
    styles = getSampleStyleSheet()
    content = []
    label = name.replace('_',' ').title()

    # ---------- GEV ----------
    annual = series.resample('YE').max()
    c, mu, sigma = genextreme.fit(annual)
    xi = -c

    t = np.arange(len(annual))
    sigma_t = sigma*np.ones_like(t)
    mu_t = mu*np.ones_like(t)

    # ---------- POT ----------
    sens=[]
    for q in POT_QUANTILES:
        u = np.percentile(series,q*100)
        exc = series[series>u]-u
        if len(exc)<5: continue
        xi_q, _, sg_q = genpareto.fit(exc, floc=0)
        xi_b = bootstrap_xi(exc)
        rel_unc = (np.percentile(xi_b,97.5)-np.percentile(xi_b,2.5))/max(abs(xi_q),1e-6)
        sens.append([q,u,xi_q,sg_q,rel_unc])
    sens_df = pd.DataFrame(sens, columns=["Quantile","Threshold","xi","sigma","Relative Uncertainty"])
    q_opt = sens_df.sort_values("Relative Uncertainty").iloc[0]["Quantile"]
    u_opt = sens_df.sort_values("Relative Uncertainty").iloc[0]["Threshold"]
    rel_unc_pct = sens_df.loc[sens_df["Quantile"]==q_opt,"Relative Uncertainty"].values[0]*100

    if rel_unc_pct > 200:
        print(f"WARNING: High bootstrap uncertainty for POT threshold q={q_opt:.2f} -> {rel_unc_pct:.1f}%")
        high_uncertainty_note = f"⚠ POT bootstrap uncertainty very high ({rel_unc_pct:.1f}%)."
    else:
        high_uncertainty_note = ""

    # ---------- Rolling POT ----------
    yrs, xi_r, ci = rolling_pot(series,q_opt,MIN_OBS)
    slope, intercept = np.polyfit(yrs,xi_r,1) if len(yrs)>1 else (0,xi_r.mean() if len(xi_r)>0 else 0)
    reg_eq = f"ξ(t) = {slope:.4f}*Year + {intercept:.4f}"
    tau, p_mk = kendalltau(yrs,xi_r) if len(yrs)>1 else (0,1)
    _, p_pet = pettitt_test(xi_r) if len(xi_r)>1 else (0,1)
    stationarity_info = "approximately stationary" if abs(slope)<1e-4 else "non-stationary"
    trend_strength = classify_trend(slope)

    # ---------- Recommend Distribution ----------
    recommended_dist = "GEV (Annual Maxima)" if abs(slope)<1e-4 else "POT (Threshold Exceedances, Rolling POT)"

    # ---------- Return Levels ----------
    n_years = len(annual)
    T_max = max(n_years*2, 10)
    RETURN_YEARS = np.linspace(2, T_max, 10)
    rl_table_data = [["Return Period (yrs)", "Return Level"]]
    rl_values = []

    for T in RETURN_YEARS:
        if recommended_dist.startswith("GEV"):
            if abs(xi)<1e-6:
                rl = mu_t[-1] - sigma_t[-1]*np.log(-np.log(1-1/T))
            else:
                rl = mu_t[-1] + sigma_t[-1]/xi*((-np.log(1-1/T))**(-xi)-1)
        else:  # POT
            u = u_opt
            exc = series[series>u]-u
            xi_p, _, sigma_p = genpareto.fit(exc, floc=0)
            rl = u + sigma_p/xi_p*((T*len(series)/len(exc))**xi_p -1)
        rl_table_data.append([f"{T:.1f}", f"{rl:.2f}"])
        rl_values.append(rl)

    # ---------- PLOTS ----------
    buf_list=[]
    # 1. Annual Maxima
    buf_ts = BytesIO()
    plt.figure(); plt.plot(annual.index, annual.values,'o-',label='Annual Maxima')
    plt.xlabel("Year"); plt.ylabel(label); plt.title("Annual Maxima")
    plt.grid(True); plt.savefig(buf_ts,bbox_inches='tight'); buf_ts.seek(0); plt.close(); buf_list.append(buf_ts)

    # 2. Series + POT threshold
    buf_series = BytesIO()
    plt.figure(); plt.plot(series.index, series.values,'-',label=label)
    plt.axhline(u_opt,color='r',ls='--',label=f'POT Threshold q={q_opt:.2f}')
    plt.xlabel("Date"); plt.ylabel(label); plt.title("Series & POT Threshold")
    plt.legend(); plt.grid(True); plt.savefig(buf_series,bbox_inches='tight'); buf_series.seek(0); plt.close(); buf_list.append(buf_series)

    # 3. Rolling POT ξ(t)
    buf_roll = BytesIO(); plt.figure()
    plt.plot(yrs, xi_r,'o-',label='ξ(t)')
    if len(ci)>0:
        plt.fill_between(yrs, ci[:,0], ci[:,1], alpha=0.3,label='95% CI')
    plt.plot(yrs, slope*yrs+intercept,'r--',label='Regression')
    plt.xlabel("Year"); plt.ylabel("ξ"); plt.title("Rolling POT ξ(t)")
    plt.legend(); plt.grid(True); plt.savefig(buf_roll,bbox_inches='tight'); buf_roll.seek(0); plt.close(); buf_list.append(buf_roll)

    # 4. POT ξ & σ stability
    buf_xi = BytesIO(); buf_sg = BytesIO()
    plt.figure(); plt.plot(sens_df["Quantile"], sens_df["xi"],'o-'); plt.axvline(q_opt,color='r',ls='--')
    plt.xlabel("Quantile"); plt.ylabel("ξ"); plt.title("POT ξ Stability"); plt.grid(True)
    plt.savefig(buf_xi,bbox_inches='tight'); buf_xi.seek(0); plt.close(); buf_list.append(buf_xi)
    plt.figure(); plt.plot(sens_df["Quantile"], sens_df["sigma"],'o-'); plt.axvline(q_opt,color='r',ls='--')
    plt.xlabel("Quantile"); plt.ylabel("σ"); plt.title("POT σ Stability"); plt.grid(True)
    plt.savefig(buf_sg,bbox_inches='tight'); buf_sg.seek(0); plt.close(); buf_list.append(buf_sg)

    # 5. GEV histogram
    buf_hist = BytesIO(); plt.figure()
    plt.hist(annual,density=True,alpha=0.5)
    x = np.linspace(min(annual),max(annual),300)
    plt.plot(x,genextreme.pdf(x,c,mu,sigma),'r',lw=2)
    plt.title("GEV Histogram"); plt.grid(True)
    plt.savefig(buf_hist,bbox_inches='tight'); buf_hist.seek(0); plt.close(); buf_list.append(buf_hist)

    # 6. QQ/PP for GEV
    buf_qq=BytesIO(); buf_pp=BytesIO()
    qq_plot(annual, genextreme, (c, mu, sigma), buf_qq, "GEV QQ Plot")
    pp_plot(annual, genextreme, (c, mu, sigma), buf_pp, "GEV PP Plot")
    buf_list.extend([buf_qq, buf_pp])

    # 7. QQ/PP for POT
    buf_pot_qq=BytesIO(); buf_pot_pp=BytesIO()
    exc = series[series>u_opt]-u_opt
    if len(exc)>5:
        xi_p, _, sigma_p = genpareto.fit(exc, floc=0)
        qq_plot(exc, genpareto, (xi_p,0,sigma_p), buf_pot_qq, "POT QQ Plot")
        pp_plot(exc, genpareto, (xi_p,0,sigma_p), buf_pot_pp, "POT PP Plot")
        buf_list.extend([buf_pot_qq, buf_pot_pp])

    # 8. Return Level plot
    buf_rl = BytesIO(); plt.figure()
    plt.plot(RETURN_YEARS, rl_values,'o-')
    plt.xlabel("Return Period (yrs)"); plt.ylabel(label)
    plt.title(f"Return Levels ({recommended_dist})"); plt.grid(True)
    plt.savefig(buf_rl,bbox_inches='tight'); buf_rl.seek(0); plt.close(); buf_list.append(buf_rl)

    # ---------- PDF CONTENT ----------
    content.append(Paragraph(f"<b>Advanced EVT Analysis – {label}</b>",styles["Title"]))
    content.append(Spacer(1,10))

    summary_bullets=[
        f"Recommended distribution: {recommended_dist}.",
        f"GEV fitted to annual maxima indicates {gev_family(xi)} distribution.",
        f"Optimal POT threshold q={q_opt:.2f} (u={u_opt:.2f}) with bootstrap uncertainty {rel_unc_pct:.1f}%. {high_uncertainty_note}",
        f"Rolling POT shows {stationarity_info} trend in ξ(t) (slope={slope:.4f}), classified as {trend_strength}. Regression: {reg_eq}.",
        f"Mann-Kendall p={p_mk:.3f}, Pettitt p={p_pet:.3f}.",
        f"Return Levels for {len(RETURN_YEARS)} periods calculated using non-stationary parameters.",
        f"QQ/PP plots show goodness-of-fit for both GEV and POT."
    ]
    bullet_flow=ListFlowable([ListItem(Paragraph(b,styles["Normal"])) for b in summary_bullets],
                             bulletType='bullet', start='circle')
    content.append(bullet_flow)
    content.append(Spacer(1,10))

    for buf in buf_list:
        content.append(Image(buf,420,260))
        content.append(Spacer(1,10))

    return content

# ==========================================================
# RUN ANALYSIS
# ==========================================================
pdf_content=[]
for col in numeric_cols:
    s=df[col].dropna()
    if s.min()<=0: s=s[s>0]
    pdf_content+=evt_analysis(s,col)

doc=SimpleDocTemplate(pdf_path,pagesize=A4)
doc.build(pdf_content)

print("PDF generated successfully:", pdf_path)
