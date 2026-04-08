import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import combinations
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')



# ── 1. Parse SEF files ──────────────────────────────────────────────────────

def parse_sef(filepath):
    """Parse a SEF .tsv file, return metadata dict and daily-mean DataFrame."""
    meta = {}
    header_line = None
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith('Year'):
                header_line = i
                break
            parts = line.split('\t', 1)
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    
    meta['Lat'] = float(meta.get('Lat', 'nan'))
    meta['Lon'] = float(meta.get('Lon', 'nan'))
    meta['Alt'] = float(meta.get('Alt', 'nan'))
    
    df = pd.read_csv(filepath, sep='\t', skiprows=header_line,
                     na_values=['NA'], encoding='utf-8', encoding_errors='replace')
    
    # Build date
    df['date'] = pd.to_datetime(
        df[['Year','Month','Day']].rename(columns={'Year':'year','Month':'month','Day':'day'}),
        errors='coerce'
    )
    df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
    df = df.dropna(subset=['date','Value'])
    
    # Daily mean (handles both daily and subdaily)
    daily = df.groupby('date')['Value'].mean().reset_index()
    daily.columns = ['date', 'ta']
    daily = daily.set_index('date').sort_index()
    
    return meta, daily


# ── 2. Load all stations ────────────────────────────────────────────────────

target_files = [
    'Marschlins/CHIMES_Marschlins_17820101-18710731_ta_subdaily_qc.tsv',
    'Rovereto/PALAEO-RA_Rovereto_Bonfioli_17820214-18390827_ta_subdaily_qc.tsv',
    'Turin/SMI_Turin_17530101-18651130_ta_subdaily_qc.tsv',
    'Milan/Maugeri-et-al_Milan_17630101-19491231_ta_daily_qc.tsv',
    'Bologna/Camuffo_Bologna_17150101-18151231_ta_daily_qc.tsv',
    'Padua/IMPROVE_Padua_17250112-19970531_ta_daily_qc.tsv',
]
files = [Path('/scratch3/PALAEO-RA/daily_data/final/') / f for f in target_files]
stations = {}
for f in files:
    try:
        meta, daily = parse_sef(f)
        name = meta.get('Name', f.stem)
        stations[name] = {'meta': meta, 'daily': daily, 'file': f.name}
        print(f"  {name}: {len(daily)} days, {daily.index.min().year}-{daily.index.max().year}, "
              f"alt={meta['Alt']}m, ({meta['Lat']:.2f}, {meta['Lon']:.2f})")
    except Exception as e:
        print(f"  SKIP {f.name}: {e}")

# ── 3. Compute monthly anomalies ────────────────────────────────────────────

def to_monthly_anomalies(daily_df):
    """Monthly mean → remove climatology → anomalies."""
    monthly = daily_df.resample('MS').agg(['mean', 'count'])
    monthly.columns = ['ta', 'count']
    # Require at least 20 days per month
    monthly.loc[monthly['count'] < 20, 'ta'] = np.nan
    monthly = monthly[['ta']].dropna()
    monthly['month'] = monthly.index.month
    clim = monthly.groupby('month')['ta'].mean()
    monthly['anom'] = monthly['ta'] - monthly['month'].map(clim)
    return monthly[['ta', 'anom']]

for name in stations:
    stations[name]['monthly'] = to_monthly_anomalies(stations[name]['daily'])

# ── 4. Pairwise comparison ──────────────────────────────────────────────────

def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat/2)**2 + np.cos(np.radians(lat1))*np.cos(np.radians(lat2))*np.sin(dlon/2)**2
    return R * 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))

pair_results = []
for (n1, s1), (n2, s2) in combinations(stations.items(), 2):
    dist = haversine(s1['meta']['Lat'], s1['meta']['Lon'],
                     s2['meta']['Lat'], s2['meta']['Lon'])
    
    # Merge on overlapping months
    merged = s1['monthly'][['anom']].join(s2['monthly'][['anom']], 
                                          lsuffix='_1', rsuffix='_2', how='inner')
    merged = merged.dropna()
    
    if len(merged) < 24:
        continue
    
    corr = merged['anom_1'].corr(merged['anom_2'])
    rmse = np.sqrt(((merged['anom_1'] - merged['anom_2'])**2).mean())
    bias = (merged['anom_1'] - merged['anom_2']).mean()
    
    pair_results.append({
        'station1': n1, 'station2': n2,
        'dist_km': dist,
        'corr': corr, 'rmse': rmse, 'bias': bias,
        'n_months': len(merged),
        'overlap_start': merged.index.min(),
        'overlap_end': merged.index.max(),
        'lat1': s1['meta']['Lat'], 'lon1': s1['meta']['Lon'],
        'lat2': s2['meta']['Lat'], 'lon2': s2['meta']['Lon'],
        'alt1': s1['meta']['Alt'], 'alt2': s2['meta']['Alt'],
    })

pairs_df = pd.DataFrame(pair_results).sort_values('dist_km')
print("\n" + "="*80)
print("PAIRWISE RESULTS (sorted by distance)")
print("="*80)
for _, r in pairs_df.iterrows():
    print(f"\n  {r['station1']} ({r['alt1']:.0f}m) <-> {r['station2']} ({r['alt2']:.0f}m)")
    print(f"    dist={r['dist_km']:.0f}km  corr={r['corr']:.3f}  RMSE={r['rmse']:.2f}°C  "
          f"bias={r['bias']:+.2f}°C  N={r['n_months']} months")
    print(f"    overlap: {r['overlap_start'].strftime('%Y-%m')} → {r['overlap_end'].strftime('%Y-%m')}")


# ── 5. FIGURES ───────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.3)

# --- Panel A: Map with station locations and connecting lines colored by corr ---
ax_map = fig.add_subplot(gs[0, :])

# Draw lines between pairs, colored by correlation
cmap = plt.cm.RdYlGn
for _, r in pairs_df.iterrows():
    color = cmap((r['corr'] - 0.3) / 0.7)  # map 0.3-1.0 to colormap
    lw = max(1, 4 * r['corr'])
    ax_map.plot([r['lon1'], r['lon2']], [r['lat1'], r['lat2']],
                color=color, linewidth=lw, alpha=0.7)

# Plot stations
for name, s in stations.items():
    ax_map.plot(s['meta']['Lon'], s['meta']['Lat'], 'ko', markersize=8, zorder=5)
    ax_map.text(s['meta']['Lon'] + 0.15, s['meta']['Lat'] + 0.12,
                f"{name}\n({s['meta']['Alt']:.0f}m)",
                fontsize=7, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.3, 1.0))
cb = plt.colorbar(sm, ax=ax_map, orientation='horizontal', pad=0.08, shrink=0.6, aspect=30)
cb.set_label('Monthly anomaly correlation', fontsize=10)
ax_map.set_xlabel('Longitude')
ax_map.set_ylabel('Latitude')
ax_map.set_aspect('auto')
ax_map.grid(True, alpha=0.2)
ax_map.set_title('(a) Station network – pairwise anomaly correlations', fontsize=12, fontweight='bold')

# --- Panel B: Correlation vs distance ---
ax_corr = fig.add_subplot(gs[1, 0])
for _, r in pairs_df.iterrows():
    marschlins_pair = 'Marschlins' in r['station1'] or 'Marschlins' in r['station2'] or \
                      'Schloss' in r['station1'] or 'Schloss' in r['station2']
    marker = '^' if marschlins_pair else 'o'
    color = 'orange' if marschlins_pair else 'steelblue'
    ax_corr.scatter(r['dist_km'], r['corr'], s=80, c=color, marker=marker,
                    edgecolor='k', linewidth=0.5, zorder=3)
    label = f"{r['station1'][:3]}-{r['station2'][:3]}"
    ax_corr.annotate(label, (r['dist_km'], r['corr']), fontsize=6,
                     xytext=(5, 5), textcoords='offset points')

ax_corr.set_xlabel('Distance (km)')
ax_corr.set_ylabel('Correlation (monthly anomalies)')
ax_corr.set_title('(b) Correlation vs. distance', fontweight='bold')
ax_corr.set_ylim(0.4, 1.02)
ax_corr.grid(True, alpha=0.3)
# Legend for Marschlins
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='steelblue', markersize=8, label='Low-altitude pairs'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='orange', markersize=8, label='Marschlins pairs (562m)'),
]
ax_corr.legend(handles=legend_elements, fontsize=8)

# --- Panel C: Example time series comparison (closest pair) ---
ax_ts = fig.add_subplot(gs[1, 1])
# Pick the pair with highest correlation
best = pairs_df.sort_values('corr', ascending=False).iloc[0]
s1 = stations[best['station1']]['monthly']
s2 = stations[best['station2']]['monthly']
merged = s1[['anom']].join(s2[['anom']], lsuffix='_1', rsuffix='_2', how='inner').dropna()

# Plot a 10-year window
mid = merged.index[len(merged)//2]
window = merged.loc[str(mid.year - 5):str(mid.year + 5)]

ax_ts.plot(window.index, window['anom_1'], '-', color='steelblue', alpha=0.8,
           label=f"{best['station1']}", linewidth=1)
ax_ts.plot(window.index, window['anom_2'], '-', color='firebrick', alpha=0.8,
           label=f"{best['station2']}", linewidth=1)
ax_ts.set_ylabel('Monthly anomaly (°C)')
ax_ts.set_title(f"(c) Best pair: {best['station1']}–{best['station2']} "
                f"(r={best['corr']:.3f}, {best['dist_km']:.0f} km)", fontweight='bold')
ax_ts.legend(fontsize=8)
ax_ts.grid(True, alpha=0.3)
fig.autofmt_xdate()

plt.savefig('image/neighbor_validation.png', dpi=200, bbox_inches='tight')
plt.savefig('image/neighbor_validation.pdf', bbox_inches='tight')
print("\n✓ Figures saved")
