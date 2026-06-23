import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.metrics import RocCurveDisplay
import matplotlib.pyplot as plt
import shap

np.random.seed(42)
N = 500

variants = {
    'rs2229774':  {'gene': 'RARG',   'maf': 0.08,  'effect': 'risk',       'or': 3.8},
    'rs7853758':  {'gene': 'SLC28A3','maf': 0.20,  'effect': 'protective', 'or': 0.5},
    'rs17863783': {'gene': 'UGT1A6', 'maf': 0.12,  'effect': 'risk',       'or': 2.1},
    'rs1137233':  {'gene': 'RAC2',   'maf': 0.15,  'effect': 'risk',       'or': 1.8},
    'rs9024':     {'gene': 'CBR1',   'maf': 0.30,  'effect': 'risk',       'or': 1.4},
    'rs2232228':  {'gene': 'HAS3',   'maf': 0.25,  'effect': 'risk',       'or': 1.5},
    'rs17222723': {'gene': 'ABCC2',  'maf': 0.10,  'effect': 'risk',       'or': 1.3},
    'rs45511401': {'gene': 'ABCC1',  'maf': 0.18,  'effect': 'protective', 'or': 0.7},
    'rs1695':     {'gene': 'GSTP1',  'maf': 0.35,  'effect': 'risk',       'or': 1.2},
    'rs20572':    {'gene': 'CBR1',   'maf': 0.22,  'effect': 'risk',       'or': 1.3},
    'rs1533682':  {'gene': 'ABCC5',  'maf': 0.28,  'effect': 'protective', 'or': 0.8},
    'rs1128503':  {'gene': 'ABCB1',  'maf': 0.45,  'effect': 'risk',       'or': 1.2},
    'rs2032582':  {'gene': 'ABCB1',  'maf': 0.42,  'effect': 'risk',       'or': 1.2},
    'rs1045642':  {'gene': 'ABCB1',  'maf': 0.48,  'effect': 'protective', 'or': 0.9},
}

def simulate_genotype(maf, n):
    p, q = 1 - maf, maf
    return np.random.choice([0, 1, 2], size=n, p=[p**2, 2*p*q, q**2])

geno = {rsid: simulate_genotype(v['maf'], N) for rsid, v in variants.items()}
df = pd.DataFrame(geno)

df['age'] = np.random.normal(45, 12, N).clip(18, 80).astype(int)
df['cumulative_dose'] = np.random.normal(300, 80, N).clip(100, 550)
df['baseline_lvef'] = np.random.normal(62, 5, N).clip(45, 75)

log_odds = np.zeros(N)
for rsid, v in variants.items():
    beta = np.log(v['or'])
    if v['effect'] == 'protective':
        beta = -beta
    log_odds += df[rsid] * beta

log_odds += (df['cumulative_dose'] - 300) / 80 * 0.6
log_odds += (df['age'] - 45) / 12 * 0.4
log_odds += (df['baseline_lvef'] - 62) / 5 * (-0.3)

intercept = -1.8
prob = 1 / (1 + np.exp(-(log_odds + intercept)))
df['cardiotoxicity'] = (np.random.random(N) < prob).astype(int)

cases = df['cardiotoxicity'].sum()
print(f"Patients: {N}")
print(f"Cardiotoxicity cases: {cases} ({cases/N:.1%})")

df.to_csv(r'D:\College\Project\Project 2\data\synthetic_patients_v2.csv', index=False)

features = list(variants.keys()) + ['age', 'cumulative_dose', 'baseline_lvef']
X = df[features]
y = df['cardiotoxicity']

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

rf = RandomForestClassifier(
    n_estimators=300,
    max_depth=5,
    min_samples_leaf=5,
    class_weight='balanced',
    random_state=42
)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
auc_scores = cross_val_score(rf, X, y, cv=cv, scoring='roc_auc')
print(f"5-Fold CV ROC-AUC: {auc_scores.mean():.3f} +/- {auc_scores.std():.3f}")

rf.fit(X_train, y_train)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

importances = pd.Series(rf.feature_importances_, index=features).sort_values(ascending=True)
colors = []
for f in importances.index:
    if f not in variants:
        colors.append('gray')
    elif variants[f]['effect'] == 'protective':
        colors.append('steelblue')
    else:
        colors.append('coral')

gene_labels = [f"{f}\n({variants[f]['gene']})" if f in variants else f for f in importances.index]
importances.plot(kind='barh', ax=axes[0], color=colors)
axes[0].set_yticklabels(gene_labels, fontsize=8)
axes[0].set_title('Feature Importance\n(coral=risk, blue=protective, gray=clinical)', fontweight='bold')
axes[0].set_xlabel('Importance')

RocCurveDisplay.from_estimator(rf, X_test, y_test, ax=axes[1])
axes[1].set_title('ROC Curve (held-out test set)', fontweight='bold')
axes[1].plot([0,1],[0,1],'k--', alpha=0.5)

plt.tight_layout()
plt.savefig(r'D:\College\Project\Project 2\results\results_classifier_v2.png', dpi=150)
plt.show()

explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_test)
sv = shap_values[1] if isinstance(shap_values, list) else shap_values[:, :, 1]

plt.figure(figsize=(10, 6))
shap.summary_plot(sv, X_test, feature_names=features, show=False)
plt.title('SHAP Summary — Cardiotoxicity Risk Drivers', fontweight='bold')
plt.tight_layout()
plt.savefig(r'D:\College\Project\Project 2\results\shap_summary_v2.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nTop 5 genetic features by importance:")
genetic_imp = pd.Series(rf.feature_importances_, index=features)
genetic_imp = genetic_imp[list(variants.keys())].sort_values(ascending=False)
for rsid, imp in genetic_imp.head(5).items():
    print(f"  {rsid} ({variants[rsid]['gene']}): {imp:.4f}")