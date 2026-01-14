import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# ===============================
# Streamlit 設定
# ===============================
st.set_page_config(page_title="Pourbaix Diagram Fe–Cu–Zn–H2O", layout="wide")

st.markdown("""
<div style="background: linear-gradient(135deg,#020617,#1e40af);
            color:white; padding:1.4rem; border-radius:14px; text-align:center;">
<h1>Fe / Cu / Zn – H₂O Pourbaix Diagram</h1>
<p>Thermodynamic stability regions (Ψ minimization)</p>
</div>
""", unsafe_allow_html=True)

# ===============================
# サイドバー
# ===============================
with st.sidebar:
    st.header("System")
    system = st.radio("Metal system", ["Fe–H2O", "Cu–H2O", "Zn–H2O"])
    phase_type = st.radio("Select phase type", ["Oxides only", "Hydroxides only"])

    st.header("Conditions")
    temp_c = st.slider("Temperature [°C]", 0, 100, 25)
    show_boundary = st.checkbox("Show phase boundaries", True)
    show_precip = st.checkbox("Show precipitation", True)

    st.header("Activities (log10)")
    log_a_M1 = st.number_input("log a(M⁺ or M²⁺)", value=-6.0, format="%.1f")
    log_a_M2 = None
    if system == "Fe–H2O":
        log_a_M2 = st.number_input("log a(M³⁺)", value=-6.0, format="%.1f")

# ===============================
# 定数
# ===============================
F = 96485.3
R = 8.31446
T = 273.15 + temp_c
S = R * T * np.log(10) / F
G_H2O = -237130  # J/mol

act_M1 = log_a_M1 * S
act_M2 = log_a_M2 * S if log_a_M2 is not None else 0.0

# ===============================
# 熱力学データ
# ===============================
if system == "Fe–H2O":
    Gf = {
        "Fe": 0.0,
        "Fe2+": -78900,
        "Fe3+": -4700,
        "Fe(OH)2": -486500,
        "Fe(OH)3": -696500,
        "Fe3O4": -1015400,
        "Fe2O3": -742200,
        "HFeO2-": -379000
    }

elif system == "Cu–H2O":
    Gf = {
        "Cu": 0.0,
        "Cu+": 50300,
        "Cu2+": 65600,
        "Cu2O": -17000,
        "CuO": -129700,
        "Cu(OH)2": -365000,
        "HCuO2-": -300000
    }

elif system == "Zn–H2O":
    Gf = {
        "Zn": 0.0,
        "Zn2+": -147800,
        "ZnO": -318300,
        "Zn(OH)2": -611000,
        "HZnO2-": -500000
    }

# ===============================
# メッシュ
# ===============================
res = 600
ph = np.linspace(0, 14, res)
E = np.linspace(-2.5, 2.5, res)
PH, EE = np.meshgrid(ph, E)

# ===============================
# Ψ 計算
# ===============================
Psi = {}

# -------------------------------
# Fe – H2O
# -------------------------------
if system == "Fe–H2O":
    Psi["Fe"] = np.full_like(PH, Gf["Fe"] / F)
    Psi["Fe2+"] = Gf["Fe2+"] / F + act_M1 - 2 * EE
    Psi["Fe3+"] = Gf["Fe3+"] / F + act_M2 - 3 * EE
    Psi["HFeO2-"] = (Gf["HFeO2-"] - 2 * G_H2O)/F - 2*EE - 3*S*PH

    if phase_type == "Oxides only":
        Psi["Fe3O4"] = ((Gf["Fe3O4"] - 4*G_H2O)/F - 8*EE - 8*S*PH)/3
        Psi["Fe2O3"] = ((Gf["Fe2O3"] - 3*G_H2O)/F - 6*EE - 6*S*PH)/2
    elif phase_type == "Hydroxides only":
        Psi["Fe(OH)2"] = (Gf["Fe(OH)2"] - 2*G_H2O)/F - 2*EE - 2*S*PH
        Psi["Fe(OH)3"] = (Gf["Fe(OH)3"] - 3*G_H2O)/F - 3*EE - 3*S*PH

# -------------------------------
# Cu – H2O
# -------------------------------
elif system == "Cu–H2O":
    Psi["Cu"] = np.full_like(PH, Gf["Cu"] / F)
    Psi["Cu+"] = Gf["Cu+"] / F + act_M1 - EE
    Psi["Cu2+"] = Gf["Cu2+"] / F + act_M2 - 2*EE
    Psi["HCuO2-"] = (Gf["HCuO2-"] - 2*G_H2O)/F - 2*EE - 3*S*PH

    if phase_type == "Oxides only":
        Psi["Cu2O"] = (Gf["Cu2O"] - G_H2O)/F - 2*EE - 2*S*PH
        Psi["CuO"]  = (Gf["CuO"] - G_H2O)/F - 2*EE - 2*S*PH

# -------------------------------
# Zn – H2O
# -------------------------------
elif system == "Zn–H2O":
    Psi["Zn"] = np.full_like(PH, Gf["Zn"] / F)
    Psi["Zn2+"] = Gf["Zn2+"]/F + act_M1 - 2*EE
    Psi["HZnO2-"] = (Gf["HZnO2-"] - 2*G_H2O)/F - 2*EE - 3*S*PH

    if phase_type == "Oxides only":
        Psi["ZnO"] = (Gf["ZnO"] - G_H2O)/F - 2*EE - 2*S*PH

# ===============================
# 相決定
# ===============================
labels = list(Psi.keys())
Psi_stack = np.stack([Psi[k] for k in labels], axis=0)
phase_index = np.argmin(Psi_stack, axis=0)
dominant_phase = np.array(labels)[phase_index]

# ===============================
# 水酸化物沈殿判定（支配相に関係なく全領域）
# ===============================
if show_precip:
    if system == "Zn–H2O" and phase_type == "Hydroxides only":
        logKsp_ZnOH2 = -16.4
        log_a_OH = PH - 14
        log_a_Zn2_sat = logKsp_ZnOH2 - 2*log_a_OH
        precip_ZnOH2 = act_M1 > log_a_Zn2_sat

    if system == "Cu–H2O" and phase_type == "Hydroxides only":
        logKsp_CuOH2 = -19.0
        log_a_OH = PH - 14
        log_a_Cu2_sat = logKsp_CuOH2 - 2*log_a_OH
        precip_CuOH2 = act_M2 > log_a_Cu2_sat

# ===============================
# 描画
# ===============================
colors = ["#94a3b8", "#3b82f6", "#facc15", "#22c55e",
          "#f87171", "#a855f7", "#fb923c", "#0ea5e9"]

fig, ax = plt.subplots(figsize=(10, 8), dpi=120)
ax.imshow(phase_index, origin="lower", extent=[0,14,-2.5,2.5],
          aspect="auto", cmap=ListedColormap(colors[:len(labels)]))

# 水分解線
ax.plot(ph, 1.229 - S*ph, "k--", alpha=0.4)
ax.plot(ph, 0.000 - S*ph, "k--", alpha=0.4)

# 水酸化物沈殿描画
if show_precip:
    if system == "Zn–H2O" and phase_type == "Hydroxides only":
        ax.contourf(PH, EE, precip_ZnOH2, colors=['#22c55e'], alpha=0.2)
    if system == "Cu–H2O" and phase_type == "Hydroxides only":
        ax.contourf(PH, EE, precip_CuOH2, colors=['#f87171'], alpha=0.2)

# 境界線
if show_boundary:
    for i in range(len(labels)):
        for j in range(i+1, len(labels)):
            ax.contour(PH, EE, Psi_stack[i]-Psi_stack[j],
                       levels=[0], colors="white", linewidths=0.6, alpha=0.6)

# ラベル
for i, name in enumerate(labels):
    mask = phase_index == i
    if np.any(mask):
        ax.text(PH[mask].mean(), EE[mask].mean(), name,
                ha="center", va="center", fontsize=10,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

ax.set_xlabel("pH")
ax.set_ylabel("Potential E [V vs SHE]")
ax.set_title(f"{system} Pourbaix Diagram @ {temp_c} °C")
ax.set_xlim(0,14)
ax.set_ylim(-2.5,2.5)
ax.grid(alpha=0.1)

st.pyplot(fig)
