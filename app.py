import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import japanize_matplotlib

# --- Streamlit 設定 ---
st.set_page_config(page_title="🧪 Fe–H2O Pourbaix Diagram (Oxides vs Hydroxides)", layout="wide")
st.markdown("""
<div style="background: linear-gradient(135deg,#0f172a,#1e3a8a);
            color:white; padding:1.5rem; border-radius:12px; text-align:center;">
    <h1>🧪 Fe–H2O Pourbaix Diagram (酸化物 / 水酸化物 選択可)</h1>
</div>
""", unsafe_allow_html=True)

# --- サイドバー ---
with st.sidebar:
    st.header("⚙️ パラメータ")
    temp_c = st.slider("温度 [°C]", 0, 100, 25)
    log_a_fe2 = st.number_input("Fe²⁺ 活量 log₁₀(a)", value=-6.0, format="%.1f")
    log_a_fe3 = st.number_input("Fe³⁺ 活量 log₁₀(a)", value=-6.0, format="%.1f")
    phase_type = st.radio("酸化物 / 水酸化物の選択", ["酸化物のみ", "水酸化物のみ"])
    show_boundary = st.checkbox("デバッグ境界線表示", value=True)

# --- 定数 ---
F = 96485.3
R = 8.31446
T = 273.15 + temp_c
S = R*T*np.log(10)/F
act_fe2 = log_a_fe2 * S
act_fe3 = log_a_fe3 * S
G_H2O = -237130  # J/mol

# --- 標準生成ギブズ [J/mol] ---
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

# --- Meshgrid ---
res = 600
ph_vec = np.linspace(0, 14, res)
e_vec = np.linspace(-2.5, 2.5, res)
PH, E = np.meshgrid(ph_vec, e_vec)

# --- Psi 計算関数 ---
def calc_psi(PH, E, phase_type):
    Psi = {}
    # 金属・イオン共通
    Psi["Fe"] = np.full_like(PH, Gf["Fe"]/F)
    Psi["Fe2+"] = Gf["Fe2+"]/F + act_fe2 - 2*E
    Psi["Fe3+"] = Gf["Fe3+"]/F + act_fe3 - 3*E
    Psi["HFeO2-"] = (Gf["HFeO2-"] - 2*G_H2O)/F - 2*E - 3*S*PH + act_fe2

    if phase_type == "水酸化物のみ":
        Psi["Fe(OH)2"] = (Gf["Fe(OH)2"] - 2*G_H2O)/F - 2*E - 2*S*PH + act_fe2
        Psi["Fe(OH)3"] = (Gf["Fe(OH)3"] - 3*G_H2O)/F - 3*E - 3*S*PH + act_fe3
    else:  # 酸化物のみ
        Psi["Fe3O4"] = ((Gf["Fe3O4"] - 4*G_H2O)/F - 8*E - 8*S*PH)/3
        Psi["Fe2O3"] = ((Gf["Fe2O3"] - 3*G_H2O)/F - 6*E - 6*S*PH)/2

    return Psi

# --- Psi 計算 ---
Psi_dict = calc_psi(PH, E, phase_type)

# --- 使用するフェーズのキーを選択 ---
if phase_type == "水酸化物のみ":
    psi_keys = ["Fe", "Fe2+", "Fe3+", "Fe(OH)2", "Fe(OH)3", "HFeO2-"]
else:
    psi_keys = ["Fe", "Fe2+", "Fe3+", "Fe3O4", "Fe2O3", "HFeO2-"]

Psi_stack = np.stack([Psi_dict[k] for k in psi_keys], axis=0)
phase_map = np.argmin(Psi_stack, axis=0)

# --- 描画 ---
colors = ['#94a3b8','#3b82f6','#facc15','#60a5fa','#f87171','#a855f7','#22c55e','#fb923c']
labels_dict = {
    "Fe": r"$Fe$",
    "Fe2+": r"$Fe^{2+}$",
    "Fe3+": r"$Fe^{3+}$",
    "Fe(OH)2": r"$Fe(OH)_2$",
    "Fe(OH)3": r"$Fe(OH)_3$",
    "Fe3O4": r"$Fe_3O_4$",
    "Fe2O3": r"$Fe_2O_3$",
    "HFeO2-": r"$HFeO_2^-$"
}
labels = [labels_dict[k] for k in psi_keys]

fig, ax = plt.subplots(figsize=(10,8), dpi=120)
ax.imshow(
    phase_map,
    origin='lower',
    cmap=ListedColormap(colors[:len(psi_keys)]),
    extent=[0,14,-2.5,2.5],
    aspect='auto'
)

# 水の分解線
ax.plot(ph_vec, 1.229 - S*ph_vec,'k--', alpha=0.4)
ax.plot(ph_vec, 0.0 - S*ph_vec,'k--', alpha=0.4)

# ラベル
for idx, lab in enumerate(labels):
    mask = (phase_map==idx)
    if np.any(mask):
        ph_c = PH[mask].mean()
        e_c = E[mask].mean()
        ax.text(ph_c, e_c, lab, color='black', fontsize=10, ha='center', va='center',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round'))

# 境界線
if show_boundary:
    line_style = {'colors':'white','linewidths':0.7,'alpha':0.6}
    psi_list = [Psi_dict[k] for k in psi_keys]
    for i in range(len(psi_list)):
        for j in range(i+1,len(psi_list)):
            ax.contour(PH, E, psi_list[i]-psi_list[j], levels=[0], **line_style)

ax.set_xlabel("pH")
ax.set_ylabel("電位 E [V vs SHE]")
ax.set_xlim(0,14)
ax.set_ylim(-2.5,2.5)
ax.grid(alpha=0.1)
ax.set_title(f"Fe–H2O Pourbaix Diagram @ {temp_c}°C, log a(Fe²⁺)={log_a_fe2}, log a(Fe³⁺)={log_a_fe3}")

st.pyplot(fig)
