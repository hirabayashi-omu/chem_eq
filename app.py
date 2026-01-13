import streamlit as st
import numpy as np
from chempy import balance_stoichiometry
import math
import matplotlib.pyplot as plt
import japanize_matplotlib
import pandas as pd

# ページ基本設定
st.set_page_config(page_title="化学平衡・熱力学計算機", layout="centered")

# カスタムCSS
st.markdown("""
<style>
    .stButton>button {
        width: 100%;
        border-radius: 8px;
        height: 3.5em;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        font-weight: bold;
    }
    .result-card {
        padding: 20px;
        border-radius: 12px;
        background-color: white;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        margin: 15px 0;
        border-left: 5px solid #764ba2;
    }
    .step-header {
        font-size: 1.4rem;
        font-weight: bold;
        color: #1e3a8a;
        margin-top: 25px;
        border-bottom: 1px solid #ddd;
        padding-bottom: 8px;
    }
</style>
""", unsafe_allow_html=True)

st.title("🧪 化学平衡・熱力学計算機")

# セッション状態の初期化
for key in ['balanced_result', 'reac_dict', 'prod_dict', 'calc_results']:
    if key not in st.session_state:
        st.session_state[key] = None

# --- STEP 1: 反応式のバランス ---
st.markdown('<div class="step-header">STEP 1: 化学反応式のバランス</div>', unsafe_allow_html=True)

reaction_input = st.text_input(
    "反応式（またはイオン反応式）を入力してください", 
    placeholder="例: MnO4- + 8H+ + 5e- -> Mn+2 + 4H2O",
    value="MnO4- + 8H+ + 5e- -> Mn+2 + 4H2O"
)
reaction_input_clean = reaction_input.replace("=", "->").replace("→", "->")

if st.button("1. 反応式をバランスする ⚖️"):
    try:
        import re
        if "->" not in reaction_input_clean:
            st.error("エラー: '->' で 反応物 -> 生成物 を区切ってください。")
        else:
            parts = reaction_input_clean.split("->")
            react_raw, prod_raw = parts[0], parts[1]
            
            def parse_side(side_str):
                # 物質を '+' で分割（イオンの '+' を避けるための前方参照/後方参照は複雑なので、単純にスペース付き '+' で分割）
                subs = []
                for s in re.split(r'\s+\+\s+|\s*\+\s+(?=[A-Z0-9])', side_str.strip()):
                    s = s.strip()
                    if not s: continue
                    # 先頭の数字（係数）を削除
                    m = re.match(r'^(\d+)?\s*(.*)$', s)
                    formula = m.group(2).strip() if m else s
                    if formula: subs.append(formula)
                return subs

            r_list = parse_side(react_raw)
            p_list = parse_side(prod_raw)
            
            common = set(r_list) & set(p_list)
            if common:
                st.error(f"エラー: 両辺に共通の物質があります: {common}")
            else:
                reac, prod = balance_stoichiometry(set(r_list), set(p_list))
                
                def to_int_dict(d):
                    res = {}
                    for k, v in d.items():
                        if hasattr(v, 'free_symbols'):
                            res[k] = int(v.subs({s: 1 for s in v.free_symbols}))
                        else:
                            res[k] = int(v)
                    return res

                st.session_state.reac_dict = to_int_dict(reac)
                st.session_state.prod_dict = to_int_dict(prod)
                
                # フォーマット
                def fmt(d):
                    return " + ".join(f"{v if v!=1 else ''}{k}" for k,v in d.items())
                
                st.session_state.balanced_result = f"{fmt(st.session_state.reac_dict)} -> {fmt(st.session_state.prod_dict)}"
                st.session_state.calc_results = None # 新しい反応なので計算結果クリア
                st.success("バランスに成功しました！")
    except Exception as e:
        st.error(f"エラー: {e}")

# 結果表示
if st.session_state.balanced_result:
    st.markdown(f'<div class="result-card">バランス済み: <b>{st.session_state.balanced_result}</b></div>', unsafe_allow_html=True)
    
    # ここでローカル変数に代入（session_stateから確実に取得）
    r_dict = st.session_state.reac_dict
    p_dict = st.session_state.prod_dict

    # --- STEP 2: 熱力学 ---
    st.markdown('<div class="step-header">STEP 2: 自由エネルギー計算</div>', unsafe_allow_html=True)
    c1, c2 = st.columns(2)
    with c1:
        tc = st.number_input("温度 (°C)", value=25.0)
        dh = st.number_input("ΔH (kJ/mol)", value=-285.8)
    with c2:
        tk = tc + 273.15
        ds = st.number_input("ΔS (J/mol·K)", value=-163.2)
    
    if st.button("2. エネルギー状態を計算 🌡️"):
        dg_j = (dh * 1000) - (tk * ds)
        dg_kj = dg_j / 1000
        try:
            k_val = math.exp(-dg_j / (8.314 * tk))
        except:
            k_val = float('inf')
        st.session_state.calc_results = {'dg': dg_kj, 'k': k_val, 'tc': tc, 'dh': dh, 'ds': ds}

    if st.session_state.calc_results:
        res = st.session_state.calc_results
        st.write(f"**ΔG:** {res['dg']:.2f} kJ/mol | **K:** {res['k']:.2e}")
        if res['dg'] < 0: st.success("自発的に進行します")
        else: st.warning("非自発的です")

    # --- STEP 3: ダイアグラム ---
    def get_val(d, target):
        for k, v in d.items():
            if k.lower() == target.lower() or k.lower() == target.lower()+'-': return v
        return 0

    n_e = abs(get_val(r_dict, 'e') - get_val(p_dict, 'e'))
    m_h = get_val(r_dict, 'H+') - get_val(p_dict, 'H+')

    if n_e > 0:
        st.markdown('<div class="step-header">STEP 3: E-pH ダイアグラム</div>', unsafe_allow_html=True)
        e0 = st.number_input("標準電位 E0 [V]", value=1.51)
        phs = np.linspace(0, 14, 100)
        es = e0 - (0.0592 * m_h / n_e) * phs
        
        fig, ax = plt.subplots()
        ax.plot(phs, es, label='Equilibrium', color='purple', lw=2)
        # 水の領域
        ax.plot(phs, 1.23 - 0.059 * phs, 'k--', alpha=0.3, label='O2/H2O')
        ax.plot(phs, 0 - 0.059 * phs, 'k--', alpha=0.3, label='H+/H2')
        ax.set_xlabel("pH"); ax.set_ylabel("E [V]"); ax.legend(); ax.grid(True)
        st.pyplot(fig)
    elif st.session_state.calc_results:
        st.markdown('<div class="step-header">STEP 3: 温度依存性</div>', unsafe_allow_html=True)
        res = st.session_state.calc_results
        tr = np.linspace(-100, 500, 100)
        dg_v = (res['dh'] * 1000 - (tr + 273.15) * res['ds']) / 1000
        fig, ax = plt.subplots()
        ax.plot(tr, dg_v, color='orange')
        ax.axhline(0, color='black', lw=1); ax.grid(True)
        ax.set_xlabel("Temp [°C]"); ax.set_ylabel("ΔG [kJ/mol]")
        st.pyplot(fig)
