import streamlit as st
import numpy as np
from chempy import balance_stoichiometry
import math
import matplotlib.pyplot as plt
import japanize_matplotlib
import pandas as pd
import re

# Optional: Try to import thermo, but use a fallback internal DB for ions
try:
    from thermo.chemical import Chemical
    HAS_THERMO = True
except ImportError:
    HAS_THERMO = False

# ページ基本設定
st.set_page_config(page_title="🧪 多反応 E-pH & 熱力学アナライザー", layout="wide")

# --- 拡張熱力学データベース (主要なイオン・酸化物) ---
# これらは標準的なライブラリでは取得しにくいため、内蔵データベースで優先処理します
THERMO_DB = {
    # 水系
    "H2O": {"h": -285.8, "s": 69.9},
    "H+1": {"h": 0.0, "s": 0.0},
    "OH-1": {"h": -230.0, "s": -10.8},
    "O2": {"h": 0.0, "s": 205.1},
    "H2": {"h": 0.0, "s": 130.7},
    "e-1": {"h": 0.0, "s": 0.0},
    # 鉄系
    "Fe": {"h": 0.0, "s": 27.3},
    "Fe+2": {"h": -89.1, "s": -137.7},
    "Fe+3": {"h": -48.5, "s": -315.9},
    "Fe(OH)2": {"h": -569.0, "s": 88.0},
    "Fe(OH)3": {"h": -823.0, "s": 106.7},
    "Fe3O4": {"h": -1118.4, "s": 146.4},
    "Fe2O3": {"h": -824.2, "s": 87.4},
    "FeCl2": {"h": -341.8, "s": 117.9},
    "FeCl3": {"h": -399.5, "s": 142.3},
    "Cl-1": {"h": -167.2, "s": 56.5},
    # 銅・亜鉛系
    "Cu": {"h": 0.0, "s": 33.1},
    "Cu+2": {"h": 64.8, "s": -99.6},
    "CuO": {"h": -157.3, "s": 42.6},
    "Zn": {"h": 0.0, "s": 41.6},
    "Zn+2": {"h": -153.3, "s": -112.1},
    "ZnO": {"h": -350.5, "s": 43.7},
}

# カスタムCSS (プレミアムデザイン)
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        padding: 2rem;
        border-radius: 15px;
        color: white;
        text-align: center;
        box-shadow: 0 4px 20px rgba(0,0,0,0.2);
        margin-bottom: 2rem;
    }
    .stButton>button {
        border-radius: 10px;
        font-weight: bold;
        transition: all 0.3s cubic-bezier(0.175, 0.885, 0.32, 1.275);
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }
    .card {
        background: white;
        padding: 1.5rem;
        border-radius: 12px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header"><h1>🧪 Multi-Reaction E-pH & Thermodynamics</h1><p>Multiple registration • LaTeX Support • Overlaid Plots</p></div>', unsafe_allow_html=True)

# セッション状態の初期化
if 'reactions_df' not in st.session_state:
    st.session_state.reactions_df = pd.DataFrame(columns=[
        '元の式', 'バランス済み', 'E0 (V)', 'dH (kJ/mol)', 'dS (J/mol·K)', 'n', 'm'
    ])

def clean_formula(formula):
    """Fe^{3+}, [FeCl3], Cl^- などの多様な表記を正規化"""
    f = formula.strip()
    # 1. LaTeXの上付き・下付きを除去: ^{3+} -> +3, _{2} -> 2
    f = re.sub(r'\^\{?(\d*)([+-])\}?', r'\2\1', f)
    f = re.sub(r'_\{?(\d+)\}?', r'\1', f)
    # 2. 角括弧を除去 (丸括弧は維持)
    f = f.replace('[', '').replace(']', '')
    # 3. 3+ -> +3 形式への統一
    f = re.sub(r'(\d+)([+-])', r'\2\1', f)
    # 4. 電子 e- の正規化
    if f.lower() in ['e', 'e-', 'e-1', 'electron']: return 'e-1'
    # 5. 符号のみの場合 (H+, Cl-) に係数を補完
    if f.endswith('+') and not any(c.isdigit() for c in f.split('+')[-1]): f += '1'
    if f.endswith('-') and not any(c.isdigit() for c in f.split('-')[-1]): f += '1'
    return f

def parse_side(side_str):
    """辺を物質リストに分解。電荷の + と区切り文字の + を判別"""
    # 前後にスペースがある '+' を優先して分割、または大文字・数字・括弧の前の '+' で分割
    parts = re.split(r'\s+\+\s+|\s*\+\s*(?=[A-Z0-9(\[])', side_str.strip())
    subs = []
    for s in parts:
        s = s.strip()
        if not s: continue
        # 先頭の数字（係数）を分離
        m = re.match(r'^(\d+)?\s*(.*)$', s)
        formula = m.group(2).strip() if m else s
        if formula:
            subs.append(clean_formula(formula))
    return subs

def get_thermo_data(species):
    """DB または thermo ライブラリからデータを取得"""
    # 1. 内蔵DB優先 (イオン等に対応)
    if species in THERMO_DB:
        return THERMO_DB[species]['h'], THERMO_DB[species]['s']
    
    # 2. thermo ライブラリ (中性分子に対応)
    if HAS_THERMO:
        try:
            c = Chemical(species)
            return c.Hf298 / 1000.0, c.S0  # kJ/mol, J/mol·K
        except:
            pass
    return 0.0, 0.0

def estimate_reaction_thermo(r_dict, p_dict):
    """反応全体の dH, dS を推定"""
    dh, ds = 0.0, 0.0
    for s, c in r_dict.items():
        h, s_val = get_thermo_data(s)
        dh -= c * h
        ds -= c * s_val
    for s, c in p_dict.items():
        h, s_val = get_thermo_data(s)
        dh += c * h
        ds += c * s_val
    return dh, ds

# --- 反応の登録セクション ---
with st.container():
    st.subheader("➕ 反応の登録")
    col_in, col_btn = st.columns([4, 1])
    with col_in:
        new_input = st.text_input("反応式を入力 (例: Fe^{3+} + 3 H2O <=> Fe(OH)3 + 3 H+)", 
                                  placeholder="MnO4- + 8 H+ + 5 e- -> Mn2+ + 4 H2O",
                                  label_visibility="collapsed")
    with col_btn:
        if st.button("追加 📥", use_container_width=True):
            if new_input:
                try:
                    # --- 1. LaTeX特有の記法をクリーニング ---
                    s = new_input
                    # 数式デリミタの除去
                    s = s.replace(r'\(', '').replace(r'\)', '').replace(r'\[', '').replace(r'\]', '')
                    # \text{...} -> ... への変換
                    s = re.sub(r'\\text\{([^}]*)\}', r'\1', s)
                    # LaTeX矢印の変換
                    s = re.sub(r'\\(?:long)?rightarrow', ' -> ', s)
                    s = re.sub(r'\\(?:long)?Rightarrow', ' -> ', s)
                    s = re.sub(r'\\(?:long)?leftrightarrow', ' <=> ', s)
                    # バックスラッシュの残骸を除去
                    s = s.replace('\\', '')
                    
                    # 記号の正規化（既存のロジック）
                    for arrow in ["<=>", "⇌", "⇄", "⇆", "<->", "-->", "->", "==", "=", "→", "⇒"]:
                        s = s.replace(arrow, " -> ")
                    
                    if " -> " not in s:
                        st.error("区切り記号 (->, <=>, =) が見つかりません。")
                    else:
                        parts = s.split(" -> ")
                        def trim_junk(t): return re.sub(r'^[<>=+]*', '', re.sub(r'[<>=+]*$', '', t.strip()))
                        r_raw, p_raw = trim_junk(parts[0]), trim_junk(parts[1])
                        
                        r_list = parse_side(r_raw)
                        p_list = parse_side(p_raw)
                        
                        # chempy によるバランス
                        reac, prod = balance_stoichiometry(set(r_list), set(p_list))
                        
                        def to_int_d(d):
                            return {k: int(v.subs({sym: 1 for sym in v.free_symbols})) if hasattr(v, 'free_symbols') else int(v) for k, v in d.items()}
                        
                        r_dict = to_int_d(reac)
                        p_dict = to_int_d(prod)
                        
                        def fmt(d): return " + ".join(f"{(str(v) if v!=1 else '')}{k}" for k,v in d.items())
                        balanced = f"{fmt(r_dict)} -> {fmt(p_dict)}"
                        
                        # 特徴抽出 (n: 電子, m: H+)
                        def count_spec(d, t):
                            for k, v in d.items():
                                kl = k.lower()
                                if t == 'e' and kl == 'e-1': return v
                                if t == 'h' and kl in ['h+1', 'h+']: return v
                            return 0
                        
                        n_e = abs(count_spec(r_dict, 'e') - count_spec(p_dict, 'e'))
                        m_h = count_spec(r_dict, 'h') - count_spec(p_dict, 'h')
                        
                        # 熱力学推定
                        dh_est, ds_est = estimate_reaction_thermo(r_dict, p_dict)
                        e0_est = 0.0
                        if n_e > 0:
                            dg_est = (dh_est * 1000) - 298.15 * ds_est
                            e0_est = -dg_est / (n_e * 96485)
                        
                        new_row = {
                            '元の式': new_input,
                            'バランス済み': balanced,
                            'E0 (V)': round(e0_est, 3),
                            'dH (kJ/mol)': round(dh_est, 1),
                            'dS (J/mol·K)': round(ds_est, 1),
                            'n': n_e,
                            'm': m_h
                        }
                        st.session_state.reactions_df = pd.concat([st.session_state.reactions_df, pd.DataFrame([new_row])], ignore_index=True)
                        st.success(f"追加: {balanced}")
                except Exception as e:
                    st.error(f"解析エラー: {e}")

# --- リストとグラフ出力 ---
if not st.session_state.reactions_df.empty:
    st.subheader("📋 反応リストとパラメータの調整")
    st.markdown("※ 数値を書き換えると、下のグラフがリアルタイムに更新されます。")
    edited_df = st.data_editor(
        st.session_state.reactions_df,
        use_container_width=True,
        num_rows="dynamic",
        column_config={
            "バランス済み": st.column_config.Column(width="large", disabled=True),
            "E0 (V)": st.column_config.NumberColumn(format="%.3f", help="電気化学反応の標準電位"),
            "n": st.column_config.Column("e-数", width="small", disabled=True),
            "m": st.column_config.Column("H+数", width="small", disabled=True),
        }
    )
    st.session_state.reactions_df = edited_df

    c1, c2 = st.columns(2)
    with c1:
        st.subheader("⚡ E-pH ダイアグラム (Overlaid)")
        fig, ax = plt.subplots(figsize=(8, 6))
        phs = np.linspace(0, 14, 100)
        # 水の安定域
        ax.plot(phs, 1.229 - 0.0592*phs, 'k--', alpha=0.3, label='O2/H2O')
        ax.plot(phs, 0 - 0.0592*phs, 'k--', alpha=0.3, label='H+/H2')
        ax.fill_between(phs, -0.0592*phs, 1.229 - 0.0592*phs, color='blue', alpha=0.05)
        
        any_e = False
        for i, row in edited_df.iterrows():
            # n, m, E0 が有効な数値であることを確認
            n = row['n']
            m = row['m']
            e0 = row['E0 (V)']
            
            if pd.notna(n) and pd.notna(m) and pd.notna(e0) and n > 0:
                e_line = e0 - (0.0592 * m / n) * phs
                label = row['バランス済み'] if pd.notna(row['バランス済み']) else f"反応 {i+1}"
                ax.plot(phs, e_line, label=f"R{i+1}: {label}", lw=2.5)
                any_e = True
        ax.set_xlabel("pH"); ax.set_ylabel("電位 E [V vs SHE]"); ax.set_xlim(0, 14); ax.set_ylim(-1.5, 2.0)
        ax.grid(True, alpha=0.3); ax.legend(fontsize='x-small', loc='upper right')
        st.pyplot(fig)

    with c2:
        st.subheader("🌡️ 自由エネルギー変化 (ΔG vs Temp)")
        fig, ax = plt.subplots(figsize=(8, 6))
        tc = np.linspace(-50, 500, 100); tk = tc + 273.15
        any_dg = False
        for i, row in edited_df.iterrows():
            n = row['n']
            dh = row['dH (kJ/mol)']
            ds = row['dS (J/mol·K)']
            
            if pd.notna(n) and pd.notna(dh) and pd.notna(ds) and n == 0:
                dg_val = (dh * 1000 - tk * ds) / 1000
                label = row['バランス済み'] if pd.notna(row['バランス済み']) else f"反応 {i+1}"
                ax.plot(tc, dg_val, label=f"R{i+1}: {label}", lw=2)
                any_dg = True
        ax.axhline(0, color='red', alpha=0.2)
        ax.set_xlabel("Temperature [°C]"); ax.set_ylabel("ΔG [kJ/mol]"); ax.grid(True, alpha=0.3)
        if any_dg: ax.legend(fontsize='x-small', loc='upper right')
        st.pyplot(fig)

    if st.button("全データをリセット 🗑️", type="secondary"):
        st.session_state.reactions_df = pd.DataFrame(columns=['元の式', 'バランス済み', 'E0 (V)', 'dH (kJ/mol)', 'dS (J/mol·K)', 'n', 'm'])
        st.rerun()
else:
    st.info("サイドバー（または上部）から反応式を入力して開始してください。複数の反応を登録して同時に比較できます。")
    with st.expander("📚 利用可能な記法のヒント"):
        st.markdown("""
        - **イオン**: `Fe^{3+}`, `Cu2+`, `Cl^-`
        - **錯体・水酸化物**: `[FeCl3]`, `Fe(OH)3`, `Fe3O4`
        - **矢印**: `->`, `<=>`, `⇌`, `＝`
        - **自動データ**: 主要な鉄・銅・水系物質は自動で熱力学データが入力されます。
        """)
