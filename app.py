import streamlit as st
import numpy as np
from chempy import balance_stoichiometry
import math
import matplotlib.pyplot as plt
import japanize_matplotlib
import pandas as pd

# ページ基本設定
st.set_page_config(page_title="化学平衡・熱力学計算機", layout="centered")

# カスタムCSSでデザインをリッチに
st.markdown("""
<style>
    .main {
        background-color: #f0f2f6;
    }
    .stButton>button {
        width: 100%;
        border-radius: 8px;
        height: 3.5em;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        font-weight: bold;
        border: none;
        transition: 0.3s;
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(0,0,0,0.2);
    }
    .result-card {
        padding: 25px;
        border-radius: 15px;
        background-color: white;
        box-shadow: 0 10px 25px rgba(0,0,0,0.05);
        margin-top: 20px;
        margin-bottom: 20px;
        border-left: 5px solid #764ba2;
    }
    .step-header {
        font-size: 1.5rem;
        font-weight: bold;
        color: #1e3a8a;
        margin-top: 30px;
        border-bottom: 2px solid #e2e8f0;
        padding-bottom: 10px;
    }
</style>
""", unsafe_allow_html=True)

st.title("🧪 化学平衡・熱力学計算機")
st.write("「反応式の係数合わせ」から「エネルギー計算による自発性判定」までを一括で行います。")

# --- STEP 1: 反応式のバランス ---
st.markdown('<div class="step-header">STEP 1: 化学反応式のバランス</div>', unsafe_allow_html=True)

reaction_input = st.text_input(
    "反応式（またはイオン反応式）を入力してください", 
    placeholder="例: Fe+3 + e- -> Fe+2  または  MnO4- + 8H+ + 5e- -> Mn+2 + 4H2O",
    value="MnO4- + 8H+ + 5e- -> Mn+2 + 4H2O"
)
reaction_input = reaction_input.replace("=", "->").replace("→", "->")
# 全体からスペースを除くのは、H+ とかを壊すのでやめる
# reaction_input = reaction_input.replace(" ", "")

# セッション状態の初期化
if 'balanced_result' not in st.session_state:
    st.session_state.balanced_result = None
if 'reac_dict' not in st.session_state:
    st.session_state.reac_dict = None
if 'prod_dict' not in st.session_state:
    st.session_state.prod_dict = None
if 'calc_results' not in st.session_state:
    st.session_state.calc_results = None

if st.button("1. 反応式をバランスする ⚖️"):
    try:
        import re
        reactants_raw, products_raw = reaction_input.split("->")
        
        # 物質を分離する関数
        def clean_substances(side_str):
            substances = []
            # 1. 物質の間の '+' をセパレータとして扱う。
            # イオンの '+' (H+, Fe+2) と区別するため、「前後にスペースがある + 」または
            # 「前後に他の物質があることが明らかな位置」で分割を試みる
            # ここではシンプルに、前後にスペースを許容する正規表現で分割
            parts = re.split(r'\s+\+\s+|\s*\+\s+(?=[A-Z0-9])|(?<=[A-Za-z0-9+-])\s+\+\s*', side_str.strip())
            
            # もし上記で上手くいかない場合のフォールバック: 
            # 物質名の標準的なパターン（英数字と末尾の電荷）で抽出を試みるのも手だが、
            # まずは各パーツをきれいに掃除する
            for s in parts:
                s = s.strip()
                if not s: continue
                # 先頭の数字（係数）を分離して化学式部分だけを取り出す
                # 例: "2H2O" -> "H2O", "5 Fe+2" -> "Fe+2"
                match = re.match(r'^(\d+)?\s*(.*)$', s)
                formula = match.group(2).strip() if match else s
                if formula:
                    substances.append(formula)
            return substances

        reactants_list = clean_substances(reactants_raw)
        products_list = clean_substances(products_raw)
        
        # 共通する物質がないかチェック
        common = set(reactants_list) & set(products_list)
        if common:
            st.error(f"エラー: 以下の物質が両辺に含まれています: {', '.join(common)}")
            st.warning("【解決策】\n"
                       "- 反応の前後で変化していない物質は除外してください。\n"
                       "- イオン反応式の場合は、電荷（Fe+2 と Fe+3 など）が正しく入力されているか確認してください。")
            st.session_state.balanced_result = None
            st.stop()

        reac, prod = balance_stoichiometry(set(reactants_list), set(products_list))
        
        # 自由度の確認（アンダーデターミンドな反応の場合、x1 などの変数が入る）
        # これを単純な整数解に変換する
        def resolve_symbols(side_dict):
            resolved = {}
            for k, v in side_dict.items():
                if hasattr(v, 'subs'): # sympyのシンボルの場合
                    # 全ての変数に1を代入して具体的な数値にする
                    val = v
                    for symbol in v.free_symbols:
                        val = val.subs(symbol, 1)
                    resolved[k] = int(val)
                else:
                    resolved[k] = int(v)
            return resolved

        reac_resolved = resolve_symbols(reac)
        prod_resolved = resolve_symbols(prod)
        
        st.session_state.reac_dict = reac_resolved
        st.session_state.prod_dict = prod_resolved
        
        def format_side(side):
            # 係数が1の場合は表示しない
            return " + ".join(f"{v if v!=1 else ''}{k}" for k,v in side.items())
        
        st.session_state.balanced_result = f"{format_side(reac_resolved)} -> {format_side(prod_resolved)}"
        
        # 変数（x1等）が含まれていた場合の警告
        has_symbols = any(hasattr(v, 'subs') for v in list(reac.values()) + list(prod.values()))
        if has_symbols:
            st.warning("⚠️ この反応は複数の独立した反応の組み合わせであるため、解が複数存在します。代表的な整数比を表示しています。")
        
        st.success("バランス成功！")
    except Exception as e:
        st.error(f"バランスエラー: {e}")
        st.session_state.balanced_result = None

if st.session_state.balanced_result:
    st.markdown(f"""
    <div class="result-card">
        <small style="color: grey;">バランス済み反応式:</small>
        <h2 style="color: #1e3a8a; margin-top: 5px;">{st.session_state.balanced_result}</h2>
    </div>
    """, unsafe_allow_html=True)

    reac_resolved = st.session_state.reac_dict
    prod_resolved = st.session_state.prod_dict

    # --- STEP 2: 熱力学計算 ---
    st.markdown('<div class="step-header">STEP 2: 自由エネルギー計算 & 自発性判定</div>', unsafe_allow_html=True)
    st.info("この反応全体の標準エンタルピー変化(ΔH)とエントロピー変化(ΔS)を入力してください。")

    col1, col2 = st.columns(2)
    with col1:
        temp_c = st.number_input("温度 (°C)", value=25.0, step=1.0, format="%.2f")
        dH_kj = st.number_input("反応エンタルピー ΔH (kJ/mol)", value=-285.8, step=10.0, format="%.2f")
    with col2:
        temp_k = temp_c + 273.15
        st.write(f"絶対温度: **{temp_k:.2f} K**")
        dS_j = st.number_input("反応エントロピー ΔS (J/mol·K)", value=-163.2, step=1.0, format="%.2f")

    if st.button("2. 自由エネルギー(ΔG)を計算して判定 🌡️"):
        # ΔG = ΔH - TΔS
        delta_g_j = (dH_kj * 1000) - (temp_k * dS_j)
        delta_g_kj = delta_g_j / 1000
        R = 8.31446
        
        try:
            K = math.exp(-delta_g_j / (R * temp_k))
        except (OverflowError, ZeroDivisionError):
            K = float('inf')
            
        st.session_state.calc_results = {
            'delta_g_kj': delta_g_kj,
            'delta_g_j': delta_g_j,
            'K': K
        }

    if st.session_state.calc_results:
        res = st.session_state.calc_results
        st.markdown('<div class="result-card">', unsafe_allow_html=True)
        st.subheader("📊 計算結果")
        
        res_col1, res_col2 = st.columns(2)
        res_col1.metric("ギブス自由エネルギー ΔG", f"{res['delta_g_kj']:.2f} kJ/mol")
        
        if res['K'] == float('inf'):
            res_col2.metric("平衡定数 K", "極めて大きな値")
        else:
            res_col2.metric("平衡定数 K", f"{res['K']:.2e}")

        st.markdown("---")
        st.write("**判定結果:**")
        
        if res['delta_g_kj'] < -0.01:
            st.success("✅ **自発的に進行します (Spontaneous)**\n\n標準状態で反応は右方向（生成系）へ進みます。")
        elif res['delta_g_kj'] > 0.01:
            st.error("❌ **非自発的です (Non-spontaneous)**\n\n標準状態で反応は逆方向へ進みやすい、またはエネルギーの供給が必要です。")
        else:
            st.warning("⚖️ **平衡状態 (Equilibrium)**\n\n反応は平衡状態付近にあります。")
            
        st.markdown('</div>', unsafe_allow_html=True)

    # --- STEP 3: ダイアグラム表示 ---
    # 電子(e-)が含まれているかチェック
    def count_electron(side_dict):
        # キーの中に "e-" または "e" が含まれるか確認
        for k, v in side_dict.items():
            if k.lower() in ["e-", "e", "electron"]:
                return v
        return 0

    n_electron_react = count_electron(reac_resolved)
    n_electron_prod = count_electron(prod_resolved)
    n_electron = abs(n_electron_react - n_electron_prod)

    # H+ の数を取得
    def count_h_plus(side_dict):
        for k, v in side_dict.items():
            if k == "H+": return v
        return 0
    
    m_h_plus = count_h_plus(reac_resolved) - count_h_plus(prod_resolved)

    if n_electron > 0:
        st.markdown('<div class="step-header">STEP 3: E-pH (プルベー) ダイアグラム</div>', unsafe_allow_html=True)
        st.info(f"電子移動 ($n={n_electron}$) を検出しました。電気化学反応として解析します。")
        
        e0 = st.number_input("標準電位 E0 [V vs SHE]", value=1.51, step=0.01, format="%.3f")
        
        # pH範囲
        ph_range = np.linspace(0, 14, 100)
        # E = E0 - (0.0592 * m / n) * pH  (簡易ネルンスト式 at 25°C)
        # 反応式: Ox + mH+ + ne- -> Red  の場合
        e_ph = e0 - (0.0592 * m_h_plus / n_electron) * ph_range

        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 反応ラインの下側にテキストを表示（簡易的な領域判定）
        # 反応 Ox + ne- -> Red において、ラインの上はOx(イオンや酸化物)、下はRed(金属など)
        mid_ph = 7
        mid_e = e0 - (0.0592 * m_h_plus / n_electron) * mid_ph
        
        # 注釈の追加
        ax.text(mid_ph, mid_e + 0.3, "イオンまたは酸化物\n(腐食 / 不動態)", 
                horizontalalignment='center', color='#764ba2', fontsize=10, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        ax.text(mid_ph, mid_e - 0.3, "金属安定領域\n(Immunity)", 
                horizontalalignment='center', color='#1e3a8a', fontsize=10, fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

        # 反応ライン
        ax.plot(ph_range, e_ph, label='反応平衡ライン', color='#764ba2', linewidth=3)
        
        # 水の安定領域 (a: 酸素発生, b: 水素発生)
        e_a = 1.229 - 0.0592 * ph_range
        e_b = 0 - 0.0592 * ph_range
        ax.plot(ph_range, e_a, 'k--', alpha=0.3, label='H2O/O2 (a)')
        ax.plot(ph_range, e_b, 'k--', alpha=0.3, label='H+/H2 (b)')
        ax.fill_between(ph_range, e_b, e_a, color='gray', alpha=0.05, label='水の安定域')

        ax.set_title(f"E-pHダイアグラム: {st.session_state.balanced_result}", fontsize=12)
        ax.set_xlabel("pH", fontsize=10)
        ax.set_ylabel("電位 E [V vs SHE]", fontsize=10)
        ax.set_xlim(0, 14)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='upper right')
        
        st.pyplot(fig)
        st.caption(f"※ ネルンストの式 $E = E^0 - \\frac{{0.0592 m}}{{n}} pH$ に基づく計算。($n={n_electron}, m={m_h_plus}$)")

        # 領域の説明を表示
        st.markdown("""
        ### 🛡️ 領域の解説
        E-pHダイアグラム上の領域は、金属の腐食挙動によって以下のように分類されます：
        
        1. **金属安定領域 (Immunity)**
           - Eが低く、金属が酸化されずに固体のままで存在する領域。**「腐食しない」**状態です。
        2. **酸化物・水酸化物被膜領域 (Passivation / 不動態化)**
           - 金属表面に安定な膜が形成され、内部を保護する領域。**「腐食が抑制される」**状態です。
        3. **溶解領域 (Corrosion / 腐食)**
           - 金属イオンとして溶液中に溶け出す領域。**「腐食が進行する」**状態です。
        """)

    else:
        # 電子が含まれない場合は温度vsΔGを表示
        st.markdown('<div class="step-header">STEP 3: 温度 vs ΔG ダイアグラム</div>', unsafe_allow_html=True)
        # (既存の温度ダイアグラムコード)
        temp_range_c = np.linspace(-100, 1000, 100)
        temp_range_k = temp_range_c + 273.15
        delta_g_range_kj = (dH_kj * 1000 - temp_range_k * dS_j) / 1000

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(temp_range_c, delta_g_range_kj, label='ΔG', color='#764ba2', linewidth=2)
        ax.axhline(0, color='black', linestyle='--', linewidth=1)
        ax.set_title(f"反応: {st.session_state.balanced_result}")
        ax.set_xlabel("温度 [°C]")
        ax.set_ylabel("ΔG [kJ/mol]")
        ax.grid(True, linestyle=':', alpha=0.6)
        st.pyplot(fig)

    # ヘルプ・計算式
    with st.expander("計算式の詳細を表示"):
        st.latex(r"\Delta G = \Delta H - T \Delta S")
        st.latex(r"K = e^{-\frac{\Delta G}{RT}}")
        st.write(r"※ $R = 8.314 \, \mathrm{J/(mol \cdot K)}$ (気体定数)")
else:
    st.info("👆 まずは反応式をバランスしてください。その後に熱力学計算が可能になります。")
