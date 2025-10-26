import streamlit as st
import re
import time
import pandas as pd
from io import StringIO

# ============================
# ✅ PAGE CONFIG
# ============================
st.set_page_config(
    page_title="DNA Pattern Matching Tool",
    page_icon="🧬",
    layout="wide"
)

# ============================
# 🎨 DARK THEME STYLE
# ============================
st.markdown("""
    <style>
        body { background-color: #0E1117; color: #FAFAFA; }
        .stApp { background-color: #0E1117; color: #FAFAFA; }
        h1, h2, h3, h4, h5 { color: #00B4D8 !important; }
        .stTextInput > div > div > input, .stTextArea textarea {
            background-color: #1E2636;
            color: white;
            border: 1px solid #00B4D8;
        }
        .stButton > button {
            background-color: #00B4D8;
            color: white;
            border-radius: 8px;
            padding: 0.6em 1.2em;
            font-weight: bold;
            border: none;
        }
        .stButton > button:hover { background-color: #0077B6; color: white; }
        .result-box {
            background-color: #1E2636;
            padding: 15px;
            border-radius: 10px;
            border: 1px solid #00B4D8;
        }
        .fasta-summary {
            background-color: #1E2636;
            border: 1px solid #00B4D8;
            border-radius: 10px;
            padding: 15px;
            margin-top: 10px;
        }
    </style>
""", unsafe_allow_html=True)

# ============================
# 🧬 HEADER
# ============================
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Tool (FASTA + KMP + Boyer–Moore)</h1>", unsafe_allow_html=True)
st.markdown("<br>", unsafe_allow_html=True)

# ============================
# 📂 FILE UPLOAD
# ============================
uploaded_file = st.file_uploader("📁 Upload a FASTA file", type=["fasta", "fa", "txt"])

fasta_dict = {}
if uploaded_file is not None:
    fasta_content = uploaded_file.getvalue().decode("utf-8", errors="ignore")

    # Split multiple FASTA entries
    entries = fasta_content.strip().split(">")
    for entry in entries:
        if not entry.strip():
            continue
        lines = entry.strip().split("\n")
        header = lines[0].strip()
        seq = "".join(lines[1:]).upper()
        seq = re.sub(r'[^ATCGN]', '', seq)
        fasta_dict[header] = seq

    st.success(f"✅ Loaded {len(fasta_dict)} FASTA sequence(s) successfully!")

    # ============================
    # 📊 FASTA SUMMARY TABLE
    # ============================
    summary_data = []
    for header, seq in fasta_dict.items():
        gc_content = (seq.count("G") + seq.count("C")) / len(seq) * 100 if len(seq) > 0 else 0
        n_count = seq.count("N")
        summary_data.append({
            "Header": header[:80] + ("..." if len(header) > 80 else ""),
            "Length (bp)": len(seq),
            "GC%": round(gc_content, 2),
            "N (Ambiguous)": n_count
        })

    df_summary = pd.DataFrame(summary_data)
    st.markdown("<h4>📊 FASTA Summary</h4>", unsafe_allow_html=True)
    st.dataframe(df_summary, use_container_width=True)

    # Allow user to select which sequence to analyze
    selected_header = st.selectbox("🔹 Select a sequence to analyze", list(fasta_dict.keys()))
    dna_sequence = fasta_dict[selected_header]
else:
    dna_sequence = ""

# ============================
# 🧬 INPUTS
# ============================
col1, col2 = st.columns(2)

with col1:
    if not dna_sequence:
        dna_sequence = st.text_area(
            "Enter DNA Sequence (if no file uploaded):",
            placeholder="Example: ATCGGATCGATCGTACGATCGA...",
            height=150
        ).strip().upper()

with col2:
    pattern = st.text_input(
        "Enter Pattern to Search:",
        placeholder="Example: CGATCGA"
    ).strip().upper()

algorithm = st.radio(
    "Select Algorithm:",
    ["KMP (Knuth–Morris–Pratt)", "Boyer–Moore"],
    horizontal=True
)

# ============================
# 🔍 SEARCH ALGORITHMS
# ============================
def kmp_search(text, pattern):
    lps = [0] * len(pattern)
    j = 0
    for i in range(1, len(pattern)):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j - 1]
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j

    positions = []
    j = 0
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            j = lps[j - 1]
        if text[i] == pattern[j]:
            j += 1
        if j == len(pattern):
            positions.append(i - j + 1)
            j = lps[j - 1]
    return positions


def boyer_moore_search(text, pattern):
    m, n = len(pattern), len(text)
    bad_char = {pattern[i]: i for i in range(m)}
    positions = []
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            positions.append(s)
            s += m - bad_char.get(text[s + m], -1) if s + m < n else 1
        else:
            s += max(1, j - bad_char.get(text[s + j], -1))
    return positions

# ============================
# 🚀 RUN SEARCH
# ============================
if st.button("🔍 Find Pattern"):
    if not dna_sequence or not pattern:
        st.warning("⚠️ Please provide both DNA sequence and pattern.")
    elif len(pattern) > len(dna_sequence):
        st.error("❌ Pattern cannot be longer than the sequence.")
    else:
        with st.spinner("Running pattern matching..."):
            progress = st.progress(0)
            time.sleep(0.2)

            if algorithm == "KMP (Knuth–Morris–Pratt)":
                positions = kmp_search(dna_sequence, pattern)
            else:
                positions = boyer_moore_search(dna_sequence, pattern)

            progress.progress(100)
            time.sleep(0.3)

        # ============================
        # 📊 RESULTS
        # ============================
        if positions:
            total_matches = len(positions)
            percentage = (total_matches * len(pattern) / len(dna_sequence)) * 100
            gc_content = (dna_sequence.count('G') + dna_sequence.count('C')) / len(dna_sequence) * 100

            max_display = 50
            display_positions = positions[:max_display]
            if len(positions) > max_display:
                display_positions.append("... (more)")

            st.markdown(f"""
            <div class='result-box'>
                <h4>✅ Pattern Found!</h4>
                <b>Sequence:</b> {selected_header if uploaded_file else "User Input"}<br>
                <b>Matches:</b> {total_matches}<br>
                <b>Sequence Length:</b> {len(dna_sequence):,} bases<br>
                <b>Pattern Length:</b> {len(pattern)} bases<br>
                <b>Match Coverage:</b> {percentage:.3f}% of total bases<br>
                <b>GC Content:</b> {gc_content:.2f}%<br>
                <b>Positions:</b> {display_positions}
            </div>
            """, unsafe_allow_html=True)

            df = pd.DataFrame({'Match_Position': positions})
            csv = df.to_csv(index=False)
            st.download_button("📥 Download Match Positions (CSV)", csv, "match_positions.csv", "text/csv")
        else:
            st.error("❌ Pattern not found in the selected sequence.")

# ============================
# ℹ️ FOOTER
# ============================
st.markdown("""
---
<p style='text-align:center; color:gray; font-size:0.9em;'>
Developed by <b>Anas Jamshed</b> 🧠 | Advanced FASTA Pattern Matching Tool | Powered by Streamlit
</p>
""", unsafe_allow_html=True)
