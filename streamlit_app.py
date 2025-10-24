import streamlit as st
import time

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
    </style>
""", unsafe_allow_html=True)

# ============================
# 🧬 HEADER
# ============================
st.markdown("""
    <h1 style="text-align:center;">🧬 DNA Pattern Matching Tool (FASTA + KMP + Boyer–Moore)</h1>
""", unsafe_allow_html=True)
st.markdown("<br>", unsafe_allow_html=True)

# ============================
# 📂 FILE UPLOAD
# ============================
uploaded_file = st.file_uploader("📁 Upload a FASTA file", type=["fasta", "fa", "txt"])

dna_sequence = ""
if uploaded_file is not None:
    fasta_content = uploaded_file.read().decode("utf-8")
    # Remove headers starting with '>'
    dna_sequence = "".join(line.strip() for line in fasta_content.splitlines() if not line.startswith(">"))
    st.success(f"✅ File loaded successfully! Sequence length: {len(dna_sequence):,} bases")

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
        )

with col2:
    pattern = st.text_input(
        "Enter Pattern to Search:",
        placeholder="Example: CGATCGA"
    )

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
    if not dna_sequence.strip() or not pattern.strip():
        st.warning("⚠️ Please provide a DNA sequence and a pattern.")
    else:
        with st.spinner("Running pattern matching..."):
            time.sleep(0.5)
            if algorithm == "KMP (Knuth–Morris–Pratt)":
                positions = kmp_search(dna_sequence, pattern)
            else:
                positions = boyer_moore_search(dna_sequence, pattern)

        if positions:
            total_matches = len(positions)
            percentage = (total_matches * len(pattern) / len(dna_sequence)) * 100

            st.markdown(f"""
            <div class='result-box'>
                <h4>✅ Pattern Found!</h4>
                <b>Matches:</b> {total_matches}<br>
                <b>Sequence Length:</b> {len(dna_sequence):,} bases<br>
                <b>Pattern Length:</b> {len(pattern)} bases<br>
                <b>Match Coverage:</b> {percentage:.3f}% of total bases<br>
                <b>Positions:</b> {positions}
            </div>
            """, unsafe_allow_html=True)
        else:
            st.error("❌ Pattern not found in the DNA sequence.")

# ============================
# ℹ️ FOOTER
# ============================
st.markdown("""
---
<p style='text-align:center; color:gray; font-size:0.9em;'>
Developed by <b>Anas Jamshed</b> 🧠 | Powered by Streamlit
</p>
""", unsafe_allow_html=True)
