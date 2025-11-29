# ================================
# DNA Pattern Matching Analyzer🧬 
# ================================
# This application compares different string matching algorithms for DNA sequence analysis
# It allows users to upload FASTA files or input DNA sequences manually and search for patterns

# Importing Libraries
import streamlit as st
import re
import time
import io
import unicodedata
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque, defaultdict
import zipfile
from datetime import datetime

# ReportLab for PDF generation (optional)
try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
        Table, TableStyle
    )
    from reportlab.lib.styles import getSampleStyleSheet
    reportlab_available = True
except:
    reportlab_available = False

# -----------------------------
# Helper: Convert matplotlib figure -> bytes
# -----------------------------
def fig_to_bytes(fig, format="png"):
    buf = io.BytesIO()
    fig.savefig(buf, format=format, dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()

# -----------------------------
# Unicode normalization + name fixes
# -----------------------------
def normalize_text(s):
    if not isinstance(s, str):
        return s
    return unicodedata.normalize("NFC", s)

fix_map = {
    "NaÃ¯ve Search": "Naïve Search",
    "NaÃ¯ve": "Naïve",
    "Boyerâ€“Moore": "Boyer–Moore",
    "Rabinâ€“Karp": "Rabin–Karp",
    "Ahoâ€“Corasick": "Aho–Corasick",
}

def fix_garbled(s):
    if not isinstance(s, str):
        return s
    for bad, good in fix_map.items():
        s = s.replace(bad, good)
    return s

# ==========================
# Page config & styling
# ==========================
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    body, .stApp { background-color: #0E1117; color: #FAFAFA; }
    h1, h2, h3 { color: #00B4D8 !important; }
    .stButton>button {
        background-color: #00B4D8; color: white; font-weight: bold;
        border-radius: 8px; border: none; padding: 0.6em 1.2em;
    }
    .stButton>button:hover { background-color: #0077B6; }
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare multiple string-matching algorithms</p>", unsafe_allow_html=True)
st.markdown("---")

# ==========================
# File upload / manual input
# ==========================
uploaded_files = st.file_uploader("📁 Upload FASTA files", type=["fasta", "fa", "txt"], accept_multiple_files=True)
sequences = {}

if uploaded_files:
    for uploaded_file in uploaded_files:
        fasta = uploaded_file.getvalue().decode("utf-8")
        lines = fasta.strip().split("\n")
        header = lines[0] if lines and lines[0].startswith(">") else uploaded_file.name
        sequence = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
        sequence = re.sub(r"[^ATCG]", "", sequence)
        sequences[header] = sequence
        st.success(f"✅ Loaded: {header} ({len(sequence)} bp)")
else:
    seq_input = st.text_area("🧬 Enter DNA Sequence Manually", height=120).strip().upper()
    if seq_input:
        sequences["Manual Entry"] = re.sub(r"[^ATCG]", "", seq_input)

# ==========================
# Pattern input & algorithm selection
# ==========================
pattern_input = st.text_input("🔍 Enter Pattern(s) (comma separated)", placeholder="ATGC,CTAGA,TTCGA").strip().upper()
algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=algorithms)

if not selected_algos:
    st.info("Please select at least one algorithm (default: all selected).")

# ==========================
# String matching implementations
# ==========================
def naive_search(text, pattern):
    if not pattern: return [], 0
    comparisons = 0
    hits = []
    for i in range(len(text) - len(pattern) + 1):
        comparisons += 1
        match = True
        for j in range(len(pattern)):
            comparisons += 1
            if text[i+j] != pattern[j]:
                match = False
                break
        if match:
            hits.append(i)
    return hits, comparisons

def kmp_search(text, pattern):
    if not pattern: return [], 0
    comparisons = 0
    m = len(pattern)
    lps = [0]*m
    j = 0
    # build lps
    for i in range(1, m):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j-1]; comparisons += 1
        comparisons += 1
        if pattern[i] == pattern[j]:
            j+=1; lps[i] = j
    # search
    hits = []; j = 0
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            j = lps[j-1]; comparisons += 1
        comparisons += 1
        if text[i] == pattern[j]:
            j += 1
        if j == m:
            hits.append(i-j+1)
            j = lps[j-1]
    return hits, comparisons

def boyer_moore_search(text, pattern):
    if not pattern: return [], 0
    comparisons = 0
    m, n = len(pattern), len(text)
    bad = {pattern[i]: i for i in range(m)}
    hits = []
    s = 0
    while s <= n - m:
        j = m-1
        while j >= 0:
            comparisons += 1
            if pattern[j] != text[s+j]:
                break
            j -= 1
        if j < 0:
            hits.append(s)
            s += (m - bad.get(text[s+m], -1)) if s + m < n else 1
        else:
            s += max(1, j - bad.get(text[s+j], -1))
    return hits, comparisons

def rabin_karp(text, pattern, d=256, q=101):
    m, n = len(pattern), len(text)
    if m == 0 or m > n: return [], 0
    comparisons = 0
    p = t = 0
    h = pow(d, m-1) % q
    hits = []
    for i in range(m):
        p = (d*p + ord(pattern[i])) % q
        t = (d*t + ord(text[i])) % q
    for s in range(n - m + 1):
        comparisons += 1
        if p == t and text[s:s+m] == pattern:
            hits.append(s)
            comparisons += m
        if s < n - m:
            t = (d*(t - ord(text[s])*h) + ord(text[s+m])) % q
            if t < 0: t += q
    return hits, comparisons

# ==========================
# Aho–Corasick (multi-pattern)
# ==========================
class AhoNode:
    def __init__(self):
        self.children = {}
        self.fail = None
        self.output = []

class AhoCorasick:
    def __init__(self, patterns):
        self.root = AhoNode()
        self.build_trie(patterns)
        self.build_failures()
    def build_trie(self, patterns):
        for pat in patterns:
            node = self.root
            for c in pat:
                if c not in node.children:
                    node.children[c] = AhoNode()
                node = node.children[c]
            node.output.append(pat)
    def build_failures(self):
        q = deque()
        for child in self.root.children.values():
            child.fail = self.root; q.append(child)
        while q:
            cur = q.popleft()
            for ch, nxt in cur.children.items():
                f = cur.fail
                while f and ch not in f.children:
                    f = f.fail
                nxt.fail = f.children[ch] if f and ch in f.children else self.root
                nxt.output += nxt.fail.output
                q.append(nxt)
    def search(self, text):
        node = self.root
        hits = defaultdict(list)
        comparisons = 0
        for i, c in enumerate(text):
            comparisons += 1
            while node and c not in node.children:
                node = node.fail; comparisons += 1
            if not node:
                node = self.root
                continue
            node = node.children[c]
            for pat in node.output:
                hits[pat].append(i - len(pat) + 1)
        return hits, comparisons

# ==========================
# PDF report function (optional)
# ==========================
def generate_pdf_report(title, df, figures):
    if not reportlab_available:
        raise RuntimeError("ReportLab not installed")
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4)
    styles = getSampleStyleSheet()
    story = [Paragraph(title, styles["Title"]), Paragraph(str(datetime.now()), styles["Normal"]), Spacer(1,12)]
    table_data = [list(df.columns)] + [[str(x) for x in row] for _, row in df.iterrows()]
    tbl = Table(table_data)
    tbl.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#00B4D8")),("TEXTCOLOR",(0,0),(-1,0),colors.white),("GRID",(0,0),(-1,-1),0.25,colors.grey)]))
    story.append(tbl); story.append(Spacer(1,12))
    for name, fig_bytes in figures:
        story.append(Paragraph(name, styles["Heading2"]))
        img = RLImage(io.BytesIO(fig_bytes), width=380, height=250)
        story.append(img); story.append(Spacer(1,12))
    doc.build(story)
    buf.seek(0)
    return buf.getvalue()

# ==========================
# Execution & session state
# ==========================
if "results_stored" not in st.session_state:
    st.session_state.results_stored = None
if "figures" not in st.session_state:
    st.session_state.figures = []

if st.button("🔍 Search Pattern"):

    if not sequences:
        st.warning("⚠️ No sequence loaded or entered.")
        st.stop()
    if not pattern_input:
        st.warning("⚠️ No pattern provided.")
        st.stop()
    if not selected_algos:
        st.warning("⚠️ No algorithm selected.")
        st.stop()

    st.session_state.figures = []
    patterns_list = [p.strip() for p in pattern_input.split(",") if p.strip()]
    all_results = []

    for header, dna_sequence in sequences.items():
        st.markdown(f"## 🧫 Results for **{header}** ({len(dna_sequence)} bp)")
        results = []

        for algo in selected_algos:

            # Aho-Corasick (multi-pattern)
            if algo == "Aho–Corasick":
                start = time.perf_counter()
                matches, comp = AhoCorasick(patterns_list).search(dna_sequence)
                elapsed = time.perf_counter() - start
                results.append({
                    "Sequence": header,
                    "Algorithm": algo,
                    "Pattern": "ALL",
                    "Matches": sum(len(v) for v in matches.values()),
                    "Comparisons": comp,
                    "Time (s)": round(elapsed, 5)
                })
                results.append({"Sequence":"", "Algorithm":"", "Pattern":"", "Matches":"", "Comparisons":"", "Time (s)":"",})
                continue

            # Single-pattern algorithms
            total_matches_list = []
            total_comp = 0
            total_time = 0

            for pat in patterns_list:
                if len(pat) == 0 or len(pat) > len(dna_sequence):
                    results.append({
                        "Sequence": header,
                        "Algorithm": algo,
                        "Pattern": pat,
                        "Matches": 0,
                        "Comparisons": 0,
                        "Time (s)": 0.0
                    })
                    continue

                start = time.perf_counter()
                matches, comp = {
                    "Naïve Search": naive_search,
                    "KMP": kmp_search,
                    "Boyer–Moore": boyer_moore_search,
                    "Rabin–Karp": rabin_karp
                }[algo](dna_sequence, pat)
                elapsed = time.perf_counter() - start

                results.append({
                    "Sequence": header,
                    "Algorithm": algo,
                    "Pattern": pat,
                    "Matches": len(matches),
                    "Comparisons": comp,
                    "Time (s)": round(elapsed, 5)
                })

                total_matches_list.append(len(matches))
                total_comp += comp
                total_time += elapsed

            results.append({
                "Sequence": header,
                "Algorithm": algo,
                "Pattern": "TOTAL",
                "Matches": sum(total_matches_list),
                "Comparisons": total_comp,
                "Time (s)": round(total_time, 5)
            })
            results.append({"Sequence":"", "Algorithm":"", "Pattern":"", "Matches":"", "Comparisons":"", "Time (s)":"",})

        df = pd.DataFrame(results)
        # Normalize algorithm names for display
        if "Algorithm" in df.columns:
            df["Algorithm"] = df["Algorithm"].apply(normalize_text).apply(fix_garbled)

        all_results.append(df)
        st.markdown("### 📊 Performance Results")
        st.dataframe(df, use_container_width=True)

        # Performance chart (TOTAL / ALL rows)
        df_chart = df[(df["Pattern"] == "TOTAL") | (df["Pattern"] == "ALL")]
        fig, ax = plt.subplots(figsize=(10,5))
        if not df_chart.empty:
            ax.bar(df_chart["Algorithm"], df_chart["Time (s)"], color="#00B4D8")
        ax.set_title(f"Algorithm Performance — {header}")
        ax.set_ylabel("Time (s)")
        ax.set_xticklabels(df_chart["Algorithm"], rotation=45, ha="right")
        plt.tight_layout()
        st.pyplot(fig)
        fname = f"perf_{re.sub(r'\\W+','_', header)}_{int(time.time())}.png"
        st.session_state.figures.append((fname, fig_to_bytes(fig)))

        # Pattern visualization (500 bp window)
        st.markdown("## 🔬 Pattern Visualizations (500 bp window)")
        for pat in patterns_list:
            st.markdown(f"### 🧬 Pattern: **{pat}**")
            matches, _ = naive_search(dna_sequence, pat)
            total_matches = len(matches)
            st.success(f"🔹 Found {total_matches} occurrence(s)")
            if total_matches == 0:
                continue
            start_pos = max(0, matches[0] - 50)
            end_pos = min(len(dna_sequence), start_pos + 500)
            window_seq = dna_sequence[start_pos:end_pos]
            st.info(f"📍 Region {start_pos} → {end_pos}")
            highlighted = window_seq.replace(
                pat,
                f"<span style='background-color:#ff4d6d; color:white; padding:2px; border-radius:4px;'>{pat}</span>"
            )
            st.markdown(f"<div style='background:#112233; padding:12px; border-radius:10px; font-family:monospace;'>{highlighted}</div>", unsafe_allow_html=True)
            st.info(f"⭐ {window_seq.count(pat)} match(es) inside this window")
            st.markdown("---")

    # combine results
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        if "Algorithm" in combined_df.columns:
            combined_df["Algorithm"] = combined_df["Algorithm"].apply(normalize_text).apply(fix_garbled)
        st.session_state.results_stored = combined_df

# ==========================
# Download section
# ==========================
if st.session_state.results_stored is not None:
    # CSV (utf-8-sig for Excel)
    csv_buf = io.StringIO()
    st.session_state.results_stored.to_csv(csv_buf, index=False)
    csv_bytes = csv_buf.getvalue().encode("utf-8-sig")
    st.download_button("📥 Download CSV", csv_bytes, "results.csv", mime="text/csv", key="csv")

    # ZIP charts
    if st.session_state.figures:
        zip_buf = io.BytesIO()
        with zipfile.ZipFile(zip_buf, "w") as z:
            for name, data in st.session_state.figures:
                z.writestr(name, data)
        zip_buf.seek(0)
        st.download_button("📦 Download Charts (ZIP)", zip_buf.getvalue(), "charts.zip", key="zip")

    # PDF
    if reportlab_available:
        pdf_bytes = generate_pdf_report("DNA Pattern Matching Report", st.session_state.results_stored, st.session_state.figures)
        st.download_button("📄 Download PDF Report", pdf_bytes, "dna_report.pdf", key="pdf")
    else:
        st.warning("Install reportlab to enable PDF export (pip install reportlab)")
