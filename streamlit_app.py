# ================================
# DNA Pattern Matching Analyzer🧬 
# ================================
# This application compares different string matching algorithms for DNA sequence analysis
# It allows users to upload FASTA files or input DNA sequences manually and search for patterns

# Importing Libraries
import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque, defaultdict
import zipfile
from datetime import datetime

# ReportLab for PDF generation
try:
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib import colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
        Table, TableStyle
    )
    from reportlab.lib.styles import getSampleStyleSheet
    reportlab_available = True
except:
    reportlab_available = False


# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

# ==========================
# STYLING
# ==========================
st.markdown("""
    <style>
        body, .stApp { background-color: #0E1117; color: #FAFAFA; }
        h1, h2, h3 { color: #00B4D8 !important; }
        .stButton>button {
            background-color: #00B4D8; color: white; font-weight: bold;
            border-radius: 8px; border: none; padding: 0.6em 1.2em;
        }
        .stButton>button:hover { background-color: #0077B6; }
        .result-box { background-color: #1E2636; padding: 15px; border-radius: 10px;
                      border: 1px solid #00B4D8; font-family: monospace;
                      line-height: 1.6; word-wrap: break-word; }
    </style>
""", unsafe_allow_html=True)

# ==========================
# HEADER
# ==========================
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare string matching algorithms for DNA sequences</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ==========================
# FILE UPLOAD / INPUT
# ==========================
uploaded_files = st.file_uploader(
    "📁 Upload FASTA files (you can select multiple)",
    type=["fasta","fa","txt"],
    accept_multiple_files=True
)

sequences = {}

if uploaded_files:
    for uploaded_file in uploaded_files:
        fasta = uploaded_file.getvalue().decode("utf-8")
        lines = fasta.strip().split("\n")
        header = lines[0] if lines[0].startswith(">") else uploaded_file.name
        sequence = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
        sequence = re.sub(r'[^ATCG]', '', sequence)
        sequences[header] = sequence
        st.success(f"✅ Loaded: {header} ({len(sequence)} bp)")
else:
    seq_input = st.text_area("🧬 Enter DNA Sequence", placeholder="ATCGGATCGATCG...", height=120).strip().upper()
    if seq_input:
        sequences["Manual Entry"] = re.sub(r'[^ATCG]', '', seq_input)

pattern_input = st.text_input(
    "🔍 Enter Pattern(s) to Search (comma separated for multiple)",
    placeholder="CGATCGA,ATGCGT"
).strip().upper()

algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=algorithms)


# ============
# ALGORITHMS
# ============

def naive_search(text, pattern):
    comparisons = 0
    results = []
    for i in range(len(text)-len(pattern)+1):
        comparisons += 1
        match = True
        for j in range(len(pattern)):
            comparisons += 1
            if text[i+j] != pattern[j]:
                match = False
                break
        if match:
            results.append(i)
    return results, comparisons

def kmp_search(text, pattern):
    comparisons = 0
    lps=[0]*len(pattern)
    j=0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]:
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        if pattern[i]==pattern[j]:
            j+=1
            lps[i]=j
    res=[]
    j=0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]:
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        if text[i]==pattern[j]:
            j+=1
        if j==len(pattern):
            res.append(i-j+1)
            j=lps[j-1]
    return res, comparisons

def boyer_moore_search(text, pattern):
    comparisons = 0
    m,n=len(pattern),len(text)
    bad_char={pattern[i]:i for i in range(m)}
    res=[]
    s=0
    while s<=n-m:
        j=m-1
        while j>=0:
            comparisons+=1
            if pattern[j]!=text[s+j]:
                break
            j-=1
        if j<0:
            res.append(s)
            s+=(m-bad_char.get(text[s+m],-1)) if s+m<n else 1
        else:
            s+=max(1,j-bad_char.get(text[s+j],-1))
    return res, comparisons

def rabin_karp(text, pattern, d=256, q=101):
    comparisons = 0
    m,n=len(pattern),len(text)
    p=t=0
    h=pow(d,m-1)%q
    res=[]
    for i in range(m):
        p=(d*p+ord(pattern[i]))%q
        t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        comparisons+=1
        if p==t:
            if text[s:s+m]==pattern:
                res.append(s)
                comparisons+=m
        if s<n-m:
            t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q
            t+=q if t<0 else 0
    return res, comparisons


# =====================
# AHO–CORASICK
# =====================
class AhoNode:
    def __init__(self):
        self.children={}
        self.fail=None
        self.output=[]

class AhoCorasick:
    def __init__(self,patterns):
        self.root=AhoNode()
        self.build_trie(patterns)
        self.build_failure_links()

    def build_trie(self,patterns):
        for pat in patterns:
            node=self.root
            for c in pat:
                if c not in node.children:
                    node.children[c]=AhoNode()
                node=node.children[c]
            node.output.append(pat)

    def build_failure_links(self):
        queue=deque()
        for child in self.root.children.values():
            child.fail=self.root
            queue.append(child)
        while queue:
            current=queue.popleft()
            for c,child in current.children.items():
                f=current.fail
                while f and c not in f.children:
                    f=f.fail
                child.fail=f.children[c] if f and c in f.children else self.root
                child.output+=child.fail.output
                queue.append(child)

    def search(self,text):
        node=self.root
        res=defaultdict(list)
        comparisons = 0
        for i,c in enumerate(text):
            comparisons+=1
            while node and c not in node.children:
                node=node.fail
                comparisons+=1
            if not node:
                node=self.root
                continue
            node=node.children[c]
            for pat in node.output:
                res[pat].append(i-len(pat)+1)
        return res, comparisons


# ==========================
# HELPER FUNCTIONS
# ==========================
def fig_to_bytes(fig, format="png"):
    buf = io.BytesIO()
    fig.savefig(buf, format=format, dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()

def generate_pdf_report(title, results_df, figures):
    if not reportlab_available:
        raise RuntimeError("Reportlab library not installed.")
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []
    story.append(Paragraph(title, styles["Title"]))
    story.append(Paragraph(f"Generated on: {datetime.now()}", styles["Normal"]))
    story.append(Spacer(1, 12))

    table_data = [list(results_df.columns)]
    for _, row in results_df.head(50).iterrows():
        table_data.append([str(x) for x in row])

    table = Table(table_data)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#00B4D8")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
    ]))

    story.append(table)
    story.append(Spacer(1, 12))

    for name, fig_bytes in figures:
        story.append(Paragraph(f"Figure: {name}", styles["Heading2"]))
        img_buf = io.BytesIO(fig_bytes)
        img = RLImage(img_buf, width=400, height=250)
        story.append(img)
        story.append(Spacer(1, 12))

    doc.build(story)
    buffer.seek(0)
    return buffer.getvalue()


# ==========================
# RUN ANALYSIS
# ==========================
if "results_stored" not in st.session_state:
    st.session_state.results_stored=None

if "figures" not in st.session_state:
    st.session_state.figures = []

if st.button("🔍 Search Pattern"):

    if not sequences or not pattern_input:
        st.warning("⚠️ Please enter sequence(s) and pattern(s).")

    else:
        patterns_list = [p.strip() for p in pattern_input.split(",") if p.strip()]
        all_results=[]

        for header,dna_sequence in sequences.items():
            st.markdown(f"## 🧫 Results for **{header}** ({len(dna_sequence)} bp)")
            results=[]

            for algo in selected_algos:

                if algo=="Aho–Corasick":
                    if len(patterns_list)<2:
                        st.warning("⚠️ Aho–Corasick requires multiple patterns.")
                        continue
                    start=time.time()
                    matches, comparisons = AhoCorasick(patterns_list).search(dna_sequence)
                    elapsed=time.time()-start
                    total_matches = sum(len(v) for v in matches.values())
                    results.append({
                        "Algorithm": algo,
                        "Pattern": "ALL PATTERNS",
                        "Matches": total_matches,
                        "Comparisons": comparisons,
                        "Time (s)": round(elapsed, 5)
                    })
                    results.append({"Algorithm":"","Pattern":"","Matches":"","Comparisons":"","Time (s)":""})

                else:
                    pattern_results=[]
                    total_comparisons=0
                    total_time=0

                    for pat in patterns_list:
                        start_pat=time.time()
                        matches, comp = {
                            "Naïve Search": naive_search,
                            "KMP": kmp_search,
                            "Boyer–Moore": boyer_moore_search,
                            "Rabin–Karp": rabin_karp
                        }[algo](dna_sequence, pat)
                        elapsed_pat=time.time()-start_pat

                        results.append({
                            "Algorithm": algo,
                            "Pattern": pat,
                            "Matches": len(matches),
                            "Comparisons": comp,
                            "Time (s)": round(elapsed_pat, 5)
                        })
                        pattern_results.append(len(matches))
                        total_comparisons+=comp
                        total_time+=elapsed_pat

                    results.append({
                        "Algorithm": algo,
                        "Pattern": "TOTAL",
                        "Matches": sum(pattern_results),
                        "Comparisons": total_comparisons,
                        "Time (s)": round(total_time, 5)
                    })
                    results.append({"Algorithm":"","Pattern":"","Matches":"","Comparisons":"","Time (s)":""})

            df=pd.DataFrame(results)
            all_results.append(df)

            st.markdown("### 📊 Performance Results")
            st.dataframe(df, use_container_width=True)

            st.markdown("### 📈 Performance Chart")
            df_chart = df[(df["Pattern"] == "TOTAL") | (df["Pattern"] == "ALL PATTERNS")]
            df_chart = df_chart[df_chart["Algorithm"] != ""]

            fig, ax = plt.subplots(figsize=(10,5))
            ax.bar(df_chart["Algorithm"], df_chart["Time (s)"], color="#00B4D8")
            ax.set_ylabel("Time (s)")
            ax.set_title(f"Algorithm Performance — {header}")
            ax.set_xticklabels(df_chart["Algorithm"], rotation=45, ha='right')
            plt.tight_layout()

            st.pyplot(fig)

            st.session_state.figures.append((
                f"performance_chart_{header}.png",
                fig_to_bytes(fig, "png")
            ))

        combined_df=pd.concat(all_results, ignore_index=True)
        st.session_state.results_stored=combined_df


# ==========================
# DOWNLOAD CSV
# ==========================
if st.session_state.results_stored is not None:
    csv_buffer=io.StringIO()
    st.session_state.results_stored.to_csv(csv_buffer,index=False)

    st.download_button(
        label="📥 Download Results as CSV",
        data=csv_buffer.getvalue(),
        file_name="dna_pattern_results.csv",
        mime="text/csv"
    )


# ==========================
# EXPORT OPTIONS
# ==========================
if st.session_state.results_stored is not None:

    st.markdown("## 📄 Export Options")

    # ZIP Export
    if st.session_state.figures:
        zip_buf=io.BytesIO()
        with zipfile.ZipFile(zip_buf, "w") as z:
            for name, fig_bytes in st.session_state.figures:
                z.writestr(name, fig_bytes)
        zip_buf.seek(0)

        st.download_button(
            "📦 Download All Charts (ZIP)",
            data=zip_buf,
            file_name="charts_bundle.zip"
        )

    # Individual PNG Downloads
    st.markdown("### 📉 Individual Chart Downloads")
    for name, fig_bytes in st.session_state.figures:
        st.download_button(
            f"📥 Download {name}",
            data=fig_bytes,
            file_name=name,
            mime="image/png"
        )

    # PDF Export
    if reportlab_available:
        pdf_bytes = generate_pdf_report(
            "DNA Pattern Matching Analyzer — Full Report",
            st.session_state.results_stored,
            st.session_state.figures
        )

        st.download_button(
            "📄 Download Full PDF Report",
            data=pdf_bytes,
            file_name="dna_report.pdf",
            mime="application/pdf"
        )
    else:
        st.warning("⚠ Install reportlab to enable PDF export: pip install reportlab")
