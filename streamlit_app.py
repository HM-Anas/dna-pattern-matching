# ================================
# DNA Pattern Matching Analyzer🧬 
# ================================
# This application compares different string matching algorithms for DNA sequence analysis
# It allows users to upload FASTA files or input DNA sequences manually and search for patterns

#Importing Libraries
import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import deque, defaultdict
from math import log2

# ----------------------------
# Page config
# ----------------------------
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

# ==========================
# STYLING
# ==========================
# Custom CSS styling for dark theme with blue accent colors
# Makes the app visually appealing with proper color schemes for DNA analysis
st.markdown("""
    <style>
        body, .stApp { background-color: #0E1117; color: #FAFAFA; }
        h1, h2, h3 { color: #00B4D8 !important; }
        .stButton>button {
            background-color: #00B4D8; color: white; font-weight: bold;
            border-radius: 8px; border: none; padding: 0.6em 1.2em;
        }
        .stButton>button:hover { background-color: #0077B6; }
        .result-box { background-color: #1E2636; padding: 15px; border-radius: 10px; border: 1px solid #00B4D8;
                      font-family: monospace; line-height: 1.6; word-wrap: break-word; }
        .highlight { color: #FFD60A; font-weight: bold; }
    </style>
""", unsafe_allow_html=True)

st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Enhanced: reverse complement, stats, density plots, and PNG export</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ----------------------------
# Utility functions
# ----------------------------
def clean_sequence(s):
    return re.sub(r'[^ATCG]', '', s.upper())

def reverse_complement(seq):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    return "".join(comp.get(b,'N') for b in seq[::-1])

def gc_content(seq):
    if len(seq)==0: return 0.0
    g = seq.count('G'); c = seq.count('C')
    return (g + c) / len(seq) * 100

def at_gc_ratio(seq):
    a = seq.count('A'); t = seq.count('T'); g = seq.count('G'); c = seq.count('C')
    denom = (g+c) if (g+c)>0 else 1
    return (a+t)/denom

def nucleotide_freq(seq):
    L = len(seq) if len(seq)>0 else 1
    return {b: seq.count(b)/L*100 for b in ['A','T','C','G']}

def shannon_entropy(seq):
    L = len(seq)
    if L==0: return 0.0
    freqs = [seq.count(b)/L for b in ['A','T','C','G']]
    ent = -sum(p*log2(p) for p in freqs if p>0)
    return ent

# ----------------------------
# Algorithms (as in your code)
# ----------------------------
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
    if len(pattern)==0:
        return [], comparisons
    lps=[0]*len(pattern); j=0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]:
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        if pattern[i]==pattern[j]:
            j+=1; lps[i]=j
    res=[]; j=0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]:
            j=lps[j-1]; comparisons += 1
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
    if m==0: return [], comparisons
    bad_char={pattern[i]:i for i in range(m)}
    res=[]; s=0
    while s<=n-m:
        j=m-1
        while j>=0:
            comparisons += 1
            if pattern[j]!=text[s+j]:
                break
            j-=1
        if j<0:
            res.append(s)
            s += (m - bad_char.get(text[s+m], -1)) if s+m < n else 1
        else:
            s += max(1, j - bad_char.get(text[s+j], -1))
    return res, comparisons

def rabin_karp(text, pattern, d=256, q=101):
    comparisons = 0
    m,n=len(pattern),len(text)
    if m==0 or n < m:
        return [], comparisons
    p=t=0; res=[]
    h=pow(d,m-1)%q
    for i in range(m):
        p=(d*p+ord(pattern[i]))%q
        t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        comparisons += 1
        if p==t:
            if text[s:s+m]==pattern:
                res.append(s)
                comparisons += m
        if s<n-m:
            t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q
            if t<0: t+=q
    return res, comparisons

class AhoNode:
    def __init__(self):
        self.children={}
        self.fail=None
        self.output=[]

class AhoCorasick:
    def __init__(self, patterns):
        self.root=AhoNode()
        self.build_trie(patterns)
        self.build_failure_links()
    def build_trie(self, patterns):
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
            for c, child in current.children.items():
                f=current.fail
                while f and c not in f.children:
                    f=f.fail
                child.fail = f.children[c] if f and c in f.children else self.root
                child.output += child.fail.output
                queue.append(child)
    def search(self, text):
        node=self.root
        res=defaultdict(list)
        comparisons=0
        for i,c in enumerate(text):
            comparisons += 1
            while node and c not in node.children:
                node=node.fail
                comparisons += 1
            if not node:
                node=self.root
                continue
            node=node.children[c]
            for pat in node.output:
                res[pat].append(i-len(pat)+1)
        return res, comparisons

algo_funcs = {
    "Naïve Search": naive_search,
    "KMP": kmp_search,
    "Boyer–Moore": boyer_moore_search,
    "Rabin–Karp": rabin_karp,
    "Aho–Corasick": None  # handled separately
}

# ----------------------------
# Inputs: upload / manual / patterns / options
# ----------------------------
uploaded_files = st.file_uploader("📁 Upload FASTA files (multiple allowed)", type=["fasta","fa","txt"], accept_multiple_files=True)
sequences = {}
if uploaded_files:
    for uploaded_file in uploaded_files:
        fasta = uploaded_file.getvalue().decode("utf-8")
        lines = fasta.strip().splitlines()
        header = lines[0] if lines and lines[0].startswith(">") else uploaded_file.name
        seq = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
        seq = clean_sequence(seq)
        sequences[header] = seq
        st.success(f"✅ Loaded: {header} ({len(seq)} bp)")
else:
    seq_input = st.text_area("🧬 Enter DNA Sequence (manual)", placeholder="ATCGGATCGATCG...", height=120).strip().upper()
    if seq_input:
        sequences["Manual Entry"] = clean_sequence(seq_input)

pattern_input = st.text_input("🔍 Enter Pattern(s) (comma-separated)", placeholder="CGA,ATG").strip().upper()
selected_algos = st.multiselect("⚙️ Select Algorithms", ["Naïve Search","KMP","Boyer–Moore","Rabin–Karp","Aho–Corasick"], default=["KMP","Boyer–Moore","Rabin–Karp"])
revcomp_checkbox = st.checkbox("🔄 Search Reverse Complement of patterns (auto add)")

# Session state for results
if "results_stored" not in st.session_state:
    st.session_state.results_stored = None

# ----------------------------
# Tabs layout (Option C)
# ----------------------------
tab_stats, tab_results, tab_visuals, tab_downloads = st.tabs(["Sequence Stats", "Match Results", "Visualizations", "Charts & Downloads"])

# ----------------------------
# Sequence Stats Tab
# ----------------------------
with tab_stats:
    st.markdown("## 📋 Sequence Statistics")
    if not sequences:
        st.info("Upload FASTA or enter sequence to view statistics.")
    else:
        for header, seq in sequences.items():
            st.markdown(f"### 🔖 {header} (length: {len(seq)} bp)")
            col1, col2, col3 = st.columns([1,1,1])
            gc = gc_content(seq)
            atgc = at_gc_ratio(seq)
            ent = shannon_entropy(seq)
            col1.metric("Length (bp)", len(seq))
            col2.metric("GC Content (%)", f"{gc:.2f}%")
            col3.metric("Shannon Entropy", f"{ent:.3f} bits")
            st.markdown("**Nucleotide Frequencies (%)**")
            freqs = nucleotide_freq(seq)
            freq_df = pd.DataFrame.from_dict(freqs, orient='index', columns=['Percent']).reset_index().rename(columns={'index':'Base'})
            st.dataframe(freq_df, use_container_width=False)
            # Pie chart
            fig, ax = plt.subplots(figsize=(4,3))
            ax.pie([freq_df['Percent'][i] for i in range(len(freq_df))], labels=freq_df['Base'], autopct='%1.1f%%')
            ax.set_title("Nucleotide Composition")
            st.pyplot(fig)

# ----------------------------
# Run analysis when button clicked (compute results)
# ----------------------------
if st.button("🔍 Run Analysis"):
    if not sequences or not pattern_input:
        st.warning("Please provide sequence(s) and pattern(s).")
    else:
        # prepare patterns list
        user_patterns = [p.strip() for p in pattern_input.split(",") if p.strip()]
        patterns_list = list(user_patterns)  # copy
        if revcomp_checkbox:
            for p in user_patterns:
                rc = reverse_complement(p)
                if rc not in patterns_list:
                    patterns_list.append(rc)
        # store results per sequence
        all_results = []
        perf_summary_records = []  # for performance chart
        density_data_store = {}  # store densities for visuals: { (header,pattern): density_array }
        matches_store = {}  # store matches details for visualizations and download

        for header, seq in sequences.items():
            st.markdown(f"## 🧫 Results for **{header}** ({len(seq)} bp)")
            results = []
            # Aho-Corasick handled separately if selected
            if "Aho–Corasick" in selected_algos:
                if len(patterns_list) < 2:
                    st.warning("Aho–Corasick requires >=2 patterns; skipped for this run.")
                else:
                    start = time.time()
                    ac = AhoCorasick(patterns_list)
                    matches_ac, comps_ac = ac.search(seq)
                    elapsed = time.time() - start
                    total_matches_ac = sum(len(v) for v in matches_ac.values())
                    results.append({"Algorithm":"Aho–Corasick","Pattern":"ALL PATTERNS","Matches": total_matches_ac, "Comparisons": comps_ac, "Time (s)": round(elapsed,5)})
                    perf_summary_records.append({"Sequence": header, "Algorithm":"Aho–Corasick", "Time (s)": round(elapsed,5)})
                    # store detailed matches
                    matches_store[(header,"Aho–Corasick")] = matches_ac
                    # Also store density for each pattern from Aho results
                    for pat in patterns_list:
                        positions = matches_ac.get(pat, [])
                        density = np.zeros(len(seq), dtype=int)
                        for pos in positions:
                            if 0 <= pos < len(seq):
                                density[pos] += 1
                        density_data_store[(header, pat)] = density

            # Single-pattern algorithms
            for algo in ["Naïve Search","KMP","Boyer–Moore","Rabin–Karp"]:
                if algo not in selected_algos: continue
                pattern_match_counts = []
                total_comp = 0
                total_time = 0.0
                pattern_matches_detail = {}
                for pat in patterns_list:
                    start = time.time()
                    matches_list, comps = algo_funcs[algo](seq, pat)
                    elapsed = time.time() - start
                    pattern_matches_detail[pat] = matches_list
                    total_comp += comps
                    total_time += elapsed
                    results.append({"Algorithm": algo, "Pattern": pat, "Matches": len(matches_list), "Comparisons": comps, "Time (s)": round(elapsed,5)})
                    perf_summary_records.append({"Sequence": header, "Algorithm": f"{algo} ({pat})", "Time (s)": round(elapsed,5)})
                    # density for this pattern:
                    density = np.zeros(len(seq), dtype=int)
                    for pos in matches_list:
                        if 0 <= pos < len(seq):
                            density[pos] += 1
                    density_data_store[(header, pat)] = density
                # add TOTAL row for algorithm summarizing all patterns
                results.append({"Algorithm": algo, "Pattern": "TOTAL", "Matches": sum(len(v) for v in pattern_matches_detail.values()), "Comparisons": total_comp, "Time (s)": round(total_time,5)})
                # store pattern matches detail
                matches_store[(header, algo)] = pattern_matches_detail

            # show dataframe
            df_results = pd.DataFrame(results)
            all_results.append(df_results)
            st.markdown("### 📊 Results Table")
            st.dataframe(df_results, use_container_width=True)

        # combine results for later download
        combined_df = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
        st.session_state.results_stored = combined_df
        st.session_state.density_data = density_data_store
        st.session_state.matches_store = matches_store
        st.session_state.perf_summary = pd.DataFrame(perf_summary_records)

        st.success("✅ Analysis complete. Navigate to other tabs for visualizations and downloads.")

# ----------------------------
# Match Results Tab (show combined table and allow CSV download)
# ----------------------------
with tab_results:
    st.markdown("## 📋 Combined Results")
    if st.session_state.results_stored is None:
        st.info("Run analysis to see results.")
    else:
        st.dataframe(st.session_state.results_stored, use_container_width=True)
        csv_buf = io.StringIO()
        st.session_state.results_stored.to_csv(csv_buf, index=False)
        st.download_button("📥 Download Combined Results (CSV)", data=csv_buf.getvalue(), file_name="dna_pattern_results.csv", mime="text/csv")

# ----------------------------
# Visualizations Tab
# ----------------------------
with tab_visuals:
    st.markdown("## 🎨 Visualizations")
    if st.session_state.results_stored is None:
        st.info("Run analysis to generate visualizations.")
    else:
        # Performance chart (summary)
        st.markdown("### 📈 Performance Summary (per algorithm / pattern)")
        perf_df = st.session_state.perf_summary.copy() if "perf_summary" in st.session_state else pd.DataFrame()
        if not perf_df.empty:
            fig1, ax1 = plt.subplots(figsize=(10,4))
            # aggregate mean time by Algorithm label (this label already includes pattern for single algorithms)
            summary = perf_df.groupby("Algorithm")["Time (s)"].mean().sort_values(ascending=False)
            ax1.bar(summary.index, summary.values, color="#00B4D8")
            ax1.set_ylabel("Time (s)")
            ax1.set_xticklabels(summary.index, rotation=45, ha='right', fontsize=9)
            ax1.set_title("Average Execution Time")
            plt.tight_layout()
            st.pyplot(fig1)

            # PNG download for performance chart
            buf = io.BytesIO()
            fig1.savefig(buf, format="png", bbox_inches='tight')
            buf.seek(0)
            st.download_button("📥 Download Performance Chart (PNG)", data=buf.getvalue(), file_name="performance_chart.png", mime="image/png")
            plt.close(fig1)
        else:
            st.info("No performance summary available.")

        st.markdown("---")
        st.markdown("### 🔍 Match Density Plots (one per pattern)")
        # iterate over stored densities and plot separate plots per (sequence, pattern)
        density_store = st.session_state.get("density_data", {})
        if not density_store:
            st.info("No density data available. Run analysis first.")
        else:
            # group by sequence header
            headers = sorted(list({k[0] for k in density_store.keys()}))
            for header in headers:
                st.markdown(f"#### Sequence: {header}")
                # find patterns for this header
                patterns_for_header = [k[1] for k in density_store.keys() if k[0]==header]
                for pat in patterns_for_header:
                    density = density_store[(header, pat)]
                    if density.sum() == 0:
                        st.warning(f"No matches for pattern {pat} in sequence {header}.")
                        continue
                    fig2, ax2 = plt.subplots(figsize=(10,2))
                    # smooth density with rolling window (e.g., window=50) for visualization clarity
                    window = max(1, min(1000, len(density)//200))  # heuristic window
                    smoothed = np.convolve(density, np.ones(window)/window, mode='same')
                    ax2.plot(smoothed, color='#FF6B6B')
                    ax2.set_title(f"Density plot for pattern: {pat}")
                    ax2.set_xlabel("Base position")
                    ax2.set_ylabel("Match density (smoothed)")
                    plt.tight_layout()
                    st.pyplot(fig2)

                    # PNG download for this density plot
                    buf2 = io.BytesIO()
                    fig2.savefig(buf2, format="png", bbox_inches='tight')
                    buf2.seek(0)
                    safe_name = f"{header.replace(' ','_')}_{pat}_density.png"
                    st.download_button(f"📥 Download '{pat}' density (PNG)", data=buf2.getvalue(), file_name=safe_name, mime="image/png")
                    plt.close(fig2)

# ----------------------------
# Charts & Downloads Tab
# ----------------------------
with tab_downloads:
    st.markdown("## 📦 Charts & Downloads")
    if st.session_state.results_stored is None:
        st.info("Run analysis to enable downloads.")
    else:
        # Download combined CSV again (handy)
        csv_buf = io.StringIO()
        st.session_state.results_stored.to_csv(csv_buf, index=False)
        st.download_button("📥 Download Results CSV", data=csv_buf.getvalue(), file_name="dna_pattern_results.csv", mime="text/csv")

        # Option: export highlighted sequences as simple annotated FASTA (basic)
        if st.button("📄 Export annotated FASTA (basic)"):
            fasta_out = []
            # create a simple annotated FASTA: header -> positions per pattern
            for (header, algo), matchdict in st.session_state.matches_store.items():
                # for Aho-Corasick matchdict is dict of patterns->positions
                fasta_out.append(f">{header} | algorithm: {algo}")
                # write a small annotation line per pattern (positions)
                if isinstance(matchdict, dict):
                    for pat, positions in matchdict.items():
                        fasta_out.append(f";{pat}:" + ",".join(map(str, positions)))
                else:
                    # Aho-case: matchdict might already be mapping
                    fasta_out.append(str(matchdict))
            fasta_text = "\n".join(fasta_out)
            st.download_button("📥 Download annotated FASTA (txt)", data=fasta_text, file_name="annotated_matches.txt", mime="text/plain")

        st.info("You can download individual charts from the Visualizations tab as PNG files.")

## End of app
