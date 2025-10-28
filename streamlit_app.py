import streamlit as st
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from io import StringIO
import base64

st.set_page_config(page_title="🧬 DNA Pattern Matching Tool", layout="wide")

st.title("🧬 DNA Pattern Matching Tool")
st.write("Upload one or more FASTA files to analyze DNA sequences, visualize length distribution, and export results.")

uploaded_files = st.file_uploader("Upload FASTA files", type=["fasta", "fa"], accept_multiple_files=True)

all_results = []

if uploaded_files:
    for file in uploaded_files:
        records = list(SeqIO.parse(file, "fasta"))
        st.subheader(f"📁 File: {file.name}")
        st.write(f"Total sequences: {len(records)}")

        # Extract details
        data = []
        for record in records:
            seq = str(record.seq).upper()
            data.append({
                "File": file.name,
                "Sequence_ID": record.id,
                "Length": len(seq),
                "GC_Content(%)": round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2) if len(seq) > 0 else 0
            })

        df = pd.DataFrame(data)
        all_results.append(df)

    # Combine all results safely
    combined_df = pd.concat(all_results, ignore_index=True)

    # Drop duplicates based on Sequence_ID + File
    combined_df = combined_df.drop_duplicates(subset=["Sequence_ID", "File"])

    st.subheader("📊 Combined Summary")
    st.dataframe(combined_df)

    # Visualization: Sequence length distribution
    st.subheader("📈 Sequence Length Distribution")
    fig, ax = plt.subplots(figsize=(8, 4))
    combined_df['Length'].hist(bins=20, ax=ax)
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Sequence Lengths")
    st.pyplot(fig)

    # Visualization: GC content per file
    st.subheader("🧫 Average GC Content per File")
    avg_gc = combined_df.groupby("File")["GC_Content(%)"].mean()
    st.bar_chart(avg_gc)

    # Export single CSV (no duplication)
    csv = combined_df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    st.markdown(
        f'<a href="data:file/csv;base64,{b64}" download="combined_sequences.csv">📥 Download Results CSV</a>',
        unsafe_allow_html=True
    )
else:
    st.info("Please upload one or more FASTA files to start the analysis.")
