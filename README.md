# 🧬 DNA Pattern Matching using Multiple String-Matching Algorithms


An interactive bioinformatics web tool built with Streamlit that enables users to perform DNA sequence pattern matching using multiple string matching algorithms — including Naïve, KMP, Boyer–Moore, Rabin–Karp, and Aho–Corasick.

This application supports multiple FASTA uploads, provides runtime comparisons, and visually highlights matching motifs within each sequence in an elegant dark mode interface.

⚙️ Features

✅ Upload multiple FASTA files simultaneously
✅ Or manually enter a custom DNA sequence
✅ Search for custom DNA motifs (e.g., ATGCGT, CGTAA)
✅ Choose one or more algorithms to compare
✅ View detailed results:

Number of matches

Runtime (execution time in seconds)

Positions of matches
✅ Highlighted sequence visualization showing match regions
✅ Bar chart comparison of algorithm performance
✅ Sleek dark mode UI with interactive elements

🧠 Algorithms Implemented
🔹 Naïve Search

A straightforward approach that checks for the pattern at every position of the sequence.

🔹 Knuth–Morris–Pratt (KMP)

Linear-time pattern searching algorithm

Uses the Longest Prefix Suffix (LPS) preprocessing table to skip redundant comparisons

🔹 Boyer–Moore

Efficient for large sequences

Uses the bad character heuristic to skip unnecessary comparisons

🔹 Rabin–Karp

Hash-based string matching algorithm

Compares hash values instead of characters for faster pattern search

🔹 Aho–Corasick

Multi-pattern search algorithm

Currently simplified here for single-pattern demonstration

🚀 How to Run Locally
1️⃣ Clone or Download the Repository
git clone https://github.com/HM-Anas/dna-pattern-matching.git
cd dna-pattern-matching

2️⃣ Install Dependencies

Make sure you have Python 3.9+ and pip installed.

pip install streamlit pandas matplotlib

3️⃣ Run the Streamlit App
streamlit run app.py

4️⃣ Open in Browser

The app will open automatically at:
👉 http://localhost:8501

🌐 Live Demo

🔗 Launch the App on Streamlit Cloud

🧩 Example Usage

Upload one or more .fasta or .fa files

Enter a pattern like ATGCGT

Select algorithms (e.g., KMP, Boyer–Moore)

Click Search Pattern

View highlighted matches, runtime comparison table, and performance chart
