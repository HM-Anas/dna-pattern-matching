# 🧬 DNA Pattern Matching Tool (KMP + Boyer–Moore)

An interactive web application that allows users to perform **DNA sequence pattern matching** using two classic algorithms — **Knuth–Morris–Pratt (KMP)** and **Boyer–Moore**.  

The app is built with **Streamlit** and provides an intuitive interface for searching specific DNA motifs or patterns within large sequences.  

---

## ⚙️ Features

✅ Upload a FASTA file or paste a DNA sequence manually  
✅ Search for custom DNA patterns (e.g., `ATG`, `CGTAA`)  
✅ Choose between **KMP** or **Boyer–Moore** algorithms  
✅ View **match counts**, **percent of total bases matched**, and **highlighted results**  
✅ Elegant **dark mode interface** with dynamic highlighting  

---

## 🧠 Algorithms Implemented

### 🔹 Knuth–Morris–Pratt (KMP)
- Efficient, linear-time pattern searching algorithm  
- Preprocesses the pattern using an LPS (Longest Prefix Suffix) table  

### 🔹 Boyer–Moore
- Searches patterns from right to left using bad-character heuristic  
- Often faster in practice for large text sequences  

---

## 🚀 How to Run Locally

### 1️⃣ Clone or Download the Repository
```bash
git clone https://github.com/HM-Anas/dna-pattern-matching.git
cd dna-pattern-matching
