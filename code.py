import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier

# ----------------------- Helper Functions -----------------------

def calculate_perplexity(sequence, k=3):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    kmer_counts = pd.Series(kmers).value_counts()
    total = sum(kmer_counts)
    probs = kmer_counts / total
    entropy = -sum(p * math.log2(p) for p in probs)
    perplexity = 2 ** entropy
    return perplexity

def detect_g_quadruplex(seq):
    pattern = r'(G{3,}\w{1,7}){3,}G{3,}'
    return len(re.findall(pattern, seq))

def detect_z_dna(seq):
    pattern = r'(CG){6,}'
    return len(re.findall(pattern, seq))

def detect_cruciform(seq):
    pattern = r'(.{4,6})[ATGC]{0,10}\1[::-1]'
    return len(re.findall(pattern, seq))

def detect_tata_box(seq):
    return len(re.findall(r'TATA[AT]A[AT]', seq))

def detect_direct_repeats(seq):
    return len(re.findall(r'(.{3,6})\1+', seq))

def extract_features(seq):
    seq = seq.upper()
    return {
        'perplexity': calculate_perplexity(seq),
        'g_quadruplex': detect_g_quadruplex(seq),
        'z_dna': detect_z_dna(seq),
        'cruciform': detect_cruciform(seq),
        'tata_box': detect_tata_box(seq),
        'direct_repeats': detect_direct_repeats(seq),
        'length': len(seq)
    }

# ----------------------- Streamlit App -----------------------

st.title("📄 DNA Motif & Perplexity Analyzer from TXT File")

uploaded_file = st.file_uploader("📂 Upload a .txt file with DNA sequences", type=["txt"])

if uploaded_file:
    st.success(f"✅ Uploaded: {uploaded_file.name}")

    content = uploaded_file.read().decode("utf-8").strip().splitlines()

    records = []
    for i, line in enumerate(content):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) == 2:
            seq_id, sequence = parts
        else:
            seq_id = f"Seq_{i+1}"
            sequence = parts[0]
        records.append((seq_id, sequence))

    st.info(f"🔍 Found {len(records)} sequences")

    feature_rows = []
    for seq_id, sequence in records:
        features = extract_features(sequence)
        features["ID"] = seq_id
        feature_rows.append(features)

    df = pd.DataFrame(feature_rows)

    clf = RandomForestClassifier()
    dummy_data = df.drop(columns=["ID"])
    clf.fit(dummy_data, [0] * len(df))

    predictions = clf.predict(dummy_data)
    df["Predicted_Region"] = predictions

    st.subheader("🔬 Analysis Results")
    st.dataframe(df)

    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button("📥 Download CSV", data=csv, file_name="txt_sequence_results.csv", mime="text/csv")
