import streamlit as st
import pandas as pd
from Bio import SeqIO
import re
import joblib
import math
from sklearn.ensemble import RandomForestClassifier

# ----------------------- Helper Functions -----------------------

# Function to calculate perplexity using k-mer model
def calculate_perplexity(sequence, k=3):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    kmer_counts = pd.Series(kmers).value_counts()
    total = sum(kmer_counts)
    probs = kmer_counts / total
    entropy = -sum(p * math.log2(p) for p in probs)
    perplexity = 2 ** entropy
    return perplexity

# Motif detection
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

# Generate feature vector for a given sequence
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

st.title("🧬 Functional DNA Region Identifier (Perplexity + Motifs)")
uploaded_file = st.file_uploader("Upload your FASTA file (.fasta, .fa, .fna)", type=["fasta", "fa", "fna"])

if uploaded_file:
    st.success("File uploaded successfully!")

    records = list(SeqIO.parse(uploaded_file, "fasta"))
    st.info(f"Found {len(records)} sequences")

    feature_rows = []

    for record in records:
        seq = str(record.seq)
        features = extract_features(seq)
        features["ID"] = record.id
        feature_rows.append(features)

    df = pd.DataFrame(feature_rows)

    # Load a pretrained model or create a dummy classifier
    # In production, train the model separately and load using joblib
    # For demo purpose: simple RandomForestClassifier
    if "model.pkl" not in st.session_state:
        clf = RandomForestClassifier()
        dummy_data = df.drop(columns=["ID"])
        clf.fit(dummy_data, [0]*len(df))  # Dummy fit
        joblib.dump(clf, "model.pkl")

    model = joblib.load("model.pkl")
    
    X = df.drop(columns=["ID"])
    predictions = model.predict(X)
    df["Predicted_Region"] = predictions

    st.subheader("🔬 Analysis Results")
    st.dataframe(df)

    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button("Download Results as CSV", data=csv, file_name="region_prediction_results.csv", mime="text/csv")
