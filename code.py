import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier

# ----------------------- Styling -----------------------
st.set_page_config(page_title="DNA Sequence Analyzer", layout="wide")
st.markdown("""
    <style>
    .main {
        background-color: #f0f2f6;
        padding: 1rem;
    }
    .title {
        text-align: center;
        font-size: 2.5rem;
        color: #2c3e50;
        padding-bottom: 10px;
    }
    .section {
        background-color: #ffffff;
        padding: 20px;
        border-radius: 15px;
        box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.1);
        margin-top: 20px;
    }
    .stButton>button {
        background-color: #4CAF50;
        color: white;
        border-radius: 10px;
        padding: 0.5em 1em;
        font-size: 1em;
        border: none;
    }
    .stDownloadButton>button {
        background-color: #e67e22;
        color: white;
        font-weight: bold;
        border-radius: 10px;
    }
    </style>
""", unsafe_allow_html=True)

st.markdown('<div class="title">🧬 Functional DNA Region Analyzer</div>', unsafe_allow_html=True)

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

def reverse_complement(seq):
    complement = str.maketrans('ATGCatgc', 'TACGtacg')
    return seq.translate(complement)[::-1]

def detect_cruciform(seq, min_len=4, max_len=6, spacer=10):
    count = 0
    seq = seq.upper()
    for size in range(min_len, max_len + 1):
        for i in range(len(seq) - 2 * size - spacer + 1):
            left = seq[i:i + size]
            right = seq[i + size + spacer:i + 2 * size + spacer]
            if reverse_complement(left) == right:
                count += 1
    return count

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

# ----------------------- File Upload -----------------------

with st.container():
    st.markdown('<div class="section">', unsafe_allow_html=True)

    uploaded_file = st.file_uploader("📄 Upload a `.txt` file with DNA sequences (each line = 1 sequence)", type=["txt"])

    if uploaded_file:
        st.success(f"✅ File uploaded: `{uploaded_file.name}`")

        content = uploaded_file.read().decode("utf-8").strip().splitlines()

        records = []
        for i, line in enumerate(content):
            parts = line.strip().split()
            if len(parts) >= 2:
                seq_id = parts[0]
                sequence = "".join(parts[1:])
            else:
                seq_id = f"Seq_{i+1}"
                sequence = parts[0]
            records.append((seq_id, sequence))

        st.info(f"🔍 {len(records)} sequences detected")

        feature_rows = []
        for seq_id, sequence in records:
            try:
                features = extract_features(sequence)
                features["ID"] = seq_id
                feature_rows.append(features)
            except Exception as e:
                st.warning(f"⚠️ Error processing sequence {seq_id}: {e}")

        if feature_rows:
            df = pd.DataFrame(feature_rows)

            clf = RandomForestClassifier()
            dummy_data = df.drop(columns=["ID"])
            clf.fit(dummy_data, [0] * len(df))

            predictions = clf.predict(dummy_data)
            df["Predicted_Region"] = predictions

            st.subheader("🔬 Analysis Results")
            st.dataframe(df, use_container_width=True)

            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download Results as CSV", data=csv,
                               file_name="DNA_Analysis_Results.csv", mime="text/csv")
        else:
            st.warning("❗ No valid sequences found.")

    else:
        st.info("💡 Please upload a `.txt` file to begin.")

    st.markdown('</div>', unsafe_allow_html=True)

