import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier

# ----------------------- Page Configuration -----------------------
st.set_page_config(page_title="DNA Analyzer", layout="wide")

# ----------------------- Custom CSS Styling -----------------------
st.markdown("""
    <style>
    body {
        background: linear-gradient(to bottom right, #ffecd2, #fcb69f);
    }
    .main {
        background-color: #fff5f0;
    }
    .title {
        color: #6a1b9a;
        text-align: center;
        font-size: 3em;
        font-weight: bold;
        margin-bottom: 20px;
    }
    .section {
        background-color: #ffffffdd;
        padding: 20px;
        border-radius: 15px;
        box-shadow: 2px 2px 10px rgba(0,0,0,0.1);
        margin: 20px 0;
    }
    .stButton>button {
        background-color: #43a047;
        color: white;
        font-weight: bold;
        border-radius: 10px;
        padding: 0.6em 1em;
        border: none;
    }
    .stDownloadButton>button {
        background-color: #e65100;
        color: white;
        font-weight: bold;
        border-radius: 10px;
    }
    </style>
""", unsafe_allow_html=True)

st.markdown('<div class="title">🧬 Colorful DNA Motif & Perplexity Analyzer</div>', unsafe_allow_html=True)

# ----------------------- Feature Functions -----------------------
def calculate_perplexity(sequence, k=3):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    kmer_counts = pd.Series(kmers).value_counts()
    total = sum(kmer_counts)
    probs = kmer_counts / total
    entropy = -sum(p * math.log2(p) for p in probs)
    perplexity = 2 ** entropy
    return perplexity

def detect_g_quadruplex(seq):
    return len(re.findall(r'(G{3,}\w{1,7}){3,}G{3,}', seq))

def detect_z_dna(seq):
    return len(re.findall(r'(CG){6,}', seq))

def detect_tata_box(seq):
    return len(re.findall(r'TATA[AT]A[AT]', seq))

def detect_direct_repeats(seq):
    return len(re.findall(r'(.{3,6})\1+', seq))

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def detect_cruciform(seq, min_len=4, max_len=6, spacer=10):
    count = 0
    for size in range(min_len, max_len + 1):
        for i in range(len(seq) - 2 * size - spacer + 1):
            left = seq[i:i + size]
            right = seq[i + size + spacer:i + 2 * size + spacer]
            if reverse_complement(left) == right:
                count += 1
    return count

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

# ----------------------- App Layout -----------------------

with st.container():
    st.markdown('<div class="section">', unsafe_allow_html=True)
    
    uploaded_file = st.file_uploader("📄 Upload a `.txt` file with DNA sequences (1 per line or with IDs)", type=["txt"])

    if uploaded_file:
        st.success(f"✅ File uploaded: `{uploaded_file.name}`")

        content = uploaded_file.read().decode("utf-8").strip().splitlines()

        records = []
        for i, line in enumerate(content):
            parts = line.strip().split()
            if len(parts) >= 2:
                seq_id, sequence = parts[0], ''.join(parts[1:])
            else:
                seq_id, sequence = f"Seq_{i+1}", parts[0]
            records.append((seq_id, sequence))

        feature_rows = []
        for seq_id, sequence in records:
            try:
                features = extract_features(sequence)
                features["ID"] = seq_id
                feature_rows.append(features)
            except Exception as e:
                st.warning(f"⚠️ Error in {seq_id}: {e}")

        if feature_rows:
            df = pd.DataFrame(feature_rows)
            clf = RandomForestClassifier()
            dummy_X = df.drop(columns=["ID"])
            clf.fit(dummy_X, [0] * len(df))
            df["Predicted_Region"] = clf.predict(dummy_X)

            st.subheader("📊 Motif & Perplexity Results")
            st.dataframe(df, use_container_width=True)

            st.download_button("📥 Download CSV", data=df.to_csv(index=False).encode("utf-8"),
                               file_name="dna_results.csv", mime="text/csv")
        else:
            st.warning("❗ No valid sequences found.")

    else:
        st.info("📝 Please upload a .txt file with DNA sequences to start.")
    
    st.markdown('</div>', unsafe_allow_html=True)
