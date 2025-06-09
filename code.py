import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier

# ---------- Page Setup ----------
st.set_page_config(page_title="Colorful DNA Analyzer", layout="wide")

# ---------- Full Color CSS ----------
st.markdown("""
    <style>
    body {
        background: linear-gradient(135deg, #f6d365 0%, #fda085 100%);
        font-family: 'Segoe UI', sans-serif;
    }
    .main {
        background-color: rgba(255, 255, 255, 0.9);
        padding: 2rem;
        border-radius: 12px;
        box-shadow: 0 0 25px rgba(0,0,0,0.1);
        margin-top: 2rem;
    }
    h1, h2, h3 {
        color: #4a148c;
        text-align: center;
    }
    .stTextInput>div>div>input {
        background-color: #fff8e1;
        color: #000;
    }
    .stDownloadButton>button, .stButton>button {
        background-color: #00897b;
        color: white;
        border-radius: 12px;
        font-weight: bold;
        transition: 0.3s;
    }
    .stDownloadButton>button:hover, .stButton>button:hover {
        background-color: #00695c;
    }
    .css-1cpxqw2, .css-ffhzg2 {
        background-color: #ffffffaa !important;
        padding: 10px;
        border-radius: 10px;
    }
    </style>
""", unsafe_allow_html=True)

# ---------- Helper Functions ----------
def calculate_perplexity(sequence, k=3):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    kmer_counts = pd.Series(kmers).value_counts()
    total = sum(kmer_counts)
    probs = kmer_counts / total
    entropy = -sum(p * math.log2(p) for p in probs)
    return 2 ** entropy

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

# ---------- App UI ----------
st.markdown('<div class="main">', unsafe_allow_html=True)

st.title("🧬 Colorful DNA Motif & Perplexity Analyzer")
uploaded_file = st.file_uploader("📄 Upload your .txt file (1 DNA sequence per line or ID + sequence)", type=["txt"])

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

        st.download_button("📥 Download Results CSV", data=df.to_csv(index=False).encode("utf-8"),
                           file_name="dna_results.csv", mime="text/csv")
    else:
        st.warning("❗ No valid sequences found.")

else:
    st.info("📝 Upload a `.txt` file containing DNA sequences.")

st.markdown('</div>', unsafe_allow_html=True)
