import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier

# ----------------------- Page Configuration -----------------------
st.set_page_config(page_title="DNA Motif & Perplexity Analyzer", layout="wide")

# ----------------------- Custom Styling -----------------------
st.markdown("""
    <style>
    body {
        background-color: #f0f9ff;
        font-family: 'Segoe UI', sans-serif;
    }
    .main-box {
        background-color: #ffffff;
        padding: 2rem;
        border-radius: 15px;
        box-shadow: 0px 4px 10px rgba(0,0,0,0.1);
        margin-top: 20px;
    }
    h1, .stFileUploaderLabel {
        color: #004d40;
        text-align: center;
    }
    .stButton>button, .stDownloadButton>button {
        background-color: #00796b;
        color: white;
        font-weight: bold;
        border-radius: 10px;
        padding: 8px 20px;
        margin-top: 10px;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover, .stDownloadButton>button:hover {
        background-color: #004d40;
    }
    .stDataFrame {
        background-color: #e0f2f1;
        border-radius: 10px;
        padding: 10px;
    }
    </style>
""", unsafe_allow_html=True)

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
            left = seq[i:i+size]
            right = seq[i+size+spacer:i+2*size+spacer]
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

# ----------------------- Streamlit App UI -----------------------
st.markdown('<div class="main-box">', unsafe_allow_html=True)

st.title("📄 DNA Motif & Perplexity Analyzer from TXT File")

uploaded_file = st.file_uploader("📂 Upload a .txt file with DNA sequences", type=["txt"])

if uploaded_file:
    st.success(f"✅ Uploaded: {uploaded_file.name}")

    try:
        content = uploaded_file.read().decode("utf-8").strip().splitlines()
    except Exception as e:
        st.error(f"Error reading file: {e}")
        st.stop()

    records = []
    for i, line in enumerate(content):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            seq_id = parts[0]
            sequence = "".join(parts[1:])
        else:
            seq_id = f"Seq_{i+1}"
            sequence = parts[0]
        records.append((seq_id, sequence))

    st.info(f"🔍 Found {len(records)} sequences")

    feature_rows = []
    for seq_id, sequence in records:
        try:
            features = extract_features(sequence)
            features["ID"] = seq_id
            feature_rows.append(features)
        except Exception as e:
            st.error(f"❌ Error in sequence {seq_id}: {e}")

    if feature_rows:
        df = pd.DataFrame(feature_rows)

        # Dummy Classifier
        clf = RandomForestClassifier()
        dummy_data = df.drop(columns=["ID"])
        clf.fit(dummy_data, [0] * len(df))
        df["Predicted_Region"] = clf.predict(dummy_data)

        st.subheader("🔬 Analysis Results")
        st.dataframe(df, use_container_width=True)

        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download Results as CSV", data=csv, file_name="txt_sequence_results.csv", mime="text/csv")
    else:
        st.warning("⚠️ No valid sequence records found.")

else:
    st.info("📌 Please upload a `.txt` file to start the analysis.")

st.markdown('</div>', unsafe_allow_html=True)
