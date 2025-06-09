import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier
import altair as alt
import pickle

# ------------------------- PAGE CONFIG --------------------------
st.set_page_config(
    page_title="Welcome to the DNA Motif & Perplexity Analyzer",
    layout="wide",
    page_icon="🧬"
)

# ------------------------- CSS Styling --------------------------
st.markdown("""
    <style>
    .main {
        background-color: #F0F2F6;
    }
    .css-1rs6os.edgvbvh3 {
        font-size: 18px;
    }
    .block-container {
        padding-top: 2rem;
    }
    .stDownloadButton {
        background-color: #4CAF50;
        color: white;
        padding: 0.5em 1em;
        border-radius: 8px;
        text-align: center;
    }
    </style>
""", unsafe_allow_html=True)

# ------------------------- LOGIN --------------------------
with st.sidebar:
    st.markdown("### 🔐 Login")
    username = st.text_input("Username")
    password = st.text_input("Password", type="password")
    login_button = st.button("Login")

if not (username == "admin" and password == "1234"):
    st.warning("Please enter valid login credentials to access the tool.")
    st.stop()

# ------------------------- TABS --------------------------
tabs = st.tabs(["🏠 Home", "📂 Analyze", "📊 Results", "📈 Insights", "💾 Save/Load", "📥 Download", "📞 Contact"])

# ------------------------- HELPER FUNCTIONS --------------------------

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

def extract_features(seq, k):
    seq = seq.upper()
    return {
        'Perplexity': calculate_perplexity(seq, k=k),
        'G-Quadruplex': detect_g_quadruplex(seq),
        'Z-DNA': detect_z_dna(seq),
        'Cruciform': detect_cruciform(seq),
        'TATA-Box': detect_tata_box(seq),
        'Direct Repeats': detect_direct_repeats(seq),
        'Sequence Length': len(seq)
    }

# ------------------------- HOME --------------------------
with tabs[0]:
    st.markdown("## 🧬 Welcome to the DNA Motif & Perplexity Analyzer")
    st.markdown("""
    This tool helps detect:
    - 🔁 Direct Repeats
    - 🧷 G-Quadruplexes
    - 🔀 Z-DNA
    - 🎯 TATA Boxes
    - 🧬 Cruciform Structures
    Use the Analyze tab to start!
    """)

# ------------------------- ANALYZE --------------------------
with tabs[1]:
    st.markdown("## 📂 Upload & Analyze DNA Sequences")
    k_value = st.slider("🔢 Choose k-mer size for Perplexity", min_value=2, max_value=6, value=3)
    uploaded_file = st.file_uploader("📄 Upload `.txt` or `.fasta` files", type=["txt", "fasta"])

    if uploaded_file:
        st.success(f"📁 File uploaded: `{uploaded_file.name}`")
        content = uploaded_file.read().decode("utf-8").strip().splitlines()
        records = []
        for line in content:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts = line.split()
            seq_id = parts[0]
            sequence = "".join(parts[1:]) if len(parts) > 1 else parts[0]
            records.append((seq_id, sequence))

        with st.spinner("🔍 Extracting features..."):
            feature_rows = []
            for seq_id, sequence in records:
                try:
                    features = extract_features(sequence, k_value)
                    features["ID"] = seq_id
                    feature_rows.append(features)
                except Exception as e:
                    st.error(f"Error processing {seq_id}: {e}")

            df = pd.DataFrame(feature_rows)
            df = df[["ID"] + [col for col in df.columns if col != "ID"]]

            clf = RandomForestClassifier()
            clf.fit(df.drop(columns=["ID"]), [0]*len(df))
            df["Prediction"] = clf.predict(df.drop(columns=["ID"]))

            st.session_state["results_df"] = df
            st.success("✅ Analysis complete! Check the 📊 Results tab.")

# ------------------------- RESULTS --------------------------
with tabs[2]:
    st.markdown("## 📊 Results Summary")
    if "results_df" in st.session_state:
        df = st.session_state["results_df"]
        selected_motif = st.selectbox("🔍 Filter by motif type", ["All"] + list(df.columns[1:-2]))
        if selected_motif != "All":
            df = df[df[selected_motif] > 0]

        def color_motifs(val):
            motif_icons = {
                "G-Quadruplex": "🧷", "Z-DNA": "🔀", "TATA-Box": "🎯",
                "Direct Repeats": "🔁", "Cruciform": "🧬"
            }
            return motif_icons.get(val, val)

        styled_df = df.style.applymap(lambda val: 'color: red' if isinstance(val, (int, float)) and val > 5 else '')
        st.dataframe(styled_df, use_container_width=True)
    else:
        st.warning("📂 Please upload and analyze sequences first.")

# ------------------------- INSIGHTS --------------------------
with tabs[3]:
    st.markdown("## 📈 Motif & Perplexity Insights")
    if "results_df" in st.session_state:
        df = st.session_state["results_df"]

        col1, col2 = st.columns(2)
        with col1:
            st.markdown("### 🧬 Motif Frequency Chart")
            melted = df.melt(id_vars=["ID", "Prediction"], value_vars=["G-Quadruplex", "Z-DNA", "TATA-Box", "Direct Repeats", "Cruciform"])
            chart = alt.Chart(melted).mark_bar().encode(
                x=alt.X("variable", title="Motif Type"),
                y=alt.Y("sum(value)", title="Total Count"),
                color="variable"
            ).properties(height=300)
            st.altair_chart(chart, use_container_width=True)

        with col2:
            st.markdown("### 📏 Sequence Length vs. Perplexity")
            scatter = alt.Chart(df).mark_circle(size=60).encode(
                x="Sequence Length",
                y="Perplexity",
                tooltip=["ID", "Perplexity", "Sequence Length"]
            ).interactive()
            st.altair_chart(scatter, use_container_width=True)
    else:
        st.warning("📂 Please upload and analyze sequences first.")

# ------------------------- SAVE/LOAD --------------------------
with tabs[4]:
    st.markdown("## 💾 Save or Load Session")
    if "results_df" in st.session_state:
        buffer = pickle.dumps(st.session_state["results_df"])
        st.download_button("💾 Save Analysis", data=buffer, file_name="session.pkl")
        uploaded_pickle = st.file_uploader("📤 Upload Saved Session (.pkl)", type=["pkl"])
        if uploaded_pickle:
            st.session_state["results_df"] = pickle.load(uploaded_pickle)
            st.success("🔄 Session restored! Check Results tab.")
    else:
        st.info("📂 Analyze data before saving or loading session.")

# ------------------------- DOWNLOAD --------------------------
with tabs[5]:
    st.markdown("## 📥 Download Your Results")
    if "results_df" in st.session_state:
        csv = st.session_state["results_df"].to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download as CSV", csv, file_name="motif_results.csv", mime="text/csv")
    else:
        st.warning("📂 No data available to download yet.")

# ------------------------- CONTACT --------------------------
with tabs[6]:
    st.markdown("## 📞 Contact Information")
    st.markdown("""
    - 👨‍🔬 **Dr. Y V Rajesh**  
      📧 yvrajesh_bt@kluniversity.in

    - 👩‍🔬 **G. Aruna Sesha Chandrika**  
      📧 chandrikagummadi1@gmail.com
    """)
