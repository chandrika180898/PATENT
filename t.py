import streamlit as st
import pandas as pd
import re
import math
import hashlib
import os
from sklearn.ensemble import RandomForestClassifier

# ------------------------- PAGE CONFIG --------------------------
st.set_page_config(
    page_title="DNA Motif Analyzer",
    layout="wide",
    page_icon="🧬"
)

# ------------------------- CSS Styling --------------------------
st.markdown("""
    <style>
    .main {
        background-color: #F0F2F6;
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

# ------------------------- AUTHENTICATION --------------------------
USER_FILE = "users.csv"

def hash_password(password):
    return hashlib.sha256(password.encode()).hexdigest()

def load_users():
    if not os.path.exists(USER_FILE):
        return {}
    df = pd.read_csv(USER_FILE)
    return dict(zip(df.username, df.password))

def save_user(username, password):
    hashed = hash_password(password)
    if os.path.exists(USER_FILE):
        df = pd.read_csv(USER_FILE)
        df = pd.concat([df, pd.DataFrame({"username": [username], "password": [hashed]})])
    else:
        df = pd.DataFrame({"username": [username], "password": [hashed]})
    df.to_csv(USER_FILE, index=False)

def authenticate(username, password, users):
    return users.get(username) == hash_password(password)

# ------------------------- SIDEBAR NAVIGATION --------------------------
if "logged_in" not in st.session_state:
    st.session_state.logged_in = False

if not st.session_state.logged_in:
    page = st.sidebar.radio("Choose Page", ["🔐 Login", "📝 Register"])
else:
    page = st.sidebar.radio("Choose Page", ["🏠 Home", "📂 Upload & Analyze", "📊 Results", "📥 Download Report", "📞 Contact", "🚪 Logout"])

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

def extract_features(seq):
    seq = seq.upper()
    return {
        'Perplexity': calculate_perplexity(seq),
        'G-Quadruplex': detect_g_quadruplex(seq),
        'Z-DNA': detect_z_dna(seq),
        'Cruciform': detect_cruciform(seq),
        'TATA-Box': detect_tata_box(seq),
        'Direct Repeats': detect_direct_repeats(seq),
        'Sequence Length': len(seq)
    }

# ------------------------- PAGES --------------------------

if page == "🔐 Login":
    st.title("🔐 Login")
    username = st.text_input("Username")
    password = st.text_input("Password", type="password")
    if st.button("Login"):
        users = load_users()
        if authenticate(username, password, users):
            st.success("✅ Logged in successfully!")
            st.session_state.logged_in = True
            st.session_state.username = username
            st.rerun()
        else:
            st.error("❌ Invalid credentials.")

elif page == "📝 Register":
    st.title("📝 Register")
    username = st.text_input("Choose a Username")
    password = st.text_input("Choose a Password", type="password")
    confirm = st.text_input("Confirm Password", type="password")
    if st.button("Register"):
        if password != confirm:
            st.error("❌ Passwords do not match.")
        else:
            users = load_users()
            if username in users:
                st.error("❌ Username already exists.")
            else:
                save_user(username, password)
                st.success("✅ Registered successfully! You can now login.")

elif page == "🚪 Logout":
    st.session_state.logged_in = False
    st.success("👋 You have been logged out.")
    st.rerun()

elif page == "🏠 Home":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()
    st.markdown("## 🧬 Welcome to the **DNA Motif & Perplexity Analyzer**")
    st.markdown("""
    This web-based tool is built for analyzing DNA sequences and detecting **non-B DNA structures**, 
    calculating **perplexity**, and identifying key **regulatory motifs** such as:
    
    - 🔁 **Direct Repeats**
    - 🧷 **G-Quadruplexes**
    - 🔀 **Z-DNA**
    - 🎯 **TATA Boxes**
    - 🧬 **Cruciform Structures**

    Go to **📂 Upload & Analyze** to begin!
    """)

elif page == "📂 Upload & Analyze":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()

    st.markdown("## 📂 Upload & Analyze DNA Sequences")
    uploaded_file = st.file_uploader("📄 Upload a `.txt` file with sequences (one per line)", type=["txt"])

    if uploaded_file:
        st.success(f"📁 File uploaded: `{uploaded_file.name}`")

        content = uploaded_file.read().decode("utf-8").strip().splitlines()
        records = []
        for i, line in enumerate(content):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            seq_id = parts[0]
            sequence = "".join(parts[1:]) if len(parts) > 1 else parts[0]
            records.append((seq_id, sequence))

        with st.spinner("🔍 Extracting features..."):
            feature_rows = []
            for seq_id, sequence in records:
                try:
                    features = extract_features(sequence)
                    features["ID"] = seq_id
                    feature_rows.append(features)
                except Exception as e:
                    st.error(f"Error processing {seq_id}: {e}")

            df = pd.DataFrame(feature_rows)
            df = df[["ID"] + [col for col in df.columns if col != "ID"]]  # Reorder

            # Dummy ML prediction
            clf = RandomForestClassifier()
            clf.fit(df.drop(columns=["ID"]), [0]*len(df))
            df["Prediction"] = clf.predict(df.drop(columns=["ID"]))

            st.session_state["results_df"] = df
            st.success("✅ Analysis complete! Check the 📊 Results tab.")

elif page == "📊 Results":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()
    st.markdown("## 📊 Results Summary")
    if "results_df" in st.session_state:
        st.dataframe(st.session_state["results_df"], use_container_width=True)
    else:
        st.warning("📂 Please upload and analyze sequences first.")

elif page == "📥 Download Report":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()
    st.markdown("## 📥 Download Your Results")
    if "results_df" in st.session_state:
        csv = st.session_state["results_df"].to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download as CSV", csv, file_name="motif_results.csv", mime="text/csv")
    else:
        st.warning("📂 No data available to download yet.")

elif page == "📞 Contact":
    st.markdown("## 📞 Contact Information")
    st.markdown("""
    👨‍🔬 **Dr. Y V Rajesh**  
      📧 yvrajesh_bt@kluniversity.in

    - 👩‍🔬 **G. Aruna Sesha Chandrika**  
      📧 chandrikagummadi1@gmail.com
    """)
