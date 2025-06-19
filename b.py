import streamlit as st
import pandas as pd
import re
import math
import hashlib
import os
from sklearn.ensemble import RandomForestClassifier
from datetime import datetime

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
HISTORY_DIR = "history"

def hash_password(password):
    return hashlib.sha256(password.encode()).hexdigest()

def load_users():
    if not os.path.exists(USER_FILE):
        return {}
    df = pd.read_csv(USER_FILE)
    return dict(zip(df.username, df.password))

def save_user(username, password):
    hashed = hash_password(password)
    df = pd.DataFrame({"username": [username], "password": [hashed]})
    if os.path.exists(USER_FILE):
        df_existing = pd.read_csv(USER_FILE)
        df = pd.concat([df_existing, df], ignore_index=True)
    df.to_csv(USER_FILE, index=False)

def authenticate(username, password, users):
    return users.get(username) == hash_password(password)

def log_user_history(username, filename, df):
    os.makedirs(HISTORY_DIR, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    history_file = os.path.join(HISTORY_DIR, f"{username}_history.csv")
    df_copy = df.copy()
    df_copy.insert(0, "Filename", filename)
    df_copy.insert(1, "Timestamp", timestamp)

    if os.path.exists(history_file):
        df_existing = pd.read_csv(history_file)
        df_copy = pd.concat([df_existing, df_copy], ignore_index=True)

    df_copy.to_csv(history_file, index=False)

def load_user_history(username):
    history_file = os.path.join(HISTORY_DIR, f"{username}_history.csv")
    if os.path.exists(history_file):
        return pd.read_csv(history_file)
    else:
        return pd.DataFrame()

# ------------------------- SIDEBAR NAVIGATION --------------------------
if "logged_in" not in st.session_state:
    st.session_state.logged_in = False

if not st.session_state.logged_in:
    page = st.sidebar.radio("Choose Page", ["🔐 Login", "📝 Register"])
else:
    page = st.sidebar.radio("Choose Page", ["🏠 Home", "📂 Upload & Analyze", "📊 Results", "📅 Download Report", "📂 History", "📞 Contact", "🚪 Logout"])

# ------------------------- FEATURE EXTRACTION --------------------------
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
    Analyze DNA sequences for **non-B DNA structures**, calculate **perplexity**, and identify motifs such as:

    - 🔁 **Direct Repeats**
    - 🧷 **G-Quadruplexes**
    - 🔀 **Z-DNA**
    - 🎯 **TATA Boxes**
    - 🧬 **Cruciform Structures**
    """)

elif page == "📂 Upload & Analyze":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()

    st.markdown("## 📂 Upload & Analyze DNA Sequences")
    uploaded_file = st.file_uploader("📄 Upload a `.txt`, `.fna`, or `.fasta` file with sequences", type=["txt", "fna", "fasta"])

    if uploaded_file:
        st.success(f"📁 File uploaded: `{uploaded_file.name}`")
        content = uploaded_file.read().decode("utf-8").strip().splitlines()
        records = []
        current_seq = []

        for line in content:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    records.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line)

        if current_seq:
            records.append("".join(current_seq))

        st.info(f"🔍 Found {len(records)} sequences")

        with st.spinner("🔬 Analyzing sequences..."):
            feature_rows = []
            for sequence in records:
                try:
                    features = extract_features(sequence)
                    feature_rows.append(features)
                except Exception as e:
                    st.error(f"Error processing a sequence: {e}")

            df = pd.DataFrame(feature_rows)

            clf = RandomForestClassifier()
            clf.fit(df, [0]*len(df))
            df["Prediction"] = clf.predict(df)

            st.session_state["results_df"] = df

            log_user_history(st.session_state.username, uploaded_file.name, df)

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

elif page == "📅 Download Report":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()
    st.markdown("## 📅 Download Your Results")
    if "results_df" in st.session_state:
        csv = st.session_state["results_df"].to_csv(index=False).encode("utf-8")
        st.download_button("📅 Download as CSV", csv, file_name="motif_results.csv", mime="text/csv")
    else:
        st.warning("📂 No data available to download yet.")

elif page == "📂 History":
    if not st.session_state.get("logged_in", False):
        st.warning("🔒 Please log in to access this page.")
        st.stop()
    st.markdown("## 📂 Past Upload History")
    df_hist = load_user_history(st.session_state.username)
    if df_hist.empty:
        st.info("🕒 No history available yet.")
    else:
        st.dataframe(df_hist, use_container_width=True)

elif page == "📞 Contact":
    st.markdown("## 📞 Contact Information")
    st.markdown("""
    👨‍🔬 **Dr. Y V Rajesh**  
      📧 yvrajesh_bt@kluniversity.in

    👩‍🔬 **G. Aruna Sesha Chandrika**  
      📧 chandrikagummadi1@gmail.com
    """)
