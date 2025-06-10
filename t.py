import subprocess
import sys

# Auto-install streamlit-authenticator if not installed
try:
    import streamlit_authenticator
except ModuleNotFoundError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "streamlit-authenticator"])
    import streamlit_authenticator

import streamlit as st
import pandas as pd
import re
import math
import smtplib
import shelve
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from sklearn.ensemble import RandomForestClassifier
from streamlit_authenticator import Authenticate

# ------------------------- CONFIG --------------------------
st.set_page_config(page_title="DNA Motif Analyzer", layout="wide", page_icon="🧬")
st.markdown("""<style>.main {background-color: #F0F2F6;}</style>""", unsafe_allow_html=True)

# ------------------------- AUTH SETUP --------------------------
with shelve.open("user_data") as db:
    if "credentials" not in db:
        db["credentials"] = {
            "usernames": {
                "admin": {
                    "name": "Admin",
                    "password": "admin123",  # Use hashed passwords in production
                    "email": "admin@example.com"
                }
            }
        }

credentials = shelve.open("user_data")["credentials"]
authenticator = Authenticate(credentials, "motif_app", "abcdef", cookie_expiry_days=1)

# ✅ Use positional "main" argument to avoid login type error
name, authentication_status, username = authenticator.login("Login", location="main")


# ------------------------- EMAIL FUNCTION --------------------------
def send_email(to_email, subject, content):
    from_email = "your_email@example.com"
    password = "your_app_password"
    msg = MIMEMultipart()
    msg["From"] = from_email
    msg["To"] = to_email
    msg["Subject"] = subject
    msg.attach(MIMEText(content, "plain"))
    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(from_email, password)
        server.send_message(msg)
        server.quit()
    except Exception as e:
        st.error(f"Email error: {e}")

# ------------------------- REGISTRATION --------------------------
if not authentication_status:
    with st.expander("Register Here"):
        new_user = st.text_input("Username")
        new_pass = st.text_input("Password", type="password")
        new_email = st.text_input("Email")
        if st.button("Register"):
            if new_user and new_pass and new_email:
                credentials["usernames"][new_user] = {
                    "name": new_user,
                    "password": new_pass,
                    "email": new_email
                }
                with shelve.open("user_data") as db:
                    db["credentials"] = credentials
                send_email(new_email, "DNA Motif Analyzer Registration", f"Welcome {new_user}! You’ve registered successfully.")
                st.success("✅ Registered successfully. Please log in.")
            else:
                st.warning("All fields are required.")

# ------------------------- HELPER FUNCTIONS --------------------------
def calculate_perplexity(sequence, k=3):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    kmer_counts = pd.Series(kmers).value_counts()
    total = sum(kmer_counts)
    probs = kmer_counts / total
    entropy = -sum(p * math.log2(p) for p in probs)
    return 2 ** entropy

def detect_g_quadruplex(seq): return len(re.findall(r'(G{3,}\w{1,7}){3,}G{3,}', seq))
def detect_z_dna(seq): return len(re.findall(r'(CG){6,}', seq))

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

def detect_tata_box(seq): return len(re.findall(r'TATA[AT]A[AT]', seq))
def detect_direct_repeats(seq): return len(re.findall(r'(.{3,6})\1+', seq))

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

# ------------------------- MAIN APP --------------------------
if authentication_status:
    st.sidebar.success(f"Welcome {name}!")
    page = st.sidebar.radio("📁 Navigation", ["🏠 Home", "📂 Upload & Analyze", "📊 Results", "📜 History", "📞 Contact"])

    if page == "🏠 Home":
        st.title("🧬 DNA Motif & Perplexity Analyzer")
        st.markdown("""
        Analyze DNA sequences to detect non-B DNA structures and regulatory motifs:
        - 🔁 Direct Repeats
        - 🧷 G-Quadruplexes
        - 🔀 Z-DNA
        - 🎯 TATA Boxes
        - 🧬 Cruciforms
        """)

    elif page == "📂 Upload & Analyze":
        uploaded_file = st.file_uploader("Upload `.txt` file with sequences (ID + sequence)", type=["txt"])
        if uploaded_file:
            content = uploaded_file.read().decode("utf-8").strip().splitlines()
            records = []
            for line in content:
                if line.strip():
                    parts = line.split()
                    seq_id = parts[0]
                    sequence = "".join(parts[1:]) if len(parts) > 1 else parts[0]
                    records.append((seq_id, sequence))

            feature_rows = []
            for seq_id, sequence in records:
                features = extract_features(sequence)
                features["ID"] = seq_id
                feature_rows.append(features)

            df = pd.DataFrame(feature_rows)
            clf = RandomForestClassifier()
            clf.fit(df.drop(columns=["ID"]), [0]*len(df))
            df["Prediction"] = clf.predict(df.drop(columns=["ID"]))

            # Store history
            with shelve.open("user_history") as db:
                if username not in db:
                    db[username] = []
                db[username].append(df.to_dict())

            st.session_state["results_df"] = df
            st.success("✅ Analysis complete. See results tab.")

    elif page == "📊 Results":
        if "results_df" in st.session_state:
            df = st.session_state["results_df"]
            st.dataframe(df, use_container_width=True)
            st.markdown("### 🔢 Motif Count Summary:")
            st.write(df.drop(columns=["ID", "Perplexity", "Sequence Length", "Prediction"]).sum())
        else:
            st.warning("No data available. Upload sequences first.")

    elif page == "📜 History":
        with shelve.open("user_history") as db:
            history = db.get(username, [])
        if history:
            st.markdown("### 📂 Your Analysis History")
            for idx, past in enumerate(history):
                st.markdown(f"#### 🧾 Analysis #{idx + 1}")
                st.dataframe(pd.DataFrame(past))
        else:
            st.warning("No history found.")

    elif page == "📞 Contact":
        st.markdown("## 📞 Contact Information")
        st.markdown("""
        - 👨‍🔬 **Dr. Y V Rajesh**  
          📧 yvrajesh_bt@kluniversity.in  
        - 👩‍🔬 **G. Aruna Sesha Chandrika**  
          📧 chandrikagummadi1@gmail.com
        """)
else:
    st.warning("Please login or register to continue.")
