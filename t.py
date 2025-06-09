import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier
from email.message import EmailMessage
import json
import os
from datetime import datetime

# ----- File paths for persistence -----
USERS_FILE = "registered_users.json"
HISTORY_FILE = "user_history.json"

# ----- Helper functions to load/save JSON -----
def load_json(file_path):
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            return json.load(f)
    return {}

def save_json(file_path, data):
    with open(file_path, "w") as f:
        json.dump(data, f, indent=2)

# ----- Load registered users and histories -----
if "registered_users" not in st.session_state:
    st.session_state["registered_users"] = load_json(USERS_FILE)

if "user_history" not in st.session_state:
    st.session_state["user_history"] = load_json(HISTORY_FILE)

# ----- Page Config & CSS -----
st.set_page_config(page_title="DNA Motif & Perplexity Analyzer", layout="wide", page_icon="🧬")
st.markdown("""
<style>
    .main { background-color: #F0F2F6; }
    .css-1rs6os.edgvbvh3 { font-size: 18px; }
    .block-container { padding-top: 2rem; }
    .stDownloadButton {
        background-color: #4CAF50;
        color: white;
        padding: 0.5em 1em;
        border-radius: 8px;
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)

# ----- Registration/Login Sidebar -----
with st.sidebar:
    st.title("👤 Account")

    if "logged_in" not in st.session_state:
        st.session_state["logged_in"] = False

    if not st.session_state["logged_in"]:
        mode = st.radio("Select Mode", ["Register", "Login"])

        if mode == "Register":
            st.subheader("📝 Register")
            reg_username = st.text_input("Choose a Username", key="reg_user")
            reg_password = st.text_input("Choose a Password", type="password", key="reg_pass")
            reg_email = st.text_input("Email", key="reg_email")

            if st.button("Register"):
                if not reg_username or not reg_password or not reg_email:
                    st.warning("Please fill all fields")
                elif reg_username in st.session_state["registered_users"]:
                    st.info("This username is already registered. Try Login or pick another username.")
                else:
                    # Save new user
                    st.session_state["registered_users"][reg_username] = {
                        "password": reg_password,
                        "email": reg_email
                    }
                    save_json(USERS_FILE, st.session_state["registered_users"])

                    st.success(f"User '{reg_username}' registered successfully! Please switch to Login tab.")
                    st.experimental_rerun()

        elif mode == "Login":
            st.subheader("🔐 Login")
            username = st.text_input("Username", key="login_user")
            password = st.text_input("Password", type="password", key="login_pass")

            if st.button("Login"):
                user_data = st.session_state["registered_users"].get(username)
                if user_data and user_data["password"] == password:
                    st.session_state["logged_in"] = True
                    st.session_state["current_user"] = username
                    st.success(f"Welcome back, {username}!")
                    st.experimental_rerun()
                else:
                    st.info("Please check your username and password. If new, register first.")

        st.stop()

    else:
        st.write(f"👋 Logged in as: **{st.session_state['current_user']}**")
        if st.button("Logout"):
            st.session_state["logged_in"] = False
            st.session_state.pop("current_user", None)
            st.experimental_rerun()

# ----- Helper Functions for DNA analysis -----

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
    return
