"""Entry point for the DRUtES configuration GUI."""

from pathlib import Path
import sys

import streamlit as st

# Running ``streamlit run drutes_gui/app.py`` puts drutes_gui on sys.path.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pages.global_configuration import GlobalConfigurationPage  # noqa: E402


def apply_drutes_theme() -> None:
    """Apply a Streamlit-compatible interpretation of the drutes.org theme."""
    st.markdown(
        """
        <style>
        :root {
            --drutes-blue: #0a84fe;
            --drutes-green: #85ca7b;
            --drutes-text: #555555;
            --drutes-dark: #444444;
            --drutes-muted: #8c8c8c;
            --drutes-border: #e5e5e5;
        }
        .stApp { background: #fafafa; }
        [data-testid="stHeader"] {
            background: rgba(255, 255, 255, .96);
            border-bottom: 1px solid var(--drutes-border);
        }
        .block-container {
            max-width: 1180px;
            padding-top: 2.2rem;
            padding-bottom: 4rem;
        }
        h1, h2, h3 { color: var(--drutes-dark); }
        h1 {
            font-size: clamp(2rem, 4vw, 3rem) !important;
            letter-spacing: -.035em;
            margin-bottom: .25rem !important;
        }
        h2, h3 {
            letter-spacing: -.02em;
        }
        div[data-testid="stVerticalBlockBorderWrapper"] {
            background: #fff;
            border: 1px solid var(--drutes-border);
            border-radius: 4px;
            box-shadow: 0 8px 30px rgba(25, 25, 25, .045);
        }
        div[data-testid="stButton"] button[kind="primary"] {
            min-height: 3.25rem;
            border-radius: 3px;
            font-weight: 700;
            letter-spacing: .01em;
            box-shadow: 0 8px 20px rgba(10, 132, 254, .18);
        }
        div[data-testid="stButton"] button[kind="primary"]:hover {
            background: #0073e6;
            border-color: #0073e6;
        }
        div[data-baseweb="select"] > div,
        div[data-testid="stNumberInputContainer"],
        div[data-testid="stTextInputRootElement"] {
            border-radius: 3px;
        }
        .drutes-kicker {
            color: var(--drutes-blue);
            font-size: .78rem;
            font-weight: 800;
            letter-spacing: .14em;
            margin: 1.5rem 0 .35rem;
            text-transform: uppercase;
        }
        .drutes-subtitle {
            color: var(--drutes-muted);
            font-size: 1.08rem;
            margin-bottom: 1.8rem;
        }
        .drutes-rule {
            background: linear-gradient(90deg, var(--drutes-blue), var(--drutes-green));
            border-radius: 2px;
            height: 4px;
            margin: .8rem 0 2rem;
            width: 88px;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def main() -> None:
    logo_path = PROJECT_ROOT / "logo.png"
    st.set_page_config(
        page_title="DRUtES GUI", page_icon=str(logo_path), layout="wide"
    )
    apply_drutes_theme()
    GlobalConfigurationPage(
        PROJECT_ROOT / "drutes.conf" / "global.conf", logo_path
    ).render()


if __name__ == "__main__":
    main()
