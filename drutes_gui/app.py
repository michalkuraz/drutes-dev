"""Entry point for the DRUtES configuration GUI."""

import html
from pathlib import Path
import re
import subprocess
import sys
import threading

import streamlit as st

# Running ``streamlit run drutes_gui/app.py`` puts drutes_gui on sys.path.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pages.global_configuration import GlobalConfigurationPage  # noqa: E402
from pages.mesh_configuration import MeshConfigurationPage  # noqa: E402
from pages.model_configuration import (  # noqa: E402
    HeatConfigurationPage,
    MatrixConfigurationPage,
    RootUptakeConfigurationPage,
)
from pages.solver_configuration import SolverConfigurationPage  # noqa: E402
from pages.richards_outputs import RichardsOutputPage  # noqa: E402
from drutespy.config.global_config import GlobalConfigFile  # noqa: E402
from drutespy.config.heat_config import HeatConfigFile  # noqa: E402
from drutespy.config.matrix_config import MatrixConfigFile  # noqa: E402
from drutespy.config.mesh_config import Mesh1DConfigFile  # noqa: E402
from drutespy.config.solver_config import SolverConfigFile  # noqa: E402
from drutespy.config.root_uptake_config import RootUptakeConfigFile  # noqa: E402


def terminal_document(output: str) -> str:
    """Build an isolated, fixed-height terminal display for model output."""
    escaped_output = html.escape(output or "Waiting for terminal output…")
    return f"""
    <!doctype html>
    <html>
      <head>
        <style>
          html, body {{ margin: 0; background: transparent; }}
          .terminal {{
            background: #111820;
            border: 1px solid #303b48;
            border-radius: 8px;
            box-shadow: 0 10px 30px rgba(0, 0, 0, .16);
            color: #d8e2ec;
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas,
              "Liberation Mono", monospace;
            overflow: hidden;
          }}
          .titlebar {{
            align-items: center;
            background: #202a35;
            border-bottom: 1px solid #303b48;
            color: #aebccc;
            display: flex;
            font: 12px -apple-system, BlinkMacSystemFont, sans-serif;
            height: 38px;
            padding: 0 14px;
          }}
          .lights {{ display: flex; gap: 7px; margin-right: 14px; }}
          .light {{ border-radius: 50%; height: 11px; width: 11px; }}
          .red {{ background: #ff5f57; }}
          .yellow {{ background: #febc2e; }}
          .green {{ background: #28c840; }}
          #output {{
            box-sizing: border-box;
            height: 400px;
            margin: 0;
            overflow: auto;
            padding: 16px 18px;
            scrollbar-color: #536274 #111820;
            white-space: pre-wrap;
            word-break: break-word;
          }}
        </style>
      </head>
      <body>
        <div class="terminal">
          <div class="titlebar">
            <span class="lights">
              <span class="light red"></span>
              <span class="light yellow"></span>
              <span class="light green"></span>
            </span>
            DRUtES terminal output
          </div>
          <pre id="output">{escaped_output}</pre>
        </div>
        <script>
          const output = document.getElementById("output");
          output.scrollTop = output.scrollHeight;
        </script>
      </body>
    </html>
    """


ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def capture_model_output(
    process: subprocess.Popen[str], output_lines: list[str]
) -> None:
    """Continuously collect merged DRUtES output without blocking Streamlit."""
    if process.stdout is None:
        return
    for line in process.stdout:
        output_lines.append(ANSI_ESCAPE.sub("", line))
    process.stdout.close()


def stop_model_process(process: subprocess.Popen[str]) -> None:
    """Terminate a running model and escalate to kill if it does not stop."""
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=2)


@st.fragment(run_every=0.5)
def render_model_process_monitor() -> None:
    """Draw the live terminal and the transient kill control."""
    process = st.session_state.get("model_process")
    output_lines = st.session_state.get("model_output_lines", [])
    if process is None:
        output = str(st.session_state.get("model_terminal_output", ""))
    else:
        output = "".join(output_lines)

    st.iframe(terminal_document(output), height=460)
    if process is None:
        return

    return_code = process.poll()
    if return_code is None:
        st.info("DRUtES is running…")
        if st.button(
            "Kill simulation",
            type="primary",
            width="stretch",
            key="kill_model_process",
        ):
            st.session_state.model_was_killed = True
            stop_model_process(process)
            st.rerun(scope="fragment")
        return

    output_thread = st.session_state.get("model_output_thread")
    if output_thread is not None:
        output_thread.join(timeout=1)
    output = "".join(output_lines)
    st.session_state.model_terminal_output = output

    if st.session_state.get("model_was_killed", False):
        st.warning("DRUtES simulation was killed by the user.")
    elif return_code == 0:
        st.success("DRUtES finished successfully.")
    else:
        st.error(f"DRUtES exited with status {return_code}.")

    if not st.session_state.get("model_completion_handled", False):
        st.session_state.model_completion_handled = True
        st.rerun()


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

    if "configuration_page" not in st.session_state:
        st.session_state.configuration_page = "global"
    if "visited_configuration_pages" not in st.session_state:
        st.session_state.visited_configuration_pages = {"global"}

    def navigate(page: str) -> None:
        st.session_state.configuration_page = page
        st.session_state.visited_configuration_pages.add(page)

    def save_all(model_type: str) -> None:
        configs = [
            GlobalConfigFile(PROJECT_ROOT / "drutes.conf" / "global.conf"),
            Mesh1DConfigFile(
                PROJECT_ROOT / "drutes.conf" / "mesh" / "drumesh1d.conf"
            ),
            SolverConfigFile(PROJECT_ROOT / "drutes.conf" / "solver.conf"),
        ]
        if model_type == "heat":
            configs.append(
                HeatConfigFile(PROJECT_ROOT / "drutes.conf" / "heat" / "heat.conf")
            )
        else:
            matrix = MatrixConfigFile(
                PROJECT_ROOT / "drutes.conf" / "water.conf" / "matrix.conf"
            )
            configs.append(matrix)
            matrix.load()
            if bool(matrix["root_water_uptake"].value):
                mesh = Mesh1DConfigFile(
                    PROJECT_ROOT / "drutes.conf" / "mesh" / "drumesh1d.conf"
                ).load()
                configs.append(
                    RootUptakeConfigFile(
                        PROJECT_ROOT
                        / "drutes.conf"
                        / "water.conf"
                        / "root4uptake.conf",
                        int(mesh["layer_count"].value),
                    )
                )
        for config in configs:
            config.load().save()
        st.session_state.all_configurations_saved = True
        st.session_state.configuration_page = "run"
        st.rerun()

    def render_model_runner() -> None:
        if not st.session_state.get("all_configurations_saved", False):
            return

        st.divider()
        st.subheader("Model terminal")
        st.caption(
            "Runs bin/drutes from the repository root. Standard output and "
            "errors are displayed below."
        )
        process = st.session_state.get("model_process")
        running = process is not None and process.poll() is None
        render_model_process_monitor()

        if running or not st.button(
            "Run model", type="primary", use_container_width=True
        ):
            return

        executable = PROJECT_ROOT / "bin" / "drutes"
        if not executable.is_file():
            st.error(f"Model executable was not found: {executable}")
            return

        output_lines: list[str] = []
        try:
            process = subprocess.Popen(
                [str(executable)],
                cwd=PROJECT_ROOT,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError as error:
            st.error(f"Unable to run DRUtES: {error}")
            return

        output_thread = threading.Thread(
            target=capture_model_output,
            args=(process, output_lines),
            daemon=True,
            name="drutes-output-reader",
        )
        st.session_state.model_process = process
        st.session_state.model_output_lines = output_lines
        st.session_state.model_output_thread = output_thread
        st.session_state.model_terminal_output = ""
        st.session_state.model_was_killed = False
        st.session_state.model_completion_handled = False
        output_thread.start()
        st.rerun()

    current_page = st.session_state.configuration_page
    if current_page != "global" and logo_path.exists():
        st.image(str(logo_path), width=420)

    if current_page == "mesh":
        mesh_page = MeshConfigurationPage(
            PROJECT_ROOT / "drutes.conf" / "mesh" / "drumesh1d.conf",
            on_return=lambda: navigate("global"),
        )
        mesh_page.on_navigate = navigate
        mesh_page.model_type = str(
            st.session_state.get("selected_problem_type", "RE")
        )
        mesh_page.visited_pages = st.session_state.visited_configuration_pages
        mesh_page.render()
    elif current_page == "solver":
        solver_page = SolverConfigurationPage(
            PROJECT_ROOT / "drutes.conf" / "solver.conf",
            on_return=lambda: navigate("global"),
        )
        solver_page.on_navigate = navigate
        solver_page.model_type = str(
            st.session_state.get("selected_problem_type", "RE")
        )
        solver_page.visited_pages = st.session_state.visited_configuration_pages
        solver_page.render()
    elif current_page == "heat":
        heat_page = HeatConfigurationPage(
            PROJECT_ROOT / "drutes.conf" / "heat" / "heat.conf",
            on_return=lambda: navigate("global"),
        )
        heat_page.on_navigate = navigate
        heat_page.on_save_all = save_all
        heat_page.model_type = "heat"
        heat_page.visited_pages = st.session_state.visited_configuration_pages
        heat_page.render()
    elif current_page == "matrix":
        matrix_page = MatrixConfigurationPage(
            PROJECT_ROOT / "drutes.conf" / "water.conf" / "matrix.conf",
            on_return=lambda: navigate("global"),
        )
        matrix_page.on_navigate = navigate
        matrix_page.on_save_all = save_all
        matrix_page.model_type = "RE"
        matrix_page.visited_pages = st.session_state.visited_configuration_pages
        matrix_page.render()
    elif current_page == "root_uptake":
        root_page = RootUptakeConfigurationPage(
            PROJECT_ROOT / "drutes.conf" / "water.conf" / "root4uptake.conf",
            on_return=lambda: navigate("global"),
        )
        root_page.on_navigate = navigate
        root_page.on_save_all = save_all
        root_page.model_type = "RE"
        root_page.visited_pages = st.session_state.visited_configuration_pages
        root_page.render()
    elif current_page == "run":
        st.markdown(
            '<p class="drutes-kicker">Simulation execution</p>',
            unsafe_allow_html=True,
        )
        st.title("Run DRUtES model")
        st.markdown(
            '<p class="drutes-subtitle">Review or re-edit configuration files, '
            'then run the configured model.</p><div class="drutes-rule"></div>',
            unsafe_allow_html=True,
        )
        st.subheader("Re-edit configurations")
        active_process = st.session_state.get("model_process")
        model_running = (
            active_process is not None and active_process.poll() is None
        )
        selected_model = str(
            st.session_state.get("selected_problem_type", "RE")
        )
        edit_pages = [
            ("global", "global.conf"),
            ("mesh", "drumesh1D.conf"),
            ("solver", "solver.conf"),
        ]
        if selected_model == "heat":
            edit_pages.append(("heat", "heat.conf"))
        else:
            edit_pages.append(("matrix", "matrix.conf"))
            if st.session_state.get("root_uptake_enabled", False):
                edit_pages.append(("root_uptake", "root4uptake.conf"))

        edit_columns = st.columns(len(edit_pages))
        for column, (page, filename) in zip(
            edit_columns, edit_pages, strict=True
        ):
            with column:
                if st.button(
                    f"Re-edit {filename}",
                    key=f"run_reedit_{page}",
                    use_container_width=True,
                    disabled=model_running,
                ):
                    navigate(page)
                    st.rerun()
        render_model_runner()
        if selected_model != "heat" and not model_running:
            RichardsOutputPage(
                PROJECT_ROOT / "out",
                PROJECT_ROOT / "drutes.conf" / "global.conf",
            ).render()
    else:
        GlobalConfigurationPage(
            PROJECT_ROOT / "drutes.conf" / "global.conf",
            logo_path,
            on_edit_mesh=lambda: navigate("mesh"),
            on_edit_solver=lambda: navigate("solver"),
            on_edit_model=navigate,
            on_save_all=save_all,
            visited_pages=st.session_state.visited_configuration_pages,
        ).render()


if __name__ == "__main__":
    main()
