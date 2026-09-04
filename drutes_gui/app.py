"""Entry point for the DRUtES configuration GUI."""

import html
from io import BytesIO
from pathlib import Path
import re
import subprocess
import sys
import threading
from zipfile import ZIP_DEFLATED, ZipFile

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
from pages.heat_outputs import HeatOutputPage  # noqa: E402
from pages.solver_time_output import SolverTimeOutputPage  # noqa: E402
from pages.simulation_log import SimulationLogPage  # noqa: E402
from drutespy.config.global_config import GlobalConfigFile  # noqa: E402
from drutespy.config.heat_config import HeatConfigFile  # noqa: E402
from drutespy.config.matrix_config import MatrixConfigFile  # noqa: E402
from drutespy.config.mesh_config import Mesh1DConfigFile  # noqa: E402
from drutespy.config.solver_config import SolverConfigFile  # noqa: E402
from drutespy.config.root_uptake_config import RootUptakeConfigFile  # noqa: E402
from drutes_gui.projects import (  # noqa: E402
    ProjectStore,
    build_configuration_archive,
)


def build_output_archive(output_directory: Path) -> bytes:
    """Return the complete output directory as an in-memory ZIP archive."""
    archive_buffer = BytesIO()
    with ZipFile(archive_buffer, mode="w", compression=ZIP_DEFLATED) as archive:
        archive.writestr(f"{output_directory.name}/", b"")
        for path in sorted(output_directory.rglob("*")):
            if not path.is_file() or path.is_symlink():
                continue
            archive.write(
                path,
                arcname=str(Path(output_directory.name) / path.relative_to(output_directory)),
            )
    return archive_buffer.getvalue()


def authentication_is_configured() -> bool:
    """Return whether all required Google OIDC settings are present."""
    try:
        auth = st.secrets["auth"]
    except (FileNotFoundError, KeyError):
        return False
    required = (
        "redirect_uri",
        "cookie_secret",
        "client_id",
        "client_secret",
        "server_metadata_url",
    )
    return all(str(auth.get(key, "")).strip() for key in required)


def activate_project(project_name: str) -> None:
    """Switch projects without retaining editor widgets from another project."""
    for key in list(st.session_state):
        del st.session_state[key]
    st.session_state.active_project = project_name
    st.session_state.configuration_page = "project_home"
    st.session_state.visited_configuration_pages = set()


def configuration_signature(project: Path) -> tuple[int, int, int]:
    """Summarize configuration state so prepared downloads cannot go stale."""
    files = [
        path
        for root_name in ("drutes.conf", "bin")
        for path in (project / root_name).rglob("*")
        if path.is_file() and not path.is_symlink()
    ]
    if not files:
        return 0, 0, 0
    statistics = [path.stat() for path in files]
    return (
        len(files),
        sum(item.st_size for item in statistics),
        max(item.st_mtime_ns for item in statistics),
    )


def render_project_gateway(store: ProjectStore, email: str) -> Path | None:
    """Select or create the authenticated user's persistent project."""
    active_name = st.session_state.get("active_project")
    if active_name:
        try:
            project = store.resolve_project(email, str(active_name))
        except (FileNotFoundError, ValueError):
            st.session_state.pop("active_project", None)
        else:
            with st.sidebar:
                st.markdown(f"**{st.user.get('name', email)}**")
                st.caption(email)
                st.markdown(f"Project: **{project.name}**")
                if st.button("Switch project", width="stretch"):
                    process = st.session_state.get("model_process")
                    if process is not None and process.poll() is None:
                        st.error("Stop the running simulation before switching projects.")
                    else:
                        st.session_state.pop("active_project", None)
                        st.rerun()
                if st.button("Log out", width="stretch"):
                    process = st.session_state.get("model_process")
                    if process is not None and process.poll() is None:
                        stop_model_process(process)
                    st.logout()
                if st.button("Prepare configuration ZIP", width="stretch"):
                    with st.spinner("Creating configuration archive…"):
                        st.session_state.configuration_archive = (
                            build_configuration_archive(project)
                        )
                        st.session_state.configuration_archive_signature = (
                            configuration_signature(project)
                        )
                archive = st.session_state.get("configuration_archive")
                archive_is_current = st.session_state.get(
                    "configuration_archive_signature"
                ) == configuration_signature(project)
                if archive is not None and archive_is_current:
                    st.download_button(
                        "Download project configuration",
                        data=archive,
                        file_name=f"{project.name}-configuration.zip",
                        mime="application/zip",
                        width="stretch",
                    )
                elif archive is not None:
                    st.session_state.pop("configuration_archive", None)
                    st.session_state.pop("configuration_archive_signature", None)
            return project

    st.markdown('<p class="drutes-kicker">Secure workspace</p>', unsafe_allow_html=True)
    st.title("Your DRUtES projects")
    st.markdown(
        '<p class="drutes-subtitle">Create a new simulation project or reopen '
        'one of your existing projects.</p><div class="drutes-rule"></div>',
        unsafe_allow_html=True,
    )
    account_column, logout_column = st.columns([4, 1])
    with account_column:
        st.write(f"Signed in as **{st.user.get('name', email)}** ({email})")
    with logout_column:
        if st.button("Log out", width="stretch"):
            st.logout()

    with st.container(border=True):
        st.subheader("Create a project")
        project_name = st.text_input(
            "Project name",
            placeholder="My simulation",
            max_chars=64,
        )
        if st.button("Create project", type="primary", width="stretch"):
            try:
                with st.spinner("Copying the default DRUtES configuration…"):
                    project = store.create_project(email, project_name)
            except (OSError, ValueError) as error:
                st.error(f"Project was not created: {error}")
            else:
                activate_project(project.name)
                st.rerun()

    projects = store.list_projects(email)
    if projects:
        with st.container(border=True):
            st.subheader("Existing projects")
            selected = st.selectbox("Project", projects)
            if st.button("Open project", width="stretch"):
                activate_project(selected)
                st.rerun()
    else:
        st.info("You do not have any projects yet.")
    return None


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

    if not authentication_is_configured():
        st.markdown('<p class="drutes-kicker">Authentication setup</p>', unsafe_allow_html=True)
        st.title("Google sign-in is not configured")
        st.info(
            "Copy .streamlit/secrets.toml.example to .streamlit/secrets.toml, "
            "then add your Google OAuth client ID, client secret, and a strong "
            "cookie secret. Restart Streamlit after saving the file."
        )
        return
    if not st.user.is_logged_in:
        st.markdown('<p class="drutes-kicker">Secure workspace</p>', unsafe_allow_html=True)
        st.title("Sign in to DRUtES")
        st.markdown(
            '<p class="drutes-subtitle">Use your Google account to access '
            'your simulation projects.</p><div class="drutes-rule"></div>',
            unsafe_allow_html=True,
        )
        if st.button("Continue with Google", type="primary", width="stretch"):
            st.login()
        return

    email = str(st.user.get("email", "")).strip().lower()
    if not email:
        st.error("Google did not provide an email address for this account.")
        if st.button("Log out"):
            st.logout()
        return
    project_store = ProjectStore(
        PROJECT_ROOT / "user",
        PROJECT_ROOT / "drutes.conf",
        PROJECT_ROOT / "bin" / "drutes",
    )
    workspace_root = render_project_gateway(project_store, email)
    if workspace_root is None:
        return

    if "configuration_page" not in st.session_state:
        st.session_state.configuration_page = "project_home"
    if "visited_configuration_pages" not in st.session_state:
        st.session_state.visited_configuration_pages = {"global"}

    def navigate(page: str) -> None:
        st.session_state.configuration_page = page
        st.session_state.visited_configuration_pages.add(page)

    def save_all(model_type: str) -> None:
        configs = [
            GlobalConfigFile(workspace_root / "drutes.conf" / "global.conf"),
            Mesh1DConfigFile(
                workspace_root / "drutes.conf" / "mesh" / "drumesh1d.conf"
            ),
            SolverConfigFile(workspace_root / "drutes.conf" / "solver.conf"),
        ]
        if model_type == "heat":
            heat = HeatConfigFile(
                workspace_root / "drutes.conf" / "heat" / "heat.conf"
            ).load()
            configs.append(heat)
            if bool(heat["couple_with_richards"].value):
                configs.append(
                    MatrixConfigFile(
                        workspace_root / "drutes.conf" / "water.conf" / "matrix.conf"
                    )
                )
        else:
            matrix = MatrixConfigFile(
                workspace_root / "drutes.conf" / "water.conf" / "matrix.conf"
            )
            configs.append(matrix)
            matrix.load()
            if bool(matrix["root_water_uptake"].value):
                mesh = Mesh1DConfigFile(
                    workspace_root / "drutes.conf" / "mesh" / "drumesh1d.conf"
                ).load()
                configs.append(
                    RootUptakeConfigFile(
                        workspace_root
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
            "Runs this project's bin/drutes from the project directory. Standard output and "
            "errors are displayed below."
        )
        process = st.session_state.get("model_process")
        running = process is not None and process.poll() is None
        render_model_process_monitor()

        if running or not st.button(
            "Run model", type="primary", use_container_width=True
        ):
            return

        executable = workspace_root / "bin" / "drutes"
        if not executable.is_file():
            st.error(f"Model executable was not found: {executable}")
            return

        output_lines: list[str] = []
        try:
            process = subprocess.Popen(
                [str(executable)],
                cwd=workspace_root,
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

    if current_page == "project_home":
        st.markdown('<p class="drutes-kicker">Project workspace</p>', unsafe_allow_html=True)
        st.title(workspace_root.name)
        st.markdown(
            '<p class="drutes-subtitle">Run the saved model again or open '
            'the visual configuration workflow.</p><div class="drutes-rule"></div>',
            unsafe_allow_html=True,
        )
        configure_column, run_column = st.columns(2)
        with configure_column:
            if st.button(
                "Reconfigure project",
                type="primary",
                width="stretch",
            ):
                st.session_state.visited_configuration_pages = {"global"}
                navigate("global")
                st.rerun()
        with run_column:
            if st.button("Run saved configuration", width="stretch"):
                try:
                    saved_model = str(
                        GlobalConfigFile(
                            workspace_root / "drutes.conf" / "global.conf"
                        ).load()["model_type"].value
                    )
                except (OSError, ValueError, KeyError) as error:
                    st.error(f"The saved configuration cannot be opened: {error}")
                else:
                    st.session_state.selected_problem_type = saved_model
                    st.session_state.all_configurations_saved = True
                    navigate("run")
                    st.rerun()
        project_output = workspace_root / "out"
        if project_output.is_dir() and any(project_output.iterdir()):
            st.info(
                "This project contains results from a previous run. Running "
                "again may replace files in its out directory."
            )
    elif current_page == "mesh":
        mesh_page = MeshConfigurationPage(
            workspace_root / "drutes.conf" / "mesh" / "drumesh1d.conf",
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
            workspace_root / "drutes.conf" / "solver.conf",
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
            workspace_root / "drutes.conf" / "heat" / "heat.conf",
            on_return=lambda: navigate("global"),
        )
        heat_page.on_navigate = navigate
        heat_page.on_save_all = save_all
        heat_page.model_type = "heat"
        heat_page.visited_pages = st.session_state.visited_configuration_pages
        heat_page.render()
    elif current_page == "matrix":
        matrix_page = MatrixConfigurationPage(
            workspace_root / "drutes.conf" / "water.conf" / "matrix.conf",
            on_return=lambda: navigate("global"),
        )
        matrix_page.on_navigate = navigate
        matrix_page.on_save_all = save_all
        matrix_page.model_type = str(
            st.session_state.get("selected_problem_type", "RE")
        )
        matrix_page.visited_pages = st.session_state.visited_configuration_pages
        matrix_page.render()
    elif current_page == "root_uptake":
        root_page = RootUptakeConfigurationPage(
            workspace_root / "drutes.conf" / "water.conf" / "root4uptake.conf",
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
        heat_coupled_with_richards = False
        if selected_model == "heat":
            try:
                heat_coupled_with_richards = bool(
                    HeatConfigFile(
                        workspace_root / "drutes.conf" / "heat" / "heat.conf"
                    ).load()["couple_with_richards"].value
                )
            except (OSError, ValueError, KeyError):
                heat_coupled_with_richards = bool(
                    st.session_state.get("heat_coupled_with_richards", False)
                )
            st.session_state.heat_coupled_with_richards = (
                heat_coupled_with_richards
            )
        edit_pages = [
            ("global", "global.conf"),
            ("mesh", "drumesh1D.conf"),
            ("solver", "solver.conf"),
        ]
        if selected_model == "heat":
            edit_pages.append(("heat", "heat.conf"))
            if heat_coupled_with_richards:
                edit_pages.append(("matrix", "matrix.conf"))
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
        if not model_running:
            SimulationLogPage(workspace_root / "out" / "DRUtES.log").render()
            output_directory = workspace_root / "out"
            if output_directory.is_dir():
                st.download_button(
                    "Download all outputs as ZIP",
                    data=build_output_archive(output_directory),
                    file_name="drutes-outputs.zip",
                    mime="application/zip",
                    key="download_all_outputs",
                    width="stretch",
                )
        try:
            record_solver_time = bool(
                SolverConfigFile(
                    workspace_root / "drutes.conf" / "solver.conf"
                ).load()["record_solver_time"].value
            )
        except (OSError, ValueError, KeyError):
            record_solver_time = False
        if record_solver_time and not model_running:
            SolverTimeOutputPage(
                workspace_root / "out" / "solver.time",
                workspace_root / "drutes.conf" / "global.conf",
            ).render()
        if not model_running:
            output_directory = workspace_root / "out"
            global_config_path = workspace_root / "drutes.conf" / "global.conf"
            if selected_model == "heat" and heat_coupled_with_richards:
                heat_tab, richards_tab = st.tabs(
                    ["Heat equation solution", "Richards equation solution"]
                )
                with heat_tab:
                    HeatOutputPage(
                        output_directory, global_config_path
                    ).render()
                with richards_tab:
                    RichardsOutputPage(
                        output_directory, global_config_path
                    ).render()
            elif selected_model == "heat":
                HeatOutputPage(
                    output_directory, global_config_path
                ).render()
            else:
                RichardsOutputPage(
                    output_directory, global_config_path
                ).render()
    else:
        GlobalConfigurationPage(
            workspace_root / "drutes.conf" / "global.conf",
            logo_path,
            on_edit_mesh=lambda: navigate("mesh"),
            on_edit_solver=lambda: navigate("solver"),
            on_edit_model=navigate,
            on_save_all=save_all,
            visited_pages=st.session_state.visited_configuration_pages,
        ).render()


if __name__ == "__main__":
    main()
