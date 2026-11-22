# DRUtES GUI

The GUI uses Google OpenID Connect authentication. Each Google account owns
persistent projects beneath `user/<google-email>/<project-name>/`. A new
project receives a private copy of the default `drutes.conf` tree and
`bin/drutes`; configuration, model execution, results, and downloads all use
that project directory.

## Install

```bash
cd /Users/miguel/drutes-dev
source /Users/miguel/.venv/bin/activate
python -m pip install -e .
```

## Configure Google authentication

1. In Google Cloud Console, create an OAuth 2.0 Client ID with application type
   **Web application**.
2. Add this authorized redirect URI for local use:
   `http://localhost:8501/oauth2callback`.
3. Copy the provided template:

   ```bash
   cp .streamlit/secrets.toml.example .streamlit/secrets.toml
   ```

4. Put the Google client ID and client secret into
   `.streamlit/secrets.toml`. Generate a strong cookie secret with:

   ```bash
   python -c "import secrets; print(secrets.token_urlsafe(32))"
   ```

`secrets.toml` and the complete `user/` workspace are ignored by Git. For a
deployed server, replace `redirect_uri` with the public HTTPS application URL
ending in `/oauth2callback`, and register that exact URI with Google.

## Run

```bash
streamlit run drutes_gui/app.py
```

After Google sign-in, create or reopen a project. The sidebar provides project
switching, logout, and a ZIP download containing the complete `drutes.conf`
tree, the project's `bin/drutes`, and project metadata.

The reusable formatting-preserving configuration API is:

```python
from drutespy.config.global_config import GlobalConfigFile

config = GlobalConfigFile("drutes.conf/global.conf").load()
print(config["dt"].value)
config["dt"].value = 50.0
config.save()
```

## Plotting large results

Select observation points and their individual properties, then choose **From**,
**To**, and **Maximum points in the entire chart** and click **Plot**. The default
and maximum budget is 10,000 samples across all selected curves, not per curve;
lower budgets are available for phones. Spatial profiles use a coordinate range,
and solver timings use the simulation-time range. Empty or reversed ranges are
reported without plotting misleading data.

Only a bounded preview is sent to the browser. Reduction retains original first
and last samples and the minima/maxima of each interior bucket in their original
order. It does not average values, but cannot retain every local feature. For
finer detail, narrow the range and click **Plot** again. Browser zoom alone only
enlarges the existing preview; it does not request more data. Previews are cleared
when the selected curves or their output-file signatures change. Source files
and downloaded data retain their full resolution. Solver timing tables, when
enabled, display the same bounded preview.
