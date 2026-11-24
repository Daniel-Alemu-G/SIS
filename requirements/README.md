# Pipeline environment

Use Python 3.12 in a project-local virtual environment. System site packages
are enabled only to reuse the installed Python bindings that match GDAL 3.8.4.
Packages missing from the system interpreter are installed into `.venv`.

```bash
python3.12 -m venv --system-site-packages .venv
.venv/bin/python -m pip install -r requirements/pipeline.txt
```

Do not install these requirements into the system interpreter.
