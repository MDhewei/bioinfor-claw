"""
bioinfor-claw  —  FastAPI Backend Server
=========================================
This server bridges the web UI with your local Python skill scripts.

Start it:
    cd /Users/whe3/Documents/bioinfor-claw
    pip install fastapi uvicorn python-multipart anthropic openai
    python server.py

Then open:
    http://localhost:8000
"""

import os
import json
import subprocess
import tempfile
import shutil
import asyncio
from pathlib import Path
from datetime import datetime
from typing import Optional

from fastapi import FastAPI, UploadFile, File, Form, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# ── Config ────────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).parent          # /Users/whe3/Documents/bioinfor-claw
RESULTS_DIR = REPO_ROOT / "web_results"      # where outputs are saved
RESULTS_DIR.mkdir(exist_ok=True)

app = FastAPI(title="bioinfor-claw API", version="1.0")

# Allow the frontend to call this server
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve the frontend static files (index.html, skills-bundle.js)
STATIC_DIR = REPO_ROOT / "web"
if STATIC_DIR.exists():
    app.mount("/app", StaticFiles(directory=str(STATIC_DIR), html=True), name="frontend")


# ── Models ────────────────────────────────────────────────────────────────────
class RunSkillRequest(BaseModel):
    skill_path: str          # e.g. "gene-list-analysis/function-annotation-for-gene-list"
    script_name: str         # e.g. "annotation_for_gene_list.py"
    args: list[str] = []     # extra CLI args
    input_data: Optional[str] = None   # inline text content (e.g. gene list)
    input_filename: Optional[str] = None  # filename hint for inline data

class ChatRequest(BaseModel):
    provider: str            # "anthropic" | "openai" | "google" | "mistral"
    model: str
    api_key: str
    system_prompt: str
    messages: list[dict]


# ── Health check ──────────────────────────────────────────────────────────────
@app.get("/api/health")
async def health():
    return {
        "status": "ok",
        "repo_root": str(REPO_ROOT),
        "python": shutil.which("python3") or shutil.which("python"),
        "timestamp": datetime.now().isoformat(),
    }


# ── List available skills ──────────────────────────────────────────────────────
@app.get("/api/skills")
async def list_skills():
    """Return all skill sets and their sub-skills found in the repo."""
    skills = {}
    skip = {"Assets", "assets", ".git", "__pycache__", "node_modules", "web", "web_results"}

    for set_dir in sorted(REPO_ROOT.iterdir()):
        if not set_dir.is_dir() or set_dir.name in skip or set_dir.name.startswith("."):
            continue
        sub = {}
        for skill_dir in sorted(set_dir.iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith("."):
                continue
            md = skill_dir / "SKILL.md"
            scripts = list((skill_dir / "scripts").glob("*.py")) if (skill_dir / "scripts").exists() else []
            sub[skill_dir.name] = {
                "has_skill_md": md.exists(),
                "scripts": [s.name for s in scripts],
                "has_requirements": (skill_dir / "requirements.txt").exists(),
            }
        if sub:
            skills[set_dir.name] = sub

    return {"skills": skills, "repo_root": str(REPO_ROOT)}


# ── Run a skill script ─────────────────────────────────────────────────────────
@app.post("/api/run")
async def run_skill(req: RunSkillRequest):
    """
    Execute a skill's Python script and return stdout/stderr + output files.
    """
    skill_dir = REPO_ROOT / req.skill_path
    if not skill_dir.exists():
        raise HTTPException(404, f"Skill not found: {req.skill_path}")

    script_path = skill_dir / "scripts" / req.script_name
    if not script_path.exists():
        # Try to find any .py file in scripts/
        scripts = list((skill_dir / "scripts").glob("*.py")) if (skill_dir / "scripts").exists() else []
        if not scripts:
            raise HTTPException(404, f"Script not found: {req.script_name}")
        script_path = scripts[0]

    # Create a per-run temp output directory
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = RESULTS_DIR / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # If inline data was provided, write it to a temp input file
    input_file = None
    if req.input_data:
        fname = req.input_filename or "input.txt"
        input_file = out_dir / fname
        input_file.write_text(req.input_data, encoding="utf-8")

    # Build the command
    python = shutil.which("python3") or shutil.which("python") or "python3"
    cmd = [python, str(script_path)]

    # Auto-inject --input and --output if the script likely supports them
    if input_file:
        cmd += ["--input", str(input_file)]
    cmd += ["--output", str(out_dir)]
    cmd += req.args  # any extra args from the UI

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,          # 5 min max
            cwd=str(skill_dir),   # run from inside the skill folder
        )
    except subprocess.TimeoutExpired:
        return JSONResponse({"error": "Script timed out (>5 min)", "run_id": run_id}, status_code=408)
    except Exception as e:
        return JSONResponse({"error": str(e), "run_id": run_id}, status_code=500)

    # Collect output files
    output_files = []
    for f in out_dir.iterdir():
        if f.name == (req.input_filename or "input.txt"):
            continue
        output_files.append({
            "name": f.name,
            "size": f.stat().st_size,
            "url":  f"/api/results/{run_id}/{f.name}",
            "type": f.suffix.lstrip(".").upper(),
        })

    return {
        "run_id":       run_id,
        "script":       str(script_path.relative_to(REPO_ROOT)),
        "returncode":   result.returncode,
        "stdout":       result.stdout[-8000:],   # last 8K chars
        "stderr":       result.stderr[-3000:],
        "output_files": output_files,
        "success":      result.returncode == 0,
    }


# ── Upload a file for analysis ─────────────────────────────────────────────────
@app.post("/api/upload")
async def upload_file(file: UploadFile = File(...)):
    """Upload a data file; returns a temp path to reference in /api/run."""
    upload_dir = RESULTS_DIR / "uploads"
    upload_dir.mkdir(exist_ok=True)
    dest = upload_dir / file.filename
    content = await file.read()
    dest.write_bytes(content)
    return {
        "filename": file.filename,
        "path":     str(dest),
        "size":     len(content),
        "url":      f"/api/uploads/{file.filename}",
    }


# ── Serve result files ─────────────────────────────────────────────────────────
@app.get("/api/results/{run_id}/{filename}")
async def get_result(run_id: str, filename: str):
    path = RESULTS_DIR / run_id / filename
    if not path.exists():
        raise HTTPException(404, "File not found")
    return FileResponse(str(path), filename=filename)


# ── LLM proxy (keeps API keys server-side) ────────────────────────────────────
@app.post("/api/chat")
async def chat_proxy(req: ChatRequest):
    """
    Optional: proxy LLM calls through the server so API keys stay server-side.
    Set ANTHROPIC_API_KEY / OPENAI_API_KEY as environment variables and
    remove the key field from frontend requests.
    """
    # Use env var if no key provided in request
    api_key = req.api_key or os.environ.get(
        "ANTHROPIC_API_KEY" if req.provider == "anthropic"
        else "OPENAI_API_KEY" if req.provider == "openai"
        else "GOOGLE_API_KEY" if req.provider == "google"
        else "MISTRAL_API_KEY", ""
    )
    if not api_key:
        raise HTTPException(401, "No API key provided")

    try:
        if req.provider == "anthropic":
            import anthropic
            client = anthropic.Anthropic(api_key=api_key)
            response = client.messages.create(
                model=req.model,
                max_tokens=4096,
                system=req.system_prompt,
                messages=req.messages,
            )
            return {"content": response.content[0].text}

        elif req.provider in ("openai", "mistral", "custom"):
            import openai
            base_url = {
                "openai":  "https://api.openai.com/v1",
                "mistral": "https://api.mistral.ai/v1",
            }.get(req.provider)
            client = openai.OpenAI(api_key=api_key, base_url=base_url)
            response = client.chat.completions.create(
                model=req.model,
                max_tokens=4096,
                messages=[{"role":"system","content":req.system_prompt}] + req.messages,
            )
            return {"content": response.choices[0].message.content}

        elif req.provider == "google":
            import google.generativeai as genai
            genai.configure(api_key=api_key)
            model = genai.GenerativeModel(req.model, system_instruction=req.system_prompt)
            history = [{"role": m["role"] if m["role"]!="assistant" else "model",
                        "parts": [m["content"]]} for m in req.messages[:-1]]
            chat = model.start_chat(history=history)
            response = chat.send_message(req.messages[-1]["content"])
            return {"content": response.text}

    except Exception as e:
        raise HTTPException(500, str(e))


# ── Serve frontend at root ─────────────────────────────────────────────────────
@app.get("/")
async def root():
    index = REPO_ROOT / "web" / "index.html"
    if index.exists():
        return FileResponse(str(index))
    return JSONResponse({"message": "bioinfor-claw API running. Place index.html in /web/ folder."})


if __name__ == "__main__":
    import uvicorn
    print("""
╔══════════════════════════════════════════════════╗
║   bioinfor-claw  Backend Server                 ║
║   http://localhost:8000                         ║
║                                                  ║
║   API docs: http://localhost:8000/docs          ║
╚══════════════════════════════════════════════════╝
    """)
    uvicorn.run("server:app", host="0.0.0.0", port=8000, reload=True)
