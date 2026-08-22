from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def test_required_files():
    required = [
        "main.nf",
        "nextflow.config",
        "conf/base.config",
        "conf/local.config",
        "conf/docker.config",
        "conf/slurm.config",
        "conf/azure.config",
        "containers/Dockerfile",
        "assets/samplesheet.csv",
    ]

    for file in required:
        assert (ROOT / file).exists(), f"Missing: {file}"


def test_workflow_modules():
    required = [
        "modules/qc/fastqc.nf",
        "modules/preprocessing/fastp.nf",
        "modules/quantification/salmon.nf",
        "modules/counting/tximport.nf",
        "modules/differential_expression/deseq2.nf",
        "modules/reporting/multiqc.nf",
    ]

    for file in required:
        assert (ROOT / file).exists(), f"Missing module: {file}"


def test_nextflow_profiles():
    for profile in ["local", "docker", "slurm", "azure"]:
        result = subprocess.run(
            ["nextflow", "config", "-profile", profile],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr


def test_provenance_config():
    config = (ROOT / "nextflow.config").read_text()

    for item in [
        "execution_report.html",
        "execution_timeline.html",
        "execution_trace.txt",
        "workflow_dag.html",
    ]:
        assert item in config
