from __future__ import annotations

import tempfile
from pathlib import Path

from scinterop import provenance


class TestProvenance:
    def test_write_conversion_record(self):
        record = provenance.write_conversion_record(
            input_path="/input/test.h5ad",
            input_format="h5ad",
            output_path="/output/test.rds",
            output_format="rds",
            success=True,
            runtime_s=1.5,
            log_dir=tempfile.mkdtemp(),
        )
        assert record["success"] is True
        assert record["input_format"] == "h5ad"
        assert record["output_format"] == "rds"
        assert "timestamp" in record
        assert "version" in record

    def test_write_conversion_record_with_error(self):
        record = provenance.write_conversion_record(
            input_path="/input/test.h5ad",
            input_format="h5ad",
            output_path="/output/test.rds",
            output_format="rds",
            success=False,
            runtime_s=0.5,
            error="Something went wrong",
            log_dir=tempfile.mkdtemp(),
        )
        assert record["success"] is False
        assert record["error"] == "Something went wrong"

    def test_read_provenance_log(self):
        log_dir = tempfile.mkdtemp()
        record1 = provenance.write_conversion_record(
            input_path="/input/1.h5ad",
            input_format="h5ad",
            output_path="/output/1.rds",
            output_format="rds",
            success=True,
            runtime_s=1.0,
            log_dir=log_dir,
        )
        record2 = provenance.write_conversion_record(
            input_path="/input/2.h5ad",
            input_format="h5ad",
            output_path="/output/2.rds",
            output_format="rds",
            success=False,
            runtime_s=0.5,
            error="Failed",
            log_dir=log_dir,
        )

        log_path = Path(log_dir) / "provenance.jsonl"
        records = provenance.read_provenance_log(log_path)

        assert len(records) == 2
        assert records[0]["success"] is True
        assert records[1]["success"] is False

    def test_read_nonexistent_log(self):
        try:
            provenance.read_provenance_log("/nonexistent/provenance.jsonl")
            assert False, "Should have raised"
        except Exception:
            pass
