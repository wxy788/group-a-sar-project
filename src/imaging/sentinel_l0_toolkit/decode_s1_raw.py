#!/usr/bin/env python3
"""Stream Sentinel-1 RAW ECHO chunks onto an image PRI axis in MAT v7.3.

Requires sentinel1decoder==2.1.0, numpy, pandas, h5py and hdf5storage. The .dat
file is opened read-only. This program uses Level0Decoder directly; it does not
use Level0File's metadata or I/Q caches. No output is written into the SAFE.

MATLAB iq has shape (PRI, range_sample), complex single. Calibration, noise and
missing PRI rows inside the first-to-last ECHO span are zero. StartPulse is
one-based from the first ECHO PRI. Decoding is batched; no full-scene I/Q array is held
in Python. The output is instrument I/Q before focusing/radiometric calibration.
Width is the longest ECHO in the entire measurement; shorter ECHO rows are
right-padded with zeros. timeline.echo_num_samples retains each row's original
complex-sample count (zero for non-ECHO/missing rows), before SWST alignment.
CSV sidecars retain all packet metadata and the acquisition-chunk inventory.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from enum import Enum
from importlib.metadata import version
from pathlib import Path
import platform
import re
import sys
import uuid

import numpy as np
import pandas as pd
import h5py
import hdf5storage
from sentinel1decoder import Level0Decoder
from sentinel1decoder.enums import SignalType
from sentinel1decoder.utilities import range_dec_to_sample_rate
from decode_s1_navigation import decode_navigation, as_matlab_navigation


REQUIRED_DECODER_VERSION = "2.1.0"
DEFAULT_CARRIER_FREQUENCY_HZ = 5.405e9
SPEED_OF_LIGHT_M_S = 299792458.0


def pulse_count_argument(value: str) -> int | None:
    """Allow the MATLAB Inf/end convention as well as a positive row count."""
    if value.lower() in {"inf", "+inf", "infinity", "end", "all"}:
        return None
    try:
        result = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--pulses must be a positive integer or Inf/end") from exc
    if result <= 0:
        raise argparse.ArgumentTypeError("--pulses must be a positive integer or Inf/end")
    return result


def validate_carrier_frequency_hz(value: float) -> float:
    frequency = float(value)
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("--carrier-frequency-hz must be positive and finite")
    wavelength = SPEED_OF_LIGHT_M_S / frequency
    if not np.isfinite(wavelength) or wavelength <= 0:
        raise ValueError("--carrier-frequency-hz must produce a positive finite wavelength")
    return frequency


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="Full path to a RAW measurement .dat; not -annot/-index.dat")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output .mat path outside the input SAFE directory")
    parser.add_argument("--pulses", type=pulse_count_argument, default=None,
                        help="Maximum PRI rows including internal zero rows; Inf/end or omission stops at the last ECHO")
    parser.add_argument("--start-pulse", type=int, default=1,
                        help="One-based PRI index from the first ECHO of the entire measurement (default: 1)")
    parser.add_argument("--carrier-frequency-hz", type=float, default=DEFAULT_CARRIER_FREQUENCY_HZ,
                        help="Positive finite carrier frequency used for radar.fc/lambda (default: 5.405e9 Hz)")
    parser.add_argument("--chunk", type=int, default=None,
                        help="Unsupported legacy option: chunk selection would discard the complete PRI timeline")
    parser.add_argument("--batch-pulses", type=int, default=256,
                        help="Maximum ECHO packets decoded per batch (default: 256)")
    return parser.parse_args()


def enum_code(value):
    return value.value if isinstance(value, Enum) else value


def text_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, Enum):
        return f"{value.name} ({value.value})"
    return str(value)


def find_column(frame: pd.DataFrame, *names: str) -> str | None:
    def normalize(name):
        return re.sub(r"[^a-z0-9]", "", str(name).lower())
    lookup = {normalize(column): column for column in frame.columns}
    for name in names:
        if normalize(name) in lookup:
            return lookup[normalize(name)]
    return None


def scalar_number(frame: pd.DataFrame, *names: str) -> tuple[float, str]:
    column = find_column(frame, *names)
    if column is None:
        return float("nan"), "unavailable; inspect packet_metadata"
    values = pd.to_numeric(frame[column].map(enum_code), errors="coerce").to_numpy(dtype=float)
    if len(values) and np.isfinite(values).all() and np.all(values == values[0]):
        return float(values[0]), str(column)
    return float("nan"), f"{column}: not constant or not numeric; inspect packet_metadata"


def matlab_metadata(frame: pd.DataFrame) -> dict:
    """Export a DataFrame as a struct of columns, retaining the column-name map."""
    flat = frame.reset_index()
    values = {}
    labels = {}
    names = []
    fields = []
    for column in flat.columns:
        field = re.sub(r"[^a-zA-Z0-9_]", "_", str(column)).strip("_")
        if not field or not field[0].isalpha():
            field = "field_" + field
        field = field[:60]
        root = field
        suffix = 1
        while field in values:
            field = root[:55] + "_" + str(suffix)
            suffix += 1
        series = flat[column]
        converted = series.map(enum_code)
        members = {int(x.value): text_value(x) for x in series if isinstance(x, Enum)}
        if members:
            # Optional calibration enums contain pd.NA, so map alone can leave
            # an object column. Store codes + NaN, not one string ref per PRI.
            converted = pd.to_numeric(converted, errors="raise")
        if pd.api.types.is_numeric_dtype(converted.dtype):
            if converted.isna().any():
                values[field] = converted.to_numpy(dtype=float, na_value=np.nan)
            else:
                values[field] = converted.to_numpy()
        else:
            values[field] = np.asarray([text_value(x) for x in series], dtype=object)
        if members:
            # A small codebook avoids millions of repeated MATLAB cell strings.
            labels[field] = {
                "codes": np.asarray(list(members), dtype=np.int64),
                "labels": np.asarray(list(members.values()), dtype=object),
            }
        names.append(str(column))
        fields.append(field)
    return {
        "column_names": np.asarray(names, dtype=object),
        "field_names": np.asarray(fields, dtype=object),
        "values": values,
        "enum_labels": labels,
        "note": "values is in packet order; index fields are zero-based; enum_labels contains codebooks",
    }


def csv_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.reset_index().copy()
    for column in result.columns:
        if any(isinstance(value, Enum) for value in result[column]):
            result[column + " Code"] = result[column].map(enum_code)
            result[column] = result[column].map(text_value)
    return result


def radar_parameters(selected: pd.DataFrame, num_samples: int,
                     carrier_frequency_hz: float = DEFAULT_CARRIER_FREQUENCY_HZ) -> dict:
    carrier_frequency_hz = validate_carrier_frequency_hz(carrier_frequency_hz)
    pri, pri_source = scalar_number(selected, "PRI")
    swst, swst_source = scalar_number(selected, "SWST")
    swl, swl_source = scalar_number(selected, "SWL")
    fs = float("nan")
    fs_source = "unavailable; inspect packet_metadata"
    rg_column = find_column(selected, "Range Decimation", "RGDEC")
    rg_code = float("nan")
    if rg_column is not None:
        rg_values = list(selected[rg_column])
        if len(rg_values) and all(enum_code(x) == enum_code(rg_values[0]) for x in rg_values):
            rg_code = int(enum_code(rg_values[0]))
            fs = float(range_dec_to_sample_rate(rg_values[0]))
            fs_source = f"sentinel1decoder.utilities.range_dec_to_sample_rate({rg_column})"
        else:
            fs_source = f"{rg_column}: varies within selected packets; inspect packet_metadata"
    result = {
        "fc": carrier_frequency_hz,
        "lambda": SPEED_OF_LIGHT_M_S / carrier_frequency_hz,
        "pri_s": pri,
        "prf_hz": 1.0 / pri if np.isfinite(pri) and pri > 0 else np.nan,
        "range_sample_rate_hz": fs,
        "range_decimation_code": rg_code,
        "swst_s": swst,
        "swl_s": swl,
        "parameter_sources": {
            "fc": "Configured carrier frequency in Hz; not inferred from Level-0 packet headers",
            "lambda": "299792458 m/s divided by fc in Hz",
            "pri_s": pri_source,
            "prf_hz": "1 / parsed PRI in seconds",
            "range_sample_rate_hz": fs_source,
            "swst_s": swst_source,
            "swl_s": swl_source,
        },
        "units_note": "fc/PRF/range sample rate: Hz; lambda: metres; PRI/SWST/SWL: seconds",
        "range_axis_note": "Relative sample times only; absolute slant-range timing needs instrument delay corrections",
    }
    if np.isfinite(fs) and fs > 0:
        result["range_sample_time_relative_s"] = np.arange(num_samples, dtype=float) / fs
    if np.isfinite(pri) and pri > 0:
        result["pulse_time_relative_s"] = np.arange(len(selected), dtype=float) * pri
    return result


def numeric_column(frame: pd.DataFrame, *names: str) -> np.ndarray:
    column = find_column(frame, *names)
    if column is None:
        raise ValueError(f"Required metadata field is unavailable: {names[0]}")
    values = pd.to_numeric(frame[column].map(enum_code), errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Metadata {column} contains missing or non-numeric values")
    return values


def echo_sample_counts(echo: pd.DataFrame) -> np.ndarray:
    """Original complex-sample counts from the ECHO packets' uint16 quad counts."""
    nq = numeric_column(echo, "Number of Quads")
    if np.any((nq <= 0) | (nq != np.floor(nq)) | (nq > np.iinfo(np.uint16).max)):
        raise ValueError("ECHO Number of Quads must be a positive uint16 integer")
    return (2 * nq).astype(np.uint32)


def validate_echo_grid(metadata: pd.DataFrame,
                       carrier_frequency_hz: float = DEFAULT_CARRIER_FREQUENCY_HZ) -> tuple[pd.DataFrame, int, dict]:
    """Require one sample spacing/swath/PRI, allowing variable window lengths/origins."""
    signals = numeric_column(metadata, "Signal Type")
    known_signals = {member.value for member in SignalType}
    if not set(signals).issubset(known_signals):
        raise ValueError("Unknown Signal Type code; cannot classify zero rows safely")
    echo = metadata.loc[signals == SignalType.ECHO.value]
    if echo.empty:
        raise ValueError("No ECHO acquisition chunks were found")
    for aliases in (("Range Decimation", "RGDEC"), ("Swath Number", "Swath"), ("PRI",)):
        values = numeric_column(echo, *aliases)
        if not np.all(values == values[0]):
            raise ValueError(f"Incompatible ECHO sampling grids: {aliases[0]} varies across the scene; "
                             "a single raw matrix cannot represent these grids")
    sample_counts = echo_sample_counts(echo)
    swl = numeric_column(echo, "SWL")
    if np.any(swl <= 0):
        raise ValueError("ECHO SWL must be positive and finite")
    pri = numeric_column(metadata, "PRI")
    if pri[0] <= 0 or not np.all(pri == pri[0]):
        raise ValueError("PRI duration varies across packets; a uniform PRI time axis is not valid")
    num_samples = int(sample_counts.max())
    radar = radar_parameters(echo, num_samples, carrier_frequency_hz)
    radar.pop("pulse_time_relative_s", None)
    swst = numeric_column(echo, "SWST")
    radar["swst_reference_s"] = float(swst.min())
    radar["shared_absolute_range_origin"] = bool(np.all(swst == swst[0]))
    radar["range_axis_note"] = (
        "Columns retain each packet's original receive-window samples, without shifting/resampling. "
        "Width is the longest ECHO in the entire measurement; shorter rows have trailing zeros. "
        "Original complex-sample counts are in timeline.echo_num_samples. "
        "Per-row SWST is in timeline.swst_s; differing SWST means different range origins. "
        "Absolute slant-range timing still needs instrument delay corrections."
    )
    return echo, num_samples, radar


def build_timeline(metadata: pd.DataFrame, start: int = 0, pulses: int | None = None,
                   pri_s: float = 1.0, *, image_origin: bool = False) -> tuple[pd.DataFrame, dict]:
    """Unwrap uint32 PRI Count in packet order, then select a complete PRI interval.

    Modulo increments must be in [1, 2**31): duplicate counts, backward resets,
    and ambiguous half-cycle-or-longer jumps are rejected instead of silently
    assigning multiple packets to one row or inventing a huge backwards gap.

    By default start is source-zero-based, preserving the navigation helper API.
    With image_origin=True, start is one-based from the first ECHO and the
    selection ends no later than the last ECHO. Counts include internal gaps.
    """
    if pulses is not None and np.isposinf(pulses):
        pulses = None
    if not isinstance(start, (int, np.integer)) or start < (1 if image_origin else 0):
        origin = "one-based from the first ECHO and must be >= 1" if image_origin else "source-zero-based and must be >= 0"
        raise ValueError(f"start is {origin}")
    if pulses is not None and (not isinstance(pulses, (int, np.integer)) or pulses <= 0):
        raise ValueError("--pulses must be a positive integer or Inf/end")
    if not isinstance(metadata.index, pd.MultiIndex) or metadata.index.nlevels != 2:
        raise ValueError("Unexpected metadata index; expected acquisition chunk and packet number")
    ordered = metadata.sort_index(level=1).copy()
    packet_ids = ordered.index.get_level_values(1).to_numpy(dtype=np.int64)
    if len(ordered) == 0 or np.any(np.diff(packet_ids) <= 0):
        raise ValueError("Packet numbers must be unique and increasing")
    counts_float = numeric_column(ordered, "PRI Count", "PRI_Count", "PRICT")
    if np.any(counts_float < 0) or np.any(counts_float >= 2**32) or np.any(counts_float != np.floor(counts_float)):
        raise ValueError("PRI Count must contain unsigned 32-bit integers")
    counts = counts_float.astype(np.uint64)
    increments = (counts[1:] - counts[:-1]) & np.uint64(0xFFFFFFFF)
    bad = np.flatnonzero((increments == 0) | (increments >= 2**31))
    if len(bad):
        j = int(bad[0])
        raise ValueError(f"PRI Count duplicate/reset or ambiguous jump at packets {packet_ids[j]} -> "
                         f"{packet_ids[j + 1]}: {counts[j]} -> {counts[j + 1]}")
    global_rows = np.concatenate((np.zeros(1, dtype=np.uint64), np.cumsum(increments, dtype=np.uint64)))
    total = int(global_rows[-1]) + 1
    signals = numeric_column(ordered, "Signal Type").astype(np.int16)
    echo_positions = np.flatnonzero(signals == SignalType.ECHO.value)
    image_first = int(global_rows[echo_positions[0]]) if len(echo_positions) else None
    image_last = int(global_rows[echo_positions[-1]]) if len(echo_positions) else None
    if image_origin:
        if image_first is None:
            raise ValueError("No ECHO packets; an image PRI origin cannot be defined")
        image_length = image_last - image_first + 1
        if start > image_length:
            raise ValueError(f"--start-pulse {start} exceeds image.Na_all {image_length}")
        start = image_first + start - 1
        stop = image_last + 1
    else:
        if start >= total:
            raise ValueError(f"source-zero start {start} exceeds the full source PRI axis length {total}")
        stop = total
    end = stop if pulses is None else min(stop, start + pulses)
    global_axis = np.arange(start, end, dtype=np.int64)
    unwrapped_axis = global_axis.astype(np.uint64) + counts[0]
    selected = (global_rows >= start) & (global_rows < end)
    selected_rows = global_rows[selected].astype(np.int64) - start
    ordered["global_pri_index_zero_based"] = global_rows
    ordered["selected_row_zero_based"] = np.int64(-1)
    ordered.loc[selected, "selected_row_zero_based"] = selected_rows
    ordered["pri_count_unwrapped_uint64"] = global_rows + counts[0]
    timeline = {
        "global_pri_index_zero_based": global_axis,
        "pri_count_uint32": (unwrapped_axis & np.uint64(0xFFFFFFFF)).astype(np.uint32),
        "pri_count_unwrapped_uint64": unwrapped_axis,
        "time_relative_s": global_axis.astype(float) * pri_s,
        "time_from_selection_start_s": np.arange(end - start, dtype=float) * pri_s,
        "packet_index_zero_based": np.full(end - start, -1, dtype=np.int64),
        "chunk_id": np.full(end - start, -1, dtype=np.int64),
        "signal_type_code": np.full(end - start, -1, dtype=np.int16),
        "echo_mask": np.zeros(end - start, dtype=bool),
        "calibration_mask": np.zeros(end - start, dtype=bool),
        "noise_mask": np.zeros(end - start, dtype=bool),
        "missing_pri_mask": np.ones(end - start, dtype=bool),
        "swst_s": np.full(end - start, np.nan),
        "echo_num_samples": np.zeros(end - start, dtype=np.uint32),
        "echo_num_samples_note": "Original ECHO complex-sample count (2 * Number of Quads) before "
                                 "tail padding and SWST alignment; zero for calibration/noise/missing PRI. "
                                 "Before alignment, valid samples occupy columns 1:echo_num_samples(row). "
                                 "After SWST alignment this remains the original count, not the last valid column.",
        "total_pris": total,
        "total_pris_note": "Full first-to-last SOURCE packet PRI span, including leading/trailing non-ECHO; see image.Na_all for image span",
        "selected_pris": end - start,
        "note": "One row per selected PRI. global_pri_index_zero_based and time_relative_s remain anchored "
                "to the FIRST SOURCE packet (not UTC), including in image-origin selections. "
                "selection time starts at zero. Missing packet/chunk/signal IDs are -1. "
                "Non-ECHO and missing-PRI I/Q rows are zero.",
        "unwrap_note": "uint32 serial-count differences in [1, 2^31); duplicate/reset/ambiguous jumps rejected",
    }
    if image_first is not None:
        timeline["image_pri_index_one_based"] = global_axis - image_first + 1
        timeline["time_from_image_start_s"] = (global_axis - image_first).astype(float) * pri_s
        ordered["image_pri_index_one_based"] = global_rows.astype(np.int64) - image_first + 1
    timeline["packet_index_zero_based"][selected_rows] = packet_ids[selected]
    timeline["chunk_id"][selected_rows] = ordered.index.get_level_values(0).to_numpy(dtype=np.int64)[selected]
    timeline["signal_type_code"][selected_rows] = signals[selected]
    timeline["missing_pri_mask"][selected_rows] = False
    timeline["echo_mask"][selected_rows] = signals[selected] == SignalType.ECHO.value
    timeline["noise_mask"][selected_rows] = signals[selected] == SignalType.NOISE.value
    calibration_codes = [member.value for member in SignalType if member.value >= 8]
    timeline["calibration_mask"][selected_rows] = np.isin(signals[selected], calibration_codes)
    selected_echo = selected & (signals == SignalType.ECHO.value)
    echo_rows = global_rows[selected_echo].astype(np.int64) - start
    timeline["echo_num_samples"][echo_rows] = echo_sample_counts(ordered.loc[selected_echo])
    swst_column = find_column(ordered, "SWST")
    if swst_column is not None:
        timeline["swst_s"][selected_rows] = numeric_column(ordered, "SWST")[selected]
    return ordered, timeline


def image_extent(metadata: pd.DataFrame, selected_pris: int | None = None) -> dict:
    """Describe the full image using source PRI indexes from build_timeline.

    This is independent of selected_row_zero_based and therefore of the slice.
    Na_all includes interior calibration/noise packets and PRIs without an ISP.
    """
    signals = numeric_column(metadata, "Signal Type")
    echo = metadata.loc[signals == SignalType.ECHO.value].sort_index(level=1)
    if echo.empty:
        raise ValueError("No ECHO packets; image extent is unavailable")
    first = int(echo["global_pri_index_zero_based"].iloc[0])
    last = int(echo["global_pri_index_zero_based"].iloc[-1])
    counts = numeric_column(echo, "PRI Count", "PRI_Count", "PRICT")
    result = {
        "Na_all": last - first + 1,
        "first_source_pri_index_zero_based": first,
        "last_source_pri_index_zero_based": last,
        "first_echo_packet_index_zero_based": int(echo.index.get_level_values(1)[0]),
        "last_echo_packet_index_zero_based": int(echo.index.get_level_values(1)[-1]),
        "first_pri_count_uint32": np.uint32(counts[0]),
        "last_pri_count_uint32": np.uint32(counts[-1]),
        "note": "Na_all is the inclusive unfolded PRI span from the FIRST to LAST ECHO in the entire file, "
                "independent of selection. Internal calibration/noise/missing-ISP rows count; leading/trailing "
                "non-ECHO source PRIs do not. StartPulse=1 denotes the first ECHO PRI.",
    }
    if selected_pris is not None:
        result["Na"] = int(selected_pris)
    return result


def chunk_inventory(metadata: pd.DataFrame) -> pd.DataFrame:
    signal_column = find_column(metadata, "Signal Type")
    rows = []
    for chunk_id, group in metadata.groupby(level=0, sort=True):
        row = {
            "chunk_id": int(chunk_id),
            "num_packets": len(group),
            "selected_packets": int((group["selected_row_zero_based"] >= 0).sum()),
            "first_packet_index": int(group.index.get_level_values(1)[0]),
            "last_packet_index": int(group.index.get_level_values(1)[-1]),
            "first_global_pri_index_zero_based": int(group["global_pri_index_zero_based"].iloc[0]),
            "last_global_pri_index_zero_based": int(group["global_pri_index_zero_based"].iloc[-1]),
            "signal_type": text_value(group[signal_column].iloc[0]),
            "signal_type_code": int(enum_code(group[signal_column].iloc[0])),
        }
        for label, aliases in {
            "swath_number": ("Swath Number", "Swath"),
            "num_quads": ("Number of Quads",), "baq_mode": ("BAQ Mode",),
            "pri_s": ("PRI",), "swst_s": ("SWST",),
        }.items():
            column = find_column(group, *aliases)
            row[label] = enum_code(group[column].iloc[0]) if column is not None else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def write_streamed_iq(output: Path, decoder, metadata: pd.DataFrame,
                      num_pris: int, num_samples: int, batch_pulses: int = 256) -> None:
    """Write original ECHO samples; leave shorter-row tails and non-ECHO rows zero."""
    if batch_pulses <= 0:
        raise ValueError("--batch-pulses must be positive")
    signal_column = find_column(metadata, "Signal Type")
    selected_echo = metadata.loc[(metadata["selected_row_zero_based"] >= 0) &
                                 (metadata[signal_column].map(enum_code) == SignalType.ECHO.value)]
    complex_dtype = np.dtype([("real", "<f4"), ("imag", "<f4")])
    with h5py.File(output, "r+") as handle:
        iq_dataset = handle.create_dataset(
            "iq", shape=(num_samples, num_pris), dtype=complex_dtype,
            chunks=(min(num_samples, 2048), min(num_pris, 64)),
            fillvalue=np.zeros((), dtype=complex_dtype),
        )
        iq_dataset.attrs["MATLAB_class"] = np.bytes_("single")
        for chunk_id, chunk in selected_echo.groupby(level=0, sort=True):
            # Level0Decoder starts a new acquisition chunk whenever NQ changes.
            chunk_lengths = echo_sample_counts(chunk)
            chunk_samples = int(chunk_lengths[0])
            if not np.all(chunk_lengths == chunk_samples) or chunk_samples > num_samples:
                raise ValueError(f"ECHO chunk {chunk_id} has inconsistent sample counts or exceeds matrix width")
            print(f"Decoding ECHO chunk {chunk_id}: {len(chunk)} selected packets", flush=True)
            for offset in range(0, len(chunk), batch_pulses):
                batch = chunk.iloc[offset:offset + batch_pulses].sort_index(level=1)
                iq = decoder.decode_packets(batch, batch_size=batch_pulses)
                if iq.shape != (len(batch), chunk_samples) or iq.dtype != np.complex64:
                    raise ValueError(f"Unexpected decoded output: shape={iq.shape}, dtype={iq.dtype}")
                if not np.isfinite(iq).all():
                    raise ValueError("Decoded I/Q contains non-finite samples")
                rows = batch["selected_row_zero_based"].to_numpy(dtype=np.int64)
                # Split sparse PRI positions into runs so gaps remain fill-value zero.
                breaks = np.r_[0, np.flatnonzero(np.diff(rows) != 1) + 1, len(rows)]
                for first, last in zip(breaks[:-1], breaks[1:]):
                    block = np.ascontiguousarray(iq[first:last].T)
                    iq_dataset[:chunk_samples, int(rows[first]):int(rows[last - 1]) + 1] = block.view(complex_dtype)
        handle.flush()


def run(args: argparse.Namespace) -> Path:
    decoder_version = version("sentinel1decoder")
    if decoder_version != REQUIRED_DECODER_VERSION:
        raise ValueError(f"Expected sentinel1decoder=={REQUIRED_DECODER_VERSION}; found {decoder_version}")
    source = args.input.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve()
    if not source.is_file() or source.suffix.lower() != ".dat":
        raise ValueError("--input must be a measurement .dat file")
    if source.name.lower().endswith(("-annot.dat", "-index.dat")):
        raise ValueError("Use the measurement .dat, not annotation/index .dat")
    if output.suffix.lower() != ".mat":
        raise ValueError("--output must end in .mat")
    for ancestor in source.parents:
        if ancestor.name.lower().endswith(".safe") and output.is_relative_to(ancestor):
            raise ValueError("Choose an output directory outside the source SAFE product")
    if output.parent == source.parent:
        raise ValueError("Choose an output directory separate from the measurement directory")
    if args.chunk is not None:
        raise ValueError("--chunk is no longer supported: all ECHO chunks are combined on the complete PRI axis. "
                         "Use --start-pulse and --pulses to select one-based image PRI rows.")
    if (args.pulses is not None and args.pulses <= 0) or args.start_pulse < 1:
        raise ValueError("--pulses must be positive and --start-pulse must be >= 1 from the first ECHO")
    carrier_frequency_hz = validate_carrier_frequency_hz(args.carrier_frequency_hz)
    if args.batch_pulses <= 0:
        raise ValueError("--batch-pulses must be positive")

    initial_stat = source.stat()
    print("Reading packet headers (the library temporarily reads the compressed .dat into memory)...", flush=True)
    decoder = Level0Decoder(str(source))
    metadata = decoder.decode_metadata(return_raw=False)
    _, num_samples, radar = validate_echo_grid(metadata, carrier_frequency_hz)
    metadata, timeline = build_timeline(metadata, args.start_pulse, args.pulses, radar["pri_s"], image_origin=True)
    image = image_extent(metadata, timeline["selected_pris"])
    # Navigation is decoded from all packet headers, independently of the I/Q slice.
    orbit, attitude = decode_navigation(metadata)
    num_pris = timeline["selected_pris"]
    timeline["range_window_offset_samples"] = (
        (timeline["swst_s"] - radar["swst_reference_s"]) * radar["range_sample_rate_hz"]
    )
    timeline["range_window_offset_note"] = (
        "Receive-window origin offset relative to the minimum ECHO SWST, in samples; "
        "raw I/Q columns are not shifted. Missing PRI rows have NaN."
    )
    radar["pulse_time_relative_s"] = timeline["time_relative_s"]
    summary = chunk_inventory(metadata)
    echoes = summary[summary["signal_type_code"] == SignalType.ECHO.value]
    estimated_bytes = num_pris * num_samples * np.dtype(np.complex64).itemsize
    print(f"Source: {source}\nSource PRI span: {timeline['total_pris']} rows; image.Na_all: {image['Na_all']}; "
          f"selected image PRI {args.start_pulse} through {args.start_pulse + num_pris - 1} (one-based, inclusive); "
          f"I/Q shape: ({num_pris}, {num_samples}); logical size {estimated_bytes / 2**30:.3f} GiB", flush=True)

    match = re.search(r"-raw-[a-z]-([hv]{2})-", source.name, re.IGNORECASE)
    polarisation = match.group(1).upper() if match else "unknown; inspect packet_metadata"
    packet_csv = output.with_name(output.stem + "_packets.csv")
    chunk_csv = output.with_name(output.stem + "_chunks.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    processing_info = {
        "source_file": str(source),
        "source_size_bytes": initial_stat.st_size,
        "source_mtime_ns": np.int64(initial_stat.st_mtime_ns),
        "polarisation": polarisation,
        "polarisation_source": "measurement filename",
        "start_pulse_one_based": args.start_pulse,
        "start_source_pri_index_zero_based": int(timeline["global_pri_index_zero_based"][0]),
        "carrier_frequency_hz": carrier_frequency_hz,
        "requested_pulses": args.pulses if args.pulses is not None else -1,
        "total_pris": timeline["total_pris"],
        "selected_pris": num_pris,
        "total_packets": len(metadata),
        "selected_packets": int((metadata["selected_row_zero_based"] >= 0).sum()),
        "echo_chunks_total": len(echoes),
        "echo_chunks_selected": int((echoes["selected_packets"] > 0).sum()),
        "decoded_pulses": int(timeline["echo_mask"].sum()),
        "selected_echo_pris": int(timeline["echo_mask"].sum()),
        "selected_calibration_pris": int(timeline["calibration_mask"].sum()),
        "selected_noise_pris": int(timeline["noise_mask"].sum()),
        "selected_missing_pris": int(timeline["missing_pri_mask"].sum()),
        "selected_zero_pris": int((~timeline["echo_mask"]).sum()),
        "range_samples": num_samples,
        "range_padding_note": "Matrix width is the maximum original ECHO sample count across the entire "
                              "measurement, independent of pulse selection. Shorter ECHO rows are padded "
                              "only at the end; timeline.echo_num_samples records their original lengths.",
        "array_order": "iq(PRI_index, range_sample_index), MATLAB complex single",
        "selection_note": "start_pulse_one_based is relative to the first ECHO PRI; requested_pulses=-1 means "
                          "through the LAST ECHO, inclusive. total_pris describes the source, image.Na_all the image. "
                          "start_source_pri_index_zero_based preserves the source origin used for navigation.",
        "iq_units": "decoded instrument I/Q units; uncalibrated",
        "processing_level": "Level-0 decoded I/Q; no range/azimuth compression or radiometric calibration",
        "decoder_version": decoder_version,
        "numpy_version": version("numpy"),
        "h5py_version": version("h5py"),
        "hdf5storage_version": version("hdf5storage"),
        "mat_format": "7.3 (HDF5), streamed complex single",
        "complete": False,
        "python_version": platform.python_version(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "packet_metadata_csv": str(packet_csv),
        "chunk_summary_csv": str(chunk_csv),
        "decoder_documentation": "https://sentinel1decoder.readthedocs.io/en/latest/api/classes/",
    }
    print(f"ECHO chunks: {processing_info['echo_chunks_selected']}/{len(echoes)}; "
          f"rows: ECHO={processing_info['selected_echo_pris']}, "
          f"calibration={processing_info['selected_calibration_pris']}, "
          f"noise={processing_info['selected_noise_pris']}, missing={processing_info['selected_missing_pris']}", flush=True)
    token = uuid.uuid4().hex
    temporary = [path.with_name(f".{path.stem}.{token}.partial{path.suffix}")
                 for path in (output, packet_csv, chunk_csv)]
    try:
        hdf5storage.savemat(str(temporary[0]), {
            "radar": radar,
            "image": image,
            "packet_metadata": matlab_metadata(metadata),
            "chunk_summary": matlab_metadata(summary),
            "processing_info": processing_info,
            "timeline": timeline,
            "orbit": as_matlab_navigation(orbit),
            "attitude": as_matlab_navigation(attitude),
        }, appendmat=False, fmt="7.3", oned_as="column", store_python_metadata=False)
        write_streamed_iq(temporary[0], decoder, metadata, num_pris, num_samples, args.batch_pulses)
        csv_metadata(metadata).to_csv(temporary[1], index=False, encoding="utf-8-sig")
        summary.to_csv(temporary[2], index=False, encoding="utf-8-sig")
        final_stat = source.stat()
        if (initial_stat.st_size, initial_stat.st_mtime_ns) != (final_stat.st_size, final_stat.st_mtime_ns):
            raise RuntimeError("The source file changed while decoding; output was not committed")
        with h5py.File(temporary[0], "r+") as handle:
            handle["processing_info/complete"][...] = np.uint8(1)
            handle.flush()
        # Each replacement is atomic. Publish MAT last, after both sidecars exist.
        temporary[1].replace(packet_csv)
        temporary[2].replace(chunk_csv)
        temporary[0].replace(output)
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)
    print(f"Saved {output}\nI/Q shape: ({num_pris}, {num_samples}); complex single; "
          f"polarisation: {polarisation}; selected image PRI interval with internal zero rows retained", flush=True)
    return output


if __name__ == "__main__":
    try:
        run(parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        sys.exit(1)
