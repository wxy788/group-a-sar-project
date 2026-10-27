"""Recover navigation/attitude from packet auxiliary words without reading I/Q.

Field conventions follow Sentinel-1 Space Packet Protocol Data Unit
S1-IF-ASD-PL-0007, section 3.2.3. No interpolation, extrapolation, quaternion
normalisation or rotation-direction convention is applied here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sentinel1decoder import _field_names as fn
from sentinel1decoder.utilities import read_subcommed_data


SOURCE = "Sentinel-1 RAW packet sub-commutated ancillary words; sentinel1decoder.utilities.read_subcommed_data"
POSITION_FIELDS = [fn.X_POS_DECODED, fn.Y_POS_DECODED, fn.Z_POS_DECODED]
VELOCITY_FIELDS = [fn.X_VEL_DECODED, fn.Y_VEL_DECODED, fn.Z_VEL_DECODED]
QUATERNION_FIELDS = [fn.Q0_DECODED, fn.Q1_DECODED, fn.Q2_DECODED, fn.Q3_DECODED]
ANGULAR_RATE_FIELDS = [fn.X_ANG_RATE_DECODED, fn.Y_ANG_RATE_DECODED, fn.Z_ANG_RATE_DECODED]


def as_matlab_navigation(navigation: dict) -> dict:
    """Keep text as MATLAB row chars while oned_as='column' handles time vectors."""
    return {key: np.asarray([[value]]) if isinstance(value, str) and value else value
            for key, value in navigation.items()}


def _column(metadata: pd.DataFrame, name: str) -> str | None:
    """Match parsed names and underscore-separated MATLAB-exported names."""
    normalize = lambda value: "".join(c for c in str(value).lower() if c.isalnum())
    wanted = normalize(name)
    return next((column for column in metadata.columns if normalize(column) == wanted), None)


def _origin(metadata: pd.DataFrame) -> float:
    coarse = _column(metadata, "Coarse Time")
    fine = _column(metadata, "Fine Time")
    if metadata.empty or coarse is None or fine is None:
        return float("nan")
    try:
        origin = float(metadata[coarse].iloc[0]) + float(metadata[fine].iloc[0])
        return origin if np.isfinite(origin) else float("nan")
    except (ValueError, TypeError):
        return float("nan")


def _unsigned_values(series: pd.Series, bits: int, name: str) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float, na_value=np.nan)
    if (not np.isfinite(values).all() or np.any(values < 0) or
            np.any(values >= 2**bits) or np.any(values != np.floor(values))):
        raise ValueError(f"{name} must contain finite unsigned {bits}-bit integers")
    return values.astype(np.uint64)


def _empty_output(origin: float) -> tuple[dict, dict]:
    common = {
        "available": False,
        "record_count": 0,
        "time_gps_s": np.empty(0, dtype=np.float64),
        "time_relative_s": np.empty(0, dtype=np.float64),
        "time_origin_gps_s": origin,
        "time_origin_note": "First source packet Coarse Time + Fine Time" if np.isfinite(origin) else
                            "First packet time unavailable/non-finite; origin and relative times remain NaN",
        "time_system": "GPS",
        "time_epoch": "1980-01-06T00:00:00 GPS",
        "source": SOURCE,
        "source_note": "Decoded onboard auxiliary records; not an external precise-orbit product.",
        "packet_count": 0,
        "subcomm_start_count": 0,
        "complete_subcomm_record_count": 0,
        "incomplete_subcomm_record_count": 0,
        "incomplete_word_sequence_count": 0,
        "discontinuous_packet_sequence_count": 0,
        "duplicate_record_count": 0,
        "packet_continuity_checked": False,
        "reason": "",
        "processing_note": "Complete 1..64 word sequences with consecutive uint32 Space Packet Count only. "
                           "Independent time sorting and exact duplicate removal; no interpolation or extrapolation. "
                           "Relative times are measured from the first source packet (Coarse Time + Fine Time).",
    }
    orbit = dict(common)
    orbit.update({
        "position_ecef_m": np.empty((0, 3), dtype=np.float64),
        "velocity_ecef_mps": np.empty((0, 3), dtype=np.float32),
        "reference_frame": "ECEF",
        "units_note": "position_ecef_m: metres; velocity_ecef_mps: metres/second; times: GPS seconds",
    })
    attitude = dict(common)
    attitude.update({
        "quaternion": np.empty((0, 4), dtype=np.float32),
        "quaternion_reference_frame": "ECI/J2000",
        "quaternion_order": "Q0, Q1, Q2, Q3; Q0 is the scalar component",
        "quaternion_rotation_direction": "Not asserted; preserve the source convention",
        "angular_velocity_body_rad_s": np.empty((0, 3), dtype=np.float32),
        "angular_velocity_reference_frame": "spacecraft body",
        "quaternion_note": "Source components are preserved without normalization or sign changes.",
        "units_note": "quaternion: dimensionless; angular_velocity_body_rad_s: radians/second; times: GPS seconds",
    })
    return orbit, attitude


def _unique_records(decoded: pd.DataFrame, time_field: str, value_fields: list[str],
                    label: str) -> pd.DataFrame:
    records = decoded[[time_field, *value_fields]].sort_values(time_field, kind="stable")
    if not np.isfinite(records.to_numpy()).all():
        raise ValueError(f"Decoded {label} contains non-finite timestamps or values")
    times = records[time_field].to_numpy()
    values = records[value_fields].to_numpy()
    repeated = times[1:] == times[:-1]
    conflicting = repeated & np.any(values[1:] != values[:-1], axis=1)
    if np.any(conflicting):
        index = int(np.flatnonzero(conflicting)[0]) + 1
        raise ValueError(f"Conflicting {label} records at GPS timestamp {times[index]}")
    return records.loc[~records[time_field].duplicated()].reset_index(drop=True)


def decode_navigation(metadata: pd.DataFrame) -> tuple[dict, dict]:
    """Return MATLAB-friendly orbit and attitude structs from full packet metadata.

    Input row order must be the original file packet order. PRI Count gaps do
    not invalidate an auxiliary sequence; only the actual Space Packet Count
    is used for packet continuity. Missing auxiliary inputs/no complete records
    return available=False with fixed-shape empty arrays and a reason.
    """
    if not isinstance(metadata, pd.DataFrame):
        raise TypeError("metadata must be a pandas DataFrame")
    metadata = metadata.reset_index(drop=True)
    origin = _origin(metadata)
    orbit, attitude = _empty_output(origin)
    outputs = (orbit, attitude)
    for output in outputs:
        output["packet_count"] = len(metadata)
    index_column = _column(metadata, fn.SUBCOM_ANC_DATA_WORD_INDEX_DECODED)
    word_column = _column(metadata, fn.SUBCOM_ANC_DATA_WORD_DECODED)
    packet_column = _column(metadata, "Space Packet Count")
    reason = ""
    if metadata.empty:
        reason = "Packet metadata is empty"
    elif index_column is None or word_column is None:
        reason = "Sub-commutated ancillary word/index columns are unavailable"
    elif packet_column is None:
        reason = "Space Packet Count is unavailable; complete packet continuity cannot be verified"
    if reason:
        for output in outputs:
            output["reason"] = reason
        return outputs

    indices = _unsigned_values(metadata[index_column], 8, str(index_column))
    if np.any(indices > 64):
        raise ValueError("Sub-commutated ancillary word indices must be 0..64")
    words = _unsigned_values(metadata[word_column], 16, str(word_column))
    packet_counts = _unsigned_values(metadata[packet_column], 32, str(packet_column))
    starts = np.flatnonzero(indices == 1)
    valid_starts = []
    invalid_words = 0
    invalid_packets = 0
    expected_indices = np.arange(1, 65, dtype=np.uint64)
    for start in starts:
        stop = int(start) + 64
        if stop > len(indices) or not np.array_equal(indices[start:stop], expected_indices):
            invalid_words += 1
            continue
        block_counts = packet_counts[start:stop]
        increments = (block_counts[1:] - block_counts[:-1]) & np.uint64(0xFFFFFFFF)
        if not np.all(increments == 1):
            invalid_packets += 1
            continue
        valid_starts.append(int(start))
    for output in outputs:
        output.update({
            "subcomm_start_count": len(starts),
            "complete_subcomm_record_count": len(valid_starts),
            "incomplete_subcomm_record_count": len(starts) - len(valid_starts),
            "incomplete_word_sequence_count": invalid_words,
            "discontinuous_packet_sequence_count": invalid_packets,
            "packet_continuity_checked": True,
        })
    if not valid_starts:
        for output in outputs:
            output["reason"] = "No complete ancillary 1..64 sequence with consecutive Space Packet Count"
        return outputs
    # The library assumes a zero-based contiguous RangeIndex. Only complete,
    # validated blocks enter it, so no packet gap can be bridged accidentally.
    positions = np.concatenate([np.arange(start, start + 64) for start in valid_starts])
    validated = pd.DataFrame({
        fn.SUBCOM_ANC_DATA_WORD_INDEX_DECODED: indices[positions].astype(np.int64),
        fn.SUBCOM_ANC_DATA_WORD_DECODED: words[positions].astype(np.int64),
    })
    decoded = read_subcommed_data(validated)
    if len(decoded) != len(valid_starts):
        raise RuntimeError("Navigation decoder returned an unexpected number of auxiliary records")
    pvt = _unique_records(decoded, fn.POD_SOLN_DATA_TIMESTAMP_DECODED,
                          POSITION_FIELDS + VELOCITY_FIELDS, "PVT")
    att = _unique_records(decoded, fn.ATTITUDE_DATA_TIMESTAMP_DECODED,
                          QUATERNION_FIELDS + ANGULAR_RATE_FIELDS, "attitude")
    for output, records, time_field in (
        (orbit, pvt, fn.POD_SOLN_DATA_TIMESTAMP_DECODED),
        (attitude, att, fn.ATTITUDE_DATA_TIMESTAMP_DECODED),
    ):
        times = records[time_field].to_numpy(copy=True)
        output.update({
            "available": True,
            "record_count": len(records),
            "time_gps_s": times,
            "time_relative_s": times - origin,
            "duplicate_record_count": len(decoded) - len(records),
        })
    orbit["position_ecef_m"] = pvt[POSITION_FIELDS].to_numpy(copy=True)
    orbit["velocity_ecef_mps"] = pvt[VELOCITY_FIELDS].to_numpy(copy=True)
    attitude["quaternion"] = att[QUATERNION_FIELDS].to_numpy(copy=True)
    attitude["angular_velocity_body_rad_s"] = att[ANGULAR_RATE_FIELDS].to_numpy(copy=True)
    return orbit, attitude
