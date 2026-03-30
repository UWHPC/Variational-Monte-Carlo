#!/usr/bin/env python3
"""Read VMC binary output and render with Rerun."""

import argparse
import struct
import sys

import numpy as np
import rerun as rr

# -- Binary format ----------------------------------------------------------
# Header (24 bytes):
#   uint64  num_particles
#   float64 box_length
#   uint64  measure_steps
#
# Per frame (32 + num_particles*24 bytes):
#   float64 local_energy
#   float64 mean_energy
#   float64 standard_error
#   float64 acceptance_rate
#   float64 positions[num_particles * 3]   (x0,y0,z0,x1,y1,z1,...)

HEADER_SIZE = 24  # bytes
SCALARS_PER_FRAME = 4  # local_energy, mean_energy, se, acceptance_rate


def read_header(path: str):
    with open(path, "rb") as f:
        raw = f.read(HEADER_SIZE)
    if len(raw) < HEADER_SIZE:
        sys.exit(f"File too small for header: {path}")
    num_particles, box_length, measure_steps = struct.unpack("<QdQ", raw)
    return num_particles, box_length, measure_steps


def frame_bytes(num_particles: int) -> int:
    return (SCALARS_PER_FRAME + num_particles * 3) * 8


def read_frames(path: str, num_particles: int):
    """Yield (local_energy, mean_energy, se, acceptance_rate, positions) per frame."""
    fb = frame_bytes(num_particles)
    pos_count = num_particles * 3
    with open(path, "rb") as f:
        f.seek(HEADER_SIZE)
        while True:
            raw = f.read(fb)
            if len(raw) < fb:
                break
            vals = np.frombuffer(raw, dtype=np.float64)
            local_e = vals[0]
            mean_e = vals[1]
            se = vals[2]
            ar = vals[3]
            positions = vals[4 : 4 + pos_count].reshape(-1, 3)
            yield local_e, mean_e, se, ar, positions


def prescan(path: str, num_particles: int, stride: int):
    """Pre-scan for global energy extrema and frame count."""
    e_min = float("inf")
    e_max = float("-inf")
    count = 0
    for local_e, *_ in read_frames(path, num_particles):
        if count % stride == 0:
            if local_e < e_min:
                e_min = local_e
            if local_e > e_max:
                e_max = local_e
        count += 1
    return count, e_min, e_max


def energy_to_color(e: float, e_min: float, e_range: float):
    """Map energy to a blue-white-red colormap."""
    if e_range == 0.0:
        t = 0.5
    else:
        t = np.clip((e - e_min) / e_range, 0.0, 1.0)
    # blue (low) -> white (mid) -> red (high)
    if t < 0.5:
        s = t * 2.0
        return (int(s * 255), int(s * 255), 255, 200)
    else:
        s = (t - 0.5) * 2.0
        return (255, int((1.0 - s) * 255), int((1.0 - s) * 255), 200)


def main():
    parser = argparse.ArgumentParser(description="Render VMC binary output with Rerun")
    parser.add_argument("--input", default="output/vmc.bin", help="Path to binary file")
    parser.add_argument("--stride", type=int, default=1, help="Render every Nth frame")
    args = parser.parse_args()

    path = args.input
    stride = max(1, args.stride)

    # -- Read header --------------------------------------------------------
    num_particles, box_length, measure_steps = read_header(path)
    print(f"Particles: {num_particles}  Box: {box_length}  Steps: {measure_steps}")

    fb = frame_bytes(num_particles)
    import os

    file_size = os.path.getsize(path)
    total_frames = (file_size - HEADER_SIZE) // fb
    rendered_frames = (total_frames + stride - 1) // stride
    print(f"Frames: {total_frames}  Stride: {stride}  Rendering: {rendered_frames}")

    # -- Pre-scan for energy normalization ----------------------------------
    print("Pre-scanning for energy range...")
    _, e_min, e_max = prescan(path, num_particles, stride)
    e_range = e_max - e_min
    print(f"Energy range: [{e_min:.4f}, {e_max:.4f}]")

    # -- Initialize Rerun ---------------------------------------------------
    rr.init("VMC_Simulation", spawn=True)

    # Static: coordinate system and bounding box
    rr.log("world", rr.ViewCoordinates.RIGHT_HAND_Z_UP, static=True)
    half = box_length / 2.0
    rr.log(
        "world/box",
        rr.Boxes3D(centers=[half, half, half], sizes=[box_length, box_length, box_length]),
        static=True,
    )

    # -- Per-frame rendering ------------------------------------------------
    print("Rendering frames...")
    frame_idx = 0
    rendered = 0
    for local_e, mean_e, se, ar, positions in read_frames(path, num_particles):
        if frame_idx % stride != 0:
            frame_idx += 1
            continue

        rr.set_time("frame", sequence=frame_idx)

        # Particle positions as 3D points
        # Color by local energy: blue (low) -> white -> red (high)
        r, g, b, a = energy_to_color(local_e, e_min, e_range)
        colors = np.full((num_particles, 4), [r, g, b, a], dtype=np.uint8)
        rr.log(
            "world/electrons",
            rr.Points3D(
                positions=positions,
                colors=colors,
                radii=np.full(num_particles, box_length * 0.004),
            ),
        )

        # Scalars
        rr.log("scalars/local_energy", rr.Scalars(local_e))
        rr.log("scalars/mean_energy", rr.Scalars(mean_e))
        if se != 0.0:
            rr.log("scalars/standard_error", rr.Scalars(se))
        rr.log("scalars/acceptance_rate", rr.Scalars(ar))

        rendered += 1
        if rendered % 500 == 0:
            print(f"  {rendered}/{rendered_frames} frames")

        frame_idx += 1

    print(f"Done. {rendered} frames logged to Rerun.")


if __name__ == "__main__":
    main()
