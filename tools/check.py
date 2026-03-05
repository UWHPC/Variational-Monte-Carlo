#!/usr/bin/env python3
from __future__ import annotations

import argparse
import difflib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def capture(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(cmd, cwd=str(cwd) if cwd else None, text=True)


def tool_path(name: str) -> str:
    env_key = name.upper().replace("-", "_")
    v = os.environ.get(env_key, "").strip()
    if v:
        return v
    p = shutil.which(name)
    if not p:
        raise FileNotFoundError(f"Required tool not found on PATH: {name}")
    return p


def git_ls_files(patterns: Iterable[str]) -> list[str]:
    out = capture(["git", "ls-files", *patterns], cwd=REPO_ROOT)
    return [line.strip() for line in out.splitlines() if line.strip()]


def ensure_build_dir(build_dir: Path, testing: bool) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    args = [
        "cmake",
        "-S",
        str(REPO_ROOT),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    if testing:
        args.append("-DBUILD_TESTING=ON")

    # Prefer Ninja when present; otherwise let CMake pick a default generator.
    if shutil.which("ninja"):
        args.extend(["-G", "Ninja"])

    run(args, cwd=REPO_ROOT)


def cmake_build(build_dir: Path, target: str | None = None) -> None:
    cmd = ["cmake", "--build", str(build_dir)]
    if target:
        cmd += ["--target", target]
    run(cmd, cwd=REPO_ROOT)


def ctest(build_dir: Path) -> None:
    run(["ctest", "--test-dir", str(build_dir), "--output-on-failure"], cwd=REPO_ROOT)


def format_apply() -> None:
    clang_format = tool_path("clang-format")
    patterns = [
        "*.c",
        "*.cc",
        "*.cpp",
        "*.cxx",
        "*.h",
        "*.hh",
        "*.hpp",
        "*.hxx",
        "*.ipp",
        "*.inl",
    ]
    files = git_ls_files(patterns)
    if not files:
        print("No files to format.")
        return
    run([clang_format, "-i", *files], cwd=REPO_ROOT)
    print(f"clang-format applied to {len(files)} file(s).")


def format_check() -> None:
    clang_format = tool_path("clang-format")
    patterns = [
        "*.c",
        "*.cc",
        "*.cpp",
        "*.cxx",
        "*.h",
        "*.hh",
        "*.hpp",
        "*.hxx",
        "*.ipp",
        "*.inl",
    ]
    files = git_ls_files(patterns)
    if not files:
        print("No files to check.")
        return

    failed = False
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)

        for rel in files:
            src = REPO_ROOT / rel
            dst = tmpdir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

            run([clang_format, "-i", str(dst)], cwd=REPO_ROOT)

            a = src.read_text(encoding="utf-8", errors="replace").splitlines(
                keepends=True
            )
            b = dst.read_text(encoding="utf-8", errors="replace").splitlines(
                keepends=True
            )
            if a != b:
                failed = True
                print(f"Formatting required: {rel}")
                diff = difflib.unified_diff(a, b, fromfile=rel, tofile=rel)
                sys.stdout.writelines(diff)
                print()

    if failed:
        raise SystemExit("clang-format check failed. Run: python tools/check.py format")

    print("clang-format check passed.")


def pragma_once_check() -> None:
    headers = git_ls_files(["*.h", "*.hh", "*.hpp", "*.hxx"])
    failed = False

    for rel in headers:
        p = REPO_ROOT / rel
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()

        # Look for #pragma once near the top, allowing leading blanks and comments.
        found = False
        for line in lines[:80]:
            s = line.strip()
            if not s:
                continue
            if s.startswith("//") or s.startswith("/*") or s.startswith("*"):
                continue
            if s.startswith("#pragma once"):
                found = True
            break

        if not found:
            print(f"Missing #pragma once: {rel}")
            failed = True

    if failed:
        raise SystemExit("pragma once check failed.")

    print("pragma once check passed.")


def tidy_check(build_dir: Path) -> None:
    clang_tidy = tool_path("clang-tidy")

    ccdb = build_dir / "compile_commands.json"
    if not ccdb.exists():
        raise SystemExit(
            f"Missing {ccdb}. Configure first: python tools/check.py configure --build-dir {build_dir}"
        )

    files = git_ls_files(["*.cc", "*.cpp", "*.cxx", "*.h", "*.hpp", "*.hxx"])
    if not files:
        print("No files to lint.")
        return

    failed = False
    for rel in files:
        print(f"clang-tidy: {rel}")
        try:
            run([clang_tidy, "-p", str(build_dir), str(REPO_ROOT / rel)], cwd=REPO_ROOT)
        except subprocess.CalledProcessError:
            failed = True

    if failed:
        raise SystemExit("clang-tidy failed.")

    print("clang-tidy passed.")


def cmd_all(build_dir: Path) -> None:
    format_check()
    pragma_once_check()
    ensure_build_dir(build_dir, testing=False)
    cmake_build(build_dir, target=None)
    test_build_dir = REPO_ROOT / "build-tests"
    ensure_build_dir(test_build_dir, testing=True)
    cmake_build(test_build_dir, target=None)
    ctest(test_build_dir)
    tidy_check(build_dir)


def main() -> int:
    ap = argparse.ArgumentParser(description="Cross-platform repo checks")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("format")
    sub.add_parser("format-check")
    sub.add_parser("pragma-once-check")

    ap_cfg = sub.add_parser("configure")
    ap_cfg.add_argument("--build-dir", default="build", help="Build directory")
    ap_cfg.add_argument("--testing", action="store_true", help="Enable BUILD_TESTING")

    ap_build = sub.add_parser("build")
    ap_build.add_argument("--build-dir", default="build", help="Build directory")
    ap_build.add_argument("--target", default=None, help="Optional CMake target")

    ap_test = sub.add_parser("test")
    ap_test.add_argument(
        "--build-dir", default="build-tests", help="Test build directory"
    )

    ap_tidy = sub.add_parser("tidy-check")
    ap_tidy.add_argument(
        "--build-dir",
        default="build",
        help="Build dir containing compile_commands.json",
    )

    ap_all = sub.add_parser("all")
    ap_all.add_argument(
        "--build-dir", default="build", help="Build dir for compile + tidy"
    )

    args = ap.parse_args()

    try:
        if args.cmd == "format":
            format_apply()
        elif args.cmd == "format-check":
            format_check()
        elif args.cmd == "pragma-once-check":
            pragma_once_check()
        elif args.cmd == "configure":
            ensure_build_dir(REPO_ROOT / args.build_dir, testing=bool(args.testing))
        elif args.cmd == "build":
            cmake_build(REPO_ROOT / args.build_dir, target=args.target)
        elif args.cmd == "test":
            ensure_build_dir(REPO_ROOT / args.build_dir, testing=True)
            cmake_build(REPO_ROOT / args.build_dir, target=None)
            ctest(REPO_ROOT / args.build_dir)
        elif args.cmd == "tidy-check":
            ensure_build_dir(REPO_ROOT / args.build_dir, testing=False)
            tidy_check(REPO_ROOT / args.build_dir)
        elif args.cmd == "all":
            cmd_all(REPO_ROOT / args.build_dir)
        else:
            ap.error("Unknown command")
        return 0
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
