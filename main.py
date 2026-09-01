#!/usr/bin/env python3
"""
Host-side capture for tc_logger.ino

Prints every line from the Mega. Press ENTER to start recording to a
timestamped CSV; press ENTER again to stop and close the file. Ctrl-C to quit.

Usage:
    pip install pyserial
    python tc_capture.py COM5            # Windows
    python tc_capture.py /dev/ttyACM0    # Linux
    python tc_capture.py /dev/cu.usbmodem14201 -b 115200 -o ./logs
"""


import argparse
import datetime as dt
import os
import sys
import threading

import serial

HEADER = "millis,sec,tempC,tempF"


class Recorder:
    def __init__(self, outdir):
        self.outdir = outdir
        self.fh = None
        self.lock = threading.Lock()
        self.rows = 0

    @property
    def active(self):
        return self.fh is not None

    def toggle(self):
        with self.lock:
            if self.fh is None:
                os.makedirs(self.outdir, exist_ok=True)
                stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
                path = os.path.join(self.outdir, f"tc_{stamp}.csv")
                self.fh = open(path, "w", newline="", encoding="utf-8")
                self.fh.write("host_time," + HEADER + "\n")
                self.rows = 0
                print(f"\n>>> RECORDING -> {path}\n", flush=True)
            else:
                path = self.fh.name
                self.fh.close()
                self.fh = None
                print(f"\n>>> STOPPED ({self.rows} rows) -> {path}\n", flush=True)

    def write(self, line):
        with self.lock:
            if self.fh is None:
                return
            host = dt.datetime.now().isoformat(timespec="milliseconds")
            self.fh.write(f"{host},{line}\n")
            self.fh.flush()
            self.rows += 1

    def close(self):
        with self.lock:
            if self.fh is not None:
                print(f"\n>>> closed {self.fh.name} ({self.rows} rows)")
                self.fh.close()
                self.fh = None


def key_watcher(rec, stop_evt):
    """Blocking stdin reader — ENTER toggles recording."""
    while not stop_evt.is_set():
        try:
            sys.stdin.readline()
        except (EOFError, ValueError):
            return
        if stop_evt.is_set():
            return
        rec.toggle()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("port", nargs="?", default="COM10",
                    help="serial port (default: COM10)")
    ap.add_argument("-b", "--baud", type=int, default=115200)
    ap.add_argument("-o", "--outdir", default="logs")
    args = ap.parse_args()

    rec = Recorder(args.outdir)
    stop_evt = threading.Event()

    with serial.Serial(args.port, args.baud, timeout=1) as ser:
        print(f"Connected to {args.port} @ {args.baud}. ENTER toggles CSV recording, Ctrl-C quits.")

        t = threading.Thread(target=key_watcher, args=(rec, stop_evt), daemon=True)
        t.start()

        try:
            while True:



                raw = ser.readline()
                if not raw:
                    continue
                line = raw.decode("utf-8", errors="replace").strip()
                if not line:
                    continue

                tag = "REC " if rec.active else "    "
                print(tag + line, flush=True)

                # Comments and the repeated header never go into the CSV.
                if line.startswith("#") or line == HEADER:
                    continue
                rec.write(line)
        except KeyboardInterrupt:
            print("\ninterrupted")
        finally:
            stop_evt.set()
            rec.close()


if __name__ == "__main__":
    main()