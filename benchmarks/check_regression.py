"""Retired historical threshold gate; retained to explain old workflow references."""
import sys

if __name__ == "__main__":
    sys.exit("Historical absolute throughput thresholds are not valid for current builds. "
             "Use the repeated-study methodology in PERFORMANCE.md; "
             "tests/test_remote_profile.py checks benchmark correctness without a speed threshold.")
