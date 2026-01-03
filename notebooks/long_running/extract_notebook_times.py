#!/usr/bin/env python3
"""Extract notebook execution times from GitHub Actions logs."""
import sys
import re
import subprocess

def extract_times(run_id):
    """Extract notebook execution times from GitHub Actions logs."""
    cmd = ['gh', 'run', 'view', str(run_id), '--log']
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error fetching logs: {result.stderr}", file=sys.stderr)
        return []
    
    logs = result.stdout
    times = []
    
    # Pattern to match: path = "notebooks/..." followed by t_elapsed = number
    # They appear on separate lines after "### ✓"
    lines = logs.split('\n')
    current_path = None
    
    for line in lines:
        # Look for path line
        path_match = re.search(r'path = "([^"]+)"', line)
        if path_match:
            current_path = path_match.group(1)
        
        # Look for t_elapsed line
        time_match = re.search(r't_elapsed = ([0-9.]+)', line)
        if time_match and current_path:
            time = float(time_match.group(1))
            # Extract just the filename
            filename = current_path.split('/')[-1] if '/' in current_path else current_path
            times.append((filename, time))
            current_path = None  # Reset after finding the pair
    
    return sorted(times, key=lambda x: x[1], reverse=True)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: extract_notebook_times.py <run_id>", file=sys.stderr)
        sys.exit(1)
    
    run_id = sys.argv[1]
    times = extract_times(run_id)
    
    print(f"{'Notebook':<50} {'Time (s)':>12}")
    print("-" * 64)
    total = 0
    for filename, time in times:
        print(f"{filename:<50} {time:>12.2f}")
        total += time
    print("-" * 64)
    print(f"{'TOTAL':<50} {total:>12.2f}")
    print(f"\n{'Notebooks processed':<50} {len(times):>12}")

