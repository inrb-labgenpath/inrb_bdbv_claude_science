#!/usr/bin/env python3
"""make_continuation_xml.py - XML file of a continuation: the XML of the chain with ONLY the chain length changed to the states that
remain.
usage: python3 make_continuation_xml.py <XML of the chain> <XML of the continuation> <length of the chain> <state of the checkpoint>
The script compares the two files line by line and STOPS (exit code 1, no XML left) unless exactly one line differs and the two lines
are equal apart from the value of chainLength. It writes <XML of the continuation>.comparison.txt."""
import hashlib, os, re, sys


def main():
    src, dst, steps, ck = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    remain = steps - ck
    if not (0 < remain < steps):
        sys.exit(f"STOP: {remain} states remain (length {steps}, checkpoint {ck})")
    a = open(src, encoding="utf-8").read()
    pat = re.compile(r'(<mcmc\b[^>]*\bchainLength=")(\d+)(")')
    m = pat.findall(a)
    if len(m) != 1 or int(m[0][1]) != steps:
        sys.exit(f"STOP: the XML of the chain does not hold exactly one mcmc element with chainLength=\"{steps}\"")
    b = pat.sub(lambda x: x.group(1) + str(remain) + x.group(3), a, count=1)
    la, lb = a.split("\n"), b.split("\n")
    diff = [(i + 1, x, y) for i, (x, y) in enumerate(zip(la, lb)) if x != y]
    ok = (len(la) == len(lb) and len(diff) == 1 and
          re.sub(r'chainLength="\d+"', 'chainLength=""', diff[0][1]) == re.sub(r'chainLength="\d+"', 'chainLength=""', diff[0][2]) and
          re.findall(r'chainLength="(\d+)"', diff[0][1]) == [str(steps)] and re.findall(r'chainLength="(\d+)"', diff[0][2]) == [str(remain)])
    if not ok:
        sys.exit("STOP: the two XML files do not differ in exactly one line, in the chain length only")
    tmp = dst + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(b)
    # second comparison, of the files as they are on disk
    da, db = open(src, "rb").read().split(b"\n"), open(tmp, "rb").read().split(b"\n")
    d2 = [i + 1 for i, (x, y) in enumerate(zip(da, db)) if x != y]
    if len(da) != len(db) or d2 != [diff[0][0]]:
        os.remove(tmp); sys.exit("STOP: the files on disk do not differ in exactly one line")
    os.replace(tmp, dst)
    with open(dst + ".comparison.txt", "w") as fh:
        fh.write(f"XML of the chain: {os.path.basename(src)}, SHA-256 {hashlib.sha256(open(src, 'rb').read()).hexdigest()}, {len(la)} lines\n"
                 f"XML of the continuation: {os.path.basename(dst)}, SHA-256 {hashlib.sha256(open(dst, 'rb').read()).hexdigest()}, {len(lb)} lines\n"
                 f"lines that differ: 1 (line {diff[0][0]})\n  chain:        {diff[0][1].strip()}\n  continuation: {diff[0][2].strip()}\n"
                 f"the two lines are equal apart from the value of chainLength: yes\n"
                 f"length of the chain {steps}; state of the checkpoint {ck}; chain length of the XML of the continuation {remain}\n")
    print(remain)


if __name__ == "__main__":
    main()
