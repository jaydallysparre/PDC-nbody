import sys
import re

f1, f2 = sys.argv[1], sys.argv[2]
eps = float(sys.argv[3]) # epsilon to use

f1_text = []
f2_text = []

with open(f1) as f:
    # ignoring last two lines
    f1_text = f.read().splitlines()[:-2]

with open(f2) as f:
    f2_text = f.read().splitlines()[:-2]

# just get the numbers which are formatted as we expect in the columns
entry_regex = r"([-+]?[0-9]\.[0-9]+[e][-+]?[0-9]+)|[+-]?nan"

diff_found = False

if len(f1_text) != len(f2_text):
    print(f"FILES NOT THE SAME LENGTH: len({f1})={len(f1_text)}, len({f2})={len(f2_text)}")
    diff_found = True

for i in range(len(f1_text)):
    nums1 = re.findall(entry_regex, f1_text[i])
    nums2 = re.findall(entry_regex, f2_text[i])

    for j in range(len(nums1)):
        nans = ("nan", "-nan", "+nan")
        if nums1[j] in nans or nums2[j] in nans:
            if nums1[j] != nums2[j]:
                print(f"DIFF AT LINE {i}")
                diff_found = True
                break
        elif abs(float(nums1[j])-float(nums2[j])) > eps:
            print(f"DIFF AT LINE {i}")
            diff_found = True
            break

sys.exit(diff_found)