import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statistics import mode
from collections import Counter

df = pd.read_csv("mask_ious.csv")
data = df["p_t"]

fig, ax = plt.subplots(figsize=(8, 4))

bin_width = 0.001
bins = np.arange(0.984, 1.0, bin_width)
counts, edges = np.histogram(data, bins=bins)

for i, count in enumerate(counts):
    x = (edges[i] + edges[i+1]) / 2
    for j in range(count):
        ax.plot(x, j + 0.5, 'o', color='steelblue', markersize=8)

ax.set_ylim(0, max(counts) + 1)
ax.set_xticks(bins)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticks([])
ax.set_xlabel("IoU: Comparing delineator's watersheds to TauDEM")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_visible(False)

plt.tight_layout()
plt.show()
#plt.savefig("dotplot.png")

def ascii_dotplot(data, bin_width=1, dot_char='*', label_fmt='{:.1f}'):
    # Bin the data
    binned = [round(x / bin_width) * bin_width for x in data]
    counts = Counter(binned)

    for value in sorted(counts):
        label = label_fmt.format(value)
        print(f"{label:>8} | {dot_char * counts[value]}")

ascii_dotplot(data, bin_width=0.001, label_fmt='{:.6f}')

# Calculate some summary statistics
mean_iou = np.mean(data)
median_iou = np.median(data)
mode_iou = mode(data)
print(f"Mean IoU: {mean_iou:.6f}")
print(f"Median IoU: {median_iou:.6f}")
print(f"Mode IoU: {mode_iou:.6f}")