import sys
from Bio import SeqIO
import re
import time

start_time = time.time()
fa = sys.argv[1]

total_lengths = 0
seq_id = []
chr_length = {}
total_gaps = 0
Gaps = 0

for seq in SeqIO.parse(fa, 'fasta'):
    seq_len = len(seq.seq)
    total_lengths += seq_len
    seq_id.append(seq.id)
    chr_length[seq.id] = seq_len
    n_count = seq.seq.upper().count("N")
    total_gaps += n_count
    if n_count > 0:
        # 计算连续N的数量
        gaps = re.split("N+", str(seq.seq).upper())
        Gaps += len(gaps) - 1

n = len(seq_id)
ave = total_lengths / n if n else 0
sorted_chr = sorted(chr_length.items(), key=lambda item: item[1], reverse=True)
chr_length_sorted = dict(sorted_chr)
largest_id, largest = sorted_chr[0] if sorted_chr else ("", 0)

# 计算累积长度
cumulative_lengths = []
cumsum = 0
for _, length in sorted_chr:
    cumsum += length
    cumulative_lengths.append(cumsum)

def get_nx0(x):
    target = total_lengths * (x / 10)
    for i, cumsum in enumerate(cumulative_lengths):
        if cumsum >= target:
            return sorted_chr[i][1], i + 1
    return 0, 0

N_values = {}
for x in range(5, 11):
    N_values[x] = get_nx0(x)

print(f"stats for {fa}")
print(f"sum = {total_lengths}, n = {n}, ave = {ave:.2f}, largest = {largest}, largest_id = {largest_id}")
for x in range(5, 11):
    nx0, number = N_values[x]
    # 查找对应的chr_id
    chr_id = next((k for k, v in sorted_chr if v == nx0), "Unknown")
    print(f"N{x*10} = {nx0}, number = {number}, Chr name is {chr_id}")
print(f"N_count = {total_gaps}")
print(f"Gaps = {Gaps}")
end_time = time.time()
print(f"Total Running Time is {end_time - start_time:.1f} seconds!")
