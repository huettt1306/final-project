import matplotlib.pyplot as plt

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
               'chr20', 'chr21', 'chr22', 'chrX']
counts = [42046, 43789, 35237, 35074, 31966, 30111, 27301, 
                          25637, 20109, 24502, 24008, 23772, 19853, 18263, 
                          17639, 17876, 15560, 15759, 11884, 13054, 8270, 
                          8738, 20035]

plt.figure(figsize=(13, 6))  

bars = plt.bar(chromosomes, counts, color='#66b3ff', edgecolor='#1a66ff', linewidth=1.2)

plt.title('Số SNP sau khi lọc', fontsize=15, pad=20, fontweight='bold')
plt.xlabel('Nhiễm sắc thể', fontsize=12, labelpad=10)
plt.ylabel('Số SNP', fontsize=12, labelpad=10)

plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)

for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height + 0.5, 
             f'{height}', 
             ha='center', va='bottom', fontsize=9)

plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.gca().set_axisbelow(True)  
plt.tight_layout()

plt.savefig('snp_stats_phased.png', dpi=300, bbox_inches='tight')
plt.close() 

print("Đã lưu biểu đồ vào 'chromosome_stats_with_chr_prefix.png'")