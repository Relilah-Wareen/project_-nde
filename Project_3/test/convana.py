import pandas as pd
import numpy as np

# 读取数据
df = pd.read_csv('../data/convergence_results.csv')
groups = df.groupby(['Method', 'p'])
results = []

for (method, p), group in groups:
    group = group.sort_values('N')
    errors = group['Error'].values
        
    # 计算收敛阶
    ratios = []
    for i in range(len(errors)-1):
        if errors[i+1] == 0 or errors[i] == 0:
            continue
        ratio = errors[i] / errors[i+1]
        order = np.log2(ratio)
        ratios.append(order)
    
    if ratios:
        avg_order = np.mean(ratios)
        results.append({
            'Method': method,
            'p': p,
            'Observed Order': round(avg_order, 2),
            'Theory Order': p  # 理论阶
        })

result_df = pd.DataFrame(results)

# 按方法和p值排序
result_df = result_df.sort_values(['Method', 'p'])
print("收敛阶分析结果：")
print(result_df.to_string(index=False))
