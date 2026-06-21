import pandas as pd

# 1. 매핑 정보 읽기 (ENSP ID -> Gene Symbol)
info = pd.read_csv('9606.protein.info.v12.0.txt', sep='\t', usecols=['#string_protein_id', 'preferred_name'])
mapping = dict(zip(info['#string_protein_id'], info['preferred_name']))

# 2. 네트워크 링크 읽기
# 링크 파일은 'protein1 protein2 combined_score' 구조입니다.
print("Converting STRING ID to Gene Symbol...")
with open('9606.protein.links.v12.0.txt', 'r') as f_in, open('string_network.txt', 'w') as f_out:
    header = f_in.readline() # 헤더 스킵
    for line in f_in:
        p1, p2, score = line.strip().split()
        score = int(score)
        
        # 신뢰도 필터링 (700점 이상 추천, 너무 낮으면 노이즈가 심함)
        if score >= 700:
            gene1 = mapping.get(p1)
            gene2 = mapping.get(p2)
            
            if gene1 and gene2:
                # 분석 코드 형식에 맞게 [Gene1] [Tab] [Gene2] [Tab] [Score/1000]
                f_out.write(f"{gene1}\t{gene2}\t{score/1000}\n")

print("Done! 'string_network.txt' has been created.")