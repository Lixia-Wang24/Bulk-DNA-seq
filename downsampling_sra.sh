#!/bin/bash
# 不下载完整sra，直接获取部分数据
#fasterq-dump SRR33264154 --split-files -X 1000 --outdir ./downsampling
#!/bin/bash
# 不下载完整sra，直接获取部分数据

fasterq-dump SRR33264154 \
  --split-files \
  --maxSpotId 1000 \
  --outdir ./downsampling \
  -e 2 \
  --mem 16G \
  -p

:<<!
fasterq-dump SRR33264154 \
  --split-files \       # 分离双端数据
  -X 1000 \             # 提取前1000个spots
  --outdir ./downsampling \  # 正确输出目录参数
  -e 2 \                # 减少线程数（降低内存）
  -m 16G \              # 显式设置内存上限
  -t /tmp \             # 指定临时目录（使用大容量空间）
  -p \                  # 显示进度
  --skip-technical      # 跳过技术reads

fasterq-dump SRR33264154 \
    --split-files \
    -X 1000 \
    --outdir ./downsampling \
    -e 2 \
    -m 16G \
    -t /tmp \
    -p \
    --skip-technical
!

